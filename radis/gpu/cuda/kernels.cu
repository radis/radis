
#include "gpu_cpu_agnostic.h"

extern "C"{

struct initData {
    float v_min;
    float dv;
    int N_v;
    int N_v_FT;
    int N_x_FT;
    float dxG;
    float dxL;
    int N_lines;
    int N_collision_partners;
    float log_c2Mm[16];
};


struct iterData {
    float p;
    float log_2p;
    float hlog_T;
    float log_rT;
    float c2T;
    float N;
    float x[16];
    float l;
    float slit_FWHM;
    float log_wG_min;
    float log_wL_min;
    int N_G;
    int N_L;
    float Q[16];
};

__device__ __constant__ struct initData init_d;
__device__ __constant__ struct iterData iter_d;


__global__ void fillLDM(
    unsigned char* iso,
    float* v0,
    float* da,  // pressure shift  in cm-1/atm
    float* S0,  // initial linestrength
    float* El,
    float* gamma_arr,
    float* na,
    float* S_klm
    ) {

    int N_G = iter_d.N_G;
    int N_L = iter_d.N_L;

    LOOP(blockIdx.x, gridDim.x){
        LOOP(threadIdx.x, blockDim.x){

            int i = threadIdx.x + blockDim.x * blockIdx.x;
            if (i >= init_d.N_lines) CONTINUE;

            //Calc v
            // ... pressure-shift
            float vi = v0[i] + iter_d.p * da[i];
            float ki = (vi - init_d.v_min) / init_d.dv;
            int k0i = (int)ki;
            int k1i = k0i + 1  ;

            if ((k0i < 0) || (k1i >= init_d.N_v)) CONTINUE;

            //Calc wG
            float log_wGi = logf(v0[i]) + init_d.log_c2Mm[iso[i]] + iter_d.hlog_T;
            float li = (log_wGi - iter_d.log_wG_min) / init_d.dxG;
            int l0i = (int)li;
            int l1i = l0i + 1;

            //Calc wL
            float gamma = 0.0;
            for (int j=0; j<init_d.N_collision_partners; j++){
                gamma += iter_d.x[j] * gamma_arr[i + j * init_d.N_lines];
            }
            float log_wLi = logf(gamma) + iter_d.log_2p + na[i] * iter_d.log_rT;
            float mi = (log_wLi - iter_d.log_wL_min) / init_d.dxL;
            int m0i = (int)mi;
            int m1i = m0i + 1;

            //Calc I
            // ... scale linestrengths under equilibrium
            float Si = iter_d.N * iter_d.x[0] * S0[i] * (expf(iter_d.c2T * El[i]) - expf(iter_d.c2T * (El[i] + v0[i]))) / iter_d.Q[iso[i]];

            float avi = ki - (float)k0i;
            float aGi = li - (float)l0i;
            float aLi = mi - (float)m0i;

            float aV00i = (1 - aGi) * (1 - aLi);
            float aV01i = (1 - aGi) * aLi;
            float aV10i = aGi * (1 - aLi);
            float aV11i = aGi * aLi;

            float Sv0i = Si * (1 - avi);
            float Sv1i = Si * avi;

            ADD(&S_klm[k0i * N_G * N_L + l0i * N_L + m0i], Sv0i * aV00i);
            ADD(&S_klm[k0i * N_G * N_L + l0i * N_L + m1i], Sv0i * aV01i);
            ADD(&S_klm[k0i * N_G * N_L + l1i * N_L + m0i], Sv0i * aV10i);
            ADD(&S_klm[k0i * N_G * N_L + l1i * N_L + m1i], Sv0i * aV11i);
            ADD(&S_klm[k1i * N_G * N_L + l0i * N_L + m0i], Sv1i * aV00i);
            ADD(&S_klm[k1i * N_G * N_L + l0i * N_L + m1i], Sv1i * aV01i);
            ADD(&S_klm[k1i * N_G * N_L + l1i * N_L + m0i], Sv1i * aV10i);
            ADD(&S_klm[k1i * N_G * N_L + l1i * N_L + m1i], Sv1i * aV11i);

        }
    }
}


__global__ void applyLineshapes(complex<float>* S_klm_FT, complex<float>* abscoeff) {

    const float pi = 3.141592653589793f;
    const float r4log2 = 0.36067376022224085f; // = 1 / (4 * ln(2))

    LOOP(blockIdx.x, gridDim.x){
        LOOP(threadIdx.x, blockDim.x){
            int k = threadIdx.x + blockDim.x * blockIdx.x;
            if (k >= init_d.N_x_FT) CONTINUE;

            float x = k / (init_d.N_v_FT * init_d.dv);
            float mul = 0.0;
            complex<float> out_complex = 0;

            for (int l = 0; l < iter_d.N_G; l++) {
                float wG = expf(iter_d.log_wG_min + l * init_d.dxG);
                for (int m = 0; m < iter_d.N_L; m++) {
                    int index = k * iter_d.N_G * iter_d.N_L + l * iter_d.N_L + m;
                    float wL = expf(iter_d.log_wL_min + m * init_d.dxL);
                    mul = expf(-r4log2 * powf(pi * x * wG, 2) - pi * x * wL) / init_d.dv / init_d.N_v_FT;
                    out_complex += mul * S_klm_FT[index];
                }
            }
            abscoeff[k].real(out_complex.real());
            abscoeff[k].imag(out_complex.imag());

        }
    }
}

__global__ void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit)  {

    LOOP(blockIdx.x, gridDim.x){
        LOOP(threadIdx.x, blockDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            if (iv < init_d.N_v) {
                transmittance_noslit[iv] = expf(-iter_d.l * abscoeff[iv]);
            }
            else {
                if (iv >= init_d.N_v_FT) CONTINUE;
                transmittance_noslit[iv] = 1.0;
            }
        }
    }
}

__global__ void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT){

    const float pi = 3.141592653589793f;
    const float r4log2 = 0.36067376022224085f; // = 1 / (4 * ln(2))

    LOOP(blockIdx.x, gridDim.x){
        LOOP(threadIdx.x, blockDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            float x = iv / (init_d.N_v_FT * init_d.dv);
            float window = expf(-r4log2 * powf(pi * x * iter_d.slit_FWHM, 2)) / init_d.N_v_FT;

            if (iv >= init_d.N_x_FT) CONTINUE;
            transmittance_FT[iv] = transmittance_noslit_FT[iv] * window;

        }
    }
}

}

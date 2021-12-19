
#ifdef __CUDACC__

#include<cupy/complex.cuh>
#define agnostic_loop(i, max_i)
#define agnostic_add(addr, val) atomicAdd((addr),(val))

#else

#include <complex>
#include <cmath>

#define agnostic_loop(i, max_i) for (i = 0; i < (max_i); i++)
#define agnostic_add(addr, val) *(addr) += (val)

#define __global__
#define __device__
#define __constant__

struct threadIdx_t {int x; int y; int z;} threadIdx;
struct blockDim_t {int x; int y; int z;} blockDim;
struct blockIdx_t {int x; int y; int z;} blockIdx;
struct gridDim_t {int x; int y; int z;} gridDim;

void set_dims(int blockDim_x, int gridDim_x){
    blockDim.x = blockDim_x;
    gridDim.x = gridDim_x;
}
#endif

#ifdef __CUDACC__
extern "C"{
#else
using namespace std;
#endif

const float pi = 3.141592653589793f;
const float r4log2 = 0.36067376022224085f; // = 1 / (4 * ln(2))

//TO-DO: These should really be in gpu.h but cupy fails to load this file if it's included.
struct initData {
    float v_min;
    float v_max;
    float dv;
    int N_v;
    float dxG;
    float dxL;
    int N_total;
    int Max_lines;
    int N_lines;
    int N_points_per_block;
    int N_threads_per_block;
    int N_blocks_per_grid;
    int N_points_per_thread;
    int	N_iterations_per_thread;
    int shared_size_floats;
    float log_c2Mm[16];
};


struct iterData {
    float p;
    float log_2p;
    float hlog_T;
    float log_rT;
    float c2T;
    float N;
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


#ifndef __CUDACC__
void* get_init_ptr(void){
    return &init_d;
}

void* get_iter_ptr(void){
    return &iter_d;
}

float get_T(){
    return expf(2*iter_d.hlog_T);
}
#endif

__global__ void fillLDM(
    unsigned char* iso,
    float* v0,
    float* da,  // pressure shift  in cm-1/atm
    float* S0,  // initial linestrength
    float* El,
    float* gamma,
    float* na,
    float* S_klm
    ) {

    int N_G = iter_d.N_G;
    int N_L = iter_d.N_L;

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            for (int n = 0; n < init_d.N_iterations_per_thread; n++) {

                int i = threadIdx.x + blockDim.x * (n + blockIdx.x * init_d.N_iterations_per_thread);

                if (i < init_d.N_lines) {
                    //Calc v
                    // ... pressure-shift
                    float vi = v0[i] + iter_d.p * da[i];
                    float ki = (vi - init_d.v_min) / init_d.dv;
                    int k0i = (int)ki;
                    int k1i = k0i + 1  ;

                    if ((k0i >= 0) && (k1i < init_d.N_v)) {

                        //Calc wG
                        float log_wGi = logf(v0[i]) + init_d.log_c2Mm[iso[i]] + iter_d.hlog_T;
                        float li = (log_wGi - iter_d.log_wG_min) / init_d.dxG;
                        int l0i = (int)li;
                        int l1i = l0i + 1;

                        //Calc wL
                        float log_wLi = logf(gamma[i]) + iter_d.log_2p + na[i] * iter_d.log_rT;
                        float mi = (log_wLi - iter_d.log_wL_min) / init_d.dxL;
                        int m0i = (int)mi;
                        int m1i = m0i + 1;

                        //Calc I
                        // ... scale linestrengths under equilibrium
                        float Si = iter_d.N * S0[i] * (expf(iter_d.c2T * El[i]) - expf(iter_d.c2T * (El[i] + v0[i]))) / iter_d.Q[iso[i]];

                        float avi = ki - k0i;
                        float aGi = li - l0i;
                        float aLi = mi - m0i;

                        float aV00i = (1 - aGi) * (1 - aLi);
                        float aV01i = (1 - aGi) * aLi;
                        float aV10i = aGi * (1 - aLi);
                        float aV11i = aGi * aLi;

                        float Sv0i = Si * (1 - avi);
                        float Sv1i = Si * avi;

                        agnostic_add(&S_klm[m0i + l0i * N_L + k0i * N_G * N_L], aV00i * Sv0i);
                        agnostic_add(&S_klm[m0i + l0i * N_L + k1i * N_G * N_L], aV00i * Sv1i);
                        agnostic_add(&S_klm[m0i + l1i * N_L + k0i * N_G * N_L], aV01i * Sv0i);
                        agnostic_add(&S_klm[m0i + l1i * N_L + k1i * N_G * N_L], aV01i * Sv1i);
                        agnostic_add(&S_klm[m1i + l0i * N_L + k0i * N_G * N_L], aV10i * Sv0i);
                        agnostic_add(&S_klm[m1i + l0i * N_L + k1i * N_G * N_L], aV10i * Sv1i);
                        agnostic_add(&S_klm[m1i + l1i * N_L + k0i * N_G * N_L], aV11i * Sv0i);
                        agnostic_add(&S_klm[m1i + l1i * N_L + k1i * N_G * N_L], aV11i * Sv1i);
                    }
                }
            }
        }
    }
}


__global__ void applyLineshapes(complex<float>* S_klm_FT, complex<float>* abscoeff) {

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            int k = threadIdx.x + blockDim.x * blockIdx.x;
            if (k < init_d.N_v + 1) {
                float x = k / (2 * init_d.N_v * init_d.dv);
                float mul = 0.0;
                complex<float> out_complex = 0;
                // float out_re = 0.0;
                // float out_im = 0.0;
                float wG, wL;
                int index;

                for (int l = 0; l < iter_d.N_G; l++) {
                    wG = expf(iter_d.log_wG_min + l * init_d.dxG);
                    for (int m = 0; m < iter_d.N_L; m++) {
                        index = k + l * (init_d.N_v+1) + m * iter_d.N_G * (init_d.N_v+1);
                        wL = expf(iter_d.log_wL_min + m * init_d.dxL);
                        mul = expf(-r4log2 * powf(pi * x * wG, 2) - pi * x * wL) / init_d.dv;
                        out_complex += mul * S_klm_FT[index];
                        //out_complex += LDM[index];
                    }
                }
                abscoeff[k].real(out_complex.real());
                abscoeff[k].imag(out_complex.imag());
            }
        }
    }
}

__global__ void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit)  {

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            if (iv < init_d.N_v) {
                transmittance_noslit[iv] = expf(-iter_d.l * abscoeff[iv]);
            }
            else if (iv < init_d.N_v*2){
                transmittance_noslit[iv] = 1.0;
            }
        }
    }
}

__global__ void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT){

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            float x = iv / (2 * init_d.N_v * init_d.dv);
            float window = expf(-r4log2 * powf(pi * x * iter_d.slit_FWHM, 2));

            if (iv < init_d.N_v + 1) {
                transmittance_FT[iv] = transmittance_noslit_FT[iv] * window;
            }
        }
    }
}


#ifdef __CUDACC__
}
#endif

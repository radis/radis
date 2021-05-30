
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
    int N_wG;
    int N_wL;
    int N_wG_x_N_wL;
    int N_total;
    int Max_lines;
    int N_lines;
    int N_points_per_block;
    int N_threads_per_block;
    int N_blocks_per_grid;
    int N_points_per_thread;
    int	N_iterations_per_thread;
    int shared_size_floats;
};


struct iterData {
    float p;
    float log_p;
    float hlog_T;
    float log_rT;
    float c2T;
    float N;
    float l;
    float slit_FWHM;
    float log_wG_min;
    float log_wL_min;
    float log_dwG;
    float log_dwL;
    float log_c2Mm[16];
};

__device__ __constant__ struct initData init_params_d;
__device__ __constant__ struct iterData iter_params_d;


#ifndef __CUDACC__
void set_init_params(void* host_ptr){
    init_params_d = *reinterpret_cast<initData*> (host_ptr);
}

void set_iter_params(void* host_ptr){
    iter_params_d = *reinterpret_cast<iterData*> (host_ptr);
}

float get_T(){
    return expf(2*iter_params_d.hlog_T);
}
#endif

__global__ void fillDLM(
    unsigned char* iso,
    float* v0,
    float* da,
    float* S0,
    float* El,
    float* log_2gs,
    float* na,
    float* DLM,
    float* Q // Q is an array of size max(isotopes_id) + 1 (should move to iter_params_d)
    ) {

    int NwL = init_params_d.N_wL;
    int NwGxNwL = init_params_d.N_wG_x_N_wL;

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            for (int n = 0; n < init_params_d.N_iterations_per_thread; n++) {

                int i = threadIdx.x + blockDim.x * (n + blockIdx.x * init_params_d.N_iterations_per_thread);

                if (i < init_params_d.N_lines) {
                    //Calc v
                    float v_dat = v0[i] + iter_params_d.p * da[i];
                    float iv = (v_dat - init_params_d.v_min) / init_params_d.dv;
                    int iv0 = (int)iv;
                    int iv1 = iv0 + 1  ;

                    if ((iv0 >= 0) && (iv1 < init_params_d.N_v)) {

                        //Calc wG
                        float log_2vMm = logf(v0[i]) + iter_params_d.log_c2Mm[iso[i]];
                        float log_wG_dat = log_2vMm + iter_params_d.hlog_T;
                        float iwG = (log_wG_dat - iter_params_d.log_wG_min) / iter_params_d.log_dwG;
                        int iwG0 = (int)iwG;
                        int iwG1 = iwG0 + 1;

                        //Calc wL
                        float log_wL_dat = log_2gs[i] + iter_params_d.log_p + na[i] * iter_params_d.log_rT;
                        float iwL = (log_wL_dat - iter_params_d.log_wL_min) / iter_params_d.log_dwL;
                        int iwL0 = (int)iwL;
                        int iwL1 = iwL0 + 1;

                        //Calc I
                        float I_add = iter_params_d.N * S0[i] * (expf(iter_params_d.c2T * El[i]) - expf(iter_params_d.c2T * (El[i] + v0[i]))) / Q[iso[i]];

                        float av = iv - iv0;
                        float awG = (iwG - iwG0) * expf((iwG1 - iwG) * iter_params_d.log_dwG);
                        float awL = (iwL - iwL0) * expf((iwL1 - iwL) * iter_params_d.log_dwL);

                        float aV00 = (1 - awG) * (1 - awL);
                        float aV01 = (1 - awG) * awL;
                        float aV10 = awG * (1 - awL);
                        float aV11 = awG * awL;

                        float Iv0 = I_add * (1 - av);
                        float Iv1 = I_add * av;

                        agnostic_add(&DLM[iwL0 + iwG0 * NwL + iv0 * NwGxNwL], aV00 * Iv0);
                        agnostic_add(&DLM[iwL0 + iwG0 * NwL + iv1 * NwGxNwL], aV00 * Iv1);
                        agnostic_add(&DLM[iwL0 + iwG1 * NwL + iv0 * NwGxNwL], aV01 * Iv0);
                        agnostic_add(&DLM[iwL0 + iwG1 * NwL + iv1 * NwGxNwL], aV01 * Iv1);
                        agnostic_add(&DLM[iwL1 + iwG0 * NwL + iv0 * NwGxNwL], aV10 * Iv0);
                        agnostic_add(&DLM[iwL1 + iwG0 * NwL + iv1 * NwGxNwL], aV10 * Iv1);
                        agnostic_add(&DLM[iwL1 + iwG1 * NwL + iv0 * NwGxNwL], aV11 * Iv0);
                        agnostic_add(&DLM[iwL1 + iwG1 * NwL + iv1 * NwGxNwL], aV11 * Iv1);
                    }
                }
            }
        }
    }
}


__global__ void applyLineshapes(complex<float>* DLM, complex<float>* abscoeff) {

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            if (iv < init_params_d.N_v + 1) {
                float x = iv / (2 * init_params_d.N_v * init_params_d.dv);
                float mul = 0.0;
                complex<float> out_complex = 0;
                // float out_re = 0.0;
                // float out_im = 0.0;
                float wG, wL;
                int index;

                for (int iwG = 0; iwG < init_params_d.N_wG; iwG++) {
                    wG = expf(iter_params_d.log_wG_min + iwG * iter_params_d.log_dwG);
                    for (int iwL = 0; iwL < init_params_d.N_wL; iwL++) {
                        //index = iwG + iwL * init_params_d.N_wG + iv * init_params_d.N_wG_x_N_wL;
                        index = iv + iwG * (init_params_d.N_v+1) + iwL * init_params_d.N_wG * (init_params_d.N_v+1);
                        wL = expf(iter_params_d.log_wL_min + iwL * iter_params_d.log_dwL);
                        mul = expf(-r4log2 * powf(pi * x * wG, 2) - pi * x * wL) / init_params_d.dv;
                        out_complex += mul* DLM[index];
                        //out_complex += DLM[index];
                    }
                }
                abscoeff[iv].real(out_complex.real());
                abscoeff[iv].imag(out_complex.imag());
            }
        }
    }
}

__global__ void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit)  {

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            if (iv < init_params_d.N_v) {
                transmittance_noslit[iv] = expf(-iter_params_d.l * abscoeff[iv]);
            }
            else if (iv < init_params_d.N_v*2){
                transmittance_noslit[iv] = 1.0;
            }
        }
    }
}

__global__ void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT){

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
            int iv = threadIdx.x + blockDim.x * blockIdx.x;
            float x = iv / (2 * init_params_d.N_v * init_params_d.dv);
            float window = expf(-r4log2 * powf(pi * x * iter_params_d.slit_FWHM, 2));

            if (iv < init_params_d.N_v + 1) {
                transmittance_FT[iv] = transmittance_noslit_FT[iv] * window;
            }
        }
    }
}


#ifdef __CUDACC__
}
#endif
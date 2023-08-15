
#include "gpu_cpu_agnostic.h"

extern "C"{

struct myStruct {
    int n;
    float x;
    float y;
};


__device__ __constant__ myStruct struct_obj;
__device__ __constant__ int Nt;
__device__ __constant__ int Nf;
__device__ __constant__ float dt;
__device__ __constant__ float wL;


__global__ void add_ints(int* a, int* b, int* c, int N){
    for (int i=0; i<N; i++){
        c[i] = a[i] + b[i];
    }
}


__global__ void applyLineshapes(complex<float>* data) {

    const float pi = 3.141592653589793f;

    LOOP(threadIdx.x, blockDim.x){
        LOOP(blockIdx.x, gridDim.x){
            int k = threadIdx.x + blockDim.x * blockIdx.x;
            if (k < Nf) {
                float x = k / (Nt * dt);
                float mul = 0.0;
                complex<float> out_complex = 0;

                mul = expf(- pi * x * wL);
                out_complex += mul * data[k];

                data[k].real(out_complex.real());
                data[k].imag(out_complex.imag());
            }
        }
    }
}


}

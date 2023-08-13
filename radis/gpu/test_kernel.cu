


#include <cuda/std/complex>
using namespace cuda::std;



#define agnostic_loop(i, max_i)
#define agnostic_add(addr, val) atomicAdd((addr),(val))




extern "C"{


__device__ __constant__ int Nt;
__device__ __constant__ int Nf;
__device__ __constant__ float dt;
__device__ __constant__ float wL;

__global__ void applyLineshapes(complex<float>* data) {

    const float pi = 3.141592653589793f;

    agnostic_loop(threadIdx.x, blockDim.x){
        agnostic_loop(blockIdx.x, gridDim.x){
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

#ifdef __CUDACC__

#include <cuda/std/complex>

using namespace cuda::std;

#define LOOP(i, max_i)
#define ADD(addr, val) atomicAdd((addr),(val))

#else

#include <complex>
#include <cmath>

#define LOOP(i, max_i) for (i = 0; i < (max_i); i++)
#define ADD(addr, val) *(addr) += (val)

#define __global__ __declspec(dllexport)
#define __device__ __declspec(dllexport)
#define __constant__

using namespace std;

extern "C"{

struct threadIdx_t {size_t x; size_t y; size_t z;} threadIdx;
struct blockIdx_t {size_t x; size_t y; size_t z;} blockIdx;

__declspec(dllexport) struct blockDim_t {size_t x; size_t y; size_t z;} blockDim;
__declspec(dllexport) struct gridDim_t {size_t x; size_t y; size_t z;} gridDim;

}

#endif

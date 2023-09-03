#ifdef __CUDACC__

#include <cuda/std/complex>

using namespace cuda::std;

#define LOOP(i, max_i)
#define ADD(addr, val) atomicAdd((addr),(val))
#define CONTINUE return

#else

#include <complex>
#include <cmath>


#if defined _WIN32
#define LIBRARY_API __declspec(dllexport)

#else
#define LIBRARY_API

#endif



#define LOOP(i, max_i) for (i = 0; i < (max_i); i++)
#define ADD(addr, val) *(addr) += (val)
#define CONTINUE continue

#define __global__ LIBRARY_API
#define __device__ LIBRARY_API
#define __constant__

using namespace std;

extern "C"{

struct threadIdx_t {size_t x; size_t y; size_t z;} threadIdx;
struct blockIdx_t {size_t x; size_t y; size_t z;} blockIdx;

LIBRARY_API struct blockDim_t {size_t x; size_t y; size_t z;} blockDim;
LIBRARY_API struct gridDim_t {size_t x; size_t y; size_t z;} gridDim;

}

#endif



#ifdef __CUDACC__

#include<cuda/std/complex>
using namespace cuda::std;

#else

#include<complex>
#include <cmath>

#define __global__
#define __device__
#define __constant__

#endif

using namespace std;

// GPU/CPU functions:
__global__ void fillLDM(unsigned char* iso, float* v0, float* da, float* S0, float* El, float* gamma, float* na, float* S_klm);
__global__ void applyLineshapes(complex<float>* S_klm_FT, complex<float>* abscoeff);
__global__ void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
__global__ void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT);

// CPU only:
#ifndef __CUDACC__

void* get_init_ptr(void);
void* get_iter_ptr(void);
void set_dims(int blockDim_x, int gridDim_x);
float get_T();

#endif


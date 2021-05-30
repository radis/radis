
#include<complex>
using namespace std;

// GPU/CPU functions:
void fillDLM(unsigned char* iso,float* v0,float* da,float* S0, float* El, float* log_2gs,float* na, float* DLM, float* Q);
void applyLineshapes(complex<float>* DLM, complex<float>* abscoeff);
void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT);

// CPU only:
void set_init_params(void* host_ptr);
void set_iter_params(void* host_ptr);
void set_dims(int blockDim_x, int gridDim_x);
float get_T();



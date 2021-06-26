
#include<complex>
using namespace std;

// GPU/CPU functions:
void fillLDM(unsigned char* iso, float* v0, float* da, float* S0, float* El, float* gamma, float* na, float* S_klm);
void applyLineshapes(complex<float>* S_klm_FT, complex<float>* abscoeff);
void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT);

// CPU only:
void* get_init_ptr(void);
void* get_iter_ptr(void);
void set_dims(int blockDim_x, int gridDim_x);
float get_T();



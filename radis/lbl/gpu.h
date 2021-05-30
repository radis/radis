#include<complex>
using namespace std;

void fillDLM(unsigned char* iso,float* v0,float* da,float* S0, float* El, float* log_2gs,float* na, float* DLM, float* Q);
void applyLineshapes(complex<float>* DLM, complex<float>* abscoeff);
void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT);


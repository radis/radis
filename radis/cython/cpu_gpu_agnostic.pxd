from numpy cimport complex64_t
cdef extern from "../gpu.h":
    void fillDLM(unsigned char* iso,float* v0,float* da,float* S0, float* El, float* log_2gs,float* na, float* DLM, float* Q);
    void applyLineshapes(complex64_t* DLM, complex64_t* abscoeff);
    void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
    void applyGaussianSlit(complex64_t* transmittance_noslit_FT, complex64_t* transmittance_FT);

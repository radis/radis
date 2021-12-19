from numpy cimport complex64_t

cdef extern from "../lbl/gpu.h":

    void fillLDM(unsigned char* iso, float* v0, float* da, float* S0, float* El, float* gamma, float* na, float* S_klm);
    void applyLineshapes(complex64_t* S_klm_FT, complex64_t* abscoeff);
    void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
    void applyGaussianSlit(complex64_t* transmittance_noslit_FT, complex64_t* transmittance_FT);
    void* get_init_ptr();
    void* get_iter_ptr();
    void set_dims(int blockDim_x, int gridDim_x);
    float get_T();

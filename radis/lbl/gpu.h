
#include<complex>
using namespace std;

struct initData {
    float v_min;
    float v_max;
    float dv;
    int N_v;
    int N_wG;
    int N_wL;
    int N_wG_x_N_wL;
    int N_total;
    int Max_lines;
    int N_lines;
    int N_points_per_block;
    int N_threads_per_block;
    int N_blocks_per_grid;
    int N_points_per_thread;
    int	N_iterations_per_thread;
    int shared_size_floats;
};


struct iterData {
    float p;
    float log_p;
    float hlog_T;
    float log_rT;
    float c2T;
    float N;
    float l;
    float slit_FWHM;
    float log_wG_min;
    float log_wL_min;
    float log_dwG;
    float log_dwL;
    float log_c2Mm[16];
};

// GPU/CPU functions:
void fillDLM(unsigned char* iso,float* v0,float* da,float* S0, float* El, float* log_2gs,float* na, float* DLM, float* Q);
void applyLineshapes(complex<float>* DLM, complex<float>* abscoeff);
void calcTransmittanceNoslit(float* abscoeff, float* transmittance_noslit);
void applyGaussianSlit(complex<float>* transmittance_noslit_FT, complex<float>* transmittance_FT);

// CPU only:
void set_init_params(initData params);
void set_iter_params(iterData params);
void set_dims(int blockDim_x, int gridDim_x);

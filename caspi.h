#ifndef caspi
#define caspi

//MACROS
#define square(x) ((x)*(x))
#define DEBUG 0
#define PRINT_SIM_MAPS 0
#define WISDOM

//Globals (scenes)
extern  double t_bin;
extern  int N_bin;
extern  int N_bin_half;
extern  int N_sim;
extern  int s_patch;
extern  int r_intra;
extern  int r_inter;
extern  double FWHM;
extern  double d_min;
extern  double d_max;
//#define cmap_name "custom2";

//params
#define SAVE_INTERM 0 //1: save all intermediate intensities, 0: don't
#define PI 3.14159265358979323846
#define EPS 2.220446049250313e-16

//filtering
extern  int skip; //filtering skip size, 1 means dense filtering
extern  int s_intra; //search size within a frame

//itensity-based filtering
extern  int pseudo_int; //1: pseudo intensity, 0: true insensity
extern  int IBF; //1: yes, 0: no
extern  double th_IBF3D; //default
extern  double th_IBF4D; //default
extern  int  max_int_prct; //100: normalize intensity with max

//system
extern int N_cycle; //laser cycle numbers
extern  double c; //light speed

//misc
extern  double d_range; //measurable depth range
extern  double sigma_t; //std for the given FWHM in time doamin
extern  double sigma_f; //stf for the given FWHM in freq domain
extern  double bin_f; //bin size in freq domain
extern  int N_sig_f; //number of singla bins freq domain

extern  double sigma_bin;
extern  int bin_max; //bin size

// sizes
extern int size_x;
extern int size_y;
extern int N_frame;

//Function prototypes
double* normpdf();
double* pseudo_intensity(double *hst_map_set);
double* sum_dim(double *hst_map, int bin_size);
double prctile(double *i_map_hst);
void divide(double *i_map_hst, double prctile);
double *coates_correction(double *hst_map_set);
double *coates_grp(double *hst_map);
double sum_1d(double *hst);
double *coates(double *hst);
double *flux_estimate(double *hst);
double *getNthFrame(double *hst_map_set, int frame);
fftw_complex* fft(double *flux_map_set);
fftw_complex* fft_halfbin(fftw_complex *FLUX_map_set_noisy1D);
void intensities(double *i_map_set_wiener3D, fftw_complex *FLUX_map_set_noisy1D);
fftw_complex* getNthFrame_fftw(fftw_complex *FLUX_map_set_noisy1D, int frame);
double* ifft(fftw_complex *FLUX_map_noisy1D);
void _init();

void ht3D(fftw_complex *FLUX_map_noisy1D, double *i_map_1D, double *MSE_map, double *flux_map_ht3D, double *IBFmask3D);
void extract_patch(fftw_complex* FLUX_map_noisy1D, fftw_complex* FLUX_patch_noisy1D, int start_x, int start_y);
void fft2_patch(fftw_complex *FLUX_map_1D, fftw_complex *FLUX_patch_3D, fftw_complex *SIG_patch_noisy3D, int start_x, int start_y, double *ME, double *MSE, double *S_power);
void sig_and_noisy(fftw_complex* FLUX_patch_noisy3D, fftw_complex *SIG_patch_noisy3D, fftw_complex *NOISY_patch_noisy3D);
void mean_abs(double *ME, double *ME_sqr, double *S_power, fftw_complex *NOISY_patch_noisy3D, fftw_complex *SIG_patch_noisy3D);
void extract_i_patch(int start_y, int start_x, double *i_patch, double *i_map_1D);
void filtering_fftn(double *i_map_1D, fftw_complex *SIG_patch_noisy3D, fftw_complex *SIG_patch_ht3D, int start_x, int start_y);
void filtering(fftw_complex *SIG_patch_ht3D, fftw_complex *SIG_patch_noisy3D, double *I_patch);
void shrink_hard_sig(double noise_th, int *coeff_hard, fftw_complex *SIG_patch_ht3D, fftw_complex *SIG_patch_noisy3D);
void iffn_flux_patch(double *flux_map_ht3D, double *weight_sum_map, fftw_complex *SIG_patch_ht3D, double weight, int start_x, int start_y);
void weight_sum_map_patch(double weight, int start_x, int start_y, double *weight_sum_map);
void flux_map_3D_patch(double weight, int start_x, int start_y, double *flux_map_ht3D, double *flux_patch);
void filtered_flux(double *flux_map_ht3D, double *weight_sum_map);

void wf3D(double *flux_map_wiener3D, double *flux_map_ht3D, fftw_complex *FLUX_map_noisy1D, double *MSE_map);

fftw_complex* fft_wf3D(double *flux_map_ht3D);
void sig_patches(fftw_complex *FLUX_patch_noisy3D, fftw_complex *SIG_patch_noisy3D, fftw_complex *FLUX_patch_ht3D, fftw_complex *SIG_patch_ht3D);
void make_s_power(double *S_power, fftw_complex *SIG_patch_ht3D);

void wiener_srinkage(fftw_complex *SIG_patch_wiener3D, fftw_complex *SIG_patch_noisy3D, double N_power, double *S_power);
void i_map_wiener_frame(int frame, double *i_map_set_wiener3D, double *i_map_wiener3D);

void sim_patches(double* i_map_set_wiener3D, double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set);

void psuedo_patch(int start_x, int start_y, int frame, double *patch, double *i_map_set_wiener3D);

double dist_calc(double *T, double *patch);
void xydistance(double *xyDist_buffer, int x, double *x_buffer, int y, double *y_buffer);
void framedistance(double *frameDist_buffer, int frame, double *frame_buffer);
void distDistance(double *dist_buffer);
int *sub2ind();

void swap(double *a, double *b);
void sort(double array[], double x_buffer[], double y_buffer[], double frame_buffer[], double dist_buffer[], int low, int high);
double partition(double array[], double x_buffer[], double y_buffer[], double frame_buffer[], double dist_buffer[], int low, int high);
void printArray(double *array, int size);

void sim_map(double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set, int x, int y, int frame, double *x_buffer, double *y_buffer, double *frame_buffer);

void ht4D(double *flux_map_set_ht4D, fftw_complex *FLUX_map_set_noisy1D, double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set, double *i_map_set, double *MSE_map_set, double *IBFmask4D_set);

void fft2_ht4D(fftw_complex *FLUX_patchset_noisy3D, fftw_complex *FLUX_map_set_noisy1D, int x_sim, int y_sim, int frame_sim, int sim_idx);

void fft_ht4D(fftw_complex *FLUX_patchset_noisy4D, fftw_complex *FLUX_patchset_noisy3D);
void sig_and_noisy4D(fftw_complex* FLUX_patchset_noisy4D, fftw_complex *SIG_patchset_noisy4D, fftw_complex *NOISY_patchset_noisy4D);
void mean_abs4D(double *ME, double *MSE, double *S_power, fftw_complex *NOISY_patchset_noisy4D, fftw_complex *SIG_patchset_noisy4D);

void extract_i_patch4D(int start_x, int start_y, int frame, double *i_patch, double *i_map_set);

void filtering4D(fftw_complex *SIG_patch_ht4D, fftw_complex *SIG_patchset_noisy4D, double *I_patch);

void shrink_hard_sig4D(double noise_th, int *coeff_hard, fftw_complex *SIG_patch_ht4D, fftw_complex *SIG_patchset_noisy4D);

void iffn_flux_patchset(double *flux_patchset, fftw_complex *SIG_patch_4D);

void weight_sum_map_patch4D(double weight, int sim_idx, int frame, int start_x, int start_y, double *weight_sum_map, double *flux_map_set_ht4D, double *flux_patchset);

void filtered_flux4D(double *flux_map_set_ht4D, double *weight_sum_map_set);


void wf4D(double *flux_map_set_wiener4D, double *flux_map_set_ht4D, fftw_complex *FLUX_map_set_noisy1D, double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set, double *MSE_map_set);

fftw_complex* fft_wf4D(double *flux_map_set_ht4D);
void sig_patches4D(fftw_complex *FLUX_patchset_noisy4D, fftw_complex *SIG_patchset_noisy4D, fftw_complex *FLUX_patchset_ht4D, fftw_complex *SIG_patchset_ht4D);
void make_s_power4D(double *S_power, fftw_complex *SIG_patchset_ht4D);
void wiener_srinkage4D(fftw_complex *SIG_patch_wiener4D, fftw_complex *SIG_patchset_noisy4D, double N_power, double *S_power);

#endif

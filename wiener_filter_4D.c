#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"
#include <time.h>
#include <string.h>
#include <math.h>

void wf4D(double *flux_map_set_wiener4D, double *flux_map_set_ht4D, fftw_complex *FLUX_map_set_noisy1D, double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set, double *MSE_map_set) {
    
    double y_sim;
    double x_sim;
    double frame_sim;
    double N_power;
    double weight;
    
    double* weight_sum_map_set = (double*)calloc((size_x*size_y*N_frame), sizeof(double));
    
    fftw_complex* FLUX_patchset_noisy3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch*s_patch*N_bin_half*N_sim));
    
    fftw_complex* FLUX_patchset_ht3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch*s_patch*N_bin_half*N_sim));
    
    fftw_complex* FLUX_patchset_noisy4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch*s_patch*N_bin_half*N_sim));
    
    fftw_complex* FLUX_patchset_ht4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch*s_patch*N_bin_half*N_sim));
    
    fftw_complex* SIG_patchset_noisy4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f*N_sim));
    
    fftw_complex* SIG_patchset_ht4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f*N_sim));
    
    fftw_complex* SIG_patch_wiener4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f*N_sim));
    
    double *S_power = malloc(sizeof(double) * (s_patch*s_patch*N_sig_f*N_sim));
    
    double* flux_wiener_patchset = malloc(sizeof(double) * (s_patch*s_patch*N_bin*N_sim));
    
    fftw_complex *FLUX_map_set_ht1D = fft_wf4D(flux_map_set_ht4D);
    
    for (int frame=0; frame<N_frame; frame++) {
        for (int x=0; x< (size_x-(s_patch-1)); x=x+skip) {
            for (int y=0; y< (size_y-(s_patch-1)); y=y+skip) {
                
                memset(FLUX_patchset_noisy3D, 0, (s_patch*s_patch*N_bin_half*N_sim) * sizeof(fftw_complex));
                
                memset(FLUX_patchset_ht3D, 0, (s_patch*s_patch*N_bin_half*N_sim) * sizeof(fftw_complex));
                
                for (int sim_idx=0; sim_idx<N_sim; sim_idx++) {
                    
                    y_sim = y_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    x_sim = x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    frame_sim = frame_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    
                    fft2_ht4D(FLUX_patchset_noisy3D, FLUX_map_set_noisy1D, (int)x_sim, (int)y_sim, (int)frame_sim, sim_idx);
                    
                    fft2_ht4D(FLUX_patchset_ht3D, FLUX_map_set_ht1D, (int)x_sim, (int)y_sim, (int)frame_sim, sim_idx);
                }
                
                fft_ht4D(FLUX_patchset_noisy4D, FLUX_patchset_noisy3D);
                
                fft_ht4D(FLUX_patchset_ht4D, FLUX_patchset_ht3D);
                
                sig_patches4D(FLUX_patchset_noisy4D, SIG_patchset_noisy4D, FLUX_patchset_ht4D, SIG_patchset_ht4D);
                
                N_power = MSE_map_set[(x*size_y*N_frame) + (y*N_frame) + frame];
                
                make_s_power4D(S_power, SIG_patchset_ht4D);
                
                wiener_srinkage4D(SIG_patch_wiener4D, SIG_patchset_noisy4D, N_power, S_power);
                
                iffn_flux_patchset(flux_wiener_patchset, SIG_patch_wiener4D);
                
                weight = 1/(N_power + EPS);
                
                for (int sim_idx=0; sim_idx<N_sim; sim_idx++) {
                    y_sim = y_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    x_sim = x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    frame_sim = frame_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    
                    
                    weight_sum_map_patch4D(weight, sim_idx, frame, x, y, weight_sum_map_set, flux_map_set_wiener4D, flux_wiener_patchset);
                }
            }
        }
    }
    filtered_flux4D(flux_map_set_wiener4D, weight_sum_map_set);
}

void wiener_srinkage4D(fftw_complex *SIG_patch_wiener4D, fftw_complex *SIG_patchset_noisy4D, double N_power, double *S_power) {
    
    double wiener_coeff = 0.0;
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                for (int z=0; z<N_sig_f; z++) {
                    wiener_coeff = 1/(1 + N_power/S_power[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim]);
                    
                    SIG_patch_wiener4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = wiener_coeff * SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0];
                    
                    SIG_patch_wiener4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = wiener_coeff * SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1];
                    
                }
            }
        }
    }
}

void make_s_power4D(double *S_power, fftw_complex *SIG_patchset_ht4D) {
    double complex_mag = 0.0;
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                for (int z=0; z<N_sig_f; z++) {
                    complex_mag = pow(sqrt(square(SIG_patchset_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0]) + square(SIG_patchset_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1])), 2);
                    
                    S_power[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim] = complex_mag;
                    if (complex_mag == 0.0) {
                        S_power[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim] = EPS;
                    }
                    
                }
            }
        }
    }
    
}


void sig_patches4D(fftw_complex *FLUX_patchset_noisy4D, fftw_complex *SIG_patchset_noisy4D, fftw_complex *FLUX_patchset_ht4D, fftw_complex *SIG_patchset_ht4D) {
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                for (int z=0; z<N_sig_f; z++) {
                    SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = FLUX_patchset_noisy4D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][0];
                    SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = FLUX_patchset_noisy4D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][1];
                    
                    SIG_patchset_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = FLUX_patchset_ht4D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][0];
                    SIG_patchset_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = FLUX_patchset_ht4D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][1];
                }
            }
        }
    }
}

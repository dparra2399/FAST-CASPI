#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"
#include <time.h>
#include <string.h>
#include <math.h>

void ht4D(double *flux_map_set_ht4D, fftw_complex *FLUX_map_set_noisy1D, double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set, double *i_map_set, double *MSE_map_set, double *IBFmask4D_set) {
    
    double y_sim;
    double x_sim;
    double frame_sim;
    
    double ME;
    double MSE;
    double noise_sigma;
    double noise_scale;
    double noise_th;
    double S_power;
    double N_power;
    double weight;
    
    double *i_patch = malloc(sizeof(double) * s_patch * s_patch);
    
    double *I_patch = malloc(sizeof(double) * s_patch * s_patch);
    
    fftw_complex* FLUX_map_set_ht4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(size_y*size_x*N_bin*N_frame));
    
    double* weight_sum_map_set = (double*)calloc((size_x*size_y*N_frame), sizeof(double)); //map of weight sum
    
    fftw_complex* FLUX_patchset_noisy3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch*s_patch*N_bin_half*N_sim));
    
    fftw_complex* FLUX_patchset_noisy4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch*s_patch*N_bin_half*N_sim));
    
    fftw_complex* SIG_patchset_noisy4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f*N_sim));
    
    fftw_complex* NOISY_patchset_noisy4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*(N_bin_half - N_sig_f)*N_sim));
    
    fftw_complex* SIG_patch_ht4D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f*N_sim));
    
    int* coeff_hard = malloc(sizeof(int)*(s_patch*s_patch*N_sig_f*N_sim));
    
    double* flux_patchset = malloc(sizeof(double) * (s_patch*s_patch*N_bin*N_sim));
    
    for (int frame = 0; frame < N_frame; frame++) {
        for (int x=0; x< (size_x-(s_patch-1)); x=x+skip) {
            for (int y=0; y< (size_y-(s_patch-1)); y=y+skip) {
                
                memset(FLUX_patchset_noisy3D, 0, (s_patch*s_patch*N_bin_half*N_sim) * sizeof(fftw_complex));
                
                for (int sim_idx=0; sim_idx<N_sim; sim_idx++) {
                    
                    y_sim = y_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    x_sim = x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    frame_sim = frame_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    
                    //printf("y_sim = %.1f\t x_sim = %.1f\t frame_sim %.1f\n", y_sim, x_sim, frame_sim);
                    
                    fft2_ht4D(FLUX_patchset_noisy3D, FLUX_map_set_noisy1D, (int)x_sim, (int)y_sim, (int)frame_sim, sim_idx);
                    
                }
                
                fft_ht4D(FLUX_patchset_noisy4D, FLUX_patchset_noisy3D);
                
                
                sig_and_noisy4D(FLUX_patchset_noisy4D, SIG_patchset_noisy4D, NOISY_patchset_noisy4D);
                
                mean_abs4D(&ME, &MSE, &S_power, NOISY_patchset_noisy4D, SIG_patchset_noisy4D);
                
                MSE_map_set[(x*size_y*N_frame) + (y*N_frame) + frame] = MSE;
                
                if (N_frame == 0 && N_sim > 0) {
                    noise_sigma = 6;
                } else {
                    noise_sigma = 3;
                }
                noise_scale = 1 + noise_sigma*sqrt(square(tgamma(1)/tgamma(1.5))- 1);
                noise_th = noise_scale*ME;
                
                //printf("ME = %.25f\n\n", ME);
                //printf("MSE = %.25f\n\n", MSE);
                if (IBF) {
                    N_power = MSE;
                    
                    if (N_power/S_power >= th_IBF3D) {
                        IBFmask4D_set[(x*size_y*N_frame) + (y*N_frame) + frame] = 1.0;
                        
                        extract_i_patch4D(x, y, frame, i_patch, i_map_set);
                        
                        //fftn_i_patch(I_patch, i_patch);
                        
                        filtering4D(SIG_patch_ht4D, SIG_patchset_noisy4D, I_patch);
                    } else {
                        shrink_hard_sig4D(noise_th, coeff_hard, SIG_patch_ht4D, SIG_patchset_noisy4D);
                    }
                } else {
                    shrink_hard_sig4D(noise_th, coeff_hard, SIG_patch_ht4D, SIG_patchset_noisy4D);
                }
                
                iffn_flux_patchset(flux_patchset, SIG_patch_ht4D);
                
                weight = 1/(MSE + EPS);
                
                for (int sim_idx=0; sim_idx<N_sim; sim_idx++) {
                    y_sim = y_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    x_sim = x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    frame_sim = frame_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (sim_idx*N_frame) + frame];
                    
                    weight_sum_map_patch4D(weight, sim_idx, frame, x, y, weight_sum_map_set, flux_map_set_ht4D, flux_patchset);
                }
            }
        }
    }
    filtered_flux4D(flux_map_set_ht4D, weight_sum_map_set);
    
}

void filtered_flux4D(double *flux_map_set_ht4D, double *weight_sum_map_set) {
    
    for (int frame=0; frame<N_frame; frame++) {
        for (int z=0; z<N_bin; z++) {
            //printf("\n\n BIN NUMER %d \n\n", z);
            for (int x=0; x<size_x; x++) {
                for (int y=0; y<size_y; y++) {
                    flux_map_set_ht4D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame] = fmax(( flux_map_set_ht4D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame] / weight_sum_map_set[(x*size_y*N_frame) + (y*N_frame) + frame]), 0);

                    //printf("flux_map_set_ht4D[%d][%d][%d][%d] = %.10f\n", x, y, z, frame, flux_map_set_ht4D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame]);
                }
            }
        }
    }
}

void weight_sum_map_patch4D(double weight, int sim_idx, int frame, int start_x, int start_y, double *weight_sum_map, double *flux_map_set_ht4D, double *flux_patchset) {
    
    int count_x = 0;
    int count_y = 0;
    
    for (int x=start_x; x<(start_x+s_patch); x++) {
        for (int y=start_y; y<(start_y+s_patch); y++) {
            weight_sum_map[(x*size_y*N_frame) + (y*N_frame) + frame] = weight_sum_map[(x*size_y*N_frame) + (y*N_frame) + frame] + weight;
           
            for (int z=0; z<N_bin; z++) {
                flux_map_set_ht4D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame]  = flux_map_set_ht4D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame] + weight* flux_patchset[(count_x*s_patch*N_bin*N_sim) + (count_y*N_bin*N_sim) + (z*N_sim) + sim_idx];
            }
            count_y++;
        }
        count_x++;
    }
}

void shrink_hard_sig4D(double noise_th, int *coeff_hard, fftw_complex *SIG_patch_ht4D, fftw_complex *SIG_patchset_noisy4D) {
    double complex_mag;
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int z=0; z<N_sig_f; z++) {
            //printf("BIN NUMER %d\n\n\n", z);
            for (int x=0; x<s_patch; x++) {
                for (int y=0; y<s_patch; y++) {
                
                    complex_mag =  sqrt(square(SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0]) + square(SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1]) );
                    
                    if (complex_mag <= noise_th) {
                        SIG_patch_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = 0.0;
                        SIG_patch_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = 0.0;
                        
                        coeff_hard[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim] = 0;
                    } else {
                        SIG_patch_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0];
                        SIG_patch_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1];
                        
                        coeff_hard[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim] = 1;
                    }
                    
                    //printf("shrink_hard_sig: SIG_patch_ht3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0],   SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
                    
                    //printf("shrink_hard_sig: coeff_hard[%d][%d][%d] = %d\n", x, y, z, coeff_hard[(x*s_patch*N_sig_f) + (y*N_sig_f) + z]);
                }
            }
        }
    }
    
}


void filtering4D(fftw_complex *SIG_patch_ht4D, fftw_complex *SIG_patchset_noisy4D, double *I_patch) {
    for (int sim=0; sim<N_sim; sim++) {
        for (int z=0; z<N_sig_f; z++) {
            //printf("BIN NUMER %d\n\n\n", z);
            for (int x=0; x<s_patch; x++) {
                for (int y=0; y<s_patch; y++) {
                    SIG_patch_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] * I_patch[(x*size_y) + y];
                    
                    SIG_patch_ht4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] * I_patch[(x*size_y) + y];
                    
                    //printf("SIG_patch_ht3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0],   SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
                }
            }
        }
    }
}

void extract_i_patch4D(int start_x, int start_y, int frame, double *i_patch, double *i_map_set) {
    int count_x = 0;
    int count_y = 0;
    double sum = 0.0;
    double tmp;
    //printf("\n\n\n\n\n\n\n\n\n");
    count_x = 0;
    for (int x=start_x; x<(start_x+s_patch); x++) {
        count_y = 0;
        for (int y=start_y; y<(start_y+s_patch); y++) {
            tmp = i_map_set[(x*size_y*N_frame) + (y*N_frame) + frame];
            i_patch[(count_x*s_patch) + count_y] = tmp;
            sum = sum + tmp;
            count_y++;
        }
        count_x++;
    }
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            i_patch[(x*s_patch) + y] /= sum;
            //printf("i_patch[%d][%d] = %.3f\n", x, y, i_patch[(x*s_patch) + y]);
        }
    }
}


void sig_and_noisy4D(fftw_complex* FLUX_patchset_noisy4D, fftw_complex *SIG_patchset_noisy4D, fftw_complex *NOISY_patchset_noisy4D) {
    
    int count_z = 0;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            for (int sim=0; sim<N_sim; sim++) {
                
                count_z = 0;
                for (int z=0; z<N_bin_half; z++) {
                    if (z<N_sig_f) {
                        SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0] = FLUX_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0];
                        SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1] = FLUX_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1];
                        
                        //printf("SIG_patchset_noisy4D[%d][%d][%d][%d] = %.15f + %.15fi\n", x, y, z, sim, SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0], SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1]);
                        
                    } else {
                        NOISY_patchset_noisy4D[(x*s_patch*(N_bin_half-N_sig_f)*N_sim) + (y*N_sig_f*N_sim) + (count_z*N_sim) + sim][0] = FLUX_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0];
                        NOISY_patchset_noisy4D[(x*s_patch*(N_bin_half-N_sig_f)*N_sim) + (y*N_sig_f*N_sim) + (count_z*N_sim) + sim][1] = FLUX_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1];
                        
                        //printf("NOISY_patchset_noisy4D[%d][%d][%d][%d] = %.15f + %.15fi\n", x, y, count_z, sim,  NOISY_patchset_noisy4D[(x*s_patch*(N_bin_half-N_sig_f)*N_sim) + (y*N_sig_f*N_sim) + (count_z*N_sim) + sim][0] , NOISY_patchset_noisy4D[(x*s_patch*(N_bin_half-N_sig_f)*N_sim) + (y*N_sig_f*N_sim) + (count_z*N_sim) + sim][1]);
                        
                        count_z++;
                    }
                }
            }
        }
    }
    
}

void mean_abs4D(double *ME, double *MSE, double *S_power, fftw_complex *NOISY_patchset_noisy4D, fftw_complex *SIG_patchset_noisy4D) {
    double sum_ME = 0.0;
    double sum_MSE = 0.0;
    double num_NOISY = s_patch*s_patch*(N_bin_half-N_sig_f)*N_sim;
    double complex_mag_NOISY = 0.0;
    double complex_mag_SIG = 0.0;
    double sum_S_power = 0.0;
    double num_SIG = s_patch*s_patch*N_sig_f*N_sim;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            for (int sim=0; sim<N_sim; sim++) {
                for (int z=0; z<(N_bin_half-N_sig_f); z++) {
                    complex_mag_NOISY = sqrt(square(NOISY_patchset_noisy4D[(x*s_patch*(N_bin_half-N_sig_f)*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0]) + square(NOISY_patchset_noisy4D[(x*s_patch*(N_bin_half-N_sig_f)*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1]) );
                    sum_ME = sum_ME + complex_mag_NOISY;
                    sum_MSE = sum_MSE + square(complex_mag_NOISY);
                }
                for(int z=0; z<N_sig_f; z++) {
                    complex_mag_SIG = sqrt(square(SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0]) + square(SIG_patchset_noisy4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1]));
                    sum_S_power = sum_S_power + square(complex_mag_SIG);
                }
                
            }
        }
    }
    *ME = sum_ME / num_NOISY;
    *MSE = sum_MSE / num_NOISY;
    *S_power = sum_S_power / num_SIG;
}


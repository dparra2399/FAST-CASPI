#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"


/*
 RETURN PARAMS: flux_map_ht3D, IBFmask3D, MSE_map
 
 */
void ht3D(fftw_complex *FLUX_map_noisy1D, double *i_map_1D, double *MSE_map, double *flux_map_ht3D, double *IBFmask3D) {
    
    double* weight_sum_map = (double*)calloc((size_x*size_y), sizeof(double)); //map of weight sum
    
    int* coeff_hard = malloc(sizeof(int)*(s_patch*s_patch*N_sig_f));
    
    fftw_complex* FLUX_patch_noisy1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_bin_half));
    
    fftw_complex* FLUX_patch_noisy3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_bin_half));
    
    fftw_complex* SIG_patch_noisy3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f));
    
    fftw_complex* SIG_patch_ht3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f));
    
    
    double ME;
    double MSE;
    double noise_scale;
    double noise_th;
    double S_power;
    double N_power;
    double weight;
    
    //highSNR:(y,x)=(58,70), medimumSNR:(138,99), lowSNR:(28,103)
    for (int x=0; x< (size_x-(s_patch-1)); x=x+skip) {
        for (int y=0; y< (size_y-(s_patch-1)); y=y+skip) {
            
            fft2_patch(FLUX_map_noisy1D, FLUX_patch_noisy3D, SIG_patch_noisy3D, x, y, &ME, &MSE, &S_power);
            //printf("ME = %.25f\n\n", ME);
            //printf("MSE = %.25f\n\n", MSE);
            MSE_map[(x*size_y) + y] = MSE;

            noise_scale = 1 + 3*sqrt(square(tgamma(1)/tgamma(1.5))- 1);
            noise_th = noise_scale*ME;

            //Itensity based-filtering or hard thresholding
            if (IBF) {
                N_power = MSE; //pseudo signal power before hard-thresholding

               if (N_power/S_power >= th_IBF3D) {
                   IBFmask3D[(x*size_y) + y] = 1.0;

                   filtering_fftn(i_map_1D, SIG_patch_noisy3D, SIG_patch_ht3D, x, y);

               } else {
                   shrink_hard_sig(noise_th, coeff_hard, SIG_patch_ht3D, SIG_patch_noisy3D);
               }
            } else {
                shrink_hard_sig(noise_th, coeff_hard, SIG_patch_ht3D, SIG_patch_noisy3D);
            }
            
            weight = 1/(MSE + EPS);
            
            iffn_flux_patch(flux_map_ht3D, weight_sum_map, SIG_patch_ht3D, weight, x, y);

        }
    }
    
    fftw_export_wisdom_to_filename("wisdom.txt");
    
    filtered_flux(flux_map_ht3D, weight_sum_map);
    
    free(weight_sum_map);
    free(SIG_patch_ht3D);
    free(SIG_patch_noisy3D);
    free(FLUX_patch_noisy3D);
    free(FLUX_patch_noisy1D);
    free(coeff_hard);
    
    
}

void filtered_flux(double *flux_map_3D, double *weight_sum_map) {
    
    for (int z=0; z<N_bin; z++) {
        //printf("\n\n BIN NUMER %d \n\n", z);
        for (int x=0; x<size_x; x++) {
            for (int y=0; y<size_y; y++) {
                flux_map_3D[(x*size_y*N_bin) + (y*N_bin) + z] = fmax((flux_map_3D[(x*size_y*N_bin) + (y*N_bin) + z] / weight_sum_map[(x*size_y) + y]), 0);
                
                //printf("flux_map_3D[%d][%d][%d] = %.10f\n", x, y, z, flux_map_3D[(x*size_y*N_bin) + (y*N_bin) + z]);
            }
        }
    }
}

void shrink_hard_sig(double noise_th, int *coeff_hard, fftw_complex *SIG_patch_ht3D, fftw_complex *SIG_patch_noisy3D) {
    double complex_mag;
    
    for (int z=0; z<N_sig_f; z++) {
        //printf("BIN NUMER %d\n\n\n", z);
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
            
                complex_mag =  sqrt(square(SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0]) + square(SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]) );
                
                if (complex_mag <= noise_th) {
                    SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = 0.0;
                    SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = 0.0;
                    
                    coeff_hard[(x*s_patch*N_sig_f) + (y*N_sig_f) + z] = 0;
                } else {
                    SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0];
                    SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1];
                    
                    coeff_hard[(x*s_patch*N_sig_f) + (y*N_sig_f) + z] = 1;
                }
                
                //printf("shrink_hard_sig: SIG_patch_ht3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0],   SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
                
                //printf("shrink_hard_sig: coeff_hard[%d][%d][%d] = %d\n", x, y, z, coeff_hard[(x*s_patch*N_sig_f) + (y*N_sig_f) + z]);
            }
        }
    }
    
}

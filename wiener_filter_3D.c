#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"

void wf3D(double *flux_map_wiener3D, double *flux_map_ht3D, fftw_complex *FLUX_map_noisy1D, double *MSE_map) {
    
    double *weight_sum_map = (double*)calloc((size_y*size_x), sizeof(double));
    
    fftw_complex* FLUX_patch_ht1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_bin_half));
    
    fftw_complex* FLUX_patch_ht3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_bin_half));
    
    fftw_complex* FLUX_patch_noisy1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_bin_half));
    
    fftw_complex* FLUX_patch_noisy3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_bin_half));
    
    fftw_complex* SIG_patch_noisy3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f));
    
    fftw_complex* SIG_patch_ht3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f));
    
    fftw_complex* SIG_patch_wiener3D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(s_patch*s_patch*N_sig_f));
    
    fftw_complex *FLUX_map_ht1D = fft_wf3D(flux_map_ht3D);
    
    double *S_power = malloc(sizeof(double) * (s_patch*s_patch*N_sig_f));
    
    double N_power;
    double S_power_single;
    double MSE;
    double ME;
    
    double weight;
    
    for (int x=0; x< (size_x-(s_patch-1)); x=x+skip) {
        for (int y=0; y< (size_y-(s_patch-1)); y=y+skip) {
            fft2_patch(FLUX_map_ht1D, FLUX_patch_ht3D, SIG_patch_ht3D, x, y, &ME, &MSE, &S_power_single);
            
            fft2_patch(FLUX_map_noisy1D, FLUX_patch_noisy3D, SIG_patch_noisy3D, x, y, &ME, &MSE, &S_power_single);
            
            N_power = MSE_map[(x*size_y) + y];
            
            wiener_srinkage(SIG_patch_wiener3D, SIG_patch_noisy3D, N_power, S_power);
            
            weight = 1/(N_power + EPS);
            
            iffn_flux_patch(flux_map_wiener3D, weight_sum_map, SIG_patch_wiener3D, weight, x, y);
            
        }
    }
    filtered_flux(flux_map_wiener3D, weight_sum_map);
}

void wiener_srinkage(fftw_complex *SIG_patch_wiener3D, fftw_complex *SIG_patch_noisy3D, double N_power, double *S_power) {
    double wiener_coeff = 0.0;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            for (int z=0; z<N_sig_f; z++) {
                wiener_coeff = 1/(1 + N_power/S_power[(x*s_patch*N_sig_f) + (y*N_sig_f) + z]);
                
                SIG_patch_wiener3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = wiener_coeff * SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0];
                
                SIG_patch_wiener3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = wiener_coeff * SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1];
                
                //printf("SIG_patch_wiener3D[%d][%d][%d] = %.15f + %.15f\n", x, y, z,  SIG_patch_wiener3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0],  SIG_patch_wiener3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
            }
        }
    }
}
void make_s_power(double *S_power, fftw_complex *SIG_patch_ht3D) {
    double complex_mag = 0.0;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            for (int z=0; z<N_sig_f; z++) {
                complex_mag = pow(sqrt(square(SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0]) + square(SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1])), 2);
                
                S_power[(x*s_patch*N_sig_f) + (y*N_sig_f) + z] = complex_mag;
                if (complex_mag == 0.0) {
                    S_power[(x*s_patch*N_sig_f) + (y*N_sig_f) + z] = EPS;
                }
                
                //printf("S_power[%d][%d][%d] = .%.30f\n", x, y, z, S_power[(x*s_patch*N_sig_f) + (y*N_sig_f) + z]);
            }
        }
    }
    
}

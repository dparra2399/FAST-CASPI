#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"


void fft2_patch(fftw_complex *FLUX_map_1D, fftw_complex *FLUX_patch_3D, fftw_complex *SIG_patch_noisy3D, int start_x, int start_y, double *ME, double *MSE, double *S_power) {

    int count_x = 0;
    int count_y = 0;
    int count_z = 0;
    
    double sum_ME = 0.0;
    double sum_MSE = 0.0;
    double complex_mag_NOISY = 0.0;
    double complex_mag_SIG = 0.0;
    double sum_S_power = 0.0;
    double num_NOISY = s_patch*s_patch*(N_bin_half-N_sig_f);
    double num_SIG = s_patch*s_patch*N_sig_f;
    
    fftw_init_threads();

    const int dims = s_patch * s_patch * N_bin_half;
    fftw_complex *input = fftw_alloc_complex(dims);
    const int outdims = s_patch * s_patch * N_bin_half;
    fftw_complex *output = fftw_alloc_complex(outdims);


    // setup "plans" for forward and backward transforms
    const int rank = 2;
    const int howmany = N_bin_half;
    const int istride = N_bin_half;
    const int ostride = N_bin_half;
    const int idist = 1;
    const int odist = 1;
    int n[] = {s_patch, s_patch};
    int *inembed = NULL, *onembed = NULL;

    fftw_plan fp = fftw_plan_many_dft(rank, n, howmany, input, inembed, istride, idist, output, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    

    for (int z=0; z<N_bin_half; z++) {
        count_x = 0;
        for (int x=start_x; x<(start_x+s_patch); x++) {
            count_y = 0;
            for (int y=start_y; y<(start_y+s_patch); y++) {
                input[(count_x*s_patch*N_bin_half) + (count_y*N_bin_half) + z][0] = FLUX_map_1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0];
                input[(count_x*s_patch*N_bin_half) + (count_y*N_bin_half) + z][1] = FLUX_map_1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1];
                
                count_y++;
            }
            count_x++;
        }
    }
    
    fftw_execute(fp);
    
    
    /*PUT INTO NEW FUNCTION*/
    for (int y=0; y<s_patch; y++) {
        for (int x=0; x<s_patch; x++) {
            for (int z=0; z<N_bin_half; z++) {
                if (z<N_sig_f) {
                    SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = output[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
                    SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = output[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
                    
                    complex_mag_SIG = sqrt(square(output[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0]) + square(output[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1]));
                    sum_S_power = sum_S_power + square(complex_mag_SIG);
                } else {
                    
                    complex_mag_NOISY = sqrt(square(output[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0]) + square(output[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1]) );
                    sum_ME = sum_ME + complex_mag_NOISY;
                    sum_MSE = sum_MSE + square(complex_mag_NOISY);
                }
            }
        }
    }
    
    *ME = sum_ME / num_NOISY;
    *MSE = sum_MSE / num_NOISY;
    *S_power = sum_S_power / num_SIG;
    
    
    fftw_free(input);
    fftw_free(output);
    fftw_cleanup_threads();
    
    //fftw_destroy_plan(fp);
}


void filtering_fftn(double *i_map_1D, fftw_complex *SIG_patch_noisy3D, fftw_complex *SIG_patch_ht3D, int start_x, int start_y) {
    
    double sum = 0.0;
    int count_y = 0;
    int count_x = 0;
    double mag;
    
    /*SUM ALL ELEMENTS*/
    for (int x=start_x; x<(start_x+s_patch); x++) {
        for (int y=start_y; y<(start_y+s_patch); y++) {
            sum = sum + i_map_1D[(x*size_y) + y];
        }
    }
    /*SUM ALL ELEMENTS*/
    
    const int dims = s_patch * s_patch;
    fftw_complex *input = fftw_alloc_complex(dims);
    const int outdims = s_patch * s_patch;
    fftw_complex *output = fftw_alloc_complex(outdims);


    fftw_init_threads();
    
    
    // setup "plans" for forward and backward transforms
    const int rank = 2;
    const int howmany = 1;
    int idist = 0;
    int odist = 0; /* unused because howmany = 1 */
    int istride = 1;
    int ostride = 1; /* array is contiguous in memory */
    int n[] = {s_patch, s_patch};
    int *inembed = NULL, *onembed = NULL;
    

    fftw_plan fp = fftw_plan_many_dft(rank, n, howmany, input, inembed, istride, idist, output, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    
    
    for (int x=start_x; x<(start_x+s_patch); x++) {
        count_y = 0;
        for (int y=start_y; y<(start_y+s_patch); y++) {
            input[(count_x*s_patch) + count_y][0] = (double)i_map_1D[(x*size_y) + y] / (double)sum;
            input[(count_x*s_patch) + count_y][1] = 0.0;
            count_y++;
        }
        count_x++;
    }
    
    fftw_execute(fp);
    
    
    for (int y=0; y<s_patch; y++) {
        for (int x=0; x<s_patch; x++) {
            mag = sqrt( square(output[(x*s_patch) + y][1]) + square(output[(x*s_patch) + y][0]) );
            for (int z=0; z<N_sig_f; z++) {
                SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] * mag;
                
                SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] * mag;
            }
        }
    }
    
    fftw_free(input);
    fftw_free(output);
    fftw_cleanup_threads();
    
    fftw_destroy_plan(fp);
}

void iffn_flux_patch(double *flux_map_ht3D, double *weight_sum_map, fftw_complex *SIG_patch_ht3D, double weight, int start_x, int start_y) {

    int count_x = start_x;
    int count_y = start_y;
    
    fftw_init_threads();
    
    const int dims = s_patch * s_patch * N_bin;
    fftw_complex *input = fftw_alloc_complex(dims);
    const int outdims = s_patch * s_patch * N_bin;
    double *output = fftw_alloc_real(outdims);
    
    
    // setup "plans" for forward and backward transforms
    const int rank = 3;
    const int howmany = 1;
    int idist = 0;
    int odist = 0; /* unused because howmany = 1 */
    int istride = 1;
    int ostride = 1; /* array is contiguous in memory */
    int n[] = {s_patch, s_patch, N_bin};
    int *inembed = NULL, *onembed = NULL;

    /*#ifdef WISDOM
        fftw_import_wisdom_from_filename("my_wisdom_file");
    
        fftw_plan bp = fftw_plan_many_dft_c2r(rank, n, howmany, input, inembed, istride, idist, output, onembed, ostride, odist, FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    #else*/
    
    fftw_plan bp = fftw_plan_many_dft_c2r(rank, n, howmany, input, inembed, istride, idist, output, onembed, ostride, odist, FFTW_ESTIMATE);
    
    //#endif
    
    for (int x=0; x<s_patch;x++) {
        for (int y=0; y<s_patch; y++) {
            for (int z=0; z<N_bin; z++) {
                if (z<N_sig_f) {
                    input[(x*s_patch*N_bin) + (y*N_bin) + z][0] = SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0];
                    input[(x*s_patch*N_bin) + (y*N_bin) + z][1] = SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1];
                } else {
                    input[(x*s_patch*N_bin) + (y*N_bin) + z][0] = 0.0;
                    input[(x*s_patch*N_bin) + (y*N_bin) + z][1] = 0.0;
                }
            }
        }
    }
    
    fftw_execute(bp);
    
    for (int x=0; x<s_patch; x++) {
        count_y = start_y;
        for (int y=0; y<s_patch; y++) {
            for (int z=0; z<N_bin; z++) {
                flux_map_ht3D[(count_x*size_y*N_bin) + (count_y*N_bin) + z] = flux_map_ht3D[(count_x*size_y*N_bin) + (count_y*N_bin) + z] + weight * (output[(x*s_patch*N_bin) + (y*N_bin) + z]/(s_patch*s_patch));
            }
            
            weight_sum_map[(count_x*size_y) + count_y] = weight_sum_map[(count_x*size_y) + count_y] + weight;
            count_y++;
        }
        count_x++;
    }

    /*#ifndef WISDOM
        fftw_export_wisdom_to_filename("my_wisdom_file");
    #endif*/
    
    fftw_free(input);
    fftw_free(output);
    
    fftw_cleanup_threads();
    fftw_destroy_plan(bp);
}



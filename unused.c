#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"
#include <time.h>

void extract_patch(fftw_complex* FLUX_map_1D, fftw_complex* FLUX_patch_1D, int start_x, int start_y) {
    int count_x = 0;
    int count_y = 0;
    //printf("\n\n\n\n\n\n\n\n\n");
    for (int z=0; z<N_bin_half; z++) {
        count_x = 0;
        //printf("\nBin Number %d\n", z);
        for (int x=start_x; x<(start_x+s_patch); x++) {
            count_y = 0;
            for (int y=start_y; y<(start_y+s_patch); y++) {
                FLUX_patch_1D[(count_x*s_patch*N_bin_half) + (count_y*N_bin_half) + z][0] = FLUX_map_1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0];
                FLUX_patch_1D[(count_x*s_patch*N_bin_half) + (count_y*N_bin_half) + z][1] = FLUX_map_1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1];
                
                //printf("FLUX_patch_1D[%d][%d][%d] = %.30f + %.30fi\n", count_x, count_y, z, FLUX_patch_1D[(count_x*s_patch*N_bin_half) + (count_y*N_bin_half) + z][0], FLUX_patch_1D[(count_x*s_patch*N_bin_half) + (count_y*N_bin_half) + z][1]);
                count_y++;
            }
            count_x++;
        }
    }
}

void mean_abs(double *ME, double *ME_sqr, double *S_power, fftw_complex *NOISY_patch_noisy3D, fftw_complex *SIG_patch_noisy3D) {
    double sum_ME = 0.0;
    double sum_ME_sqr = 0.0;
    double num_NOISY = s_patch*s_patch*(N_bin_half-N_sig_f);
    double complex_mag_NOISY = 0.0;
    double complex_mag_SIG = 0.0;
    double sum_S_power = 0.0;
    double num_SIG = s_patch*s_patch*N_sig_f;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            for (int z=0; z<(N_bin_half-N_sig_f); z++) {
                complex_mag_NOISY = sqrt(square(NOISY_patch_noisy3D[(x*s_patch*(N_bin_half-N_sig_f)) + (y*(N_bin_half-N_sig_f)) + z][0]) + square(NOISY_patch_noisy3D[(x*s_patch*(N_bin_half-N_sig_f)) + (y*(N_bin_half-N_sig_f)) + z][1]) );
                sum_ME = sum_ME + complex_mag_NOISY;
                sum_ME_sqr = sum_ME_sqr + square(complex_mag_NOISY);
            }
            for(int z=0; z<N_sig_f; z++) {
                complex_mag_SIG = sqrt(square(SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0]) + square(SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]) );
                sum_S_power = sum_S_power + square(complex_mag_SIG);
            }
            
        }
    }
    *ME = sum_ME / num_NOISY;
    *ME_sqr = sum_ME_sqr / num_NOISY;
    *S_power = sum_S_power / num_SIG;
}

void sig_and_noisy(fftw_complex* FLUX_patch_noisy3D, fftw_complex *SIG_patch_noisy3D, fftw_complex *NOISY_patch_noisy3D) {
    
    int count_z = 0;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            count_z = 0;
            for (int z=0; z<N_bin_half; z++) {
                if (z<N_sig_f) {
                    SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = FLUX_patch_noisy3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
                    SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = FLUX_patch_noisy3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
                    
                } else {
                    NOISY_patch_noisy3D[(x*s_patch*(N_bin_half-N_sig_f)) + (y*(N_bin_half-N_sig_f)) + count_z][0] = FLUX_patch_noisy3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
                    NOISY_patch_noisy3D[(x*s_patch*(N_bin_half-N_sig_f)) + (y*(N_bin_half-N_sig_f)) + count_z][1] = FLUX_patch_noisy3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
                    
                    count_z++;
                }
            }
        }
    }
    
   /*for (int z=0; z<N_sig_f; z++) {
        printf("\n\n\nBIN NUMBER %d\n", z);
        for (int x=0; x<(s_patch); x++) {
            for (int y=0; y<(s_patch); y++) {
                printf("SIG_patch_noisy3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0], SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
            }
        }
    }
    printf("\n\n\n\n\n\n\n\n\n\n\n\n");
    for (int z=0; z<(N_bin_half-N_sig_f); z++) {
        printf("\n\n\nBIN NUMBER %d\n", z);
        for (int x=0; x<(s_patch); x++) {
            for (int y=0; y<(s_patch); y++) {
                printf("NOISY_patch_noisy3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, NOISY_patch_noisy3D[(x*s_patch*(N_bin_half-N_sig_f)) + (y*(N_bin_half-N_sig_f)) + z][0], NOISY_patch_noisy3D[(x*s_patch*(N_bin_half-N_sig_f)) + (y*(N_bin_half-N_sig_f)) + z][1]);
            }
        }
    }*/
    
}
void filtering(fftw_complex *SIG_patch_ht3D, fftw_complex *SIG_patch_noisy3D, double *I_patch) {
    
    for (int z=0; z<N_sig_f; z++) {
        //printf("BIN NUMER %d\n\n\n", z);
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] * I_patch[(x*size_y) + y];
                
                SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] * I_patch[(x*size_y) + y];
                
                //printf("SIG_patch_ht3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0],   SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
            }
        }
    }
}



void extract_i_patch(int start_x, int start_y, double *i_patch, double *i_map_1D) {
    int count_x = 0;
    int count_y = 0;
    double sum = 0.0;
    //printf("\n\n\n\n\n\n\n\n\n");
    count_x = 0;
    for (int x=start_x; x<(start_x+s_patch); x++) {
        count_y = 0;
        for (int y=start_y; y<(start_y+s_patch); y++) {
            i_patch[(count_x*s_patch) + count_y] = i_map_1D[(x*size_y) + y];
            sum = sum + i_patch[(count_x*s_patch) + count_y];
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

void flux_map_3D_patch(double weight, int start_x, int start_y, double *flux_map_3D, double *flux_patch) {
    
    int count_x = 0;
    int count_y = 0;
    
    for (int z=0; z<N_bin; z++) {
        count_x = 0;
        for (int x=start_x; x<(start_x+s_patch); x++) {
            count_y = 0;
            for (int y=start_y; y<(start_y+s_patch); y++) {
                flux_map_3D[(x*size_y*N_bin) + (y*N_bin) + z] = flux_map_3D[(x*size_y*N_bin) + (y*N_bin) + z] + weight* flux_patch[(count_x*s_patch*N_bin) + (count_y*N_bin) + z];
                
                //printf("flux_map_ht3d[%d][%d][%d] = %.10f\n", x, y, z, flux_map_3D[(x*size_y*N_bin) + (y*N_bin) + z]);
                
                //printf("flux_patch[%d][%d][%d] = %.15f\n", count_x, count_y, z, flux_patch[(count_x*s_patch*N_bin) + (count_y*N_bin) + z]);
                count_y++;
            }
            count_x++;
        }
    }
}

void weight_sum_map_patch(double weight, int start_x, int start_y, double *weight_sum_map) {
    for (int x=start_x; x<(start_x+s_patch); x++) {
        for (int y=start_y; y<(start_y+s_patch); y++) {
            weight_sum_map[(x*size_y) + y] = weight_sum_map[(x*size_y) + y] + weight;
            //printf("weight_sum_map[%d][%d] = %.3f\n", x, y, weight_sum_map[(x*size_y) + y]);
        }
    }
}



/*
 ALL THE FFTW STUFF THAT WAS USED
 */
/*
for (int z=0; z<N_bin_half; z++) {
    count_x = 0;
    for (int x=start_x; x<(start_x+s_patch); x++) {
        count_y = 0;
        for (int y=start_y; y<(start_y+s_patch); y++) {
            in[count_y][0] = FLUX_map_1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0];
            in[count_y][1] = FLUX_map_1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1];
            count_y++;
           
        }
        
        p = fftw_plan_dft_1d(s_patch, in, out, FFTW_FORWARD, FFTW_MEASURE);
        fftw_execute(p);
        
        for (int y=0; y<s_patch; y++) {
            FLUX_patch_3D[(count_x*s_patch*N_bin_half) + (y*N_bin_half) + z][0] = out[y][0];
            FLUX_patch_3D[(count_x*s_patch*N_bin_half) + (y*N_bin_half) + z][1] = out[y][1];
        }
        
        fftw_destroy_plan(p);
        count_x++;
    }
}

count_z = 0;
for (int z=0; z<N_bin_half; z++) {
    for (int y=0; y<s_patch; y++) {
        for (int x=0; x<s_patch; x++) {
            in[x][0] = FLUX_patch_3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
            in[x][1] = FLUX_patch_3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
        }
        
        p = fftw_plan_dft_1d(s_patch, in, out, FFTW_FORWARD, FFTW_MEASURE);
        fftw_execute(p);
        
        for (int x=0; x<s_patch; x++) {
            FLUX_patch_3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0] = out[x][0];
            FLUX_patch_3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1] = out[x][1];
            
            if (z<N_sig_f) {
                SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = out[x][0];
                SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = out[x][1];
                
                complex_mag_SIG = sqrt(square(out[x][0]) + square(out[x][1]));
                sum_S_power = sum_S_power + square(complex_mag_SIG);
            } else {
                
                complex_mag_NOISY = sqrt(square(out[x][0]) + square(out[x][1]) );
                sum_ME = sum_ME + complex_mag_NOISY;
                sum_MSE = sum_MSE + square(complex_mag_NOISY);
                
                count_z++;
            }
        }
        fftw_destroy_plan(p);
    }
}*/

/*
for (int x=start_x; x<(start_x+s_patch); x++) {
    count_y = 0;
    for (int y=start_y; y<(start_y+s_patch); y++) {
        in[count_y][0] = (double)i_map_1D[(x*size_y) + y] / (double)sum;
        in[count_y][1] = 0.0;
        //printf("in[%d] = %.3f + %.3fi\n", count_y, in[count_y][0], in[count_y][1]);
        count_y++;
    }
    
    p = fftw_plan_dft_1d(s_patch, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(p);
    
    for (int y=0; y<s_patch; y++) {
        I_patch_tmp[(count_x*s_patch) + y][0] = out[y][0];
        I_patch_tmp[(count_x*s_patch) + y][1] = out[y][1];
        //printf("out[%d] = %.15f + %.15fi\n", y, out[y][0], out[y][1]);
    }
    count_x++;
    fftw_destroy_plan(p);
}

for (int y=0; y<s_patch; y++) {
    for (int x=0; x<s_patch; x++) {
        in[x][0] = I_patch_tmp[(x*s_patch) + y][0];
        in[x][1] = I_patch_tmp[(x*s_patch) + y][1];
        //printf("in[%d] = %.15f + %.15fi\n", x, in[x][0], in[x][1]);
    }
    
    p = fftw_plan_dft_1d(s_patch, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(p);
    
    for (int x=0; x<s_patch; x++) {
        //printf("out[%d] = %.15f + %.15fi\n", x, out[x][0], out[x][1]);
        mag = sqrt( square(out[x][1]) + square(out[x][0]) );
        for (int z=0; z<N_sig_f; z++) {
            SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] * mag;
            
            SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] * mag;
        }
    }
    fftw_destroy_plan(p);
}*/




/*
for (int z=0; z<N_bin; z++) {
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            if (z<N_sig_f) {
                in_patch[y][0] = SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0];p
                in_patch[y][1] = SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1];
            } else {
                in_patch[y][0] = 0.0;
                in_patch[y][1] = 0.0;
            }
        }
        p = fftw_plan_dft_1d(s_patch, in_patch, out_patch, FFTW_BACKWARD, FFTW_MEASURE);
        fftw_execute(p);

        for (int y=0; y<s_patch; y++) {
            out_patch[y][0] *= 1./s_patch;
            out_patch[y][1] *= 1./s_patch;
            temp[(x*s_patch*N_bin) + (y*N_bin) + z][0] = out_patch[y][0];
            temp[(x*s_patch*N_bin) + (y*N_bin) + z][1] = out_patch[y][1];
        }
        fftw_destroy_plan(p);
    }
}
for (int z=0; z<N_bin; z++) {
    for (int y=0; y<s_patch; y++) {
        for (int x=0; x<s_patch; x++) {
            in_patch[x][0] = temp[(x*s_patch*N_bin) + (y*N_bin) + z][0];
            in_patch[x][1] = temp[(x*s_patch*N_bin) + (y*N_bin) + z][1];
        }
        p = fftw_plan_dft_1d(s_patch, in_patch, out_patch, FFTW_BACKWARD, FFTW_MEASURE);
        fftw_execute(p);
        for (int x=0; x<s_patch; x++) {
            out_patch[x][0] *= 1./s_patch;
            out_patch[x][1] *= 1./s_patch;
            temp[(x*s_patch*N_bin) + (y*N_bin) + z][0] = out_patch[x][0];
            temp[(x*s_patch*N_bin) + (y*N_bin) + z][1] = out_patch[x][1];
        }
        fftw_destroy_plan(p);
    }
}

for (int x=0; x<s_patch; x++) {
    count_y = start_y;
    for (int y=0; y<s_patch; y++) {
        for (int z=0; z<N_bin; z++) {
            in_bin[z][0] = temp[(x*s_patch*N_bin) + (y*N_bin) + z][0];
            in_bin[z][1] = temp[(x*s_patch*N_bin) + (y*N_bin) + z][1];
        }
        
        p = fftw_plan_dft_c2r_1d(N_bin, in_bin, out_bin, FFTW_MEASURE);
        fftw_execute(p);
        
        for (int z=0; z<N_bin; z++) {
            out_bin[z] *= 1./N_bin;
            
            flux_map_ht3D[(count_x*size_y*N_bin) + (count_y*N_bin) + z] = flux_map_ht3D[(count_x*size_y*N_bin) + (count_y*N_bin) + z] + weight * out_bin[z];
        }
        
        fftw_destroy_plan(p);
        
        weight_sum_map[(count_x*size_y) + count_y] = weight_sum_map[(count_x*size_y) + count_y] + weight;
        count_y++;
    }
    count_x++;
}*/

void sig_patches(fftw_complex *FLUX_patch_noisy3D, fftw_complex *SIG_patch_noisy3D, fftw_complex *FLUX_patch_ht3D, fftw_complex *SIG_patch_ht3D) {
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            for (int z=0; z<N_sig_f; z++) {
                SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = FLUX_patch_noisy3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
                SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = FLUX_patch_noisy3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
                
                SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0] = FLUX_patch_ht3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
                SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1] = FLUX_patch_ht3D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
            }
        }
    }
    
    /*for (int z=0; z<N_sig_f; z++) {
         printf("\n\n\nBIN NUMBER %d\n", z);
         for (int x=0; x<(s_patch); x++) {
             for (int y=0; y<(s_patch); y++) {
                 printf("SIG_patch_noisy3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0], SIG_patch_noisy3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
             }
         }
     }
     printf("\n\n\n\n\n\n\n\n\n\n\n\n");
     for (int z=0; z<N_sig_f; z++) {
         printf("\n\n\nBIN NUMBER %d\n", z);
         for (int x=0; x<(s_patch); x++) {
             for (int y=0; y<(s_patch); y++) {
                 printf("SIG_patch_ht3D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][0], SIG_patch_ht3D[(x*s_patch*N_sig_f) + (y*N_sig_f) + z][1]);
             }
         }
     }*/
}

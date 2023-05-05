#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"

void fft_ht4D(fftw_complex *FLUX_patchset_noisy4D, fftw_complex *FLUX_patchset_noisy3D) {
    fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_sim);
    fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_sim);
    
    fftw_plan p;
    
    fftw_init_threads();
    
    for (int z=0; z<N_bin_half; z++) {
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                for (int sim=0; sim<N_sim; sim++) {
                    in[sim][0] = FLUX_patchset_noisy3D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][0];
                    in[sim][1] = FLUX_patchset_noisy3D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][1];
                }
                
                p = fftw_plan_dft_1d(N_sim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(p);
                
                for (int sim=0; sim<N_sim; sim++) {
                    FLUX_patchset_noisy4D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][0] = out[sim][0];
                    FLUX_patchset_noisy4D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim][1] = out[sim][1];
                }
            }
        }
    }
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup_threads();
}

void fft2_ht4D(fftw_complex *FLUX_patchset_noisy3D, fftw_complex *FLUX_map_set_noisy1D, int x_sim, int y_sim, int frame_sim, int sim_idx) {
    
    int count_x = 0;
    int count_y = 0;
    
    fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s_patch);
    fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s_patch);
    fftw_complex *FLUX_patch_noisy1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s_patch*s_patch*N_bin_half);
    
    fftw_plan p;
    
    fftw_init_threads();

    //printf("---- fft2_ht4d ---\n");
    for (int z=0; z<N_bin_half; z++) {
        //printf("\n");
        count_x = 0;
        for (int x=x_sim; x<x_sim+s_patch; x++) {
            count_y = 0;
            for (int y=y_sim; y<y_sim+s_patch; y++) {
                
                in[count_y][0] = FLUX_map_set_noisy1D[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame_sim][0];
                in[count_y][1] = FLUX_map_set_noisy1D[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame_sim][1];
                
                //printf("%.15f + %.15fi\t", in[count_y][0], in[count_y][1]);
                
                count_y++;
            }
            p = fftw_plan_dft_1d(s_patch, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            
            for (int y=0; y<s_patch; y++) {
                FLUX_patch_noisy1D[(count_x*s_patch*N_bin_half) + (y*N_bin_half) + z][0] = out[y][0];
                FLUX_patch_noisy1D[(count_x*s_patch*N_bin_half) + (y*N_bin_half) + z][1] = out[y][1];
            }
            fftw_destroy_plan(p);
            count_x++;
        }
    }
    
    for (int z=0; z<N_bin_half; z++) {
        for (int y=0; y<s_patch; y++) {
            for (int x=0; x<s_patch; x++) {
                in[x][0] = FLUX_patch_noisy1D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][0];
                in[x][1] = FLUX_patch_noisy1D[(x*s_patch*N_bin_half) + (y*N_bin_half) + z][1];
            }
            
            p = fftw_plan_dft_1d(s_patch, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            
            for (int x=0; x<s_patch; x++) {
                FLUX_patchset_noisy3D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim_idx][0] = out[x][0];
                FLUX_patchset_noisy3D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim_idx][1] = out[x][1];
                
                  //printf("%.10f + %.10fi \t", FLUX_patchset_noisy3D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim_idx][0], FLUX_patchset_noisy3D[(x*s_patch*N_bin_half*N_sim) + (y*N_bin_half*N_sim) + (z*N_sim) + sim_idx][1]);
            }
            //printf("\n");
            fftw_destroy_plan(p);
        }
    }
    
    fftw_free(in);
    fftw_free(out);
    fftw_free(FLUX_patch_noisy1D);
    fftw_cleanup_threads();
}


fftw_complex* fft_wf3D(double *flux_map_ht3D) {
    fftw_complex *out;
    fftw_complex *in;
    
    fftw_plan p;
    fftw_complex *FLUX_map_ht1D;
    
    FLUX_map_ht1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N_bin_half * size_x * size_y));
    
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    
    fftw_init_threads();
    
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            //printf("\n\n\n");
            for (int z=0; z<N_bin; z++) {
                in[z][0] = flux_map_ht3D[(x*size_y*N_bin) + (y*N_bin) + z];
                in[z][1] = 0.0;
                //printf("in[%d] = %.15f\n", z, in[z][0]);
            }
            
            p = fftw_plan_dft_1d(N_bin, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            
            for (int z=0; z<N_bin_half; z++) {
                
                FLUX_map_ht1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0] = out[z][0];
                FLUX_map_ht1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1] = out[z][1];
                //printf("FLUX_map_ht1D[%d][%d][%d]= %.30f + %.30fi\n", x, y, z, FLUX_map_ht1D[(x*size_y*N_bin) + (y*N_bin_half) + z][0], FLUX_map_ht1D[(x*size_y*N_bin) + (y*N_bin_half) + z][1]);
                    
            }
            fftw_destroy_plan(p);
        }
    }
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup_threads();
    return FLUX_map_ht1D;
}

fftw_complex* fft_wf4D(double *flux_map_set_ht4D) {
    fftw_complex *out;
    fftw_complex *in;
    
    fftw_plan p;
    fftw_complex *FLUX_map_set_ht1D;
    
    FLUX_map_set_ht1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N_bin_half * size_x * size_y * N_frame));
    
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    
    fftw_init_threads();
    
    for (int frame=0; frame<N_frame; frame++) {
        for (int x=0; x<size_x; x++) {
            for (int y=0; y<size_y; y++) {
                //printf("\n\n\n");
                for (int z=0; z<N_bin; z++) {
                    in[z][0] = flux_map_set_ht4D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame];
                    in[z][1] = 0.0;
                    //printf("in[%d] = %.15f\n", z, in[z][0]);
                }
                
                p = fftw_plan_dft_1d(N_bin, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(p);
                
                for (int z=0; z<N_bin_half; z++) {
                    
                    FLUX_map_set_ht1D[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][0] = out[z][0];
                    FLUX_map_set_ht1D[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][1] = out[z][1];
                    //printf("FLUX_map_ht1D[%d][%d][%d]= %.30f + %.30fi\n", x, y, z, FLUX_map_ht1D[(x*size_y*N_bin) + (y*N_bin_half) + z][0], FLUX_map_ht1D[(x*size_y*N_bin) + (y*N_bin_half) + z][1]);
                        
                }
                fftw_destroy_plan(p);
            }
        }
    }
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup_threads();
    return FLUX_map_set_ht1D;
}

void iffn_flux_patchset(double *flux_patchset, fftw_complex *SIG_patch_4D) {
    fftw_complex *in_patch;
    fftw_complex *in_bin;
    fftw_complex *in_sim;
    
    fftw_complex *out_patch;
    fftw_complex *out_bin;
    double *out_sim;
    
    fftw_complex *temp;
    fftw_plan p;
    
    
    in_patch = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s_patch);
    out_patch = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s_patch);
    
    in_bin = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    out_bin = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    
    in_sim = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_sim);
    out_sim = malloc(sizeof(double) * N_sim);
    
    temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s_patch * s_patch * N_bin * N_sim));
    
    fftw_init_threads();
    
    //printf("Checkpoint 0\n");
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int z=0; z<N_bin; z++) {
            for (int x=0; x<s_patch; x++) {
                for (int y=0; y<s_patch; y++) {
                    if (z<N_sig_f) {
                        in_patch[y][0] = SIG_patch_4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][0];
                        in_patch[y][1] = SIG_patch_4D[(x*s_patch*N_sig_f*N_sim) + (y*N_sig_f*N_sim) + (z*N_sim) + sim][1];
                    } else {
                        in_patch[y][0] = 0.0;
                        in_patch[y][1] = 0.0;
                    }
                }
                
                p = fftw_plan_dft_1d(s_patch, in_patch, out_patch, FFTW_BACKWARD, FFTW_ESTIMATE);
                fftw_execute(p);

                for (int y=0; y<s_patch; y++) {
                    out_patch[y][0] *= 1./s_patch;
                    out_patch[y][1] *= 1./s_patch;
                    temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][0] = out_patch[y][0];
                    temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][1] = out_patch[y][1];
                }
                fftw_destroy_plan(p);
            }
        }
    }
    
    //printf("Checkpoint 1\n");
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int z=0; z<N_bin; z++) {
            for (int y=0; y<s_patch; y++) {
                for (int x=0; x<s_patch; x++) {
                    in_patch[x][0] = temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][0];
                    in_patch[x][1] = temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][1];
                }
                
                p = fftw_plan_dft_1d(s_patch, in_patch, out_patch, FFTW_BACKWARD, FFTW_ESTIMATE);
                fftw_execute(p);
                
                for (int x=0; x<s_patch; x++) {
                    out_patch[x][0] *= 1./s_patch;
                    out_patch[x][1] *= 1./s_patch;
                    temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][0] = out_patch[x][0];
                    temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][1] = out_patch[x][1];
                }
            fftw_destroy_plan(p);
            }
        }
    }

    //printf("Checkpoint 2\n");
    
    for (int sim=0; sim<N_sim; sim++) {
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                for (int z=0; z<N_bin; z++) {
                    in_bin[z][0] = temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][0];
                    in_bin[z][1] = temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][1];
                    //printf("in_bin[%d] = %.15f + %.15fi\n", z, in_bin[z][0], in_bin[z][1]);
                }
                
                p = fftw_plan_dft_1d(N_bin, in_bin, out_bin, FFTW_ESTIMATE, FFTW_BACKWARD);
                fftw_execute(p);
                
                for (int z=0; z<N_bin; z++) {
                    out_bin[z][0] *= 1./N_bin;
                    out_bin[z][1] *= 1./N_bin;
                    temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][0] = out_bin[z][0];
                    temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][1] = out_bin[z][1];
                }
            fftw_destroy_plan(p);
            }
        }
    }
    
    //printf("Checkpoint 3\n");
    
    for (int z=0; z<N_bin; z++) {
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                for (int sim=0; sim<N_sim; sim++) {
                    in_sim[sim][0] = temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][0];
                    in_sim[sim][1] = temp[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim][1];
                    //printf("in_bin[%d] = %.15f + %.15fi\n", z, in_bin[z][0], in_bin[z][1]);
                }
                
                p = fftw_plan_dft_c2r_1d(N_sim, in_sim, out_sim, FFTW_ESTIMATE);
                fftw_execute(p);
                
                for (int sim=0; sim<N_sim; sim++) {
                    out_sim[sim] *= 1./N_sim;
                    flux_patchset[(x*s_patch*N_bin*N_sim) + (y*N_bin*N_sim) + (z*N_sim) + sim] = out_sim[sim];
                }
            fftw_destroy_plan(p);
            }
        }
    }
    
    fftw_free(in_patch);
    fftw_free(out_patch);
    fftw_free(in_bin);
    fftw_free(out_bin);
    fftw_free(in_sim);
    fftw_free(out_sim);
    fftw_cleanup_threads();
    
    /*for (int z=0; z<N_bin; z++) {
        printf("NEW BIN %d\n", z);
        for (int x=0; x<s_patch; x++) {
            for (int y=0; y<s_patch; y++) {
                printf("flux_patch[%d][%d][%d] = %.15f\n", x, y, z, flux_patch[(x*s_patch*N_bin) + (y*N_bin) + z]);
            }
        
        }
    }*/
}


double* ifft(fftw_complex *FLUX_map_noisy1D) {
    fftw_complex *in;
    double *out;
    fftw_plan p;
    
    double *FLUX_map_1D = malloc(sizeof(double) * (N_bin * size_x * size_y));
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    out = malloc(sizeof(double) * N_bin);
    fftw_init_threads();
    
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            for (int z=0; z<N_bin; z++) {
                if (z<N_sig_f) {
                    in[z][0] = FLUX_map_noisy1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0];
                    in[z][1] = FLUX_map_noisy1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1];
                } else {
                    in[z][0] = 0.0;
                    in[z][1] = 0.0;
                }
                //printf("in[%d] = %.5f + %.5fi\n", z, in[z][0], in[z][1]);
            }
            p = fftw_plan_dft_c2r_1d(N_bin, in, out, FFTW_ESTIMATE);
            fftw_execute(p);
            for (int z=0; z<N_bin; z++) {
                out[z] *= 1./N_bin;
                if (out[z] < 0) {
                    FLUX_map_1D[(x*size_y*N_bin) + (y*N_bin) + z] = 0;
                } else {
                    FLUX_map_1D[(x*size_y*N_bin) + (y*N_bin) + z] = out[z]; // N_bin;
                }
                //printf("out[%d][0] = %.15f + %.15fi\n", z, out[z][0], out[z][1]);
            //printf("FLUX_map_1D[%d][%d][%d] = %.15f\n", x, y, z, FLUX_map_1D[(x*size_y*N_bin) + (y*N_bin) + z]);
            }
            fftw_destroy_plan(p);
        }
    }
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup_threads();
    return FLUX_map_1D;
}


fftw_complex* fft(double *flux_map_set) {
    fftw_complex *out;
    fftw_complex *in;
    fftw_plan p;
    fftw_complex *FLUX_map_set_noisy1D;
    FLUX_map_set_noisy1D = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N_bin * size_x * size_y * N_frame));
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N_bin);
    
    fftw_init_threads();
    
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            for (int frame=0; frame<N_frame; frame++) {
                //printf("\n\n\n");
                for (int z=0; z<N_bin; z++) {
                    in[z][0] = flux_map_set[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame];
                    in[z][1] = 0.0;
                    //printf("in[%d] = %.15f\n", z, in[z]);
                }
                p = fftw_plan_dft_1d(N_bin, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(p);
                for (int z=0; z<N_bin; z++) {
                    FLUX_map_set_noisy1D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame][0] = out[z][0];
                    FLUX_map_set_noisy1D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame][1] = out[z][1];
                    //printf("FLUX_map_set_noisy1D_full[%d][%d][%d][%d] = %.15f + %.15fi\n", x, y, z, frame, FLUX_map_set_noisy1D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame][0], FLUX_map_set_noisy1D[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame][1]);
                    
                }
                fftw_destroy_plan(p);
            }
        }
    }
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup_threads();
    return FLUX_map_set_noisy1D;
}


fftw_complex* fft_halfbin(fftw_complex *FLUX_map_set_noisy1D_full) {
    fftw_complex *FLUX_map_set_noisy1D_halfbin = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N_bin_half * size_x * size_y * N_frame));
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            for (int z=0; z<N_bin_half; z++) {
                for (int frame=0; frame<N_frame; frame++) {
                    FLUX_map_set_noisy1D_halfbin[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][0] = FLUX_map_set_noisy1D_full[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame][0];
                    
                    FLUX_map_set_noisy1D_halfbin[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][1] = FLUX_map_set_noisy1D_full[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame][1];
                    
                    //printf("FLUX_map_set_noisy1D_halfbin[%d][%d][%d][%d] = %.15f + %.15fi\n", x, y, z, frame, FLUX_map_set_noisy1D_halfbin[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][0], FLUX_map_set_noisy1D_halfbin[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][1]);
                }
            }
        }
    }
    return FLUX_map_set_noisy1D_halfbin;
}


fftw_complex* getNthFrame_fftw(fftw_complex *FLUX_map_set_noisy1D, int frame) {
    fftw_complex *FLUX_map_noisy1D = (fftw_complex*)fftw_malloc((size_x*size_y*N_bin_half) * sizeof(fftw_complex));
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            for (int z=0; z<N_bin_half; z++) {
                FLUX_map_noisy1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0] = FLUX_map_set_noisy1D[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][0];
                FLUX_map_noisy1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1] = FLUX_map_set_noisy1D[(x*size_y*N_bin_half*N_frame) + (y*N_bin_half*N_frame) + (z*N_frame) + frame][1];
            //printf("FLUX_map_noisy1D[%d][%d][%d] = %.15f + %.15fi\n", x, y, z, FLUX_map_noisy1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][0], FLUX_map_noisy1D[(x*size_y*N_bin_half) + (y*N_bin_half) + z][1]);
            }
        }
    }
    return FLUX_map_noisy1D;
}

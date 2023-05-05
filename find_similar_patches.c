#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"
#include <time.h>
#include <string.h>
#include <math.h>

int s_intra;
int s_inter;


void sim_patches(double* i_map_set_wiener3D, double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set) {
    
    double *patch = malloc(sizeof(double)*(s_patch*s_patch));
    double *T = malloc(sizeof(double)* (s_patch*s_patch));
    
    int frame_start;
    int frame_end;
    
    int y_start;
    int y_end;
    int x_start;
    int x_end;
    
    int x_ind;
    int y_ind;
    int z_ind;
    int ind_srch;
    
    double dist;
    s_intra = 2*r_intra + 1;  // intra-frame search size
    s_inter = 2*r_inter + 1;  // inter-frame search size
    int n = (s_intra*s_intra*s_inter);
    
    int *ind_srch_matrix = sub2ind();
    double *xyDist_buffer = (double*)malloc(sizeof(double) * (s_intra*s_intra*s_inter));
    double *frameDist_buffer = (double*)malloc(sizeof(double) * (s_intra*s_intra*s_inter));
    
    double *x_buffer = malloc(sizeof(double) * s_intra*s_intra*s_inter);
    double *y_buffer = malloc(sizeof(double) * s_intra*s_intra*s_inter);
    double *frame_buffer = malloc(sizeof(double) * s_intra*s_intra*s_inter);
    double *dist_buffer = malloc(sizeof(double) * s_intra*s_intra*s_inter);
    
    //memset(xyDist_buffer, INFINITY, sizeof(double) * (s_intra*s_intra*s_inter));
    for (int i=0; i<(s_intra*s_intra*s_inter); i++) {
        frameDist_buffer[i] = INFINITY;
        xyDist_buffer[i] = INFINITY;
    }
    //memset(frameDist_buffer, INFINITY, sizeof(double) * (s_intra*s_intra*s_inter));
    
    for (int frame=0; frame<N_frame; frame++) {
        for (int y=0; y<= (size_y-(s_patch-1)); y=y+skip) {
            for (int x=0; x<= (size_x-(s_patch-1)); x=x+skip) {
                psuedo_patch(x, y, frame, patch, i_map_set_wiener3D);
                
                frame_start = (int) fmax(frame - r_inter, 0);
                frame_end = (int) fmin(frame + r_inter, N_frame);
                
                y_start = (int) fmax(y - r_intra, 0);
                y_end = (int) fmin(y + r_intra, size_y - (s_patch-1));

                x_start = (int) fmax(x - r_intra, 0);
                x_end = (int) fmin(x + r_intra, size_x - (s_patch-1));
                
                //absolute x position
                memset(x_buffer, -1.0, sizeof(double)*s_intra*s_intra*s_inter);
                
                memset(y_buffer, -1.0, sizeof(double) * (s_intra*s_intra*s_inter));
                //absolute y position
                
                //frame position
                memset(frame_buffer, -1.0, sizeof(double) * (s_intra*s_intra*s_inter));
                
                // distance between patches
                memset(dist_buffer, -1.0, sizeof(double) * (s_intra*s_intra*s_inter));
                
                //printf("frame_start: %d, frame_end: %d, y_start: %d, y_end: %d, x_start: %d, x_end: %d\n", frame_start, frame_end, y_start, y_end, x_start, x_end);
                
                for (int frame_srch=frame_start; frame_srch <= frame_end; frame_srch++) {
                    for (int y_srch=y_start; y_srch <= y_end; y_srch++) {
                        for (int x_srch=x_start; x_srch <= x_end; x_srch++) {
                            //printf("PATCH SEARCH\n\n");
                            psuedo_patch(x_srch, y_srch, frame_srch, T, i_map_set_wiener3D);
                            dist = dist_calc(T, patch);
                            
                            //printf("dist = %.5f\n", dist);
                            
                            y_ind = (r_intra)+(y_srch-y);
                            x_ind = (r_intra)+(x_srch-x);
                            z_ind = (r_inter)+(frame_srch-frame);
                            //printf("x_ind = %d, y_ind = %d, z_ind = %d\n", x_ind, y_ind, z_ind);
                            ind_srch = ind_srch_matrix[(x_ind*s_intra*s_inter) + (y_ind*s_inter) + z_ind];
                            
                            dist_buffer[ind_srch] = dist;
                            x_buffer[ind_srch] = x_srch;
                            y_buffer[ind_srch] = y_srch;
                            frame_buffer[ind_srch] = frame_srch;
                            
                            //printf("ind_srch = %d, x_srch = %d\n", ind_srch, x_srch);
                        }
                    }
                }
                
                //FIND SIMILAR PATCHES BABBBY
                xydistance(xyDist_buffer, x, x_buffer, y, y_buffer);
                sort(xyDist_buffer, x_buffer, y_buffer, frame_buffer, dist_buffer, 0, n-1);
                
                //printArray(x_buffer, n);
                //printArray(xyDist_buffer, n);
                
                framedistance(frameDist_buffer, frame, frame_buffer);
                sort(frameDist_buffer, x_buffer, y_buffer, frame_buffer, dist_buffer, 0, n-1);
                
                distDistance(dist_buffer);
                sort(dist_buffer, x_buffer, y_buffer, frame_buffer, dist_buffer, 0, n-1);
                
                sim_map(x_sim_map_set, y_sim_map_set, frame_sim_map_set, x, y, frame, x_buffer, y_buffer, frame_buffer);
                
            }
        }
    }
    
    free(patch);
    free(T);
    free(ind_srch_matrix);
    free(x_buffer);
    free(y_buffer);
    free(frame_buffer);
    free(dist_buffer);
    free(xyDist_buffer);
    free(frameDist_buffer);
    
    if (PRINT_SIM_MAPS) {
        printf("\nx_sim_map\n");
        for (int frame=0; frame<N_frame; frame++) {
            for (int z=0; z<N_sim; z++) {
                printf("\n\n----BIN NUMBER %d----\n\n", z);
                for (int y=0; y<size_y; y++) {
                    for (int x=0; x<size_x; x++) {
                        printf("%.0f\t",x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame]);
                    }
                    printf("\n");
                }
            }
        }
       
        printf("\ny_sim_map\n");
        for (int frame=0; frame<N_frame; frame++) {
            for (int z=0; z<N_sim; z++) {
                printf("\n\n----BIN NUMBER %d----\n\n", z);
                for (int y=0; y<size_y; y++) {
                    for (int x=0; x<size_y; x++) {
                        printf("%.0f\t",y_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame]);
                    }
                    printf("\n");
                }
            }
        }
        
        for (int frame=0; frame<N_frame; frame++) {
            for (int z=0; z<N_sim; z++) {
                //printf("BIN NUMBER %d\n\n\n", z);
                for (int x=0; x<size_x; x++) {
                    for (int y=0; y<size_y; y++) {
                        //printf("frame_sim_map_set[%d][%d][%d][%d] = %.1f\n", frame, x, x, y, frame_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame]);
                    }
                }
            }
        }

    }
}

void distDistance(double *dist_buffer) {

    for (int i=0; i<(s_intra*s_intra*s_inter); i++) {
        if (isnan(dist_buffer[i])) {
            dist_buffer[i] = INFINITY;
        }
    }
}

void sim_map(double *x_sim_map_set, double *y_sim_map_set, double *frame_sim_map_set, int x, int y, int frame, double *x_buffer, double *y_buffer, double *frame_buffer) {
    
    for (int z=0; z<N_sim; z++) {
        x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame] = x_buffer[z];
        y_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame] = y_buffer[z];
        frame_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame] = frame_buffer[z];
        
        //printf("x_sim_map_set[%d][%d][%d][%d] = %.1f\n", frame, x, x, y, x_sim_map_set[(x*size_y*N_sim*N_frame) + (y*N_sim*N_frame) + (z*N_frame) + frame]);
        
    }
    //printf("\n\n\n\n\n\n\n\n\n\n");
}

void xydistance(double *xyDist_buffer, int x, double *x_buffer, int y, double *y_buffer) {
    int size = s_intra*s_intra*s_inter;
    
    for (int i=0; i<size; i++) {
        if (!isnan(x_buffer[i])) {
            xyDist_buffer[i] = square(x-x_buffer[i]) + square(y-y_buffer[i]);
        } else {
            xyDist_buffer[i] = INFINITY;
        }
        //printf("xyDist_buffer[%d] = %.0f\n", i, xyDist_buffer[i]);
    }
    //printf("\n");
}

void framedistance(double *frameDist_buffer, int frame, double *frame_buffer) {
    int size = s_intra*s_intra*s_inter;
    
    for (int i=0; i<size; i++) {
        if (!isnan(frame_buffer[i])) {
            frameDist_buffer[i] = square(frame-frame_buffer[i]);
        } else {
            frameDist_buffer[i] = INFINITY;
        }
    }
}

int *sub2ind() {
    int *ind_srch_matrix = malloc(sizeof(int) * (s_intra*s_intra*s_inter));
    
    int count = 1;
    
    for (int i=0; i<s_intra; i++) {
        for (int j=0; j<s_intra; j++) {
            for (int k=0; k<s_inter; k++) {
                ind_srch_matrix[(i*s_intra*s_inter) + (j*s_inter) + k] = count;
                count++;
                
                //printf("ind_srch_matrix[%d][%d][%d] = %d\n", j, i, k,ind_srch_matrix[(j*s_intra*s_inter) + (i*s_inter) + k]);
            }
        }
    }
    
    return ind_srch_matrix;
    
}

double dist_calc(double *T, double *patch) {
    double dist = 0.0;
    
    for (int x=0; x<s_patch; x++) {
        for (int y=0; y<s_patch; y++) {
            dist = dist + square(patch[(x*s_patch) + y] - T[(x*s_patch) + y]);
            //printf("patch[%d][%d] = %.5f, T = %.5f\n", x, y, patch[(x*s_patch) + y], T[(x*s_patch) + y]);
        }
    }
    
    return dist;
}

void psuedo_patch(int start_x, int start_y, int frame, double *patch, double *i_map_set_wiener3D) {
    int count_x = 0;
    int count_y = 0;
    for (int x=start_x; x<(start_x+s_patch); x++) {
        count_y = 0;
        for (int y=start_y; y<(start_y+s_patch); y++) {
            patch[(count_x*s_patch) + count_y] = i_map_set_wiener3D[(x*s_patch) + y];
                
                //printf("patch[%d][%d] = %.30f\n", count_x, count_y, patch[(count_x*s_patch) + count_y]);
                count_y++;
        }
        count_x++;
    }
}

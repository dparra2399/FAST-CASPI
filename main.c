#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"
#include <time.h>
//#include"mat.h"
//#include"matrix.h"

//Globals (scenes)
 double t_bin;
 int N_bin; //default 1024
int N_bin_half;
 int N_sim;
 int s_patch;
 int r_intra;
 int r_inter;
 double FWHM;
 double d_min;
 double d_max;

//filtering
 int skip; //filtering skip size, 1 means dense filtering
 int s_intra; //search size within a frame

//itensity-based filtering
 int pseudo_int; //1: pseudo intensity, 0: true insensity
 int IBF; //1: yes, 0: no
 double th_IBF3D; //default
 double th_IBF4D; //default
 int max_int_prct; //100: normalize intensity with max

//system
int N_cycle; //laser cycle number
double c; //light speed

//misc
 double d_range; //measurable depth range
 double sigma_t; //std for the given FWHM in time doamin
 double sigma_f;//stf for the given FWHM in freq domain
 double bin_f; //bin size in freq domain
 int N_sig_f; //number of singla bins freq domain

 double sigma_bin;
 int bin_max;

int size_y; //209
int size_x; //167
int N_frame;



int main()
{	
	double *hst_map_set; //contiguous array
    clock_t t;
    _init();
    //mxArray *temp_set;
	//mxArray *cell;
	//MATFile *pmat = matOpen("Art_hst.mat","r");
	
	
	//allocate memory for hst_map_set
	hst_map_set = malloc((size_y*size_x*N_bin*N_frame)*sizeof(double));
	//TEMPORARY HSt_MAP_SET
	for (int x=0; x<size_x; x++) {
		for (int y=0; y<size_y; y++) {
			for (int z=0; z<N_bin; z++) {
				for (int n=0; n<N_frame; n++) {
                    hst_map_set[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + n] = 1.0; //(double) ((x+1)*(y+1)*(z+1)*(n+1) /  (double)(size_x*size_y*N_bin*N_frame));
                   // printf("hst_map_set[%d][%d][%d][%d] = %.15f\n", x, y, z, n, hst_map_set[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + n]);
				}
			}
		}
	}
	/*
	if (pmat == NULL) {
		printf("Error creating MATFile\n");
		return(EXIT_FAILURE);
	} else {
		printf("Sucessful file read\n");
	}
	temp_set = matGetVariable(pmat, "hst_map_set");
	if (temp_set == NULL) {
		printf("Error using matGetVariable\n");
		return(EXIT_FAILURE);
	} else {
		printf("Sucessful matGetVariable\n");
	}
	
	double *ptr =  mxGetDoubles(temp_set);
	if (ptr == NULL) {
		printf("Error using mxGetDoubles\n");
		return(EXIT_FAILURE);
	} else {
		printf("Sucessful matGetDoubles\n");
	}*/
    
	double *gauss = normpdf();
	
    if (pseudo_int) {
        double *i_map_hst = pseudo_intensity(hst_map_set);
    }
    t = clock();
	double *flux_map_noisy = coates_correction(hst_map_set);
    t = clock() - t;
    double time_taken_coates = ((double)t)/CLOCKS_PER_SEC;

    printf("coates correction took %f seconds to execute \n", time_taken_coates);
    
    
    fftw_complex *FLUX_map_set_noisy1D_full = fft(flux_map_noisy);

    fftw_complex *FLUX_map_set_noisy1D = fft_halfbin(FLUX_map_set_noisy1D_full);
    
    double* i_map_set_wiener3D = (double*)calloc((size_y*size_x*N_frame), sizeof(double));
    
    t = clock();
    if (pseudo_int) {
        intensities(i_map_set_wiener3D, FLUX_map_set_noisy1D);
    }
    t = clock() - t;
    double time_taken_intent = ((double)t)/CLOCKS_PER_SEC;
    printf("intensity filtering took %f seconds to execute \n", time_taken_intent);
    
    double *x_sim_map_set = (double*)calloc((size_y*size_x*N_sim*N_frame), sizeof(double));
    double *y_sim_map_set = (double*)calloc((size_y*size_x*N_sim*N_frame), sizeof(double));
    double *frame_sim_map_set = (double*)calloc((size_y*size_x*N_sim*N_frame), sizeof(double));
    
    t = clock();
    sim_patches(i_map_set_wiener3D, x_sim_map_set, y_sim_map_set, frame_sim_map_set);
    t = clock() - t;
    double time_taken_sim_patches = ((double)t)/CLOCKS_PER_SEC;
    printf("similar patched took %f seconds to execute \n", time_taken_sim_patches);
    
    
    double* MSE_map_set = malloc(sizeof(double)*(size_x*size_y*N_frame));
    double* IBFmask4D_set = (double*)calloc((size_y*size_x*N_frame), sizeof(double));
    double* flux_map_set_ht4D = (double*)calloc((size_y*size_x*N_bin*N_frame), sizeof(double));
    t=clock();
    
    ht4D(flux_map_set_ht4D, FLUX_map_set_noisy1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, i_map_set_wiener3D, MSE_map_set, IBFmask4D_set);
    
    if (pseudo_int) {
        double *i_map_ht4D = pseudo_intensity(flux_map_set_ht4D);
    }
    
    t= clock() - t;
    double time_taken_ht4D_patches = ((double)t)/CLOCKS_PER_SEC;
    printf("hard thresholding 4D took %f seconds to execute \n", time_taken_ht4D_patches);
    
    double* flux_map_set_wiener4D = (double*)calloc((size_y*size_x*N_bin*N_frame), sizeof(double));
    
    t = clock();
    
    wf4D(flux_map_set_wiener4D, flux_map_set_ht4D, FLUX_map_set_noisy1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, MSE_map_set);
    
    if (pseudo_int) {
        double *i_map_wiener4D = pseudo_intensity(flux_map_set_wiener4D);
    }
    
    t= clock() - t;
    double time_taken_wf4D_patches = ((double)t)/CLOCKS_PER_SEC;
    printf("wiener filtering 4D took %f seconds to execute \n", time_taken_wf4D_patches);
    
    double total_time = time_taken_coates + time_taken_intent + time_taken_sim_patches + time_taken_ht4D_patches + time_taken_wf4D_patches;
    printf("total time took = %f seconds \n", total_time);
    
	printf("Sucessful End\n");
	return 0;
}

void intensities(double *i_map_set_wiener3D, fftw_complex *FLUX_map_set_noisy1D) {
    fftw_complex *FLUX_map_noisy1D;
    double *FLUX_map_1D;
    double *i_map_1D;
    double *i_map_ht3D;
    double *i_map_wiener3D;
    
    double* MSE_map = malloc(sizeof(double)*(size_x*size_y)); //map of mean squared noise
    double* flux_map_ht3D = (double*)calloc((size_y*size_x*N_bin), sizeof(double)); //recovered flux map
    double* IBFmask3D = (double*)calloc((size_y*size_x), sizeof(double));
    double* flux_map_wiener3D = (double*)calloc((size_y*size_x*N_bin), sizeof(double));
    
    int* coeff_hard = malloc(sizeof(int)*(s_patch*s_patch*N_sig_f));
    for (int frame=0; frame<N_frame; frame++) {
        FLUX_map_noisy1D = getNthFrame_fftw(FLUX_map_set_noisy1D, frame);
        FLUX_map_1D = ifft(FLUX_map_noisy1D);
        
        i_map_1D = sum_dim(FLUX_map_1D, N_bin);
        divide(i_map_1D, prctile(i_map_1D));
        
        
        clock_t t = clock();
        
        ht3D(FLUX_map_noisy1D, i_map_1D, MSE_map, flux_map_ht3D, IBFmask3D);
        
        t= clock() - t;
        double time_taken_ht3D_patches = ((double)t)/CLOCKS_PER_SEC;
        printf("hard thresholding 3D took %f seconds to execute \n", time_taken_ht3D_patches);
        
        i_map_ht3D = sum_dim(flux_map_ht3D, N_bin);
        divide(i_map_ht3D, prctile(i_map_ht3D));
        
        wf3D(flux_map_wiener3D, flux_map_ht3D, FLUX_map_noisy1D, MSE_map);
        
        i_map_wiener3D = sum_dim(flux_map_wiener3D, N_bin);
        divide(i_map_wiener3D, prctile(i_map_wiener3D));
        i_map_wiener_frame(frame, i_map_set_wiener3D, i_map_wiener3D);
        
        if (DEBUG) {
            /*for (int z=0; z<N_bin; z++) {
                printf("\n\n BIN NUMER %d \n\n", z);
                for (int x=0; x<size_x; x++) {
                    for (int y=0; y<size_y; y++) {
                        printf("flux_map_ht3d[%d][%d][%d] = %.15f\n", x, y, z, flux_map_ht3D[(x*size_y*N_bin) + (y*N_bin) + z]);
                    }
                }
            }
            
            for (int x=0; x<size_x; x++) {
                for (int y=0; y<size_y; y++) {
                    printf("IBFmask3D[%d][%d] = %.15f\n", x, y, IBFmask3D[(x*size_y) + y]);
                }
            }
            
            for (int x=0; x<size_x; x++) {
                for (int y=0; y<size_y; y++) {
                    printf("MSE_map[%d][%d] = %.15f\n", x, y, MSE_map[(x*size_y) + y]);
                }
            }*/
            
            for (int z=0; z<N_bin; z++) {
                printf("\n\n BIN NUMER %d \n\n", z);
                for (int x=0; x<size_x; x++) {
                    for (int y=0; y<size_y; y++) {
                        printf("flux_map_wiener3d[%d][%d][%d] = %.15f\n", x, y, z, flux_map_wiener3D[(x*size_y*N_bin) + (y*N_bin) + z]);
                    }
                }
            }
            
        }
    }
}
void i_map_wiener_frame(int frame, double *i_map_set_wiener3D, double *i_map_wiener3D) {
    
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            i_map_set_wiener3D[(x*size_y*N_frame) + (y*N_frame) + frame] = i_map_wiener3D[(x*size_y) + y];
        }
    }
}

double sum_1d(double *hst) {
	double sum = 0.0;
	for (int i = 0; i<N_bin; i++) {
		sum = hst[i] + sum;
	}
	return sum;
}

double* normpdf() {
	int size = bin_max;
	double* pdf = malloc(size*sizeof(double));
	double m = bin_max/2;
	double s = sigma_bin;
	for (int i=1; i<size+1; i++) {
		pdf[i-1] =(1/(s*sqrt(2* PI))) * exp(-0.5*pow((i-m)/s,2.0));
		//calculate pdf using mean and sd
	}
	pdf[size] = 0;
	return pdf;
}

double* pseudo_intensity(double *hst_map_set) {
	double *hst_map;
	double *i_map_hst;
    for (int i=0; i<N_frame; i++) {
        hst_map = getNthFrame(hst_map_set, i);
        i_map_hst = sum_dim(hst_map, N_bin);
        divide(i_map_hst, prctile(i_map_hst));
    }
	//printf("i_map_hst[5][5] = %.1f\n", i_map_hst[5*size_y + 5]);
	return i_map_hst;
}

double *getNthFrame(double *hst_map_set, int frame) {
	double *hst_map = malloc((size_x*size_y*N_bin) * sizeof(double));
	for (int x=0; x<size_x; x++) {
		for (int y=0; y<size_y; y++) {
			for (int z=0; z<N_bin; z++) {
				hst_map[(x*size_y*N_bin) + (y*N_bin) + z] = hst_map_set[(x*size_y*N_bin*N_frame) + (y*N_bin*N_frame) + (z*N_frame) + frame];
			}
		}
	}
	return hst_map;
}

void divide(double *i_map_hst, double prctile) {
	for (int i=0; i<size_x; i++) {
		for (int j=0; j<size_y; j++) {
			i_map_hst[(i*size_y) + j] = i_map_hst[(i*size_y) + j] / prctile;
			//printf("i_map_whatever[%d][%d] = %.5f\n", i, j, i_map_hst[(i*size_y) + j]);
		}
	}
}

double prctile(double *sum_map_hst) {
	double max = sum_map_hst[0];
	if (max_int_prct == 100) {
		//get max value of hst_map
		//TEMPORARY AS WELL
		for (int i=0; i<size_x; i++) {
			for (int j=0; j<size_y; j++) {
				if (sum_map_hst[(i*size_y) + j] > max) {
				max = sum_map_hst[(i*size_y) + j];
				}
			}
		}
	}
	return max;
}


double* sum_dim(double *hst_map, int bin_size)  {
	double *i_map_hst = malloc((size_x*size_y)*sizeof(double));
	double sum;
	
	
	for (int x=0; x<size_x; x++) {
		for (int y=0; y<size_y; y++) {
			sum = 0.0;
            for (int z=0; z<bin_size; z++) {
				sum = sum + hst_map[(x*size_y*bin_size) + (y*bin_size) + z];
			}
			i_map_hst[(x*size_y) + y]	= sum;
			//printf("i_map_hst_whatever[%d][%d] = %.5f\n", x, y, i_map_hst[x*size_y + y]);
		}
	}
	return i_map_hst;
}
	

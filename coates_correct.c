#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"

double *coates_correction(double *hst_map_set) {
    double *flux_map_noisy; //(
    double *hst_map;
    for (int i=0; i<N_frame; i++) {
        hst_map = getNthFrame(hst_map_set, i);
        flux_map_noisy = coates_grp(hst_map);
        if (pseudo_int) {
            double *i_map_coates = sum_dim(flux_map_noisy, N_bin);
            divide(i_map_coates, prctile(i_map_coates));
        }
    }
    return flux_map_noisy;
}

double *coates_grp(double *hst_map) {
    double *hst = malloc(N_bin*sizeof(double));
    double *flux_hat_grp = malloc((size_y*size_x*N_bin)*sizeof(double));
    double *flux_hat;
    for (int x=0; x<size_x; x++) {
        for (int y=0; y<size_y; y++) {
            //squeeze hst_map[x][y]
            for (int z=0; z<N_bin; z++) {
                hst[z] = hst_map[(x*size_y*N_bin) + (y*N_bin) + z];
                //printf("hst_tmp[%d] = %.5f\n", z, hst[z]);
            }
            flux_hat = coates(hst);
            for (int z=0; z<N_bin; z++) {
                flux_hat_grp[(x*size_y*N_bin) + (y*N_bin) + z] = flux_hat[z];
                //printf("Flux_hat_grp[%d] = %.5f\n", ((x*size_y*N_bin) + (y*N_bin) + z), flux_hat_grp[(x*size_y*N_bin) + (y*N_bin) + z]);
            }
        }
    }
    return flux_hat_grp;
}

double *coates(double *hst) {
    double *flux_hat;
    double sum = sum_1d(hst);
    //printf("Sum should be 1024: %.5f\n", sum);
    if (sum >= N_cycle) {
        N_cycle = sum + 1;
    }
    
    flux_hat = flux_estimate(hst);
    return flux_hat;
}

/**CAN BE OPTIMIZED**/
double *flux_estimate(double *hst) {
    double *numer = malloc(N_bin * sizeof(double));
    double *denom = malloc(N_bin * sizeof(double));
    double *flux_hat = malloc(N_bin * sizeof(double));
    double cumsum = 0.0;
    numer[0] = 0.0;
    
    for (int i=0; i<N_bin; i++) {
        cumsum = hst[i] + cumsum;
        denom[i] = cumsum;
    }
    
    for (int i=1; i<N_bin; i++) {
        numer[i] = denom[i-1];
    }
    
    for (int i=0; i<N_bin; i++) {
        flux_hat[i] = log((N_cycle - numer[i])/(N_cycle - denom[i]));
        if (flux_hat[i] < 0) {
            flux_hat[i] = 0;
        }
        //printf("flux_hat[%d] = %.5f\n", i, flux_hat[i]);
    }
    
    free(numer);
    free(denom);
    return flux_hat;
}

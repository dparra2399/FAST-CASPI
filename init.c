#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"


//Initialize all global vars
void _init() {
    t_bin = 80E-12;
    N_bin = 1024; //default 1024
    N_bin_half = round(((double)N_bin-1)/2)+1;
    N_sim = 10;
    s_patch = 2; //default 8
    r_intra = 10;
    r_inter = 0;
    FWHM = 400E-12;
    d_min = 1.40;
    d_max = 2.17;

   //filtering
    skip = 1; //filtering skip size, 1 means dense filtering
    s_intra = 2*r_intra + 1; //search size within a frame

   //itensity-based filtering
    pseudo_int = 1; //1: pseudo intensity, 0: true insensity
    IBF = 1; //1: yes, 0: no
    th_IBF3D = 0.8; //default
    th_IBF4D = 0.9; //default
    max_int_prct = 100; //100: normalize intensity with max

   //system
   N_cycle = 1000; //laser cycle number
   c = 3E+8; //light speed
    
   size_y = 206; //209
   size_x = 167; //167
   N_frame = 1;

    
    d_range = c*N_bin*t_bin/2; //measurable depth range
    sigma_t = FWHM/(2*sqrt(2*log(2))); //std for the given FWHM in time doamin
    sigma_f = 1/(2*M_PI*sigma_t) ;//stf for the given FWHM in freq domain
    bin_f = 1/(t_bin*N_bin); //bin size in freq domain
    N_sig_f = ceil(3*sigma_f/bin_f); //number of singla bins freq domain

    sigma_bin = sigma_t/t_bin;
    bin_max = 10*round(sigma_bin);
}


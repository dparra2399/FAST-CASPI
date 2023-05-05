#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "caspi.h"
#include <time.h>
#include <string.h>
#include <math.h>


// Quick sort in C
//NOT DAVID PARRAS CODE THIS IS FROM THUS LINK:
//https://www.programiz.com/dsa/quick-sort

// function to swap elements
void swap(double *a, double *b) {
  double t = *a;
  *a = *b;
  *b = t;
}

// function to find the partition position
double partition(double array[], double x_buffer[], double y_buffer[], double frame_buffer[], double dist_buffer[], int low, int high) {
  
  // select the rightmost element as pivot
  double pivot = array[high];
  
  // pointer for greater element
  int i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swap(&array[i], &array[j]);
      swap(&x_buffer[i], &x_buffer[j]);
      swap(&y_buffer[i], &y_buffer[j]);
      swap(&frame_buffer[i], &frame_buffer[j]);
      swap(&dist_buffer[i], &dist_buffer[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swap(&array[i + 1], &array[high]);
  swap(&x_buffer[i + 1], &x_buffer[high]);
  swap(&y_buffer[i + 1], &y_buffer[high]);
  swap(&frame_buffer[i + 1], &frame_buffer[high]);
  swap(&dist_buffer[i + 1], &dist_buffer[high]);
  
  // return the partition point
  return (i + 1);
}

void sort(double array[], double x_buffer[], double y_buffer[], double frame_buffer[], double dist_buffer[], int low, int high) {
    if (low < high) {
    
        // find the pivot element such that
        // elements smaller than pivot are on left of pivot
        // elements greater than pivot are on right of pivot
        double pi = partition(array, x_buffer, y_buffer, frame_buffer, dist_buffer, low, high);
    
        // recursive call on the left of pivot
        sort(array, x_buffer, y_buffer, frame_buffer, dist_buffer, low, pi - 1);
    
        // recursive call on the right of pivot
        sort(array, x_buffer, y_buffer, frame_buffer, dist_buffer, pi + 1, high);
    }
}

// function to print array elements
void printArray(double *array, int size) {
  for (int i = 0; i < size; i++) {
    printf("array[%d] = %.0f\n", i, array[i]);
  }
  printf("\n");
}

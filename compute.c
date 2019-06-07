#include <unistd.h>
#include <math.h>

#include "compute.h"

void compute_edge_A(unsigned int tea) {
  usleep(tea);
}

void compute_edge_B(unsigned int teb) {
  usleep(teb);
}

void compute_middle_A(unsigned int tma) {
  usleep(tma);
}

void compute_middle_B(unsigned int tmb) {
  usleep(tmb);
}


void compute_edgeA(double* bufA, int bx, int by) {
  double tmp = 0.0;

  // compute full columns on the left and right
  for(int j = 1; j < by+1; ++j) {
    for(int i = 1; i < bx+1; i+=bx-1) {
      tmp  += bufA[indx(i,j)]/2.0 + (bufA[indx(i-1,j)] +
                                     bufA[indx(i+1,j)] +
                                     bufA[indx(i,j-1)] +
                                     bufA[indx(i,j+1)])/4.0/2.0;
    }
  }
}


void compute_middleA(double* bufA, int bx, int by) {
  double tmp = 0.0;

  // compute full columns on the left and right
  for(int j = 2; j < by; ++j) {
    for(int i = 2; i < bx; ++i) {
      for(int k=0; k< 1; k++){

        bufA[indx(i,j)] += bufA[indx(i,j)]/2.0 + (double)k*(bufA[indx(i-1,j)] +
                                                            bufA[indx(i+1,j)] +
                                                            bufA[indx(i,j-1)] +
                                                            bufA[indx(i,j+1)])/4.0/2.0;

        //bufA[indx(i,j)]= pow(bufA[indx(i,j)], 0.245);

      }
      //bufA[indx(i,j)]=  bufA[indx(i,j)]/1000.0;
    }
    //bufA[indx(bx/3,j)] =  bufA[indx(bx/2,j)];
  }
}


void compute_edgeB(double* bufB, int bx, int by) {
  double tmp = 0.0;

  // compute full columns on the left and right
  for(int j = 1; j < by+1; ++j) {
    for(int i = 1; i < bx+1; i+=bx-1) {
      tmp  += bufB[indx(i,j)]/2.0 + (bufB[indx(i-1,j)] +
                                     bufB[indx(i+1,j)] +
                                     bufB[indx(i,j-1)] +
                                     bufB[indx(i,j+1)])/4.0/2.0;
    }
  }
}

void compute_middleB(double* bufB, int bx, int by) {
  double tmp = 0.0;

  // compute full columns on the left and right
  for(int j = 2; j < by; ++j) {
    for(int i = 2; i < bx; ++i) {
      for(int k=0; k< 1; k++){

        bufB[indx(i,j)] += bufB[indx(i,j)]/2.0 + (double)k*(bufB[indx(i-1,j)] +
                                                            bufB[indx(i+1,j)] +
                                                            bufB[indx(i,j-1)] +
                                                            bufB[indx(i,j+1)])/4.0/2.0;

        //bufB[indx(i,j)]= pow(bufB[indx(i,j)], 0.245);

      }
      //bufB[indx(i,j)] =  bufB[indx(i,j)]/1000.0;
    }
    //bufB[indx(bx/3,j)] =  bufB[indx(bx/2,j)];
  }
}

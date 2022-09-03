// -*- C++ -*-

#include "pomp.h"

// stratified resampling
void systematic_resampling (int *n, double *w, int *perm) 
{
  double u, du;
  int i, j;

  for (j = 1; j < *n; j++) w[j] += w[j-1];
  for (j = 0; j < *n; j++) w[j] /= w[*n-1];

  GetRNGstate();
  du = 1.0 / ((double) *n);
  u = runif(-du,0);
  PutRNGstate();

  for (i = 0, j = 0; j < *n; j++) {
    u += du;
    while (u > w[i]) i++;
    perm[j] = i + 1; // must use 1-based indexing for compatibility with R!
  }
}


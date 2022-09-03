// -*- C++ -*-

#include "pomp.h"

// B-spline basis
void bspline (double *y, double *x, int *nx, int *i, int *p, double *knots, int *nknots)
{
  int j;
  double a, b;
  double y1[*nx], y2[*nx];
  int i2, p2;

  if ((*i < 0) || (*p < 0) || (*i + *p >= *nknots-1)) 
    error("bad arguments in bspline");

  if (*p == 0) {
    for (j = 0; j < *nx; j++) {
      y[j] = ((knots[*i] <= x[j]) && (x[j] < knots[*i+1]));
    }
  } else {
    i2 = *i+1;
    p2 = *p-1;
    bspline(y1,x,nx,i,&p2,knots,nknots);
    bspline(y2,x,nx,&i2,&p2,knots,nknots);
    for (j = 0; j < *nx; j++) {
      a = (x[j]-knots[*i]) / (knots[*i+*p]-knots[*i]);
      b = (knots[i2+*p]-x[j]) / (knots[i2+*p]-knots[i2]);
      y[j] = a * y1[j] + b * y2[j];
    }
  }
}

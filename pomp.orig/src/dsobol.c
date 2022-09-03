// -*- C++ -*-

#include <R.h>
#include <Rmath.h>

void F77_NAME(insobl)(int *, int *, int *, int *);
void F77_NAME(gosobl)(double *);

void dsobol (double *data, int *dim, int *n)
{
  int flag[2], taus = 0;
  int k;

  F77_CALL(insobl)(flag, dim, n, &taus);

  if (!flag[0])
    error("dimension is not OK in sobol");

  if (!flag[1]) 
    error("number of points requested is not OK in sobol");

  for (k = 0; k < *n; k++)
    F77_CALL(gosobl)(data + k*(*dim));

}

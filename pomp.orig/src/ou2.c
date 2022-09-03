// -*- C++ -*-

#include <R.h>
#include <Rmath.h>

struct lookup_table {
  int length, width;
  int index;
  double *x;
  double **y;
};

// prototypes
double expit (double x);
double logit (double x);
double transmission (int dim, double *basis, double *coef, double trend, double t);
void normal_dmeasure (int *n, int *index, double *X, double *y, double *f);
void normal_rmeasure (int *n, int *index, double *X, double *obs);
void ou2 (double *x, 
	  double alpha1, double alpha2, double alpha3, double alpha4, 
	  double sigma1, double sigma2, double sigma3);
void ou2_adv (double *X, int *nvar, int *np, int *stateindex, int *parindex);

// simple 2D Ornstein-Uhlenbeck process
void ou2 (double *x, 
	  double alpha1, double alpha2, double alpha3, double alpha4, 
	  double sigma1, double sigma2, double sigma3)
{
  double eps[2], y[2];

  if (!(R_FINITE(x[0]))) return;
  if (!(R_FINITE(x[1]))) return;
  if (!(R_FINITE(alpha1))) return;
  if (!(R_FINITE(alpha2))) return;
  if (!(R_FINITE(alpha3))) return;
  if (!(R_FINITE(alpha4))) return;
  if (!(R_FINITE(sigma1))) return;
  if (!(R_FINITE(sigma2))) return;
  if (!(R_FINITE(sigma3))) return;

  eps[0] = rnorm(0,1);
  eps[1] = rnorm(0,1);

  y[0] = alpha1*x[0]+alpha3*x[1]+sigma1*eps[0];
  y[1] = alpha2*x[0]+alpha4*x[1]+sigma2*eps[0]+sigma3*eps[1];

  x[0] = y[0];
  x[1] = y[1];
}

#define STATE      (x[stateindex[0]])
#define ALPHA1     (x[parindex[0]])
#define ALPHA2     (x[parindex[1]])
#define ALPHA3     (x[parindex[2]])
#define ALPHA4     (x[parindex[3]])
#define SIGMA1     (x[parindex[4]])
#define SIGMA2     (x[parindex[5]])
#define SIGMA3     (x[parindex[6]])

// advance the matrix of particles one time unit
void ou2_adv (double *X, int *nvar, int *np, int *stateindex, int *parindex)
{
  double *x;
  int k;
  GetRNGstate();		// initialize R's pseudorandom number generator
  for (k = 0; k < *np; k++) {
    x = &X[*nvar*k];		// get address of k-th particle
    ou2(&STATE,-exp(ALPHA1),exp(ALPHA2),-exp(ALPHA3),-exp(ALPHA4),exp(SIGMA1),SIGMA2,exp(SIGMA3)); // advance particle
  }
  PutRNGstate();		// finished with R's random number generator
}

#undef STATE
#undef ALPHA1
#undef ALPHA2
#undef ALPHA3
#undef ALPHA4
#undef SIGMA1
#undef SIGMA2
#undef SIGMA3

#define STATE (x[index[0]])
#define TAU   (x[index[1]])

void normal_dmeasure (int *n, int *index, double *X, double *y, double *f) {
  int p, nv = n[0], np = n[1];
  double *x, v, tol = 1.0e-18;
  for (p = 0; p < np; p++) {
    x = &X[nv*p];
    v = exp(TAU);
    if (R_FINITE(v)) {
      f[p] = dnorm(*y,STATE,v+tol,0)+tol;
    } else {
      f[p] = tol;
    }
  }
}

void normal_rmeasure (int *n, int *index, double *X, double *obs) {
  int p, nv = n[0], np = n[1];
  double *x, v, tol = 1.0e-18;
  GetRNGstate();
  for (p = 0; p < np; p++) {
    x = &X[nv*p];
    v = exp(TAU);
    if (R_FINITE(v)) {
      obs[p] = rnorm(STATE,v+tol);
    } else {
      obs[p] = R_NaReal;
    }
  }
  PutRNGstate();
}

#undef STATE
#undef TAU

// linear interpolation on a lookup table
void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt)
{
  int flag = 0;
  int j;
  double e;
  tab->index = findInterval(tab->x,tab->length,x,TRUE,TRUE,tab->index,&flag);
  if (flag != 0)              // we are extrapolating
    warning("table_lookup: extrapolating (flag %d) at %le", flag, x);
  e = (x - tab->x[tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
  for (j = 0; j < tab->width; j++) {
    y[j] = tab->y[j][tab->index] * e + tab->y[j][tab->index-1] * (1-e);
  }
  if (dydt != 0) {
    for (j = 0; j < tab->width; j++) {
      dydt[j] = (tab->y[j][tab->index] - tab->y[j][tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
    }
  }
}

double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

double logit (double x) {
  return log(x/(1-x));
}


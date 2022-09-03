// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

# One may find these two functions from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);


// prototypes
double expit (double x);
double logit (double x);
void euler_multinomial (int ntrans, double n, double *rate, double *trans, double dt);
void negbin_rmeasure (int *n, int *index, double *X, double *cases);
void negbin_dmeasure (int *n, int *index, double *X, double *y1, double *y2, double *f);
void basic_seir_b (double *Sh, double *Eh, double *Ih, double *Th, double *Dh, double *Rh, double *Rh1, double *Rh2, double *Th1, double *Th2,  double *C,  double *CC,
                   double beta, double sigmah, double theta, double gammah, double r, double kappah, double psi,
                   double dt, double v1, double eta);
void basic_sirs_pois (double *X, double *t1, double *t2, int *nstep,
			double *t_start, double *t_end,
		      int *nvar, int *np,
		      int *stateindex, int *parindex,  int *seas_dim, double *vac, double *impt);

//expp<-function(gamma,n1,dt){365*dt*n1/(1.0-exp(-n1*gamma*dt))}
//logg<-function(gamma,n1,dt){-log(1.0-(365*dt*n1)/gamma)/(n1*dt)}

double expp(double gamma, double n1, double dt)   //gamma to infectious period in days
{
return(365.0*dt*n1/(1.0-exp(-n1*gamma*dt)));
}
double logg(double gamma, double n1, double dt)  //infectious period in days to gamma
{
return(-log(1.0-(365.0*dt*n1)/gamma)/(n1*dt));
}

void basic_seir_b (double *Sh, double *Eh, double *Ih, double *Th, double *Dh, double *Rh,  double *Rh1, double *Rh2, double *Th1, double *Th2, double *C, double *CC,
                   double beta, double sigmah, double theta, double gammah, double r, double kappah,
                   double psi, double dt, double v1, double eta)

{
  double rate[8], trans[11], n;
  if (!(R_FINITE(Sh[0]))) return;
  if (!(R_FINITE(Eh[0]))) return;
  if (!(R_FINITE(Ih[0]))) return;
  if (!(R_FINITE(Th[0]))) return;
  if (!(R_FINITE(Rh[0]))) return;
  if (!(R_FINITE(Rh1[0]))) return;
  if (!(R_FINITE(Rh2[0]))) return;
  if (!(R_FINITE(Th1[0]))) return;
  if (!(R_FINITE(Th2[0]))) return;
  if (!(R_FINITE(Dh[0]))) return;
  if (!(R_FINITE(C[0]))) return;
  if (!(R_FINITE(CC[0]))) return;

  if (!(R_FINITE(beta))) return;
  if (!(R_FINITE(eta))) return;
  if (!(R_FINITE(kappah))) return;
  if (!(R_FINITE(sigmah))) return;
  if (!(R_FINITE(theta))) return;
  if (!(R_FINITE(gammah))) return;
  if (!(R_FINITE(v1))) return;

// if rounding is not done elsewhere, we must do rounding here!
// the state variables are assumed to be positive integers!
   Rh[0] += rint(Sh[0]*v1);
   Sh[0] -= rint(Sh[0]*v1);
    if(Sh[0]<0)Sh[0]=0;

	n = Sh[0]+Eh[0]+Ih[0]+Th[0]+Rh[0]+Rh1[0]+Rh2[0]+Th1[0]+Th2[0];
  rate[0] = beta*Ih[0]/n; 			// infection rate 
  rate[1] = sigmah; 			// from Eh to Ah 
  rate[2] = theta*gammah;      		// from Ih to Th 
  rate[3] = (1-theta)*gammah;   	// from Ih to Rh
  rate[4] = (1-theta*r)*3*kappah;   	         	// from Th to recovery 
  rate[5] = r*theta*3*kappah;   	         	// from Th to death 
  rate[6] = 3*psi;
  rate[7] = rate[4]+rate[5];
  euler_multinomial(1,Sh[0],&rate[0],&trans[0],dt);
  euler_multinomial(1,Eh[0],&rate[1],&trans[1],dt);
  euler_multinomial(2,Ih[0],&rate[2],&trans[2],dt);
  euler_multinomial(1,Th[0],&rate[7],&trans[9],dt);
  euler_multinomial(1,Th1[0],&rate[7],&trans[10],dt);
  euler_multinomial(2,Th2[0],&rate[4],&trans[4],dt);
  euler_multinomial(1,Rh[0],&rate[6],&trans[6],dt);
  euler_multinomial(1,Rh1[0],&rate[6],&trans[7],dt);
  if(Rh1[0]<1)trans[7]=0;
  euler_multinomial(1,Rh2[0],&rate[6],&trans[8],dt);
  if(Rh2[0]<1)trans[8]=0;

  Sh[0]  +=  (trans[8] - trans[0]);
  Eh[0]  +=  (trans[0] - trans[1]);
  Ih[0]  +=  (trans[1] - trans[2] - trans[3]);
  Th[0]  += (trans[2] - trans[9]);
  Th1[0]  += (trans[9] - trans[10]);
  Th2[0]  +=  (trans[10] - trans[4] - trans[5]);
  Rh[0]  +=  (trans[3] + trans[4] - trans[6]);
  Rh1[0]  += (trans[6] - trans[7]);
  Rh2[0]  += (trans[7] - trans[8]);
  Dh[0]  +=  trans[5] + Rf_rpois(eta);
  C[0]   +=  (trans[2]+trans[3]);
  CC[0]  +=  (trans[8]);
 while(Eh[0]<1){Eh[0]+=1; Sh[0]-=1;}
 while(Ih[0]<1){Ih[0]+=1; Sh[0]-=1;}

}


// void basic_seir_b (double *Sh, double *Eh, double *Ah, double *Ih, double *Th, double *Dh, double *Rh,  double *C,
//                   double beta, double sigmah, double theta, double delta, double gammah, double psi, double kappah,
//                   double dt)

#define SIGMAH   (x[parindex[0]])
#define THETA    (x[parindex[1]])
#define P0       (x[parindex[2]])
#define PHI      (x[parindex[3]])
#define GAMMAH   (x[parindex[4]])
#define PSI      (x[parindex[5]])
#define KAPPAH   (x[parindex[6]])
#define ETA      (x[parindex[7]])
#define LOGBETA  (x[parindex[8]])
#define BPOP     (x[parindex[9]])
#define ALPHA    (x[parindex[10]])

#define BSH      (x[stateindex[0]])
#define BEH      (x[stateindex[1]])
#define BIH      (x[stateindex[2]])
#define BTH      (x[stateindex[3]])
#define BRH      (x[stateindex[4]])
#define BRH1      (x[stateindex[5]])
#define BRH2      (x[stateindex[6]])
#define TH1      (x[stateindex[7]])
#define TH2      (x[stateindex[8]])
#define BDH      (x[stateindex[9]])
#define CC       (x[stateindex[10]])
#define BCA      (x[stateindex[11]])


void basic_sirs_pois (double *X, double *t1, double *t2, int *nstep,
		      double *t_start, double *t_end,
		      int *nvar, int *np,
		      int *stateindex, int *parindex, int *seas_dim, double *vac, double *impt)
{
  double dt = (*t2 - *t1) / ((double) *nstep);
  float beta,r,cut;
  double t, *x;
  float xa[*seas_dim/2+1] ,ya[*seas_dim/2+1], y2a[*seas_dim/2+1], za[*seas_dim/2+1], z2a[*seas_dim/2+1], slope;
  int j, p, step, day;
  
  GetRNGstate();		// initialize R's pseudorandom number generator


  for (p = 0; p < *np; p++) {	// set the number of deaths to zero
    x = &X[*nvar * p];
    BDH = 0.0;
  }

   for (p = 0; p < *np; p++) {

   x = &X[*nvar * p];

	for(j = 0; j < *seas_dim/2; j++) {
	xa[j+1] = 0.0 + (1.0)/((float)*seas_dim/2-1.0)*j ;
	ya[j+1] = (float)(&LOGBETA)[j];
	za[j+1] = (float)(&LOGBETA)[j+(*seas_dim/2)];
	}

//  slope = ya[*seas_dim];
//  ya[*seas_dim] = ya[1]; // to make the end and 1st equal !!
  cut =2021-(3.36-exp(ETA))/33.6;
  slope = 0;
  spline(xa, ya, *seas_dim/2, slope, slope, y2a);
  spline(xa, za, *seas_dim/2, slope, slope, z2a);
	//if(p<0)
	//for(j=0;j< *seas_dim;j++)
	//Rprintf("%10.5f %10.5f %10.5f\n",xa[j+1],ya[j+1],y2a[j+1]);

  for (step = 0, t = *t1; step < *nstep; step++, t += dt) {
  day=(int)floor(366*(t - *t_start));
 // splint(xa, ya, y2a,*seas_dim,t-floor(t), &beta);
  if(t<cut)
  splint(xa, ya, y2a,*seas_dim/2,(t-*t_start)/(cut-*t_start), &beta);
  if(t>cut)
  splint(xa, za, z2a,*seas_dim/2,(t-cut)/(*t_end-cut), &beta);
	if(p<0)
	Rprintf("%10.5f %10.5f\n",t,beta);


//void basic_seir_b (double *Sh, double *Eh, double *Ih, double *Th, double *Dh, double *Rh,  double *C,
//                   double beta, double sigmah, double p1, double theta, double gammah, double psi, double kappah,
//                   double dt)
r=1;
	beta=exp(beta);
	if(t>cut){r=exp(ALPHA);}
	if( (t>cut) & (t<cut+1.0/365)){BSH+=round((BRH+BRH1+BRH2)*expit(PSI));
		BRH=round(BRH*(1-expit(PSI)));
		BRH1=round(BRH1*(1-expit(PSI)));
		BRH2=round(BRH2*(1-expit(PSI)));
	}
basic_seir_b(&BSH,  &BEH, &BIH, &BTH, &BDH, &BRH, &BRH1, &BRH2, &TH1, &TH2, &BCA, &CC,
   beta, exp(SIGMAH), expit(THETA), exp(GAMMAH), r, exp(KAPPAH), 0.5, dt, expit(P0)*vac[day], exp(ETA)); 
    }
  }

  PutRNGstate();		// finished with R's random number generator

}


// simulate Euler-multinomial transitions
void euler_multinomial (int ntrans, double n, double *rate, double *trans, double dt) {
  double p = 0.0;
  int k;
  for (k = 0; k < ntrans; k++) p += rate[k]; // total event rate
  n = rbinom(n,1-exp(-p*dt));	// total number of events
  ntrans -= 1;
  for (k = 0; k < ntrans; k++) {
    trans[k] = ((n > 0) && (p > 0)) ? rbinom(n,rate[k]/p) : 0;
    n -= trans[k];
    p -= rate[k];
  }
  trans[ntrans] = n;
}


double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

double logit (double x) {
  return log(x/(1-x));
}

#define SC   (x[index[0]])
#define GC   (x[index[1]])
#define TAU  (x[index[2]])

void negbin_dmeasure (int *n, int *index, double *X, double *y1, double *y2, double *f) {
  int p, nv = n[0], np = n[1];
  double *x,  size, prob, f1, f2, tau, tol = 1.0e-18;
  for (p = 0; p < np; p++) {
    x = &X[nv*p];
    tau = exp(TAU);
    size = 1.0/tau;
    prob = 1.0/(1.0+SC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
     if (ISNA(*y1))
        f1 = 1.0;
        else
      f1 = dnbinom(*y1,size,prob,0)+tol;
    } else {
      f1 = tol;
    }
   prob = 1.0/(1.0+GC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
     if (ISNA(*y2))
        f2 = 1.0;
         else
      f2 = dnbinom(*y2,size,prob,0)+tol;
    } else {
      f2 = tol;
    }
        f[p]=f2;
  }
}


void negbin_rmeasure (int *n, int *index, double *X, double *cases) {
  int p, nv = n[0], np = n[1];
  double *x, tau, size, prob, tol = 1.0e-18;
  GetRNGstate();
  for (p = 0; p < np; p++) {
    x = &X[nv*p];
    tau = exp(TAU);
    size = 1.0/tau;
    prob = 1.0/(1.0+SC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
      cases[2*p] = rnbinom(size,prob+tol);
    } else {
      cases[2*p] = R_NaReal;
    }
        prob = 1.0/(1.0+GC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
      cases[2*p+1] = rnbinom(size,prob+tol);
    } else {
      cases[2*p+1] = R_NaReal;
    }
  }
  PutRNGstate();
}

#undef SC
#undef GC
#undef TAU




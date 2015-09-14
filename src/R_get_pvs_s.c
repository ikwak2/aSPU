#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"

void get_pvs_s(double *Zs, double *Ts,
	       double *npow, double *rnms, int n_pow, int n_Zs,
	       int n_perm, int Ps, double *pvs);

void R_get_pvs_s(double *Zs, double *Ts,
		 double *npow, double *rnms, int *n_pow, int *n_Zs,
		 int *n_perm, int *Ps, double *pvs) {
  get_pvs_s(Zs, Ts, npow, rnms, *n_pow, *n_Zs, *n_perm, *Ps, pvs);
}

void get_pvs_s(double *Zs, double *Ts,
	       double *npow, double *rnms, int n_pow, int n_Zs,
	       int n_perm, int Ps, double *pvs)
{
  int i, j, b, rr, cc, k;
  double *pPerm0, *T0s, *U0, ss, *P0s, *minP0s, minpPerm0, minp;

  pPerm0 = (double *) R_alloc ( n_pow, sizeof(double) ) ;
  T0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  U0 = (double *) R_alloc ( n_Zs, sizeof(double) ) ;
  P0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  minP0s = (double *) R_alloc ( n_perm, sizeof(double) ) ;


  for(i = 0 ; i < n_pow; i++ ) {
    pvs[i] = i;
  }
  pvs[n_pow] = minp;
}



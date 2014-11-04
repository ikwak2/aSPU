#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"

void get_pvs3(double *XUs, double T,
	      double *r,
	      int nr_XUs, int nc_XUs, int n_perm, int n_r, double *pv);

void R_get_pvs3(double *XUs, double *T,
		double *r,
		int *nr_XUs, int *nc_XUs, int *n_perm, int *n_r, double *pv) {
  get_pvs3(XUs, *T, r, *nr_XUs, *nc_XUs, *n_perm, *n_r, pv);
}

void get_pvs3(double *XUs, double T,
	      double *r,
	      int nr_XUs, int nc_XUs, int n_perm, int n_r, double *pv)
{
  int i, b, rr, cc;
  double *T0s, *U0, ss;
  T0s = (double *) R_alloc ( n_perm, sizeof(double) ) ;
  U0 = (double *) R_alloc ( nc_XUs, sizeof(double) ) ;


  for( i = 0 ; i < n_perm ; i++ ) {
    // r <- sample(r, length(r))
    double_permute(r, n_r);
    // U0 <- as.vector(t(XUs) %*% r0)
    for(cc = 0 ; cc < nc_XUs ; cc++) {
      ss = 0;
      for(rr = 0 ; rr < nr_XUs ; rr++) {
	ss += XUs[ cc*nr_XUs + rr ] * r[rr] ;
      }
      U0[cc] = ss;
    }

    //    for (j in 1:length(pow))
    //      if (pow[j] < Inf)
    //        T0s[b, j] = sum(U0^pow[j])
    //      else T0s[b, j] = max(abs(U0))
    ss = 0;
    for( b = 0 ; b < nc_XUs ; b ++) {
      ss += U0[b];
    }
    T0s[i] = ss;
  }


  //  pPerm0[j] = sum(abs(Ts[j]) <= abs(T0s[, j]))/n.perm
  ss = 0;
  for( i = 0 ; i < n_perm ; i++) {
    if ( abss(T) <= abss(T0s[ i ]) )
      ss += 1;
  }
  pv[0] = T;
  pv[1] = ss / n_perm;
}


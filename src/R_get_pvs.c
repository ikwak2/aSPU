#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "R_get_pvs.h"

void get_pvs(double *XUs, double *Ts,
	     double *npow, double *r, int n_pow,
	     int nr_XUs, int nc_XUs, int n_perm, int n_r, double *pvs);

void R_get_pvs(double *XUs, double *Ts,
	       double *npow, double *r, int *n_pow,
	       int *nr_XUs, int *nc_XUs, int *n_perm, int *n_r, double *pvs) {
  get_pvs(XUs, Ts, npow, r, *n_pow, *nr_XUs, *nc_XUs, *n_perm, *n_r, pvs);
}

void get_pvs(double *XUs, double *Ts,
	     double *npow, double *r, int n_pow,
	     int nr_XUs, int nc_XUs, int n_perm, int n_r, double *pvs)
{
  int i, j, b, rr, cc, one, k;
  double *pPerm0, *T0s, *r0, *U0, ss, *P0s, *minP0s, minpPerm0, minp;
  int *bb;

  pPerm0 = (double *) R_alloc ( n_pow, sizeof(double) ) ;
  T0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  U0 = (double *) R_alloc ( nc_XUs, sizeof(double) ) ;
  bb =  (int *) R_alloc ( n_r, sizeof(int) ) ;
  P0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  minP0s = (double *) R_alloc ( n_perm, sizeof(double) ) ;


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
    for( j = 0 ; j < (n_pow) ; j++ ) {

      ss = 0;
      if(npow[j] != 0) {
	for( b = 0 ; b < nc_XUs ; b ++) {
	  ss += pow( U0[b], npow[j] );
	}
      } else {
	ss = 0;
	for( b = 0 ; b < nc_XUs ; b ++) {
	  if( ss <  abss(U0[b]) )
	    ss = abss(U0[b]);
	}
      }
      T0s[i*n_pow + j] = ss;
    }
  }


  //  pPerm0[j] = sum(abs(Ts[j]) <= abs(T0s[, j]))/n.perm
  for( j = 0 ; j < n_pow ; j++ ) {
    ss = 0;
    for( i = 0 ; i < n_perm ; i++) {
      if ( abss(Ts[j]) < abss(T0s[ i*n_pow + j ]) )
	ss += 1;
    }
    pPerm0[j] = ss / n_perm;
  }

  //	P0s = (n.perm-rank(abs(T0s)))/(n.perm - 1)

  for( j = 0 ; j < n_pow ; j++ ) {
    for( i = 0 ; i < n_perm ; i++) {
      ss = 0;
      for( k = 0 ; k < n_perm ; k++) {
	if(k != i) {
	  if ( abss(T0s[ i*n_pow + j]) <= abss(T0s[ k*n_pow + j]) )
	    ss++;
	}
      }
      P0s[ i*n_pow + j] = (ss + 1) / n_perm ;
    }
  }
  //  if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
  for(i = 0 ; i < n_perm ; i++ ) {
    minP0s[i] = P0s[i*n_pow + 0]  ;
  }
  for(j = 1 ; j < n_pow ; j ++) {
    for(i = 0 ; i < n_perm ; i++ ) {
      if(minP0s[i] > P0s[i*n_pow + j])
	minP0s[i] = P0s[i*n_pow + j];
    }
  }

  // minp = sum( min(pPerm0) > minP0s ) / B
  ss = 1;
  for(i=0 ; i < n_pow; i++) {
    if( ss > pPerm0[i] )
      ss = pPerm0[i];
  }
  minpPerm0 = ss;

  ss = 0;
  for(i = 0; i < n_perm; i++ ) {
    if(minpPerm0 > minP0s[i])
      ss += 1;
  }
  minp = ss / n_perm;

  // pvs <- c(pPerm0, minP)
  for(i = 0 ; i < n_pow; i++ ) {
    pvs[i] = pPerm0[i];
  }
  pvs[n_pow] = minp;
}

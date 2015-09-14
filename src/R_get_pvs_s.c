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
  T0s = (double *) R_alloc ( n_pow, sizeof(double) ) ;
  U0 = (double *) R_alloc ( n_Zs, sizeof(double) ) ;
  P0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  minP0s = (double *) R_alloc ( n_perm, sizeof(double) ) ;


  for( i = 0 ; i < n_perm ; i++ ) {

    for(rr = 0 ; rr < n_Zs ; rr++) {
      if( Ps == 1 )
	U0[rr] = fabs(rnms[ n_Zs * i + rr ] );
      else
	U0[rr] = rnms[ n_Zs * i + rr ] ;
    }

    for( j = 0 ; j < (n_pow) ; j++ ) {

      ss = 0;
      if(npow[j] != 0) {
	for( b = 0 ; b < n_Zs ; b ++) {
	  ss += pow( U0[b], npow[j] );
	}
      } else {
	for( b = 0 ; b < n_Zs ; b ++) {
	  if( ss <  fabs(U0[b]) )
	    ss = fabs(U0[b]);
	}
      }
      T0s[i*n_pow + j] = ss;
    }
  }
  //  printf("\n");


  //  for( j = 0 ; j < n_pow ; j++) {
  //    printf("%f ", Ts[j]);
  //  }
  //  printf("\n");


  for( j = 0 ; j < n_pow ; j++ ) {
    ss = 0;
    for( i = 0 ; i < n_perm ; i++) {
      if ( fabs(Ts[j]) <= fabs(T0s[   j ]) )
	ss += 1;

      //      printf("%f  ", T0s[i*n_pow + j ] );
    }
    pPerm0[j] = ss / n_perm;
    //    printf("%f  ", pPerm0[j]);
    //    printf("\n");
  }
  //  printf("\n");

  //	P0s = (n.perm-rank(abs(T0s)))/(n.perm - 1)

  for( j = 0 ; j < n_pow ; j++ ) {
    for( i = 0 ; i < n_perm ; i++) {
      ss = 0;
      for( k = 0 ; k < n_perm ; k++) {
	if(k != i) {
	  if ( fabs(T0s[  j]) <= fabs(T0s[  j]) )
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

  // minp = (sum( min(pPerm0) >= minP0s )+1) / (B+1)
  ss = 1;
  for(i=0 ; i < n_pow; i++) {
    if( ss > pPerm0[i] )
      ss = pPerm0[i];
  }
  minpPerm0 = ss;

  ss = 0;
  for(i = 0; i < n_perm; i++ ) {
    if(minpPerm0 >= minP0s[i])
      ss += 1;
  }
  minp = (ss+1) / (n_perm+1);

  // pvs <- c(pPerm0, minP)
  for(i = 0 ; i < n_pow; i++ ) {
    pvs[i] = pPerm0[i];
  }
  pvs[n_pow] = minp;
}

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}


void get_pvs_s(double *Zs, double *CovSsqrt, double *Ts,
	       double *npow, int *seeds, int n_pow, int n_Zs,
	       int n_perm, int Ps, double *pvs);

void R_get_pvs_s(double *Zs, double *CovSsqrt, double *Ts,
		 double *npow, int *seeds, int *n_pow, int *n_Zs,
		 int *n_perm, int *Ps, double *pvs) {
  get_pvs_s(Zs, CovSsqrt, Ts, npow, seeds, *n_pow, *n_Zs, *n_perm, *Ps, pvs);
}

void get_pvs_s(double *Zs, double *CovSsqrt, double *Ts,
	       double *npow, int *seeds, int n_pow, int n_Zs,
	       int n_perm, int Ps, double *pvs)
{
  int i, j, b, rr, cc, k;
  double *pPerm0, *T0s, *U0, *nR, ss, *P0s, *minP0s, minpPerm0, minp;
  //  int *bb;

  nR = (double *) R_alloc ( n_Zs, sizeof(double) ) ;
  pPerm0 = (double *) R_alloc ( n_pow, sizeof(double) ) ;
  T0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  U0 = (double *) R_alloc ( n_Zs, sizeof(double) ) ;
  //  bb =  (int *) R_alloc ( n_r, sizeof(int) ) ;
  P0s = (double *) R_alloc ( n_perm * n_pow, sizeof(double) ) ;
  minP0s = (double *) R_alloc ( n_perm, sizeof(double) ) ;


  for( i = 0 ; i < n_perm ; i++ ) {
    //        U00<-rnorm(k, 0, 1)
    //    printf("Random smp :");

    for(j = 0 ; j < n_Zs ; j++ ) {
      nR[j] = randn(0,1);
      //      printf("%f  ", nR[j]);
    }
    //    printf("\n");

    // 	U0<-CovSsqrt %*% U00

    //    printf("Random smp2 :");

    for(cc = 0 ; cc < n_Zs ; cc++) {
      ss = 0;
      for(rr = 0 ; rr < n_Zs ; rr++) {
	ss += CovSsqrt[ rr*n_Zs + cc ] * nR[rr] ;

	//	printf("%f ", CovSsqrt[rr*n_Zs + cc]);
      }

      if( Ps == 1 )
	U0[cc] = abss(ss);
      else
	U0[cc] = ss;
      //      printf("%f  ", U0[cc]);
    }
    //printf("\n");
//printf("\n");

    //    for(j = 0 ; j < n_Zs ; j++ ) {
    // printf("%f  ", U0[j]);
    //}
    //printf("\n");


    for( j = 0 ; j < (n_pow) ; j++ ) {

      ss = 0;
      if(npow[j] != 0) {
	for( b = 0 ; b < n_Zs ; b ++) {
	  ss += pow( U0[b], npow[j] );
	}
      } else {
	for( b = 0 ; b < n_Zs ; b ++) {
	  if( ss <  abss(U0[b]) )
	    ss = abss(U0[b]);
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
      if ( abss(Ts[j]) <= abss(T0s[ i*n_pow + j ]) )
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



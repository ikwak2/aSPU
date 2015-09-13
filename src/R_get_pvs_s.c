#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"

# include <complex.h>
# include <time.h>


/******************************************************************************/

float r4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a unit pseudorandom R4.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r4_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R4_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  const int i4_huge = 2147483647;
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}



/******************************************************************************/

float r4_normal_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_NORMAL_01 returns a unit pseudonormal R4.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    The Box-Muller method is used, which is efficient, but
    generates two values at a time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2013

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R4_NORMAL_01, a normally distributed random value.
*/
{
  float r1;
  float r2;
  const double r4_pi = 3.141592653589793;
  float x;

  r1 = r4_uniform_01 ( seed );
  r2 = r4_uniform_01 ( seed );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r4_pi * r2 );

  return x;
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
      nR[j] = r4_normal_01(&seeds[i*n_Zs + j]);
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



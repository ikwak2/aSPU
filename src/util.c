#include <R.h>
#include "util.h"

/**********************************************************************
 *
 * random_int
 *
 * Generates a random int integer between "low" and "high", inclusive.
 *
 *  Input:
 *
 *    low
 *
 *    high
 *
 **********************************************************************/
int random_int(int low, int high)
{
  return((int)(unif_rand()*(double)(high - low + 1)) + low);
}

/**********************************************************************
 *
 * double_permute
 *
 *   This function randomly permutes a vector of doubles
 *
 * Input:
 *
 *   array = vector of doubles; on output, it contains a random
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/
void double_permute(double *array, int len)
{
  int i, which;
  double tmp;

  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}


double abss(double a) {
  if (a < 0) {
    return(-a);
  } else {
    return(a);
  }
}



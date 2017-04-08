#include <armadillo>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0 (arma::mat tXUs, arma::mat yresid, arma::mat powV, int nperm) {
  //const int n = X1.n_rows;
  const int npow = powV.n_rows;
  // containers
  arma::mat T0s1(nperm,npow);
  T0s1.fill(0);

  for (int i = 0; i < nperm; i++) {
    arma::mat r10 = shuffle(yresid);
    arma::mat U01 = tXUs * r10;

    for (int j = 0; j < npow; j++) {
      if (powV(j,0) == 0) {
        arma::mat tmpU01 = abs(U01);
        T0s1(i,j) = tmpU01.max();
      } else {
        T0s1(i,j) = accu(pow(U01,powV(j,0)));
      }
    }
  }
  Rcpp::List res;
  res["T0"] =T0s1;

  return(res);
}

// [[Rcpp::export()]]
Rcpp::List calcT0sim (arma::mat CvSqrt, arma::mat powV, int nperm) {
  //const int n = X1.n_rows;
  const int npow = powV.n_rows;
  const int k = CvSqrt.n_rows;
  // containers
  arma::mat T0s1(nperm,npow);
  T0s1.fill(0);

  for (int i = 0; i < nperm; i++) {
    arma::mat U00 = arma::randn(k,1);
    arma::mat U0 = CvSqrt * U00;

    for (int j = 0; j < npow; j++) {
      if (powV(j,0) == 0) {
        arma::mat tmpU01 = abs(U0);
        T0s1(i,j) = tmpU01.max();
      } else {
        T0s1(i,j) = accu(pow(U0,powV(j,0)));
      }
    }
  }
  Rcpp::List res;
  res["T0"] =T0s1;

  return(res);
}


// [[Rcpp::export()]]
Rcpp::List calcT0simM2 (arma::mat A, arma::mat B, arma::mat powV, arma::mat pow2V, int nperm, int Pval, arma::mat Ts) {
  //const int n = X1.n_rows;
  const int npow = powV.n_rows;
  const int npow2 = pow2V.n_rows;
  const int nphe = A.n_rows;
  const int nsnp = B.n_rows;
  // containers
  int np = npow * npow2;
  arma::mat T0s1(nperm, np);
  T0s1.fill(0);

  for (int i = 0; i < nperm; i++) {
    arma::mat Z00 = arma::randn(nphe,nsnp);
    arma::mat Z0 = arma::trans(A * Z00 * B);

    if(Pval == 1) {
      Z0 = abs(Z0);
    }

    for (int j = 0; j < npow; j++) {
      //      printf("j = %d  ", j);


      for( int j2 = 0; j2 < npow2; j2++) {
        //      printf("j2 = %d  ", j2);

        if (powV(j,0) == 0) {
          arma::mat tmpZ01 = abs(Z0);
          arma::mat Tmp(1, nphe);
          Tmp.fill(0);
          for(int t1 = 0; t1 < nphe; t1 ++) {
            Tmp(0,t1) = max(tmpZ01.col(t1)) ;
          }
          T0s1(i, j2 + j*npow2 ) = accu(pow(Tmp,pow2V(j2,0)));
        } else {
          arma::mat Tmp(1, nphe);
          Tmp.fill(0);
          for(int t1 = 0; t1 < nphe; t1 ++) {
            //            Tmp(0,t1) = accu(pow(Z0.col(t1), powV(j1,0)) ) ;
            arma::mat Tmp2 = pow(Z0.col(t1), powV(j,0));

            Tmp(0,t1) = accu(Tmp2) ;
          }

          T0s1(i, j2 + j*npow2 ) = accu(pow(Tmp,pow2V(j2,0)));


        }
      }
    }
  }

  arma::mat Pperm(1, npow*npow2);
  Pperm.fill(0);

  for(int i = 0; i < npow*npow2; i++ ) {
    for(int j = 0; j < nperm ; j++ ) {
      if( std::abs(T0s1(j,i)) > std::abs(Ts(0,i)) ) {
        Pperm(0,i) = Pperm(0,i) + 1;
      }
    }
  }

  arma::mat Pvs(1, npow*npow2);
  for(int i = 0; i < npow*npow2; i++ ) {
    Pvs(0,i) = Pperm(0,i)/nperm;
  }

  double minPperm = 1;
  minPperm = Pvs.min();

  arma::mat PPvs(nperm, np);
  PPvs.fill(0);

  for(int j = 0; j < nperm; j++) {
    for(int i = 0; i < np; i++) {
      for(int k = 0; k < nperm; k++) {
        if(k != j) {
          if( std::abs(T0s1(k,i)) > std::abs(T0s1(j,i)) ) {
            PPvs(j, i) = PPvs(j,i) + 1;
          }
        }
      }
    }
  }

  arma::mat mPs(1,nperm);
  mPs.fill(nperm);

  for(int j = 0; j < nperm ; j++) {
    for(int i = 0; i < np ; i++) {
      if( mPs(0,j) > PPvs(j,i)) {
        mPs(0,j) = PPvs(j,i);
      }
    }
  }

  int minp = 0;
  for(int j=0; j < nperm; j++) {
    if( minPperm > mPs(0,j)/nperm ) {
      minp = minp + 1;
    }
  }

  Rcpp::List res;
  res["Pvs"] = Pvs;
  res["minp"] = minp;

  return(res);
}



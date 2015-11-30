# include <RcppArmadillo.h>
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


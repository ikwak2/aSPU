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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List aSPUsPathEngine(Rcpp::List CH, Rcpp::List CHcovSq, arma::vec pow1, arma::vec pow2, int nGenes, int n_perm, int k, int Ps, arma::vec nSNPs0, arma::vec StdTs) {
  const int n_ch = CH.size();
  const int n_pow1 = pow1.size();
  const int n_pow2 = pow2.size();
  
  arma::mat T0(n_perm, n_pow1*nGenes);

  // iterate for n_perm
  for(int b=0; b < n_perm; b++) {
    Rcpp::NumericVector U00 = rnorm(k,0,1);
    arma::vec u0 = as<arma::vec>(U00);
    arma::vec U0(1);

    // iterate for chromosome
    for(int b2 = 0; b2 < n_ch ; b2++) {
      NumericVector CH1 = CH[b2] ;
      uvec idxu = as<uvec>(CH1) - 1 ;
      
      arma::vec TT = as<arma::mat>(CHcovSq[b2])*u0.elem(idxu) ;
      U0 = join_cols(U0,TT) ;
    }
    arma::vec UU2 = U0.subvec(1,U0.size()-1);

    // Use absolute values when inputs are p-values
    if( Ps == 1 ) {
      UU2 = abs(UU2);
    }

    // iterate for pow1
    for( int j = 0 ; j < pow1.size() ; j++) {
      int SNPstart = 0;

      // itreate for genes
      for(int iGene = 0 ; iGene < nGenes ; iGene++ ) {
        if( iGene != 0) {
          SNPstart = sum( nSNPs0.subvec(0,iGene-1) ) ; 
        }
        int idx1 = SNPstart;
        int idx2 = SNPstart+nSNPs0(iGene)-1;

        if( pow1[j] > 0) {
          arma::vec tmp1 = pow(UU2.subvec(idx1, idx2),pow1[j]);
          double tmp2 = sum(tmp1);

          // Calculate first level statistics
          if( tmp2 > 0 ) {
            T0(b, j*nGenes+iGene) = pow(std::abs(tmp2)/nSNPs0[iGene] , 1/pow1[j]);
          } else {
            T0(b, j*nGenes+iGene) =  -pow(std::abs(tmp2)/nSNPs0[iGene] , 1/pow1[j]);
          }
          
        } else {
          arma::vec T0tp = abs(UU2.subvec(idx1, idx2));
          T0(b, j*nGenes+iGene) = max(abs(UU2.subvec(idx1, idx2)));
        }
        
      }
    }
  }

  // Calculate 2nd level Statistics 
  
  arma::vec Ts2(n_pow1*n_pow2);
  arma::mat T0s2(n_perm, n_pow1*n_pow2);
  
  for(int j2=0 ; j2 < n_pow2; j2++) {
    for(int j=0 ; j < n_pow1; j++) {
      if(pow2[j2] > 0) {

        // Calculate 2nd level statistics
        arma::vec tmp3 = pow(StdTs.subvec(j*nGenes,(j+1)*nGenes-1), pow2[j2]);
        double tmp4 = sum(tmp3);
        Ts2(j2*n_pow1 +j) = tmp4;
        
        for(int b = 0 ; b < n_perm ; b++ ) {
          arma::rowvec T0tmp = T0.row(b);
          arma::rowvec tmp3 = pow(T0tmp.cols(j*nGenes,(j+1)*nGenes-1), pow2[j2]);

          float tmp4 = sum(tmp3);
          T0s2(b, j2*n_pow1 +j) = tmp4;
          
        }
      } else {
        Ts2(j2*n_pow1 + j) = max( arma::abs(StdTs.subvec(j*nGenes,(j+1)*nGenes-1)) ) ;
        
        for(int b = 0 ; b < n_perm ; b++ ) {
          arma::rowvec T0tmp = T0.row(b);
          T0s2(b, j2*n_pow1 + j) = max( arma::abs(T0tmp.cols(j*nGenes,(j+1)*nGenes-1)) ) ;
        }
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("T0") = T0,
                            Rcpp::Named("Ts2") = Ts2,
                            Rcpp::Named("T0s2") = T0s2);
}


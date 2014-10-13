#' SUM test (residual permutation, permutation part coded in C)
#'
#' It gives the p-values of the SUM test based on the permutation of residuals.
#'
#' @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#'     or It can be any quantitative traits.
#'
#' @param X genotype data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele.
#'
#' @param cov covariates
#'
#' @param model Use "gaussian" for quantitative trait (Default)
#'    , and Use "binomial" for binary trait.
#'
#' @param n.perm number of permutation
#'
#' @export
#' @return SUM Test Statistic and p-value for SUM test.
#'
#' @examples
#'
#' data(exdat)
#' out <- SUMtestC(exdat$Y, exdat$X, cov = NULL, model = "binomial", n.perm = 1000)
#' out
#'
#' @seealso \code{\link{aSPU}}

SUMtest <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), n.perm=1000){

    model = match.arg(model)

    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)

#### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        XUs<-Xg <- X
        r<-Y-mean(Y)
        U<-as.vector(t(Xg) %*% r)
   } else {
       tdat1<-data.frame(trait=Y, cov)
       fit1<-glm(trait~.,family=model,data=tdat1)
       pis<-fitted.values(fit1)
       XUs<-matrix(0, nrow=n, ncol=k)
       Xmus = X
       for(i in 1:k){
           tdat2<-data.frame(X1=X[,i], cov)
           fit2<-glm(X1~.,data=tdat2)
           Xmus[,i]<-fitted.values(fit2)
           XUs[, i]<-(X[,i] - Xmus[,i])
       }
       r<-Y - pis
       U<-t(XUs) %*% r
   }

    ##observed statistics
    T = sum(U)

    ## residual permutation
    nr_XUs = nrow(XUs)
    nc_XUs = ncol(XUs)
    n_perm = n.perm

    output <- .C("R_get_pv",
                 as.double(XUs),
                 as.double(T),
                 as.double(r),
                 as.integer(nr_XUs),
                 as.integer(nc_XUs),
                 as.integer(n_perm),
                 pv = as.double( 0 ),
                 PACKAGE="aSPU")

    pv <- output$pv

    list(T = T, pv = output$pv)
}

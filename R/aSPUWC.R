#' Variance-weighted Sum of powered score (SPUw) test; using permutations to get the p-values (NO adjustment for covariates yet, C version).
#'
#' It gives the p-values of the SPUw test and aSPUw test based on based on the permutation of residuals.  (NO adjustment for covariates yet, C version)
#'
#' @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#'     or It can be any quantitative traits.
#'
#' @param X genotype data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele.
#'
#' @param pow power used in SPU test.
#'
#' @param n.perm number of permutation
#'
#' @param userank similar code with the original version if TRUE. by definition if FALSE
#'
#' @export
#' @return Test Statistics and p-values for SPU tests and aSPU test.
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPUWC(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
#' out
#'
#' @seealso \code{\link{aSPU}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPUWC <- function(Y, X, pow=c(1:8, Inf), n.perm=1000, userank = T){

    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)

    Xg <- X

############SSU#####################
    Xbar<-apply(Xg, 2, mean)
    #   Xgb<-Xg
    #   for(i in 1:nrow(Xg))
    #      Xgb[i,]<-Xg[i,]-Xbar
    # faster than doing above:
    subtract<-function(x, y) { x - y }
    Xgb=t(apply(Xg, 1, subtract, Xbar))

    r=Y-mean(Y)

#########Score vector:
    U<-as.vector( t(Xg) %*% r)
    #Var(U):
    v0=mean(Y)*(1-mean(Y))
    p=length(U)
    diagCovS=rep(0, p)
    for(i in 1:p)
        diagCovS[i] = v0 * sum(Xgb[,i]^2)
    diagCovS<-ifelse(diagCovS>1e-20, diagCovS, 1e-20)
    diagSDs<-sqrt(diagCovS)

   # test stat's:
    Ts<-rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
            Ts[j] = sum((U/diagSDs)^pow[j]) else Ts[j] = max(abs(U/diagSDs))
    #     VarTs[j] = var(Upow)
    }

   # permutations:
    XUs = Xg
    n_pow = length(pow)
    nr_XUs = nrow(XUs)
    nc_XUs = ncol(XUs)
    n_perm = n.perm
    n_r = length(r)
    if(userank)
        r <- jitter(r, amount = 0.0001)

    output <- .C("R_get_pvs2",
                 as.double(XUs),
                 as.double(Ts),
                 as.double(npow),
                 as.double(r),
                 as.double(diagSDs),
                 as.integer(n_pow),
                 as.integer(nr_XUs),
                 as.integer(nc_XUs),
                 as.integer(n_perm),
                 as.integer(n_r),
                 pvs = as.double( rep(0,n_pow + 1) ),
                 PACKAGE="aSPU")

    pvs <- output$pvs

    Ts <- c(Ts, min(pvs[1:n_pow]) )
    names(Ts) <- c(paste("SPUw", pow, sep=""), "aSPUw")
    names(pvs) = names(Ts)

    list(Ts = Ts, pvs = pvs)
}

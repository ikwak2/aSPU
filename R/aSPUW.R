#' Variance-weighted Sum of powered score (SPUw) test; using permutations to get the p-values (NO adjustment for covariates yet).
#'
#' It gives the p-values of the SPUw test and aSPUw test based on based on the permutation of residuals.  (NO adjustment for covariates yet)
#'
#' @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#'     or It can be any quantitative traits. Vector with length n (number of observations)
#'
#' @param X genotype data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele. Matrix with dimension n by g (n : number of observation, p : number of genotype data)
#'
#' @param pow power used in SPU test. Vector of g number of power.
#'
#' @param n.perm number of permutation
#'
#' @export
#' @return Test Statistics and p-values for SPU tests and aSPU test.
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPUW(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
#' out
#'
#' @seealso \code{\link{aSPU}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPUW <- function(Y, X, pow=c(1:8, Inf), n.perm=1000){

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
    pPerm0<-rep(0, length(pow))
    T0s = numeric(n.perm)
    s <- sample(1:10^5,1)
    Y0 <- numeric(n)
    for (j in 1:length(pow)){
        set.seed(s) # to ensure the same samples are drawn for each pow
        for(b in 1:n.perm){
            r0 <- sample(r, length(r))
            U0<-as.vector( t(Xg) %*% r0)
     # test stat's:
            if (pow[j] < Inf) {T0s[b] = round( sum((U0/diagSDs)^pow[j]), digits=8) }
            if (pow[j] == Inf) {T0s[b] = round( max(abs(U0/diagSDs)), digits=8)}
        }

        pPerm0[j] = round((sum(abs(Ts[j]) <= abs(T0s)))/(n.perm), digits=8)
        P0s = (n.perm-rank(abs(T0s))+1)/(n.perm)
        if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
                                        # cat("g=",g,"\n")
    }
    cat("P0s caculated","\n")
    Paspu<-(sum(minp0<=min(pPerm0))+1)/(n.perm+1)
    pvs <- round(c(pPerm0, Paspu), digits = 3)

    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPUw", pow, sep=""), "aSPUw")
    names(pvs) = names(Ts)

    list(Ts = Ts, pvs = pvs)
}


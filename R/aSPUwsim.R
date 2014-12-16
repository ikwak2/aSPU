#' Variance-weighted adaptive Sum of powered score (SPUw) test; using permutations to get the p-values (NO adjustment for covariates yet).
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
#' @param cov covariates. Matrix with dimension n by k (n :number of observation, k : number of covariates)
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


aSPUwsim <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), pow=c(1:8, Inf), n.perm=1000){

    model <- match.arg(model)

    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)

    if (is.null(cov)){
        ## NO nuisance parameters:
        Xg <- X
        Xbar<-apply(Xg, 2, mean)
        subtract<-function(x, y) { x - y }
        Xgb=t(apply(Xg, 1, subtract, Xbar))

        r=Y-mean(Y)

        U<-as.vector( t(Xg) %*% r)

        ##cov of the score stats:
        CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)

    } else {
        ## with nuisance parameters:
        tdat1<-data.frame(trait=Y, cov)
        fit1<-glm(trait~.,family=model,data=tdat1)
        pis<-fitted.values(fit1)
        Us<-matrix(0, nrow=n, ncol=k)
        for(i in 1:k){
            tdat2<-data.frame(X1=X[,i], cov)
            fit2<-glm(X1~.,data=tdat2)
            X1mus<-fitted.values(fit2)
            r <- Y - pis
            Us[, i]<-(Y - pis)*(X[,i] - X1mus)
        }
        U<-apply(Us, 2, sum)
        CovS<-matrix(0, nrow=k, ncol=k)
        for(i in 1:n)
            CovS<-CovS + Us[i,] %*% t(Us[i,])
    }



    Vs<-diag(CovS)
    diagSDs<-ifelse(Vs>1e-20, sqrt(Vs), 1e-10)

#    svd.CovS<-svd(CovS)
#    CovSsqrt<-svd.CovS$u %*% diag(sqrt(svd.CovS$d))


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




#calculate a permuttaion p-value matrix based on a permuted stat matrix:
PermPvs<-function(T0s){
    B=nrow(T0s); n=ncol(T0s)
    P0s<-matrix(1, nrow=B, ncol=n)
    for(j in 1:n)
        for(b in 1:B)
            P0s[b,j] = sum( abs(T0s[b,j]) < abs(T0s[-b,j]) )/(B-1)
    return(P0s)
}

# Main function for the SPU tests and the adaptive SPU (aSPU) test:

SPUwtest<-function(Y, X, Z=NULL, pow=c(2, 4, 16), B=100, NormalApprox=T,
                   PermRawP=F ){

    n<-length(Y)
    k<-ncol(X)
    k2<-ncol(Z)

    #######construction of the score vector and its cov matrix:
    if (is.null(Z)){
        ## NO nuisance parameters:
        Xg <- X
        Xbar<-apply(Xg, 2, mean)
        Xgb<-Xg
        for(i in 1:nrow(Xg))
            Xgb[i,]<-Xg[i,]-Xbar
        ##score vector:
        U<-t(Xg) %*% (Y-mean(Y))
        ##cov of the score stats:
        CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)

    } else {
        ## with nuisance parameters:
        tdat1<-data.frame(trait=Y, Z)
        fit1<-glm(trait~.,family="binomial",data=tdat1)
        pis<-fitted.values(fit1)
        Us<-matrix(0, nrow=n, ncol=k)
        for(i in 1:k){
            tdat2<-data.frame(X1=X[,i], Z)
            fit2<-glm(X1~.,data=tdat2)
            X1mus<-fitted.values(fit2)
            Us[, i]<-(Y - pis)*(X[,i] - X1mus)
        }
        U<-apply(Us, 2, sum)
        CovS<-matrix(0, nrow=k, ncol=k)
        for(i in 1:n)
            CovS<-CovS + Us[i,] %*% t(Us[i,])
    }

    Vs<-diag(CovS)
    SDs<-ifelse(Vs>1e-20, sqrt(Vs), 1e-10)
    Uw = U/SDs

    svd.CovS<-svd(CovS)
    CovSsqrt<-svd.CovS$u %*% diag(sqrt(svd.CovS$d))

    # test stat's:
    Ts<-VarTs<-rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
            Ts[j] = sum(Uw^pow[j]) else Ts[j] = max(abs(Uw))
    #     VarTs[j] = var(Upow)
    }

    # simulations: simulate the score vector from its null distr,
    T0s = matrix(0, nrow=B, ncol=length(pow))
    for(b in 1:B){
        U00<-rnorm(k, 0, 1)
        U0<-CovSsqrt %*% U00

        #    alternatively, by permutations:
        #     Y0 <- sample(Y, length(Y))
        #     U0<-as.vector( t(Xg) %*% (Y0-meanY))

        # test stat's:
        Uw0 = U0/SDs
        for(j in 1:length(pow))
            if (pow[j] < Inf)
                T0s[b, j] = sum(Uw0^pow[j]) else T0s[b, j] = max(abs(Uw0))

    }

    # p-values based on the indepdence assumption on SNPs (i.e. components of U):
    # need to correct the mean of the moment of a normal variate:
    # if Z ~U(mu, sigma), then E((Z-mu)^p) = sigma^p * (p-1)!!
    #   pIndep = 1 - pchisq(((Ts-meanTs)^2)/(length(U)*VarTs), df=1)

    # permutation-based p-values:
    pPerm = pPerm0 = rep(NA, length(pow));
    pvs = NULL;
    if (NormalApprox){
        for(j in 1:length(pow)) {
            mu = mean(T0s[,j]); sigma = sd(T0s[,j])
            pPerm[j] = 1 - pchisq(((Ts[j]-mu)/sigma)^2, df=1)
        }
        pvs<-c(pvs, pPerm)
    }
    if (PermRawP) {
        for(j in 1:length(pow)) {
            pPerm0[j] = sum( abs(Ts[j]) < abs(T0s[,j]))/B
        }
        P0s = PermPvs(T0s)
        minP0s = apply(P0s, 1, min)
        minP =  sum( min(pPerm0) > minP0s )/B
        pvs<-c(pvs, pPerm0, minP)
    }

    Ts <- c(Ts, min(pPerm0) )
                                        # return(c(pIndep, pPerm))
    return(pvs)
}

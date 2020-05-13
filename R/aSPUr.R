#' Robust Sum of powered score (SPU) tests and aSPU test for a quantitative trait
#'
#' The test is based on the Huber loss function and using the parametric
#' bootstrap for inference (i.e. bootstrapping residuals).
#'
#' @param Y a vector of quantitative traits (QTs).
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP (or a predictor). The value of each SNP is the # of the copies
#'     for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ).
#'
#' @param cov Covariates. A matrix with dimension n by k2 (n :number of observation, k2 : number of covariates).
#'
#' @param pow power used in SPUr test. A vector of the powers.
#'
#' @param B number of bootstraps.
#'
#' @param C Constant in huber loss function. C = 1.345 is chosen to maintain a high efficiency for a Normal error.
#'
#' @return p-values of the SPUr tests in the order of supplied pow values; finally, the p-value of the aSPUr test (that combines the SPUs tests with pow by taking their min P-value and adjust for multiple testing).
#'
#' @author Yiwei Zhang and Wei Pan
#'
#' @references
#' Peng Wei, Ying Cao, Yiwei Zhang, Zhiyuan Xu, Il-Youp Kwak, Eric Boerwinkle, Wei Pan (2016) On Robust Association Testing for Quantitative Traits and Rare Variants, G3, 6(12) 3941-3950.
#'
#' @examples
#'
#'
#' data(exdat)
#'
#' ## example analysis using aSPU test on exdat data.
#'
#' QT <- jitter(exdat$Y)
#'
#' out <- aSPUr(Y = QT, X = exdat$X, cov = NULL, B = 100)
#' out
#' ## This is a vector of p-values for SPUr and aSPUr tests.
#' ## SPU1 to SPUInf corresponds with the option pow=c(1:8, Inf)
#' ## They are p-values for corresponding SPUr tests.
#' ## The last element is p-value of aSPUr test.
#'
#' @seealso \code{\link{aSPU}}
#'
#' @export

aSPUr<-function(Y, X, cov=NULL, pow=c(1:8, Inf), B=1000,C=1.345){
    
    n<-length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k<-ncol(X)
    
    if (is.null(cov)){
        ## NO nuisance parameters:
        Xg <- X
        
        ##score vector:
        ##   U<-t(Xg) %*% (Y-mean(Y))
        yresids<-Y-mean(Y)
        ##sigma0=sqrt(sum(yresids^2)/(n-1))
        ##yfits<-rep(mean(Y), n)
################################### this is using huber's loss
        sigma.mad=mad(Y-mean(Y))
        
        ##yy=ifelse(abs(Y-mean(Y))<=C*sigma.mad, (Y-mean(Y))/sigma.mad, sign(Y-mean(Y))*C)
        ##above changed to the below 6/5/15 by wp:
        yy=ifelse(abs(Y-mean(Y))<=C*sigma.mad, (Y-mean(Y)), sign(Y-mean(Y))*C*sigma.mad)
        U<-t(Xg)%*%yy
        
        ##cov of the score stats:
                                        #   CovS<- (sigma0^2)*(t(Xgb) %*% Xgb)
    } else {
        ## with nuisance parameters:
        tdat1<-data.frame(trait=Y, cov)
        fit1<-glm(trait~., data=tdat1)
        yfits<-fitted.values(fit1)
        yresids<-fit1$residuals
        fit1res1<-summary(fit1)
        sigma0<-sqrt(fit1res1$dispersion)
        
        Us<-XUs<-matrix(0, nrow=n, ncol=k)
        Xmus = X
        for(i in 1:k){
            tdat2<-data.frame(X1=X[,i], cov)
            fit2<-glm(X1~., data=tdat2)
            Xmus[,i]<-fitted.values(fit2)
            XUs[, i]<-(X[,i] - Xmus[,i])
        }
################################### this is using huber's loss
        sigma.mad=mad(Y-yfits)
        
        ##yy=ifelse(abs(Y-yfits)<=C*sigma.mad, (Y-yfits)/sigma.mad, sign(Y-yfits)*C)
        ##above changed to the below 6/3/15 by wp:
        yy=ifelse(abs(Y-yfits)<=C*sigma.mad, (Y-yfits), sign(Y-yfits)*C*sigma.mad)
        U<-t(XUs)%*%yy
        
        ##U<-t(XUs) %*% (Y - yfits)
        ##CovS<-matrix(0, nrow=k, ncol=k)
        ##for(i in 1:n)
        ##  CovS<-CovS + (sigma0^2)*XUs[i,] %*% t(XUs[i,])
    }
    ## test stat's:
    Ts<-rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf){
            Ts[j] = sum(U^pow[j]) 
        } else { Ts[j] = max(abs(U))
        }
    }
    
    ## bootstrap: 
    T0s = T02s = matrix(0, nrow=B, ncol=length(pow))
    Y0=Y
    for(b in 1:B){
        if (is.null(cov)) {
            Y0 <- sample(Y, length(Y))
#########Null score vector:
##U0<-t(Xg) %*% (Y0-mean(Y0))
            
            ##sigma.mad0=mad(Y0-mean(Y0))
            ##yy=ifelse(abs(Y0-mean(Y0))<=C*sigma.mad0, (Y0-mean(Y0))/sigma.mad0, sign(Y0-mean(Y0))*C)
            ##above changed to the below 6/3/15 by wp:
            yy=ifelse(abs(Y0-mean(Y0))<=C*sigma.mad, (Y0-mean(Y0)), sign(Y0-mean(Y0))*C*sigma.mad)
            U0<-t(Xg)%*%yy
        }
     else{
         ## with nuisance parameters:
         ##Y0 <- yfits + rnorm(n, mean=0, sd=sigma0 )
         Y0 <- yfits + sample(yresids, n, replace = T )
         tdat0<-data.frame(trait=Y0, cov)
         fit0<-glm(trait~., data=tdat0)
         yfits0<-fitted.values(fit0)
         ##U0<-t(XUs) %*% (Y0 - yfits0)
         
         ##sigma.mad0=mad(Y0-yfits0)
         ##yy=ifelse(abs(Y0-yfits0)<=C*sigma.mad0, (Y0-yfits0)/sigma.mad0, sign(Y0-yfits0)*C)
         ##above changed to the below 6/3/15 by wp:
         yy=ifelse(abs(Y0-yfits0)<=C*sigma.mad, (Y0-yfits0), sign(Y0-yfits0)*C*sigma.mad)
         
         U0<-t(XUs)%*%yy
     }
        
     # test stat's:
        for(j in 1:length(pow))
            if (pow[j] < Inf){
                T0s[b, j] = sum(U0^pow[j]) 
            } else { T0s[b, j] = max(abs(U0))
            }
        
    }
    
    ## bootstrap-based p-values:
    pPerm = pPerm0 = rep(NA, length(pow));
    pvs = NULL;
    
    for(j in 1:length(pow)) {
        pPerm0[j] = sum( abs(Ts[j]) <= abs(T0s[,j]))/B
    }
    P0s = PermPvs(T0s)
    minP0s = apply(P0s, 1, min)
    minP =  sum( min(pPerm0) >= minP0s )/B
    pvs<-c(pPerm0, minP)
    names(pvs) <- c(paste("SPUr", pow, sep=""), "aSPUr")
    return(pvs) 
}  


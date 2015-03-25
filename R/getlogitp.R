#' Get logistic p-value
#'
#' Get p-value using logistic regresion for each SNPs
#'
#' @param Y Response or phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele. Matrix with dimension n by g (n : number of observation, p : number of genotype data)
#'
#' @export
#' @return p-values for each SNPs.
#'
#' @examples
#'
#' simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
#'                            LOR=.2, rholim=c(0,0),
#'                            n=100, MAFlim=c(0.05, 0.4), p0=0.05)
#' logitp <- getlogitp(simula$Y, simula$X)
#'
#'
#' @seealso \code{\link{aSPUw}}


getlogitp <- function(Y, X) {
    n.Y <- length(Y)
    nc.X <- ncol(X)

    logitPs <- NULL
    for(i in 1:nc.X) { # i  =1
        tdat1<-data.frame(trait=Y, X = X[,i])
        fit1<-glm(trait~.,family="binomial",data=tdat1)
        logitPs <- c(logitPs, summary(fit1)$coeff[2,4])
    }
    logitPs
}

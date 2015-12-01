#' An Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics.(SNP based)
#'
#' SNP based Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics.
#'
#' @param Z matrix of summary Z-scores, SNPs in rows and traits in columns. Or a vector of summary Z-scores for a single snp
#'
#' @param v estimated covariance matrix based on the summary Z-scores (output of estcov)
#'
#' @param pow power used in SPU test. A vector of the powers.
#'
#' @param tranform if TRUE, the inference is made on transformed Z
#' 
#' @param B number of Monte Carlo samples simulated to compute p-values
#'
#' @return compute p-values for SPU(gamma) i.e. pow=1:8, and infinity
#'             aSPU, based on the minimum p-values over SPU(power)
#'             each row for single SNP
#'

#' @author Junghi Kim, Yun Bai and Wei Pan 
#'
#' @references
#' Junghi Kim, Yun Bai and Wei Pan (2016) An Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics, Genetic Epidemiology DOI 10.1002/gepi.21931
#'
#' @examples
#'
#' # -- n.snp: number of SNPs
#' # -- n.trait: number of traits
#' # -- n.subject: number of subjects
#'
#' n.snp <- 100
#' n.traits <- 10
#' n.subjects <- 1000
#' traits <- matrix(rnorm(n.subjects*n.traits), n.subjects, n.traits)
#' v <- cov(traits)
#' allZ <- rmvnorm(n.snp, sigma=v)
#' colnames(allZ) <- paste("trait", 1:n.traits, sep="")
#' rownames(allZ) <- paste("snp", 1:n.snp, sep="")
#'
#' 
#' v <- estcov(allZ)
#' MTaSPUs(Z = allZ, v = v, B = 100, pow = c(1:4, Inf), tranform = FALSE)
#' MTaSPUs(Z = allZ[1,], v = v, B = 100, pow = c(1:4, Inf), tranform = FALSE)
#' minP(Zi= allZ[1,], v = v)
#'
#' @seealso \code{\link{minP}} \code{\link{estcov}}

MTaSPUs <- function(Z, v, B, pow, tranform = FALSE){
    if (tranform){
        v <- ginv(v)
        Z <- tcrossprod(Z, v)
    }

    if (is.vector(Z)) {
        N <- 1; K <- length(Z)
    }else {
        N <- dim(Z)[1]; K <- dim(Z)[2]
    }

    set.seed(1000)
    Z0 <- rmvnorm(B, mean = rep(0, nrow(v)), sigma = v)
    pval <- matrix(0, N, length(pow)+1)

    ponum <- pow[pow < Inf]
    ## SPU for integer power
    for(k0 in 1:length(ponum)){
        k <- ponum[k0]
        if (N == 1) {
            z1 <- abs(sum(Z^k))
        } else {
            z1 <- abs(rowSums(Z^k))
        }
        z0b <- abs(rowSums(Z0^k))
        for(i in 1:N){
            pval[i,k0] <- (1+sum(z0b>z1[i]))/(B+1)
        }
    }

    ## SPU(max)
    if (Inf %in% pow){
        if(N == 1) {
            z1 <- max(abs(Z))
        } else {
            z1 <- rowMaxs(abs(Z))
        }
        z0 <- rowMaxs(abs(Z0))
        for(i in 1:N){
            pval[i,length(pow)] = (1+sum(z0>z1[i]))/(B+1)
        }
    }

    ## aSPU
    if (N == 1) {
        p1m <- min(pval[,1:length(pow)])
    } else {
        p1m <- rowMins(pval[,1:length(pow)])
    }
    p0 <- matrix(NA, B, length(pow))
    for(k0 in 1:length(ponum)){
        k <- ponum[k0]
        zb <- abs(rowSums(Z0^k))
        p0[,k0] <- (1+B-rank(abs(zb)))/B
    }


    if (Inf %in% pow){
        zb <- rowMaxs(abs(Z0))
        p0[,length(pow)] = (1+B-rank(zb))/B
    }
    p0m = rowMins(p0)
    for(i in 1:N){
        pval[i,length(pow)+1] = (1+sum(p0m<p1m[i]))/(B+1)
    }

    if(Inf %in% pow) s <- c(paste("SPU", ponum, sep=""),"SPUInf","aSPU") else {
                                                                             s <- c(paste("SPU", ponum, sep=""),"aSPU")}
    colnames(pval) <- s
    rownames(pval) <- rownames(Z)
    return(pval)
}

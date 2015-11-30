#' An Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics.(SNP based)
#'
#' SNP based Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics.
#'
#' @param Z matrix of summary Z-scores, SNPs in rows and traits in columns. Or a vector of summary Z-scores for a single snp
#'
#' @param v estimated covariance matrix based on the summary Z-scores (output of estcov)
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
#' sumaspu(Z = allZ, v = v, B = 100, pow = c(1:4, Inf), tranform = FALSE)
#' sumaspu(Z = allZ[1,], v = v, B = 100, pow = c(1:4, Inf), tranform = FALSE)
#' minP(Zi= allZ[1,], v = v)
#'
#' @seealso \code{\link{minP}} \code{\link{estcov}}

sumaspu = function(Z, v, B){

    N = dim(Z)[1]; K = dim(Z)[2]
    set.seed(1000)
    Z0 = rmvnorm(B, sigma=v)
    pval = matrix(NA, N,10)

    ## SPU(1:8)
    for(k in 1:8){
        zb = abs(rowSums(Z^k))
        z0b = abs(rowSums(Z0^k))
        for(i in 1:N){
            pval[i,k] = (1+sum(z0b>zb[i]))/(B+1)
        }
    }

    ## SPU(max)
    z1 = apply(abs(Z), 1, max)
    z0 = apply(abs(Z0), 1, max)
    for(i in 1:N){
        pval[i,9] = (1+sum(z0>z1[i]))/(B+1)
    }

    ## aSPU
    p1m = apply(pval[,1:9], 1, min)
    p0 = matrix(NA, B,9)
    for(k in 1:8){
        zb = abs(rowSums(Z0^k))
        p0[,k] = (1+B-rank(abs(zb)))/B
    }
    zb = apply(abs(Z0), 1, max)
    p0[,9] = (1+B-rank(zb))/B
    p0m = apply(p0, 1, min)
    for(i in 1:N){
        pval[i,10] = (1+sum(p0m<p1m[i]))/(B+1)
    }
    colnames(pval) <- c(paste("SPU", 1:9, sep=""), "aSPU")
    rownames(pval) <- rownames(Z)
    return(pval)
}

#' estcov
#'
#' Estimated covariance matrix based on the summary Z-scores
#'
#' @param allZ matrix of summary Z-scores for all SNP. each row for SNP; each column for single trait
#'
#' @return estimated correlation matrix
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
#' @seealso \code{\link{sumaspu}} \code{\link{minP}}

estcov <- function(allZ){
    n.snp <- dim(allZ)[1]
    n.traits <- dim(allZ)[2]
    minZ = rowMins(2*pnorm(-abs(allZ)))
    nullSNPs = which(minZ>(0.05/(n.snp/n.traits)))
    return( cor(allZ[nullSNPs,]))
}

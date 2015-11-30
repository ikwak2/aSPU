#' minP test.
#'
#' Return exact minP test. 
#'
#' @param Zi a vector of summary Z-scores for single SNP
#'
#' @param v estimated covariance matrix based on the summary Z-scores (output of estcov)
#'
#' @return return exact minP test
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
#' @seealso \code{\link{estcov}} \code{\link{sumaspu}}

minP <- function(Zi, v){
    n <- dim(v)[1]
    x <- as.numeric(max(abs(Zi)))
    return(as.numeric(1 - pmvnorm(lower=c(rep(-x,n)), upper=c(rep(x,n)), mean=c(rep(0,n)), sigma=v)))
}


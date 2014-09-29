#' Sum of powered score (SPU) test
#'
#' It gives the p-values of the SPS test and aSPU test.
#'
#' @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#' or It can be any quantitative traits.
#'
#' @param X genotype data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele.
#'
#' @param cov covariates
#'
#' @param resample Use "perm" for residual permutations and "boot" for parametric bootstrap
#'
#' @param model Use "gaussian" for quantitative trait, and use "binomial" for binary trait.
#'
#' @param n.perm number of permutation
#'
#' @param version "1" for the use of vector T0. "2" for the use of matrix T0. Generally matrix version is faster but when n.perm is so big "2" does not work. This case we should use "1"
#'
#' @export
#' @return Test Statistics and p-values for SPS tests and aSPU test.
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot", model = "binary", pow = c(1:8, Inf), n.perm = 1000, version = 2)
#' out
#'
#' @seealso \code{\link{aSPUperm}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPU <- function(Y, X, cov=NULL, resample = c("boot", "perm"), model=c("gaussian", "binary"), pow = c(1:8, Inf), n.perm = 1000, version = 2:1 ) {

    model <- match.arg(model)
    resample <- match.arg(resample)
    version <- match.arg(version)

    if(resample == "boot") {
        if(version == 2) {
            aSPUboot2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm)
        } else {
            aSPUboot(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm)
        }
    } else {
        if(version == 2) {
            aSPUperm2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm)
        } else {
            aSPUperm(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm)
        }
    }

}


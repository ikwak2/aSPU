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
#' @param version "vec" for the use of vector in permutation. "mat" for the use of matrix in permutation. Generally matrix version is faster but when n.perm is so big "mat" version does not work. This case we should use "vec" version.
#'
#' @export
#' @return Test Statistics and p-values for SPS tests and aSPU test.
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot", model = "binomial", pow = c(1:8, Inf), n.perm = 1000, version = "mat")
#' out
#'
#' @seealso \code{\link{aSPUperm}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPU <- function(Y, X, cov=NULL, resample = c("boot", "perm"), model=c("gaussian", "binomial"), pow = c(1:8, Inf), n.perm = 1000, version = c("mat","vec") ) {

    model <- match.arg(model)
    resample <- match.arg(resample)
    version <- match.arg(version)

    if(resample == "boot") {
        if(version == "mat") {
            aSPUboot2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        } else {
            aSPUboot(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        }
    } else {
        if(version == "mat") {
            aSPUperm2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        } else {
            aSPUperm(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        }
    }

}

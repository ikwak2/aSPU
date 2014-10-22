#' Sum of powered score (SPU) test
#'
#' It gives the p-values of the SPU test and aSPU test.
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
#' @param pow power used in SPU test.
#'
#' @param n.perm number of permutation
#'
#' @param userank use similar code with original version if T, by definition if F
#'
#' @export
#' @return Test Statistics and p-values for SPU tests and aSPU test.
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot", model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
#' out
#'
#' @seealso \code{\link{aSPUperm}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPU <- function(Y, X, cov=NULL, resample = c("perm", "boot"), model=c("gaussian", "binomial"), pow = c(1:8, Inf), n.perm = 1000, userank = T ) {

    model <- match.arg(model)
    resample <- match.arg(resample)

    if(resample == "boot") {
        if(n.perm > 10^5) {
            aSPUboot(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        } else {
            aSPUboot2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        }
    } else {
        aSPUpermC(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model, userank = userank)
    }

}

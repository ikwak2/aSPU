#' Variance-weighted adaptive Sum of powered score tests. (SPUw and aSPUw)
#'
#' It gives the p-values of the SPUw test and aSPUw test based
#' on based on the permutation of residuals or simulation from its distripution.
#'
#' @param Y Response or phenotype data. It can be disease lables; =0 for controls, =1 for cases.
#'     or It can be any quantitative traits. Vector with length n (number of observations)
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP. The value of each element is the # of the copies
#'     for an allele. Matrix with dimension n by g (n : number of observation, p : number of genotype data)
#'
#' @param cov covariates. Matrix with dimension n by k (n :number of observation, k : number of covariates)
#'
#' @param resample Use "perm" for residual permutations and "sim" for simulations from the distribution.
#'
#' @param model Use "gaussian" for quantitative trait, and use "binomial" for binary trait.
#'
#' @param pow power used in SPUw test. Vector of g number of power.
#'
#' @param n.perm number of permutations or bootstraps
#'
#' @param userank use similar code with original version if T, by definition if F
#'
#' @export
#' @return Test Statistics and p-values for SPUw tests and aSPUw test.
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPUW(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
#' out
#'
#' @seealso \code{\link{aSPU}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPUw <- function(Y, X, cov=NULL, resample = c("perm", "asymp"), model=c("gaussian", "binomial"), pow = c(1:8, Inf), n.perm = 1000, userank = T ) {

    model <- match.arg(model)
    resample <- match.arg(resample)

    if(resample == "asymp") {
        aSPUwsim(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
    } else {
        aSPUwpermC(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model, userank = userank)
    }
}

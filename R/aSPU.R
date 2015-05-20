#' Sum of Powered Score (SPU) tests and adaptive SPU (aSPU) test.
#'
#' It gives p-values of the SPU tests and aSPU test.
#'
#' @param Y Response or phenotype data. It can be disease indicators; =0 for controls, =1 for cases.
#' Or it can be a quantitative trait. A vector with length n (number of observations)
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP (or a predictor). The value of each SNP is the # of the copies
#'     for an allele. A matrix with dimension n by p (n : number of observation, p : number of genotype data)
#'
#' @param cov Covariates. A matrix with dimension n by k (n :number of observation, k : number of covariates)
#'
#' @param resample Use "perm" for residual permutations, "sim" for simulations from the null distribution, and "boot" for parametric bootstrap.
#'
#' @param model Use "gaussian" for quantitative trait, and use "binomial" for binary trait.
#'
#' @param pow power used in SPU test. A vector of the powers.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @return A list object, Ts : Test Statistics for the SPU and aSPU test.
#'         pvs : p-values for the SPU and aSPU test.
#'
#' @author Il-Youp Kwak, Junghi Kim, Yiwei Zhang and Wei Pan
#'
#' @references
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
#' A powerful and adaptive association test for rare variants,
#' Genetics, 197(4), 1081-95
#'
#' @examples
#'
#' data(exdat)
#'
#' ## example analysis using aSPU test on exdat data.
#' out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot",
#'            model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
#'
#' out$Ts
#' # These are list of Test Statistics for SPU and aSPU tests.
#' out$pvs
#' # These are p-values of Test Statistics for SPU and aSPU tests.
#'
#' @seealso \code{\link{aSPUw}}


aSPU <- function(Y, X, cov=NULL, resample = c("perm", "sim", "boot"), model=c("gaussian", "binomial"), pow = c(1:8, Inf), n.perm = 1000) {

    model <- match.arg(model)
    resample <- match.arg(resample)

    if(resample == "boot") {
        if(n.perm > 10^5) {
            aSPUboot(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        } else {
            aSPUboot2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        }
    } else {

        if(resample == "sim") {
            if(n.perm > 10^5) {
                aSPUsim1(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
            } else {
                aSPUsim2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
            }
        } else {
            userank = T;
            aSPUpermC(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model, userank = userank)
        }
    }
}

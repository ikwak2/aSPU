#' Sum of Powered Score (SPUs) tests and adaptive SPU (aSPUs) test with GWAS summary statistics. (C coded)
#'
#' It gives p-values of the SPUs tests and aSPUs test with GWAS summary statistics.
#'
#' @param Zs Z-scores for each SNPs. It could be P-values if the Ps option is TRUE.
#'
#' @param corrSNP Correaltion matirx of SNPs. Estimated from the reference population.
#'
#' @param pow power used in SPU test. A vector of the powers.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @param Ps TRUE if input is p-value, FALSE if input is Z-scores. The default is FALSE.
#'
#' @return A list object, Ts : test statistics for the SPU tests (in the order of the specified pow) and finally for the aSPU test.
#'         pvs : p-values for the SPU and aSPU tests.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
#' A powerful and adaptive association test for rare variants,
#' Genetics, 197(4), 1081-95
#'
#' Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014) Comparison of statistical tests for group differences in brain functional networks, NeuroImage, 1;101:681-694
#'
#' @examples
#'
#' data(kegg9)
#' ## example analysis using aSPUM test.
#' g <- kegg9$gene.info[1,1]  #  SOAT1
#' ## Take snps mapped on gene "SOAT1" from the information of gene.info and snp.info.
#' snps <- which( ( kegg9$snp.info[,2] == kegg9$gene.info[kegg9$gene.info[,1] == g, 2] ) &
#'                  (kegg9$snp.info[,3] > kegg9$gene.info[kegg9$gene.info[,1] == g, 3] ) &
#'                  (kegg9$snp.info[,3] < kegg9$gene.info[kegg9$gene.info[,1] == g, 4] )  )
#' ## Take subsets
#' newP <- kegg9$nP[snps] ;
#' ldsub <- kegg9$ldmatrix[snps, snps];
#' ## Get p-value for gene SOAT1. Read vignette for details.
#' out <- aSPUsC(newP, corrSNP=ldsub , pow=c(1,2,4,8, Inf), n.perm=100, Ps=TRUE)
#'
#' out$Ts
#' # This is a vector of Test Statistics for SPUM and aSPUM tests.
#' # SPU1 to SPUInf corresponds with the option pow=c(1:8, Inf)
#' # They are SPU test statistics.
#' # The last element aSPUM is minimum of them, aSPUM statistic.
#'
#' out$pvs
#' # This is a vector of p-values for SPUM and aSPUM tests.
#' # SPUM1 to SPUMInf corresponds with the option pow=c(1:8, Inf)
#' # They are p-values for corresponding SPUM tests.
#' # The last element is p-value of aSPUM test.
#'
#' @seealso \code{\link{aSPUw}} \code{\link{aSPU}}

aSPUsC <- function(Zs, corrSNP, pow = c(1,2,4,8, Inf), n.perm = 1000, Ps = FALSE)
{
    n <- length(Zs)
    k <- n
    if(Ps == TRUE)
      Zs <- qnorm(1 - Zs/2)

    U <- Zs
    CovS <- corrSNP

    eS <- eigen(CovS, symmetric = TRUE)
    ev <- eS$values

    CovSsqrt <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), k)

#    svd.CovS <- svd(CovS)
#    CovSsqrt <- svd.CovS$u %*% diag(sqrt(svd.CovS$d))

    Ts = rep(NA, length(pow))
    npow = pow
    for (j in 1:length(pow)) {
        if (pow[j] < Inf)
            Ts[j] = sum(U^pow[j])
        else {
            Ts[j] = max(abs(U))
            npow[j] = 0
        }
    }

    n_pow = length(pow)
    n_Zs = length(Zs)
    n_perm = n.perm


    CovSsqrt %*% rnorm(n_Zs)

    A <- matrix(rnorm(n_Zs*n_perm), n_Zs, n_perm)


    rnms <- c(apply(A, 2, function(x) CovSsqrt %*% x))
    ##    rnorm(n_Zs*n_perm)
    ##    seeds = sample(1:10000000,n_Zs*n_perm)
    Pss <- ifelse(Ps, 1,0)

#    dyn.load("R_get_pvs_s.so")
    output <- .C("R_get_pvs_s",
                 as.double(Zs),
                 as.double(Ts),
                 as.double(npow),
                 as.double(rnms),
 #                as.integer(seeds),
                 as.integer(n_pow),
                 as.integer(n_Zs),
                 as.integer(n_perm),
                 as.integer(Pss),
                 pvs = as.double( rep(0,n_pow + 1) ),
                 PACKAGE="aSPU")

    pvs <- output$pvs

    Ts <- c(Ts, min(pvs[1:length(pow)]))
    names(Ts) <- c(paste("SPUs", pow, sep = ""), "aSPUs")
    names(pvs) = names(Ts)
    list(Ts = Ts, pvs = pvs)
}



#' Sum of Powered Score (SPU) tests and adaptive SPU (aSPU) test for meta-analyzed data.
#'
#' It gives p-values of the SPU tests and aSPU test for meta-analyzed data.
#'
#' @param Zs Z-scores for each SNPs. It could be P-values if the Ps option is TRUE. 
#'
#' @param CorrSNP Correaltion matirx of SNPs. Estimated from the reference population.
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
#' g <- gene.info[1,1]  #  SOAT1 
#' ## Take snps mapped on gene "SOAT1" from the information of gene.info and snp.info. 
#' snps <- which( ( snp.info2[,2] == gene.info[gene.info[,1] == g, 2] ) &
#'                  (snp.info2[,3] > gene.info[gene.info[,1] == g, 3] ) &
#'                  (snp.info2[,3] < gene.info[gene.info[,1] == g, 4] )  )
#' ## Take subsets
#' newP <- nZ[snps] ;
#' ldsub <- ldmatrix[snps, snps];
#' ## Get p-value for gene SOAT1. Read vignette for details.
#' out <- aSPUM(newP, corrSNP=ldsub , pow=c(1,2,4,8, Inf), n.perm=100, Ps=TRUE)
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

aSPUM <- function(Zs, corrSNP, pow = c(1,2,4,8, Inf), n.perm = 1000, Ps = FALSE)
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
    for (j in 1:length(pow)) {
        if (pow[j] < Inf)
            Ts[j] = sum(U^pow[j])
        else Ts[j] = max(abs(U))
    }
    pPerm0 = rep(NA, length(pow))
    T0s = numeric(n.perm)
    s <- sample(1:10^5, 1)
    for (j in 1:length(pow)) {
        set.seed(s)
        for (b in 1:n.perm) {
            U00 <- rnorm(k, 0, 1)
            U0 <- CovSsqrt %*% U00
            if (Ps == TRUE)
              U0 <- abs(U0)
            if (pow[j] < Inf) {
	       T0s[b] = round(sum(U0^pow[j]), digits = 8)
            }
            if (pow[j] == Inf) {
                T0s[b] = round(max(abs(U0)), digits = 8)
            }
        }
        pPerm0[j] = round(sum(abs(Ts[j]) <= abs(T0s))/n.perm,
            digits = 8)
        P0s = ((n.perm - rank(abs(T0s))) + 1)/(n.perm)
        if (j == 1)
            minp0 = P0s
        else minp0[which(minp0 > P0s)] = P0s[which(minp0 > P0s)]
    }
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPUstat", pow, sep = ""), "aSPUstat")
    names(pvs) = names(Ts)
    list(Ts = Ts, pvs = pvs)
}










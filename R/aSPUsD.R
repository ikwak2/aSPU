#' Sum of Powered Score (SPUs) tests and adaptive SPU (aSPUs) test for single trait - SNP set association with GWAS summary statistics (distribution based).
#'
#' It gives p-values of the SPUs tests and aSPUs test with GWAS summary statistics.
#'
#' @param Zs Z-scores for each SNPs. It could be P-values if the Ps option is TRUE. 
#'
#' @param corrSNP Correlation matirx of the SNPs to be tested; estimated from a
#' reference panel (based on the same set of the reference alleles as
#' used in calculating Z-scores).
#'
#' @param Ps TRUE if input is p-value, FALSE if input is Z-scores. The default is FALSE.
#'
#' @return pvs : p-values for the SPUsD and aSPUsD tests.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#'
#' Gongjun Xu, Lifeng Lin, Peng Wei and Wei Pan (2016)
#' An adaptive two-sample test for high-dimensional means, Biometrika (2016) 103 (3): 609-624.
#' Il-Youp Kwak, Wei Pan (2015) Adaptive Gene- and Pathway-Trait Association Testing with GWAS Summary Statistics,
#' Bioinformatics, 32(8), 1178-1184
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
#' out <- aSPUsD(newP, corrSNP=ldsub, Ps=TRUE)
#'
#' out
#'
#' @seealso \code{\link{aSPUs}} \code{\link{aSPU}}
#'
#' @export

aSPUsD <- function(Zs, corrSNP, Ps = FALSE)
{
    n <- length(Zs)
    k <- n
    if(Ps == TRUE)
      Zs <- qnorm(1 - Zs/2)

    U <- Zs
    CovS <- corrSNP

    p1 = SumSqU(U, CovS)
    p2 = Sum(U, CovS)
    p3 = UminPd(U, CovS)
    mp = min(p1,p2,p3)
    p = 1 - (1-mp)^3

    pvs = c(p1,p2,p3,p)
    names(pvs) = c("Psum","Pssu","PminP","PaSPU")
    pvs
}

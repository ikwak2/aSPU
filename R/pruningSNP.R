#' Prune SNPs (related to aSPUs, aSPUsPath, MTaSPUsSet)
#'
#' When correlation matirx have highly correlated SNPs, the performance of aSPUs, aSPUsPath and MTaSPUsSet are not very good. Do pruning using this function.
#'
#' @param corSNP Correlation matirx of the SNPs to be tested; estimated from a
#' reference panel (based on the same set of the reference alleles as
#' used in calculating Z-scores).
#'
#' @param rup pruning criteria, erase one SNP when correlation between two SNPs are larger than this value.
#'
#' @param rdown pruning criteria, erase one SNP when correlation between two SNPs are smaller than this value.
#' 
#' @return a list object pruned.corSNP and to.erase. pruned.corSNP is pruned correlation matrix. to.erase is SNP index to erase to get purned object. (i.e. corSNP[-to.erase, -to.erase] = pruned.corSNP )
#'
#' @author Il Youp Kwak
#'
#' @examples
#'
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
#'
#' prSNP <- pruneSNP(ldsub) 
#' newP <- newP[-prSNP$to.erase]
#' ldsub <- ldsub[-prSNP$to.erase, -prSNP$to.erase ]
#' 
#' ## Get p-value for gene SOAT1. Read vignette for details.
#' out <- aSPUs(newP, corSNP=ldsub , pow=c(1:8, Inf), n.perm=100, Ps=TRUE)
#'
#' out$Ts
#' # This is a vector of Test Statistics for SPUM and aSPUM tests.
#' # SPUs1 to SPUsInf corresponds with the option pow=c(1:8, Inf)
#' # They are SPUs test statistics.
#' # The last element aSPUs is minimum of them, aSPUs statistic.
#'
#' out$pvs
#' # This is a vector of p-values for SPUs and aSPUs tests.
#' # SPUs1 to SPUsInf corresponds with the option pow=c(1:8, Inf)
#' # They are p-values for corresponding SPUs tests.
#' # The last element is p-value of aSPUs test.
#'
#' @seealso \code{\link{aSPUs}} \code{\link{aSPUsPath}} \code{\link{MTaSPUsSet}}
#'
#' @export

pruneSNP <- function(corSNP, rup = .95, rdown = -.95) {

    if(!is.matrix(corSNP)) {
        stop("corSNP should be a matrix")
    }
    if(ncol(corSNP) != nrow(corSNP) )
        stop("corSNP is a correlation matirx. the number of rows and columns should match")

    if(ncol(corSNP) <= 1 )
        stop(" only 1 SNP, nothing to prune.")
    
    ers <- NULL
    while( sum( corSNP < -0.95 ) > 0 ) {
        toerase <- which( apply( corSNP < rdown, 1, sum) == max( apply( corSNP < rdown, 1, sum) ) )
        corSNP <- corSNP[-toerase[1], -toerase[1]]
        ers <- c(ers, toerase[1])
        if( is.matrix(corSNP) == FALSE) break
    }
        
    while( sum( corSNP > 0.95 ) > ncol(corSNP) ) {
        toerase <- which( apply( corSNP > rup, 1, sum) == max( apply( corSNP > rup, 1, sum) ) )
        corSNP <- corSNP[-toerase[1], -toerase[1]]
        ers <- c(ers, toerase[1])
        if( is.matrix(corSNP) == FALSE) break
    }

    list(corSNP = corSNP, to.erase = ers)
}


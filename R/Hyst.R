#' HYST
#'
#' Get p-value using GATES-Simes method.
#'
#' @param pvec p-values for each SNP.
#'
#' @param ldmatrix numeric. A correlation matrix of SNPs, dimensions matching the p and snps arguments.
#'
#' @param snp.info A matrix with snp data. Data are coded as 0, 1, 2,
#' corresponding to homozygotes for the major allele,
#' heterozygotes, and homozygotes for the minor allele,
#' respectively. Each row of the matrix is one SNP and each
#' column represents one subject.
#'
#' @param gene.info A matrix with four columns with the order of gene.Name,
#' chr, Start, and End. The name of the gene is under gene.Name.
#' The chromosome is chr, and the start and end coordinate
#' (start < end) are in the Start and End columns. It is suggested
#' that the genes are sorted by their locations in the genome.
#'
#' @export
#' @return p-values for each SNPs.
#'
#' @author Il-Youp Kwak
#'
#' @references
#' Miao-Xin Li, Johnny S.H. Kwan and Pak C. Sham (2012)
#' HYST: A Hybrid Set-Based Test for Genome-wide Association Studies, with Application to Protein-Protein Interaction-Based Association Analysis
#'
#' @examples
#'
#' simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
#'                            LOR=.2, rholim=c(0,0),
#'                            n=100, MAFlim=c(0.05, 0.4), p0=0.05)
#' logitp <- getlogitp(simula$Y, simula$X)
#'
#' ## get correlation of SNPs using controls
#' ldmat <- cor(simula$X[ simula$Y == 0, ])
#' out <- Hyst(pvec = logitp, ldmatrix = ldmat, snp.info = simula$snp.info,
#'             gene.info = simula$gene.info)
#'
#' @seealso \code{\link{GatesSimes}} \code{\link{GATES2}}


Hyst <- function(pvec, ldmatrix, snp.info, gene.info) {
    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) { #g = 6
        if( sum(snp.info[,2]==gene.info[g,1]) == 1 ) {
            Pg <- pvec[snp.info[,2]==gene.info[g,1]]
            keyG <- 1
            M <- 1
        } else {
            pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
            M <- length(pvecG)
            ldmat <- ldmatrix[snp.info[,2]==gene.info[g,1], snp.info[,2]==gene.info[g,1]]
#            Xn <- X[controls, snp.info[,2]==gene.info[g,1]]
            o.pv <- order(pvecG)
            ldmat2 <- ldmat[o.pv,o.pv]
            out <- GATES2(ldmatrix = ldmat2, p = sort(pvecG) )
            keyG <- out[2]
            Pg <- out[1]
        }

        PGs <- c(PGs, Pg)
        keyGs <- c(keyGs, gpos + keyG)
        gpos = gpos + M

    }
#    cat("keyGs :", keyGs);
#    cat("\n");
#    cat("PGs :", PGs);

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(ldmatrix[keyGs[i], keyGs[i]])
#            r = abs(cor(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)
    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)
#    Physt2 = 1 - pchisq( Hyst , df = 2*n.gene)
}

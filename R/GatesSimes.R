#' GATES-Simes method
#'
#' Get p-value using GATES-Simes method.
#'
#' @param pvec p-values for each SNP.
#'
#' @param ldmatrix The correlation matrix among SNPs,
#'
#' @param snp.info SNP information matrix, The 1st column is SNP ids, 2nd column is SNP chromosome, 3rd column indicate SNP location.
#'
#' @param gene.info GENE information matrix, The 1st column is SNP ids, 2nd column is SNP chromosome, 3rd column indicate SNP location.
#'
#' @export
#' @return p-values for each SNPs.
#'
#' @author Il-Youp Kwak
#'
#' @references
#' Hongsheng Gui, Miaoxin Li, Pak C Sham and Stacey S Cherny (2011)
#' Comparisons of seven algorithms for pathway analysis using the WTCCC Crohn's Disease
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
#' out <- GatesSimes(pvec = logitp, ldmatrix = ldmat, snp.info = simula$snp.info,
#'                   gene.info = simula$gene.info)
#'
#'
#' @seealso \code{\link{Hyst}}  \code{\link{GATES2}}


GatesSimes <- function(pvec, ldmatrix, snp.info, gene.info) {

    gene.info <- gene.info[ gene.info[,4] != 0, ]
    n.gene <- nrow(gene.info)

    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0
    for(g in 1:n.gene) { # g = 1

        if( sum(snp.info[,2]==gene.info[g,1]) == 1 ) {
            Pg <- pvec[snp.info[,2]==gene.info[g,1]]
        } else {
            pvecG <- pvec[snp.info[,2]==gene.info[g,1]]
            ldmat <- ldmatrix[snp.info[,2]==gene.info[g,1], snp.info[,2]==gene.info[g,1]]
#            Xn <- X[controls, snp.info[,2]==gene.info[g,1]]
            o.pv <- order(pvecG)
            ldmat2 <- ldmat[o.pv,o.pv]
#            ldmat <- cor(Xn[,o.pv])
            Pg <- GATES2(ldmatrix = ldmat2, p = sort(pvecG) )[1]
        }

        PGs <- c(PGs, Pg)
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )
    PgatesSime
}

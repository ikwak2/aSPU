#' GATES-Simes method
#'
#' Get p-value using GATES-Simes method.
#'
#' @param pvec p-values for each SNP.
#'
#' @param ldmatrix numeric. A correlation matrix of SNPs, dimensions matching the p and snps arguments.
#'
#' @param snp.info SNP information matrix, The 1st column is SNP ids, 2nd column is SNP chromosome, 3rd column indicate SNP location.
#'
#' @param gene.info GENE information matrix, The 1st column is GENE ids, 2nd column is GENE chromosome, 3rd and 4th column indicate where the gene location starts and ends.
#'
#' @export
#' @return p-value based on Gates-Simes method.
#'
#' @author Il-Youp Kwak, Peng Wei and Wei Pan
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

    n.gene <- nrow(gene.info)

    GL <- NULL;
    for(g in 1:n.gene) { # g = 2

        snpTF <- ( snp.info[,2] == gene.info[g,2] &
                      gene.info[g,3] <= as.numeric(snp.info[,3]) &
                          gene.info[g,4] >= as.numeric(snp.info[,3]) )

        if( sum(snpTF) != 0)
            GL[[g]] <- which(snpTF)
    }

    PGs <- NULL;
    for ( i in 1:length(GL) ) { #i = 26
        if ( length(GL[[i]]) == 1 ) {
            Pg <- pvec[ GL[[i]] ]
        } else {
            pvecG <- pvec[ GL[[i]] ]
            ldmat <- ldmatrix[GL[[i]], GL[[i]]]
            o.pv <- order(pvecG)
            ldmat2 <- ldmat[o.pv,o.pv]
            Pg <- GATES2(ldmatrix = ldmat2, p = sort(pvecG) )[1]
        }
        PGs <- c(PGs, Pg)
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )
    PgatesSime
}


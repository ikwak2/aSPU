#' Gene-based Association Test that uses Extended Simes procedure (GATES)
#'
#' Get p-value using GATES method. Usually it used to get genomewise p-values.
#' This function is taken from postgwas package.
#' There is a little modification of the code GATES in postgwas package.
#' 1) The approximated matrix may have negative eigen value, we modified it not to have negative value 2) we add one more return (the key gene location) for Hyst method.
#'
#' @param ldmatrix numeric. A correlation matrix of SNPs,
#'          dimensions matching the p and snps arguments.
#'
#' @param p p-value for each SNPs.
#'
#' @export
#' @return Combined p-values using GATES method and the key gene location.
#'
#' @references Miao-Xin Li, Hong-Sheng Gui, Johnny S.H. Kwan and Pak C. Sham (2011)
#' GATES: A Rapid and Powerful Gene-Based Association Test Using Extended Simes Procedure
#' The American Journal of Human Genetics 88, 283-293
#'
#' @author Milan Hiersche(taken from pastgwas package), Il-Youp Kwak(modified a little)
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
#' (gatesp <- GATES2(ldmat, logitp)[1] )
#'
#'
#' @seealso \code{\link{Hyst}} \code{\link{GatesSimes}}


GATES2 <- function (ldmatrix, p)
{
    snpcount <- length(p)
    if (!all(dim(ldmatrix) == snpcount))
        stop("function GATES: Argument 'ldmatrix' is not rectangular or does not match the length of argument vector 'p'.\n")
    if (any(is.na(ldmatrix)))
        stop("function GATES: Argument 'ldmatrix' may not contain NA values.\n")
    if (snpcount < 1)
        stop("function GATES: No SNP provided.\n")
    if (snpcount == 1)
        return(p)
    ldmatrix <- 0.7723 * ldmatrix^6 - 1.5659 * ldmatrix^5 + 1.201 *
        ldmatrix^4 - 0.2355 * ldmatrix^3 + 0.2184 * ldmatrix^2 +
        0.6086 * ldmatrix
    eff.snpcount.fun <- function(ldmat) {
        ldmat <- as.matrix(ldmat)
        snpcount.local <- dim(ldmat)[1]
        if (snpcount.local <= 1)
            return(1)
        ev <- eigen(ldmat, only.values = TRUE)$values

        if( sum(ev < 0) !=0 ) {
            ev <- ev[ev > 0]
            ev <- ev / sum(ev) * snpcount.local
        }

        ev <- ev[ev > 1]
        snpcount.local - sum(ev - 1)
    }
    eff.snpcount.global <- eff.snpcount.fun(ldmatrix)

    candid <- sapply(1:snpcount, function(i) {
        (eff.snpcount.global * p[i]) / eff.snpcount.fun(ldmatrix[1:i,1:i])
    }
                     )

    Pg <- min(candid)
    keyGloc <- which(min(candid) == candid)[1]
    out <- c(Pg, keyGloc)
    names(out) <- c("Pg", "keyGeneLoc")

    out
}

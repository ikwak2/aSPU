#' Pathway based Sum of Powered Score tests (SPUsPath) and adaptive SPUpath (aSPUsPath) test for single trait - pathway association with GWAS summary statistics. (vector version, fast with large n), Cpp version
#'
#' It gives p-values of the SPUsPath tests and aSPUsPath test with GWAS summary statistics.
#'
#' @param Zs Z-scores for each SNPs. It could be P-values if the Ps option is TRUE. 
#'
#' @param corSNP Correlation matirx of the SNPs to be tested; estimated from a
#' reference panel (based on the same set of the reference alleles as
#' used in calculating Z-scores).
#'
#' @param snp.info SNP information matrix, the 1st column is SNP id, 2nd column is chromosome #, 3rd column indicates SNP location.
#'
#' @param gene.info GENE information matrix, The 1st column is GENE id, 2nd column is chromosome #, 3rd and 4th column indicate start and end positions of the gene.
#'
#' @param pow SNP specific power(gamma values) used in SPUpath test.
#'
#' @param pow2 GENE specific power(gamma values) used in SPUpath test.
#'
#' @param n.perm number of permutations.
#'
#' @param Ps TRUE if input is p-value, FALSE if input is Z-scores. The default is FALSE.
#'
#' @param prune if it is TRUE, do pruing before the test using pruneSNP function. 
#'
#' @export
#' @return P-values for SPUMpath tests and aSPUMpath test.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#' Il-Youp Kwak, Wei Pan (2015)
#' Adaptive Gene- and Pathway-Trait Association Testing with GWAS Summary Statistics,
#' Bioinformatics, 32(8):1178-1184 
#'
#' @examples
#' data(kegg9)
#'
#' # p-values of SPUpath and aSPUpath tests.
#' out.a <- aSPUsPath(kegg9$nP, corSNP = kegg9$ldmatrix, pow=c(1:8, Inf),
#'                   pow2 = c(1,2,4,8), 
#'                   snp.info=kegg9$snp.info, gene.info = kegg9$gene.info,
#'                   n.perm=10, Ps = TRUE)
#'
#' out.a
#'
#' @seealso \code{\link{aSPUs}}


aSPUsPath2C <- function (Zs, corSNP, pow = c(1, 2, 4, 8, Inf), pow2 = c(1, 2, 4, 8),
                         snp.info, gene.info, n.perm = 1000, Ps = FALSE, prune = TRUE) {

    ## some input checking stuff in here
    if( any(is.na(X)) ) {
        stop("NA exist in gene matrix. ")
    }

    if( any(is.na(Y)) ) {
        stop("NA exist in phenotype. ")
    }

    
    
    if (prune == TRUE) {
        pr <- pruneSNP(corSNP)
        if (length(pr$to.erase) > 0) {
            Zs <- Zs[-pr$to.erase]
            corSNP <- corSNP[-pr$to.erase, -pr$to.erase]
            snp.info <- snp.info[snp.info[, 1] %in% names(Zs), 
                ]
        }
    }
    if (length(Zs) <= 1) {
        stop("less than 1 SNP.")
    }
    if (!(length(Zs) == dim(corSNP)[1] & dim(corSNP)[2] == dim(corSNP)[1])) {
        stop("dimension do not match. Check dimension of Zs and corSNP.")
    }
    ko <- length(Zs)
    n.gene <- nrow(gene.info)
    GL <- list(0)
    GLch <- NULL
    i = 1
    for (g in 1:n.gene) {
        snpTF <- (snp.info[, 2] == gene.info[g, 2] & gene.info[g, 
            3] <= as.numeric(snp.info[, 3]) & gene.info[g, 4] >= 
            as.numeric(snp.info[, 3]))
        if (sum(snpTF) != 0) {
            GL[[i]] <- which(snpTF)
            GLch <- c(GLch, gene.info[g, 2])
            i = i + 1
        }
    }
    chrs <- unique(GLch)
    CH <- list(0)
    CH.CovSsqrt <- list(0)
    for (i in 1:length(chrs)) {
        c = chrs[i]
        CH[[i]] <- unlist(GL[which(GLch == c)])
        Covtemp <- corSNP[CH[[i]], CH[[i]]]
        eS <- eigen(Covtemp, symmetric = TRUE)
        ev <- eS$values
        k1 <- length(ev)
        CH.CovSsqrt[[i]] <- eS$vectors %*% diag(sqrt(pmax(ev, 
            0)), k1)
    }
    Zs = Zs[unlist(GL)]
    nSNPs0 = unlist(lapply(GL, length))
    k <- length(Zs)
    if (Ps == TRUE) 
        Zs <- qnorm(1 - Zs/2)
    U <- Zs
    nGenes = length(nSNPs0)
    StdTs <- rep(0, length(pow) * nGenes)
    for (j in 1:length(pow)) {
        for (iGene in 1:nGenes) {
            if (iGene == 1) 
                SNPstart = 1
            else SNPstart = sum(nSNPs0[1:(iGene - 1)]) + 1
            indx = (SNPstart:(SNPstart + nSNPs0[iGene] - 1))
            if (pow[j] < Inf) {
                a = (sum(U[indx]^pow[j]))
                StdTs[(j - 1) * nGenes + iGene] = sign(a) * ((abs(a)/nSNPs0[iGene])^(1/pow[j]))
            } else {
                StdTs[(j - 1) * nGenes + iGene] = max(abs(U[indx]))
            }
        }
        
    }

    Ts2 <- rep(0, length(pow) * length(pow2))
    for (j2 in 1:length(pow2)) {
        for (j in 1:length(pow)) {
            if (pow2[j2] < Inf) {
                Ts2[(j2 - 1) * length(pow) + j] = sum(StdTs[((j - 
                  1) * nGenes + 1):(j * nGenes)]^pow2[j2])
            }
            else {
                Ts2[(j2 - 1) * length(pow) + j] = max(StdTs[((j - 
                  1) * nGenes + 1):(j * nGenes)])
            }
        }
    }

    if(max(pow) == Inf) {
        pow[which(pow ==Inf)] = -1
    }

    if(max(pow2) == Inf) {
        pow2[which(pow2 ==Inf)] = -1
    }
    
    s <- sample(1:10^5, 1)    
    Results = aSPUsPathEngine2(CH, CH.CovSsqrt, pow, pow2, nGenes, n.perm , k, Ps, nSNPs0, Ts2, s)

    minp0 <- Results$minp0
    pPerm0 <- Results$pPerm0
    P0s <- Results$P0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    Ts <- c(Ts2, min(pPerm0))

    if(min(pow) == -1) {
        pow[which(pow == -1 )] = Inf
    }

    if(min(pow2) == -1) {
        pow2[which(pow2 == -1)] = Inf
    }

    nmvec <- NULL
    for (ii in pow2) {
        for (jj in pow) {
            nmvec <- c(nmvec, paste("SPUsPath", jj, ",", ii, 
                sep = ""))
        }
    }
    nmvec <- c(nmvec, "aSPUsPath")
    names(Ts) = nmvec
    names(pvs) = nmvec
    list(Ts = Ts, pvs = pvs)
}


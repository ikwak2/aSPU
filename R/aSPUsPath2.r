#' Pathway based Sum of Powered Score tests (SPUsPath) and adaptive SPUpath (aSPUsPath) test for single trait - pathway association with GWAS summary statistics. (vector version, fast when n is large)
#'
#' It gives p-values of the SPUsPath tests and aSPUsPath test with GWAS summary statistics. Faster than aSPUsPath function when n is large (N > 10^4).
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
#' \dontrun{out.a <- aSPUsPath2(kegg9$nP, corSNP = kegg9$ldmatrix, pow=c(1:8, Inf),
#'                  pow2 = c(1,2,4,8), 
#'                  snp.info=kegg9$snp.info, gene.info = kegg9$gene.info,
#'                  n.perm=1000, Ps = TRUE) }
#'
#' #out.a
#'
#' @seealso \code{\link{aSPUs}}
#'
#' @export

aSPUsPath2 <- function (Zs, corSNP, pow = c(1, 2, 4, 8, Inf), pow2 = c(1, 2, 4, 8),
                        snp.info, gene.info, n.perm = 1000, Ps = FALSE, prune = TRUE) {
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

    ##    s <- sample(1:10^8, n.perm)
    s <- sample(1:10^5)    
    pPerm0 = rep(NA, length(pow)*length(pow2))
    T0s = numeric(n.perm)    
    Ts2t <- rep(0, length(pow) * length(pow2))
    T0st = rep(0, nGenes*length(pow))
    for (j2 in 1:length(pow2)) { # j2 = 1 
        for (j in 1:length(pow)) { # j = 1
            set.seed(s)
            for (b in 1:n.perm) {
##                set.seed(s[b])
                U00 <- rnorm(k, 0, 1)
                U0 <- NULL
                for (ss in 1:length(CH)) {
                    U0 <- c(U0, CH.CovSsqrt[[ss]] %*% U00[CH[[ss]]])
                }
                if (Ps == TRUE) 
                    U0 <- abs(U0)
                
                for (iGene in 1:nGenes) { #iGene = 1
                    
                    if (iGene == 1) {
                        SNPstart = 1
                    } else { SNPstart = sum(nSNPs0[1:(iGene - 1)]) + 1 }
                    indx = (SNPstart:(SNPstart + nSNPs0[iGene] - 1))
                    if (pow[j] < Inf) {
                        a = (sum(U0[indx]^pow[j]))
                        T0st[(j - 1) * nGenes + iGene] = sign(a) * 
                            ((abs(a)/nSNPs0[iGene])^(1/pow[j]))
                    } else T0st[(j - 1) * nGenes + iGene] = max(abs(U0[indx]))
                }
                
                if (pow2[j2] < Inf) {
                    T0s[b] =
                        sum(T0st[((j-1) * nGenes + 1):(j * nGenes)]^pow2[j2])
                } else {
                    T0s[b] =
                        max(T0st[((j-1) * nGenes + 1):(j * nGenes)])
                }
            }
            
            pPerm0[(j2 - 1) * length(pow) + j] = sum(abs(Ts2[(j2 - 1) * length(pow) + j]) <= abs(T0s))/n.perm
            P0s = ((n.perm - rank(abs(T0s))) + 1)/(n.perm)
            if (j == 1 && j2 ==1) 
                minp0 = P0s
            else minp0[which(minp0 > P0s)] = P0s[which(minp0 > P0s)]
        }
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    Ts <- c(Ts2, min(pPerm0))
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


#' gene-Multitrait Sum of Powered Score (MTSPUsSet) tests and adaptive MTSPUsSet (MTaSPUsSet) test for multi trait - SNP set association with GWAS summary statistics.
#'
#' It gives p-values of the MTSPUsSet tests and MTaSPUsSet test with GWAS summary statistics.
#'
#' @param Zs Z-score matrix. row represent SNPs and column represent traits. It could be P-values if the Ps option is TRUE. 
#'
#' @param corSNP Correlation matirx of the SNPs to be tested; estimated from a
#' reference panel (based on the same set of the reference alleles as
#' used in calculating Z-scores).
#'
#' @param corPhe Correlation matirx of phenotypes to be tested; Estimated from Z-scores.
#' 
#' @param pow SNP specific power(gamma values) used in MTSPUsSet test.
#'
#' @param pow2 GENE specific power(gamma values) used in MTSPUsSet test.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @param Ps TRUE if input is p-value, FALSE if input is Z-scores. The default is FALSE.
#'
#' @param prune if it is TRUE, do pruing before the test using pruneSNP function. 
#'
#' @return A vector object, MTSPUsSet test P values and MTaSPUsSet P value.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#'
#' Il-Youp Kwak, Wei Pan (2016)
#' Gene- and pathway-based association tests for multiple
#'      traits with GWAS summary statistics, Bioinformatics, doi:10.1093/bioinformatics/btw577
#'
#' @examples
#'
#' data(SAMD11)
#' attach(SAMD11)
#' ## example analysis using MTaSPUsSet test.
#' (outFZ <- MTaSPUsSet(ZsF, corSNP=corSNPF, corPhe = corPheF,
#'       pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=10, Ps=FALSE))
#'
#' 
#' @seealso \code{\link{MTaSPUsSetC}} 


MTaSPUsSet <- function(Zs, corSNP, corPhe, pow=c(1,2,4,8),
                       pow2 = c(1,2,4,8), n.perm=5000, Ps = FALSE, prune = TRUE) {

    if( dim(corSNP)[1] != dim(corSNP)[2] ){
        stop(" corSNP should be correlation matrix with same number of rows and columns. ")
    }

    if( dim(corPhe)[1] != dim(corPhe)[2] ){
        stop(" corSNP should be correlation matrix with same number of rows and columns. ")
    }

    if( dim(Zs)[1] != ncol(corSNP) ) {
        stop(" The number of rows (SNPs) do not match with correlation matrix among SNPs (corSNP) ")
    }

    if( dim(Zs)[1] != ncol(corSNP) ) {
        stop(" The number of columns (traits) do not match with correlation matrix among traits (corPhe) ")
    }

    if(prune== TRUE) {
        pr <- pruneSNP(corSNP)
        if( length(pr$to.erase) > 0 ) {
            Zs <- as.matrix(Zs[-pr$to.erase,])
            corSNP <- corSNP[-pr$to.erase, -pr$to.erase]
        }
    }
    
    if( dim(Zs)[1] <= 1 ) {
        stop("less than 1 SNP.")
    }

    nsnp <- dim(Zs)[1]
    nphe <- dim(Zs)[2]
    
    V <- corSNP
    U <- corPhe
    
    e.U <- eigen(U)
    e.U$values[ e.U$values < 0 ] = 0
    A <- e.U$vectors %*% diag(sqrt(e.U$values)) %*% t(e.U$vectors)
    e.V <- eigen(V)
    e.V$values[ e.V$values < 0 ] = 0
    B <- e.V$vectors %*% diag(sqrt(e.V$values)) %*% t(e.V$vectors)
    
    Zs <- as.matrix(Zs)
    if(Ps == TRUE)
        Zs <- qnorm(1 - Zs/2)
    
    Ts <- rep(0, length(pow)*length(pow2))

    for(p1 in 1:length(pow)) {
        for(p2 in 1:length(pow2)) {

            if(pow[p1] < Inf) {
                Ts[p2 + (p1-1) * length(pow2)] = sum(apply(Zs, 2, function(x) sum(x^pow[p1]) )^pow2[p2])
            } else {
                Ts[p2 + (p1-1) * length(pow2)] = sum(apply(Zs, 2, max )^pow2[p2])
            }
        }
    }

####
## residual permutation
    pPerm0 = rep(NA,length(pow)*length(pow2))
                                        #    T0s1 = numeric(length(pow)*length(pow2) )
    T0s = numeric(n.perm)
    s <- sample(1:10^5,1)

    for (p1 in 1:length(pow)){
        for (p2 in 1:length(pow2)){
            set.seed(s) # to ensure the same samples are drawn for each pow
            for (b in 1:n.perm){
                Z <- matrix(rnorm(n=nsnp*nphe, 0,1), nphe, nsnp)
                Z0 <- t(A %*% Z %*% B)

                if(Ps == TRUE)
                    Z0 <- abs(Z0)

                if(pow[p1] < Inf) {
                    T0s[b] = sum(apply(Z0, 2, function(x) sum(x^pow[p1]) )^pow2[p2])
                } else {
                    T0s[b] = sum(apply(Z0, 2, max )^pow2[p2])
                }
            }

            pPerm0[p2 + (p1-1) * length(pow2)] = sum(abs(Ts[p2 + (p1-1) * length(pow2)])<=abs(T0s)) / n.perm
            P0s = ( (n.perm-rank(abs(T0s))) + 1 )/(n.perm)
            if (p2 + (p1-1) * length(pow2) == 1 ) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
        }

    }

    Paspu<-(sum(minp0<=min(pPerm0))+1)/(n.perm+1)
    pvs <- c(pPerm0, Paspu)

    nmvec <- NULL;
    for(nm in paste("MTSPUsSet",pow,",", sep=""))
        nmvec <- c(nmvec, paste(nm, pow2, sep="") )

    nmvec <- c(nmvec, "MTaSPUsSet")
    names(pvs) <- nmvec
    pvs
}

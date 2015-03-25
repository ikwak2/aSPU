
###extract the first fw PCs such that they explain at least cutoff*100%
### of total variation from X;
### Input: X: nObs by nSNPs;
###        cutoff: threshold to determine how many first few PCs to use;
### Output: a matrix consisting of the first few PCs;

extractPCs<-function(X, cutoff=0.95){

    if (is.null(ncol(X)) || ncol(X)<2) Xpcs<-X
    else{
        X.pca<-prcomp(X, center=FALSE)
        ##find the min k s.t. the first k PCs explain >= cutoff variations:
        Xevs<-X.pca$sdev^2
        for(k in 1:ncol(X))
            if (sum(Xevs[1:k])/sum(Xevs) >= cutoff) break;
        Xpcs<-X %*% X.pca$rot[,1:k]
    }

    Xpcs
}

#calculate a permuttaion p-value matrix based on a permuted stat matrix:
PermPvs<-function(T0s){
    B=nrow(T0s); n=ncol(T0s)
    P0s<-matrix(1, nrow=B, ncol=n)
    for(j in 1:n)
        for(b in 1:B)
            P0s[b,j] = sum( abs(T0s[b,j]) < abs(T0s[-b,j]) )/(B-1)
    return(P0s)
}

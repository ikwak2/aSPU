#' Image plot of P-value matrix.
#'
#' It gives a P-value image for a given P-value matrix.
#'
#' @param Ps P-value matrix. Row represent traits, column represent SNPs
#'
#' @param zlim -log 10 transformed p value range. 
#'
#' @param main main title.
#'
#' @param xlab lable for x axis. Default is "SNPs".
#'
#' @param thresh Default is set to widely used genome wise threshold -log(5e-8,10).
#'
#' @param trait.names A vector of trait names. 
#'
#' @param yt it controls the location of trait names. 
#'
#' @return Image of P-values matrix.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#' Il-Youp Kwak, Wei Pan (2017)
#' Gene- and pathway-based association tests for multiple traits with GWAS summary statistics, Bioinformatics. 33(1), 64-71
#'
#' @examples
#' ## Say we have 3 traits and their p-values at 5 SNPs. 
#' Ps <- rbind( c(0.001, 0.4, 0.5, 0.00000001, .1),
#'            c(0.03, 0.3, 0.3, 0.00001, .2),
#'            c(0.01, 0.2, 0.4, 0.001, .0001) )
#'
#' ## We can visualize it using plotPmat function.
#' plotPmat(Ps)
#' 
#' 
#' @seealso \code{\link{aSPUs}}

plotPmat <- function(Ps, zlim=NULL, main = NULL, yt = NULL, xlab = "SNPs", thresh=-log(5e-8,10), trait.names=NULL  ) {

    if(!is.matrix(Ps))
        stop("Ps need to be a matrix");

    if(is.null(trait.names))
        trait.names = paste("tr", 1:dim(Ps)[1], sep = "")

    log10P <- -log(Ps,10)
    
    if(is.null(zlim)) {
        maxP <- max(log10P, na.rm=TRUE)
        zlim <- c(-maxP, maxP)
    }

    if(zlim[1] != zlim[2]) {
        va <- max(abs(zlim[1]), abs(zlim[2]) )
        zlim = c(-va, va)
    }

    pos = 1:ncol(log10P)
    y = 1:nrow(log10P)
    log10P <- log10P
    log10P[log10P>thresh] <- log10P[log10P>thresh] * -1

    val <- sqrt(seq(0, 1, len=251))
    col <- c(rgb(val, val, 1), rgb(1, rev(val), rev(val))[-1])

    if(is.null(yt)) {
        yt = - 0.09655 * length(pos) + 0.4931034
    }

    par(mar=c(3.1, 3.6, 3.6, 0.6), las=1)
    layout(cbind(1, 2), widths=c(4.72, 1.28))
    image(pos, y, t(log10P), xaxt="n", yaxt="n",  ylab="", xlab="",
          zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", main = main)
    title(xlab=xlab, mgp=c(2, 0, 0))

    for(i in 1:dim(Ps)[1] ) {
       text(yt,i, trait.names[i], xpd = TRUE)
    }

    par(mar=c(4.1, 3.6, 3.6, 3.1))
    lodscale <- seq(0, zlim[2], length=length(col))

    ls2 <- lodscale
    ls2[ls2 > thresh] = ls2[ls2 > thresh]*-1
    par(mar=c(8.1, 0.6, 4.6, 4.1))
    image(0, lodscale, t(cbind(ls2)), xaxt="n", xlab="", yaxt="n",
          col=col, zlim=zlim)
    axis(side=4, at=pretty(lodscale), labels=abs(pretty(lodscale)))
#    title(main = expression(paste("-", log[10], "P") ))
    text(1, zlim[2]+1, expression(paste("-", log[10], "P") ), xpd=TRUE)
}


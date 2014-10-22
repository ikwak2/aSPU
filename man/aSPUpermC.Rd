\name{aSPUpermC}
\alias{aSPUpermC}
\title{Sum of powered score (SPU) test (perm, permutation part coded in C)}
\usage{
aSPUpermC(Y, X, cov = NULL, model = c("gaussian", "binomial"),
  pow = c(1:8, Inf), n.perm = 1000, userank = T)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases.  or It can be any quantitative
  traits.}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele.}

  \item{cov}{covariates}

  \item{model}{Use "gaussian" for quantitative trait
  (Default) , and Use "binomial" for binary trait.}

  \item{pow}{power used in SPU test.}

  \item{n.perm}{number of permutation}

  \item{userank}{similar to the original code if TRUE, by
  definition if FALSE}
}
\value{
Test Statistics and p-values for SPU tests and aSPU test.
}
\description{
It gives the p-values of the SPU test and aSPU test based
on based on the permutation of residuals.
}
\examples{
data(exdat)
out <- aSPUpermC(exdat$Y, exdat$X, cov = NULL, model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
out
}
\seealso{
\code{\link{aSPU}}, \code{\link{aSPUperm2}},
\code{\link{aSPUboot}}, \code{\link{aSPUboot2}}
}


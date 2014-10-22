\name{aSPUWC}
\alias{aSPUWC}
\title{Variance-weighted Sum of powered score (SPUw) test; using permutations to get the p-values (NO adjustment for covariates yet, C version).}
\usage{
aSPUWC(Y, X, pow = c(1:8, Inf), n.perm = 1000, rank = T)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases.  or It can be any quantitative
  traits.}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele.}

  \item{pow}{power used in SPU test.}

  \item{n.perm}{number of permutation}

  \item{userank}{similar code with the original version if
  TRUE. by definition if FALSE}
}
\value{
Test Statistics and p-values for SPU tests and aSPU test.
}
\description{
It gives the p-values of the SPUw test and aSPUw test based
on based on the permutation of residuals.  (NO adjustment
for covariates yet, C version)
}
\examples{
data(exdat)
out <- aSPUWC(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
out
}
\seealso{
\code{\link{aSPU}}, \code{\link{aSPUperm2}},
\code{\link{aSPUboot}}, \code{\link{aSPUboot2}}
}


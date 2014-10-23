\name{aSPUW}
\alias{aSPUW}
\title{Variance-weighted Sum of powered score (SPUw) test; using permutations to get the p-values (NO adjustment for covariates yet).}
\usage{
aSPUW(Y, X, pow = c(1:8, Inf), n.perm = 1000)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases.  or It can be any quantitative
  traits. Vector with length n (number of observations)}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele. Matrix with dimension n by g (n
  : number of observation, p : number of genotype data)}

  \item{pow}{power used in SPU test. Vector of g number of
  power.}

  \item{n.perm}{number of permutation}
}
\value{
Test Statistics and p-values for SPU tests and aSPU test.
}
\description{
It gives the p-values of the SPUw test and aSPUw test based
on based on the permutation of residuals.  (NO adjustment
for covariates yet)
}
\examples{
data(exdat)
out <- aSPUW(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
out
}
\seealso{
\code{\link{aSPU}}, \code{\link{aSPUperm2}},
\code{\link{aSPUboot}}, \code{\link{aSPUboot2}}
}


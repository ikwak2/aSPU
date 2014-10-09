\name{aSPU}
\alias{aSPU}
\title{Sum of powered score (SPU) test}
\usage{
aSPU(Y, X, cov = NULL, resample = c("perm", "boot"), model = c("gaussian",
  "binomial"), pow = c(1:8, Inf), n.perm = 1000, version = c("C", "mat",
  "vec"))
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases. or It can be any quantitative
  traits.}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele.}

  \item{cov}{covariates haha random stuff}

  \item{resample}{Use "perm" for residual permutations and
  "boot" for parametric bootstrap}

  \item{model}{Use "gaussian" for quantitative trait, and
  use "binomial" for binary trait.}

  \item{pow}{power used in SPU test.}

  \item{n.perm}{number of permutation}

  \item{version}{"vec" for the use of vector in
  permutation. "mat" for the use of matrix in permutation.
  Generally matrix version is faster but when n.perm is so
  big "mat" version does not work. This case we should use
  "vec" version. We have "C" version for residual
  permutation method. It coded in C and much faster than
  "mat" and "vec" method.}
}
\value{
Test Statistics and p-values for SPS tests and aSPU test.
}
\description{
It gives the p-values of the SPS test and aSPU test.
}
\examples{
data(exdat)
out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot", model = "binomial", pow = c(1:8, Inf), n.perm = 1000, version = "mat")
out
}
\seealso{
\code{\link{aSPUperm}}, \code{\link{aSPUperm2}},
\code{\link{aSPUboot}}, \code{\link{aSPUboot2}}
}


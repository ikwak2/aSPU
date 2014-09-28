\name{aSPU}
\alias{aSPU}
\title{Sum of powered score (SPU) test}
\usage{
aSPU(Y, X, cov = NULL, resample = c("boot", "perm"), pow = c(1:8, Inf),
  n.perm = 1000, version = 2:1)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases. or It can be any quantitative
  traits.}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele.}

  \item{cov}{covariates}

  \item{resample}{Use "boot" for parametric bootstrap and
  "perm" for residual permutations}

  \item{model}{Use "gaussian" for quantitative trait, and
  use "binomial" for binary trait.}

  \item{n.perm}{number of permutation}

  \item{version}{"1" for the use of vector T0. "2" for the
  use of matrix T0. Generally matrix version is faster but
  when n.perm is so big "2" does not work. This case we
  should use "1"}
}
\value{
Test Statistics and p-values for SPS tests and aSPU test.
}
\description{
It gives the p-values of the SPS test and aSPU test.
}
\examples{
data(exdat)
out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot", pow = c(1:8, Inf), n.perm = 1000, version = 2)
out
}
\seealso{
\code{\link{aSPUperm}}, \code{\link{aSPUperm2}},
\code{\link{aSPUboot}}, \code{\link{aSPUboot2}}
}


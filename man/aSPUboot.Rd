\name{aSPUboot}
\alias{aSPUboot}
\title{Sum of powered score (SPU) test (boot, version 1, vector used in permutation)}
\usage{
aSPUboot(Y, X, cov = NULL, model = (c("gaussian", "binary")), pow = c(1:8,
  Inf), n.perm = 1000)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases. or it can be any quantitative
  traits.}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele.}

  \item{cov}{covariates}

  \item{model}{Use "gaussian" for quantitative trait, and
  use "binomial" for binary trait.}

  \item{pow}{power used in SPU test.}

  \item{n.perm}{number of permutation}
}
\value{
Test Statistics and p-values for SPS tests and aSPU test.
}
\description{
It gives the p-values of the SPS test and aSPU test based
on the parametric bootstrap. (This is version 1, matrix
version is faster but if it doesn't work, we should use
version 1, vector version)
}
\examples{
data(exdat)
out <- aSPUboot(exdat$Y, exdat$X, cov = NULL, model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
out
}
\seealso{
\code{\link{aSPU}}, \code{\link{aSPUperm}},
\code{\link{aSPUperm2}}, \code{\link{aSPUboot2}}
}


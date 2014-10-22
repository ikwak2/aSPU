\name{aSPU}
\alias{aSPU}
\title{Sum of powered score (SPU) test}
\usage{
aSPU(Y, X, cov = NULL, resample = c("perm", "boot"), model = c("gaussian",
  "binomial"), pow = c(1:8, Inf), n.perm = 1000, userank = T)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases. or It can be any quantitative
  traits.}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele.}

  \item{cov}{covariates}

  \item{resample}{Use "perm" for residual permutations and
  "boot" for parametric bootstrap}

  \item{model}{Use "gaussian" for quantitative trait, and
  use "binomial" for binary trait.}

  \item{pow}{power used in SPU test.}

  \item{n.perm}{number of permutation}

  \item{userank}{use similar code with original version if
  T, by definition if F}
}
\value{
Test Statistics and p-values for SPU tests and aSPU test.
}
\description{
It gives the p-values of the SPU test and aSPU test.
}
\examples{
data(exdat)
out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot", model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
out
}
\seealso{
\code{\link{aSPUperm}}, \code{\link{aSPUperm2}},
\code{\link{aSPUboot}}, \code{\link{aSPUboot2}}
}


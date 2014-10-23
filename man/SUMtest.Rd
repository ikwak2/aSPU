\name{SUMtest}
\alias{SUMtest}
\title{SUM test (residual permutation, permutation part coded in C)}
\usage{
SUMtest(Y, X, cov = NULL, model = c("gaussian", "binomial"),
  n.perm = 1000)
}
\arguments{
  \item{Y}{phenotype data. It can be disease lables; =0 for
  controls, =1 for cases.  or It can be any quantitative
  traits. Vector with length n (number of observations)}

  \item{X}{genotype data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele. Matrix with dimension n by g (n
  : number of observation, p : number of genotype data)}

  \item{cov}{covariates. Matrix with dimension n by k (n
  :number of observation, k : number of covariates)}

  \item{model}{Use "gaussian" for quantitative trait
  (Default) , and Use "binomial" for binary trait.}

  \item{n.perm}{number of permutation}
}
\value{
SUM Test Statistic and p-value for SUM test.
}
\description{
It gives the p-values of the SUM test based on the
permutation of residuals.
}
\examples{
data(exdat)
out <- SUMtest(exdat$Y, exdat$X, cov = NULL, model = "binomial", n.perm = 1000)
out
}
\seealso{
\code{\link{aSPU}}
}


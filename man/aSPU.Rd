\name{aSPU}
\alias{aSPU}
\title{Adaptive Sum of powered score (SPU) tests. (SPU and aSPU)}
\usage{
aSPU(Y, X, cov = NULL, resample = c("perm", "sim", "boot"), 
   model = c("gaussian", "binomial"), pow = c(1:8, Inf), 
   n.perm = 1000, userank = T)
}
\arguments{
  \item{Y}{Response or phenotype data. It can be disease lables; =0 for
  controls, =1 for cases. or It can be any quantitative
  traits. Vector with length n (number of observations)}

  \item{X}{Genotype or other data; each row for a subject, and each
  column for an SNP. The value of each element is the # of
  the copies for an allele. Matrix with dimension n by g (n
  : number of observation, p : number of genotype data)}

  \item{cov}{covariates. Matrix with dimension n by k (n
  :number of observation, k : number of covariates)}

  \item{resample}{Use "perm" for residual permutations, "sim" for
    simulations from the distribution and "boot" for parametric
    bootstrap}

  \item{model}{Use "gaussian" for quantitative trait, and
  use "binomial" for binary trait.}

  \item{pow}{power used in SPU test. Vector of g number of
  power.}

  \item{n.perm}{number of permutation or bootstraps}

  \item{userank}{use similar code with original version if
  T, by definition if F}
}
\value{
Test Statistics and p-values for SPU tests and aSPU test.
}

\references{
  Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen, Peng Wei (2014) A powerful and adaptive association test for rare variants
  \emph{Genetics} 114.165035.
}


\description{
It gives the p-values of the SPU test and aSPU test.
}
\examples{
data(exdat)
out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot",
         model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
out
}
\seealso{
\code{\link{aSPUw}}
}


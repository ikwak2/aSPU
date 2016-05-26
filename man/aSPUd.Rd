\name{aSPUd}
\alias{aSPUd}
\title{Adaptive Sum of powered score (SPU) tests (SPU and aSPU) (distribution based).}
\usage{
aSPUd(Y, X, cov = NULL, model = c("gaussian", "binomial"))
}
\arguments{
\item{Y}{Response or phenotype data. It can be a disease indicator; =0 for controls, =1 for cases.
Or it can be a quantitative trait. A vector with length n (number of observations).}

\item{X}{Genotype or other data; each row for a subject, and each column
for an SNP (or a predictor). The value of each SNP is the # of the copies
for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ).}

\item{cov}{Covariates. A matrix with dimension n by p (n :number of observation, p : number of covariates).}

\item{model}{Use "gaussian" for a quantitative trait, and use "binomial" for a binary trait.}

}
\value{
p-values for SPU(1), SPU(2), minP tests and aSPU test.
}
\description{
It gives the p-values of the SPU(1), SPU(2) and minP tests and aSPU test based on the distribution.
}
\examples{

data(exdat)

## example analysis using aSPU test on exdat data.
out <- aSPUd(exdat$Y, exdat$X, cov = NULL, model = "binomial")

out
}
\author{
Il-Youp Kwak and Wei Pan
}
\references{
Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
A powerful and adaptive association test for rare variants,
Genetics, 197(4), 1081-95

Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014) Comparison of statistical tests for group differences in brain functional networks, NeuroImage, 1;101:681-694
}
\seealso{
\code{\link{aSPU}}
}


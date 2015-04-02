\name{getlogitp}
\alias{getlogitp}
\title{Get logistic p-value}
\usage{
getlogitp(Y, X)
}
\arguments{
  \item{Y}{Response or phenotype data. It can be disease
  lables; =0 for controls, =1 for cases.}

  \item{X}{Genotype or other data; each row for a subject,
  and each column for an SNP. The value of each element is
  the # of the copies for an allele. Matrix with dimension
  n by g (n : number of observation, p : number of genotype
  data)}
}
\value{
p-values for each SNPs.
}
\description{
Get p-value using logistic regresion for each SNPs
}
\examples{
simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
                           LOR=.2, rholim=c(0,0),
                           n=100, MAFlim=c(0.05, 0.4), p0=0.05)
logitp <- getlogitp(simula$Y, simula$X)
}
\seealso{
\code{\link{GATES2}} \code{\link{GatesSimes}}
\code{\link{Hyst}}
}


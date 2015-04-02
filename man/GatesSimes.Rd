\name{GatesSimes}
\alias{GatesSimes}
\title{GATES-Simes method}
\usage{
GatesSimes(pvec, ldmatrix, snp.info, gene.info)
}
\arguments{
  \item{pvec}{p-values for each SNP.}

  \item{ldmatrix}{numeric. A correlation matrix of SNPs,
  dimensions matching the p and snps arguments.}

  \item{snp.info}{SNP information matrix, The 1st column is
  SNP ids, 2nd column is SNP chromosome, 3rd column
  indicate SNP location.}

  \item{gene.info}{GENE information matrix, The 1st column
  is GENE ids, 2nd column is GENE chromosome, 3rd and 4th
  column indicate where the gene location starts and ends.}
}
\value{
p-value based on Gates-Simes method.
}
\description{
Get p-value using GATES-Simes method.
}
\examples{
simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
                           LOR=.2, rholim=c(0,0),
                           n=100, MAFlim=c(0.05, 0.4), p0=0.05)
logitp <- getlogitp(simula$Y, simula$X)

## get correlation of SNPs using controls
ldmat <- cor(simula$X[ simula$Y == 0, ])
out <- GatesSimes(pvec = logitp, ldmatrix = ldmat, snp.info = simula$snp.info,
                  gene.info = simula$gene.info)
}
\author{
Il-Youp Kwak
}
\references{
Hongsheng Gui, Miaoxin Li, Pak C Sham and Stacey S Cherny
(2011) Comparisons of seven algorithms for pathway analysis
using the WTCCC Crohn's Disease
}
\seealso{
\code{\link{Hyst}} \code{\link{GATES2}}
}


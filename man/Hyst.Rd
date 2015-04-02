\name{Hyst}
\alias{Hyst}
\title{HYST}
\usage{
Hyst(pvec, ldmatrix, snp.info, gene.info)
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
p-values for each SNPs.
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
out <- Hyst(pvec = logitp, ldmatrix = ldmat, snp.info = simula$snp.info,
            gene.info = simula$gene.info)
}
\author{
Il-Youp Kwak
}
\references{
Miao-Xin Li, Johnny S.H. Kwan and Pak C. Sham (2012) HYST:
A Hybrid Set-Based Test for Genome-wide Association
Studies, with Application to Protein-Protein
Interaction-Based Association Analysis
}
\seealso{
\code{\link{GatesSimes}} \code{\link{GATES2}}
}


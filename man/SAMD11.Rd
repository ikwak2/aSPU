\docType{data}
\name{SAMD11}
\alias{SAMD11}
\title{An example data for MTaSPUsSet test}
\usage{
data(kegg9)
}
\description{
  Genetic Investigation of ANthropometric Traits (GIANT) consortium
  data contain P-values of 2.7 million SNPs with six anthropometric
  traits that are well established to represent body size and shape:
  height, weight, BMI, waist circumference (WC), hip
  circumference(HIP). The data is downloadable from
  https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files

  We mapped SNPs on gene SAMD11. This subdata contains P values, Z
  scores of SNPs mapped on gene SAMD11. It also contains correlation
  among SNPs and correlation among phenotypes for demostrated example of
  MTaSPUsSet test. 
  
}
\examples{
data(SAMD11)

## Z-scores for man
SAMD11$ZsM

## P-values for man
SAMD11$PsM

## correlation among SNPs for man
SAMD11$corSNPM

## correlation among phenotypes for man
SAMD11$corPheM

}
\references{
Shungin D, Winkler TW, Croteau-Chonka DC, Ferreira T, Locke AE, Magi R,
Strawbridge R, Pers TH, Fischer K, Justice AE, Workalemahu T, Wu JM, et
al. (2015). New genetic loci link adipose and insulin biology to body
fat distribution. Nature 518: 187-196.
}


\keyword{datasets}

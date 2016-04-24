\docType{data}
\name{someGs}
\alias{someGs}
\title{An example data for MTaSPUsSet test}
\usage{
data(someGs)
}
\description{
  Genetic Investigation of ANthropometric Traits (GIANT) consortium
  data contain P-values of 2.7 million SNPs with six anthropometric
  traits that are well established to represent body size and shape:
  height, weight, BMI, waist circumference (WC), hip
  circumference(HIP). The data is downloadable from
  https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files

  We mapped SNPs on gene LCORL, RASA2, STK33 and RPGRIP1L. This subdata
  contains P values, correlation matrices for demostrated example of
  MTaSPUsSet test. 
  
}
\examples{
data(someGs)

## P-values for LCORL
someGs$LCORL[[1]]

## correlation matrix for LCORL
someGs$LCORL[[2]]

## P-values for RASA2
someGs$RASA2[[1]]

## correlation matrix for LCORL
someGs$RASA2[[2]]

# See vignettes for more details.

}
\references{
Shungin D, Winkler TW, Croteau-Chonka DC, Ferreira T, Locke AE, Magi R,
Strawbridge R, Pers TH, Fischer K, Justice AE, Workalemahu T, Wu JM, et
al. (2015). New genetic loci link adipose and insulin biology to body
fat distribution. Nature 518: 187-196.
}


\keyword{datasets}

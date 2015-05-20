\docType{data}
\name{exdat}
\alias{exdat}
\title{An example data set}
\usage{
data(exdat)
}
\description{
  The exdat data set is list of three objects. exdat$Y is a vector of length
  1000, 500 0s and 500 1s. exdat$X is a matrix of 1000 by 10. This
  simulated X matrix is assumed to be rare variants. All
  eliments are 0, 1 or 2. exdat$SNP0indx is a vector of length 10. The
  values are 1 or 0, 0 indicate corresponding column of X matrix have no
  association with Y vector.
}
\examples{
data(exdat)
exdat$X[c(1:10, 501:510),]
exdat$Y[c(1:10, 501:510)]

# 3rd and 10th column of X have no association with Y
exdat$SNP0indx

}
\keyword{datasets}


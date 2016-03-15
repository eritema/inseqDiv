\name{pmdSamp}
\alias{pmdSamp}

\title{
Function for resampling and calculating standard deviation, ploting of error plot for PMD index
}
\description{
The function can calculate the stablity of PMD index by calculating standard deviation by resamping the insertion site data. 
}
\usage{
pmdSamp(x, ntime, frac)
}

\arguments{
  \item{x}{
Insertion Sequence data matrix
}
  \item{ntime}{
number of times we want resampling e.g 1000
}
  \item{frac}{
fraction of sample we want resample e.g 0.75 
}
}
\details{
In order to calculate the standard deviation of PMD index, performed analysis by resampling Insertion site data. Resulting out put is in the form of matrix which contain Richness, Evenness and PMD index values and their respective standard deviation for defined fraction of sampling, n times of sampling. These informtions are calculated for every sample of matrix. 
}

\author{
shahzad 
}

\examples{
data(insertionSeq)
resultMat<-pmdSamp(insertionSeq,100,0.75)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

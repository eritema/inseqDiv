\name{expMatrix}
\alias{expMatrix}
\docType{data}
\title{
Experimental matrix for the insertion Sequence detection experiments for given insertion Sequence data. 
}
\description{
Experimental matrix contain information for various samples of experiment, their corresponding experimental factors. 
}
\usage{data(expMatrix)}
\format{
  A data frame with 17 samples of insertionSeq matrix on 4 experimental factors. Samples are taken on rows and experimental factors on columns. 
  
}
\details{
Data frame containing categorical and numerical experimental condition variables in columns for various samples on provided in rows. Sample names should be identical with sample names in Insertion Sequence data matrix and should be as row names. Experimental variables should be represented in columns as colmns names. 
}

\examples{
data(expMatrix)
}
\keyword{datasets}
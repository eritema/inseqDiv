\name{expDiversity}
\alias{expDiversity}

\title{
Function for calculating diversity of insertion sites as per various experimental factors
}
\description{
It provides the PMD index for the Insertion sites belong to various experimental factors in various samples. 
}
\usage{
expDiversity(expMat, factor, insertionSeq)
}

\arguments{
  \item{expMat}{
data frame matrix which is provided by user and contain experimental factors as columns. Experimental factors are either categorical or numerical variables. Experimental factors are taken as columns and sample names as row names. 
}
  \item{factor}{
is Experimental factor which is to be taken for calculating the pmd index of the samples belonging to variable experimental factor. 
}
  \item{insertionSeq}{
is insertion Sequence data frame matrix. It contains columns as samples while rows as insertion sites. 
}
}

\author{
Shahzad Ahmad 
}

\examples{
data(insertionSeq)
data(expMatrix)
expDiversity(expMatrix, "Enzyme", insertionSeq) 
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

\name{pmd}
\alias{pmd}

\title{
Function for calculating the PMD index of the insertion sequence data 
}
\description{
Provides sample PMD index, Evenness, Shannon, Simpson and Richness values for each of sample in data frame
}
\usage{
pmd(x,thr=1.1)
}

\arguments{
  \item{x}{
Insertion Sequence data matrix. x should be of type data frame or matrix. further detailed properties of insertion site data set are given in insertionSeq package data "insertionSeq" documentation.
}
}
\details{
In order to implement the PMD index for the insertion site data, we formulated pmd() function. It calculated the pmd index by taking account of evenness and richness of data samples. As we provide input data as mtrix of various samples, this returns the pmd index value, richness and evenness for each of sample in form of matrix.
}

\author{
shahzad
}

\examples{
data(insertionSeq)
pmd(insertionSeq)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

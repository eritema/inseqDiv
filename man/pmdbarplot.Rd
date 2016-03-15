\name{pmdbarplot}
\alias{pmdbarplot}

\title{
Function for creating barplot of PMD data of pmd() function of data Matrix. 
}
\description{
This function plots the barplot for results of pmd() function. It deals with pmd values of various samples as length of bar either towards clonal or monoclonal space of plot. 
}
\usage{
pmdbarplot(pmdres)
}

\arguments{
  \item{pmdres}{
  result matrix output by pmd() function. Argument should be of type matrix and specific formate as of output my pmd( ) function 
}
}
\details{
Representation of PMD index is important to show the variation of clonality among various samples. we have used barplot for the visualization of pmd index results for various samples. we are representing the samples being more clonal as green bars as containing pmd index value greater than zero while sample shown to be clonal express less pmd index value on scale and shown as red bars.  
}

\author{
shahzad ahmad
}

\examples{
data(insertionSeq)
pmdresult<-pmd(insertionSeq)
pmdbarplot(pmdresult)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
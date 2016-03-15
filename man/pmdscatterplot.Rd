\name{pmdscatterplot}
\alias{pmdscatterplot}

\title{
Function for ploting results of function pmd() in scatter plot 
}
\description{
Provides a scatter plot for result Matrix of PMD function. It plots the Richness of component of matrix on x exis and Evenness component on Y axix. Each dot on plot represent a sample in experiment and it can also be labelled with its name. Color of plot is given according to the value of PMD index. Red color indicated sample being more monoclonal and blue being more polyclonal.  
}
\usage{
pmdscatterplot(pmdres, nameofplot)
}

\arguments{
  \item{pmdres}{
is the Matrix provided by function pmd() 
}
  \item{nameofplot}{
argument is a character type argument proving name of the resulting plot.
}
}
\details{
Provide scatter plot for the pmd function results 
}

\author{
shahzad Ahmad
}

\examples{
data(insertionSeq)
pmdRes<-pmd(insertionSeq)
pmdscatterplot(pmdRes,"test Data 1")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

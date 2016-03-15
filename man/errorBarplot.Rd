\name{errorBarplot}
\alias{errorBarplot}

\title{
Function for ploting the results of pmdsamp function. 
}
\description{
Function for ploting the results of pmdsamp function.
}
\usage{
errorBarplot(X, labelling = TRUE)
}

\arguments{
  \item{X}{
resulting result matrix of pmdSamp
}
  \item{labelling}{
String name of the graph
}
}
\details{
it describes the error defined by the plots 
}


\author{
shahzad
}

\examples{
data(insertionSeq)
pmdsampResult<-pmdSamp(insertionSeq,100,0.75)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

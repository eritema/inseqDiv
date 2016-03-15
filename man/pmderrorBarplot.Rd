\name{pmderrorBarplot}
\alias{pmderrorBarplot}

\title{
Function creats error barplots for the pmdSamp results 
}
\description{
Function creats error barplots for the pmdSamp results 
}
\usage{
pmderrorBarplot(pmdres)
}

\arguments{
  \item{pmdres}{
pmd results of pmdSamp function 
}
}
\details{
provides standard deviation in the form of error barplot 
}

\examples{
pmdres<-pmdSamp(insertionSeq,100,0.7)
pmderrorBarplot(pmdres)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

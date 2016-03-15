\name{reconstruct}
\alias{reconstruct}

\title{
Supportive function for the pmdSamp
}
\description{
Supportive function for the pmdSamp function to convert the integration site data into base sampling data 
}
\usage{
reconstruct(abun.vec)
}

\arguments{
  \item{abun.vec}{
vector if untegration sites called by function pmdSamp
}
}
\details{
pmdSamp requires a function to convert the 
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (abun.vec) 
{
    len.abun <- length(abun.vec)
    numVec <- c(1:len.abun)
    finraw <- numeric()
    for (i in 1:length(numVec)) {
        temp <- rep(numVec[i], abun.vec[i])
        finraw <- c(finraw, temp)
    }
    return(finraw)
  }
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

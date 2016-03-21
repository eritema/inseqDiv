\name{pmdscatterplot}
\alias{pmdscatterplot}

\title{
Function for ploting results of function pmd() in scatter plot 
}
\description{
Provides a scatter plot on the clonal plane for the samples analyzed via pmd(). It plots the Richness on x axis and Evenness component on Y axix. Each dot on plot represent a experimental sample and can be labelled with its name. Color of plot is given according to the value of PMD index. Red color indicated sample being more monoclonal and blue being more polyclonal.  
}
\usage{
pmdscatterplot(pmdres, nameofplot="",labelling=T,labelsize=3,spotRadius=7)
}

\arguments{
  \item{pmdres}{
    is the Matrix provided by function pmd() 
  }
  \item{nameofplot}{
    title of the plot
  }
  \item{labelling}{
    flag that controls the presence/absence of sample labels in the plot
  }
  \item{labelsize}{
    dimension of the labels
  }
  \item{spotRadius}{
    dimension of the data spot
  }
}
\details{
Provide a scatter plot for the pmd function results.  
}

\author{
Shahzad Ahmad
}

\examples{
data(insertionSeq)
pmdRes<-pmd(insertionSeq)
pmdscatterplot(pmdRes,"test Data 1")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

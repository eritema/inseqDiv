\name{heatmapMatrix}
\alias{heatmapMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to creat heatmap data matrix
}
\description{
This function creats a data matrix ultimately used for the creation of heatmap. It takes the Insertion Sequence data as input and output the resulting matrix containing the clonality Values for the desired number most clonal Insertion Sites.
}
\usage{
heatmapMatrix(data, n)
}

\arguments{
  \item{data}{
Insertion Sequence data matrix
}
  \item{n}{
Threshold value of numeric type for most clonal Insertion Sites in Insertion Site data e.g 20, 30
}
}
\details{
For following the clonality of insertion Sites and Samples during various sampling time periods during gene therapy clinical trials, we choose a heatmap method. for producing a heatmap we requires a data matrix from which we can creat the heatmap. In this function we input the insetion sequence data, threshold for most clonal genes and output the desired heatmap matrix.

}

\examples{
data(insertionSeq)
heatmapMatrix<-heatmapMatrix(insertionSeq,10)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

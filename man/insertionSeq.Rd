\name{insertionSeq}
\alias{insertionSeq}
\docType{data}
\title{
Insertion Sequence data matrix of one patient for various samples. 
}
\description{
This matrix contain detected insertion sites at various sampling times for one patient 
}
\usage{data(insertionSeq)}
\format{
  A data frame with 12041 insertion sites in rows on the 17 time periods (samples) as columns.  
  
}
\details{
Insertion site matrix contain data frame which represent insertion sites as rows of matrix and columns as samples taken at various time periods or various conditions. Insertion sites should be represented as integer value while empty spaces in columns as zero. row names whould be insertion site names while columns names as various time periods or sample names. 
}

\references{
cartier, Nathalie et.al.(2009).Hematopoietic Stem Cell Gene Therapy with a Lentiviral Vector in X-Linked Adrenoleukodystrophy.Science,Vol. 326 no. 5954 pp. 818-823.
}
\examples{
data(insertionSeq)
}
\keyword{datasets}

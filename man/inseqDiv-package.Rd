\name{inseqDiv-package}
\alias{inseqDiv-package}
\alias{inseqDiv}
\docType{package}
\title{
Insertion Sequence Diversity analysis. 
}
\description{
inseqDiv include functions that can perform diversity analysis of insertion sites in cancer Gene therapy patients. we focuses on the qutitative isolation of the mono, oligo and polyclonality of the gene corrected clones by using new index PMD index. It also include funtion for determining the effect of experimental factors on the diversity of transgenic clones.
}
\details{
\tabular{ll}{
Package: \tab inseqDiv\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-04\cr
License: \tab GPL-3\cr
}
Using our package you can calculate pmd index of samples and various visualization schemes to view the status of samples in large set of experiments. We can calculate PMD index by pmd() function. Results of this function can be translated into various plots. 
We also need to follow the clonality of insertion sites during various time periods in Gene therapy trials. we have developed a function which could give a mtrix for creation of heatmap of clonality results. 
}
\author{
Fronza, Raffaele and Ahmad, shahzad

Maintainer:Shahzad Ahmad
<shahzad.gcu@hotmail.com>
}

\keyword{ package }

\examples{
data(insertionSeq)
pmdresult<-pmd(insertionSeq)
barplot(pmdresult)
}

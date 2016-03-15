reduceMatrix <- function(sorta,n=2) {
  j=numeric(0)
  #test<-matrix(0,nrow=dim(sorta)[1],ncol=dim(sorta)[2])
  for (i in 1:dim(sorta)[1]) {
    tab<-table(sorta[i,])
    if (!length(tab[names(tab)==-100])) {
      j<-c(j,i)
    }
    else if(!(dim(sorta)[2]-(tab[names(tab)==-100])>n-1)){
      j<-c(j,i)
    }
  }  
  return(sorta[-j,])
}
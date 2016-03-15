simpson<-function(x,index="invsimps") {
  simpson<-numeric(ncol(x))
  for(i in 1:dim(x)[2]) {
    freq<-prop.table(t(x)[i,])
    if(index=="simpson")
      simpson[i]<-sum(freq*freq)
    else 
      simpson[i]<-1/sum(freq*freq)
  }
  return(simpson)
}
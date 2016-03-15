shannon<-function(x,base=2) {
  shann<-numeric(ncol(x))
  for(i in 1:dim(x)[2]) {
    freq<-prop.table(t(x)[i,])
    shann[i]<--sum(freq[freq>0]*log(freq[freq>0],base))
  }
  return(shann)
}
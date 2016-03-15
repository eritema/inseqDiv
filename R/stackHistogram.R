stackHistogram<-function(newdata) {
  for(i in 1:dim(newdata)[2]) {
    if(prec!=newdata[i,5]) {
      prec=newdata[i,5]
      j=0
      print(prec)
    }
    else {
      if(j<10){
        print(prec)
      }
      j=j+1}
  }
}
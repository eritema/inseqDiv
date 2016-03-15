reconstruct <-
function(abun.vec){
    
    len.abun<-length(abun.vec)
    numVec<-c(1:len.abun)
    finraw<-numeric()
    
    for(i in 1:length(numVec)){
        
        temp<-rep(numVec[i],abun.vec[i])
        finraw<-c(finraw,temp)
        
    }
    
    return(finraw)
}

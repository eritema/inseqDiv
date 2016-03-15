experimentalDiversity<-function(data,ISfactor=c(1,2,4),factorList=c(14,13,12,11),thr=3){
  ### collapse code 
  CollapseData<-dataTransformation(data,ISfactor,factorList,thr)
  pmdres<-pmd(CollapseData)
  ## combining factors required in experimental matrix
  factor_vec<-row.names(pmdres)
  ## splitting the factors according to delimeter into list
  expFactor<-strsplit(factor_vec, "_")
  ## converting list into matrix and then dataframe 
  coln<-length(factorList)
  matexpFactor<- matrix(unlist(expFactor), ncol =coln, byrow = TRUE)
  datafram.expfac<-as.data.frame(matexpFactor)
  colnames(datafram.expfac)<-colnames(data[factorList])
  ###removing duplicates by first row 
  dataFram<-cbind(datafram.expfac,pmdres)
  #colnames(datafram.expfac)<-colnames(data[factorList])
  row.names(dataFram)<-factor_vec
  #datafram.expfac<-datafram.expfac[,-1]
  return(dataFram)
}
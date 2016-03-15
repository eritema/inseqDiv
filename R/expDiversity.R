expDiversity <-
function(expMat,factor,insertionSeq,boxplot=TRUE){
    
   variable <- as.vector(unique(expMat[,factor]))
   
    finalMat<-matrix()
 
    for(i in 1:length(variable)){
        
        temp<-names(which(apply(expMat==variable[i], 1, any)))
	
        mat<-insertionSeq[,temp]
	
        pmdVec<-pmd(mat)
	t<-dim(pmdVec)
	category<-rep(variable[i],t[1])
	
	pmdMat<-as.data.frame(cbind(pmdVec,category))

	if(i==1){
	finalMat<-pmdMat
	}
	else{

	finalMat<-rbind(finalMat,pmdMat)

          }
    }
	finalMat<-transform(finalMat, pmdIndex = as.numeric(as.character((pmdIndex))))
	finalMat<-transform(finalMat, Richness = as.numeric(as.character((Richness))))
	finalMat<-transform(finalMat, Evenness = as.numeric(as.character((Evenness))))
	print(summary(aov(pmdIndex~category, finalMat)))
if(boxplot==TRUE){
        return(finalMat)
	}else {
finalMat<-aggregate(.~category, data=finalMat, mean)
row.names(finalMat)<-finalMat[,1]
finalMat<-finalMat[,-1]
return(finalMat)
}
}
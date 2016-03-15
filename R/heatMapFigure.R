reduceMatrix <- function(sorta,n=2) {
  j=numeric(0)
  #test<-matrix(0,nrow=dim(sorta)[1],ncol=dim(sorta)[2])
  for (i in 1:dim(sorta)[1]) {
    tab<-table(sorta[i,])
    if(!is.null(tab[names(tab)==-100]) && !(dim(sorta)[2]-(tab[names(tab)==-100])>n-1)){
      j<-c(j,i)
    }
  } 
  
  return(sorta[-j,])
}

heatMapFigure<-function (data, rowC=F, colC=F, time=F, ncom=1, g100=F, autoRange=T, minimo=-6, key=F, cexRowValue=0,cexColValue=0) {
  print(ncom)
	library(gplots)
	a<-data.matrix(data)
	if(ncom>1) {
		a<-reduceMatrix(a,ncom)
    print(head(a))
	}
	if(autoRange) {
	  minimo<-as.integer(names(table(a)[2]))*2
    if ((minimo<=0)&&(minimo>-1)) minimo=-1
    else if ((minimo>-4)&&(minimo<-1)) minimo=-4
	}
		
	if(g100) {
		g100.palette <- colorRampPalette(c("white","red","orange","yellow","green","blue","gray"),space = "Lab", bias=0.4)
		minRange<-minimo
		maxRange<-0
	}
	else {
		g100.palette <- colorRampPalette(c("white","red","gray","blue"),space = "Lab", bias=1.7)
		minRange<-minimo*2
		maxRange<-abs(minRange)
	}
	if(time) {
		sorta<-a[do.call(order,as.data.frame(a)),]
		a<-sorta
		tree="none"
	}
	else if ((!rowC)&&(colC)) 
		tree="col"
	else if ((rowC)&&(!colC)) 
		tree="row"
	else if ((rowC)&&(colC))
		tree="both"
	else
		tree="none"
  
  if(key==T) 
    keysizeValue=0.85
  else keysizeValue=0.2
  
  if ((cexRowValue==0)&&(cexColValue==0))
	  heatmap.2(a,Rowv=rowC,Colv=colC,col = g100.palette(256),dendrogram=tree,trace="none",density.info="none",symbreaks=T,breaks=seq(minRange,maxRange,(abs(minRange)+maxRange)/256),key=key,keysize=keysizeValue)  # cexRow=0.2,cexCol=0.2, lhei=c(0.1,1),lwid=c(0.1,1))
  else if ((cexRowValue!=0)&&(cexColValue==0))
    heatmap.2(a,Rowv=rowC,Colv=colC,col = g100.palette(256),dendrogram=tree,trace="none",density.info="none",symbreaks=T,breaks=seq(minRange,maxRange,(abs(minRange)+maxRange)/256),key=key, keysize=keysizeValue,cexRow=cexRowValue)
	else if ((cexRowValue==0)&&(cexColValue!=0))
	  heatmap.2(a,Rowv=rowC,Colv=colC,col = g100.palette(256),dendrogram=tree,trace="none",density.info="none",symbreaks=T,breaks=seq(minRange,maxRange,(abs(minRange)+maxRange)/256),key=key, keysize=keysizeValue,cexCol=cexColValue)
  else 
    heatmap.2(a,Rowv=rowC,Colv=colC,col = g100.palette(256),dendrogram=tree,trace="none",density.info="none",symbreaks=T,breaks=seq(minRange,maxRange,(abs(minRange)+maxRange)/256),key=key, keysize=keysizeValue,cexRow=cexRowValue,cexCol=cexColValue)
}

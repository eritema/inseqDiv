## data= data matrix from HISAP where first three colums should be
## 1= Chromosome 2=IntegrationLocus 3=Number of Reads. Other columns can be any experimental factors
## Function "dataTransformation" By default use 1,2,4 columns as Insertion site name. (Here i have also defined default for experimental)
## any range can be provided 1 to 6 in experimental factors
## Final Matrix can be used in any function of inseqDiv package

dataTransformation<-function(data,ISrange=c(1,2,4),factorRange=c(12,13,14),thr=3,sub,fraction){
    # Sort data by (chromosome,position)
    if (sub==TRUE) {
      data<-subSample(data,fraction)
      print(length(data[,3]<3))
    }
    sortdata<-data[order(as.numeric(data[,1]),as.numeric(data[,2]),na.last=NA),]
    prec<--10
    cycle<-1
    cycDisp<-0

    # Collapse IS that are closer than (thr) bases
    print("Starting collapsing")
    #mat<-data.matrix(sortdata[,1:3])

    mat<-sortdata[,1:3]
    #print(dim(mat))
    #print(head(mat,3))
    #for (i in 1:dim(sortdata)[1]) {
    for (i in 1:nrow(mat)) {
        if(cycle>=1000) {cycDisp<-cycle+cycDisp;print(cycDisp);cycle=0}   # Job advancement
        cycle<-cycle+1
        # if the location is within thr (3) bases change to the previous
        position<-mat[i,2]
        #if (abs(sortdata[i,2]-prec)<=thr) {prec<-sortdata[i,2];sortdata[i,2]<-first;}
        if (abs(position-prec)<=thr) {prec<-position;mat[i,2]<-first;}
          # else is a new IS
        #else {first<-sortdata[i,2];prec<-first}
        else {first<-position;prec<-first}
    }
    print("End collapsing")

    # substitute the original IS list with the collapsed
    sortdata[,1:3]<-as.data.frame(mat)
    CollapseData<-sortdata
    #print(dim(CollapseData))
    #print(head(CollapseData,3))
    IS_vec<-do.call(paste, c(CollapseData[ISrange], sep = "_"))
    factor_vec<-do.call(paste, c(CollapseData[factorRange], sep = "_"))
    numberOfread_vec<-CollapseData[,3]
    dataMatrix<-cbind(factor_vec,IS_vec,numberOfread_vec)
    #print(dim(dataMatrix))
    #print(head(dataMatrix,3))
    dataMatrix<-transform(dataMatrix, numberOfread_vec = as.numeric(as.character((numberOfread_vec))))

    dataMatrix<-as.data.frame(dataMatrix)
    

    abundance <- as.numeric(dataMatrix[,3])

    temp<- data.frame(rep(NA, sum(abundance)))
    print("start Data Construction")
    for (i in 1:(ncol(dataMatrix)-1)) {
            temp[,i] <- rep(as.character(dataMatrix[,i]), abundance)
        }
    print("End Data Construction")
    #print(dim(temp))
    #print(head(temp))
    finalMat<- as.data.frame.matrix(table(temp))
    finalMat<-t(finalMat)
    #print(head(finalMat))
    return(finalMat)
}

subSample<-function(data,fraction) {
  print("start data subsampling:")
  print(dim(data))
  matApp<-data[rep(seq_len(nrow(data)), data[,3]), ]
  matApp[,3]<-1
  #print("Dimension before Subsampling:")
  #print(dim(matApp)[1])
  matApp<-matApp[sample(nrow(matApp),as.integer(sum(matApp$V3)*(1-fraction)),replace=T ), ]
#  index<-sample(1:nrow(matApp),as.integer(sum(matApp$V3)*fraction),replace=T)
#  matApp<-matApp[-index,]
  #print("Dimension after Subsampling:")
  #print(dim(matApp)[1])
  matApp<-aggregate(V3~., data=matApp, FUN=sum)
  
  # column 3 should be the ones that contains readsNumber
  app<-matApp[,ncol(matApp)]
  matApp[,ncol(matApp)]<-matApp[,3]
  matApp[,3]<-app
  
  print(dim(matApp))
  print("End data subsampling")
  return(matApp)
}

pmd <-
function(data,ISrange=c(1,2,4),factorRange=c(12,13,14),thr=3,threshold = 1.1,noNaN=TRUE,chao=FALSE,sub=FALSE,fraction=0.5){
  
  
  pmd<-list(data=data,ISrange=ISrange,factorRange=factorRange,thr=thr,threshold=threshold,noNaN=noNaN) #build the object
    class(pmd) <- "pmd"
    x<-dataTransformation(pmd$data,pmd$ISrange,pmd$factorRange,pmd$thr,sub,fraction)
    speciesData2<-t(x)
    pmd$speciesData<-speciesData2
    ## calculation of the richness of species data
    pmd$shann<-shannon(x)
    pmd$simpson<-simpson(x)
    # if chao=FALSE then the richnes is computed without correction
    if(!chao) {
      for(i in 1:nrow(speciesData2)){
        tem<-0
        for(j in 1:ncol(speciesData2)){
          if(speciesData2[i,j]>0){
            tem=tem+1
          }
          pmd$richness[i]<-log(tem)
        }
      }
    # if chau=TRUE we use chao as abundance richness estimator 
    # if the number of specie is low is possible that the chao estimator return NaN or Inf
    } else {
      #print(nrow(speciesData2))
      #print(head(speciesData2),1)
      for(i in 1:nrow(speciesData2)){
        tem<-0
        singleton<-1
        doubleton<-1
        for(j in 1:ncol(speciesData2)){
          if(speciesData2[i,j]>0){
            tem=tem+1
          }
          if (speciesData2[i,j]==1){
            singleton=singleton+1
          } else if (speciesData2[i,j]==2){
            doubleton=doubleton+1
          }
        }
        #print(tem)
        pmd$richness[i]<-log(tem+(singleton*singleton/(2*doubleton)))
        #print (c(tem,singleton,doubleton,pmd$richness[i]))
      }
    }
    ## calculation of eveness of the species
    abun.vec<-rowSums(speciesData2)
    maxrow<-apply(speciesData2,1,max)
    vec.eve<-numeric(0)
    for(i in 1:length(maxrow)){
        evenness<--log(maxrow[i]/abun.vec[i])
        vec.eve[i]<-evenness
    }
    ## calcultion of the pmd index of the species
    pmd.index<-numeric(0)
    for(i in 1:length(pmd$richness)){
        dist<-pmd$richness[i]-vec.eve[i]
        pmd.sup<-log2(vec.eve[i]/dist)
        pmd.index<-c(pmd.index,pmd.sup)
    }
    pmd$pmd<-pmd.index
    pmd$eveness<-vec.eve
    ## Creation of Matrix
    eulerMat <- matrix(c(pmd$richness,pmd$shann,pmd$simpson,pmd$eveness,pmd$pmd), nrow = nrow(speciesData2), ncol = 5,dimnames=list(row.names(speciesData2), c("Richness","Shannon","Simpson","Evenness","pmdIndex")))
    ## removal of NaN and Inf Value
    #eulerMat<-eulerMat[!rowSums(!is.finite(eulerMat)),]
    print(eulerMat)
    if (noNaN){
      #eulerMat<-eulerMat[!rowSums(is.nan(eulerMat)),]
      eulerMat<-eulerMat[!rowSums(!is.finite(eulerMat)),]
      eulerMat<-eulerMat[eulerMat[,1]> threshold,]
    }
    ## removal of row having richness less than adjusted threshold
      #eulerMat<-eulerMat[eulerMat[,1]> threshold,]
    pmd$eulerMat<-data.frame(eulerMat)
    return(pmd)
}

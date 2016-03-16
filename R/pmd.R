## data= data matrix from HISAP where first three colums should be
## 1= Chromosome 2=IntegrationLocus 3=Number of Reads. Other columns can be any experimental factors
## Function "dataTransformation" By default use 1,2,4 columns as Insertion site name. (Here i have also defined default for experimental)
## any range can be provided 1 to 6 in experimental factors
## Final Matrix can be used in any function of inseqDiv package

dataTransformation<-function(data,ISrange=c(1,2,4),factorRange=c(12,13,14),thr=3){
    # Sort data by (chromosome,position)
    sortdata<-data[order(as.numeric(data[,1]),as.numeric(data[,2]),na.last=NA),]

    prec<--10
    cycle<-1
    cycDisp<-0

    # Collapse IS that are closer than (thr) bases
    print("Starting collapsing")
    #mat<-data.matrix(sortdata[,1:3])
    mat<-sortdata[,1:3]
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

    IS_vec<-do.call(paste, c(CollapseData[ISrange], sep = "_"))
    factor_vec<-do.call(paste, c(CollapseData[factorRange], sep = "_"))
    numberOfread_vec<-CollapseData[,3]
    dataMatrix<-cbind(factor_vec,IS_vec,numberOfread_vec)

    dataMatrix<-transform(dataMatrix, numberOfread_vec = as.numeric(as.character((numberOfread_vec))))

    dataMatrix<-as.data.frame(dataMatrix)

    abundance <- as.numeric(dataMatrix[,3])

    temp<- data.frame(rep(NA, sum(abundance)))
    print("start Data Construction")
    for (i in 1:(ncol(dataMatrix)-1)) {
            temp[,i] <- rep(as.character(dataMatrix[,i]), abundance)
        }
    print("End Data Construction")
    finalMat<- table(temp)
    finalMat<-t(as.data.frame.matrix(finalMat))
    return(finalMat)
}

pmd <-
function(data,ISrange=c(1,2,4),factorRange=c(12,13,14),thr=3,threshold = 1.1,noNaN=TRUE){
    pmd<-list(data=data,ISrange=ISrange,factorRange=factorRange,thr=thr,threshold=threshold,noNaN=noNaN) #build the object
    class(pmd) <- "pmd"
    x<-dataTransformation(pmd$data,pmd$ISrange,pmd$factorRange,pmd$thr)
    speciesData2<-t(x)
    ## calculation of the richness of species data
    pmd$shann<-shannon(x)
    pmd$simpson<-simpson(x)
    pmd$richness<-numeric(nrow(speciesData2))
    for(i in 1:nrow(speciesData2)){
        tem<-0
        for(j in 1:ncol(speciesData2)){
            if(speciesData2[i,j]!=0){
                tem=tem+1
            }
        }
        pmd$richness[i]<-log(tem)
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

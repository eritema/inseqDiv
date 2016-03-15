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

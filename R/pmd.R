## data= data matrix from HISAP where first three colums should be
## 1= Chromosome 2=IntegrationLocus 3=Number of Reads. Other columns can be any experimental factors
## Function "dataTransformation" By default use 1,2,4 columns as Insertion site name. (Here i have also defined default for experimental)
## any range can be provided 1 to 6 in experimental factors
## Final Matrix can be used in any function of inseqDiv package

dataTransformation<-function(data,ISrange=c(1,2,4),factorRange=c(12,13,14),thr=3,sub,fraction){
    
  
    # Check if the user asked for subsampling
    if (sub==TRUE) {
      data<-subSample(data,fraction)
      print(length(data[,3]<3))
    }
  
  
  # Sort data by (chromosome,position)
    sortdata<-data[order(as.numeric(data[,1]),as.numeric(data[,2]),na.last=NA),]

    # Set some flag variables
    prec<--10
    cycle<-1
    cycDisp<-0

    # Collapse IS that are closer than (thr) bases
    print("Starting collapsing")
    mat<-sortdata[,1:3]
    # for each integration in the IS table
    for (i in 1:nrow(mat)) {
      # check if the progression bar should be updated every 10000 cycles
        if(cycle>=1000) {cycDisp<-cycle+cycDisp;print(cycDisp);cycle=0}   # Job advancement
        cycle<-cycle+1
        # if the location is within thr (3) bases change the current position to the first
        position<-mat[i,2]
        if (abs(position-prec)<=thr) {prec<-position;mat[i,2]<-first;}
        # else is a new IS and the position is correct
        else {first<-position;prec<-first}
    }
    print("End collapsing")

    # substitute the original IS list with the collapsed
    sortdata[,1:3]<-as.data.frame(mat)
    CollapseData<-sortdata
    
    # construct names with colums in ISrange
    IS_vec<-do.call(paste, c(CollapseData[ISrange], sep = "_"))
    
    # consruct factor names with colums in factorRange
    factor_vec<-do.call(paste, c(CollapseData[factorRange], sep = "_"))
    
    # Extract the reads count
    numberOfread_vec<-CollapseData[,3]
    
    # construct the datamatrix
    dataMatrix<-cbind(factor_vec,IS_vec,numberOfread_vec)
    
    # to refactorize
    dataMatrix<-transform(dataMatrix, numberOfread_vec = as.numeric(as.character((numberOfread_vec))))
    
    # transform into dataframe (redundant with the operation before?)
    dataMatrix<-as.data.frame(dataMatrix)
    
    # abundance == numberOfread_vec. Refactorize
    abundance <- as.numeric(dataMatrix[,3])
    
    # initiaize a dataframe with as many rows as the total number of reads
    temp<- data.frame(rep(NA, sum(abundance)))
    print("start Data Construction")
    
    # for each element repeat its name as many time as is his abbundance (each name is a level)
    for (i in 1:(ncol(dataMatrix)-1)) {
            temp[,i] <- rep(as.character(dataMatrix[,i]), abundance)
        }
    print("End Data Construction")
    
    # Construct a contingency table using levels
    finalMat<- as.data.frame.matrix(table(temp))
    finalMat<-t(finalMat)
    return(finalMat)
}


# Subsample the original IS table (data) removing randomly (uniformly) n*fraction
# of the reads
subSample<-function(data,fraction) {
  print("start data subsampling:")
  
  # Replicate each row #nreads times
  # seq_len(num) return a regular sequence from 1 to num
  # rep(seq[],times[]) replicate the element seq[i] times[i] times
  # Combining the two return an index that is used to replicate the rows of mapApp
  matApp<-data[rep(seq_len(nrow(data)), data[,3]), ]
  
  # Each row has to count 1 in clonality
  matApp[,3]<-1
  
  # sample the rows of matApp 1-fraction times
  matApp<-matApp[sample(nrow(matApp),as.integer(sum(matApp$V3)*(1-fraction)),replace=T ), ]

  # aggregate split the matrix (matApp) on the basis of the factors on the right side of the formula (all: .)
  # operate the numerical function (sum) of the factor(s) on the left side (V3)
  matApp<-aggregate(V3~., data=matApp, FUN=sum)
  
  # column 3 should be the ones that contains readsNumber
  app<-matApp[,ncol(matApp)]
  matApp[,ncol(matApp)]<-matApp[,3]
  matApp[,3]<-app
  
  print("End data subsampling")
  return(matApp)
}

pmd <-
function(data,ISrange=c(1,2,4),factorRange=c(12,13,14),thr=3,threshold = 1.1,noNaN=TRUE,chao=FALSE,sub=FALSE,fraction=0.5){
  
  # define the arguments
  pmd <-
    list(
      data = data,ISrange = ISrange,factorRange = factorRange,thr = thr,threshold =
        threshold,noNaN = noNaN
    ) #build the object
  # instantiate the class
  class(pmd) <- "pmd"
  
  # Reformat the IS table into a dataframe that contains as rows the
  # unique IS and as column the factors (see ISOT@DKFZ.G100:MS)
  x <-
    dataTransformation(pmd$data,pmd$ISrange,pmd$factorRange,pmd$thr,sub,fraction)
  speciesData2 <- t(x)
  pmd$speciesData <- speciesData2
  ## calculation of the richness of species data
  pmd$shann <- shannon(x)
  pmd$simpson <- simpson(x)
  # if chao=FALSE then the richnes is computed without correction
  if (!chao) {
    for (i in 1:nrow(speciesData2)) {
      tem <- 0
      for (j in 1:ncol(speciesData2)) {
        if (speciesData2[i,j] > 0) {
          tem = tem + 1
        }
        pmd$richness[i] <- log(tem)
      }
    }
    # if chau=TRUE we use chao as abundance richness estimator
    # if the number of specie is low is possible that the chao estimator return NaN or Inf
    # chow estimator simply correct the richness on the basis of the distribution
    # of single and doubletones: Rich_corr=R_obs+{n_Single^2}\over{2 n_double}
  } else {
    for (i in 1:nrow(speciesData2)) {
      tem <- 0
      singleton <- 1
      doubleton <- 1
      for (j in 1:ncol(speciesData2)) {
        if (speciesData2[i,j] > 0) {
          tem = tem + 1
        }
        if (speciesData2[i,j] == 1) {
          singleton = singleton + 1
        } else if (speciesData2[i,j] == 2) {
          doubleton = doubleton + 1
        }
      }
      pmd$richness[i] <- log(tem + (singleton * singleton / (2 * doubleton)))
    }
  }
  ## calculation of eveness of the species
  abun.vec <- rowSums(speciesData2)
  maxrow <- apply(speciesData2,1,max)
  vec.eve <- numeric(0)
  vec.eve<--log(maxrow/abun.vec)
#   for (i in 1:length(maxrow)) {
#     evenness <- -log(maxrow[i] / abun.vec[i])
#     vec.eve[i] <- evenness
#   }
  
  ## calcultion of the pmd index of the species
  pmd.index <- numeric(0)
  for (i in 1:length(pmd$richness)) {
    dist <- pmd$richness[i] - vec.eve[i]
    pmd.sup <- log2(vec.eve[i] / dist)
    pmd.index <- c(pmd.index,pmd.sup)
  }
  pmd$pmd <- pmd.index
  pmd$eveness <- vec.eve
  ## Creation of Matrix
  eulerMat <-
    matrix(
      c(pmd$richness,pmd$shann,pmd$simpson,pmd$eveness,pmd$pmd), nrow = nrow(speciesData2), ncol = 5,dimnames =
        list(
          row.names(speciesData2), c("Richness","Shannon","Simpson","Evenness","pmdIndex")
        )
    )

  #eulerMat<-eulerMat[!rowSums(!is.finite(eulerMat)),]
  if (noNaN) {
    # removal of NaN and Inf Value 
    eulerMat <- eulerMat[!rowSums(!is.finite(eulerMat)),]
    # removal of row having richness less than adjusted threshold 
    eulerMat <- eulerMat[eulerMat[,1] > threshold,]
  }
  ## removal of row having richness less than adjusted threshold
  #eulerMat<-eulerMat[eulerMat[,1]> threshold,]
  pmd$eulerMat <- data.frame(eulerMat)
  return(pmd)
}

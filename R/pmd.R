pmd <-
function(x,threshold = 1.1,noNaN=TRUE){
  
        speciesData2<-t(x)
        ## calculation of the richness of species data
        shann<-shannon(x)
        simpson<-simpson(x)
        richness<-numeric(nrow(speciesData2))
        for(i in 1:nrow(speciesData2)){
            tem<-0
            for(j in 1:ncol(speciesData2)){
                if(speciesData2[i,j]!=0){
                    tem=tem+1
                }                
            }            
            richness[i]<-log(tem)
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
        for(i in 1:length(richness)){
            dist<-richness[i]-vec.eve[i]
            pmd<-log2(vec.eve[i]/dist)
            pmd.index<-c(pmd.index,pmd)
        }
        ## Creation of Matrix 
        eulerMat <- matrix(c(richness,shann,simpson,vec.eve,pmd.index), nrow = nrow(speciesData2), ncol = 5,dimnames=list(row.names(speciesData2), c("Richness","Shannon","Simpson","Evenness","pmdIndex")))
        ## removal of NaN and Inf Value 
        #eulerMat<-eulerMat[!rowSums(!is.finite(eulerMat)),]
        if (noNaN){
          #eulerMat<-eulerMat[!rowSums(is.nan(eulerMat)),]
          eulerMat<-eulerMat[!rowSums(!is.finite(eulerMat)),]
          eulerMat<-eulerMat[eulerMat[,1]> threshold,]
        }
	      ## removal of row having richness less than adjusted threshold 
	        #eulerMat<-eulerMat[eulerMat[,1]> threshold,]
        return(data.frame(eulerMat))
    }

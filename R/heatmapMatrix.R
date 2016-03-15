heatmapMatrix <-
function(data,n,index="PMD"){
        
        ###calculate the values of pmd index for the data to insert in matrix
        print("Started processing")
        res.fram<-as.data.frame(pmd(data))
        print("Stop pmd")
        #data<-data.matrix(data)
        pmdVec<-res.fram$pmdIndex
        richness<-res.fram$Richness
        ## select the n number of max insertion site counts and making other as zero
        
        data<-t(data)
        print("Stop data transformation")
        abun.vec<-rowSums(data)
        print("Stop summation")
        cycle<-0
        cycDisp<-0
        for(i in 1:nrow(data)){
            if(cycle>=2) {cycDisp<-cycle+cycDisp;print(cycDisp);cycle=0}   
            cycle<-cycle+1 
            vec<-tail(sort(data[i,]),n)
            vec<-vec[ vec != 0 ]
            minVec<-min(vec)
            for(j in 1:ncol(data)){
                
                if(data[i,j]<minVec)
                {
                    data[i,j]<- -100
                }
            }
        }
        print("Stop cycle1")
        ### replacing the max number values with pmd values of their respective experiment 
        cycle<-0
        cycDisp<-0
        for(i in 1:nrow(data)){
            if(cycle>=2) {cycDisp<-cycle+cycDisp;print(cycDisp);cycle=0}   
            cycle<-cycle+1
            maxVec<-max(data[i,])
            for(j in 1:ncol(data)){
                
                if(data[i,j] > 0)
                    
                {
                    
                    if (index=="PMD"){
                        
                        data[i,j]<-data[i,j]/maxVec*pmdVec[i]
                    }else if (index=="Evenness"){
                        data[i,j]<-log(abun.vec[i]/data[i,j])
                    }else{
                        data[i,j]<-(data[i,j]/maxVec*richness[i])
                        
                    }
                    
                }
            }
        }
        ## returing the matrix
        print("Stop cycle 2")
        data1<-t(data)
        data1[data1 == -100] <- NA
        findata1<-data1[rowSums(is.na(data1)) < ncol(data1),]
        findata1[is.na(findata1)] <- -100
        print("Stop processing")
        return(findata1)
    }
pmdSamp <-
function(x,ntime=30,frac=0.5){
    
    ##exporting and transposing data 
    speciesData2<-t(x)
    ## Final vectors for storing values 
    mean.rich<-numeric()
    sd.rich<-numeric()
    mean.shann<-numeric()
    sd.shann<-numeric()
    mean.simp<-numeric()
    sd.simp<-numeric()
    mean.even<-numeric()
    sd.even<-numeric()
    mean.pmd<-numeric()
    sd.pmd<-numeric()
    ## looping over the matrix
    for(i in 1:nrow(speciesData2)){
        ## taking from row only non zero values and length of vector is used for selecting sampling number from same vector
        
        temp<-subset(speciesData2[i,],speciesData2[i,] !=0)
        temp2<-reconstruct(temp) 
        len<-length(temp2)
        ### temporary vectors for saving each of the row variables 
        
        even.vec<-numeric()
        shan.vec<-numeric()
        simp.vec<-numeric()
        richn.vec<-numeric()
        pmd.vec<-numeric()
        ### sampling decided by the user 
        y<-0
        while(y<ntime){
            y<-y+1
            samp<-sample(temp2,size=frac*len)
            
            uniqu_species<-unique(samp)
            
            vec.abundance <- numeric(0)
            
            #for(i in 1:length(uniqu_species)){
            #    vec.abundance<-c(vec.abundance,length(samp[samp==uniqu_species[i]]))
            #}
            vec.abundance<-table(samp)
            
            richness<-log(length(uniqu_species))
            evenness<--log(max(vec.abundance)/round(frac*len))
            shan.vec<-c(shan.vec,shannon(as.matrix(vec.abundance)))
            simp.vec<-c(simp.vec,simpson(as.matrix(vec.abundance)))
            #calculation of pmd
            dist<-richness-evenness
            pmd<-log2(evenness/dist)
            #saving in vectors 
            even.vec<-c(even.vec,evenness)
            richn.vec<-c(richn.vec,richness)
            pmd.vec<-c(pmd.vec,pmd)
            
        }
        ### filling of the final vectors 
        mean.rich<-c(mean.rich,mean(richn.vec))
        sd.rich<-c(sd.rich,sd(richn.vec))
        mean.shann<-c(mean.shann,mean(shan.vec))
        sd.shann<-c(sd.shann,sd(shan.vec))
        mean.simp<-c(mean.simp,mean(simp.vec))
        sd.simp<-c(sd.simp,sd(simp.vec))
        mean.even<-c(mean.even,mean(even.vec))
        sd.even<-c(sd.even,sd(even.vec))
        mean.pmd<-c(mean.pmd,mean(pmd.vec))
        sd.pmd<-c(sd.pmd,sd(pmd.vec))
    }
    ## filling of matrix for the display 
    mat<-matrix(c(mean.rich,sd.rich,mean.shann,sd.shann,mean.simp,sd.simp,mean.even,sd.even,mean.pmd,sd.pmd),ncol=10,nrow=nrow(speciesData2),dimnames=list(row.names(speciesData2), c("Richness","sd.Richness","Shannon","sd.Shannon","Simpson","sd.Simpson","Evenness","sd.Evenness","pmdIndex","sd.pmdIndex")))
    return(data.frame(mat))
}

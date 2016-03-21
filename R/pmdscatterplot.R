pmdscatterplot <-
function(pmdres,nameofplot="",labelling=T,labelsize=3,spotRadius=7){

    x<-as.data.frame(pmdres)
    pmd_Index<-x[,5]


  d <- qplot(x[,1], x[,4], data=x, 
               colour= pmd_Index,
               xlab="Richness",
               ylab="Evenness",
               main=nameofplot)
              
  if(labelling==TRUE){
         #d + scale_colour_gradient2(low="red", high="red") +
    d+geom_point(size=spotRadius)+
    geom_abline(slope=0.5,colour="red",linetype=3) + 
    geom_abline(slope=1,colour="black",linetype=1) + 
    geom_text(data = x, aes(Richness,Evenness, label = row.names(x)),
              colour="black", 
              vjust=0,
              hjust = 1,
              size=labelsize) + 
    geom_abline(slope=0,colour="black",linetype=1) + 
    ylim(0,max(x[,1])) + 
    xlim(0,max(x[,1])) +theme(legend.position="none")
  }
  else{

          #d + scale_colour_gradient2(low="red", high="blue") + 
    d+geom_point(size=spotRadius)+
    geom_abline(slope=0.5,colour="red",linetype=3) + 
    geom_abline(slope=1,colour="black",linetype=1) + 
    geom_abline(slope=0,colour="black",linetype=1) + 
    ylim(0,max(x[,1])) + 
    xlim(0,max(x[,1]))        

  }      
}
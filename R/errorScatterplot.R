
errorScatterplot <-
function(pmdRes,labelling=TRUE){
##matrix to data frame
    res.fram<-as.data.frame(pmdRes)
## measure to define dimention of plot
    maxrich<-max(ceiling(res.fram$Richness))
## Vectors to plot
    Richness<-res.fram$Richness
    Evenness<-res.fram$Evenness
## first layer of ggplot2
    p<-ggplot(data = res.fram,aes(x =Richness ,y = Evenness))+ xlim(0,maxrich) + ylim(0,maxrich) + geom_errorbar(aes(ymin = Evenness-res.fram$sd.Evenness,ymax =Evenness+res.fram$sd.Evenness)) + geom_errorbarh(aes(xmin = Richness-res.fram$sd.richness,xmax =Richness+res.fram$sd.richness )) + geom_abline(slope=c(0.5,1,0),colour=c("red","black","black"),linetype=c(3,1,1))
## Condition if true will add the labels of the sample name or otherwise
    if(labelling==TRUE){
        p + geom_point() + geom_text(data = res.fram, aes(Richness,Evenness, label = row.names(res.fram)), hjust = 1,size=3)
    }
    else{
        p + geom_point()
    }

}

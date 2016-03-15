pmdbarplot <-function(pmdres,index='pmd',labAngle=45,heig=-1,fontSize=1,color=NULL) {

  #     x<-as.data.frame(pmdres)
#     x$sample<-row.names(x)
#     x$colour <- ifelse(x$pmdIndex < 0, "green","red")
#     x$hjust <- ifelse(x$pmdIndex > 0, 1.3, -0.3)
#     
#     
#     
#     ggplot(x, aes_string(sample, col, label = sample,size=1)) + 
#           geom_text(aes(y = 0,colour = colour)) + 
#           geom_bar(stat = "identity",aes(fill = colour)) 
#     
#     last_plot() + coord_flip() + scale_x_discrete(breaks = "cut") + theme_bw() +
#         theme(legend.position = "none")

  if (index=="shannon")
    #barplot(pmdres[,2],names.arg=rownames(pmdres))
    x<-barplot(pmdres[,2],xaxt="n",col=color)
  else if (index=="simpson")
    x<-barplot(pmdres[,3],xaxt="n",col=color)
  else if (index=="rich")
    x<-barplot(pmdres[,1],xaxt="n",col=color)
  else if (index=="even")
    x<-barplot(pmdres[,4],xaxt="n")
  else x<-barplot(pmdres[,5],xaxt="n",col=color)
  
  labs <- rownames(pmdres)
  text(x=x-.25, y=heig, labs, xpd=TRUE, srt=labAngle,cex=fontSize)
}

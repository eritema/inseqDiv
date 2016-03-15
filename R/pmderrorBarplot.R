pmderrorBarplot<-function(pmdres,index="pmd") 
{
  
    x <- as.data.frame(pmdres)

    x$sample <- row.names(x)

if (index=="pmd") {
    x$colour <- ifelse(x$pmdIndex < 0, "green", "red")
  
    x$hjust <- ifelse(x$pmdIndex > 0, 1.3, -0.3)
    ggplot(x, aes(sample, pmdIndex, label = sample, hjust = hjust
              , size = 1)) + geom_text(aes(y = 0, colour = colour)) + 
    geom_bar(stat = "identity", aes(fill = colour)) + 
    geom_errorbar(aes(ymax=pmdIndex+sd.pmdIndex, ymin=pmdIndex-sd.pmdIndex), position="identity",size=0.4)

    last_plot() + coord_flip() + scale_x_discrete(breaks = "cut") + 
    theme_bw() + theme(legend.position = "none")
  }
else if (index=="shannon") {
    x$colour <- ifelse(x$Shannon < 0, "green", "red")
  
    x$hjust <- ifelse(x$Shannon > 0, 1.3, -0.3)
    ggplot(x, aes(sample, Shannon, label = sample, hjust = hjust
                , size = 1)) + geom_text(aes(y = 0, colour = colour)) + 
      geom_bar(stat = "identity", aes(fill = colour)) + 
      geom_errorbar(aes(ymax=Shannon+sd.Shannon, ymin=Shannon-sd.Shannon), position="identity",size=0.4)
  
    last_plot() + coord_flip() + scale_x_discrete(breaks = "cut") + 
      theme_bw() + theme(legend.position = "none")
}
    else if (index=="simpson") {
      x$colour <- ifelse(x$Simpson < 0, "green", "red")
      
      x$hjust <- ifelse(x$Simpson > 0, 1.3, -0.3)
      ggplot(x, aes(sample, Simpson, label = sample, hjust = hjust
                    , size = 1)) + geom_text(aes(y = 0, colour = colour)) + 
        geom_bar(stat = "identity", aes(fill = colour)) + 
        geom_errorbar(aes(ymax=Simpson+sd.Simpson, ymin=Simpson-sd.Simpson), position="identity",size=0.4)
      
      last_plot() + coord_flip() + scale_x_discrete(breaks = "cut") + 
        theme_bw() + theme(legend.position = "none")
    }
}

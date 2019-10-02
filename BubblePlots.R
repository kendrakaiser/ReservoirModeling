load("ResModLHS5000run.RData")
load("maxStorage.RData")

exceeds <- function(cleanedData){
  modEval<-lapply(1:21, matrix, data= NA, nrow=n, ncol=4)
  
  for (i in 1:21){
    for (j in 1:n){
      Q<-cleanedData$Q[[i]][,j]
      S<-cleanedData$stor[[i]][,j]
      maxS<-maxStorage[[i]]
      colnames(modEval[[i]])<-c('DaysStor', 'VolStor', 'DaysQlim', 'VolQlim')
      
      daysSover<- which(S > maxS)
      Sover=length(daysSover)
      modEval[[i]][j,'DaysStor']<-Sover
      if(Sover>0){
        modEval[[i]][j,'VolStor']<-sum(S[daysSover]-maxS[daysSover])
      } else{modEval[[i]][j,'VolStor']<-0}
      
      daysQover <- which(Q >qlim[1:196,2])
      Qover<-length(daysQover)
      modEval[[i]][j,'DaysQlim']<-Qover
      if (Qover>0){
        modEval[[i]][j,'VolQlim']<-sum(Q[daysQover]-qlim[daysQover,2])
      } else{modEval[[i]][j,'VolQlim']<-0}
    }
  }
  return(modEval)
}

modEval<-exceeds(both)

plotBubbles<- function(wy, Vollim){
  wydata <- as.data.frame(cbind(modEval[[wy]], bothLHS$data))
  
  s<- ggplot(wydata, aes(x=s, y=m, size = DaysStor, fill = VolStor)) +
    geom_point(shape=21) +
    scale_fill_continuous(low = "plum1", high = "purple4", limits=c(1, Vollim))+ #, breaks = c(50000, 2500000, 5000000, 10000000, 20000000)) + # range = c(0, 200000))+
    labs(size = "Days over Storage", fill = "Vol over Storage")+
    scale_size_continuous(range = c(0, 5), limits=c(1,100), breaks= c(0,25,50,75,100))+
    theme_classic() +
    theme(axis.line = element_line(color = "grey85"), axis.ticks = element_line(color = "grey85"))
  
  v<- ggplot(wydata, aes(x=s, y=m, size = DaysQlim, fill = VolQlim)) +
    geom_point(shape=21) +
    labs(size = "Days over Q Limits", fill = "Volume Over Q Limits")+
    scale_fill_continuous(low = "lightpink1", high = "red3", limits=c(1, 100000))+ #, breaks = c(50000, 2500000, 5000000)) + #, range = c(0, 5000000)
    scale_size_continuous(range = c(0, 5), limits=c(1,100), breaks= c(0,25,50,75,100))+
    theme_classic() +
    theme(axis.line = element_line(color = "grey85"), axis.ticks = element_line(color = "grey85"))
  
  
  gs<-ggplotGrob(s)
  gv<-ggplotGrob(v)
  
  nam<-paste("g",wy, sep="") 
  assign(nam, rbind(gs, gv, size = "first"), envir = .GlobalEnv)
  gg<-get(nam)
  
  gg$widths <- unit.pmax(gs$widths, gv$widths)
  grid.newpage()
  grid.draw(gg)
}

low=c(4,16)
medlow=c(10,8,18,7)
av=c(5,13,3,17,6,21,12,11)
high=c(19,1,15,2,14,9,20)

for (wy in medlow){
  plotBubbles(wy, 80000)
}
for (wy in av){
  plotBubbles(wy, 2000000)
}
for (wy in 1:21){
  plotBubbles(wy, 10000000)
}

grid.arrange(g10, g8, g18, g7, nrow=1)
grid.arrange(g5, g13, g3, g17, nrow=1)
grid.arrange(g6, g21, g12, g11, nrow=1)
grid.arrange(g19, g1, g15, g2, nrow=1)
grid.arrange(g14, g9, g20, nrow=1)
#final set of figs to show in poster
grid.arrange(g8, g18, g5, g13, g3, g21, g12, g15, g2, nrow=3)


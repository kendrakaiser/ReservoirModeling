library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(grid)


fill<- rep(c("Modeled", "Observed", "Max"), each=196)
datelab<-seq(as.Date("1997-01-01"), as.Date("1997-07-15"), by="1 day")
x<- datelab
x2<- rep(datelab, times=4)
qs<-rep(c("Inflow", "Limit","Modeled Q", "Observed Q"), each=196)
gnames<- c("g98", "g99", "g00", 'g01', "g02", "g03", "g04","g05", 'g06', "g07", "g08", "g09", "g10", "g11", "g12", "g13", "g14", "g15", "g16", "g17")

plotfn<- function(wy){
  vol<-(c(FC$storN[FC$WY == yrs[wy]],FC$AF[FC$WY == yrs[wy]],FC$maxS[FC$WY == yrs[wy]]))/10000
  df <- data.frame(fill,x,vol)
  
  q<- c(FC$Q[FC$WY == yrs[wy]], qlim[1:196,2],FC$qoN[FC$WY == yrs[wy]], FC$Qo[FC$WY == yrs[wy]])/1000
  dfQ<-data.frame(q, x2, qs)
  
  ps<-ggplot(df, aes(x=x, y=vol, fill=fill)) +
    geom_hline(yintercept=maxAF/10000,linetype="dashed")+
    geom_area(position = "identity")+
    scale_fill_manual(values = alpha(c("#7fcdbb", "#ce1256", "#15D1FF"), c(0.65,1,0.70)))+ #green red blue "#7fcdbb", "#ce1256", "#2c7fb8" #abdda4
    scale_x_date(date_breaks = "1 month", date_labels =  "%b", expand=c(0,0))+
    xlab(NULL)+
    theme_classic()+
    theme(axis.text.x = element_blank())+
    theme(legend.position="none")+
    ylab("Storage (10,000 ac-ft)")+
    coord_cartesian(ylim=c(minS/10000,maxAF/10000))
  
  pq<-ggplot(dfQ, aes(x=x2, y=q, colour=qs, group= qs, linetype = qs))+
    geom_line() +
    scale_color_manual(values=c("#1018A4","#969696","#ce1256", "#15B9FF"))+
    scale_linetype_manual(values=c("solid","dashed","solid","solid")) +
    scale_x_date(date_breaks = "1 month", date_labels =  "%b", expand=c(0,0))+
    theme_classic()+
    theme(legend.position="none") +
    scale_y_continuous(breaks=seq(0,12,2))+
    ylab("Discharge (1,000 cfs)")+
    xlab("Date")
  
  gs<-ggplotGrob(ps)
  gq<-ggplotGrob(pq)
  
  nam<-paste("g",wy, sep="") 
  assign(nam, rbind(gs, gq, size = "first"), envir = .GlobalEnv)
  gg<-get(nam)
  
  gg$widths <- unit.pmax(gs$widths, gq$widths)
  grid.newpage()
  grid.draw(gg)
}
#plotfn(20)

for (wy in 1:20){
  plotfn(wy)
}


grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, nrow=1)
grid.arrange(g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, nrow=1)


  vol<-(c(stor[1:195,wy],afmg[1:195,wy],maxS[1:195,wy]))/10000
  df <- data.frame(fill,x,vol)

  
ps<-ggplot(df, aes(x=x, y=vol, fill=fill)) +
    geom_hline(yintercept=26.4900,linetype="dashed")+
    geom_area(position = "identity")+
    scale_fill_manual(values = alpha(c("#7fcdbb", "#ce1256", "#15D1FF"), c(0.65,1,0.70)))+ #green red blue "#7fcdbb", "#ce1256", "#2c7fb8" #abdda4
    scale_x_date(date_breaks = "1 month", date_labels =  "%b", expand=c(0,0))+
    xlab(NULL)+
    theme_classic()+
    theme(axis.text.x = element_blank())+
    theme(legend.position="none")+
    ylab("Storage (10,000 ac-ft)")+
    coord_cartesian(ylim=c(0,27), expand=c(0,0))


min(stor[,20])

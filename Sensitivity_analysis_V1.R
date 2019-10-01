#s: number of days to discharge water over
#m: number of prior days to estimate change in storage

library(pse)
library(DescTools)
library(ggplot2)
library(hexbin)
library(parallel)
library(grid)
library(gridExtra)

#----------------------------------------------------------------------------
#set parameters
maxDays=28
q.arg<- list(list("min"=1, "max"=maxDays), list("min"=1, "max"=maxDays))
names(q.arg)<-c("s", "m")
factors<-c("s", "m")
n=5000


#define probability function - discrete uniform density function for integer parameters
qdunif<-function(p, min, max){
  floor(qunif(p, min, max))}
#wrap model
modelRun<-function(params){
  return(mapply(outflowStor, params[,1], params[,2]))
}

#----------------------------------------------------------------------------
# Re-organize output for analysis - all runs for each wy in one matrix
#----------------------------------------------------------------------------
cleanData<-function(LHSrun){
  wy_Q<-lapply(1:21, matrix, data= NA, nrow=196, ncol=n)
  wy_stor<-lapply(1:21, matrix, data= NA, nrow=196, ncol=n)
  maxStor<-matrix(data= NA, nrow=196, ncol=21)
  reps = seq(1, length(LHSrun), 21) #this is whats breaking it 
  
  count=0
  for (i in 1:21){
    for (j in 1:n){
      wy_Q[[i]][,j]<- LHSrun[[reps[j]+count]][,5]
      wy_stor[[i]][,j]<- LHSrun[[reps[j]+count]][,3]
    }
    count=count+1
    maxStor[,i]<- LHSrun[[i]][,1]
  }
  out<-list(wy_Q, wy_stor, maxStor)
  names(out)<-(c("Q", "stor", "maxStor"))
  return(out)
}

#----------------------------------------------------------------------------
#calculates the number of days model exceeded max storage or discharge limits 
# and total wy volume of exceedance
#----------------------------------------------------------------------------

exceeds <- function(cleanedData){
  modEval<-lapply(1:21, matrix, data= NA, nrow=n, ncol=4)
  
  for (i in 1:21){
    for (j in 1:n){
      Q<-cleanedData$Q[[i]][,j]
      S<-cleanedData$stor[[i]][,j]
      maxS<-cleanedData$maxStor[,i]
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


#create hypercubes with different parameter sets
bothLHS <-LHS(model = NULL, factors, N=n, q='qdunif', q.arg, nboot=1)

#new.cluster <- parallel::makePSOCKcluster(c("localhost", "localhost"))
#clusterLHS <- LHS(modelRun(bothLHS$data), factors, N=50, cl = new.cluster)
#stopCluster(new.cluster)

outB<-modelRun(bothLHS$data)
Sys.time()

both<-cleanData(outB)
modEval<-exceeds(both)

maxLimDays <- matrix(data=NA, nrow=21, ncol=2)
for (wy in 1:21){
  wydat<-modEval[[wy]]
  maxLimDays[wy, 1] <-max(wydat[,1]) #DaysStor
  maxLimDays[wy, 2] <-max(wydat[,3]) #DaysQlim
}

#----------------------------------------------------------------------------
#plot
## add ability to plot this for all years on one figure with equal scales -- then subset by below, average, or above average
wy=2

plotBubbles<- function(wy){
  wydata <- as.data.frame(cbind(modEval[[wy]], bothLHS$data))
  
  d <- ggplot(data=wydata, aes(x=s, y=m))
  d + geom_hex(bins=15)
  
  s<- ggplot(wydata, aes(x=s, y=m, size = DaysStor, fill = VolStor)) +
    geom_point(shape=21) +
    scale_fill_continuous(low = "plum1", high = "purple4", limits=c(1,20000000))+ #, breaks = c(50000, 2500000, 5000000, 10000000, 20000000)) + # range = c(0, 200000))+
    labs(size = "Days over Storage", fill = "Vol over Storage")+
    scale_size_continuous(range = c(0, 5), limits=c(1,100), breaks= c(0,25,50,75,100))+
    theme_classic() +
    theme(axis.line = element_line(color = "grey85"), axis.ticks = element_line(color = "grey85"))
  
  v<- ggplot(wydata, aes(x=s, y=m, size = DaysQlim, fill = VolQlim)) +
    geom_point(shape=21) +
    labs(size = "Days over Q Limits", fill = "Volume Over Q Limits")+
    scale_fill_continuous(low = "lightpink1", high = "red3", limits=c(1,5000000))+ #, breaks = c(50000, 2500000, 5000000)) + #, range = c(0, 5000000)
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

for (wy in 1:20){
  plotBubbles(wy)
  }

grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, nrow=3)
grid.arrange(g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, nrow=3)

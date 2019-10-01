#s: number of days to discharge water over
#m: number of prior days to estimate change in storage

library(pse)
library(DescTools)
library(ggplot2)
library(hexbin)

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
  reps = seq(1, length(LHSrun), 21)
  
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

#----------------------------------------------------------------------------
#set parameters
maxDays=28
q.arg<- list(list("min"=1, "max"=maxDays), list("min"=1, "max"=maxDays))
names(q.arg)<-c("s", "m")
factors<-c("s", "m")
n=500

#create hypercubes with different parameter sets
bothLHS <-LHS(model = NULL, factors, N=n, q='qdunif', q.arg, nboot=1)
outB<-modelRun(bothLHS$data)

both<-cleanData(outB)
modEval<-exceeds(both)

#----------------------------------------------------------------------------
#plot
## add ability to plot this for all years on one figure with equal scales -- then subset by below, average, or above average
wy=1
for (wy in 1:21){
  wydata <- as.data.frame(cbind(modEval[[wy]], bothLHS$data))
  d <- ggplot(data=wydata, aes(x=s, y=m))
  d + geom_hex()
  
  pl<- ggplot(wydata, aes(x=s, y=m, size = DaysStor, fill = VolStor)) +
    geom_point(shape=21)+
    scale_fill_continuous(low = "plum1", high = "purple4")+
    labs(size = "Days over Storage Limits", fill = "Volume Over Storage Limits")
  pl
  
  ql<- ggplot(wydata, aes(x=s, y=m, size = DaysQlim, fill = VolQlim)) +
    geom_point(shape=21)+
    labs(size = "Days over Q Limits", fill = "Volume Over Q Limits")+
    scale_fill_continuous(low = "lightpink1", high = "red3")
  ql
}


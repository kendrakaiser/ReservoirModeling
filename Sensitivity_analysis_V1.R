#s: number of days to discharge water over
#m: number of prior days to estimate change in storage

library(pse)
library(DescTools)

#define probability function - discrete uniform density function for integer parameters
qdunif<-function(p, min, max){
  floor(qunif(p, min, max))}
#wrap model
modelRun<-function(params){
  return(mapply(outflowStor, params[,1], params[,2]))
}

#set parameters
maxDays=28
q.arg<- list(list("min"=1, "max"=maxDays), list("min"=1, "max"=maxDays))
names(q.arg)<-c("s", "m")
factors<-c("s", "m")


#myLHS<-LHS(model=modelRun, factors, N=100, q='qdunif', q.arg, nboot=4)
n=5
#create hypercubes with different parameter sets
bothLHS <-LHS(model = NULL, factors, N=n, q='qdunif', q.arg, nboot=1)
outB<-modelRun(bothLHS$data)

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

both<-cleanData(outB)

#turn this into a function that reads in LHS output
#Days managers or model exceeded max storage or discharge limits

exceeds <- function(cleanedData){
  modEval<-lapply(1:21, matrix, data= NA, nrow=4, ncol=n)
  
  for (i in 1:21){
    for (j in 1:n){
      Q<-cleanedData$Q[[i]][,j]
      S<-cleanedData$stor[[i]][,j]
      maxS<-cleanedData$maxStor[,i]
      rownames(modEval[[i]])<-c('DaysStor', 'VolStor', 'DaysQlim', 'VolQlim')
      
      daysSover<- which(S > maxS)
      Sover=length(stor_exceed)
      modEval[[i]]['DaysStor',j]<-Sover
      if(Sover>0){
        modEval[[i]]['VolStor',j]<-sum(S[daysSover]-maxS[daysSover])
      }

      daysQover <- which(Q >qlim[1:196,2])
      Qover<-length(daysQover)
      modEval[[i]]['DaysQlim',j]<-Qover
      if (Qover>0){
        modEval[[i]]['VolQlim',j]<-sum(Q[daysQover]-qlim[daysQover,2])
      } else{modEval[[i]]['VolQlim',j]<-0}
    }
  }
  return(modEval)
}

modEval<-exceeds(both)
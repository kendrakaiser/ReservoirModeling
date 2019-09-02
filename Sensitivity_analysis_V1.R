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
q.arg<- list(list("min"=1, "max"=maxDays), list("min"=1, "max"=maxDays))
names(q.arg)<-c("s", "m")
factors<-c("s", "m")
maxDays=28


#myLHS<-LHS(model=modelRun, factors, N=100, q='qdunif', q.arg, nboot=4)
n=500
#create hypercubes with different parameter sets
bothLHS <-LHS(model = NULL, factors, N=n, q='qdunif', q.arg, nboot=1)
outB<-modelRun(bothLHS$data)




#turn this into a functino that reads in LHS output
#Days managers or model exceeded max storage or discharge limits
Mod_exceedDate=matrix(data=NA, nrow = 60, ncol = 21)
Mod_overQlim=matrix(data=NA, nrow = 110, ncol = 21)
DaysOver_Mod=matrix(data=NA, nrow = 21, ncol = 3)
colnames(DaysOver_Mod)<-c("ModOverS", "ModOverQ", "ModVolQ")

for (wy in 1:21){
  MaxStor<- results[[wy]][,1]
  ModStor<- results[[wy]][,3]
  ModQ<- results[[wy]][,5]
  
  Mod_exceed <- which(ModStor > MaxStor)
  ll=length(Mod_exceed)
  if(ll>0){
    Mod_exceedDate[1:ll, wy]<-Mod_exceed
  }
  topped <- which(ModStor > maxAF)
  lll<-length(topped)
  if (lll>0){
    days_topped[1:lll, wy]<-topped
  }
 
  Mod_Q <- which(ModQ >qlim[1:196,2])
  mll<-length(Mod_Q)
  if (mll>0){
    Mod_overQlim[1:mll, wy]<-Mod_Q
    Mod_Vol_over_Qlim<-sum(ModQ[Mod_Q]-qlim[Mod_Q,2])
  } else{Mod_Vol_over_Qlim<-0}
  
  DaysOver_Mod[wy,]<-c(ll,mll, Mod_Vol_over_Qlim)
}
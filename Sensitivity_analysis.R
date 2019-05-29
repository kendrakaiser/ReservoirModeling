#s: number of prior days to estimate change in storage
#m: planning window

### Need to figure out how to run the model and save all output so we can do sensitivity analysis on storage:maxStorage and discharge

library(pse)

#define probability function - discrete uniform density function for integer parameters
qdunif<-function(p, min, max){
  floor(qunif(p, min, max))}
#wrap model
modelRun<-function(params){
  return(mapply(outflowStor, params[,1], params[,2]))
}

#set parameters
q.arg<- list(list("min"=1, "max"=10), list("min"=1, "max"=15))
names(q.arg)<-c("s", "m")
factors<-c("s", "m")


#myLHS<-LHS(model=modelRun, factors, N=100, q='qdunif', q.arg, nboot=4)
n=100
#create hypercube 
bothLHS <-LHS(model = NULL, factors, N=n, q='qdunif', q.arg, nboot=1)
#bothLHS<-tell(bothLHS, bothLHS$data) #res<-get.results(bothLHS)
outB<-modelRun(bothLHS$data)

#pse plots -- not sure how helpful
#plotscatter(bothLHS,index.res=c(5, 8, 10),  add.lm=FALSE)#stack=TRUE, index.res=c(250, 255, 260)
#plotecdf(bothLHS, stack=TRUE)
#plotprcc(bothLHS, stack=TRUE)

#set only S to change
q.argS<- list(list("min"=1, "max"=10), list("min"=6, "max"=8))
sLHS<-LHS(model = NULL, factors, N=n, q='qdunif', q.argS, nboot=1)
out_S<-modelRun(sLHS$data)
#set only M to change
q.argM<- list(list("min"=2, "max"=6), list("min"=1, "max"=15))
mLHS<-LHS(model = NULL, factors, N=n, q='qdunif', q.argM, nboot=4)
out_M<-modelRun(mLHS$data)

#--------------------------------
# Re-organize output for analysis
#--------------------------------
wy_Q<-lapply(1:21, matrix, data= NA, nrow=196, ncol=n)
count=0
for (i in 1:21){
  for (j in 1:r){
    count=count+1
    wy_Q[[i]][,j]<- outB[[i+(i*r)]][,5] #I dont think this is right 
  }
}

matplot(wy_Q[[7]], type='l')


#--------------------------------
#Mean discharge for each day under all model runs with confidence intervals
#--------------------------------

OutMeans<-matrix(data=NA, nrow = 196, ncol = 21)
z<-lapply(1:21, matrix, data=NA, nrow = 196, ncol =2)
sd_doy<-matrix(data=NA, nrow = 196, ncol = 21)
for (i in 1:21){
  OutMeans[,i]<-rowMeans(wy_Q[[i]], na.rm = FALSE, dims = 1)
  z[[i]]<-apply(wy_Q[[i]], 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
  sd_doy[,i]<- apply(wy_Q[[i]], 1, sd)
}

matplot(OutMeans, type='l')

OutMeans[,2]<-rowMeans(out_S, na.rm = FALSE, dims = 1)
OutMeans[,3]<-rowMeans(out_M, na.rm = FALSE, dims = 1)
zS<-apply(out_S, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
zM<-apply(out_M, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
sd_S<- apply(out_S, 1, sd)
sd_M<- apply(out_M, 1, sd)

plot(sd_all, type='l')
lines(1:4116, sd_S, type='l', col='green')
lines(1:4116, sd_M, type='l', col='blue')

#--------------------------------
#calculate days over the discharge limits, and total volume over discharge limits as a function of s/m

days_over<-matrix(data=NA, nrow = 21, ncol = 100)
vol_over<-matrix(data=NA, nrow = 21, ncol = 100)
for (i in 1: 100){
  for (wy in 1:21){
    wy_out= out[FC$WY == yrs[wy],i]
    ids<-which(wy_out > 10000)
    days_over[wy, i] = length(ids)
    vol_over[wy, i] =sum(wy_out[ids] -10000)
  }
}

plot(bothparams$m, vol_over[1,])
plot(bothparams$s, vol_over[1,])

###s
days_over_S<-matrix(data=NA, nrow = 21, ncol = 50)
vol_over_S<-matrix(data=NA, nrow = 21, ncol = 50)
for (i in 1: 50){
  for (wy in 1:21){
    wy_out= out_S[FC$WY == yrs[wy],i]
    ids<-which(wy_out > 10000)
    days_over_S[wy, i] = length(ids)
    vol_over_S[wy, i] =sum(wy_out[ids] -10000)
  }
}

plot(Sparams$m, vol_over_S[1,])
plot(Sparams$s, vol_over_S[1,])

###M
days_over_M<-matrix(data=NA, nrow = 21, ncol = 100)
vol_over_M<-matrix(data=NA, nrow = 21, ncol = 100)
for (i in 1: 100){
  for (wy in 1:21){
    wy_out= out_M[FC$WY == yrs[wy],i]
    ids<-which(wy_out > 10000)
    days_over_M[wy, i] = length(ids)
    vol_over_M[wy, i] =sum(wy_out[ids] -10000)
  }
}

plot(Mparams$m, vol_over_M[1,])
plot(Mparams$s, vol_over_M[1,])

par(mfrow=c(3,1)) 

plot(OutMeans[,1], type='l', lwd='2')
lines(1:4116, z[2,], type='l', col='blue')
lines(1:4116, z[1,], type='l', col='blue' )

plot(OutMeans[,2], type='l', lwd='2')
lines(1:4116, zS[2,], type='l', col='blue')
lines(1:4116, zS[1,], type='l', col='blue')

plot(OutMeans[,3], type='l', lwd='2')
lines(1:4116, zM[2,], type='l', col='blue')
lines(1:4116, zM[1,], type='l', col='blue')

par(mfrow=c(1,1)) 
plot(OutMeans[,1], type='l', lwd='2')
lines(1:4116, OutMeans[,2], type='l', lwd='1', col='blue')
lines(1:4116, OutMeans[,3], type='l', lwd='1', col='green')
lines(1:4116, rep(10000, 4116), lty=3, col="grey")

OutSum<-matrix(data=NA, nrow = jul, ncol = 100)
obsSum<-matrix(data=NA, nrow = jul, ncol = 21)

for (wy in 1:21){
  obsSum[,wy]<-cumsum(FC$Qo[FC$WY == wy])
  
  q=out[FC$WY == wy]
  for (i in 1:100){
    OutSum[,i]<-cumsum(q[,i])
  }
  
}







#cumulative sum of discahrge

matplot(OutSum, type='l', col='blue')
lines(FC$doy[1:jul], obsSum, type='l', col='black', lwd='2')

matplot(out, type='l', col='blue')
lines(OutMean, type='l', col="black", lwd=2)

myLHS <-LHS(model = modelRun, factors, N=8, q='qdunif', q.arg, nboot=4)



print.LHS <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("Model:\n"); print (x$model);
  cat("Factors:\n"); print (x$factors);
  cat("Results:\n"); print (x$res.names);
  cat("PRCC:\n"); print (x$prcc);
}

tstr<-print.LHS(myLHS)
tstr[[1]][["y"]]

data<-myLHS$data
index.res<-1:get.noutputs(myLHS)
index.data <- 1:get.ninputs(myLHS)
res<-myLHS$res

dat <- as.vector(get.results(myLHS)[,index.res[1]])
g <- rep(index.res, each=dim(obj$res)[1])
Ecdf(dat, group=g, col=col, xlab=xlab, ...)

plotprcc(myLHS, stack=TRUE)
plotscatter(myLHS, index.res=c(250, 255, 260), add.lm=FALSE)
plotprcc(myLHS, index.res=c(15, 30, 60, 90, 120, 150))
#s: number of prior days to estimate change in storage
#m: planning window

### Need to figure out how to run the model and save all output so we can do sensitivity analysis on storage:maxStorage and discharge

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
q.arg<- list(list("min"=1, "max"=14), list("min"=1, "max"=14))
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
q.argS<- list(list("min"=1, "max"=14), list("min"=6, "max"=8))
sLHS<-LHS(model = NULL, factors, N=n, q='qdunif', q.argS, nboot=1)
out_S<-modelRun(sLHS$data)
#set only M to change
q.argM<- list(list("min"=6, "max"=8), list("min"=1, "max"=14))
mLHS<-LHS(model = NULL, factors, N=n, q='qdunif', q.argM, nboot=4)
out_M<-modelRun(mLHS$data)

#------------------------------------------------------------------------------------------------
# Re-organize output for analysis - all runs for each wy in one matrix
#------------------------------------------------------------------------------------------------
cleanData<-function(LHSrun){
  wy_Q<-lapply(1:21, matrix, data= NA, nrow=196, ncol=n)
  wy_stor<-lapply(1:21, matrix, data= NA, nrow=196, ncol=n)
  reps = seq(1, length(LHSrun), 21)

  count=0
  for (i in 1:21){
    for (j in 1:n){
      wy_Q[[i]][,j]<- LHSrun[[reps[j]+count]][,5]
      wy_stor[[i]][,j]<- LHSrun[[reps[j]+count]][,3]
    }
    count=count+1
  }
  out<-list(wy_Q, wy_stor)
  names(out)<-(c("Q", "stor"))
  return(out)
}

both<-cleanData(outB)
sOut<-cleanData(out_S)
mOut<-cleanData(out_M)

#plot all outputs from the model runs

pal <- colorRampPalette(c("yellow", "green", "blue"))
cols<-pal(21)

StorSensPlot<-function(LHS, vals, wy, sm, dataOut){
  paramVal<-data.matrix(LHS$data[sm])
  colRamp<-pal(vals)[as.numeric(cut(paramVal,breaks = vals))]
  plt<-matplot(dataOut$stor[[wy]], type='l', col = colRamp)+
    ColorLegend(x='topleft', cols = pal(vals), labels=1:vals, cntrlbl=TRUE)+
    lines(dataOut$maxS[[wy]][,1], type='l') #max storage
  return(plt) ### color coded to m/s values
}

QSensPlot<-function(LHS, vals, wy, sm, dataOut){
  paramVal<-data.matrix(LHS$data[sm])
  colRamp<-pal(vals)[as.numeric(cut(paramVal,breaks = vals))]
  plt<- matplot(dataOut$Q[[wy]], type='l', col = colRamp) +
        lines(qlim[,2], type='l', lty=3, col='grey17')
  return(plt) ### color coded to m/s values
}

par(mfrow=c(2,2))
wy=10
StorSensPlot(bothLHS, 14, wy, 'm', both)
StorSensPlot(bothLHS, 14, wy, 's', both)
QSensPlot(bothLHS, 14, wy, 'm', both)
QSensPlot(bothLHS, 14, wy, 's', both)
mtext("Model with both variables changing", outer = TRUE, cex = 1.5)

par(mfrow=c(2,2))

StorSensPlot(mLHS, 14, wy, 'm', mOut)
StorSensPlot(mLHS, 2, wy, 's', sOut)
QSensPlot(mLHS, 14, wy, 'm', mOut)
QSensPlot(mLHS, 2, wy, 's', sOut)+


par(mfrow=c(2,2))
StorSensPlot(sLHS, 2, wy, 'm', mOut)
StorSensPlot(sLHS, 14, wy, 's', sOut)
QSensPlot(sLHS, 2, wy, 'm', mOut)
QSensPlot(sLHS, 14, wy, 's', sOut)
mtext("Model Variability of S", outer = TRUE, cex = 1.5)

par(mfrow=c(4,1))
StorSensPlot(bothLHS, 14, wy, 'm', both)
StorSensPlot(bothLHS, 14, wy, 's', both)
StorSensPlot(mLHS, 14, wy, 'm', mOut)
StorSensPlot(sLHS, 14, wy, 's', sOut)

for (wy in 1:21){
  par(mfrow=c(4,1))
  StorSensPlot(bothLHS, 14, wy, 'm', both)
  StorSensPlot(bothLHS, 14, wy, 's', both)
  StorSensPlot(mLHS, 14, wy, 'm', mOut)
  StorSensPlot(sLHS, 14, wy, 's', sOut)
}
#------------------------------------------------------------------------------------------------
# calculate days over the discharge limits, and total volume over discharge limits as a function of s/m
#------------------------------------------------------------------------------------------------
yrs<- 1998:2018
inflSum<-matrix(data=NA, nrow = jul, ncol = 21)
for (wy in 1:21){
  inflSum[,wy]<-cumsum(FC$Q[FC$WY == yrs[wy]])
}
totQ<-inflSum[jul,]

colRamp<-pal(21)[as.numeric(cut(totQ,breaks = 21))]


overage<-function(LHSout){
  days_over<-matrix(data=NA, nrow = n, ncol = 21)
  vol_over<-matrix(data=NA, nrow = n, ncol = 21)
  for (i in 1:n){
    for (wy in 1:21){
      wy_dat= LHSout$Q[[wy]][,i]
      ids<-which(wy_dat > 10000)
      days_over[i, wy] = length(ids)
      vol_over[i, wy] =sum(wy_dat[ids]-10000)
    }
  }
  out<-list(days_over, vol_over)
  names(out)<-c("days_over", 'vol_over')
  return(out)
}

cols<-par(21)

plt_overage<-function(LHSout, params, sm){
  Qover<-overage(LHSout)
  a<-Qover$days_over
  y<-as.data.frame(a[,1])
  x<-bothLHS$data[sm]
  plt<- plot(x[,],y[,], ylab="Days Over Discharge Limit", xlab="Discharge Planning Days", ylim=c(0.5, 35), pch=19, col=cols[1])
    for (i in 2:21){
     y<-as.data.frame(a[,i])
     x<-bothLHS$data[sm]
     points(x[,],y[,], pch=19, col=cols[i]) #color ramp isnt quite right
    } 
  return(plt)
}

plt_overage(mOut, mLHS$data, 'm')

#------------------------------------------------------------------------------------------------
# Mean discharge for each day of each water year under all model runs with confidence intervals and standard deviations
#------------------------------------------------------------------------------------------------
# need to group water years by flow - either cluster analysis to group, of high/low/average

Qstats<-function(dataOut){
  
  OutMeans<-matrix(data=NA, nrow = 196, ncol = 21)
  z<-lapply(1:21, matrix, data=NA, nrow = 196, ncol =2)
  sd_doy<-matrix(data=NA, nrow = 196, ncol = 21)
  for (i in 1:21){
    wy_Q<-dataOut$Q[[i]]
    OutMeans[,i]<-rowMeans(wy_Q, na.rm = FALSE, dims = 1)
    z[[i]]<-apply(wy_Q, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
    sd_doy[,i]<- apply(wy_Q, 1, sd)
  }
  out<-list(OutMeans, z, sd_doy)
  names(out)<-c('mean', 'ci', 'sd')
  return(out)
}
storstats<-function(dataOut){
  
  OutMeans<-matrix(data=NA, nrow = 196, ncol = 21)
  z<-lapply(1:21, matrix, data=NA, nrow = 196, ncol =2)
  sd_doy<-matrix(data=NA, nrow = 196, ncol = 21)
  for (i in 1:21){
    wy_Q<-dataOut$stor[[i]]
    OutMeans[,i]<-rowMeans(wy_Q, na.rm = FALSE, dims = 1)
    z[[i]]<-apply(wy_Q, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
    sd_doy[,i]<- apply(wy_Q, 1, sd)
  }
  out<-list(OutMeans, z, sd_doy)
  names(out)<-c('mean', 'ci', 'sd')
  return(out)
}

bothQstats<-Qstats(both)
sQstats<-Qstats(sOut)
mQstats<-Qstats(mOut)

bothSstats<-storstats(both)
sSstats<-storstats(sOut)
mSstats<-storstats(mOut)


matplot(bothSstats$sd, type='l', ylim=c(0,84000))
matplot(mSstats$sd, type='l', ylim=c(0,84000))
matplot(sSstats$sd, type='l', ylim=c(0,84000))

matplot(bothQstats$sd, type='l')
matplot(mQstats$sd, type='l')
matplot(sQstats$sd, type='l')


#------------------------------------------------------------------------------------------------
# cumulative sum of discahrge - need to update
#------------------------------------------------------------------------------------------------

OutSum<-matrix(data=NA, nrow = jul, ncol = 100)
obsSum<-matrix(data=NA, nrow = jul, ncol = 21)

for (wy in 1:21){
  q=out[FC$WY == wy]
  for (i in 1:100){
    OutSum[,i]<-cumsum(q[,i])
  }
  
}









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
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

#set only S to change - M still veries between 13/14 - can give a single value for the qdunif..
q.argS<- list(list("min"=1, "max"=maxDays), list("min"=13, "max"=15))
sLHS<-LHS(model = NULL, factors, N=n, q='qdunif', q.argS, nboot=1)
out_S<-modelRun(sLHS$data)
#set only M to change
q.argM<- list(list("min"=13, "max"=15), list("min"=1, "max"=maxDays))
mLHS<-LHS(model = NULL, factors, N=n, q='qdunif', q.argM, nboot=4)
out_M<-modelRun(mLHS$data)

#set only M to change with a lower S
q.argM<- list(list("min"=6, "max"=8), list("min"=1, "max"=maxDays))
mLHS_low<-LHS(model = NULL, factors, N=n, q='qdunif', q.argM, nboot=4)
out_M_lowS<-modelRun(mLHS_low$data)


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
sOut<-cleanData(out_S)
mOut<-cleanData(out_M)
mOut_sl<-cleanData(out_M_lowS)

#plot all outputs from the model runs

pal <- colorRampPalette(c("yellow", "green", "blue"))
cols<-pal(21)

StorSensPlot<-function(LHS, vals, wy, sm, dataOut){
  paramVal<-data.matrix(LHS$data[sm])
  colRamp<-pal(vals)[as.numeric(cut(paramVal,breaks = vals))]
  plt<-matplot(dataOut$stor[[wy]], type='l', col = colRamp, ylim=c(300000, 1100000))+
    ColorLegend(x='topleft', cols = pal(vals), labels=1:vals, cntrlbl=TRUE) +
    lines(dataOut$maxStor[,wy], type='l') #max storage
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
StorSensPlot(bothLHS, maxDays, wy, 'm', both)
StorSensPlot(bothLHS, maxDays, wy, 's', both)
QSensPlot(bothLHS, maxDays, wy, 'm', both)
QSensPlot(bothLHS, maxDays, wy, 's', both)
mtext("Model with both variables changing", outer = TRUE, cex = 1.5)

par(mfrow=c(4,1))
StorSensPlot(bothLHS, maxDays, wy, 'm', both)
StorSensPlot(bothLHS, maxDays, wy, 's', both)
StorSensPlot(mLHS, maxDays, wy, 'm', mOut)
StorSensPlot(sLHS, maxDays, wy, 's', sOut)


par(mfrow=c(2,2)) #show difference in sensitivity between m and s

StorSensPlot(mLHS, maxDays, wy, 'm', mOut)
StorSensPlot(sLHS, maxDays, wy, 's', sOut)
QSensPlot(mLHS, maxDays, wy, 'm', mOut)
QSensPlot(sLHS, maxDays, wy, 's', sOut)

wy=10
par(mfrow=c(2,2)) #show difference in sensitivity between m when s is set at 7 or 14

StorSensPlot(mLHS, maxDays, wy, 'm', mOut) #s=14
StorSensPlot(mLHS_low, maxDays, wy, 'm', mOut_sl) #s=7
QSensPlot(mLHS, maxDays, wy, 'm', mOut)
QSensPlot(mLHS_low, maxDays, wy, 'm', mOut_sl)

hml<-c(4,3,20)

par(mfcol=c(2,3))
for (wy in hml){
  StorSensPlot(bothLHS, maxDays, wy, 'm', both)
  QSensPlot(bothLHS, maxDays, wy, 'm', both)
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

overageDat<-function(LHSout, LHSmod){
  over<-overage(LHSout)
  a<<-over$days_over
  av<<-as.vector(a)

  v<<-over$vol_over
  vv<<-as.vector(v)
  
  y<<-as.matrix(LHSmod$data['s'])
  x<<-as.matrix(LHSmod$data['m'])
}

daysOverMin<-matrix(data=NA, nrow=21, ncol=3)
daysOverMin[,1]<-apply(a, 2, FUN=min)
daysOverMin[,2]<-apply(a, 2, FUN=max)
daysOverMin[,3]<-apply(a, 2, FUN=mean)

volOverMin<-matrix(data=NA, nrow=21, ncol=3)
volOverMin[,1]<-apply(v, 2, FUN=min)
volOverMin[,2]<-apply(v, 2, FUN=max)
volOverMin[,3]<-apply(v, 2, FUN=mean)

save("both", "bothLHS", "a", "v", "daysOverMin", "volOverMin", "totQ", "inflSum", file="workspace.RData")
overageDat(both, bothLHS)
par(mfrow=c(4,4))
#plot number of days over discharge limits
for (i in 1:21){
  ai<-a[,i]
  maxOver<-max(ai)
  
  if(maxOver>1){
  colRamp<<-pal(maxOver)[as.numeric(cut(ai, breaks = maxOver))]
  plot(x, y, ylab="s", xlab="m", ylim=c(1, 28), xlim=c(1, 28), pch=19, col=colRamp)+
    ColorLegend(x='top', cols = pal(maxD), labels=seq(from=0, to=82, by=10), cntrlbl=TRUE, horiz=TRUE)
  }
}
  
par(mfrow=c(4,4))
#plot volumer of water over discharge limits
for (i in 1:21){
  ai<-v[,i]
  maxOver<-max(ai)
  
  if(maxOver>1){
    colRamp<<-pal(maxOver)[as.numeric(cut(ai, breaks = maxOver))]
    plot(x, y, ylab="s", xlab="m", ylim=c(1, 28), xlim=c(1, 28), pch=19, col=colRamp)+
      ColorLegend(x='top', cols = pal(maxD), labels=seq(from=0, to=82, by=10), cntrlbl=TRUE, horiz=TRUE)
  }
}


StorOverage<-function(LHSout){
  days_over<-matrix(data=NA, nrow = n, ncol = 21)
  vol_over<-matrix(data=NA, nrow = n, ncol = 21)
  for (i in 1:n){
    for (wy in 1:21){
      wy_dat= LHSout$stor[[wy]][,i]
      ids<-which(wy_dat > LHSout$maxStor[,wy])
      days_over[i, wy] = length(ids)
      vol_over[i, wy] =sum(wy_dat[ids]-LHSout$maxStor[ids,wy])
    }
  }
  out<-list(days_over, vol_over)
  names(out)<-c("days_over_S", 'vol_over_S')
  return(out)
}

storOver<-StorOverage(both)

Sd<-storOver$days_over
dd<-as.vector(Sd)
vS<-storOver$vol_over
vv<-as.vector(vS)
maxvO<-max(vv)

par(mfrow=c(4,5))
#plot number of days or volume over storage limits
for (i in 1:21){
  ai<-Sd[,i]
  maxOver<-max(ai)
  
  if(maxOver>1){
    colRamp<<-pal(40)[as.numeric(cut(ai, breaks = 40))]
    plot(x, y, ylab="s", xlab="m", ylim=c(1, 28), xlim=c(1, 28), pch=19, col=colRamp)+
      ColorLegend(x='top', cols = pal(maxD), labels=seq(from=0, to=82, by=10), cntrlbl=TRUE, horiz=TRUE)
  }
}
#------------------------------------------------------------------------------------------------
# Mean discharge for each day of each water year under all model runs with confidence intervals and standard deviations
#------------------------------------------------------------------------------------------------
# need to group water years by flow - either cluster analysis to group, of high/low/average

Qstats<-function(dataOut){
  
  OutMeans<-matrix(data=NA, nrow = 196, ncol = 21)
  z<-lapply(1:21, matrix, data=NA, nrow = 196, ncol =2)
  sd_doy<-matrix(data=NA, nrow = 196, ncol = 21)
  cv_doy<-matrix(data=NA, nrow = 196, ncol = 21)
  for (i in 1:21){
    wy_Q<-dataOut$Q[[i]]
    OutMeans[,i]<-rowMeans(wy_Q, na.rm = FALSE, dims = 1)
    z[[i]]<-apply(wy_Q, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
    sd_doy[,i]<- apply(wy_Q, 1, sd)
    cv_doy[,i]<- apply(wy_Q, 1, cv)
  }
  out<-list(OutMeans, z, sd_doy, cv_doy)
  names(out)<-c('mean', 'ci', 'sd', 'cv')
  return(out)
}
storstats<-function(dataOut){
  
  OutMeans<-matrix(data=NA, nrow = 196, ncol = 21)
  z<-lapply(1:21, matrix, data=NA, nrow = 196, ncol =2)
  sd_doy<-matrix(data=NA, nrow = 196, ncol = 21)
  cv_doy<-matrix(data=NA, nrow = 196, ncol = 21)
  for (i in 1:21){
    wy_S<-dataOut$stor[[i]]
    OutMeans[,i]<-rowMeans(wy_S, na.rm = FALSE, dims = 1)
    z[[i]]<-apply(wy_S, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
    sd_doy[,i]<- apply(wy_S, 1, sd)
    cv_doy[,i]<- apply(wy_S, 1, cv)
  }
  out<-list(OutMeans, z, sd_doy, cv_doy)
  names(out)<-c('mean', 'ci', 'sd', 'cv')
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

high=c(1,2,9,14,15,21)
low=c(4,8)
av=c(3,5,6,7,10,11,12,13,16,17,18,19,20)

#cumulative inflow to system
Qsum<-matrix(data=NA, nrow=21, ncol=1)
for (wy in 1:21){
  csum<-cumsum(FC$Q[FC$WY == yrs[wy]])
  Qsum[wy]<-csum[196]
}

#color ramp of cumulative inflow
paramVal<-data.matrix(Qsum)
colRamp<-pal(21)[as.numeric(cut(paramVal,breaks = 21))]

par(mfcol=c(3,1))
matplot(bothQstats$cv[,high], type='l', col = colRamp)
matplot(bothQstats$cv[,av], type='l', col = colRamp)
matplot(bothQstats$cv[,low], type='l', col = colRamp)

#------------------------------------------------------------------------------------------------
# cumulative sum of discharge - need to update
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
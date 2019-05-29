#Initial analysis of reservoir model output and manager operations

m=10
s=5

results<- outflowStor(s,m)

#plot the initial results
for (wy in 1:2){
  plot(results[[wy]][,1], type='l', ylim=c(300000, 1010200)) #max storage
  lines(results[[wy]][,4], col='orange') #modeled storage
  lines(FC$AF[FC$WY == yrs[wy]], col='green') #actual storage
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(results[[wy]][,7], type='l', lty=3, col='orange', ylim=c(0,16000)) #modeled discharge
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
}


#make a few plots - e.g number, timing and volume that the managers exceeded the maxS for the day - and other plots from the nicholas data site
#determine average day they start drafting for irrigation - e.g. once dsdt is negative and never goes positive again
#--------------
#Days model exceeded max storage = 170
exceed <- which(FC$stor > FC$maxS)
exceedDate=matrix(data=NA, nrow = length(exceed), ncol = 2)
exceedDate[,1] <- FC$WY[exceed]
exceedDate[,2] <- FC$doy[exceed]
hist(exceedDate[,2])

#Days managers exceeded max storage limits = 228
Mexceed <- which(FC$AF > FC$maxS)
MexceedDate=matrix(data=NA, nrow = length(Mexceed), ncol = 2)
MexceedDate[,1] <- FC$WY[Mexceed]
MexceedDate[,2] <- FC$doy[Mexceed]
hist(MexceedDate[,2])

topped <- which(FC$stor > maxAF)
FC$topped<- FC$stor > maxAF
hist(FC$doy[FC$topped == 'TRUE'])

#direction of change in storage and last day of increasing sotrage
#coefficient of variation for each year of managed system, natural, and modeled
library(sjstats)
dervDS=matrix(data=NA, nrow = 195, ncol = 1)
TotalQ=matrix(data=NA, nrow = 21, ncol = 1)
cvQ=matrix(data=NA, nrow = 21, ncol = 1)
cvQman=matrix(data=NA, nrow = 21, ncol = 1)
cvQmod=matrix(data=NA, nrow = 21, ncol = 1)
DOY_decreasingS=matrix(data=NA, nrow = 21, ncol = 1)

for (wy in 1:21){
  stor=FC$AF[FC$WY == yrs[wy]]
  Q=FC$Q[FC$WY == yrs[wy]]
  Qo=FC$Qo[FC$WY == yrs[wy]]
  qo=FC$qo[FC$WY == yrs[wy]]
  TotalQ[wy] = sum(Q)
  cvQ[wy]= cv(Q)
  cvQman[wy]= cv(Qo)
  cvQmod[wy]= cv(qo)
  for(i in 1:195){
    dervDS[i]=(stor[i+1]/stor[i])
  }
  DOY_decreasingS[wy]=max(which(dervDS > 1))
}

#last day of increasing storage as a function of inflow
plot(TotalQ, DOY_decreasingS)


plot(cvQmod)
plot(cvQ, cvQman)





library(RColorBrewer)
cols <- brewer.pal(9, "YlGnBu")
pal <- colorRampPalette(cols[3:9])
zz=pal(28)

matplot(1:196, qoS, type='l', col=zz, xlab= 'Day of Year', ylab= 'Discharge (cfs)')
legend(1, 12000, legend=c("1", "7", "14", "28"),
       col=c(zz[1], zz[7], zz[14], zz[28]), title = "Storage Forecast Days", lty=1:2, cex=0.8,
       box.lty=0)

matplot(1:196, storS, type='l', col=pal(28), xlab= 'Day of Year', ylab= 'Storage (AF)')
legend(1, 9e+05, legend=c("1", "7", "14", "28"),
       col=c(zz[1], zz[7], zz[14], zz[28]), title = "Storage Forecast Days", lty=1:2, cex=0.8,
       box.lty=0)

matplot(1:196, storFs, type='l', col=pal(28), xlab= 'Day of Year', ylab= 'Storage Forecast (AF)')
legend(1, 9e+05, legend=c("1", "7", "14", "28"),
       col=c(zz[1], zz[7], zz[14], zz[28]), title = "Storage Forecast Days", lty=1:2, cex=0.8,
       box.lty=0)

matplot(1:196, qoM, type='l', col=pal(10))
legend(1, 12000, legend=c("1", "7", "14", "28"),
       col=c(zz[1], zz[7], zz[14], zz[28]), title = "Advanced Planning Days", lty=1:2, cex=0.8,
       box.lty=0)

matplot(1:196, storM, type='l', col=pal(10), xlab= 'Day of Year', ylab= 'Storage (AF)')
legend(1, 9e+05, legend=c("1", "7", "14", "28"),
       col=c(zz[1], zz[7], zz[14], zz[28]), title = "Advanced Planning Days", lty=1:2, cex=0.8,
       box.lty=0)

matplot(1:196, storFm, type='l', col=pal(10), xlab= 'Day of Year', ylab= 'Storage Forecast (AF)')
legend(1, 9e+05, legend=c("1", "7", "14", "28"),
       col=c(zz[1], zz[7], zz[14], zz[28]), title = "Advanced Planning Days", lty=1:2, cex=0.8,
       box.lty=0)
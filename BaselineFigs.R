# Baseline analysis of state of system as managed now

load("maxStorage.RData")
load("ResSystemData.RData")
yrs<- 1998:2018

perc<-matrix(data=NA, nrow = 196, ncol = 21)
for (wy in 1:21){
  AF=FC$AF[FC$WY == yrs[wy]]
  perc[,wy]=AF/maxStorage[[wy]]
}

plot(perc[,1], type='l', ylim=c(0.2, 1.8))
for (wy in 1:21){
lines(perc[,wy], col='grey') #modeled storage
}

percQs<-apply(perc, 1, quantile, probs = c(0.05, 0.5, 0.95),  na.rm = TRUE)

matplot(perc, col='grey', type='l') 
lines(percQs[1,], col='blue') 
lines(percQs[2,], col='black') 
lines(percQs[3,], col='blue') 

#cumulative inflow to system
Qsum<-matrix(data=NA, nrow=21, ncol=1)
for (wy in 1:21){
  csum<-cumsum(FC$Q[FC$WY == yrs[wy]])
  Qsum[wy]<-csum[196]
}
Qsum_stats<-quantile(Qsum, probs = c(0.05, .25, 0.5, .75, 0.95),  na.rm = TRUE)


low=c(4,16)
medlow=c(10,8,18,7)
av=c(5,13,3,17,6,21,12,11,19)
high=c(1,15,2,14,9,20)

#color ramp of cumulative inflow
paramVal<-data.matrix(Qsum)
colRamp<-pal(21)[as.numeric(cut(paramVal,breaks = 21))]


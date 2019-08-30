#Initial analysis of reservoir model output and manager operations

m=14
s=14

results<- outflowStor(s,m)

#plot the initial results
for (wy in 1:21){
  plot(results[[wy]][,1], type='l', ylim=c(300000, 1010200)) #max storage
  lines(results[[wy]][,3], col='orange') #modeled storage
  lines(FC$AF[FC$WY == yrs[wy]], col='green') #actual storage
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  #plot(results[[wy]][,5], type='l', lty=3, col='orange', ylim=c(0,16000)) #modeled discharge
  #lines(qlim[,2], type='l', lty=3, col='grey17')
  #lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
}

#running standard deviation of managers outflow
library(caTools)

runQsd<-runsd(x, k)


#make a few plots - e.g number, timing and volume that the managers exceeded the maxS for the day - and other plots from the nicholas data site
#determine average day they start drafting for irrigation - e.g. once dsdt is negative and never goes positive again
#--------------

#Days managers or model exceeded max storage or discharge limits
MexceedDate=matrix(data=NA, nrow = 60, ncol = 21)
Mod_exceedDate=matrix(data=NA, nrow = 60, ncol = 21)
days_topped=matrix(data=NA, nrow = 30, ncol = 21)

Mang_overQlim=matrix(data=NA, nrow = 110, ncol = 21)
Mod_overQlim=matrix(data=NA, nrow = 110, ncol = 21)

DaysOver=matrix(data=NA, nrow = 21, ncol = 6)
colnames(DaysOver)<-c("ManOverS", "ModOverS", "ManOverQ", "ModOverQ", "ManVolQ", "ModVolQ")

for (wy in 1:21){
  MaxStor<- results[[wy]][,1]
  ModStor<- results[[wy]][,3]
  ModQ<- results[[wy]][,5]
  ManQ<-FC$Qo[FC$WY == yrs[wy]]
  
  ManagedStor<- FC$AF[FC$WY == yrs[wy]]
  Mexceed <- which(ManagedStor > MaxStor)
  l=length(Mexceed)
  if(l>0){
    MexceedDate[1:l, wy] <-Mexceed
  }
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
  Mg_Q<- which(ManQ > qlim[1:196,2])
  ml<-length( Mg_Q)
  if (ml>0){
    Mang_overQlim[1:ml, wy]<-Mg_Q
    Man_Vol_over_Qlim<-sum(ManQ[Mg_Q]-qlim[Mg_Q,2])
  } else{Man_Vol_over_Qlim<-0}
  Mod_Q <- which(ModQ >qlim[1:196,2])
  mll<-length(Mod_Q)
  if (mll>0){
    Mod_overQlim[1:mll, wy]<-Mod_Q
    Mod_Vol_over_Qlim<-sum(ModQ[Mod_Q]-qlim[Mod_Q,2])
  } else{Mod_Vol_over_Qlim<-0}
  
  DaysOver[wy,]<-c(l,ll,ml,mll,Man_Vol_over_Qlim, Mod_Vol_over_Qlim)
}


#direction of change in storage and last day of increasing storage
#coefficient of variation for each year of managed system, natural, and modeled
library(sjstats)
dervDS=matrix(data=NA, nrow = 195, ncol = 1)
dervDQo=matrix(data=NA, nrow = 195, ncol = 1)
dervDQ_mod=matrix(data=NA, nrow = 195, ncol = 1)
dervD_nyc=matrix(data=NA, nrow = 195, ncol = 1)
dervD_nycMod=matrix(data=NA, nrow = 195, ncol = 1)

TotalQ=matrix(data=NA, nrow = 21, ncol = 1)
cvQin=matrix(data=NA, nrow = 21, ncol = 1)
cvQman=matrix(data=NA, nrow = 21, ncol = 1)
cvQmod=matrix(data=NA, nrow = 21, ncol = 1)

DOY_dx=matrix(data=NA, nrow = 21, ncol = 5)
colnames(DOY_dx)<-c("StorDec", "QoInc", "ModQinc", "NYCQinc", "ModNYCinc")

for (wy in 1:21){
  stor=FC$AF[FC$WY == yrs[wy]]
  Q=FC$Q[FC$WY == yrs[wy]]
  Qo=FC$Qo[FC$WY == yrs[wy]]
  Mod_Q=results[[wy]][,5]
  Modnyc=results[[wy]][,6]
  nyc=FC$NYC_cfs[FC$WY == yrs[wy]]
  TotalQ[wy] = sum(Q)
  
  cvQin[wy]= cv(Q)
  cvQman[wy]= cv(Qo)
  cvQmod[wy]= cv(Mod_Q)
  for(i in 1:195){
    dervDS[i]=(stor[i+1]/stor[i])
    dervDQo[i]=(Qo[i+1]/Qo[i])
    dervDQ_mod[i]=(Mod_Q[i+1]/Mod_Q[i])
    dervD_nyc[i]=(nyc[i+1]/nyc[i])
    dervD_nycMod[i]=(Modnyc[i+1]/Modnyc[i])
  }
  DOY_dx[wy,]=c(max(which(dervDS > 1)), min(which(dervDQo > 1)), min(which(dervDQ_mod > 1)), min(which(dervD_nyc > 1)), min(which(dervD_nycMod > 1)))
  
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
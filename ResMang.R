#import, sum storage for Anderson Ranch, Arrowrock, and Lucky Peak Reservoirs
setwd("~/Documents/GitRepos/ReservoirModeling")
#setwd("D:/Documents/GitRepos/ReservoirModeling")

#--------------------------------------
# Import reservoir data 
#--------------------------------------
# missing data was replaced by average of bracketing values
res<-read.csv("Data/BRB_reservoir_data_1997-2018_noleap.csv")
res$Date <-as.Date(res$Date, format ="%m/%d/%y")
res$totAF<- res$andAF +res$arkAF +res$lucAF
res$qoV<-res$qo*24*60*60*.0000229569 #convert from flow rate (cfs) to volume (acre-feet)
res$resid<- res$in_unreg- res$in_computed
res$Y = as.numeric(format(res$Date, format = "%Y"))
res$M = as.numeric(format(res$Date, format = "%m"))
res$WY[res$M <10]= res$Y[res$M <10]
res$WY[res$M >= 10] = res$Y[res$M >= 10]+1
lowell<-read.csv("Data/lowell_data.csv")
res$low<-lowell$low_af

# load Rule Curve data --------
fv<-read.csv("Data/ForecastVol.csv")
fcVol<-read.csv("Data/FloodControlVol.csv", header=FALSE) #DOY, date, volumes
fcVol[,3:36]<- fcVol[,3:36]*1000 #bc flood control space is in 1000 ac-ft see plate 7-1 in the WCM
fcVol$V2<- as.Date(fcVol$V2, format ="%m/%d/%Y")
prj<-read.csv("Data/Inflw_prj.csv") #coefficients for projection eqn
prjAP<-read.csv("Data/Inflw_prjAPril1.csv") #coefficients for projection eqn after April 1
wfc<-read.csv("Data/plate7-2.csv") #winter flood control space
qlim<- read.csv("Data/MaxQ.csv", header = FALSE) #discharge mins and maxes
#conversions from vol to flow (ac-ft to cfs)
v2f<-43560.000443512/(24*60*60)
f2v<-24*60*60*.0000229569 
#Reservoir Storage from WCM
minQ<-240
maxAF<-1010188
minS<-  41000 +11630+ 28767 #total inactive capacity AND, ARK, LP


#subset (DOY 1: July 31st) from full timeseries -----
jul=196 #change to 212 #july 31st once update fcVol csv
reps <-100
idx<-which(res$doy >= 1 & res$doy <= jul)
rows<- length(idx)
FC<- data.frame(res$doy[idx],res$WY[idx], res$totAF[idx], res$in_unreg[idx], res$qo[idx])
colnames(FC)<-c("doy", "WY", "AF", "Q", "Qo")
yrs<- 1998:2018

#calculate daily forecast values
forecast<-read.csv("Data/LP_coordinatedForecasts.csv")
volF<-vector(length = nrow(FC))
for (i in 1:nrow(FC)){
  if (any(forecast$wy == FC$WY[i] & forecast$doy == FC$doy[i])){
    ii=which(forecast$wy == FC$WY[i] & forecast$doy == FC$doy[i])
    volF[i] <- forecast$ForecastVol[ii]
  } else{volF[i] <- volF[i-1] - FC$Q[i-1]*f2v}
}

FC<-cbind(FC, volF)
FC$volF[FC$volF < 0] <- 0

#find the index of the first doy
doy1<-matrix(data=NA, ncol=1, nrow=21)
for (wy in 1:21){
  doy1[wy]<- which(FC$WY == yrs[wy] & FC$doy == 1)
}

##--------------------------------------
#       DEFINE FUNCTIONS 
##--------------------------------------

# REQUIRED storage - lookup day of year, inflow volume ----
# plate 7-1 and plate 7-3 - Needs sum of timeseries
reqStor<- function(sumQin,doy){ 
 
  fvol<-round(sumQin/1000000, digits=1)
  fcol<- as.numeric(fv$col[fv$fv == fvol])
  
  fcs<- as.numeric(fcVol[doy, fcol+2]) 
  #Winter flood control space (low flow years Plate 7-2) ----
  if (doy < 91 && fvol > 1.2 && fvol < 1.8){
    wcol<- as.numeric(fv$wcol[fv$wfc == fvol])
    fcs<- as.numeric(wfc[doy,wcol])
  }
  if (is.na(fcs) == TRUE){
    fcs<- as.numeric(fcVol[doy, fcol+2])
  }
  #required flood control space (storage volume)

  return(fcs)
} 

#predict max storage in the next m days
predMaxS<- function(m){
  for (day in 1:jul){
    volF <- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow
    count=0
    
    if (volF >= 0 && day < (jul-m)){
      for (it in day:(day+m-1)){
        count=count+1
        vF = volF - (volF/(jul-day))*(day-it) #predict maxS given equal distribution of inflow 
        reqS <- reqStor(vF, it) 
        maxS[day, count] <- (maxAF-reqS) #max storage today given the whole years inflow
      }
    } else {storD <- 0}
  
    maxS[is.na(maxS)]<- maxAF
  }
  return(maxS)
}
#determine minimum daily release before April 1
minRelease<- function(day, volF){
  ix= which(prj$start <= day & prj$end >= day)
  vol1 = volF*prj$b[ix] + prj$c[ix]
  if (prj$start[ix] != day){
    vol2 = volF*prj$b[ix+1] + prj$c[ix+1]
    frac = (day - prj$start[ix])/(prj$end[ix]-prj$start[ix])
    volFmar <-  frac*vol1 + (1-frac)*vol2
  } else {volFmar = vol1}
  
  volFresid <- volF - volFmar
  FCvolAP<- reqStor(volFresid, 91) #flood control space required on april 1
  minEvac<- FCvolAP - availStor[day] #minimum evaculation btw today and April 1
  minReleaseVol <- minEvac+volFmar
  Qmin <- (minReleaseVol*v2f)/(jul-day+1) #associated  qmin
  
  ##If statements that constrain for high flows and ramp rates?
}
#determine minimum daily release after April 1
minReleaseApril<-function(day, volF){
  ix30 = findInterval(day, prjAP$doy)
  ix15=ix30-1
  volF_target15 <- volF*prjAP$b[ix15] + prjAP$c[ix15]
  volF_target30 <- volF*prjAP$b[ix30] + prjAP$c[ix30]
  residual15<- volF-volF_target15
  residual30<- volF- volF_target30
  
  
  if (day < jul-30){
    FCvol30<- reqStor(residual30, day+30) #required storage space
    minEvac30 <- FCvol30 - availStor[day]
    minReleaseVol30 <- minEvac30 + volF_target30
    q30<-(minReleaseVol30*v2f)/(jul-day+1) 
  } else {q30 <- minQ}
  if (day < jul-15){
    FCvol15<-reqStor(residual15, day+15)
    minEvac15 <- FCvol15 - availStor[day]
    minReleaseVol15 <- minEvac15 + volF_target15
    q15<-(minReleaseVol15*v2f)/(jul-day+1)
  } else {q15 <- minQ}
  
  Qmin <- max(q15, q30)
  
}

#forecast what the storage would be in s days given previous ∆ in S and make those changes over m days
forecastS<-function(s,m,day){
  if (day > s+1){
    dsdt= (stor[day] - stor[day-s])/s
    storF[day] <<- (dsdt*m)+stor[day] 
  } 
}  

#evaluate change in storage in regard to the forecasted storage to prevent going over maxS
#update discharge, change in storage and day+1 storage
evalS<- function(Qin, day, stor, maxS, Qmin,s,m){
  
  if (storF[day] >= maxS[day,m] && day > s+1){ #&& day <188
    dsdtMax= (storF[day] - maxS[day,m])/s
    qo[day] <- Qmin[day] + (dsdtMax*v2f)
    flag = 'TRUE' #true we need to increase ramp rates to get rid of the water
  } else {qo[day] <- Qmin[day]
  flag='FALSE'}
  
  if (stor[day] > maxAF){
    addQ = (stor[day] - maxAF)*v2f
    minFCq[day] <- minFCq[day] + addQ
    flag= 'TRUE'
  } else {flag='FALSE'}
  
  #if the calculated discharge is greater than +/- 500 set it to +/- 500
  if (flag == 'TRUE'){ #going over maxS or maxAF
    ramp= 1000
  } else {ramp = 500}
  
  if (qo[day] > (qo[day-1] + ramp) && day > 2 ){
    qo[day] <- qo[day-1] + ramp
  } else if (qo[day] < (qo[day-1]-500) && day > 2 ){
    qo[day] <- qo[day-1]-500}
  
  if (availStor[day] <= (1000*v2f)){ #this puts a hard constraints on not topping the dam - but doesnt work 21 times
    qo[day] <- qo[day] + 1000
    availStor[day] <- availStor[day] + 1000*f2v
    stor[day] <- stor[day] - (1000*f2v)
  }
  
  dS[day] <- (Qin[day]- qo[day])*f2v
  
  if (stor[day] <= minS){ #dont let storage go below the minimum
    qo[day] <- minQ
    dS[day] <- -minQ*f2v 
  }
  
  if (day < jul){
    stor[day+1]<<-stor[day] + dS[day]  #AF in the reservoir
  }
  
  qo[day]<<-qo[day]
  dS[day]<<-dS[day]
}

#determine change in storage and outflow for a given water year and forecast window
outflowStor<-function(wy,s,m){
  #   set up blank matricies 
  #------------------------------------------
  stor<<-matrix(data=NA, nrow = jul, ncol = 1)
  maxS<<-matrix(data=NA, nrow = jul, ncol = m) 
  availStor<<-matrix(data=NA, nrow = jul, ncol = 1)
  minFCq<<-matrix(data=NA, nrow = jul, ncol = 1)
  storF<<-matrix(data=NA, nrow = jul, ncol = 1)
  qo<<-matrix(data=NA, nrow = jul, ncol = 1) #modeled outflow from reservoir
  dS<<-matrix(data=NA, nrow = jul, ncol = 1)
  resS<<-matrix(data=NA, nrow = jul, ncol = 1)
  #---------
  # initalize
  #----- 
  stor[1] <<- FC$AF[doy1[wy]] #initialize with actual storage on Jan 1
  Qin<- FC$Q[FC$WY == yrs[wy]]
  maxS <<- predMaxS(m) #vector of 198 days of max storage out to M days
  
  #----- run all the functions to get to discharge and updated storage
  for (day in 1:jul){ 
    volF<<- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow  
    availStor[day] <<- maxAF-stor[day]
    # Min flood control release for storage goals on April 1 and every 15 days after
    if (day < 91){
      minFCq[day]<<-minRelease(day, volF)
    } else {
      minFCq[day]<<-minReleaseApril(day, volF)
    }
    #minimum discharge is 240
    if (minFCq[day] < minQ){
      minFCq[day] <<- minQ
    } 
    
    forecastS(s,m,day)
    evalS(Qin, day, stor, maxS, minFCq, s, m)
  }
  
  out<<- cbind(availStor, storF, stor, dS, minFCq, qo)
  colnames(out)<-c("availStor", "storF", 'stor', 'dS', 'minFCq', 'qo')
  return(out)
}


#---------------------------------------------------------------
#   determine change in storage and outflow for any day of year
#---------------------------------------------------------------
#select Qin from matrix or array? how to save output better?
#s: number of prior days to estimate change in storage
#m: planning window
m=10
s=5
wy=2
results<-list()

params<-cbind(c(2,2,2,2,2,2,2,2,2,2),c(1,2,3,4,5,6,7,8,9,10), c(1,2,3,4,5,6,7,8,9,10))

for (wy in 1:21){
  results[[wy]]<- outflowStor(wy,5,10)
}

modelRun<-function(params){
  return(mapply(outflowStor, params[,1], params[,2], params[,3]))
}


#plot the initial results
for (wy in 1:21){
  plot(FC$maxS[FC$WY == yrs[wy]], type='l', ylim=c(300000, 1010200))
  lines(FC$stor[FC$WY == yrs[wy]], col='orange')
  lines(FC$AF[FC$WY == yrs[wy]], col='green') 
 
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(FC$qo[FC$WY == yrs[wy]], type='l', lty=3, col='orange', ylim=c(0,16000))
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
}



#make a few plots - e.g number, timing and volume that the managers exceeded the maxS for the day - and other plots from the nicholas data site
#determine average day they start drafting for irrigation - e.g. once dsdt is negative and never goes positive again

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
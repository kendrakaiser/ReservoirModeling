#import, sum storage for Anderson Ranch, Arrowrock, and Lucky Peak Reservoirs
setwd("~/Documents/GitRepos/ReservoirModeling")
#setwd("D:/Documents/GitRepos/ReservoirModel")

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

# REQUIRED storage - lookup day of year, inflow volume 
# plate 7-1 and plate 7-3 - Needs sum of timeseries
resStor<- function(sumQin,doy){ 
 
  fvol<-round(sumQin/1000000, digits=1)
  fcol<- as.numeric(fv$col[fv$fv == fvol])
  stor<- as.numeric(fcVol[doy, fcol+2])
  #Winter flood control space (low flow years Plate 7-2) ----
  if (doy < 91 && fvol > 1.2 && fvol < 1.8){
    wcol<- as.numeric(fv$wcol[fv$wfc == fvol])
    stor<- as.numeric(wfc[doy,wcol])
  }
  if (is.na(stor) == TRUE){
    stor<- as.numeric(fcVol[doy, fcol+2])
  }
  #required storage volume
  outList<-(list("stor" = stor))
  return(outList)
} 

#determine minimum daily release before April 1
minRelease<- function(){
  ix= which(prj$start <= day & prj$end >= day)
  vol1 = volF*prj$b[ix] + prj$c[ix]
  if (prj$start[ix] != day){
    vol2 = volF*prj$b[ix+1] + prj$c[ix+1]
    frac = (day - prj$start[ix])/(prj$end[ix]-prj$start[ix])
    volFmar <-  frac*vol1 + (1-frac)*vol2
  } else {volFmar = vol1}
  
  volFresid <- volF - volFmar
  FCvolAP<- resStor(volFresid, 91) #flood control space required on april 1
  minEvac<- FCvolAP$stor - stor[day] #minimum evaculation btw today and April 1
  minReleaseVol <- minEvac+volFmar
  Qmin<- (minReleaseVol*v2f)/(jul-day+1) #associated  qmin
  
  ##If statements that constrain for high flows and ramp rates
}

#determine minimum daily release after April 1
minReleaseApril<-function(){
  ix30 = findInterval(day, prjAP$doy)
  ix15=ix30-1
  volF_target15 <- volF*prjAP$b[ix15] + prjAP$c[ix15]
  volF_target30 <- volF*prjAP$b[ix30] + prjAP$c[ix30]
  residual15<- volF-volF_target15
  residual30<- volF- volF_target30
  
  if (day < jul-30){
    FCvol30<- resStor(residual30, day+30)
    minEvac30 <- FCvol30$stor - stor[day]
    minReleaseVol30 <- minEvac30 + volF_target30
    q30<-(minReleaseVol30*v2f)/(jul-day+1) 
  } else {q30 <- minQ}
  if (day < jul-15){
    FCvol15<-resStor(residual15, day+15)
    minEvac15 <- FCvol15$stor - stor[day]
    minReleaseVol15 <- minEvac15 + volF_target15
    q15<-(minReleaseVol15*v2f)/(jul-day+1)
  } else {q15 <- minQ}
  
  Qmin<- max(q15, q30)
  
}

# CHANGE STORAGE - based on current storage, minimum discharge and ramping rates
#needs whole timeseries and maxS of any given day 
changeS<- function(Qin, day, stor, maxS, Qmin){ #consider if these all need to be here

  if (stor[day] >= maxS) {#if S > max AF, change of volume in reservoir
    dS[day] <- maxS - stor[day] - (Qin[day]*f2v) #calculate ∆ volume of water in the reservoir
    qo[day] <- -dS[day]*v2f
    ##******************************
    #this is the problem - moved ramping rates here, but they need to be contingent on the previous days discharge and how close to the maxS 
    ##******************************
  #} else if (day > 21 && minFCq[day] > minFCq [day-1]) {
    #qo[day] <- qo[day-1] + 500 
    #dS[day] <- (Qin[day]- qo[day])*f2v
  #} else if (day > 21 && minFCq[day] <= minFCq[day-1]) {
   # qo[day] <- qo[day-1] - 500 
    #dS[day] <- (Qin[day]- qo[day])*f2v
  } else if (day > 1) {
    qo[day] <- qo[day-1] 
    dS[day] <- (Qin[day]- qo[day])*f2v
  } else {
    qo[day] <- minQ 
    dS[day] <- (Qin[day]- qo[day])*f2v
  }

  if (stor[day] <= minS){ #dont let storage go below the minimum
    qo[day] <- minQ
    dS[day] <- -minQ*f2v 
  }
  if (day < jul){
    stor[day+1]<-stor[day] + dS[day]  #AF in the reservoir
  }
  
  outlist<-(list("stor"=stor, "dS"=dS, "qo"=qo, "day"=day))
  return(outlist)
}


#Ramping rate is +/- 500 cfs per day  --> distributes the water over following days
ramprate <- function(qo, stor){

  for (day in 1:jul){
    if (day > 10 && qo[day] < (qo[day-1] - 500)){ #todays qo < 500cfs than yesterdays
      Q=qo[day]
      qo[day] <- qo[day-1]-500
      dS[day] <- Qin[day]- (qo[day]*f2v)
      Qd<-Q-qo[day]
      dD <- day + round((Qd)/500) 
        if (dD > jul) {dD=jul}
      qo[(day+1):dD] <- qo[(day+1):dD]+ Qd/(round((Qd)/500))
      dS[(day+1):dD] <- (Qin[(day+1):dD]-qo[(day+1):dD])*f2v 
      stor[day:dD] <- stor[day:dD] + dS[day:dD]
  }else if (day >10 && qo[day] > (qo[day-1] + 500)){
      #todays qo > 500cfs than yesterday
    #kinda works - line by line, but something is wrong w looping or something
      Q=qo[day]
      qo[day] <- qo[day-1]+500
      dS[day] <- Qin[day]- (qo[day]*f2v)
      Qd<-Q-qo[day]
      dD <- day + round((Qd)/500) 
        if (dD > jul) {dD=jul}
      qo[(day+1):dD] <- qo[(day+1):dD]+ Qd/(round((Qd)/500))
      dS[(day+1):dD] <- (Qin[(day+1):dD]-qo[(day+1):dD])*f2v 
      stor[day:dD] <- stor[day:dD] + dS[day:dD]
    } else{qo[day]<-qo[day]}
  }
  outlist<-(list("dS"=dS, "qo"=qo, "stor"=stor))
  return(outlist) 
}
#-------------------------------------------------------------
#   set up blank matricies 
#------------------------------------------------------------
resS=matrix(data=NA, nrow = jul, ncol = 1)
stor=matrix(data=NA, nrow = jul, ncol = 1)
maxS=matrix(data=NA, nrow = jul, ncol = 1)
dS=matrix(data=NA, nrow = jul, ncol = 1)
minFCq=matrix(data=NA, nrow = jul, ncol = 1)
qo=matrix(data=NA, nrow = jul, ncol = 1)


#---------------------------------------------------------------
#   determine change in storage and outflow for any day of year
#---------------------------------------------------------------
#select Qin from matrix or array? - turn into funtion

for (wy in 1){
  stor[1]<-FC$AF[doy1[wy]] #initialize with actual storage on Jan 1
  Qin<- FC$Q[FC$WY == yrs[wy]]

  for (day in 1:jul){
    volF= FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow
    if (volF >= 0){
      maxSday<- resStor(volF, day) 
      storD<-maxSday$stor
    } else {storD <- 0}
    
    maxS[day] <- maxAF-storD #max storage today given the whole years inflow
    
    #Qmin for storage goals on April 1 and every 15 days after
    if (day < 91){
      minFCq[day]<-minRelease()
    } else {
      minFCq[day]<-minReleaseApril()
    }
    #minimum discharge is 240
    if (minFCq[day] <= minQ){
      minFCq[day] <- minQ
    }
    
    #min discharge today given forecasted storage in the next n days
    #add variables for surplus storage and associated Qmin 
    
    
    resS <- changeS(Qin, day, stor, maxS[day], minFCq)
    qo[day]<-resS$qo[day]
    dS[day]<- resS$dS[day]
    
    if (day < jul){
      stor[day+1]<-resS$stor[day+1] #AF in the reservoir
    }
  }

  #rr=ramprate(qo, stor)
  
  #save intial run output ----  # change for debugging - put clear matricies at beginning
  FC$qo[FC$WY == yrs[wy]]<-qo[,]   #rr$qo
  FC$stor[FC$WY == yrs[wy]]<-stor[,]   #rr#stor
  FC$maxS[FC$WY == yrs[wy]]<-maxS[,]
  FC$Qmin[FC$WY == yrs[wy]]<-minFCq[,]
}


#plot the initial results
for (wy in 1){
  plot(FC$maxS[FC$WY == yrs[wy]], type='l', ylim=c(300000, 1010200))
  lines(FC$stor[FC$WY == yrs[wy]], col='orange')
  lines(FC$AF[FC$WY == yrs[wy]], col='green') 
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(FC$qo[FC$WY == yrs[wy]], type='l', lty=3, col='orange', ylim=c(0,16000))
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
}

#------------------------------
#enter iterative loop to impose discahrge limits and decrease daily variability until it works #------------------------------
  
spd=7 #number of days to average values over 

for (wy in 1:21){ #this is all just a mess o stuf right now...
  counter=0
  print(wy)
  
  cy<-FC[FC$WY == yrs[wy],]
  qqlim=qlim[1:196,2]
  while ((any(cy$qo > qlim[1:196,2])) && counter < 1000 | ((any(cy$stor > cy$maxS)) && counter < 1000)) {
    counter= counter+1
    print(counter)
    
  #this just forces any day where Q goes over the limit to be at the limit and updates storage 
    cy<-FC[FC$WY == yrs[wy],]
    ind<-which(cy$qo > qqlim)
    id<-ind[1]
    #qDist<-cy$qo[ind] - qqlim[ind]
    indL<-which(cy$qo< qqlim)
    for (i in id:length(indL)){
     if (cy$qo[i] > qqlim[i]){
        cy$qo[id] = cy$qo[id] - (qqlim[indL[i]]-cy$qo[indL[i]])
        cy$stor[indL[i]]= cy$stor[indL[i]] - (qqlim[indL[i]]-cy$qo[indL[i]])*f2v
        cy$qo[indL[i]] = qqlim[indL[i]]
        
    }}
    
    plot(cy$qo[cy$WY == yrs[wy]], type='l', lty=3, col='orange', ylim=c(0,16000))
    lines(qlim[,2], type='l', lty=3, col='grey17')
    lines(cy$Qo[cy$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
    
    
    cy$qo[ind] = qqlim[ind]
    dS[ind] <- (cy$Q[ind] - cy$qo[ind]) * f2v 
    cy$stor[ind]<- cy$stor[ind-1] + dS[ind] 
    
    
    
   #this is one way to keep it from dropping way below storage min, but more of a bandaid 
    ifx<- which(cy$stor <= minS)
    dS[ifx] <- (cy$Q[ifx]- minQ)*f2v
    cy$qoN[ifx]<- minQ 
    
    #trying to allow for a faster ramping rate if the storage is going to go over maxS 
    if (any(cy$stor > cy$maxS) && counter <2){
      ind<-which(cy$stor > cy$maxS)
      dS[ind]<- dS[ind] - (500*f2v)
      cy$qo[ind]<- cy$qo[ind] + 500
      cy$stor[ind]<-cy$stor[ind] - 500*f2v
    }
    
    
    #this was the first attempt at smoothing - by finding all q>qMax and distributing those flows over n prior days and update the dS and stor
    #but something is wrong in here - discharge gets smoothed but storage is unreasonable - something about indexing
    if (any(cy$qo > qlim[1:196,2])){ 
      id<-min(ind)
      if (id > spd+1){
        cy$qo[((id-spd):id)] <- sum(cy$qo[((id-spd):id)]) / (1+spd)
        dS[((id-spd):id)] <- (cy$Q[((id-spd):id)] - cy$qo[((id-spd):id)]) * f2v 
        cy$stor[((id-spd):id)]<- cy$stor[((id-spd):id)] + dS[((id-spd):id)]  
        
      } else if (id < spd){ #trying to deal with the first week of Jan
        cy$qo[(id:(id+(spd*2)))] <- sum(cy$qo[(id:(id+(spd*2)))]) / (1+(spd*2))
        dS[(id:(id+(spd*2)))] <- (cy$Q[(id:(id+(spd*2)))] - cy$qo[(id:(id+(spd*2)))]) * f2v
        cy$stor[(id:(id+(spd*2)))]<- cy$stor[(id:(id+(spd*2)))] + dS[(id:(id+(spd*2)))]  
      }
    }
  }
  
  FC$qoN[FC$WY == yrs[wy]]<- cy$qo
  FC$storN[FC$WY == yrs[wy]]<- cy$stor
  
}


datelab<-seq(as.Date("1997-01-01"), as.Date("1997-07-31"), by="1 day")

for (wy in 1:21){
  plot(FC$maxS[FC$WY == yrs[wy]], type='l', ylim=c(300000, 1010200))
  lines(FC$storN[FC$WY == yrs[wy]], col='orange')
  lines(FC$AF[FC$WY == yrs[wy]], col='green') 
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(FC$qoN[FC$WY == yrs[wy]], type='l', lty=3, col='orange', ylim=c(0,16000))
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
}


for (wy in 20){
  plot(cy$maxS[cy$WY == yrs[wy]], type='l', ylim=c(300000, 1010200))
  lines(cy$stor[cy$WY == yrs[wy]], col='orange')
  lines(cy$AF[cy$WY == yrs[wy]], col='green') 
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(cy$qo[cy$WY == yrs[wy]], type='l', lty=3, col='orange', ylim=c(0,16000))
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(cy$Qo[cy$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow
}
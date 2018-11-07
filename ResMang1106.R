#import, sum storage for Anderson Ranch, Arrowrock, and Lucky Peak Reservoirs
#setwd("~/Dropbox/BSU/R/WaterActors/ToyModels/ReservoirModeling")
setwd("D:/Dropbox/BSU/R/WaterActors/ToyModels/ReservoirModeling")

# Import reservoir data ------
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
forecast<-read.csv("Data/LP_coordinatedForecasts.csv")

#run for DOY 1: July 31st subset from full timeseries -----
jul=196 #change to 212 #july 31st once update fcVol csv
reps <-100
idx<-which(res$doy >= 1 & res$doy <= jul)
rows<- length(idx)
FC<- data.frame(res$doy[idx],res$WY[idx], res$totAF[idx], res$in_unreg[idx], res$qo[idx])
colnames(FC)<-c("doy", "WY", "AF", "Q", "Qo")
yrs<- 1998:2018

# load Rule Curve data --------
fv<-read.csv("Data/ForecastVol.csv")
fcVol<-read.csv("Data/FloodControlVol.csv", header=FALSE) #DOY, date, volumes
fcVol[,3:36]<- fcVol[,3:36]*1000 #bc flood control space is in 1000 ac-ft see plate 7-1 in the WCM
fcVol$V2<- as.Date(fcVol$V2, format ="%m/%d/%Y")
prj<-read.csv("Data/Inflw_prj.csv") #coefficients for projection eqn
wfc<-read.csv("Data/plate7-2.csv") #winter flood control space
qlim<- read.csv("Data/MaxQ.csv", header = FALSE) #discharge mins and maxes
Qpred<-read.csv("Data/Qpred.csv", header= FALSE)
#conversions from vol to flow (ac-ft to cfs)
v2f<-43560.000443512/(24*60*60)
f2v<-24*60*60*.0000229569 
#Reservoir Storage from WCM
minQ<-240
maxAF<-1010188
minS<-  41000 +11630+ 28767 #total inactive capacity AND, ARK, LP

#find the index of the first doy
doy1<-matrix(data=NA, ncol=1, nrow=21)
for (wy in 1:21){
  doy1[wy]<- which(FC$WY == yrs[wy] & FC$doy == 1)
}

#REQUIRED storage - lookup day of year, inflow volume ####
### plate 7-1 and plate 7-3 ### ----
#Needs sum of timeseries
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

#needs whole timeseries and maxS of any given day 
changeS<- function(Qin, day, stor, maxS, Qmin){ #consider if these all need to be here

  if (stor[day] >= maxS) {#if S > max AF, change of volume in reservoir
    dS[day] <- maxS - stor[day] - (Qin[day]*f2v) #calculate change of volume of water in the reservoir
    qo[day] <- -dS[day]*v2f
  } else {
    qo[day] <- Qmin ##this is tricky bc then the discharge will always be high
    dS[day] <- (Qin[day]- Qmin)*f2v
  } 
  #Ramping rate is +/- 500 cfs per day
  if (day > 1 && qo[day] < qo[day-1] - 500){
    qo[day] <- qo[day-1] - 500
    dS[day] <- (Qin[day]-qo[day])*f2v 
  } else if (day >1 && qo[day] > qo[day-1] + 500){
    qo[day] <- qo[day-1] + 500
    dS[day] <- (Qin[day]-qo[day])*f2v 
  } else{}
  if (stor[day] <= minS){
    qo[day] <- minQ
    dS[day] <- minQ*f2v 
  }
  if (day < jul){
    stor[day+1]<-stor[day] + dS[day]  #AF in the reservoir
  }
  
  outlist<-(list("stor"=stor, "dS"=dS, "qo"=qo, "day"=day))
  return(outlist)
}

#### set up blank matricies -----
resS=matrix(data=NA, nrow = jul, ncol = 1)
stor=matrix(data=NA, nrow = jul, ncol = 1)
maxS=matrix(data=NA, nrow = jul, ncol = 1)
dS=matrix(data=NA, nrow = jul, ncol = 1)
QminAP=matrix(data=NA, nrow = jul, ncol = 1)
qo=matrix(data=NA, nrow = jul, ncol = 1)

#select Qin from matrix or array? - turn into funtion

for (wy in 1:21){
  stor[1]<-FC$AF[doy1[wy]]
  Qin<- FC$Q[FC$WY == yrs[wy]]
#### determine reservoir storage and discharge for any given day of year
  for (day in 1:jul){
    # use Coordinated forecast on 1st of every month
    if (any(forecast$doy == day && forecast$wy == yrs[wy])){
      volF <- forecast[(forecast$doy == day & forecast$wy == yrs[wy]),1]
    }
    else{volF <- sum(Qin[day:jul])*f2v} #otherwise use the actual inflow data - UPDATE this w daily forecast timeseries
    
    maxSday<- resStor(volF, day) 
    maxS[day] <- maxAF-maxSday$stor #max storage today given the whole years inflow
    
    #Determine April 1 FC space and Qmin
    if (any(prj$start == day)){
      ix= which(any(prj$start == day))
      volFmar<- volF*prj$b[ix] + prj$c[ix]
      volFresid <- volF - volFmar
    
      FCvolAP<- resStor(volFresid, 91) #flood control space required on april 1
      APmaxS<- maxAF- FCvolAP$stor #max storage on april 1
      QminAP[day]<- (volF - FCvolAP$stor + volFmar)*v2f/(jul-day+1) #associated  qmin
    } else {QminAP[day] <- minQ}
    
    resS <- changeS(Qin, day, stor, maxS[day], QminAP[day])
    qo[day]<-resS$qo[day]
    
    if (day < jul){
      stor[day+1]<-resS$stor[day+1] #AF in the reservoir
    }
  }

  #save intial run output ----
  FC$qo[FC$WY == yrs[wy]]<-qo[,]
  FC$stor[FC$WY == yrs[wy]]<-stor[,]
  FC$maxS[FC$WY == yrs[wy]]<-maxS[,]
  FC$Qmin[FC$WY == yrs[wy]]<-QminAP[,]
  
  plot(maxS, type='l', ylim=c(300000, 1010200))
  lines(stor, col='orange')
  lines(FC$AF[FC$WY == yrs[wy]], col='green') 
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(qo, type='l', lty=3, col='orange', ylim=c(0,16000))
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow

}

#enter iterative loop until it works ----
spd=15 #number of days to average values over 
#bumping this up this high definitely helped mediate high flows - but the storage values are all over the place

for (wy in 1:21){
  counter=0
  print(wy)
  
  cy<-FC[FC$WY == yrs[wy],]
  
  while ((any(cy$qo > qlim[1:196,2])) && counter < 200 | ((any(cy$stor > cy$maxS)) && counter < 200)) {
    counter= counter+1
    print(counter)
    wy
    ifx<- which(cy$stor <= cy$minS)
    dS[ifx] <- (cy$Q[ifx]- minQ)*f2v
    cy$qo[ifx]<- minQ 
    
    #trying to allow for a faster ramping rate if the storage is going to go over maxS - need to do it better tho
    if (any(cy$stor > cy$maxS) && counter <2){
      ind<-which(cy$stor > cy$maxS)
      dS[ind]<- dS[ind] - (500*f2v)
      cy$qo[ind]<- cy$qo[ind] + 500
      cy$stor[ind]<-cy$stor[ind] - 500*f2v
    }
    
    if (any(cy$qo > qlim[1:196,2])){ #find all q>qMax distribute those flows over prior days and update the dS and stor
      ind<-which(cy$qo > qlim[1:196,2])
      id<-min(ind)
      
      if (id > spd+1){
        cy$qo[((id-spd):id)] <- sum(cy$qo[((id-spd):id)]) / (1+spd)
        dS[((id-spd):id)] <- (cy$Q[((id-spd):id)] - cy$qo[((id-spd):id)]) * f2v 
        cy$stor[((id-spd):id)]<- cy$stor[((id-spd):id)] + dS[((id-spd):id)]  
      } else if (id < spd){
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
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

#Relationship between change in discharge and flow -----
flow=seq(from=240,to=10740, by=450)
flow2=seq(740, 11240, 450)
perc= flow/flow2
p_dQ<-lm(perc~ log(flow))
a<- as.numeric(p_dQ$coefficients[1])
b<- as.numeric(p_dQ$coefficients[2])

p<- function(flow){
  a + log(flow)*b
}

# REQUIRED storage - lookup day of year, inflow volume ----
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

#predict max storage in the next m days
predMaxS<- function(){
  for (day in 1:jul){
    volF<<- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow
    count=0
    if (volF >= 0 && day< jul){
      for (it in day:(day+m-1)){
        count=count+1
        vF = volF - (volF/(jul-day))*(day-it) #predict maxS given equal distribution of inflow over next m days
        maxSday<- resStor(vF, it) 
        storD<-maxSday$stor
        maxS[day, count] <- (maxAF-storD) #max storage today given the whole years inflow
      }
    } else {storD <- 0}
  
    maxS[is.na(maxS)]<- maxAF
  }
  return(maxS)
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
  
  ##If statements that constrain for high flows and ramp rates?
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

#evaluate change in storage to prevent going over maxS
evalS<- function(Qin, day, stor, maxS, Qmin, n){
  
  if (day > n+1){
  dsdt= (stor[day] - stor[day-n])/n
  storF[day] = (dsdt*n)+stor[day] #forecast what the storage would be in n days given previous ∆ in S
  
    if (storF[day] >= maxS[day,n] && day <188){
      dsdtMax= (storF[day] - maxS[day,n])/n
      Qnew = Qmin[day] + (dsdtMax*v2f)
    } else {Qnew = Qmin[day]}

  qo[day] <- Qnew
  
  } else {
    qo[day] <- Qmin[day]
  }
  
  #if the calculated discharge is greater than +/- 500 set it to +/- 500
  if (qo[day] > qo[day-1]+500 && day > 2){
    qo[day] = qo[day-1]+500
  } else if (qo[day] < qo[day-1]-500 && day > 2){
    qo[day] = qo[day-1]-500}
  
  #if S > maxS, change of volume in reservoir
  #if (stor[day] >= maxS) {
   # qo[day] <- qo[day-1] + 500}
  
  dS[day] <- (Qin[day]- qo[day])*f2v
  
  if (stor[day] <= minS){ #dont let storage go below the minimum
    qo[day] <- minQ
    dS[day] <- -minQ*f2v 
  }
  
  if (day < jul){
    stor[day+1]<-stor[day] + dS[day]  #AF in the reservoir
  }
  
  outlist<-(list("stor"=stor, "dS"=dS, "qo"=qo, "storF"=storF))
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
  } else {
    qo[day] <- Qmin[day]
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

s=5 #number of prior days to estimate change in storage
m=10 #planning window
#-------------------------------------------------------------
#   set up blank matricies 
#------------------------------------------------------------
resS=matrix(data=NA, nrow = jul, ncol = 1)
stor=matrix(data=NA, nrow = jul, ncol = 1)
maxS=matrix(data=NA, nrow = jul, ncol = m)
dS=matrix(data=NA, nrow = jul, ncol = 1)
minFCq=matrix(data=NA, nrow = jul, ncol = 1)
qo=matrix(data=NA, nrow = jul, ncol = 1)
storF=matrix(data=NA, nrow = jul, ncol = 1)


#---------------------------------------------------------------
#   determine change in storage and outflow for any day of year
#---------------------------------------------------------------
#select Qin from matrix or array? - turn into funtion

#for (wy in 2){
wy=2  
  stor[1]<-FC$AF[doy1[wy]] #initialize with actual storage on Jan 1
  Qin<- FC$Q[FC$WY == yrs[wy]]

  maxS<-predMaxS()



  
  for (day in 1:jul){ 
    # Min flood control release for storage goals on April 1 and every 15 days after
    if (day < 91){
      minFCq[day]<-minRelease()
    } else {
      minFCq[day]<-minReleaseApril()
    }
    #minimum discharge is 240
    if (minFCq[day] < minQ){
      minFCq[day] <- minQ
    }
    
    
    resS <- evalS(Qin, day, stor, maxS, minFCq, s)
    qo[day]<-resS$qo[day]
    dS[day]<- resS$dS[day]
    storF[day]<-resS$storF[day]
    
    if (day < jul){
      stor[day+1]<-resS$stor[day+1] #AF in the reservoir
    }
  }


  
  #save intial run output ----  # change for debugging - put clear matricies at beginning
  FC$qo[FC$WY == yrs[wy]]<-qo[,]   #rr$qo
  FC$stor[FC$WY == yrs[wy]]<-stor[,]   #rr#stor
  FC$maxS[FC$WY == yrs[wy]]<-maxS[,1]
  FC$Qmin[FC$WY == yrs[wy]]<-minFCq[,]
  FC$storF[FC$WY == yrs[wy]]<-storF[,]
  
#}



#plot the initial results
wy=2
  plot(FC$maxS[FC$WY == yrs[wy]], type='l', ylim=c(300000, 1010200))
  lines(FC$stor[FC$WY == yrs[wy]], col='orange')
  lines(FC$AF[FC$WY == yrs[wy]], col='green') 
  lines(FC$storF[FC$WY == yrs[wy]], col='blue') 
  lines(maxS[,m], col="red")
  
  #plot(Qmin, type='l', col='blue', ylim=c(0,16000))
  plot(FC$qo[FC$WY == yrs[wy]], type='l', lty=3, col='orange', ylim=c(0,16000))
  lines(qlim[,2], type='l', lty=3, col='grey17')
  lines(FC$Qo[FC$WY == yrs[wy]], type='l', lty=5, lwd='1', col='skyblue1') #manged outflow



## Model of Reservoir Storage and Outflow in the Boise River Basin
## Code based of of the ACOE Water Control Manual
## Kendra E Kaiser

setwd("~/Documents/GitRepos/ReservoirModeling")

#--------------------------------------
# Import reservoir data 
#--------------------------------------
#sum storage for Anderson Ranch, Arrowrock, and Lucky Peak Reservoirs
# missing data was replaced by average of bracketing values
res<-read.csv("Data/BRB_reservoir_data_1997-2018_noleap.csv")
res$Date <-as.Date(res$Date, format ="%m/%d/%y")
res$totAF<- res$andAF +res$arkAF +res$lucAF
res$qoV<-res$qo*24*60*60*.0000229569 #convert from flow rate (cfs) to volume (acre-feet)
res$resid<- res$in_unreg- res$in_computed #qu is the correct unregulated flow estimate
res$Y = as.numeric(format(res$Date, format = "%Y"))
res$M = as.numeric(format(res$Date, format = "%m"))
res$WY[res$M <10]= res$Y[res$M <10]
res$WY[res$M >= 10] = res$Y[res$M >= 10]+1
lowell<-read.csv("Data/lowell_data.csv") #fb = reservoir water surface elevation in ft
res$low<-lowell$low_af
nyc<-read.csv("Data/NY_canal_data.csv")
res$NYC<-nyc$bsei_qj[1:7646]

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
FC<- data.frame(res$doy[idx],res$WY[idx], res$totAF[idx], res$in_unreg[idx], res$qo[idx], res$low[idx], res$NYC[idx])
colnames(FC)<-c("doy", "WY", "AF", "Q", "Qo", "LowellAF", "NYC_cfs")
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

##INITALIZE storF with the intial storage of the reservoir

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
  
  #if march 1st check lowell volume and if < full fill at whatever rate every day till full
  if (day<51){
      Lowell_cfs<<-0
} else if (day >= 40 && day < 91 && LowellAF[day-1] < 155237){
      low_under <- 155237 - lowell$low_af[day-1] 
      #calc how much under maximum storage the lake is
      Lowell_flow_est <- low_under/(91-day) * v2f #days between start fill date (Feb21st) to March 31st
      #if the min flow isnt enough to slowly fill Lowell then increase it 
      if (Qmin < Lowell_flow_est){ 
        Qmin <- Lowell_flow_est
      } else{
        Lowell_flow_est <- Qmin
      }
      
      Lowell_cfs[day] <<- Lowell_flow_est
      LowellAF[day] <<- lowell$low_af[day-1] + (Lowell_flow_est * f2v)
  
  } else if (day >=91){
    LowellAF[day]<<-LowellAF[day-1]
    Lowell_cfs<<-0}

  return(Qmin)
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
    dsdt= (stor[day] - stor[day-s])/s #calculate daily change in storage over past s days --> change this to 1) take average of dS over s days, 2) create weighted average
    storF[day] <<- (dsdt*m)+stor[day] #forecasted flow for today is the change in storage* m days subtracted form todays storage -- does that make sense??
  } else {storF[day]<<- stor[day]}
}  
#evaluate change in storage in regard to the forecasted storage to prevent going over maxS
#update discharge, change in storage and day+1 storage
evalS<- function(Qin, day, stor, maxS, Qmin,s,m){
  if (storF[day] >= maxS[day,m] && day > s+1){ #&& day <188 #this is max storage m days out
    dsdtMax= (storF[day] - maxS[day,m])/m #max change in storage is m days out divided by m
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
  ###UPDATE here
  #this puts a hard constraints on not topping the dam - but doesnt work 21 times
  if (availStor[day] <= (1000*v2f)){ 
    qo[day] <- qo[day] + 1000
    availStor[day] <- availStor[day] + 1000*f2v
    stor[day] <- stor[day] - (1000*f2v)
  }
  
  dS[day] <- (Qin[day]- qo[day])*f2v
  #dont let storage go below the minimum
  if (stor[day] <= minS){ 
    qo[day] <- minQ
    dS[day] <- -minQ*f2v 
  }
  #update storage for the next day
  if (day < jul){
    stor[day+1]<<-stor[day] + dS[day]  #AF in the reservoir
  }
  #write out values
  qo[day]<<-qo[day]
  dS[day]<<-dS[day]
}

#determine change in storage and outflow for a given water year and forecast window
outflowStor<-function(s,m){
  results<-list()
  discharge<-matrix(data=NA, nrow = jul, ncol=21)
  for (wy in 1:21){
    #   set up blank matricies 
    #------------------------------------------
    stor<<-matrix(data=NA, nrow = jul, ncol = 1)
    maxS<<-matrix(data=NA, nrow = jul, ncol = m) 
    availStor<<-matrix(data=NA, nrow = jul, ncol = 1)
    minFCq<<-matrix(data=NA, nrow = jul, ncol = 1)
    storF<<-matrix(data=NA, nrow = jul, ncol = 1)
    qo<<-matrix(data=NA, nrow = jul, ncol = 1) #modeled outflow from reservoir
    dS<<-matrix(data=NA, nrow = jul, ncol = 1)
    LowellAF<<-matrix(data=NA, nrow=jul, ncol = 1) #modeled Lowell Acre-feet
    Lowell_cfs<<-matrix(data=NA, nrow=jul, ncol = 1) #modeled flow being sent to Lowell
    #---------
    # initalize
    #----- 
    stor[1] <<- FC$AF[doy1[wy]] ####"this will only work if sent to global env "#initialize with actual storage on Jan 1
    Qin<- FC$Q[FC$WY == yrs[wy]]
    maxS <<- predMaxS(m) #vector of 198 days of max storage out to M days
    LowellAF<<-lowell$low_af[1:51]  ###"this needs an update to be associated with each water year - berak up the date"
    #----- run all the functions to get to discharge and updated storage
      for (day in 1:jul){ 
      volF<<- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow  
      availStor[day] <<- maxAF-stor[day]
      # Min flood control release for storage goals on April 1 and every 15 days after
      if (day < 91){
      minFCq[day]<<-minRelease(day, volF)
    ## # "Error in minFCq[day] <<- minRelease(day, volF) : 
## # replacement has length zero" #put unit test here
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
  
    out<- cbind(maxS[,1], storF, stor, minFCq, qo) #availStor dS, only use one arrow or it will overwrite in global env
    colnames(out)<-c('maxS','storF', 'stor', 'minQ', 'qo')
    results[[wy]]<-out
    #discharge[,wy]<-qo
    
  }
  #Q<-c(discharge)
  return(results)
}

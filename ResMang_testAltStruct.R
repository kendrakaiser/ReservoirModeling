predMaxS<- function(m){
  maxS<-matrix(data=NA, nrow = jul, ncol = m) 
  
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

minRelease<- function(day, volF, availStor){
  ix= which(prj$start <= day & prj$end >= day)
  vol1 = volF*prj$b[ix] + prj$c[ix]
  if (prj$start[ix] != day){
    vol2 = volF*prj$b[ix+1] + prj$c[ix+1]
    frac = (day - prj$start[ix])/(prj$end[ix]-prj$start[ix])
    volFmar <-  frac*vol1 + (1-frac)*vol2
  } else {volFmar = vol1}
  
  volFresid <- volF - volFmar
  FCvolAP<- reqStor(volFresid, 91) #flood control space required on april 1
  minEvac<- FCvolAP - availStor #minimum evaculation btw today and April 1
  minReleaseVol <- minEvac+volFmar
  Qmin <- (minReleaseVol*v2f)/(jul-day+1) #associated  qmin
  
  #if march 1st check lowell volume and if < full fill at whatever rate every day till full
  if (day<30){
    Lowell_cfs[day] <<-0
  } else if (day == 30 && LowellAF[day-1] < 155237){
    low_under <- 155237 - LowellAF[day-1]
    #calc how much under maximum storage the lake is
    Lowell_flow_est <<- low_under/(91-30) * v2f #days between start fill date to March 31st - equal discharge the whole time
  } else if (day >= 30 && LowellAF[day-1] < 155237){
    
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
    Lowell_cfs[day]<<-0}
  
  return(Qmin)
  ##If statements that constrain for high flows and ramp rates?
}
#determine minimum daily release after April 1
minReleaseApril<-function(day, volF, availStor){
  ix30 = findInterval(day, prjAP$doy)
  ix15=ix30-1
  volF_target15 <- volF*prjAP$b[ix15] + prjAP$c[ix15]
  volF_target30 <- volF*prjAP$b[ix30] + prjAP$c[ix30]
  residual15<- volF-volF_target15
  residual30<- volF- volF_target30
  
  
  if (day < jul-30){
    FCvol30<- reqStor(residual30, day+30) #required storage space
    minEvac30 <- FCvol30 - availStor
    minReleaseVol30 <- minEvac30 + volF_target30
    q30<-(minReleaseVol30*v2f)/(jul-day+1) 
  } else {q30 <- minQ}
  if (day < jul-15){
    FCvol15<-reqStor(residual15, day+15)
    minEvac15 <- FCvol15 - availStor
    minReleaseVol15 <- minEvac15 + volF_target15
    q15<-(minReleaseVol15*v2f)/(jul-day+1)
  } else {q15 <- minQ}
  
  Qmin <- max(q15, q30)
  
}

#forecast what the storage would be in s days given previous ∆ in S and make those changes over m days
forecastS<-function(s,m,day){
  if (day > s+1){
    dsdt= (stor[day] - stor[day-s])/s #calculate daily change in storage over past s days 
    
    storF[day] <- (dsdt*m)+stor[day] #forecasted flow for today is the change in storage* m days added to todays storage
  } else {storF[day]<- stor[day]}
  
  return(storF[day])
} 

evalS<- function(Qin, day, stor, maxS, Qmin,s,m, availStor){
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
  if (availStor <= (1000*v2f)){ 
    qo[day] <- qo[day] + 1000
    availStor<- availStor + 1000*f2v
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
  #qo[day]<<-qo[day]
  #dS[day]<<-dS[day]
  
  return(qo[day])
}


outflowStor<-function(s,m){
  results<-list()
  discharge<-matrix(data=NA, nrow = jul, ncol=21)
  
  for (wy in 1:21){
    #   set up blank matricies 
    #------------------------------------------
    stor<-matrix(data=NA, nrow = jul, ncol = 1)
    availStor<-matrix(data=NA, nrow = jul, ncol = 1)
    minFCq<-matrix(data=NA, nrow = jul, ncol = 1)
    storF<-matrix(data=NA, nrow = jul, ncol = 1)
    qo<-matrix(data=NA, nrow = jul, ncol = 1) #modeled outflow from reservoir
    dS<-matrix(data=NA, nrow = jul, ncol = 1)
    LowellAF<-matrix(data=NA, nrow=jul, ncol = 1) #modeled Lowell Acre-feet
    Lowell_cfs<-matrix(data=NA, nrow=jul, ncol = 1) #modeled flow being sent to Lowell
    #---------
    # initalize
    #----- 
    stor[1] <- FC$AF[doy1[wy]] ####"this will only work if sent to global env "#initialize with actual storage on Jan 1
    Qin<- FC$Q[FC$WY == yrs[wy]]
    maxS <- predMaxS(m) #vector of 198 days of max storage out to M days
    LowellAF[1:30]<-FC$LowellAF[doy1[wy]:(doy1[wy]+29)]  ###"this needs an update to be associated with each water year - break up the date"
    #----- run all the functions to get to discharge and updated storage
    for (day in 1:jul){ 
      volF<- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow  
      availStor[day] <- maxAF-stor[day]
      # Min flood control release for storage goals on April 1 and every 15 days after
      if (day < 91){
        minFCq[day]<-minRelease(day, volF, availStor[day])
        ## # "Error in minFCq[day] <<- minRelease(day, volF) : 
        # TODO: replacement has length zero" #put unit test here
      } else {
        minFCq[day]<-minReleaseApril(day, volF, availStor[day])
      }
      #minimum discharge is 240
      if (minFCq[day] < minQ){
        minFCq[day] <- minQ
      } 
      
      storF[day] = forecastS(s,m,day)
      #qo[day], dS[day] 
      qo[day]= evalS(Qin, day, stor, maxS, minFCq, s, m, availStor[day])
    }
    
    out<- cbind(maxS[,1], storF, stor, minFCq, qo, Lowell_cfs) #availStor dS, only use one arrow or it will overwrite in global env
    colnames(out)<-c('maxS','storF', 'stor', 'minQ', 'qo', 'Lowell_cfs')
    
    results[[wy]]<-out #there is something wrong with how the max Storage is being saved
    #discharge[,wy]<-qo
  }
  #Q<-c(discharge)
  return(results)
}

res1=outflowStor(14,14)

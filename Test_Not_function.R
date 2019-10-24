stor<-matrix(data=NA, nrow = jul, ncol = 1)
maxS<-matrix(data=NA, nrow = jul, ncol = m) 
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
  volF<<- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow  
  availStor[day] <- maxAF-stor[day]
  # Min flood control release for storage goals on April 1 and every 15 days after
  if (day < 91){
    minFCq[day]<-minRelease(day, volF)
    ## # "Error in minFCq[day] <<- minRelease(day, volF) : 
    # TODO: replacement has length zero" #put unit test here
  } else {
    minFCq[day]<-minReleaseApril(day, volF)
  }
  #minimum discharge is 240
  if (minFCq[day] < minQ){
    minFCq[day] <- minQ
  } 
  
  forecastS(s,m,day)
  evalS(Qin, day, stor, maxS, minFCq, s, m)
}

out<- cbind(maxS[,1], storF, stor, minFCq, qo, Lowell_cfs)

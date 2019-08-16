"Unit Testing Reservoir Code"
Error in minFCq[day] <- minRelease(day, volF) : 
  replacement has length zero
> volF
[1] 1230000
> day
[1] 1

stor<-matrix(data=NA, nrow = jul, ncol = 1)
maxS<-matrix(data=NA, nrow = jul, ncol = m) 
availStor<-matrix(data=NA, nrow = jul, ncol = 1)
minFCq<-matrix(data=NA, nrow = jul, ncol = 1)
storF<-matrix(data=NA, nrow = jul, ncol = 1)
qo<-matrix(data=NA, nrow = jul, ncol = 1) #modeled outflow from reservoir
dS<-matrix(data=NA, nrow = jul, ncol = 1)
LowellAF<-matrix(data=NA, nrow=jul, ncol = 1) #modeled Lowell Acre-feet
Lowell_cfs<-matrix(data=NA, nrow=jul, ncol = 1) #modeled flow being sent to Lowell

stor[1] <- FC$AF[doy1[wy]] ####"this will only work if sent to global env "#initialize with actual storage on Jan 1
Qin<- FC$Q[FC$WY == yrs[wy]]
maxS <- predMaxS(m) #vector of 198 days of max storage out to M days
LowellAF<-lowell$low_af[1:51] 
volF<- FC$volF[FC$WY == yrs[wy] & FC$doy == day] #todays forecasted inflow  
availStor[day] <- maxAF-stor[day]

ix= which(prj$start <= day & prj$end >= day) ## whats the difference with "=" v "<<-, <-"?
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
if (day >= 51 && day < 91 && LowellAF[day] < 155237){
  low_under <- 155237 - lowell$low_af[day] 
  #"This line 26 is where it be breaking -bc cant be the first day "
  #calc how much under maximum storage the lake is
  Lowell_cfs[day] <- low_under/39 * v2f #days between start fill date (Feb21st) to March 31st
  LowellAF <- lowell$low_af[day-1] + low_under/39
  
  if (Qmin < Lowell_cfs[day]){ #not set up right yet because this will only happen on feb21
    Qmin <- Lowell_cfs[day]
    LowellAF[day] <- lowell$low_af[day-1] + (Qmin*f2v)
  }
}
#"That all individually worked, something isnt translating to the functino"
minFCq[day]<-minRelease(day, volF)

"results<- outflowStor(s,m)
Error in if (day >= 51 && day < 91 && LowellAF[day - 1] < 155237) { : 
  missing value where TRUE/FALSE needed"

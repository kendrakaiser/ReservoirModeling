"Unit Testing Reservoir Code"
Error in minFCq[day] <- minRelease(day, volF) : 
  replacement has length zero
> volF
[1] 1230000
> day
[1] 1

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
if (day >= 51 && day < 91 && LowellAF[day-1] < 155237){
  low_under <- 155237 - lowell$low_af[day-1] 
  "This line 26 is where it be breaking -bc cant be the first day "
  #calc how much under maximum storage the lake is
  Lowell_cfs[day] <- low_under/39 * v2f #days between start fill date (Feb21st) to March 31st
  LowellAF <- lowell$low_af[day-1] + low_under/39
  
  if (Qmin < Lowell_cfs[day]){ #not set up right yet because this will only happen on feb21
    Qmin <- Lowell_cfs[day]
    LowellAF <- lowell$low_af[day-1] + (Qmin*f2v)
  }
}

minFCq[day]<-minRelease(day, volF)

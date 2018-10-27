# this function will calculate pediatric height z-score based on
# the CDC 2000 Height tables
# the function will take as inputs:
# height (in cm)
# age (in yrs)
# gender (1=male, or 0=female)
# it returns a z-score for the corresponding height
# limitations:
# will return invalid result if age < 0, and will use height z-score for 18-years when age > 18y

#setwd("~/Documents/Research/CKiDr_local/BPz")
load("~/Documents/Research/CKiDr_local/BPz/htz.RData")
# gender = 1 for male, 0 for female

htz <- function(ht,age,mf) {
  if (anyNA(c(ht,age,mf))) return(NA)
  # if arguments are outside of allowable ranges, return NA
  if (any(age<0,ht<0,!any(mf==0 | mf==1))) return(NA)
  
  # get age in months (and if >18y, use 216m)
  age<-age*12
  age<-min(216,age)
  # get the corresponding LMS, adjusted for position within interval
  agerow<-findInterval(age,htz_lms[,1])
  if (mf==1) cols<-2:4 else cols<-5:7
  LMS<-htz_lms[agerow,cols]+(age-htz_lms[agerow,1])*diff(as.matrix(htz_lms[agerow:(agerow+1),cols]))
  htz<-((ht/LMS[2])^LMS[1]-1)/(LMS[1]*LMS[3])
  return(as.numeric(htz))
}

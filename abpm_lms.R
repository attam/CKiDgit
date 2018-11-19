# function that uses LMS method to generate percentile or z-score based on provided BP, or inversely,
# provides BP values for a given percentile, if percentile=TRUE
# input must be either in the form age=age (in years) or height=height(in cm)
# bp = 1 (systolic) or 2 (diastolic)
# period=1 (24h) or =2 (day) or =3 (night)
# gender = 1 (male) or 0 (female)
lms<-function(age=NULL,height=NULL,bp=NULL, gender,period,percentile=NULL) {
  load("abpm_norm.RData")
  if (gender==0) {
    abpm.norm.age<-abpm.norm.age[,30:58]
    abpm.norm.height<-abpm.norm.height[1:12,30:58]
  } else {
    abpm.norm.age<-abpm.norm.age[,1:29]
    abpm.norm.height<-abpm.norm.height[,1:29]
  }
  cols<-c(1,2+1:3+(period-1)*3+(9*(bp-1)))
  if (!is.null(age)) {
    lms_data<-abpm.norm.age[,cols]
    t<-age
  } else {
    lms_data<-abpm.norm.height[,cols]
    t<-height
  }
  row.num<-which.min(abs(t-lms_data[,1]))
  lms<-lms_data[row.num,2:4]
  if (exists("sbp") | exists("dbp")) {
    bp<-ifelse(exists("sbp"),sbp,dbp)
    z<-(((bp/lms[2])^lms[1])-1)/(lms[1]*lms[3])
    return(pnorm(z))
  } else {
    bp<-lms[2]*(1+lms[1]*lms[3]*qnorm(percentile))^(1/lms[1])
    return (bp)
  }
}
# function to calculate blood pressure percentiles based on the old (4th report) guidelines
# reference: http://pediatrics.aappublications.org/content/114/Supplement_2/555
# inputs:
#    - age (yrs)
#    - height percentile (0-100)
#    - gender (male=1, female=0)
#    - sysdia (1 for systolic, 2 for diastolic)
#    - bp (either sbp or dbp) - if inverse function is used, the percentile argument must be set (see below)
#    - percentile (eg, 85) - inverse function to calculate BP at corresponding desired percentile
# limitations: age must be within range 1-18yrs, otherwise result will be NA
load("bp4th_coefs.RData")
# if (!exists("bp4th_coefs")) load("bp4th_coefs.RData")
# if (!all(c(11,4)==dim(bp4th_coefs))) load("bp4th_coefs.RData")
# convert height to inches and take only coefficients needed for computation
#height<-height/2.54

bpp4<- function(age,heightp,gender,sysdia,bp=NULL,percentile=NULL) {
  # check that all parameters are provided and are valid
  if (any(!between(age,1,18),anyNA(c(bp,age,gender,heightp)))) return(NA)
  # extract only the relevant coefficients based on gender and sysdia
  bp4th_coefs.1<-bp4th_coefs[c(2,4)-gender]
  bp4th_coefs.1<-bp4th_coefs.1[c(sysdia)]
  # calculation of average BP (mu) is based on reference above (appendix B)
  first<-function(x,age) bp4th_coefs.1[1+x,1]*(age-10)^x
  second<-function(x,heightz) bp4th_coefs.1[x+5,1]*heightz^x
  sum1<-sum(unlist(sapply(1:4,function(x) first(x,age))))
  sum2<-sum(unlist(sapply(1:4,function(x) second(x,qnorm(heightp/100)))))
  mu<-unlist(bp4th_coefs.1[1,1]+sum1+sum2)
  if (!is.null(percentile)) {
    bpp4<-round(qnorm(percentile/100)*unlist(bp4th_coefs.1[10,1])+mu,2)}
  else {
  bpp4<-round(pnorm((bp-mu)/unlist(bp4th_coefs.1[10,1]))*100,2)}
  return(bpp4)
  
  }

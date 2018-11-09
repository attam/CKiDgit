# this function will calculate pediatric blood pressure percentiles based on
# the AAP 2017 BP tables
# the actual source of the data is from Rosen et al, (2008)
# other resources used:
# link to sas version of this function
# https://sites.google.com/a/channing.harvard.edu/bernardrosner/pediatric-blood-press/childhood-blood-pressure/childhoodbppctsas
# and BC children's anthopometric excel file:
# http://www.bcchildrens.ca/health-professionals/clinical-resources/endocrinology-diabetes/tools-calculators#Anthro--calculators
#
# the function will take as inputs:
# either systolic or diastolic BP (must specify using 1 or 2 as the last argument, respectively)
# age in years (must be between 0 and 18 years)
# height in cm (must be > 0, and height z-score must be within the range +/- 3.09, otherwise, will return NA)
# 
# it returns a percentile for both SBP and DBP
# sbp and dbp must be > 0
# will return 0 or 100 if bp percentile is <1st or >99th, respectively
# outputstyle: percentile = 1, z-score=2, lookup=3 (*) note: z-scores for percentiles
# outside of the range 1-99 will be assigned z-scores of -/+ 2.33
# * if outputstyle = 3, the output will use first parameter as desired BP percentile and will return
# a BP corresponding with the desired percentile

library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
source('BPz/htz.R')
#setwd("")
if (!exists("BPcoef") | !exists("all_t")) load("BPz/BPcoefs.RData")
# gender = 1 for male, 0 for female
# sysdia = 1 for systolic, 2 for diastolic

bpp <- function(bp=bp,age,ht,mf,sysdia=sysdia,outputstyle=outputstyle) {
  if (anyNA(c(bp,age,ht,mf))) return(NA)
  # get height z-score, and if outside the range +/- 3.09, set it to 3.09
  htz<-htz(ht,age,mf)
  # if arguments are outside of allowable ranges, return NA
  if (any(!between(age,0,18),ht<0,abs(htz)>3.09,bp<0, is.na(bp))) return(NA)

  # set t values based on gender
  if(mf==1) {
    t<-all_t[,1:3]
    bpc<-BPcoef[,1:26]
    gen_o<-0
  }
    else {
      t<-all_t[,4:6]
      bpc<-BPcoef[,27:52]
      gen_o<-3
    }
  w<-(age-10.)*(ht-150.+gen_o)
  x2<-rbind(a=max(ht-t[1,1],0),b=max(ht-t[4,1],0),c=max(ht-t[5,1],0))
  y2<-rbind(a=max(age-t[1,2],0),b=max(age-t[4,2],0),c=max(age-t[5,2],0))
  w2<-rbind(a=max(w-t[1,3],0),b=max(w-t[4,3],0),c=max(w-t[5,3],0))
  x3a<-max(ht-t[2,1],0)
  x4a<-max(ht-t[3,1],0)
  y3a<-max(age-t[2,2],0)
  y4a<-max(age-t[3,2],0)
  w3a<-max(w-t[2,3],0)
  w4a<-max(w-t[3,3],0)
  
  x2s<-(x2[1]^3-x2[2]^3*((t[5,1]-t[1,1])/(t[5,1]-t[4,1]))+x2[3]^3*((t[4,1]-t[1,1])/(t[5,1]-t[4,1])))/100
  x3s<-(x3a^3-x2[2]^3*((t[5,1]-t[2,1])/(t[5,1]-t[4,1]))+x2[3]^3*((t[4,1]-t[2,1])/(t[5,1]-t[4,1])))/100
  x4s<-(x4a^3-x2[2]^3*((t[5,1]-t[3,1])/(t[5,1]-t[4,1]))+x2[3]^3*((t[4,1]-t[3,1])/(t[5,1]-t[4,1])))/100
  
  y2s<-(y2[1]^3-y2[2]^3*((t[5,2]-t[1,2])/(t[5,2]-t[4,2]))+y2[3]^3*((t[4,2]-t[1,2])/(t[5,2]-t[4,2])))/100
  y3s<-(y3a^3-y2[2]^3*((t[5,2]-t[2,2])/(t[5,2]-t[4,2]))+y2[3]^3*((t[4,2]-t[2,2])/(t[5,2]-t[4,2])))/100
  y4s<-(y4a^3-y2[2]^3*((t[5,2]-t[3,2])/(t[5,2]-t[4,2]))+y2[3]^3*((t[4,2]-t[3,2])/(t[5,2]-t[4,2])))/100  
  
  w2s<-(w2[1]^3-w2[2]^3*((t[5,3]-t[1,3])/(t[5,3]-t[4,3]))+w2[3]^3*((t[4,3]-t[1,3])/(t[5,3]-t[4,3])))/100^2
  w3s<-(w3a^3-w2[2]^3*((t[5,3]-t[2,3])/(t[5,3]-t[4,3]))+w2[3]^3*((t[4,3]-t[2,3])/(t[5,3]-t[4,3])))/100^2
  w4s<-(w4a^3-w2[2]^3*((t[5,3]-t[3,3])/(t[5,3]-t[4,3]))+w2[3]^3*((t[4,3]-t[3,3])/(t[5,3]-t[4,3])))/100^2
  
  if (sysdia==1) dbp_o<-0 else dbp_o<-13
  fx<-bpc[,1+dbp_o]+bpc[,6+dbp_o]*ht+bpc[,7+dbp_o]*x2s+bpc[,8+dbp_o]*x3s+
    bpc[,9+dbp_o]*x4s+bpc[,2+dbp_o]*age+bpc[,3+dbp_o]*y2s+bpc[,4+dbp_o]*y3s+bpc[,5+dbp_o]*y4s+
    bpc[,10+dbp_o]*w+bpc[,11+dbp_o]*w2s+bpc[,12+dbp_o]*w3s+bpc[,13+dbp_o]*w4s
  dif<-abs(bp-fx)
  p<-which.min(as.numeric(unlist(dif)))
    if (bp<fx[1]) p<-0.99 
    if (bp>fx[99]) p<-99.1 
  bpp<-p
  if (outputstyle==2) {
  bpp<-qnorm(bpp/100)
  }
  if (outputstyle==3){
    bpp<-fx[bp]
  }
  bpp
  }


# function to categorize patient's clinic BP based on the 4th report (1996,revised in 2004) guidelines:
# input: clinic SBP and DBP, age (yrs), gender (male=1, female=0), height (cm)
# for reference, see Table 3 in PDF at the following URL: http://pediatrics.aappublications.org/content/114/Supplement_2/555#ref-1

# output is BP status with the following codes:
# 0=normotensive
# 1=pre-hypertension
# 2=Stage 1 HTN
# 3=Stage 2 HTN

# bpp4<- function(age,heightp,gender,sysdia,bp=NULL,percentile=NULL)

bpstatus4th<-function(sbp,dbp,age,gender,heightp) {
  sbpp4<-bpp4(age,heightp,gender,1,sbp)
  dbpp4<-bpp4(age,heightp,gender,2,dbp)
  if (all(sbpp4<90,dbpp4<90)) return (0)
  sbpp<-sapply(c(90,95,99),function(x) bpp4(age,heightp,gender,1,percentile=x))
  dbpp<-sapply(c(90,95,99),function(x) bpp4(age,heightp,gender,2,percentile=x))
  if (any(sbp>sbpp[3]+5,dbp>dbpp[3]+5,sbp>=160,dbp>=100)) return(3)
  if (any(between(sbp,sbpp[2],sbpp[3]+5),between(dbp,dbpp[2],dbpp[3]+5),between(sbp,140,159.9),between(dbp,90,99.9))) return(2)
  return(1)
  
}
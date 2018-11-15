# function to categorize patient's clinic BP based on the 4th report (1996,revised in 2004) guidelines:
# input: clinic SBP and DBP, age (yrs), gender (male=1, female=0), height (cm)
# for reference, see Table 3 in PDF at the following URL: http://pediatrics.aappublications.org/content/114/Supplement_2/555#ref-1

# output is BP status with the following codes:
# 0=normotensive
# 1=pre-hypertension
# 2=Stage 1 HTN
# 3=Stage 2 HTN
# note: in cases where values are missing (eg, SBP, DBP, SBP % or DBP %), the BP status is determined based only on the
# available data unless all data is missing, in which case the result is NA
# bpp4<- function(age,heightp,gender,sysdia,bp=NULL,percentile=NULL)
# note: some discrepancies exist between calculated and reported BP percentile (most likely due to age uncertainty)
# for consistency with data, the optional arguments sbpp and dbpp will allow the user to directly input these values
# in which case the percentiles will not be calculated

bpstatus4th<-function(sbp,dbp,age,gender,heightp, sbpp4=NULL, dbpp4=NULL) {
  if (!is.na(age) & age>18) {
    if (!all(is.na(c(sbp,dbp)))) {
    if (sbp<120 & dbp<80) return (0)
    if (sbp>=160 | dbp>=100) return (3)
    if ((sbp>=140 & sbp<160) | (dbp>=90 & dbp<100)) return(2)
    return (1)
    } else return(NA)
  }
  if (is.null(sbpp4) | is.na(sbpp4)) sbpp4<-bpp4(age,heightp,gender,1,sbp)
  if (is.null(dbpp4) | is.na(dbpp4)) dbpp4<-bpp4(age,heightp,gender,2,dbp)
  if (!all(is.na(c(sbpp4,dbpp4)))) {
  if (all(sbpp4<90,dbpp4<90,sbp<120,dbp<80, na.rm=T)) return (0)
  sbpp<-bpp4(age,heightp,gender,1,percentile=c(90,95,99))
  dbpp<-bpp4(age,heightp,gender,2,percentile=c(90,95,99))
  if (any(sbp>sbpp[3]+5,dbp>dbpp[3]+5,sbp>=160,dbp>=100, na.rm=T)) return(3)
  if (any(between(sbp,sbpp[2],sbpp[3]+5),between(dbp,dbpp[2],dbpp[3]+5),between(sbp,140,159.9),between(dbp,90,99.9),na.rm=T)) return(2)
  return(1)
  } else return(NA)
  
}
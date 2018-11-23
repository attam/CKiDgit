# function to categorize patient's clinic BP based on new 2017 AAP guidelines:
# input: clinic SBP and DBP, age (yrs), gender (male=1, female=0), height (cm)
# for reference, see Table 3 in PDF at the following URL: https://solutions.aap.org/DocumentLibrary/pcowebinars/2017%20Hypertension%20Webinar.pdf

# output is BP status with the following codes:
# 0=normotensive
# 1=elevated BP
# 2=Stage 1 HTN
# 3=Stage 2 HTN

# note: this function is permissive for missing values and classification is made based on available data
# eg 1, if SBP is missing, and DBP<90th percentile and <80mmHg, then the patient is classified as normotensive
# eg 2, if DBP is missing and SBP >= 160mmHg, then patient is classified as stage 2 HTN
# if insufficient data is available, the result is NA

BPstatus2017<- function(sbp,dbp,age,gender,height) {
  # check if at least one of SBP or DBP exists and if age range is valid (note: <1y is invalid, will return NA)
  if (all(is.na(c(sbp,dbp)))) return(NA)
  if (any(age<1,is.na(age))) return(NA)
  if (all(age>=1,age<13)) {
    #reference bpp <- function(age,ht,mf,sysdia=sysdia,bp=NULL, z=NULL, percentile=NULL)
    bpp<-c(bpp(age,height,gender,1,sbp,z=F),bpp(age,height,gender,2,dbp,z=F))
    bp95<-c(bpp(age,height,gender,1,percentile=95),bpp(age,height,gender,2,percentile=95))
#    print(c(sbp,dbp,age,height,gender))
    if (all(bpp<90,sbp<120,dbp<80,na.rm = T)) {
      return(0)}
    if (any(c(sbp>=140,dbp>=90,sbp>=12+bp95[1],dbp>=12+bp95[2]), na.rm=T)) return(3)
    if (any(c(sbp>=130 & sbp<140, dbp>= 80 & dbp<90, sbp>=bp95[1] & sbp<bp95[1]+12, dbp>=bp95[2] & dbp<bp95[2]+12), na.rm=T)) return(2)
    if (any(c(bpp[1]>=90 & bpp[1]<95,bpp[2]>=90 & bpp[2]<95,sbp>=120 & sbp<bp95[1],dbp>=80 & dbp<bp95[2]))) return (1)
    return (-1)
  } 
  if (age>=13) {
    if (all(sbp<120,dbp<80, na.rm=T)) return(0)
    if (any(sbp>=140,dbp>=90, na.rm=T)) return(3)
    if (any(sbp>=130,dbp>=80, na.rm=T)) return(2)
    if (any(sbp>=120 & sbp<130, dbp<80)) return(1)
    else return(-1)
  }

}


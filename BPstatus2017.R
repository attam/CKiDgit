# function to categorize patient's clinic BP based on new 2017 AAP guidelines:
# input: clinic SBP and DBP, age (yrs), gender (male=1, female=0), height (cm)
# for reference, see Table 3 in PDF at the following URL: https://solutions.aap.org/DocumentLibrary/pcowebinars/2017%20Hypertension%20Webinar.pdf

# output is BP status with the following codes:
# 0=normotensive
# 1=elevated BP
# 2=Stage 1 HTN
# 3=Stage 2 HTN

BPstatus2017<- function(sbp,dbp,age,gender,height) {
  # check that all parameters are provided and are valid
  if (any(abs(htz(height,age,gender))>3.09,anyNA(c(sbp,dbp,age,gender,height)))) return(NA)
  # check if age range is valid (note: <1y is invalid, will return NA)
  if (age<1) return(NA)
  if (all(age>=1,age<13)) {
    bpp<-c(bpp(sbp,age,height,gender,1,1),bpp(dbp,age,height,gender,2,1))
#    print(c(sbp,dbp,age,height,gender))
    if (all(bpp<90,na.rm = T)) {
      return(0)}
    if (any(c(sbp>=140,dbp>=90,sbp>=12+bpp(95,age,height,gender,1,3),dbp>=12+bpp(95,age,height,gender,2,3)))) return(3)
    if (any(sbp>min(bpp(95,age,height,gender,1,3),120),dbp>min(bpp(95,age,height,gender,2,3),80))) return(2)
    else return(1)
  } 
  if (age>=13) {
    if (all(sbp<120,dbp<80)) return(0)
    if (any(sbp>=140,dbp>=90)) return(3)
    if (any(sbp>=130,dbp>=80)) return(2)
    else return(1)
  }

}

# function to find closest LVMI percentile based on data provided;
# inputs:
# - gender (male=1,female=0)
# - age (years)
# - LVMI (g/m2)
# reference: Khoury PR, Mitsnefes M, Daniels SR, Kimball TR. Age-Specific Reference Intervals for Indexed Left Ventricular Mass in Children. Journal of the American Society of Echocardiography. 2009;22(6):709-714. doi:10.1016/j.echo.2009.03.003)

LVMI_p<-function(gender,age,lvmi){
  load("LVMI_table.RData")
  if(is.na(age)) return(NA)
  if(age>=18) age<-17.9
  percentile<-c(10,25,50,75,90,95)
  slicenum<-ifelse(gender==0,2,1)
  rownum<-age %/% 2 + 2
  dt = LVMI_table[rownum,,slicenum]
  closestVal<- percentile[which.min(abs(dt-lvmi))]
  return(closestVal)
}
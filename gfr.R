# function to calculate estimated GFR
# based on ckidfull equation [the best estimated gfr for research based on Schwartz and Schneider 2012]
# function below calculates gfr based on full ckid equation if all data available, otherwise returns NA
# if optional switch permissive is set, will default to ckidfull, but if only height and Scr are available
# will return the bedside schwartz gfr
GFR<-function(height,cr,cystatin,bun,gender, permissive=NULL) {
  if (anyNA(c(height,cr,cystatin,bun,gender))) {
    if (!is.null(permissive)){
      if (anyNA(c(height,cr))) return(NA)
      return (((height/100)*0.413)/cr)
    } else return (NA)}
  term1=((height/100)/cr)^0.456
  term2=(1.8/cystatin)^0.418
  term3=(30/bun)^0.079
  term4=((height/100)/1.4)^0.179
  gfr<-39.8*term1*term2*term3*term4*ifelse(gender==1,1.076,1)
  return(gfr)
}
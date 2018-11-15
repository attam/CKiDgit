# function that determines a drug dose index (DDI)
# where DDI is a measure of how high a dose is based on the maximum advised dose for a given drug
# the formula is;
# DDI = dose/max dose
# max dose is the minimum of (maximum dose in mg/kg, max dose in mg)
# inputs:
# - dose (mg); drug (using the standard names in the data); weight (in kg); gfr (in mL/min/1.73m2)
# output will be DDI (unitless, percent of maximum advised dose)
# note: if gfr plays a role in maximum dose advised, it is adjusted accordingly based on references below if specific guidelines are provided on dose adjustments
# if no specific guidelines are provided (eg, avoid if GFR<30), DDI is not adjusted
# 1. Chu PY, Campbell MJ, Miller SG, Hill KD. Anti-hypertensive drugs in children and adolescents. World J Cardiol. 2014;6(5):234-244. doi:10.4330/wjc.v6.i5.234
# 2. Lexicomp Online. http://online.lexi.com.ezproxy.cmh.edu/lco/action/home. Accessed November 10, 2018.
# 3. Misurac J, Nichols KR, Wilson AC. Pharmacologic Management of Pediatric Hypertension. Pediatr Drugs. 2016;18(1):31-43. doi:10.1007/s40272-015-0151-3

# read the reference data on drug max doses
dose_ref <- read_csv("dosing_ref_antiHTN.csv")

DDI <-function(dose, drug,weight,gfr) {
  drug<-toupper(drug)
  if (dose<=0|is.na(match(drug[1],dose_ref$Drug))) return(NA)
  dose_ref<-dose_ref %>% filter(Drug==drug)
  
  # adjustment for decreased GFR
  if (!is.na(dose_ref[4])) {
    dose_adj<-1
  if (drug=="LISINOPRIL") dose_adj<-ifelse(between(gfr,10,50), 0.5,ifelse(gfr<10,0.25,1))
  if (drug=="CAPTOPRIL") dose_adj<-ifelse(between(gfr,10,50), 0.75,ifelse(gfr<10,0.5,1))
  dose_ref[2:3]<-dose_ref[2:3]*dose_adj
  }
  if (is.na(dose_ref[2])) {DDI<-dose/dose_ref[3]} else {
  DDI<-ifelse(weight*dose_ref[2]>dose_ref[3], dose/dose_ref[3],(dose/weight)/dose_ref[2])}
  names(DDI)<-"DDI"
  return (unlist(DDI))
}
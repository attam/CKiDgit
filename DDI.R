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
# modification: 
# - due to fact that drug dose-response curve flattens out after reaching at some point,
#   maximum DDI will be 1 (ie, any calculation of DDI>1 will be given value of 1)
dose_ref <- read_csv("dosing_ref_antiHTN.csv")

DDI <-function(dose, drug,age,weight,gfr,get_max_dose=NULL,renal=NULL,type=NULL) {
  drug<-toupper(drug)
  if (is.na(age)|dose<=0|is.na(match(drug[1],dose_ref$Drug))) return(NA)
  # special cases:
  special<-c("OLMESARTAN","IRBESARTAN","MINOXIDIL")
  if (drug %in% special){
    if (drug=="OLMESARTAN") {
      weight_cat<-cut(weight,c(0,20,35,Inf),right=FALSE)
      ddi<-switch(as.numeric(weight_cat),(dose/weight)/0.6,dose/20,dose/40)
    }
    if (drug=="IRBESARTAN") {
      age_cat<-cut(age,c(0,6,13,Inf),right=FALSE)
      ddi<-switch(as.numeric(age_cat),NA,dose/150,dose/300)
    }
    if (drug=="MINOXIDIL") {
      age_cat<-cut(age,c(0,12,Inf),right=FALSE)
      ddi<-switch(as.numeric(age_cat),dose/50,dose/100)
    }
  } else {
  dose_ref<-dose_ref %>% filter(Drug==drug)
  dose_adj<-1
  # adjustment for decreased GFR
  if (!is.na(dose_ref[4])) {
  if (drug=="LISINOPRIL") dose_adj<-ifelse(between(gfr,10,50), 0.5,ifelse(gfr<10,0.25,1))
  if (drug=="ENALAPRIL"|drug=="CAPTOPRIL") dose_adj<-ifelse(between(gfr,10,50), 0.75,ifelse(gfr<10,0.5,1))
  dose_ref[2:3]<-dose_ref[2:3]*dose_adj
  }
  if (is.na(dose_ref[2])) {
    ddi<-dose/dose_ref[3]
    max.dose<-dose_ref[3]
    } else {
    if(!is.null(type)){
      temp.min<-unlist(c(dose/dose_ref[3],(dose/weight)/dose_ref[2]))
      which_temp.min<-which.min(temp.min)
      ddi<-ifelse(weight*dose_ref[2]>dose_ref[3],min(temp.min),dose/dose_ref[3])
      max.dose<-ifelse(weight*dose_ref[2]>dose_ref[3],ifelse(which_temp.min==1,dose_ref[3],dose_ref[2]),dose_ref[3])
      } else {
        ddi<-ifelse(weight*dose_ref[2]>dose_ref[3], dose/dose_ref[3],(dose/weight)/dose_ref[2])
        max.dose<-ifelse(weight*dose_ref[2]>dose_ref[3],dose_ref[3],dose_ref[2])
      }
    #names(ddi)<-"DDI"
    }
  }
  ddi<-as.numeric(unlist(ddi))
  out<-NA
  if (!is.null(get_max_dose)) out<-max.dose
  if (!is.null(renal)) out<-dose_adj
  if (is.null(get_max_dose)&is.null(renal)) out<-unlist(min(c(ddi,1)))
  return (out)
          }
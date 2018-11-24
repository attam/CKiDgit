#loading libraries
library(readr)
library(corrplot)
library(dplyr)
library(reshape2)
library(histogram)
library(gplots)
library(tidyr)
library(stats)
#this change is to test git functionality...
#library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
#library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
#library("histogram", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")

wd<-if_else(Sys.info()["nodename"]=='matta', "~/Documents/Research/CKiDgit", "C:/Users/Orit/Downloads/CKiD/CKiDgit")
setwd(wd)
#setwd("C:/Users/Orit/Downloads/CKiD/CKiDgit")

# importing data
medsum_full <- read.csv("data/medsum_full.csv")
cardio<-read.csv("data/cardio.csv")
echo<-read.csv("data/echo.csv")
f13 <- read.csv("data/f13.csv")
f14 <- read.csv("data/f14.csv")
f15 <- read.csv("data/f15.csv")
socdem <- read.csv("data/socdem.csv")
growth <- read.csv("data/growth.csv")
pe <- read.csv("data/pe.csv")
gfrcalibratedsummary <- read.csv("data/gfrcalibratedsummary.csv")
advr <- read.csv("data/advr.csv")
kidhist <- read.csv("data/kidhist.csv")
l05 <- read.csv("data/l05.csv")
l02 <- read.csv("data/l02.csv")
l51 <- read.csv("data/l51.csv")
l31 <- read.csv("data/l31.csv")
gfr<-gfrcalibratedsummary
gh <- read.csv("data/gh.csv")

# data integrity issues:
# review of the raw data showed several instances where the historical data is not consistent
# eg, changes in LBW status (growth.csv) in the same participant across different visits (this should be the same)
# this is also seen for premature status (growth.csv)
# the table generated below shows the flagged BP counts grouped by visit (ie, problem measuring BP)
# consider excluding those with status '1'
# pe %>% select(VISIT, PEPRBLPM)%>% group_by(VISIT) %>% table %>% addmargins()
# growth appears to contain the best weights after comparing data from pe and growth

echo_df<-echo %>% select(CASEID,VISIT,ECHODATEY,LVMI,LVHF,LVHE,CAUTION)
# echo_df %>% group_by(VISIT) %>% filter(CAUTION==0) %>% summarise_at(vars(ECHODATEY), funs(date_mean=mean,date_min=min,date_max=max, n_pts=length))


# combining data files
# selection of candidate variables to include from each data file
# combining using inner_join


test<-f13 %>% select(CASEID,VISIT,GHVISDAT,GHHGRADM,GHHGRADF,GHTOTINC,GHBFAMHB,contains("HBP")) %>% tbl_df
test<-full_join(socdem %>% select(-VERSION),test)
test<-full_join(f14 %>% select(CASEID,VISIT,MHCHRLD,MH_BPD, MHFEMALE),test)
test<-full_join(f15 %>% select(CASEID,VISIT,contains("NSSTER"),NSSTPYBP,NSSTSEBP), test)
test<-full_join(growth %>% select(CASEID,VISIT,LBW,PREMATURE,AVHEIGHT,HTPCTAG,AVWEIGHT,WTPCTAG,BMI,BMIPCTAG),test)
test<-full_join(pe %>% select(CASEID,VISIT,PEPRBLPM,PEPRWEIT),test)
test<-full_join(gfr %>% select(CASEID,VISIT,MALE1FE0,SCR,BUN,CYC_DB,BEDGFR, BSA),test)
test<-full_join(advr %>% select(CASEID,VISIT,AENSHEAD,AENSDIZZ,AENSLWBP),test)
test<-full_join(gh %>% select(CASEID,VISIT,GHGENDER,GHWKSBDD,GHPREMIE),test)
test<-full_join(kidhist %>% select(CASEID,DOB,BSDATE, CKDONST,PRIMDX,GNGDIAG),test)
test<-full_join(l05 %>% select(CASEID,VISIT,RLSERCRE,RLURPROT,RLURMALB,RLURCREA),test)
test<-full_join(echo_df %>% select(CASEID,VISIT,LVHF,LVHE,LVMI,ECHODATEY,CAUTION),test)
test<-full_join(cardio %>% select(CASEID,VISIT,DB_DATE,SBP,SBPINDXAGH,SBPPCTAGH,SBPZAGH,DBP,DBPINDXAGH,DBPPCTAGH,DBPZAGH,ABPM_DATE,ABPMSUCCESS,SHYPAGH,DHYPAGH),test)

# this was compared with age from pierce et al and correlates well, but not perfectly
# note: age from pierce is in integers, not decimals
# due to DOB provided in integer form, the true age may be off by up to 1 year (eg, born 1/1/1994 = DOB 1994, baseline date 12/31/2000 = BSDATE 2000, calculated age = 6yrs, true age=7yrs)
test$age<-(test$BSDATE-test$DOB)+test$DB_DATE

# gender: since some values of MALE1FE0 are missing, will solve missing values by taking the mean of MALE1FE0 for each case (no change in gender over time)
test<-test %>% group_by(CASEID) %>% mutate(MALE1FE0=mean(MALE1FE0, na.rm=TRUE))

#add ckidfull [the best estimated gfr for research based on Schwartz and Schneider 2012]
# function below calculates gfr based on full ckid equation if all data available, otherwise returns NA
# if optional switch permissive is set, will provide alternative bedside schwartz gfr in case where full ckidgfr is not available
gfr<-function(height,cr,cystatin,bun,gender, permissive=NULL) {
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
# gfr will be based on ckidfull equation if the data is available, otherwise bedside schartz equation is used
test$gfr<-mapply(gfr,test$AVHEIGHT,test$SCR,test$CYC_DB,test$BUN,test$MALE1FE0, permissive=TRUE)
test$gfr[which(is.na(test$gfr))]<-test$BEDGFR[which(is.na(test$gfr))]

# add column for urine protein:creatinine ratio based on RLURPROT and RLURCREA
# add column for percentage of total urine protein that is albumin (Ualb_pct)
test<-test %>% mutate(Upc=RLURPROT/RLURCREA)
test<-test %>% mutate(Ualb_pct=RLURMALB/RLURPROT)

# add a column for CKD stage using ckidfull
test$CKD_stage <- cut(test$gfr, 
                      breaks = c(Inf, 90, 60, 30, 15, -Inf), 
                      labels=c(1:5),ordered_result=TRUE, 
                      right = T)

# create categories for proteinuria based on cutoffs noted
test$Upc.factor<-cut(test$Upc,breaks=c(-Inf,0.5,1,2,Inf),labels=c("normal","mild","moderate","severe"),ordered_result = T,right=F)

# loads the correct BP medication names into BP_medlist
# and corrects BP medication names for brand names and misspellings
# corrected names are assigned a new variable called med.corrected
# display information about timing of echocardiogram data, grouped by visits
# note: BP_medlist will display the misspelled/brand names of BP meds, frequency at visit 10?, and corrected names
# categorize the BP meds into groups based on BP_medgroups.csv
load("BP_medlist.RData")
BP_medgroups <- read.csv("BP_medgroups.csv")
medsum_full$med.corrected<-BP_medlist$Corrected.Name[match(medsum_full$MSMEDICA,BP_medlist$Var1)]
medsum_full$BPmedgroup<-BP_medgroups$bp_group[match(medsum_full$med.corrected,BP_medgroups$x)]
# remove all unintelligible medications from the data (may be viewed using BP_medlist[which(BP_medlist$Corrected.Name=="*invalid*"),c(2:3)])
medsum_full<-medsum_full %>% filter(med.corrected!="*invalid*")
# combine split dosing medications (eg, some patients taking different doses of same medication at different times of day, eg, labetalol)
medsum_full<-medsum_full %>% group_by(CASEID,VISIT,MSVISDAT, BPmedgroup,med.corrected) %>% summarise_at(vars(DLYDOSE,DLYFREQ,MSMISS7D), sum) %>% ungroup

# determine compliance as percent of doses taken, based on reported missed doses per week and daily frequency of medication
medcompliance<-function(miss7d,dlyfreq) {
  if (miss7d<0) miss7d<-0
  if (dlyfreq<0) dlyfreq<-NA
  return((dlyfreq*7-miss7d)/(dlyfreq*7))
}
medsum_full$compliance<-mapply(medcompliance,medsum_full$MSMISS7D,medsum_full$DLYFREQ)
medsum_full<-medsum_full %>% group_by(CASEID,VISIT) %>%  mutate(mean_compliance=mean(compliance,na.rm=TRUE)) %>% ungroup

# more organizing of medication data
# these lines generatet columns for std_dose (dose per kg, with AVWEIGHT from growth data), corrects
# for split dosing (eg, labetalol), and nests drug data into a 5 column list called rx_info
# the drug dose index is the 5th column, and uses gfr based on ckidfull equation if possible, otherwise bedside GFR is used
source("DDI.R")
medsum_full<-full_join(test %>% select(CASEID,VISIT,AVWEIGHT,gfr),medsum_full)
# function to calculate standard dose (in mg/kg)
std_dose<-function(dlydose,weight){
  if (anyNA(c(dlydose,weight))) return (NA)
  if(dlydose>0 & weight>0) {
    return(dlydose/weight)
  } else return(NA)
}
medsum_full$std_dose<-mapply(std_dose,medsum_full$DLYDOSE,medsum_full$AVWEIGHT)
medsum_full$DDI<-mapply(DDI,medsum_full$DLYDOSE, medsum_full$med.corrected,medsum_full$AVWEIGHT, medsum_full$gfr)
# identify duplicates (cases with visit codes that have coresponding visit dates that are not the same) and eliminate all the corresponding medication data
# justifications 1) unable to determine whether duplicate visits represent 'addition' of medications or 'switching' between meds
# 2) some cases have same visit code and corresponding visit dates separated in time by up to 4.28 years (obviously an error in recording data)
medsum_full<-medsum_full %>% group_by(CASEID,VISIT) %>% mutate(ndup=n_distinct(MSVISDAT), visitmean=mean(unique(MSVISDAT,na.rm=TRUE)), visitsd=sd(unique(MSVISDAT,na.rm=TRUE)))
medsum_full %>% group_by(VISIT) %>% filter(ndup>1) %>% summarise(npts.dup.visit=n_distinct(CASEID)) # show how many patients at each visit are duplicated
medsum_full %>% group_by(VISIT) %>% arrange(desc(visitsd)) # show the extreme values of visit date sd
medsum_full<-medsum_full %>%filter(ndup==1) # remove all duplicates from the data



medsum_full<-medsum_full %>% group_by(CASEID,VISIT) %>% mutate(n_agents=n_distinct(med.corrected, na.rm=TRUE))
medsum_full<-medsum_full %>% ungroup %>% nest(c(6,8,9,11,13,14), .key="rx_info") %>% select(-c(3,4,7))
medsum_full<-dcast(medsum_full,CASEID+VISIT+MSVISDAT+n_agents+mean_compliance~med.corrected, value.var="rx_info")

# save(medsum_full, file='medsum_full.1.RData')
load(file='medsum_full.1.RData')

test<-full_join(medsum_full,test)

# replace all NULL values of drugs with a tibble containing NA's using change_null_to_list function above
# note: the line below takes approx 5 min to complete, so will save results to medsum_full4.RData file and read it into memory
# the line will be commented to save time
change_null_to_list <- function(x) {
  # If x is a data frame, do nothing and return x
  # Otherwise, return a data frame with 1 row of NAs
  if (!is.null(x)) {return(x)}
  else {return(as_tibble(t(c(BPmedgroup=as.factor(NA), DLYFREQ=as.numeric(NA), DLYDOSE=as.numeric(NA),std_dose=as.numeric(NA),DDI=as.numeric(NA),compliance=as.numeric(NA)))))}
}

# for (i in 9:58) {test[[names(test)[i]]]<-lapply(test[[names(test)[i]]],change_null_to_list)}

# cardio data organization
# overview of data...
# display the time differences between ABPM and casual BP measurement dates
# cardio %>% filter(BPstatus!=-1) %>% mutate(date_dif=ABPM_DATE-DB_DATE) %>% group_by(VISIT) %>% summarise_at(vars(date_dif),funs(mean_dif_days=365*mean(.,na.rm=TRUE),max_dif_yr=max, sd_dif_yr=sd(.,na.rm=TRUE)))
# display the number of observations where the difference between the date of visit and the ABPM date is > 60 days (=TRUE)
# consider excluding these observations since the BP status is less reliable
# cardio %>% filter(BPstatus!=-1) %>% mutate(date_dif=abs(ABPM_DATE-DB_DATE)) %>% group_by(VISIT) %>% mutate(big_dif=date_dif>60/365) %>% select(big_dif) %>% table %>% addmargins
# display the counts for all the BP statuses grouped by visit
# patients with BP status -1 occured when casual BP was not classified (either no BP entered at that visit, or BP entered, but no percentile/z-score)
# cardio %>% group_by(VISIT) %>% filter(ABPMSUCCESS==1, BPstatus !=-1) %>% select(BPstatus) %>% table() %>% addmargins
# display information about timing of ABPM data, grouped by visits
# cardio%>% filter(ABPMSUCCESS==1) %>% group_by(VISIT) %>% summarise_at(vars(ABPM_DATE), funs(date_mean=mean,date_min=min,date_max=max,n_pts=length))

# add columns for new 2017 AAP BP percentiles and z-scores (this takes about 2 minutes to calculate)
source('BPz/bpzv2.R')
test$SBPPCTAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,1,test$SBP, z=F)
test$DBPPCTAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,2,test$DBP, z=F)
test$SBPZAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,1,test$SBP, z=T)
test$DBPZAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,2,test$DBP, z=T)

source('BPstatus4th.R')
source('BPstatus2017.R')
source('bpp4.R')
# classification of BP status based on 4th report and 2017 guidelines
# note: BP status is NA when BP was unavailable
test$BPstatus2017<-mapply(BPstatus2017,test$SBP,test$DBP,test$age,test$MALE1FE0,test$AVHEIGHT)
# fill in missing height percentiles to reduce number of missing BP status
missing_htpct<-which(is.na(test$HTPCTAG) & !is.na(test$AVHEIGHT))
test$HTPCTAG[missing_htpct]<- mapply(htz,test$AVHEIGHT[missing_htpct], test$age[missing_htpct], test$MALE1FE0[missing_htpct],p=T)
test$BPstatus4th<-mapply(bpstatus4th,test$SBP,test$DBP,test$age,test$MALE1FE0,test$HTPCTAG,test$SBPPCTAGH,test$DBPPCTAGH)

# save(test,file="test.RData")
load(file="test.RData")
# stacked barplot of BP status (in proportions) grouped by visit
test %>% filter(VISIT%%10==0) %>% group_by(VISIT) %>% select(BPstatus2017) %>% table(useNA = 'always') %>% prop.table(margin=1)

# cardio data organization
# analysis of BP status
# based on 2009 AHA consensus guidelines (old)
# BP status is coded as a variable called BPstatus in cardio
# 0 = normal
# 1 = White coat hypertension (WCH)
# 2 = masked hyeprtension (MH)
# 3 = ambulatory hypertension (AH)

# cardio$BPstatus=-1
# cardio$BPstatus[cardio$WKSYSINDX<0.95 & cardio$WKDIAINDX<0.95 & cardio$SLSYSINDX<0.95 & cardio$SLDIAINDX<0.95 & cardio$WKSYSLOAD<25 & cardio$WKDIALOAD<25 & cardio$SLSYSLOAD<25 & cardio$SLDIALOAD<25 & cardio$SBPPCTAGH<95 & cardio$DBPPCTAGH<95]<-0
# cardio$BPstatus[(cardio$WKSYSINDX<0.95 & cardio$WKDIAINDX<0.95 & cardio$SLSYSINDX<0.95 & cardio$SLDIAINDX<0.95 & cardio$WKSYSLOAD<25 & cardio$WKDIALOAD<25 & cardio$SLSYSLOAD<25 & cardio$SLDIALOAD<25) & (cardio$SBPPCTAGH>=95 | cardio$DBPPCTAGH>=95)]<-1
# cardio$BPstatus[(cardio$WKSYSINDX>=0.95 | cardio$WKDIAINDX>=0.95 | cardio$SLSYSINDX>=0.95 | cardio$SLDIAINDX>=0.95 | cardio$WKSYSLOAD>=25 | cardio$WKDIALOAD>=25 | cardio$SLSYSLOAD>=25 | cardio$SLDIALOAD>=25) & (cardio$SBPPCTAGH<95 & cardio$DBPPCTAGH<95)]<-2
# cardio$BPstatus[(cardio$WKSYSINDX>=0.95 | cardio$WKDIAINDX>=0.95 | cardio$SLSYSINDX>=0.95 | cardio$SLDIAINDX>=0.95 | cardio$WKSYSLOAD>=25 | cardio$WKDIALOAD>=25 | cardio$SLSYSLOAD>=25 | cardio$SLDIALOAD>=25) & (cardio$SBPPCTAGH>=95 | cardio$DBPPCTAGH>=95)]<-3
source('bpclass.R')
cardio<-cardio %>%rowwise()%>% mutate(noctHTN=sum(SYSPCTDIPPING<10,DIAPCTDIPPING<10))
cardio<-left_join(cardio,test %>% select(CASEID,VISIT,SBPPCTAGH2017,DBPPCTAGH2017, age,MALE1FE0,AVHEIGHT))
cardio<-cardio %>% rowwise() %>% mutate(BPclass= bpclass(WKSYSINDX,WKDIAINDX,SLSYSINDX,SLDIAINDX,WKSYSLOAD,WKDIALOAD,SLSYSLOAD, SLDIALOAD, SBP, DBP, SBPPCTAGH2017,DBPPCTAGH2017))
test<-full_join(cardio %>% select(CASEID,VISIT,ABPMSUCCESS, BPclass,noctHTN),test)
# combine BP class variants into larger groups using BP classification from 2014
test$BPclass.factor<-cut(test$BPclass,0:6,right=FALSE, labels=c("NL","WCH","PH","MH","AH","SAH"),ordered_result = TRUE)
# stacked barplot of BP class (in proportions) grouped by visit
table.a<-test %>% filter(VISIT%%20==0)%>% group_by(VISIT) %>% select(BPclass.factor) %>% table(useNA = 'no') 
table.a%>% plot(prop.table(margin=1)[,1], col=c("green","greenyellow","yellow","orange","orangered","red4"),main="Proportion of BP classes by visit", ylab="BP class (%)")
# proportion table of nocturnal HTN by BP class grouped by visit
table.b<-test %>% filter(VISIT%%20==0) %>% group_by(BPclass.factor) %>% select(noctHTN,VISIT) %>% table(useNA = 'no') %>% prop.table(margin=3)*100
# to save time, use this as starting point...
# duplicates<-test %>% group_by(CASEID,VISIT) %>% mutate(duplicates=n_distinct(MSVISDAT)) %>% filter(duplicates>1)
#save(test,duplicates, file="test.RData")
load(file="test.RData")

# table to demonstrate differences in classification of patients during visit 20
BPstatus_comparison<-test %>% filter(VISIT==20) %>% select(BPstatus4th,BPstatus2017) %>% table(.,useNA = "always") %>% addmargins()
write.table(BPstatus_comparison,row.names=T, col.names=NA,"BPstatus_comparison.csv")
BPstatus_comparison

# data overview found in file "overview of data"
# this is a test to see how changes as handled in git

# relationship between BP status and number of antihypertensive agents:
# reference: tutorial of chi square analysis in r: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
n_agents_20<-test %>% filter(VISIT==20)
# show the table
# excludes patients whose BP status unknown (either clinic BP or ABPM study results are missing)
n_agents_20 %>% select(n_agents,BPclass.factor) %>% table(., useNA = 'no') %>% addmargins()
# balloonplot
dt <- table(n_agents_20$n_agents,n_agents_20$BPclass.factor, exclude = c(-1,NA))
balloonplot(t(dt), main ="relationship between BP status and number of antihypertensive agents used at visit 20", xlab ="BP status", ylab="n_agents",label = FALSE, show.margins = FALSE)

# chi-squared test shows highly significant p-value of 8.4E-5
# suggesting n_agents and BP status are not randomly/evenly distributed
chisq <- chisq.test(dt)
chisq
library(corrplot,ggplot2)
corrplot(chisq$residuals, is.cor = FALSE)

# writing table files
#table1
# visit20_BPstatus_n_agents.csv: BP status (columns) by number of agents (rows) in VISIT 20 only
# note: includes those with unknown BP status (-1), due to unknown clinic BP percentile, or unsuccessfull ABPM study
table1<-test %>% filter(VISIT==20) %>% select(n_agents,BPclass.factor) %>% table() %>% addmargins()
write.table(table1,row.names=T, col.names=NA, "visit20_BPstatus_n_agents.csv")
table1

# n_agents_visit.csv: n_agents(columns) by visit (rows)
# note: includes those not taking any antihypertensive meds
table2<-test %>% filter(VISIT%%10==0) %>% select(VISIT,n_agents) %>% table()
write.table(table2, row.names=T, col.names=NA,"n_agents_visit.csv")
table2
# proportion table version of table 2
round(100*prop.table(table2,margin=1),1)

# BPstatus_visit.csv: BPclass (columns) by visit (rows, only even visits when ABPM obtained)
table3<-test %>% filter(VISIT%%20==0) %>% group_by(VISIT) %>% select(BPclass.factor) %>% table()
write.table(table3,row.names=T, col.names = NA, "BPstatus_visit.csv")
table3 %>% addmargins(2)
round(100*prop.table(table3,margin=1),1) %>% addmargins(2)

# Counts of patients on each antihypertensive medication grouped by visit and sorted (descending)
unsorted<-medsum_full.3 %>% filter(VISIT%%10==0) %>% group_by(VISIT,med.corrected) %>% summarise(n_rx=n_distinct(CASEID)) %>% arrange(desc(n_rx),.by_group=T) %>% xtabs(formula=n_rx~addNA(factor(med.corrected))+VISIT)
table4<-unsorted[sort(unsorted[,1], decreasing=T,index.return=T)$ix,] %>% addmargins()
write.table(table4,row.names=T, col.names=NA,"visit_medcounts.csv")
table4

# counts of missing BP measurements
# 1 = one of SBP or DBP are missing
# 2 = both SBP and DBP are missing
# 0 = SBP and DBP are measured
table(is.na(test$SBP)+is.na(test$DBP))

# parse drug info stored in test as separate variables for analysis
source("get_drug_info.R")
drugname<-names(test)[9:58]
drugname<-drugname[-c(3,5,6,11,20,26,33)]
temp <- paste(drugname[1:43], "<- get_drug_info(test,'",drugname,"', c(1:6))", sep="")
eval(parse(text=temp))

temp.1<-paste("test$",drugname[1:43],"_DDI<-",drugname[1:41],"$DDI", sep="")
temp.2<-paste("test$",drugname[1:43],"_std_dose<-",drugname[1:41],"$std_dose", sep="")
eval(parse(text=temp.1))
eval(parse(text=temp.2))

# draw histograms of DDI for each drug
temp<-paste("if (!all(is.na(",drugname[1:43],"$DDI))) hist(",drugname[1:43],"$DDI, xlim=c(0,3), breaks=200, main=paste(\"",drugname[1:43],"\"))", sep="")
par(mfrow=c(6,4))
par(mar=c(2,2,1,1))
eval(parse(text=temp))

# comparing patients based on DDI
temp_DDI<-test%>% select(ends_with('_DDI')) %>% as.matrix %>% apply(.,1,sort,decreasing=T,na.last=T)
test$DDI_all<-as_tibble(t(temp_DDI))


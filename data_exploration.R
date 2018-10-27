#loading libraries
library(readr)
library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library("histogram", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")

setwd("~/Documents/Research/CKiDgit")

# importing data
medsum_full <- read.csv("data/medsum_full.csv")
cardio<-read.csv("data/cardio.csv")
echo<-read.csv("data/echo.csv")
f13 <- read.csv("~/Documents/Research/CKiDr_local/data/f13.csv")
f14 <- read.csv("~/Documents/Research/CKiDr_local/data/f14.csv")
f15 <- read.csv("~/Documents/Research/CKiDr_local/data/f15.csv")
socdem <- read.csv("~/Documents/Research/CKiDr_local/data/socdem.csv")
growth <- read.csv("~/Documents/Research/CKiDr_local/data/growth.csv")
pe <- read.csv("~/Documents/Research/CKiDr_local/data/pe.csv")
gfrcalibratedsummary <- read.csv("~/Documents/Research/CKiDr_local/data/gfrcalibratedsummary.csv")
advr <- read.csv("~/Documents/Research/CKiDr_local/data/advr.csv")
kidhist <- read.csv("~/Documents/Research/CKiDr_local/data/kidhist.csv")
l05 <- read.csv("~/Documents/Research/CKiDr_local/data/l05.csv")
l02 <- read.csv("~/Documents/Research/CKiDr_local/data/l02.csv")
l51 <- read.csv("~/Documents/Research/CKiDr_local/data/l51.csv")
l31 <- read.csv("~/Documents/Research/CKiDr_local/data/l31.csv")
gfr<-gfrcalibratedsummary
gh <- read.csv("~/Documents/Research/CKiDr_local/data/gh.csv")

# data integrity issues:
# review of the raw data showed several instances where the historical data is not consistent
# eg, changes in LBW status (growth.csv) in the same participant across different visits (this should be the same)
# this is also seen for premature status (growth.csv)
# the table generated below shows the flagged BP counts grouped by visit (ie, problem measuring BP)
# consider excluding those with status '1'
pe %>% select(VISIT, PEPRBLPM)%>% group_by(VISIT) %>% table %>% addmargins()
# growth appears to contain the best weights after comparing data from pe and growth
#l31<-rename(l31, "SMWBAGEP.2"=SMWBAGEP)

echo_df<-echo %>% select(CASEID,VISIT,ECHODATEY,LVMI,LVHF,LVHE,CAUTION)
# loads the correct BP medication names into BP_medlist
# and corrects BP medication names for brand names and misspellings
# corrected names are assigned a new variable called med.corrected
# display information about timing of echocardiogram data, grouped by visits
echo_df %>% group_by(VISIT) %>% filter(CAUTION==0) %>% summarise_at(vars(ECHODATEY), funs(date_mean=mean,date_min=min,date_max=max, n_pts=length))

load("BP_medlist.RData")
# note: BP_medlist will display the misspelled/brand names of BP meds, frequency at visit 10?, and corrected names
medsum_full$med.corrected<-BP_medlist$Corrected.Name[match(medsum_full$MSMEDICA,BP_medlist$Var1)]
# categorize the BP meds into groups based on BP_medgroups.csv
BP_medgroups <- read.csv("~/Documents/Research/CKiDr_local/BP_medgroups.csv")
medsum_full$BPmedgroup<-BP_medgroups$bp_group[match(medsum_full$med.corrected,BP_medgroups$x)]
# add a column containing the number of distinct antihypertensive agents (n_agents) taken by a patient at that visit
medsum_full<-medsum_full %>% filter(!is.na(BPmedgroup)) %>% group_by(CASEID,VISIT) %>% mutate(n_agents=n_distinct(med.corrected))

# displays the number of patients at each visit in the medsum_full data
medsum_full %>% select(CASEID,VISIT) %>% group_by(VISIT) %>% distinct() %>% tally

# displays the number of patients taking BP meds grouped by visit and number of agents
medsum_full %>% select(CASEID,VISIT,n_agents) %>% group_by(VISIT) %>% select(n_agents) %>% table()

# generates and displays a column for the number of antihypertensive agents for each patient at each visit
medsum_full %>% filter(!is.na(BPmedgroup)) %>% group_by(CASEID,VISIT) %>% summarise(n_agents=n_distinct(med.corrected))
                                                                                    
# display information about timing of ABPM data, grouped by visits
cardio%>% filter(ABPMSUCCESS==1) %>% group_by(VISIT) %>% summarise_at(vars(ABPM_DATE), funs(date_mean=mean,date_min=min,date_max=max,n_pts=length))

# analysis of BP status
# based on 2009 AHA consensus guidelines
# BP status is coded as a variable called BPstatus in cardio
# 0 = normal
# 1 = White coat hypertension (WCH)
# 2 = masked hyeprtension (MH)
# 3 = ambulatory hypertension (AH)
cardio$BPstatus=-1
cardio$BPstatus[cardio$WKSYSINDX<0.95 & cardio$WKDIAINDX<0.95 & cardio$SLSYSINDX<0.95 & cardio$SLDIAINDX<0.95 & cardio$WKSYSLOAD<25 & cardio$WKDIALOAD<25 & cardio$SLSYSLOAD<25 & cardio$SLDIALOAD<25 & cardio$SBPPCTAGH<95 & cardio$DBPPCTAGH<95]<-0
cardio$BPstatus[(cardio$WKSYSINDX<0.95 & cardio$WKDIAINDX<0.95 & cardio$SLSYSINDX<0.95 & cardio$SLDIAINDX<0.95 & cardio$WKSYSLOAD<25 & cardio$WKDIALOAD<25 & cardio$SLSYSLOAD<25 & cardio$SLDIALOAD<25) & (cardio$SBPPCTAGH>=95 | cardio$DBPPCTAGH>=95)]<-1
cardio$BPstatus[(cardio$WKSYSINDX>=0.95 | cardio$WKDIAINDX>=0.95 | cardio$SLSYSINDX>=0.95 | cardio$SLDIAINDX>=0.95 | cardio$WKSYSLOAD>=25 | cardio$WKDIALOAD>=25 | cardio$SLSYSLOAD>=25 | cardio$SLDIALOAD>=25) & (cardio$SBPPCTAGH<95 & cardio$DBPPCTAGH<95)]<-2
cardio$BPstatus[(cardio$WKSYSINDX>=0.95 | cardio$WKDIAINDX>=0.95 | cardio$SLSYSINDX>=0.95 | cardio$SLDIAINDX>=0.95 | cardio$WKSYSLOAD>=25 | cardio$WKDIALOAD>=25 | cardio$SLSYSLOAD>=25 | cardio$SLDIALOAD>=25) & (cardio$SBPPCTAGH>=95 | cardio$DBPPCTAGH>=95)]<-3

# display the time differences between ABPM and casual BP measurement dates
cardio %>% filter(BPstatus!=-1) %>% mutate(date_dif=ABPM_DATE-DB_DATE) %>% group_by(VISIT) %>% summarise_at(vars(date_dif),funs(mean_dif_days=365*mean(.,na.rm=TRUE),max_dif_yr=max, sd_dif_yr=sd(.,na.rm=TRUE)))
# display the number of observations where the difference between the date of visit and the ABPM date is > 60 days (=TRUE)
# consider excluding these observations since the BP status is less reliable
cardio %>% filter(BPstatus!=-1) %>% mutate(date_dif=ABPM_DATE-DB_DATE) %>% group_by(VISIT) %>% mutate(big_dif=date_dif>60/365) %>% select(big_dif) %>% table %>% addmargins
# display the counts for all the BP statuses grouped by visit
# patients with BP status -1 occured when casual BP was not classified (either no BP entered at that visit, or BP entered, but no percentile/z-score)
cardio %>% group_by(VISIT) %>% filter(ABPMSUCCESS==1, BPstatus !=-1) %>% select(BPstatus) %>% table() %>% addmargins

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
test<-full_join(cardio %>% select(CASEID,VISIT,DB_DATE,SBP,SBPINDXAGH,SBPPCTAGH,SBPZAGH,DBP,DBPINDXAGH,DBPPCTAGH,DBPZAGH,ABPM_DATE,ABPMSUCCESS,BPstatus),test)

# this was compared with age from pierce et al and correlates well, but not perfectly
# note: age from pierce is in integers, not decimals
test$age<-(test$BSDATE-test$DOB)+test$DB_DATE

# gender: since some values of MALE1FE0 are missing, will solve missing values by taking the mean of MALE1FE0 for each case (no change in gender over time)
test<-test %>% group_by(CASEID) %>% mutate(male1fe0=mean(MALE1FE0, na.rm=TRUE))


# organize medication data

#old: medsum_full.2<- medsum_full.2 %>% group_by(CASEID,VISIT,MSVISDAT,med.corrected, DLYFREQ,BPmedgroup,n_agents) %>% summarise_at(vars(DLYDOSE),sum)
# old: medsum_full.2<- medsum_full.2 %>% group_by(CASEID,VISIT,MSVISDAT,med.corrected, n_agents,BPmedgroup,DLYFREQ) %>% summarise_at(vars(DLYDOSE),sum) %>% ungroup() %>% nest(c(6:8), .key="rx_info")
# old version: medsum_full.3<-dcast(medsum_full.2,CASEID+VISIT+MSVISDAT+DLYFREQ+BPmedgroup+n_agents~med.corrected, value.var="DLYDOSE")
#old: medsum_full.4<-medsum_full.3 %>% group_by(CASEID,VISIT) %>% summarise_at(c(6:57),mean, na.rm=TRUE)
# old: calculate mg/kg for each drug
#test<-test %>% mutate_at(.vars=c(5:55),funs(mg_kg=if (!is.null(.[[]])) unlist(.)[3]/AVWEIGHT else 0))

# this line generates columns for std_dose (dose per kg, with AVWEIGHT from growth data), corrects
# for split dosing (eg, labetalol), and nests drug data into a 4 column list called rx_info
medsum_full.2<-medsum_full
medsum_full.3<-full_join(test %>% select(CASEID,VISIT,AVWEIGHT),medsum_full.2) %>% select(CASEID,VISIT,MSVISDAT,med.corrected,n_agents,BPmedgroup,DLYFREQ,DLYDOSE,AVWEIGHT) %>% mutate(std_dose=DLYDOSE/AVWEIGHT) %>% group_by(CASEID,VISIT,MSVISDAT,med.corrected, n_agents,AVWEIGHT,BPmedgroup,DLYFREQ) %>% summarise_at(vars(DLYDOSE),sum) %>% ungroup() %>% mutate(std_dose=DLYDOSE/AVWEIGHT)
medsum_full.4<-medsum_full.3 %>% nest(c(7:10), .key="rx_info")
medsum_full.4<-dcast(medsum_full.4,CASEID+VISIT+MSVISDAT+n_agents~med.corrected, value.var="rx_info")

test<-full_join(medsum_full.4,test)


#add LVMI/BSA variable
#test<- test %>% mutate(LVMI_BSA=LVMI/BSA)

#add ckidfull [the best estimated gfr for research based on Schwartz and Schneider 2012]
test<-test %>% mutate(term1=((AVHEIGHT/100)/SCR)^0.456)
test<-test %>% mutate(term2=(1.8/CYC_DB)^0.418)
test<-test%>% mutate(term3=(30/BUN)^0.079)
test<-test %>% mutate(term4=((AVHEIGHT/100)/1.4)^0.179)
test<-test %>% mutate(ckidfull=39.8*term1*term2*term3*term4*if_else(male1fe0==1,1.076,1))

# add column for urine protein:creatinine ratio based on RLURPROT and RLURCREA
# add column for percentage of total urine protein that is albumin (Ualb_pct)
test %>% mutate(Upc=RLURPROT/RLURCREA)
test %>% mutate(Ualb_pct=RLURMALB/RLURPROT)

# add a column for CKD stage using ckidfull
test$CKD_stage <- cut(test$ckidfull, 
                      breaks = c(-Inf, 15, 30, 60, 90, Inf), 
                      labels=c(5:1),ordered_result=TRUE, 
                      right = FALSE)

# add columns for new 2017 AAP BP percentiles and z-scores (this takes about 2 minutes to calculate)
source('BPz/bpzv2.R')
test$SBPPCTAGH2017<-mapply(bpp,test$SBP,test$age, test$AVHEIGHT,test$male1fe0,1,1)
test$DBPPCTAGH2017<-mapply(bpp,test$DBP,test$age, test$AVHEIGHT,test$male1fe0,2,1)
test$SBPZAGH2017<-mapply(bpp,test$SBP,test$age, test$AVHEIGHT,test$male1fe0,1,2)
test$DBPZAGH2017<-mapply(bpp,test$DBP,test$age, test$AVHEIGHT,test$male1fe0,2,2)


# data overview found in file "overview of data"



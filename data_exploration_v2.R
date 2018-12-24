#loading libraries
library(readr)
library(corrplot)
library(dplyr)
library(reshape2)
library(histogram)
library(gplots)
library(tidyr)
library(stats)
library(tibble)
library(ggplot2)
detach("package:MASS", unload=TRUE)

wd<-if_else(Sys.info()["nodename"]=='matta', "~/Documents/Research/CKiDgit", "C:/Users/Orit/Downloads/CKiD/CKiDgit")
setwd(wd)

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
gfr <- read.csv("data/gfrcalibratedsummary.csv")
advr <- read.csv("data/advr.csv")
kidhist <- read.csv("data/kidhist.csv")
l05 <- read.csv("data/l05.csv")
l02 <- read.csv("data/l02.csv")
l51 <- read.csv("data/l51.csv")
l31 <- read.csv("data/l31.csv")
gh <- read.csv("data/gh.csv")

# data integrity issues:
# review of the raw data showed several instances where the historical data is not consistent
# eg, changes in LBW status (growth.csv) in the same participant across different visits (this should be the same)
# this is also seen for premature status (growth.csv)
# the table generated below shows the flagged BP counts grouped by visit (ie, problem measuring BP)
# consider excluding those with status '1'
# pe %>% select(VISIT, PEPRBLPM)%>% group_by(VISIT) %>% table %>% addmargins()
# growth appears to contain the best weights after comparing data from pe and growth

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

source("gfr.R")
# gfr will be based on ckidfull equation if the data is available, otherwise bedside schartz equation is used
test$gfr<-mapply(GFR,test$AVHEIGHT,test$SCR,test$CYC_DB,test$BUN,test$MALE1FE0, permissive=TRUE)

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
load("BP_medlist.RData",verbose = TRUE)
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
medsum_full<-left_join(test %>% select(CASEID,VISIT,AVWEIGHT,age,gfr),medsum_full)
medsum_full<-medsum_full %>% mutate(std_dose=ifelse(DLYDOSE>0,DLYDOSE/AVWEIGHT,NA))
source("DDI.R")
medsum_full$ddi<-mapply(DDI,medsum_full$DLYDOSE,medsum_full$med.corrected,medsum_full$age,medsum_full$AVWEIGHT,medsum_full$gfr)

# identify duplicates (cases with visit codes that have coresponding visit dates that are not the same) and eliminate all the corresponding medication data
# justifications 1) unable to determine whether duplicate visits represent 'addition' of medications or 'switching' between meds
# 2) some cases have same visit code and corresponding visit dates separated in time by up to 4.28 years (obviously an error in recording data)
medsum_full<-medsum_full %>% group_by(CASEID,VISIT) %>% mutate(ndup=n_distinct(MSVISDAT))
medsum_full<-medsum_full %>%filter(ndup==1) # remove all duplicates from the data
medsum_full<-medsum_full %>% group_by(CASEID,VISIT) %>% mutate(n_agents=n_distinct(med.corrected, na.rm=TRUE))
medsum_full<-medsum_full %>% group_by(CASEID,VISIT) %>% mutate(sum_DDI=sum(ddi))
medsum_full.old<-medsum_full
medsum_full<-medsum_full %>% ungroup %>% nest(.,BPmedgroup,DLYFREQ,DLYDOSE,std_dose,ddi,compliance, .key="rx_info")
medsum_full<-dcast(medsum_full,CASEID+VISIT+MSVISDAT+n_agents+mean_compliance+sum_DDI~med.corrected, value.var="rx_info")
test<-full_join(medsum_full,test)

# replace all NULL values of drugs with a tibble containing NA's using change_null_to_list function above
# note: the line below takes approx 5 min to complete, so will save results to medsum_full4.RData file and read it into memory
# the line will be commented to save time
change_null_to_list <- function(x) {
  # If x is a data frame, do nothing and return x
  # Otherwise, return a data frame with 1 row of NAs
  if (!is.null(x)) {return(x)}
  else {return(as_tibble(t(c(BPmedgroup=as.factor(NA), DLYFREQ=as.numeric(NA), DLYDOSE=as.numeric(NA),std_dose=as.numeric(NA),ddi=as.numeric(NA),compliance=as.numeric(NA)))))}
}

for (i in which(names(test)=="ACETAZOLAMIDE"):which(names(test)=="VERAPAMIL")) {test[[names(test)[i]]]<-lapply(test[[names(test)[i]]],change_null_to_list)}

# parse drug info stored in test as separate variables for analysis
source("get_drug_info.R")
drugname<-names(test)[which(names(test)=="ACETAZOLAMIDE"):which(names(test)=="VERAPAMIL")]
# exclude combo drugs (contain '/', eg "AMLODIPINE/HCTZ")
drugname<-drugname[-grep("\\/",drugname)]
temp <- paste(drugname[1:43], "<- get_drug_info(test,'",drugname,"', c(1:6))", sep="")
eval(parse(text=temp))

temp.1<-paste("test$",drugname[1:43],"_ddi<-",drugname[1:43],"$ddi", sep="")
temp.2<-paste("test$",drugname[1:43],"_std_dose<-",drugname[1:43],"$std_dose", sep="")
eval(parse(text=temp.1))
eval(parse(text=temp.2))



# add columns for new 2017 AAP BP percentiles and z-scores (this takes about 2 minutes to calculate)
source('bpp4.R')
source('BPz/bpzv2.R')
source('BPstatus4th.R')
source('BPstatus2017.R')
source('bpclass.R')
source('bpclass2.R')
test$SBPPCTAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,1,test$SBP, z=F)
test$DBPPCTAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,2,test$DBP, z=F)
test$SBPZAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,1,test$SBP, z=T)
test$DBPZAGH2017<-mapply(bpp,test$age,test$AVHEIGHT,test$MALE1FE0,2,test$DBP, z=T)

# classification of BP status based on 4th report and 2017 guidelines
# note: BP status is NA when BP was unavailable
test$BPstatus2017<-mapply(BPstatus2017,test$SBP,test$DBP,test$age,test$MALE1FE0,test$AVHEIGHT)
# fill in missing height percentiles to reduce number of missing BP status
missing_htpct<-which(is.na(test$HTPCTAG) & !is.na(test$AVHEIGHT))
test$HTPCTAG[missing_htpct]<- mapply(htz,test$AVHEIGHT[missing_htpct], test$age[missing_htpct], test$MALE1FE0[missing_htpct],p=T)
test$BPstatus4th<-mapply(bpstatus4th,test$SBP,test$DBP,test$age,test$MALE1FE0,test$HTPCTAG,test$SBPPCTAGH,test$DBPPCTAGH)

# cardio data organization
# analysis of BP status
# based on 2009 AHA consensus guidelines (old)
# BP status is coded as a variable called BPstatus in cardio
# 0 = normal
# 1 = White coat hypertension (WCH)
# 2 = masked hyeprtension (MH)
# 3 = ambulatory hypertension (AH)

# noctHTN variable: BP% dipping <10%
cardio<-cardio %>%rowwise()%>% mutate(noctHTN=sum(SYSPCTDIPPING<10,DIAPCTDIPPING<10))
cardio<-left_join(cardio,test %>% select(CASEID,VISIT,SBPPCTAGH2017,DBPPCTAGH2017, age,MALE1FE0,AVHEIGHT))
cardio<-cardio %>% rowwise() %>% mutate(BPclass= bpclass2(WKSYSINDX,WKDIAINDX,SLSYSINDX,SLDIAINDX,WKSYSLOAD,WKDIALOAD,SLSYSLOAD, SLDIALOAD, SBP, DBP, SBPPCTAGH2017,DBPPCTAGH2017))
test<-full_join(cardio %>% select(CASEID,VISIT,ABPMSUCCESS,BPclass,noctHTN),test)
# combine BP class variants into larger groups using BP classification from 2014
test$BPclass.factor<-cut(test$BPclass,0:4,right=FALSE, labels=c("NL","WCH","MH","AH"),ordered_result = TRUE)
test$BPclass.2.factor<-factor(test$BPclass.factor<="WCH", labels=c("NL+WCH","MH+AH"))
# Create variable LVMIp using LVMIp function to calculate LVMI percentiles
test<-test %>% filter(!is.na(LVMI)) %>% rowwise() %>% mutate(LVMIp=LVMI_p(MALE1FE0,age,LVMI))
test$BMIz<-qnorm(test$BMIPCTAG/100)
test$GNGDIAG.factor<-factor(test$GNGDIAG<=2,labels=c("glom","non-glom"))

# save(test,file="test2.RData")
load(file="test.RData")
### data organization complete ###


### data analysis begin ###
# define cohort
combo_cases<-medsum_full.old %>% filter(VISIT==20,length(grep("\\/",med.corrected)!=0L)) %>% ungroup() %>% select(CASEID) %>% unique() %>% unlist() %>% as.numeric()
cases<-test %>% filter(!CASEID %in% combo_cases,VISIT==20, n_agents>0,!is.na(sum_DDI),ABPMSUCCESS==1,!BPclass.factor=="WCH",!is.na(BUN), !is.na(CYC_DB),!is.na(SCR),!is.na(AVHEIGHT), !is.na(AVWEIGHT), !is.na(BPstatus2017)) %>% ungroup() %>%select(CASEID) %>% unique %>% unlist %>% as.numeric()# exclude those taking combo medications
test.cohort<-test %>% filter(VISIT==20, CASEID %in% cases)

test.cohort<-test.cohort %>% mutate(sum_DDI.ln=log(sum_DDI))


cum_DDI.normotensive<-test.cohort %>% filter(BPclass.factor<="WCH") %>% ungroup %>% select(sum_DDI.ln) %>% unlist %>% as.numeric()
cum_DDI.hypertensive<-test.cohort %>% filter(BPclass.factor>"WCH") %>% ungroup %>% select(sum_DDI.ln) %>% unlist %>% as.numeric()




# timing of PE to ABPM
# library(ggalt)
# theme_set(theme_classic())
# test.cohort$CASEID.factor=factor(test.cohort$CASEID, levels=as.char(test.cohort$CASEID))
# gg<-ggplot(test.cohort,aes(x=MSVISDAT,xend=ABPM_DATE,y=CASEID.factor,group=CASEID.factor))+
#   geom_dumbbell(color="#a3c4dc",
#                size=0.75,
#                point.colour.1="#0e668b")+
#   scale_x_continuous(label=percent)+
#   labs(x=NULL,
#        y=NULL,
#        title="Time differences: PE vs ABPM",
#        caption="CKiD")+
#   theme(plot.title=element_text(hjust=0.5,face="bold"))
# plot(gg)


library(epitab)
library(neat_table)



# stacked barplot of BP class (in proportions) grouped by visit
table.a.2<-test %>% filter(VISIT%%20==0)%>% group_by(VISIT) %>% select(BPclass2.factor) %>% table(useNA = 'no') 
par(mfrow=c(1,1))
table.a.2%>% plot(prop.table(margin=1)[,1], col=c("green","greenyellow","yellow","orange","orangered","red4"),main="Proportion of BP classes by visit", ylab="BP class (%)")
table.a.2

# proportion table of nocturnal HTN by BP class grouped by visit
table.b<-test %>% filter(VISIT%%20==0) %>% group_by(BPclass.factor) %>% select(noctHTN,VISIT) %>% table(useNA = 'no')
table.b %>% addmargins
table.b.2<-test %>% filter(VISIT%%20==0) %>% group_by(BPclass2.factor) %>% select(noctHTN,VISIT) %>% table(useNA = 'no')
table.b.2 %>% addmargins

# table to demonstrate differences in classification of patients during visit 20
BPstatus_comparison<-test.cohort %>% select(BPstatus4th,BPstatus2017) %>% table(.,useNA = "ifany") %>% addmargins()
write.table(BPstatus_comparison,row.names=T, col.names=NA,"BPstatus_comparison.csv")
BPstatus_comparison

# data overview found in file "overview of data"
# this is a test to see how changes as handled in git

# relationship between BP status and number of antihypertensive agents:
# reference: tutorial of chi square analysis in r: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
# show the table
# excludes patients whose BP status unknown (either clinic BP or ABPM study results are missing)
temp<-test.cohort %>% select(n_agents,BPclass2.factor) %>% table(., useNA = 'no') %>% addmargins()
write.table(temp,row.names=T,col.names=NA,"BPclass_by_n_agents_v20.csv")

# balloonplot (association between n_agents vs BPclass2)
dt <- table(test.cohort$n_agents,test.cohort$BPclass2.factor, exclude = c(-1,NA))
balloonplot(t(dt), main ="relationship between BP status and number of antihypertensive agents used at visit 20", xlab ="BP status", ylab="n_agents",label = FALSE, show.margins = FALSE)

# chi-squared test shows highly significant p-value of 8.4E-5
# suggesting n_agents and BP status are not randomly/evenly distributed
chisq <- chisq.test(dt)
chisq
library(corrplot,ggplot2)
corrplot(chisq$residuals, is.cor = FALSE, title="chi-squared test residuals for\n BP class vs n_agents\n (X-squared = 27.946, df = 12, p-value = 0.005634)",addCoef.col = 1,mar=c(1,0,2.5,1),number.cex=0.6)

# insights from this plot - demonstrates where n_agents and BP class deviates from the expected, if no relationship existed between these factors:
# 1. n_agents positively correlates with AH group (ie, subjects with AH tend to have more anti-HTN agents)
# 2. n_agents negatively correlated with MH group (ie, fewer than expected subjects with 3-4 agents, and more than expected subjects with 0 agents - suggests possible undertreatment - consistent with Barletta et al)
# 3. n_agents positively correlates with WCH group (more than expected subjects with 3 and 4 agents, fewer than expected subjects with 1-2 agents - suggests possible overtreatment - not previously reported, but makes sense)
# 4. n_agents does not have consistent pattern within the normotensive group (ie, more than expected on 1 agent, fewer than expected with 0,2-4 agents)

chisq.2 <- chisq.test(dt.2)
chisq.2
library(corrplot)
corrplot(chisq.2$residuals, is.cor = FALSE, title="chi-squared test residuals for\n BP class vs n_agents",addCoef.col = 1,mar=c(1,0,2.5,1),number.cex=0.6)

# writing table files
#table1
# visit20_BPstatus_n_agents.csv: BP status (columns) by number of agents (rows) in VISIT 20 only
# note: includes those with unknown BP status (-1), due to unknown clinic BP percentile, or unsuccessfull ABPM study
test$n_agents.factor<-as.factor(test$n_agents)
test$VISIT.factor<-as.factor(test$VISIT[test$VISIT%%10==0])
table1<-test %>% filter(VISIT==20) %>% select(n_agents,BPclass.factor) %>% table() %>% addmargins()
write.table(table1,row.names=T, col.names=NA, "visit20_BPstatus_n_agents.csv")
table1

test$n_agents.factor<-as.factor(test$n_agents)
test$VISIT.factor<-as.factor(test$VISIT[test$VISIT%%10==0])
table1.2<-test %>% filter(VISIT==20) %>% select(n_agents,BPclass2.factor) %>% table() %>% addmargins()
write.table(table1.2,row.names=T, col.names=NA, "visit20_BPstatus_n_agents_v2.csv")
table1.2

# # neat_tables
# 
# #remove missing values for each factor
# test.20<-test %>% filter(VISIT==20) %>% select(n_agents.factor,BPclass.factor) %>% na.omit()
# cont.table1<-contingency_table(list("n_agents"="n_agents.factor"),outcomes=list("BP class"="BPclass.factor"),crosstab_funcs=list(freq()),data=test.20)
# neat_table(cont.table1)
# 
# test.20<-test %>% filter(VISIT==20) %>% select(n_agents.factor,BPclass2.factor) %>% na.omit()
# cont.table1.2<-contingency_table(list("n_agents"="n_agents.factor"),outcomes=list("BP class2"="BPclass2.factor"),crosstab_funcs=list(freq()),data=test.20)
# neat_table(cont.table1.2)

# n_agents_visit.csv: n_agents(columns) by visit (rows)
# note: includes those not taking any antihypertensive meds
table2<-test %>% filter(VISIT%%10==0) %>% select(VISIT,n_agents) %>% table(useNA = "always") %>% addmargins()
write.table(table2, row.names=T, col.names=NA,"n_agents_visit.csv")
table2

table2.2<-test %>% filter(VISIT%%10==0, CASEID %in% cases) %>% select(VISIT,n_agents) %>% table(useNA = "always") %>% addmargins()
write.table(table2.2, row.names=T, col.names=NA,"n_agents_visit_cohort.csv")
table2.2

# BPstatus_visit.csv: BPclass (columns) by visit (rows, only even visits when ABPM obtained)
table3<-test %>% filter(VISIT%%20==0,CASEID %in% cases) %>% group_by(VISIT) %>% select(BPclass2.factor) %>% table() %>% addmargins(2)
write.table(table3,row.names=T, col.names = NA, "BPstatus_visit_cohort.csv")
table3
# round(100*prop.table(table3,margin=1),1) %>% addmargins(2)
# test.3<-test%>% filter(VISIT%%20==0)
# test.3$VISIT.factor<-as.factor(test.3$VISIT)
# test.3<-test.3 %>%  select(VISIT.factor,BPclass.factor,BPclass2.factor) %>% na.omit()
# cont.table3<-contingency_table(list("VISIT"="VISIT.factor"),outcomes=list("BP class"="BPclass.factor"),crosstab_funcs=list(freq()),data=test.3)
# neat_table(cont.table3)

# table3.2<-test %>% filter(VISIT%%20==0) %>% group_by(VISIT) %>% select(BPclass2.factor) %>% table()
# write.table(table3.2,row.names=T, col.names = NA, "BPstatus_visit_v2.csv")
# table3.2 %>% addmargins(2)
# round(100*prop.table(table3.2,margin=1),1) %>% addmargins(2)
# test.3$VISIT.factor<-as.factor(test.3$VISIT)
# cont.table3.2<-contingency_table(list("VISIT"="VISIT.factor"),outcomes=list("BP class2"="BPclass2.factor"),crosstab_funcs=list(freq()),data=test.3)
# neat_table(cont.table3.2)

# Counts of patients on each antihypertensive medication grouped by visit and sorted (descending)
unsorted<-medsum_full.3 %>% filter(VISIT%%10==0, CASEID %in% cases) %>% group_by(VISIT,med.corrected) %>% summarise(n_rx=n_distinct(CASEID)) %>% arrange(desc(n_rx),.by_group=T) %>% xtabs(formula=n_rx~addNA(factor(med.corrected))+VISIT)
table4<-unsorted[sort(unsorted[,1], decreasing=T,index.return=T)$ix,] %>% addmargins()
write.table(table4,row.names=T, col.names=NA,"visit_medcounts_cohort.csv")
table4

# overview of DDI
write_csv(medsum_full.old %>% filter(VISIT==20, CASEID %in% as.numeric(unlist(cases_high_DDI)),DDI>=2) %>% select(CASEID,age,AVWEIGHT,gfr,med.corrected,DLYDOSE,std_dose,max_dose, max_dosev2,renal_adjust,DDI,DDIv2),"high_DDI_v20.csv")
cases_high_DDI<-test.cohort %>% filter(VISIT==20, mean_DDI >=2) %>% select(CASEID)

test.20<-test %>% filter(VISIT==20)
test.20.DDI<-as_tibble(cbind(test.20$DDI_all[,1:4],test.20$n_agents,test.20$BPclass2.factor))
names(test.20.DDI)<-c(paste("DDI_",1:4,sep=''),"n_agents","BPclass2")
test.20.melted<-melt(test.20.DDI,id="n_agents")
test.20.melted.filtered<-test.20.melted %>% filter(value!=0)
test.20.melted.filtered$variable.factor<-as.factor(test.20.melted.filtered$variable)
test.20.melted.filtered$n_agents.factor<-as.factor(test.20.melted.filtered$n_agents)
ggplot(aes(y=value,fill=n_agents.factor,x=variable.factor),data=subset(test.20.melted.filtered, value<3))+labs(title="Drug Dose Index (DDI) vs n_agents",x="DDI for medications 1-4", y="DDI (dose/max dose)")+geom_boxplot()
ggplot(aes(y=value,fill=variable.factor,x=n_agents.factor),data=subset(test.20.melted.filtered, value<3))+labs(title="Drug Dose Index (DDI) vs n_agents\n(Visit 20)",x="# BP meds taken", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="DDI 1-4"))+geom_boxplot()

# DDI vs BP class
test.20<-test %>% filter(VISIT==20)
test.20.DDI<-as_tibble(cbind(test.20$DDI_all[,1:4],test.20$BPclass.factor))
names(test.20.DDI)<-c(paste("DDI_",1:4,sep=''),"BPclass")
test.20.melted.2<-melt(test.20.DDI,id="BPclass")
test.20.melted.2$BPclass.factor<-factor(test.20.melted.2$BPclass, labels=c("NL","WCH","PH","MH","AH","SAH"),ordered = TRUE)
test.20.melted.2.filtered<-test.20.melted.2 %>% filter(!is.na(BPclass.factor))
test.20.melted.2.filtered$variable.factor<-as.factor(test.20.melted.2.filtered$variable)
ggplot(aes(y=value,fill=BPclass.factor,x=variable.factor),data=subset(test.20.melted.2.filtered, value<1.5))+labs(title="Drug Dose Index (DDI) vs BP class",x="DDI for medications 1-4", y="DDI (dose/max dose)")+labs(title="Drug Dose Index (DDI) vs BP class\n(Visit 20)",x="DDI 1-4", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="BP class"))+geom_boxplot()

test.20<-test %>% filter(VISIT==20)
test.20.DDI<-as_tibble(cbind(test.20$DDI_all[,1:4],test.20$BPclass2.factor))
names(test.20.DDI)<-c(paste("DDI_",1:4,sep=''),"BPclass2")
test.20.melted.2<-melt(test.20.DDI,id="BPclass2")
test.20.melted.2$BPclass2.factor<-factor(test.20.melted.2$BPclass2, labels=c("NL","WCH","MH","AH"),ordered = TRUE)
test.20.melted.2.filtered<-test.20.melted.2 %>% filter(!is.na(BPclass2.factor))
test.20.melted.2.filtered$variable.factor<-as.factor(test.20.melted.2.filtered$variable)
ggplot(aes(y=value,fill=BPclass2.factor,x=variable.factor),data=subset(test.20.melted.2.filtered, value<1.5))+labs(title="Drug Dose Index (DDI) vs BP class",x="DDI for medications 1-4", y="DDI (dose/max dose)")+labs(title="Drug Dose Index (DDI) vs BP class\n(Visit 20)",x="DDI 1-4", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="BP class2"))+geom_boxplot()

# DDI vs nocturnal HTN
test.20<-test %>% filter(VISIT==20)
test.20.DDI<-as_tibble(cbind(test.20$DDI_all[,1:4],test.20$noctHTN))
names(test.20.DDI)<-c(paste("DDI_",1:4,sep=''),"noctHTN")
test.20.melted.6<-melt(test.20.DDI,id="noctHTN")
test.20.melted.6$noctHTN.factor<-as.factor(test.20.melted.6$noctHTN)
test.20.melted.6.filtered<-test.20.melted.6 %>% filter(!is.na(noctHTN.factor))
test.20.melted.6.filtered$variable.factor<-as.factor(test.20.melted.6.filtered$variable)
ggplot(aes(y=value,fill=noctHTN.factor,x=variable.factor),data=subset(test.20.melted.6.filtered, value<3))+labs(title="Drug Dose Index (DDI) vs nocturnal HTN status",x="DDI for medications 1-4", y="DDI (dose/max dose)")+labs(title="Drug Dose Index (DDI) vs nocturnal HTN status\n(Visit 20)",x="DDI 1-4", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="Nocturnal HTN (yes/no)"))+geom_boxplot()

# mean DDI vs BP class
test.20<-test %>% filter(VISIT==20,n_agents>=2)
test.20.DDI<-as_tibble(cbind(test.20$mean_DDI,test.20$BPclass.factor))
names(test.20.DDI)<-c("mean_DDI","BPclass")
test.20.melted.4<-melt(test.20.DDI,id="BPclass")
test.20.melted.4$BPclass.factor<-factor(test.20.melted.4$BPclass, labels=c("NL","WCH","PH","MH","AH","SAH"),ordered = TRUE)
test.20.melted.4.filtered<-test.20.melted.4 %>% filter(!is.na(BPclass.factor))
test.20.melted.4.filtered$variable.factor<-as.factor(test.20.melted.4.filtered$variable)
ggplot(aes(y=value,fill=BPclass.factor),data=subset(test.20.melted.4.filtered, value<3))+labs(title="Mean Drug Dose Index (DDI) vs BP class",x="BP class", y="mean DDI (dose/max dose)")+labs(title="Mean Drug Dose Index (DDI) vs BP class\n(Visit 20)",x="BP class", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="BP class"))+geom_boxplot()+scale_fill_brewer(palette="Spectral", direction=-1)

test.20<-test %>% filter(VISIT==20,n_agents>=2)
test.20.DDI<-as_tibble(cbind(test.20$mean_DDI,test.20$BPclass2.factor))
names(test.20.DDI)<-c("mean_DDI","BPclass2")
test.20.melted.4<-melt(test.20.DDI,id="BPclass2")
test.20.melted.4$BPclass2.factor<-factor(test.20.melted.4$BPclass2, labels=c("NL","WCH","MH","AH"),ordered = TRUE)
test.20.melted.4.filtered<-test.20.melted.4 %>% filter(!is.na(BPclass2.factor))
test.20.melted.4.filtered$variable.factor<-as.factor(test.20.melted.4.filtered$variable)
ggplot(aes(y=value,fill=BPclass2.factor),data=subset(test.20.melted.4.filtered, value<3))+labs(title="Mean Drug Dose Index (DDI) vs BP class",x="BP class", y="mean DDI (dose/max dose)")+labs(title="Mean Drug Dose Index (DDI) vs BP class\n(Visit 20)",x="BP class", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="BP class2"))+geom_boxplot()+scale_fill_brewer(palette="Spectral", direction=-1)

# mean_DDI vs CKD stage
test.20<-test %>% filter(VISIT==20, n_agents>=1)
test.20.DDI<-as_tibble(cbind(test.20$mean_DDI,test.20$CKD_stage))
names(test.20.DDI)<-c("mean_DDI","CKD_stage")
test.20.melted.5<-melt(test.20.DDI,id="CKD_stage")
test.20.melted.5$CKD_stage.factor<-factor(test.20.melted.5$CKD_stage, labels=c(1:5),ordered = TRUE)
test.20.melted.5.filtered<-test.20.melted.5 %>% filter(!is.na(CKD_stage.factor))
test.20.melted.5.filtered$variable.factor<-as.factor(test.20.melted.5.filtered$variable)
ggplot(aes(y=value,fill=CKD_stage.factor),data=subset(test.20.melted.5.filtered, value<3))+labs(title="Mean Drug Dose Index (DDI) vs CKD stage",x="CKD stage", y="mean DDI (dose/max dose)")+labs(title="Mean Drug Dose Index (DDI) vs CKD stage\n(Visit 20)",x="BP class", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="CKD stage"))+geom_boxplot()+scale_fill_brewer(palette="Spectral", direction=-1)

# DDI vs LVMIp
test.20<-test %>% filter(VISIT==20)
test.20.DDI<-as_tibble(cbind(test.20$DDI_all[,1:4],test.20$LVMIp))
names(test.20.DDI)<-c(paste("DDI_",1:4,sep=''),"LVMIp")
test.20.melted.3<-melt(test.20.DDI,id="LVMIp")
test.20.melted.3$LVMIp.factor<-factor(test.20.melted.3$LVMIp, labels=c(10,25,50,75,90,95),ordered = TRUE)
test.20.melted.3.filtered<-test.20.melted.3 %>% filter(!is.na(LVMIp.factor))
test.20.melted.3.filtered$variable.factor<-as.factor(test.20.melted.3.filtered$variable)
ggplot(aes(y=value,fill=LVMIp.factor,x=variable.factor),data=subset(test.20.melted.3.filtered, value<3))+labs(title="Drug Dose Index (DDI) vs LVMI percentile",x="DDI for medications 1-4", y="DDI (dose/max dose)")+labs(title="Drug Dose Index (DDI) vs LVMI percentile\n(Visit 20)",x="DDI 1-4", y="DDI (dose/max dose)")+guides(fill=guide_legend(title="LVMI percentile"))+geom_boxplot()

# Logistic Regression


# load visit 20 data
load(file="test.20.RData")
library(MASS)
m<-polr(test.20$BPclass2.factor~unlist(test.20$DDI_1)+unlist(test.20$DDI_2)+test.20$Upc+test.20$gfr+test.20$n_agents+test.20$BMIPCTAG+test.20$CKDONST+test.20$GHTOTINC+test.20$MALE1FE0+test.20$age,Hess=TRUE)
summary(m)

# relationship between mean_DDI and n_agents (ANA)
test.20.filtered<-test.20 %>% filter(mean_DDI<=3, !is.na(n_agents), n_agents>0)
group_by(test.20.filtered, n_agents) %>%
  summarise(
    count = n(),
    mean = mean(mean_DDI, na.rm = TRUE),
    sd = sd(mean_DDI, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(mean_DDI ~ n_agents, data = test.20.filtered)
# Summary of the analysis
summary(res.aov)
ggline(test.20.filtered, x = "n_agents", y = "mean_DDI", 
       add = c("mean_se", "jitter"), 
       order = c(1:4),
       ylab = "mean_DDI", xlab = "n_agents")

# looking at strata based on n_agents
# n_agents=1
test.20.filtered_n1<-test.20 %>% filter(mean_DDI<=3, !is.na(BPclass2.factor),n_agents==1)
group_by(test.20.filtered_n1, BPclass2.factor) %>%
  summarise(
    count = n(),
    mean = mean(mean_DDI, na.rm = TRUE),
    sd = sd(mean_DDI, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(mean_DDI ~ BPclass2.factor, data = test.20.filtered_n1)
# Summary of the analysis for n_agents=1
summary(res.aov)
ggline(test.20.filtered_n1, x = "BPclass2.factor", y = "mean_DDI", 
       add = c("mean_se", "jitter"), 
       order = c("NL","WCH","MH","AH"),
       ylab = "mean_DDI", xlab = "BPclass2.factor")

# within visit20, n_agents=1, compare mean_DDI for MH vs AH
test.20.filtered_n1.2<-test.20 %>% filter(mean_DDI<=3,n_agents==1,BPclass2.factor=="MH"|BPclass2.factor=="AH")
group_by(test.20.filtered_n1.2, BPclass2.factor) %>%
  summarise(
    count = n(),
    mean = mean(mean_DDI, na.rm = TRUE),
    sd = sd(mean_DDI, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(mean_DDI ~ BPclass2.factor, data = test.20.filtered_n1.2)
# Summary of the analysis for n_agents=1
summary(res.aov)
ggline(test.20.filtered_n1.2, x = "BPclass2.factor", y = "mean_DDI", 
       add = c("mean_se", "jitter"), 
       order = c("MH","AH"),
       ylab = "mean_DDI", xlab = "BPclass2.factor")

# ----------
# define cohort (n=497) to be used for all subsequent analysis
# inclusion criteria based on VISIT 20:
# successfull ABPM
# available data for: SCR, AVWEIGHT,AVHEIGHT

cohort<-test %>% filter(VISIT==20,ABPMSUCCESS==1, !is.na(SCR),!is.na(AVHEIGHT),!is.na(AVWEIGHT)) %>% select(CASEID) %>% unlist %>% as.numeric()
library(RColorBrewer)
# histogram by medication class and CKD stage
# rename BPmedgroup for those not taking BP meds to "none"
BPmedgroups<-levels(addNA(medsum_full.old$BPmedgroup))
BPmedgroups[14]<-"none"
#medsum_full.old$BPmedgroup.2<-factor(medsum_full.old$BPmedgroup, levels=BPmedgroups)
#medsum_full.old$BPmedgroup.2[which(is.na(medsum_full.old$BPmedgroup.2))]<-"none"
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases), aes(CKD_stage))
colourCount = 14
getPalette = colorRampPalette(brewer.pal(14, "Set2"))
g+
  #geom_bar(aes(fill=BPmedgroup.2),width=0.5)+
  scale_fill_manual(values = colorRampPalette(brewer.pal(14, "Set1"))(colourCount))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Histogram on categorical variable", subtitle="BP medication class by CKD stage (VISIT 20)")+
  guides(fill=guide_legend(title="BP med class"))

# ordered bar plot of medications by frequency at visit 20
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases,!is.na(med.corrected)), aes(x=reorder(med.corrected,med.corrected,function(x)-length(x))))
g+
  geom_bar(aes(fill=CKD_stage),width=0.5)+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(5, "Set1")))(5))+
  labs(title="Frequency of medications at visit 20", caption="source:CKiD",x="Medication name")+
  theme(axis.text.x=element_text(angle=65,vjust=0.6))

# number of agents and BP med class at V20
colourCount=14
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases,!is.na(n_agents)), aes(n_agents))
g+
  #geom_bar(aes(fill=reorder(BPmedgroup.2,BPmedgroup.2,function(x) - length(x))),width=0.5)+
  scale_fill_manual(values = colorRampPalette(brewer.pal(14, "Set1"))(colourCount))+
  labs(title="BP medication class by number of agents",subtitle="(VISIT 20)", caption="source:CKiD",x="number of BP agents")+
  guides(fill=guide_legend(title="BP medication class"))

# number of agents and BP status at V20
colourCount = 4
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases,!is.na(n_agents), !is.na(BPclass2.factor)), aes(n_agents))
g+
  geom_bar(aes(fill=BPclass2.factor),width=0.5)+
  scale_fill_manual(labels=c("Normotensive","White-coat HTN","Masked HTN","Ambulatory HTN"),values = colorRampPalette(rev(brewer.pal(4, "Spectral")))(colourCount))+
  labs(title="BP status by number of agents",subtitle="(VISIT 20)", caption="source:CKiD",x="number of BP agents")+
  guides(fill=guide_legend(title="BP status"))

# CKD stage and BP status at V20      
colourCount = 4
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases, !is.na(BPclass2.factor)), aes(CKD_stage))
g+
  geom_bar(aes(fill=BPclass2.factor),width=0.5)+
  scale_fill_manual(labels=c("Normotensive","White-coat HTN","Masked HTN","Ambulatory HTN"),values = colorRampPalette(rev(brewer.pal(4, "Spectral")))(colourCount))+
  labs(title="BP status by CKD stage",subtitle="(VISIT 20)", caption="source:CKiD",x="CKD stage")+
  guides(fill=guide_legend(title="BP status"))

# BP status and proteinuria at V20      
colourCount = 4
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases, !is.na(BPclass2.factor)), aes(Upc.factor))
g+
  geom_bar(aes(fill=BPclass2.factor),width=0.5)+
  scale_fill_manual(labels=c("Normotensive","White-coat HTN","Masked HTN","Ambulatory HTN"),values = colorRampPalette(rev(brewer.pal(4, "Spectral")))(colourCount))+
  labs(title="BP status by severity of proteinuria (VISIT 20)",subtitle="(Up:c ratio, mg/mg cr) \n[0-0.5=normal,0.5-1=mild,1-2=moderate,>2=severe]", caption="source:CKiD",x="Degree of Proteinuria")+
  guides(fill=guide_legend(title="BP status"))

# BP status and LVMIp at V20
colourCount = 6
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases, !is.na(BPclass2.factor),!is.na(LVMIp)), aes(BPclass2.factor))
g+
  geom_bar(aes(fill=LVMIp.factor),width=0.5)+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(6, "Spectral")))(colourCount))+
  labs(title="LVMI percentile by BP status (VISIT 20)",subtitle="LVMI %-ile (Khouri et al, 2009)", caption="sources: Khoury PR, Mitsnefes M, Daniels SR, Kimball TR. \n Age-Specific Reference Intervals for Indexed Left Ventricular Mass in Children. \n Journal of the American Society of Echocardiography. 2009;22(6):709-714 \n and CKiD",x="BP status")+
  guides(fill=guide_legend(title="LVMI percentile"))

# n_agents and LVMIp at V20
colourCount = 6
g<-ggplot(medsum_full.old %>% filter(VISIT==20, CASEID %in% cases, !is.na(n_agents),!is.na(LVMIp)), aes(n_agents))
g+
  geom_bar(aes(fill=LVMIp.factor),width=0.5)+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(6, "Spectral")))(colourCount))+
  labs(title="LVMI percentile by number of BP agents (VISIT 20)",subtitle="LVMI %-ile (Khouri et al, 2009)", caption="sources: Khoury PR, Mitsnefes M, Daniels SR, Kimball TR. \n Age-Specific Reference Intervals for Indexed Left Ventricular Mass in Children. \n Journal of the American Society of Echocardiography. 2009;22(6):709-714 \n and CKiD",x="# BP medications")+
  guides(fill=guide_legend(title="LVMI percentile"))

# timing of ABPM and visit number
g<-ggplot(test %>% filter(CASEID %in% cases),aes(ABPM_DATE))
g+geom_density(aes(fill=VISIT.factor),alpha=0.8)+
  labs(title="ABPM date grouped by visit number", caption="source: CKiD",x="ABPM date (yrs from entry into CKiD)",fill="VISIT")

# -----------
# correlation between n_agents and BP status
dt <- table(test.cohort$n_agents,test.cohort$BPclass2.factor, exclude = c(-1,NA))
chisq <- chisq.test(dt)
chisq
corrplot(chisq$residuals, is.cor = FALSE, title="chi-squared test residuals for\n BP class vs n_agents\n (X-squared = 27.946, df = 12, p-value = 0.005634)",addCoef.col = 1,mar=c(1,0,2.5,1),number.cex=0.6)

# standard dose by medication
temp<-medsum_full.old %>% filter(CASEID %in% cases, VISIT==20)
g<-ggplot(temp, aes(std_dose))
colourCount=length(unique(temp$med.corrected))
g+ geom_density(aes(fill=factor(med.corrected)),alpha=0.5)+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(colourCount, "Spectral")))(colourCount))+
  labs(title="Standardized dose (mg/kg) vs BP medication",
       subtitle="Std dose grouped by medication",
       caption="source:CKiD",
       x="Standardized dose (mg/kg)",
       fill="BP medication")+
  coord_cartesian(xlim = c(0,1.2))

# correlation between factors:
library(ggcorplot)
temp<-test.cohort %>% filter(VISIT==20) %>% select(n_agents,Upc,gfr,BMIPCTAG,age,MALE1FE0,BPclass2,SBPZAGH2017,DBPZAGH2017,LVMIp,CKDONST, GHTOTINC)
corr<-round(cor(temp,use="complete.obs"),2)
ggcorrplot(corr,hc.order=TRUE,
           type="lower",
           lab=TRUE,
           lab_size=3,
           method="circle",
           colors=c("tomato2","white","spring green3"),
           title="Correlogram of Visit 20 cohort",
           ggtheme=theme_bw)

# distribution of mean_DDI for visit 20
temp<-test.cohort %>% filter(VISIT==20,!is.na(mean_DDI)) %>% select(mean_DDI) %>% unlist %>% as.numeric()
hist(temp,axes=FALSE,breaks=150,xlim=c(0.1,2),main="Histogram of mean DDI (VISIT 20)", xlab="mean DDI")
lines(density(temp,adjust=0.5,weights=rep(18/length(temp),length(temp))))
axis(1, lwd = 2, xaxp = c(-0.0, 2, 20))
axis(2,lwd=2,xaxp=c(0,30,5))

# distribution of mean_DDI by n_agents (V20)
temp<-test.cohort %>% filter(VISIT==20,!is.na(mean_DDI)) %>% select(mean_DDI,n_agents)
temp$n_agents.factor<-factor(temp$n_agents)
g<-ggplot(temp,aes(mean_DDI))
g+geom_density(aes(fill=n_agents.factor),alpha=0.8)+
  labs(title="Distribution of mean_DDI grouped by n_agents",
       subtitle="(VISIT 20)",
       caption="Source:CKiD",
       x="mean DDI",
       fill="# BP meds")

g<-ggplot(temp,aes(mean_DDI))+scale_fill_brewer(palette="Spectral")
g+geom_histogram(aes(fill=n_agents.factor),
                 binwidth=0.1,
                 col="black",
                 size=0.1)+
  labs(title="Distribution of mean_DDI by n_agents",
       subtitle="(VISIT 20)")

# distribution of mean_DDI by BPstatus (V20)
temp<-test.cohort %>% filter(VISIT==20,!is.na(mean_DDI),!is.na(BPclass2.factor)) %>% select(mean_DDI,BPclass2.factor)
g<-ggplot(temp,aes(mean_DDI))
g+geom_density(aes(fill=BPclass2.factor), adjust=0.25,alpha=0.8)+
  labs(title="Distribution of mean_DDI grouped by BP status",
       subtitle="(VISIT 20)",
       caption="Source:CKiD",
       x="mean DDI",
       fill="BP status")+
  coord_cartesian(xlim = c(0,1.2))

# distribution of mean_DDI by GNGDIAG (V20)
temp<-test.cohort %>% filter(VISIT==20,!is.na(mean_DDI),!is.na(BPclass2.factor)) %>% select(mean_DDI,BPclass2.factor)
g<-ggplot(temp,aes(mean_DDI))
g+geom_density(aes(fill=BPclass2.factor), adjust=0.25,alpha=0.8)+
  labs(title="Distribution of mean_DDI grouped by BP status",
       subtitle="(VISIT 20)",
       caption="Source:CKiD",
       x="mean DDI",
       fill="BP status")+
  coord_cartesian(xlim = c(0,1.2))

write_csv(medsum_full.old %>% filter(VISIT==20, CASEID %in% as.numeric(unlist(cases_high_DDI)),DDI>=2) %>% select(CASEID,age,AVWEIGHT,gfr,med.corrected,DLYDOSE,std_dose,max_dose, max_dosev2,renal_adjust,DDI,DDIv2),"high_DDI_v20.csv")


## Linear Regression modeling for DDI
# set up: check DDI for normality
par(mfrow=c(2,1))
qqnorm(test.cohort$mean_DDI, pch = 1, frame = FALSE,main="mean DDI")
qqline(test.cohort$mean_DDI, col = "steelblue", lwd = 2)
qqnorm(log(test.cohort$mean_DDI), pch = 1, frame = FALSE,main="transformed: ln(mean DDI)")
qqline(log(test.cohort$mean_DDI), col = "steelblue", lwd = 2)

ggplot(test.cohort,aes(x=mean_DDI,y=SBPZAGH2017))+
  geom_point(aes(col=CKD_stage,size=n_agents))+
  geom_smooth(method="loess",se=F)+
  xlim(c(0,3))+
  ylim(c(-3,3))+
  geom_encircle(aes(x=mean_DDI,y=SBPZAGH2017),
                data=test.cohort %>% filter(BPclass2.factor=="MH"),
                color="red",
                size=2,
                expand=0.08)+
  labs(subtitle="test",
       y="SBPz",
       x="mean_DDI",
       title="scatterplot+encircle",
       caption="Source:CKiD")

# timing of ABPM, echo, PE, meds
# gather all times into test
test<-full_join(cardio %>% select(CASEID,VISIT,ABPM_DATE),test)
test<-full_join(pe %>% select(CASEID,VISIT,PEVISDAT),test)
test<-full_join(echo %>% select(CASEID,VISIT,ECHODATEY),test)
test<-full_join(medsum_full %>% select(CASEID,VISIT,MSVISDAT),test)

# filter only cohort patients
temp<-test %>% select(-c("DDI_all")) %>%filter(CASEID %in% cases)
# melt dates into new dataframe to display as histograms
temp<-melt(temp,id=c("CASEID","VISIT"),measure.vars=c("ECHODATEY","PEVISDAT","MSVISDAT","ABPM_DATE"))
g<-ggplot(temp %>% filter(VISIT==20),aes(value))
g+geom_density(aes(fill=factor(variable)),adjust=0.5,alpha=0.5)+
  labs(title="Timing of ABPM, EHCO, PE, MEDS data",
       subtitle="(VISIT 20)",
       caption="Source:CKiD",
       x="Date (yrs from cohort entry)",
       fill="Data type")+
  coord_cartesian(xlim = c(0.5,1.5))

# show all visits
g<-ggplot(temp,aes(value))
g+geom_density(aes(fill=factor(variable)),adjust=0.5,alpha=0.5)+
  labs(title="Timing of ABPM, EHCO, PE, MEDS data",
       subtitle="(all visits)",
       caption="Source:CKiD",
       x="Date (yrs from cohort entry)",
       fill="Data type")


# mean_DDI analysis (boxplots by different variables)
library(ggthemes)

# by n_agents
g<-ggplot(test.cohort %>% filter(between(n_agents.factor,2,5)),aes(n_agents.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(n_agents.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by number of BP agents",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="# BP meds",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,3))+
  guides(fill=guide_legend(title="# BP meds"))

# by CKD_stage
test.cohort$CKD_stage.factor<-factor(test.cohort$CKD_stage)
g<-ggplot(test.cohort ,aes(CKD_stage.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(CKD_stage.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by CKD stage",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="CKD stage",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  guides(fill=guide_legend(title="CKD stage"))

# by Upc group
g<-ggplot(test.cohort ,aes(Upc.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(Upc.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by degree of proteinuria",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="degree of proteinuria",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  guides(fill=guide_legend(title="Degree of proteinuria"))

# by BP status
g<-ggplot(test.cohort ,aes(BPclass2.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(BPclass2.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by BP status",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="BP status",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  guides(fill=guide_legend(title="BP status"))

# by primary CKD diagnosis
g<-ggplot(test.cohort ,aes(GNGDIAG.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(GNGDIAG.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by CKD diagnosis (G/NG)",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="primary CKD dx",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  guides(fill=guide_legend(title="Primary CKD Diagnosis"))

# by HTN status
test.cohort$HTN.factor<-factor(test.cohort$HTN)
g<-ggplot(test.cohort ,aes(HTN.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(HTN.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by BP status (HTN - yes/no)",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="HTN (yes/no)",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  guides(fill=guide_legend(title="HTN status (yes/no)"))

colourCount=4
g<-ggplot(test.cohort ,aes(BPclass2.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(BPclass2.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by BP status",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="BP status",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(4, "Spectral")))(4))+
  guides(fill=guide_legend(title="BP status"))

colourCount=4
g<-ggplot(test.cohort ,aes(Upc.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(Upc.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by Degree of proteinuria",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="Degree of proteinuria",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(4, "Spectral")))(4))+
  guides(fill=guide_legend(title="Degree of proteinuria"))

colourCount=6
test.cohort<-test.cohort %>% filter(!is.na(LMVIp))
test.cohort$LVMIp.factor<-factor(test.cohort$LVMIp)
g<-ggplot(test.cohort ,aes(LVMIp.factor,mean_DDI))
g+geom_boxplot(aes(fill=factor(LVMIp.factor)))+
  theme(axis.text.x=element_text(vjust=0.6))+
  labs(title="Boxplot of mean_DDI by LVMI%ile",
       subtitle="(VISIT 20)",
       caption="source: CKiD",
       x="LVMI %ile",
       y="mean_DDI")+
  coord_cartesian(ylim = c(0,2.8))+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(6, "Spectral")))(6))+
  guides(fill=guide_legend(title="LVMI %ile"))


# univariate analysis of ln(mean_DDI)
library(ggpubr)
# Compare means of mean_DDI by BP status
compare_means(mean_DDI.ln ~ BPclass2.factor,  data = test.cohort %>% filter(n_agents>=1),method="anova")
compare_means(mean_DDI.ln ~ BPclass2.factor,  data = test.cohort %>% filter(n_agents>=1),method="t.test")

# Compare means of mean_DDI by n_agents 
compare_means(mean_DDI.ln ~ n_agents,  data = test.cohort %>% filter(n_agents>=1),method="anova")
compare_means(mean_DDI.ln ~ n_agents,  data = test.cohort %>% filter(n_agents>=1),method="t.test")

my_comparisons <- list( c("1", "2") )
ggboxplot(test.cohort %>% filter(n_agents>=1), x = "n_agents.factor", y = "mean_DDI.ln",
          color = "n_agents.factor", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  coord_cartesian(ylim = c(-2,3))

compare_means(SBPZAGH2017 ~ n_agents,  data = test.cohort,method="anova")
compare_means(SBPZAGH2017 ~ n_agents,  data = test.cohort,method="t.test")

my_comparisons <- list( c("0","3"),c("1", "3"),c("2","3") )
ggboxplot(test.cohort %>% filter(!is.na(n_agents.factor)), x = "n_agents.factor", y = "SBPZAGH2017",
          color = "n_agents.factor", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  coord_cartesian(ylim = c(-3,4))

##########
# cardio data organization
# overview of data...
# display the time differences between ABPM and casual BP measurement dates
temp<-cardio %>% filter(VISIT %%20==0,CASEID %in% cases) %>% mutate(date_dif=ABPM_DATE-DB_DATE) %>% group_by(VISIT) %>% summarise_at(vars(date_dif),funs(mean_dif_days=365*mean(.,na.rm=TRUE),max_dif_yr=max(.,na.rm=TRUE), sd_dif_days=365*sd(.,na.rm=TRUE)))
write_csv(temp,"ABPM_PE_time_dif.csv")
# display the number of observations where the difference between the date of visit and the ABPM date is > 90 days (=TRUE)
# consider excluding these observations since the BP status is less reliable
temp<-cardio %>% filter(VISIT %% 20 ==0,CASEID %in% cases) %>% mutate(date_dif=abs(ABPM_DATE-DB_DATE)) %>% group_by(VISIT) %>% mutate(big_dif=date_dif>90/365) %>% select(big_dif) %>% table %>% addmargins
write.table(temp,"ABPM_PE_bigtime_dif.csv")
# display the counts for all the BP statuses grouped by visit
# patients with BP status -1 occured when casual BP was not classified (either no BP entered at that visit, or BP entered, but no percentile/z-score)
# cardio %>% group_by(VISIT) %>% filter(ABPMSUCCESS==1, BPstatus !=-1) %>% select(BPstatus) %>% table() %>% addmargins
# display information about timing of ABPM data, grouped by visits
temp<-cardio%>% filter(VISIT %% 20==0,ABPMSUCCESS==1, CASEID %in% cases) %>% group_by(VISIT) %>% summarise_at(vars(ABPM_DATE), funs(date_mean=mean,date_min=min,date_max=max,n_pts=length))
write_csv(temp,"ABPM_timing_by_visit.csv")

# stacked barplot of BP status (in proportions) grouped by visit
test %>% filter(VISIT%%10==0) %>% group_by(VISIT) %>% select(BPstatus2017) %>% table(useNA = 'no') %>%  plot()

# draw histograms of DDI for each drug
temp<-paste("if (!all(is.na(",drugname[1:43],"$DDI))) hist(",drugname[1:43],"$DDI, xlim=c(0,3), breaks=200, main=paste(\"",drugname[1:43],"\"))", sep="")
par(mfrow=c(6,4))
par(mar=c(2,2,1,1))
eval(parse(text=temp))

# comparing patients based on DDI
# create tibble containing DDI for all meds, sorted in order of descending DDI
# create separate columns for each medication (med1-med5)
# create variable mean_DDI representing mean DDI for each patient (ie, mean DDI of all meds for each patient, patients with 0 meds have mean DDI = 0)
#temp_DDI<-test%>% select(ends_with('_DDI')) %>% as.matrix %>% apply(.,1,sort,decreasing=T,na.last=T)
test$mean_DDI<-NULL
test$DDI_all<-as.data.frame(test %>% select(ends_with("_DDI")))
DDI_all.2<-apply(test$DDI_all,1,FUN=function(x) sort(x,decreasing=TRUE,na.last=TRUE))
test$DDI_all<-t(DDI_all.2)[,1:5]
#test$DDI_all[which(test$n_agents==0),1]<-0
temp<-paste("test$DDI_",1:5,"<-test$DDI_all[,",1:5,"]",sep='')
eval(parse(text=temp))
test$mean_DDI<-NA
test$mean_DDI<-apply(test$DDI_all,1,FUN=function(x) mean(x,na.rm=TRUE))

# compare LVMIp to LVHE: show table of all LVMIp percentiles and the average proportion of patients categorized as LVH (using LVHE)
test %>% select(CASEID,VISIT,LVMI,LVHF,LVHE, LVMIp) %>% group_by(LVMIp) %>% summarise(mean_LVHE=mean(LVHE, na.rm=TRUE), mean_LVHF=mean(LVHF,na.rm=TRUE))
test %>% select(CASEID,VISIT,LVMI,LVHF,LVHE, LVMIp) %>% group_by(LVHE,LVHF) %>% summarise(mean_LVMIp=mean(LVMIp, na.rm=TRUE))

library(ggplot2)
library(ggcorrplot)
library(dplyr)
test.cohort.2<-test.cohort %>% select(BPclass,Ualb_pct,GNGDIAG,BMIz,gfr,age,n_agents,LVMIp,GHTOTINC,MALE1FE0,Upc,sum_DDI)
corr<-round(cor(test.cohort.2),1)
ggcorrplot(corr,hc.order=FALSE,
           type="lower",
           lab=TRUE,
           lab_size=3,
           method="circle",
           colors=c("dark red","white","dark green"),
           title="Correlogram of study population characteristics",
           ggtheme=theme_bw)

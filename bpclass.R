# BP status classification using:
# Flynn JT, Daniels SR, Hayman LL, et al. Update: Ambulatory Blood Pressure Monitoring in Children and Adolescents. Hypertension. 2014;63(5):1116-1135.
# classification:

# Code	BP class
# 0	Normal
# 1	White coat hypertension (WCH)
# 2	Prehypertension (PH)
# 3	Masked hypertension (MH)
# 4	Ambulatory hypertension (AH)
# 5	Severe ambulatory hypertension (SAH)
# 2.1	Prehypertension variant 1
# 3.1	*Masked hypertension variant 1
# 3.2	Masked hypertension variant 2
# 3.3	Masked hypertension variant 3
# 3.4	Masked hypertension variant 4
# 3.5	Masked hypertension variant 5
# 4.1	*Ambulatory hypertension variant 1
# 4.2	Ambulatory hypertension variant 2
# 4.3 Ambulatory hypertension variant 3
# (* = defined as unclassified by UpToDate/AHA 2014)

# parameters needed to determine status include:
# index (wake and sleep, systolic and diastolic) = 4 values
# load (wake and sleep, systolic and diastolic) = 4 values
# BP and BP percentiles (systolic and diastolic) = 4 values

# wake and sleep systolic and diastolic means
# note: mean indices and loads are calculated using wake and sleep values of each
bpclass<-function(wksi,wkdi,slsi,sldi,wksl,wkdl,slsl,sldl,sbp,dbp,sbpp,dbpp) {
  # if all indices, loads, means and percentiles are missing, return NA, otherwise use the maximum index, load, mean and percentile, respectively
  maxindex<-ifelse(all(is.na(c(wksi,wkdi,slsi,sldi))),NA,max(wksi,wkdi,slsi,sldi, na.rm=TRUE))
  maxload<-ifelse(all(is.na(c(wksl,wkdl,slsl,sldl))),NA,max(wksl,wkdl,slsl,sldl, na.rm=TRUE))
  # maxsysmean<-ifelse(all(is.na(c(wksysmn,slsysmn))),NA,max(wksysmn,slsysmn,na.rm=TRUE))
  # maxdiamean<-ifelse(all(is.na(c(wkdiamn,sldiamn))),NA,max(wkdiamn,sldiamn,na.rm=TRUE))
  maxbpp<-ifelse(all(is.na(c(sbpp,dbpp))),NA, max(sbpp,dbpp, na.rm=TRUE))
  # if all maximums are missing, return NA
  if (anyNA(c(maxindex,maxload,maxbpp))) return (NA)
  if (all(maxbpp<90,maxindex<1,maxload<25)) return(0) #NL
  if (all(maxbpp<90,maxload>=25,maxindex>=1)) return (3) #MH
  if (all((maxbpp>=90|sbp>=120|dbp>=80) & maxbpp<95,maxload>=25 & maxload<50,maxindex<1)) return(2) #PH
  if (all(maxbpp>=95,maxload<25,maxindex<1)) return(1) #WCH
  if (all(maxbpp>=95,maxload>=25 & maxload<50,maxindex>=1)) return(4) #AH
  if (all(maxbpp>=95,maxload>=50,maxindex>=1)) return(5) #SAH
  if (all(maxbpp<90,maxindex<1,maxload>=25 & maxload<50)) return(3.1) # MH1
  if (all(maxbpp<90,maxindex>=1,maxload<25)) return(3.2) #MH2
  if (all((maxbpp>=90|sbp>=120|dbp>=80) & maxbpp<95,maxindex<1,maxload<25)) return(2.1) #PH1
  if (all((maxbpp>=90|sbp>=120|dbp>=80) & maxbpp<95,maxindex>=1,maxload<25)) return(2.2) #PH2
  if (all((maxbpp>=90|sbp>=120|dbp>=80) & maxbpp<95,maxindex>=1,maxload>=25 & maxload<50)) return(3.3) #MH3
  if (all((maxbpp>=90|sbp>=120|dbp>=80) & maxbpp<95,maxindex>=1,maxload>=50)) return(3.4) #MH4
  if (all(maxbpp>=95,maxindex>=1,maxload<25)) return(4.2) #AH2
  if (all(maxbpp>=95,maxindex<1,maxload>=25 & maxload<50)) return(4.1) #AH1
  if (all(maxbpp>=95,maxindex<1,maxload>=50)) return(4.3) #AH3
  if (all((maxbpp>=90|sbp>=120|dbp>=80) & maxbpp<95,maxindex<1,maxload>=50)) return(2.3) #PH3
  if (all(maxbpp<90,maxindex<1,maxload>=50)) return(3.5) #MH5

  
  # if (all(maxbpp>=95,maxindex>=1 , maxload>=50)) return(5) #SAH
  # if (all(maxbpp>=95,maxindex<1,maxload<25)) return(1) # WCH
  # if (all(maxbpp<90,maxindex>=1, between(maxload,25,50))) return(3) # MH
  # if (all(maxbpp>=95,maxindex>=1 , between(maxload,25,50))) return(4) # AH
  # if (all((maxbpp>=90 & maxbpp<95)| sbp>=120 | dbp>=80,maxindex<1,maxload>=25)) return(2) #pre-HTN
  # if (all(maxload>=25, maxindex<1, (maxbpp>=95 | maxbpp<90))) return(3.5) # unclassified (no consensus on management)
  # if (all((maxbpp>=90 & maxbpp<95),maxload<25,maxindex<1)) return(0.5) # unclassified (my own definition - less severe than WCH, but milder than pre-HTN)
  # if (all(maxbpp<90,maxload>=50,maxindex>=1)) return(3.5) # unclassified (my own definition - severe masked HTN - similar to MH, but load >50%)
  # if (all(maxbpp>=90,maxindex>=1,maxload>=25))
  
  return (-1)
}

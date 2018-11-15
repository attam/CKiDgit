# BP status classification using:
# Flynn JT, Daniels SR, Hayman LL, et al. Update: Ambulatory Blood Pressure Monitoring in Children and Adolescents. Hypertension. 2014;63(5):1116-1135.
# classification:
# normal = 0
# white-coat HTN = 1
# pre-HTN = 2
# masked HTN = 3
# ambulatory HTN = 4
# severe ambulatory HTN = 5
# unclassified (elevated BP load with normal ambulatory mean and normal or elevated clinic BP) = 3.5
# parameters needed to determine status include:
# index (wake and sleep, systolic and diastolic) = 4 values
# load (wake and sleep, systolic and diastolic) = 4 values
# BP and BP percentiles (systolic and diastolic) = 4 values
# working group systolic and diastolic 90% limit
# wake and sleep systolic and diastolic means
# note: mean indices and loads are calculated using wake and sleep values of each
bpclass<-function(wksi,wkdi,slsi,sldi,wksl,wkdl,slsl,sldl,sbp,dbp,sbpp,dbpp,s90,d90,wksysmn,slsysmn,wkdiamn,sldiamn) {
  maxindex<-ifelse(all(is.na(c(wksi,wkdi,slsi,sldi))),NA,max(wksi,wkdi,slsi,sldi, na.rm=TRUE))
  maxload<-ifelse(all(is.na(c(wksl,wkdl,slsl,sldl))),NA,max(wksl,wkdl,slsl,sldl, na.rm=TRUE))
  maxsysmean<-ifelse(all(is.na(c(wksysmn,slsysmn))),NA,max(wksysmn,slsysmn,na.rm=TRUE))
  maxdiamean<-ifelse(all(is.na(c(wkdiamn,sldiamn))),NA,max(wkdiamn,sldiamn,na.rm=TRUE))
  maxbpp<-ifelse(all(is.na(c(sbpp,dbpp))),NA, max(sbpp,dbpp, na.rm=TRUE))
  if (all(is.na(c(maxindex,maxload,maxbpp)))) return (NA)
  if (all(maxbpp<90,maxsysmean<s90,maxdiamean<d90,maxload<25,na.rm=TRUE)) return(0) #NL
  if (all(maxbpp>95,maxindex>1 , maxload>50, na.rm=TRUE)) return(5) #SAH
  if (all(maxbpp>95,maxindex<1,maxload<25, na.rm=TRUE)) return(1) # WCH
  if (all(maxbpp<95,maxindex>1 | between(maxload,25,50),na.rm=TRUE)) return(3) # MH
  if (all(maxbpp>95,maxindex>1 , between(maxload,25,50),na.rm=TRUE)) return(4) # AH
  if (all((maxbpp>=90 | sbp>120 | dbp>80),maxindex<1,maxload>=25, na.rm=TRUE)) return(2) #pre-HTN
  if (all(maxload>25, maxsysmean<s90,maxdiamean<d90, (maxbpp>=95 | maxbpp<90),na.rm=TRUE)) return(3.5) # unclassified (no consensus on management)
  return (NA)
}
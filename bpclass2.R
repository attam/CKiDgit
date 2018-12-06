# BP status classification:
# Samuels J, Ng D, Flynn JT, et al. AMBULATORY BLOOD PRESSURE PATTERNS IN CHILDREN WITH CHRONIC KIDNEY DISEASE. Hypertension. 2012;60(1):43-50. doi:10.1161/HYPERTENSIONAHA.111.189266
# classification:

# Code	BP class
# 0	Normal
# 1	White coat hypertension (WCH)
# 2	Prehypertension (PH)
# 3	Masked hypertension (MH)
# 4	Ambulatory hypertension (AH)

# parameters needed to determine status include:
# index (wake and sleep, systolic and diastolic) = 4 values
# load (wake and sleep, systolic and diastolic) = 4 values
# BP and BP percentiles (systolic and diastolic) = 4 values

# wake and sleep systolic and diastolic means
# note: mean indices and loads are calculated using wake and sleep values of each
bpclass2<-function(wksi,wkdi,slsi,sldi,wksl,wkdl,slsl,sldl,sbp,dbp,sbpp,dbpp) {
  # if all indices, loads, means and percentiles are missing, return NA, otherwise use the maximum index, load, mean and percentile, respectively
  maxindex<-ifelse(all(is.na(c(wksi,wkdi,slsi,sldi))),NA,max(wksi,wkdi,slsi,sldi, na.rm=TRUE))
  maxload<-ifelse(all(is.na(c(wksl,wkdl,slsl,sldl))),NA,max(wksl,wkdl,slsl,sldl, na.rm=TRUE))
  maxbpp<-ifelse(all(is.na(c(sbpp,dbpp))),NA, max(sbpp,dbpp, na.rm=TRUE))
  # if all maximums are missing, return NA
  if (anyNA(c(maxindex,maxload,maxbpp))) return (NA)
  if (all(maxbpp<95,maxindex<1,maxload<25)) return(0) #NL
  if (all(maxbpp<95,(maxload>=25|maxindex>=1))) return (2) #MH
  if (all(maxbpp>=95,maxload<25,maxindex<1)) return(1) #WCH
  if (all(maxbpp>=95,(maxload>=25 | maxindex>=1))) return(3) #AH
  # test to check if any uncategorized patients
  return (-1)
}

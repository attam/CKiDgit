# function to extract the drug info from data structured in the CKiDgit project
# inputs: 
# - x=data.frame
# - rx=drug name (not case sensitive, but needs to be correctly spelled, and using generic drug names)
# - params= columns wanted (eg, 1:4 or c(1,3,4)) [columns: BPmedgroup,DLYFREQ,DLYDOSE,std_dose]

# function to extract rx_info given x (a data.frame), rx (either name of the drug, eg "AMLODIPINE", or index corresponding to the column containing the drug info) and param (1 = BPmedgroup, 2=DLYFREQ, 3=DLYDOSE, 4=std_dose ie, per kg)
get_drug_info <-function(x,rx,params) {
  if (typeof(rx)=="character") rx<-which(names(x)==toupper(rx))
  nrows<-dim(x)[1]
  out<-matrix(unlist(x[,rx]),nrow=nrows,ncol=length(params),byrow=TRUE,dimnames=(list(NULL,c("BPmedgroup","DLYFREQ","DLYDOSE","std_dose","DDI"))))
  return (as_tibble(out) %>% select(params))
  # return(unlist(lapply(x[,rx],'[[',param)))}

}

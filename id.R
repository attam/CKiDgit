# function to identify CASEID and VISIT corresponding to given row number of x
# if inv parameter is set (=TRUE), will assume n is 2-element vector containing CASEID and VISIT, and returns corresponding row of x
id<-function(x,y,inverse=NULL){
  if (!is.null(inverse)) {
    return(which(x$CASEID==y[1] & x$VISIT==y[2]))    
  } else {
    return(c(x$CASEID[y],x$VISIT[y]))
  }
}
#' A function to get the location change between samples
#'
#' This function determines 
#' @title locationChange
#' @aliases locationChange
#' @keywords location change RNA
#' @usage locationChange(sample, point=rep(0,nrow(sample)), base=sample[1,], margin=1, tol=0.1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param point An optional numeric vector containing the location of the disruption (e.g. the mutation in an RNA)
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @export
#' @details This function calculates the distance of change from the disruption using the change in pattern to determine the location of change.
#' @return A numeric vector of location changes.
#' @author Chanin Tolson
#' @seealso  \code{\link{getChangeParams}}
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get location change
#' loc = locationChange(samp_nreduce)
#'
locationChange = function(sample, point=rep(0,nrow(sample)), base=sample[1,], margin=1, tol=0.1){
  
  #set optional paramater margin
  if(missing(margin)) {
    margin = 1
  } else {
    if(!(margin %in% c(1,2))){
      warning("Margin value not valid. Margin set to default.")
      margin = 1
    }
    if(margin==2){
      sample = t(sample)
    }
  }
  
  #set optional paramater base
  if(missing(base)) {
    base = sample[1,]
  } else {    
    base = base
  }
  
  #set optional paramater point
  if(missing(point)){
    point = rep(0, nrow(sample))
  } else {
    point = point
  }
  
  #set optional paramater  tol
  if(missing(tol)){
    tol = 0.1
  } else {
    if(tol < 0){
      warning("Tol value not valid. Tol set to default.")
      tol = 0.1
    }
    tol = tol
  }
  
  #values cannot be negative (raise the minimum value to at least 0)
  if(sum(sample<0, na.rm=T)>0){
    sample = sample-min(sample, na.rm=T) 
    base = base-min(sample, na.rm=T) 
  }
  if(sum(base<0, na.rm=T)>0){
    base = base-min(base, na.rm=T)  
    sample = sample-min(base, na.rm=T) 
  }
  
  #function to get the location of pattern changes
  patternLocs = function(samp, base, tol){
    #initialize
    samp = as.numeric(samp)
    base = as.numeric(base)
    patterns = rep(0, length(samp)-1)
    patternb = rep(0, length(base)-1)
    
    #get samp pattern
    pat = (c(samp, NA) - c(NA, samp))[2:(length(samp))]
    patterns[pat>tol] = 1
    patterns[pat<(-tol)] = -1
    
    #get base pattern
    pat = (c(base, NA) - c(NA, base))[2:(length(base))]
    patternb[pat>tol] = 1
    patternb[pat<(-tol)] = -1
    
    loc = which(abs(patterns-patternb)>0, arr.ind=T)
    
    return(loc)
  }
  
  if(dim(sample)[1]==1){
    loc = NULL
    options(warn=-1)
    loc[[1]] = t(apply(sample, 1, patternLocs, base=base, tol=tol))
    options(warn=0)
    loc[is.na(loc)] = 1 
  } else{
    options(warn=-1)
    loc = apply(sample, 1, patternLocs, base=base, tol=tol)
    options(warn=0)
    loc[is.na(loc)] = 1 
  }

  
  #find average distance from the disruption
  dist = rep(0, length(loc))
  for(i in 1:length(loc)){
    if(length(loc[[i]])!=0){
      dist[i] = mean(abs(loc[[i]]-point[i]), na.rm=T)
    }
  }
  sample[,is.na(base)]=NA
  len = ncol(sample)-rowSums(is.na(sample))
  dist = dist/sqrt(len)
  
  #return average distance from the disruption
  return(dist)
}

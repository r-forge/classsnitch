#' A function to get the location change between samples
#'
#' This function determines 
#' @title locationChange
#' @aliases locationChange
#' @keywords location change RNA
#' @usage locationChange(sample, point=rep(0,nrow(sample)), base=sample[1,], margin=1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param point An optional numeric vector containing the location of the disruption (e.g. the mutation in an RNA)
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @export
#' @import changepoint
#' @import dtw
#' @details This function calculates the distance of change from the disruption using changepoint analysis to determine the location of change.
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
locationChange = function(sample, point, base=sample[1,], margin=1){
  
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
  
  #calculate timewarp change
  timewarp = function(x, y){
    warp = dtw(x, y)
    
    return(warp$index1s/warp$index2s)
  }
  tw = apply(sample, 1, timewarp, y=base)
  
  #calculate change point
  cp = sapply(tw, cpt.mean, method="PELT", class=T)
  loc = sapply(cp, cpts)
  
  #find distance from the disruption
  dist = rep(0, length(loc))
  for(i in 1:length(loc)){
    if(length(loc[[i]])!=0){
      dist[i] = mean(abs(loc[[i]]-point[i]), na.rm=T)
    }
  }
  
  #return location change
  return(dist)
}

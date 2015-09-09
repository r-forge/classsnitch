#' A function to get the average trace between samples
#'
#' This function compares the average trace difference between samples.
#' @title traceChange
#' @aliases traceChange
#' @keywords trace change RNA
#' @usage traceChange(sample, base=sample[1,], margin=1, 
#'    point=rep(0,nrow(sample)), window=ncol(sample))
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param point An optional numeric vector containing the location of the disruption (e.g. the mutation in an RNA)
#' @param window An optional number indicating the number of columns around the disruption to calculate. Default is the entire trace.
#' @export
#' @details This function calculates the average trace difference between the base vector and each row (or column) in sample.
#' @return A numeric vector of trace changes.
#' @author Chanin Tolson
#' @seealso  \code{\link{getChangeParams}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get trace change
#' tc = traceChange(samp_nreduce)
#'
traceChange = function(sample, base=sample[1,], margin=1, point=rep(0,nrow(sample)), window=ncol(sample)){
  
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
  
  #set optional paramater window
  if(missing(window)){
    window = floor(ncol(sample)/2)
  } else {
    window = floor(window/2)
  }
  
  diff = function(samp, base, point){
    d = samp-base
    
    return(d)
  }
  
  #calculate trace change
  d = apply(sample, 1, diff, base=base, point=point)
  
  #narrow down by location
  tracediff = NULL
  for(i in 1:nrow(sample)){
    if((point[i]-window)<1){
      start = 1
    } else if((point[i]-window)>=1){
      start = point[i]-window 
    }
    if((point[i]+window)>ncol(sample)){
      end = ncol(sample)
    } else if((point[i]+window)<=ncol(sample)){
      end = point[i]+window 
    }
    tracediff[i] = mean(sample[i,start:end])
  }
  
  #return magnitude change
  return(tracediff)
}

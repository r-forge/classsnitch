#' A function to get the pattern change between samples
#'
#' This function compares the pattern change between samples.
#' @title patternChange
#' @aliases patternChange 
#' @keywords pattern change RNA
#' @usage patternChange(sample, base=sample[1,], margin=1, tol=0.1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the value to which the samples are to be compared (e.g. a wild type SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @export
#' @details The pattern for a single SHAPE reacitivity trace is the pattern of increase in reactivity or decrease in reactivity between nucleotides. If the change is less than the tolerance value, it is considered a none change. The pattern change value is the Pearson correlation coefficient between the base vector pattern and the pattern of each row (or column) in sample.  
#' @return A numeric vector of pattern changes.
#' @author Chanin Tolson
#' @seealso  \code{\link{getChangeParams}}
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get patterne change
#' pat = patternChange(samp_nreduce)
#'
patternChange = function(sample, base=sample[1,], margin=1, tol=0.1){
  
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
  if(sum(sample<0)>0){
    sample = sample-min(sample, na.rm=T) 
    base = base-min(sample, na.rm=T) 
  }
  if(sum(base<0)>0){
    base = base-min(base, na.rm=T)  
    sample = sample-min(base, na.rm=T) 
  }
  
  #function to get pattern change correlation
  pattern = function(samp, base, tol){
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
    
    return(cor(patterns, patternb, method="pearson"))
  }
  
  #calculate pattern change
  pat = apply(sample, 1, pattern, base=base, tol=tol)
  pat[is.na(pat)] = 0
  
  #return pattern change
  return(pat)
}

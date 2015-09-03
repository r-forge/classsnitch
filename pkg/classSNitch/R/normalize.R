#' A Between-Sample Normalization Function
#'
#' This function performs between-sample normalization.
#' @title normalize
#' @keywords normalize between-sample RNA
#' @usage normalize(sample, base=sample[1,], margin=1)
#' @param sample A numeric matrix containing values to be normalized (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the sample is to be normalized (e.g. a wild type SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if sample is organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @export
#' @details This function normalizes the average value of the base vector to 1.5. Each row (or column) in sample is then normalized by minimizing the absolute difference between the base and the sample row (or column). 
#' @return A normalized numeric matrix with the same dimensions as sample.
#' @author Chanin Tolson
#' @seealso  \code{\link{getChangeParams}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#'            
#' #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' base = sample(1:100, size=10)
#' #normalize
#' samp_norm = normalize(sample, base)
#'
normalize = function(sample, base, margin){
  
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
  
  #set average reactivity for wild-type to 1.5 
  base = (1.5*length(base)/sum(base))*base
  
  #function to optimize difference between wild-type and sample
  optimizeNorm = function (x, samp, base){
    sum(abs(base-x*samp))
  }
  
  #loop through each sample
  opt = apply(sample, 1, optimize, f=optimizeNorm, interval=c(0,1), tol=0.0001, base=base)
  samp_norm = do.call(rbind, lapply(1:length(opt), function(i){sample[i,]*opt[[i]]$minimum}))
  
  #organize output data
  if(margin==2){
    samp_norm = t(samp_norm)
  }
  samp_norm = matrix(unlist(samp_norm), nrow=nrow(samp_norm))
  
  #return normalized sample
  return(samp_norm)
}


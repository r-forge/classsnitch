#' A function to get magnitude and pattern change parameters
#'
#' This calculates the magnitude and pattern change in SHAPE reactivity.
#' @title getChangeParams
#' @aliases getChangeParams
#' @keywords change parameters magnitude pattern RNA
#' @usage getChangeParams(sample, base=NULL, margin=1, trim=0, high=NULL, tol=0.1, outfile=NULL, append=F)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wild type SHAPE trace). Default is the first trace in each file.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param trim An optional number indicating the number of nucleotides to be trimed from the ends. Default is 0.
#' @param high An optional number indicating the reactivity above which reactivities are considered high. Default is third quartile of the sample in each file.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @param point An optional numerical vector indicating the location of disruption (e.g. mutation point)
#' @param outfile An optional string indicating the name of the output file. The output file will consist of two columns (magnitude change and pattern change). Default will not output a file.
#' @param append An optional boolean to append the file if an outfile is given. Default is FALSE. 
#' @export
#' @details This function normalizes and reduces the noise in the sample. The magnitude, pattern and location change are calculated for the sample using the magnitudeChange, patternChange and locationChange functions. 
#' @return 
#' \describe{
#'  \item{"outmat"}{A three column numeric matrix for magnitude, pattern and location change.} 
#'  \item{"outfile"}{An optional output file for the matrix.}
#' }
#' @author Chanin Tolson
#' @seealso  \code{\link{magnitudeChange}} \code{\link{patternChange}} \code{\link{normalize}} \code{\link{reduceNoise}} \code{\link{classifyRNA}} \code{\link{predict.classifyRNA}} \code{\link{locationChange}}
#' @examples #input files
#' data("shape_ex")
#' #get change parameters
#' params = getChangeParams(shape_ex, trim=5, outfile="out.txt")
#'
getChangeParams = function(sample, base=NULL, margin=1, trim=0, high=NULL, tol=0.1, point, outfile=NULL, append=F){
  
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
  
  #set optional paramater trim
  if(missing(trim)){
    trim = 0
  } else {
    if(trim < 0 || trim > dim(sample)[margin]){
      warning("Trim value not valid. Trim set to default.")
      trim = 0
    }
    trim = trim-1
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
  
  #set optional paramater point
  if(missing(point)){
    point = rep(0, nrow(sample))
  } else {
    point = point
  }
  
  #set optional paramater append
  if(missing(append)){
    append = F
  } else {
    if(!(append %in% c(TRUE, FALSE))){
      warning("Append value not valid. Append set to default.")
      append = F
    }
    append = append
  }

  #normalize
  samp_norm = normalize(sample, base)
  base = (1.5*length(base)/sum(base))*base
  
  #set optional paramater high
  if(missing(high)){
    high = boxplot(as.numeric(unlist(samp_norm)), plot=F)$stats[4]
  } else {
    if(high < 0){
      warning("High value not valid. High set to default.")
      high = boxplot(as.numeric(unlist(samp_norm)), plot=F)$stats[4]
    }
    high = high
  }
      
  #high peak filter
  samp_qual = reduceNoise(samp_norm, base, trim=trim, high=high)
      
  #magnitude change
  mag = magnitudeChange(samp_qual, base)
      
  #pattern change
  pat = patternChange(samp_qual, base, tol=tol)
  
  #location change
  loc = locationChange(samp_qual, point, base)
  
  #combine parameters
  params = cbind(cbind(mag, pat), loc)
  rownames(params) = rownames(sample)
  colnames(params) = c("magnitude change", "pattern change", "location change")
  
  #write parameter outfile
  if(missing(outfile)){
  } else{
    if(append==F){
      write.table(params, outfile, quote=F, sep="\t", row.names=T, col.names=T)
    } else{
      write.table(params, outfile, quote=F, sep="\t", row.names=T, col.names=F, append=T)
    }
  }
  
  #return magnitude and pattern change
  return(params)
}
#' A function to get the features for describing RNA structure change. These features can be used in classification of RNA structure change.
#'
#' This calculates the change in SHAPE reactivity traces.
#' @title getFeatures
#' @aliases getFeatures
#' @keywords change features RNA structure
#' @usage getFeatures(sample, base=NULL, margin=1,  norm=T, noise=T, trim=0, high=NULL,
#'    tol=0.1, point=rep(0,nrow(sample)), window=ncol(sample), outfile=NULL, append=F)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wild type SHAPE trace). Default is the first trace in each file.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param norm An optional boolean to normalize the sample. Default is TRUE. 
#' @param noise An optional boolean to reduce noise in the sample. Default is TRUE. 
#' @param trim An optional number indicating the number of nucleotides to be trimed from the ends. Default is 0.
#' @param high An optional number indicating the reactivity above which reactivities are considered high. Default is third quartile of the sample in each file.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @param point An optional numerical vector indicating the location of disruption (e.g. mutation point)
#' @param outfile An optional string indicating the name of the output file. The output file will consist of two columns (magnitude change and pattern change). Default will not output a file.
#' @param append An optional boolean to append the file if an outfile is given. Default is FALSE. 
#' @param window An optional number indicating the number of columns around the disruption to calculate. Default is the entire trace.
#' @export
#' @details This function calculates the magnitude correlation coefficient, pattern  correlation coefficient, average change distance, dynamic time warping, average trace difference and rna length. 
#' @return 
#' \describe{
#'  \item{"outmat"}{A three column numeric matrix for magnitude, pattern, location and timewarp change.} 
#'  \item{"outfile"}{An optional output file for the matrix.}
#' }
#' @author Chanin Tolson
#' @seealso  \code{\link{normalize}} \code{\link{reduceNoise}} \code{\link{getMagCC}} \code{\link{getPatternCC}} \code{\link{getChangeDist}} \code{\link{getTimeWarping}} \code{\link{getTraceDiff}}
#' @examples #input files
#' data("shape_ex")
#' #get features
#' params = getFeatures(shape_ex, trim=5, outfile="out.txt")
#'
getFeatures = function(sample, base=NULL, margin=1, norm=T, noise=T, trim=0, high=NULL, tol=0.1, point=rep(0,nrow(sample)), window=ncol(sample), outfile=NULL, append=F){
  
  #set sample parameter
  sample = as.matrix(sample)
  if(dim(sample)[1]==1 || dim(sample)[2]==1){
    sample = t(sample)
  }
  
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
  
  #set optional paramater norm
  if(missing(norm)){
    norm = TRUE
  } else {
    norm = norm
  }
  
  #set optional paramater norm
  if(missing(noise)){
    noise = TRUE
  } else {
    noise = noise
  }
  
  #set optional paramater trim
  if(missing(trim)){
    trim = 0
  } else {
    if(trim < 0 || trim > dim(sample)[2]){
      warning("Trim value not valid. Trim set to default.")
      trim = 0
    } else if(trim > 0){
      trim = trim-1 
    }
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
  
  #set optional paramater window
  if(missing(window)){
    window = floor(ncol(sample)/2)
  } else {
    window = floor(window/2)
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

  #set optional paramater normalize
  if(norm == TRUE){
    #normalize
    samp_norm = normalize(sample, base)
    base = (1.5*length(base)/sum(base, na.rm=T))*base
    base[base<(-0.5)] = 0  
  } else {
    samp_norm = sample
  }
  
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
  
  #reduce noise
  if(noise == TRUE){
    samp_qual = reduceNoise(samp_norm, base, trim=trim, high=high)
  } else{
    samp_qual = samp_norm
  }
      
  #magnitude change
  mag = getMagCC(samp_qual, base)
      
  #pattern change
  pat = getPatternCC(samp_qual, base, tol=tol)
  
  #location change
  loc = getChangeDist(samp_qual, point=point, base=base, tol=tol)
  
  #timewarp change
  tw = getTimeWarping(samp_qual, base)
  
  #timewarp change
  tc = getTraceDiff(samp_qual, base, point=point, window=window)
  
  #rna length
  len = ncol(samp_qual)-sum(is.na(samp_qual))
  
  #combine features
  features = cbind(mag, pat, loc, tw, tc, len)
  rownames(features) = rownames(sample)
  colnames(features) = c("magnitude change", "pattern change", "change distance",  "time warping", "trace difference", "length")
  
  #write parameter outfile
  if(missing(outfile)){
  } else{
    if(append==F){
      write.table(features, outfile, quote=F, sep="\t", row.names=T, col.names=T)
    } else{
      write.table(features, outfile, quote=F, sep="\t", row.names=T, col.names=F, append=T)
    }
  }
  
  #return features
  return(features)
}
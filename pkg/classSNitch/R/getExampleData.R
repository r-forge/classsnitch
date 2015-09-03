#' A function to load data for use in the classifyRNA package
#'
#' This function returning the data sets available for the classifyRNA package
#' @title getExampleData
#' @keywords example data RNA
#' @usage getExampleData(name)
#' @param name A string indicating which data set to return
#' @export
#' @details This function returns data from mutmap.Rdata, magpat_ex.Rdata or shape_ex.Rdata
#' @return A numeric matrix for mutmap, magpat_ex, or shape_ex.
#' @note The data has been gathered from the RNA Mapping Database mutate and map experiments.
#' @author Chanin Tolson
#' @references \href{http://rmdb.stanford.edu/}{RNA Mapping Database}
#' @seealso \code{\link{classifyRNA}} \code{\link{getChangeParams}} \code{\link{mutmap}} \code{\link{magpat_ex}} \code{\link{shape_ex}}
#' @examples obj = getExampleData("mutmap")
#'
getExampleData = function(name){
  
  if(!(name %in% c("mutmap", "magpat_ex", "shape_ex"))){
    stop("Data file not found.")
  }
  
  #loads data
  load(paste0(paste0(paste0(find.package("classifyRNA"), "/data/"), name), ".Rdata"))
  
  #save as variables
  obj = get(name)
  
  #return data
  return(obj)
}

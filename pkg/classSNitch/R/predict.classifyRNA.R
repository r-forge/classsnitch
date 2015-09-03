#' A function to classify rna structure change
#'
#' This function predicts RNA structure change in SHAPE data
#' @title predict.classifyRNA
#' @aliases predict, predict.classifyRNA
#' @keywords predict prediction RNA structure change
#' @usage predict(object, sample, type="response")
#' @param object A classifyRNA object (see classifyRNA function).
#' @param sample A matrix of predictors for magnitude and pattern change (e.g. output from getChangeParams())
#' @param type A string indicating what type of random forest data to return ("response", "vote", "prob"). Default is "response". 
#' @export
#' @import randomForest
#' @details This function predicts RNA structure change in SHAPE data using a random forest classifier. 
#' @return A numeric matrix of "response", "vote" or "prob" predictions.
#' @note Options for type: 
#' \itemize{
#'  \item{"response"}{ Predicted classes (classes with majority vote)}
#'  \item{"vote"}{ Vote count fraction (one column for each class and one row for each input)} 
#'  \item{"class"}{ Class probabilities (one column for each class and one row for each input)}
#' }   
#' @author Chanin Tolson
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package)
#' @seealso  \code{\link{getChangeParams}} \code{\link{classifyRNA}} \code{\link{getExampleData}} 
#' @examples #input data
#' data("magpat_ex")
#' #build classifier
#' cr = classifyRNA(classes=2)
#' #get prediction
#' cr_pred = predict(cr, magpat_ex, type="vote")
#'
predict.classifyRNA = function(object, sample, type="response"){
  
  #check paramater object
  if(attr(object, "class") == "classifyRNA"){
    rf = structure(object, class = "randomForest")
  } else {
    stop("Object must be a classifyRNA object")
  }
  
  #set optional paramater type
  if(missing(type)){
    type = "response"
  } else {
    if(!(type %in% c("response", "vote", "prob"))){
      warning("Type value not valid. Type set to default.")
      type = "response"
    }
    type = type
  }
  
  #random forest prediction
  colnames(sample) = c("mag_corr100", "pat_corr100")
  rf_pred = predict(rf, sample, type=type)
  
  #return predictions
  return(rf_pred)
}

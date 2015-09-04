#' A function to classify rna structure change
#'
#' This function predicts RNA structure change in SHAPE data
#' @title predict.classifyRNA
#' @aliases predict.classifyRNA
#' @keywords predict prediction RNA structure change
#' @usage predict(object, sample=NULL)
#' @param object A classifyRNA object (see classifyRNA function).
#' @param sample A matrix of predictors for magnitude and pattern change (e.g. output from getChangeParams())
#' @export
#' @import randomForest
#' @details This function predicts RNA structure change in SHAPE data using a random forest classifier. 
#' @return A list of "response", "vote" or "prob" predictions.
#' \describe{
#'  \item{"response"}{ Predicted classes (classes with majority vote)}
#'  \item{"vote"}{ Vote count fraction (one column for each class and one row for each input)} 
#'  \item{"class"}{ Class probabilities (one column for each class and one row for each input)}
#' }
#' @author Chanin Tolson
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package)
#' @seealso  \code{\link{getChangeParams}} \code{\link{classifyRNA}}
#' @examples #input data
#' data("magpat_ex")
#' #build classifier
#' cr = classifyRNA(classes=2)
#' #get prediction
#' cr_pred = predict(cr, magpat_ex[,2:4])
#'
predict.classifyRNA = function(object, ...){
  
  #check paramater object
  if(attr(object, "class") == "classifyRNA"){
    rf = structure(object, class = "randomForest")
  } else {
    stop("Object must be a classifyRNA object")
  }
  
  #get optional parameters
  rf_pred = NULL
  opt = list(...)
  if(length(opt)>0){
    sample = opt[[1]]
    colnames(sample) = c("mag", "pat", "loc")
    rf_pred[[1]] = predict(rf, sample, type="response")
    rf_pred[[2]] = predict(rf, sample, type="vote", norm.votes=F)
    rf_pred[[3]] = predict(rf, sample, type="prob")
  } else{
    sample = NULL 
    rf_pred[[1]] = predict(rf, type="response")
    rf_pred[[2]] = predict(rf, type="vote", norm.votes=F)
    rf_pred[[3]] = predict(rf, type="prob")
  }
  names(rf_pred) = c("response", "vote", "prob")

  #return predictions
  return(rf_pred)
}

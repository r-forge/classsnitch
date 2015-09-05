#' A function to classify rna structure change
#'
#' This function predicts RNA structure change in SHAPE data
#' @title predict.classifyRNA
#' @keywords predict prediction RNA structure change
#' @usage predict(object, sample=NULL, ...)
#' @param object An object of classifyRNA (see classifyRNA function).
#' @param sample An optional matrix of predictors for magnitude, pattern, location and pattern change (e.g. output from getChangeParams())
#' @param ...	Further arguments passed to or from other methods.
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
#' cr_pred = predict(cr, magpat_ex[,2:5])
#'
predict.classifyRNA = function(object, sample=NULL, ...){
  
  #check paramater object
  if(attr(object, "class") == "classifyRNA"){
    rf = structure(object, class = "randomForest")
  } else {
    stop("Object must be a classifyRNA object")
  }
  
  #get optional parameters
  rf_pred = NULL
  if(!missing(sample)){
    colnames(sample) = c("mag", "pat", "loc", "tw")
    rf_pred[[1]] = predict(rf, sample, type="response")
    rf_pred[[2]] = predict(rf, sample, type="vote", norm.votes=F)
    rf_pred[[3]] = predict(rf, sample, type="prob")
  } else{
    rf_pred[[1]] = predict(rf, type="response")
    rf_pred[[2]] = predict(rf, type="vote", norm.votes=F)
    rf_pred[[3]] = predict(rf, type="prob")
  }
  names(rf_pred) = c("response", "vote", "prob")

  #return predictions
  return(rf_pred)
}

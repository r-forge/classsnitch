#' A function to build a classifier for RNA structure change
#'
#' This function builds a random forest classifier for RNA structure change in SHAPE data
#' @title classifyRNA
#' @aliases classifyRNA
#' @keywords classifier RNA structure change random forest
#' @usage classifyRNA(data=NULL, classes=2)
#' @param data Optional data to build the classifier. Default is pre-loaded data.
#' @param classes An optional number indicating which class style to use. Only used when data is not supplied. Default is 2.
#' @export
#' @export classify_default
#' @import randomForest
#' @details This function builds a random forest classifier for RNA structure change using the randomForest package.
#' @return A classifyRNA object, based on randomForest object (see randomForest package)
#' \describe{
#'  \item{call}{The original call to randomForest}
#'  \item{type}{One of regression, classification, or unsupervised.}
#'  \item{predicted}{The predicted values of the input data based on out-of-bag samples.}
#'  \item{importance}{A matrix with nclass + 2 columns. The first nclass columns are the class-specific measures computed as mean descrease in accuracy. The nclass + 1st column is the mean descrease in accuracy over all classes. The last column is the mean decrease in Gini index.}
#'  \item{importanceSD}{The standard errors of the permutation-based importance measure. A p by nclass + 1 matrix corresponding to the first nclass + 1 columns of the importance matrix.}
#'  \item{ntree}{Number of trees grown.}
#'  \item{mtry}{Number of predictors sampled for spliting at each node.}
#'  \item{forest}{A list that contains the entire forest}
#'  \item{err.rate}{Vector error rates of the prediction on the input data, the i-th element being the (OOB) error rate for all trees up to the i-th.}
#'  \item{confusion}{The confusion matrix of the prediction (based on OOB data).}
#'  \item{votes}{A matrix with one row for each input data point and one column for each class, giving the fraction or number of (OOB) votes from the random forest.}
#'  \item{oob.times}{Number of times cases are out-of-bag (and thus used in computing OOB error estimate)}
#'  \item{proximity}{A matrix of proximity measures among the input (based on the frequency that pairs of data points are in the same terminal nodes).}
#' }
#' @note Organization of the data file: header=TRUE, tab-delimited .txt file
#' \itemize{
#'  \item{"column 1"}{ RNA name} 
#'  \item{"column 2"}{ Classification} 
#'  \item{"column 3"}{ Magnitude predictor} 
#'  \item{"column 4"}{ Pattern predictor} 
#' }
#' Options for classes: 
#' \itemize{
#'  \item{"1"}{ none v. local v. global}
#'  \item{"2"}{ none v. local/global}
#'  \item{"3"}{ global v. none/local} 
#'  \item{"4"}{ local v. global}
#' }
#' The default data has been gathered from the RNA Mapping Database mutate and map experiments.
#' @author Chanin Tolson
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package) \cr\cr
#' \href{http://rmdb.stanford.edu/}{RNA Mapping Database}
#' @seealso  \code{\link{getChangeParams}} \code{\link{predict.classifyRNA}} \code{\link{getExampleData}}
#' @examples
#' #build classifier
#' rf = classifyRNA(classes=2)
#' #get confusion matrix
#' rf$confusion
#' 
classifyRNA = function(data=NULL, classes=2){
  
  #set optional paramater classes
  if(missing(classes)){
    classes = 2
  } else {
    if(!(classes %in% c(1,2,3,4))){
      warning("classes set to default.")
      classes = 2
    }
    classes = classes
  }
  
  #set optional paramater data
  if(missing(data)){
    data = classSNitch::classify_default
    responses = data[,classes]
    input = data[,8:9]
  } else {
    data = data
    if(ncol(data) != 3){
      stop("Incorrect data file format.")
    }
    responses = data[,1]
    input = data[,2:3]
  }
  
  #get parameters
  input = cbind(as.factor(responses), input)
  rownames(input) = rownames(data)
  colnames(input) = c("class", "mag_corr100", "pat_corr100") 
  input = input[-which(is.na(input[,1]),arr.ind=T),]
  
  #random forest classification
  rf = randomForest(class~., data=input, importance=TRUE, proximity=TRUE, ntree=5001)
  
  #convert to a classifyRNA object
  cr = structure(rf, class = "classifyRNA")
  
  #return classifyRNA object
  return(cr)
}


#' Package for the autonomous classification of RNA structure change
#'
#' Mutations in RNA will create a riboSNitch, if important structural elements are disrupted. Recent ultra-high throughput techniques, such as SHAPE-MaP and PARS, enable the collection of structural RNA information on a genome-wide scale. With the ability to gather genome-wide structural information on RNA, it is important to accurately classify these structural data in order to identify those structural changes that result in a phenotypic outcome. We have developed an automated approach to classify structure change in SHAPE data. This method utilizes random forest classification on a set of pattern and magnitude parameters from the mutate-and-map SHAPE data set (or another user specified data set) to build a classifier. The classifier is then used to identify structure change in other SHAPE traces. Enabling scientists to identify structure change may help guide experiments that examine RNA structure and its role in biological processes. 
#'
#' @docType package
#' @name classSNitch
#' @author Chanin Tolson
#' @import gplots
#' @import ROCR
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package) \cr\cr
#' \href{http://rmdb.stanford.edu/}{RNA Mapping Database}
#' @seealso  \code{\link{getFeatures}} \code{\link{classifyRNA}}
#' @examples #get change features
#' library("ROCR")
#' library("gplots")
#' 
#' data("shape_ex")
#' sample = getFeatures(shape_ex[2:nrow(shape_ex),], base=shape_ex[1,], trim=5)
#' 
#' #predict change
#' data("mutmap")
#' cr = classifyRNA(mutmap)
#' cr_pred = predict(cr, sample, type="response")
#' 
#' #plot ROC curve (no change v. local/global change)
#' data("mutmap")
#' predobj = prediction(cr$votes[,1], mutmap[,1]==1)
#' perfobj = performance(predobj, 'tpr', 'fpr')
#' aucobj = performance(predobj, 'auc')
#' plot(perfobj@@x.values[[1]], perfobj@@y.values[[1]], lwd=2, 
#'      type="l", xlab="Specificity", ylab="Sensitivity")
#' points(c(-1,2),c(-1,2), col="red", type="l")
#' text(0.8, 0.2, paste("AUC: ", format(aucobj@@y.values, digits=2), sep=""), cex=1)
#' 
NULL
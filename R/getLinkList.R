#' @title getLinkList
#'
#' @description \code{getLinkList} Converts the weight matrix, as returned by \code{\link{GENIE3}}, to a sorted list of regulatory links (most likely links first).
#'
#' @param weightMatrix Weighted adjacency matrix as returned by \code{\link{GENIE3}}.
#' @param reportMax Maximum number of links to report. The default value NULL means that all the links are reported.
#' @param threshold Only links with a weight equal or above the threshold are reported. Default: threshold = 0, i.e. all the links are reported.
#'
#' @return List of regulatory links in a data frame. Each line of the data frame corresponds to a link. The first column is the regulatory gene, the second column is the target gene, and the third column is the weight of the link.
#'
#' @seealso \code{\link{GENIE3}}
#'
#' @examples
#' ## Generate fake expression matrix
#' exprMat <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(exprMat) <- paste("Gene", 1:20, sep="")
#' colnames(exprMat) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weightMat <- GENIE3(exprMat, regulators=paste("Gene", 1:5, sep=""))
#'
#' ## Get ranking of edges
#' linkList <- getLinkList(weightMat)
#' head(linkList)
#' @export
getLinkList <- function(weightMatrix, reportMax=NULL, threshold=0) {
    if(!is.numeric(threshold)) {
    	stop("threshold must be a number.")
    }

	# Only process weights off-diagonal
  regulatorsInTargets <- rownames(weightMatrix)[rownames(weightMatrix) %in% colnames(weightMatrix)]
  if(length(regulatorsInTargets)==1) weightMatrix[regulatorsInTargets, regulatorsInTargets] <- NA
  if(length(regulatorsInTargets)>1) diag(weightMatrix[regulatorsInTargets, regulatorsInTargets]) <- NA
  
  linkList <- reshape2::melt(weightMatrix, na.rm=TRUE)
  colnames(linkList) <- c("regulatoryGene", "targetGene", "weight")
  linkList <- linkList[linkList$weight>=threshold,]
  linkList <- linkList[order(linkList$weight, decreasing=TRUE),]
  
  if(!is.null(reportMax)) {
  	linkList <- linkList[1:min(nrow(linkList), reportMax),]
  }
  rownames(linkList) <- NULL
  
  uniquePairs <- nrow(unique(linkList[,c("regulatoryGene", "targetGene")]))
  if(uniquePairs < nrow(linkList)) 
    warning("There might be duplicated regulator-target (gene id/name) pairs.")
  
  return(linkList)
}

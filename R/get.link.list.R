#' @title get.link.list
#' 
#' @description \code{get.link.list} Converts the weight matrix, as returned by \code{\link{GENIE3}} or \code{\link{tlGENIE3}}, to a sorted list of regulatory links (most likely links first).
#' 
#' @param weight.matrix Weighted adjacency matrix as returned by \code{\link{GENIE3}} or \code{\link{tlGENIE3}}.
#' @param report.max Maximum number of links to report. The default value NULL means that all the links are reported.
#' @param threshold Only links with a weight above the threshold are reported. Default: threshold = 0, i.e. all the links are reported.
#' 
#' @return List of regulatory links in a data frame. Each line of the data frame corresponds to a link. The first column is the regulatory gene, the second column is the target gene, and the third column is the weight of the link.
#'
#' @seealso \code{\link{GENIE3}}, \code{\link{tlGENIE3}}
#'
#' @examples
#' ## Generate fake expression matrix
#' expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
#' colnames(expr.matrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weight.matrix <- GENIE3(expr.matrix, regulators=paste("Gene", 1:5, sep=""))
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
get.link.list <- function(weight.matrix, report.max=NULL, threshold=0) {
    if(!is.numeric(threshold)) {
    	stop("threshold must be a number.")
    } 

	# Only process weights off-diagonal
	diag(weight.matrix) <- NA
    link.list <- melt(weight.matrix, na.rm=TRUE)
    colnames(link.list) <- c("regulatory.gene", "target.gene", "weight")
    link.list <- link.list[link.list$weight>=threshold,]
    link.list <- link.list[order(link.list$weight, decreasing=TRUE),]
  
    if(!is.null(report.max)) {
    	link.list <- link.list[1:min(nrow(link.list), report.max),]
    } 
  
    rownames(link.list) <- NULL
  
    return(link.list)
}

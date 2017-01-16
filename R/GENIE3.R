#' @title GENIE3
#'
#' @description \code{GENIE3} Infers a gene regulatory network (in the form of a weighted adjacency matrix) from expression data, using ensembles of regression trees.
#'
#' @param expr.matrix Expression matrix (genes x samples). Every row is a gene, every column is a sample.
#' @param tree.method Tree-based method used. Must be either "RF" for Random Forests (default) or "ET" for Extra-Trees.
#' @param K Number of candidate regulators randomly selected at each tree node (for the determination of the best split). Must be either "sqrt" for the square root of the total number of candidate regulators (default), "all" for the total number of candidate regulators, or a stricly positive integer.
#' @param ntrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param regulators Subset of genes used as candidate regulators. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param ncores Number of cores to use for parallel computing. Default: 1.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: FALSE.
#' @param seed Random number generator seed for replication of analyses. The default value NULL means that the seed is not reset.
#'
#' @return Weighted adjacency matrix of inferred network. Element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene j.
#'
#' @examples
#' ## Generate fake expression matrix
#' exprMatrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
#' colnames(exprMatrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weightMatrix <- GENIE3(exprMatrix, regulators=paste("Gene", 1:5, sep=""))
#'
#' ## Get ranking of edges
#' linkList <- get.link.list(weightMatrix)
#' head(linkList)
#' @export
GENIE3 <- function(expr.matrix, tree.method="RF", K="sqrt", ntrees=1000, regulators=NULL, ncores=1, verbose=FALSE, seed=NULL)
{
  ############################################################
	# check input arguments
	if (!is.matrix(expr.matrix) && !is.array(expr.matrix)) {
		stop("Parameter expr.matrix must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
	}

	if (length(dim(expr.matrix)) != 2) {
		stop("Parameter expr.matrix must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
	}

	if (is.null(rownames(expr.matrix))) {
		stop("expr.matrix must specify the names of the genes in rownames(expr.matrix).")
	}

	if (tree.method != "RF" && tree.method != "ET") {
		stop("Parameter tree.method must be \"RF\" (Random Forests) or \"ET\" (Extra-Trees).")
	}

	if (K != "sqrt" && K != "all" && !is.numeric(K)) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}

	if (is.numeric(K) && K<1) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}

	if (!is.numeric(ntrees) || ntrees<1) {
		stop("Parameter ntrees should be a stricly positive integer.")
	}

	if (!is.null(regulators)) {
	  if(length(regulators)<2) stop("Provide at least 2 potential regulators.")

		if (!is.vector(regulators)) {
			stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
		}

		if (is.character(regulators) && length(intersect(regulators,rownames(expr.matrix))) == 0) {
			stop("The genes must contain at least one candidate regulator.")
		}

		if (is.numeric(regulators) && max(regulators) > nrow(expr.matrix)) {
			stop("At least one index in regulators exceeds the number of genes.")
		}

    if(is.numeric(regulators))
    {
        regulators <- rownames(expr.matrix)[regulators]
    }
	}

	if (!is.numeric(ncores) || ncores<1) {
		stop("Parameter ncores should be a stricly positive integer.")
	}
  ############################################################


  # transpose expression matrix to (samples x genes)
  expr.matrix <- t(expr.matrix)
  num.samples <- nrow(expr.matrix)
  num.genes <- ncol(expr.matrix)
  gene.names <- colnames(expr.matrix)

  # get names of input genes
  if(is.null(regulators)) {
    input.gene.names <- gene.names
  } else {
    # input gene indices given as integers
    if (is.numeric(regulators)) {
      input.gene.names <- gene.names[regulators]
      # input gene indices given as names
    } else {
      input.gene.names <- regulators
      # for security, abort if some input gene name is not in gene names
      missing.gene.names <- setdiff(input.gene.names, gene.names)
      if (length(missing.gene.names) != 0) stop(paste("Genes missing in the expression matrix:", paste(missing.gene.names, collapse=", ")))
    }
  }
  rm(regulators)

  # set random number generator seed if seed is given
  if (!is.null(seed)) set.seed(seed)

	# tree method
	if (tree.method == 'RF') {
		RF_randomisation <- 1
		ET_randomisation <- 0
		bootstrap_sampling <- 1
	} else {
		RF_randomisation <- 0
		ET_randomisation <- 1
		bootstrap_sampling <- 0
	}

	if (verbose) message(paste("Tree method: ", tree.method,
	                           "\nK: ", K,
	                            "\nNumber of trees: ", ntrees, sep=""))
  # other default parameters
  nmin <- 1
  permutation_importance <- 0

  # setup weight matrix
  weight.matrix <- matrix(0.0, nrow=length(input.gene.names), ncol=num.genes)
  rownames(weight.matrix) <- input.gene.names
  colnames(weight.matrix) <- gene.names

  # compute importances for every target gene
  if (ncores==1) {
  	# serial computing
  	if (verbose) message("Using 1 core.")
    for (target.gene.name in gene.names) {
      if (verbose) message(paste("Computing gene ", which(gene.names == target.gene.name), "/", num.genes, ": ",target.gene.name, sep=""))

      # remove target gene from input genes
      these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
      num.input.genes <- length(these.input.gene.names)
      mtry <- setMtry(K, num.input.genes)

      x <- expr.matrix[,these.input.gene.names]
      y <- expr.matrix[,target.gene.name]

      im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
              as.single(c(x)),as.single(c(y)),as.integer(nmin),
              as.integer(ET_randomisation),as.integer(RF_randomisation),
              as.integer(mtry),as.integer(ntrees),
              as.integer(bootstrap_sampling),as.integer(permutation_importance),
              as.double(vector("double",num.input.genes)))[[12]]

      # normalize variable importances
      im <- im / sum(im)
      weight.matrix[these.input.gene.names, target.gene.name] <- im
    }
	} else {
	  #requireNamespace("foreach"); requireNamespace("doRNG"); requireNamespace("doParallel")
	  require("foreach"); require("doRNG"); require("doParallel")

		# parallel computing
	  doParallel::registerDoParallel(); options(cores=ncores)
    if (verbose) message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))

    weight.matrix.reg <- foreach::foreach(target.gene.name=gene.names, .combine=cbind) %dorng%
	  # weight.matrix.reg <- doRNG::"%dorng%"(foreach::foreach(target.gene.name=gene.names, .combine=cbind),
    {
      # remove target gene from input genes
      these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
      num.input.genes <- length(these.input.gene.names)
      mtry <- setMtry(K, num.input.genes)

      x <- expr.matrix[,these.input.gene.names]
      y <- expr.matrix[,target.gene.name]

      im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
        as.single(c(x)),as.single(c(y)),as.integer(nmin),
        as.integer(ET_randomisation),as.integer(RF_randomisation),
        as.integer(mtry),as.integer(ntrees),
        as.integer(bootstrap_sampling),as.integer(permutation_importance),
        as.double(vector("double",num.input.genes)))[[12]]

      # normalize variable importances
      im <- im / sum(im)

      c(setNames(0, target.gene.name), setNames(im, these.input.gene.names))[input.gene.names]
    }
    attr(weight.matrix.reg, "rng") <- NULL
    weight.matrix[input.gene.names,] <- weight.matrix.reg
	}
  return(weight.matrix)
}

# mtry <- setMtry(K, num.input.genes)
setMtry <- function(K, num.input.genes)
{
  # set mtry
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.input.genes))
  } else {
    mtry <- num.input.genes
  }

  return(mtry)
}


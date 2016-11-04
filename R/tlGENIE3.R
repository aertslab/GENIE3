#' @title tlGENIE3
#' 
#' @description \code{tlGENIE3} Infers a gene regulatory network (in the form of a weighted adjacency matrix) from time series of expression data, using ensembles of regression trees. The method can also be used to infer a network jointly from timeseries and steady-state expression data.
#'
#' @param TS.data List of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment. Each row of a matrix is a gene, each column is a time point. The number of time points does not need to be the same in all the experiments, but the genes must be the same.
#' @param SS.data Steady-state expression matrix (genes x samples). Every row is a gene, every column is a sample. The default value NULL means that no steady-state data are used for network inference.
#' @param tree.method Tree-based method used. Must be either "RF" for Random Forests (default) or "ET" for Extra-Trees.
#' @param K Number of candidate regulators randomly selected at each tree node (for the determination of the best split). Must be either "sqrt" for the square root of the total number of candidate regulators (default), "all" for the total number of candidate regulators, or a stricly positive integer.
#' @param ntrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param h Value of the time lag (i.e. number of time points between the expressions of the candidate regulators and the expression of the target gene). Default: 1.
#' @param regulators Subset of genes used as candidate regulators. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param ncores Number of cores to use for parallel computing. Default: 1.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: FALSE.
#' @param seed Random number generator seed for replication of analyses. The default value NULL means that the seed is not reset.
#'
#' @return Weighted adjacency matrix of inferred network. Element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene j. 
#' 
#' @examples
#' ## Generate fake time series of expression data
#' data1 <- matrix(sample(1:10, 70, replace=TRUE), nrow=10)
#' rownames(data1) <- paste("Gene", 1:10, sep="")
#' colnames(data1) <- paste("Time point", 1:7, sep="")
#'
#' data2 <- matrix(sample(1:10, 50, replace=TRUE), nrow=10)
#' rownames(data2) <- paste("Gene", 1:10, sep="")
#' colnames(data2) <- paste("Time point", 1:5, sep="")
#'
#' TS.data <- list(data1,data2) 
#'
#' ## Run tlGENIE3
#' weight.matrix <- tlGENIE3(TS.data)
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
tlGENIE3 <- function(TS.data, SS.data=NULL, tree.method="RF", K="sqrt", ntrees=1000, h=1, regulators=NULL, ncores=1, verbose=FALSE, seed=NULL) {

	# check input arguments
	
	if (!is.list(TS.data) && !is.vector(TS.data)){
		stop("Parameter TS.data must be a list of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment.")
	}
	
	for (k in seq(from=1, to=length(TS.data))) {
		if (!is.matrix(TS.data[[k]]) && !is.array(TS.data[[k]])) {
			stop("Parameter TS.data must be a list of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment.")
		}
	
		if (length(dim(TS.data[[k]])) != 2) {
			stop("Parameter TS.data must be a list of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment.")
		}
	
		if (is.null(rownames(TS.data[[k]]))) {
			stop("Each k-th matrix of TS.data must specify the names of the genes in rownames(TS.data[[k]]).")
		}
	}
	
	gene.names <- rownames(TS.data[[1]])
	num.genes <- length(gene.names)
	
	for (k in seq(from=2, to=length(TS.data))){
		if (length(union(gene.names, rownames(TS.data[[k]]))) != num.genes) {
			stop("The genes/rows must be the same in all the matrices of TS.data.")
		}
	}
	
	if (!is.null(SS.data)) {
		
		if (!is.matrix(SS.data) && !is.array(SS.data)) {
			stop("Parameter SS.data must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
		}
	
		if (length(dim(SS.data)) != 2) {
			stop("Parameter SS.data must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
		}
	
		if (is.null(rownames(SS.data))) {
			stop("SS.data must specify the names of the genes in rownames(SS.data).")
		}
		
		if (length(intersect(gene.names, rownames(SS.data))) != num.genes) {
			stop("The genes/rows of SS.data must be the same as the genes in TS.data.")
		}
		
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
	
	if (!is.numeric(h) || h<1) {
		stop("Parameter h should be a stricly positive integer.")
	}
	
	num.samples.min <- dim(TS.data[[1]])[2]
	for (k in seq(from=2, to=length(TS.data))) {
		num.samples.min <- min(num.samples.min,dim(TS.data[[k]])[2])
	}
	if (h >= num.samples.min) {
		stop("Parameter h must be lower than the lowest number of time points in the time series experiments.")
	}
	
	if (!is.null(regulators)) {
		if (!is.vector(regulators)) {
			stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
		}
		
		if (is.character(regulators) && length(intersect(regulators,gene.names)) == 0) {
			stop("The genes must contain at least one candidate regulator.")
		}
		
		if (is.numeric(regulators) && max(regulators) > num.genes) {
			stop("At least one index in regulators exceeds the number of genes.")
		}
	}
	
	if (!is.numeric(ncores) || ncores<1) {
		stop("Parameter ncores should be a stricly positive integer.")
	}
	
	
	
	# set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # transpose all expression matrices to (samples x genes)
    for (k in seq(from=1, to=length(TS.data))) {
    	TS.data[[k]] <- t(TS.data[[k]])
    }
	if (!is.null(SS.data)) {
		SS.data <- t(SS.data)
	}
	
    # setup weight matrix
    weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
    rownames(weight.matrix) <- gene.names
    colnames(weight.matrix) <- gene.names
	
    # get names of input genes
    if (is.null(regulators)) {
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
            if (length(missing.gene.names) != 0) {
                for (missing.gene.name in missing.gene.names) {
                    cat(paste("Gene ", missing.gene.name,
                              " was not in the expression matrix\n", sep=""))
                }
                stop("Aborting computation")
            }
        }
    }
	
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
	
	if (verbose) {
        cat(paste("Tree method: ", tree.method, "\nK: ", K,
	              "\nNumber of trees: ", ntrees, "\n\n",
                  sep=""))
        flush.console()
	}
	
	
	# generate time-lagged inputs and outputs
	num.exp <- length(TS.data)
	num.samples.time <- 0
	for (k in seq(from=1, to=length(TS.data))) {
		num.samples.time <- num.samples.time + dim(TS.data[[k]])[1]
	}
	
	input.matrix.time <- matrix(0.0, nrow=num.samples.time-h*num.exp, ncol=num.genes)
	output.matrix.time <- matrix(0.0, nrow=num.samples.time-h*num.exp, ncol=num.genes)
	
	num.samples.count <- 1
	
	for (k in seq(from=1, to=length(TS.data))){
		current.timeseries <- TS.data[[k]][,gene.names]
		num.points <- dim(current.timeseries)[1]
		current.timeseries.input <- current.timeseries[1:(num.points-h),]
		current.timeseries.output <- current.timeseries[(h+1):num.points,]
		num.samples.current <- dim(current.timeseries.input)[1]
		input.matrix.time[num.samples.count:(num.samples.count+num.samples.current-1),] <- current.timeseries.input
		output.matrix.time[num.samples.count:(num.samples.count+num.samples.current-1),] <- current.timeseries.output
		num.samples.count <- num.samples.count + num.samples.current
	}
	
	# add steady-state data (if any)
	if (!is.null(SS.data)) {
		input.all <- rbind(SS.data[,gene.names],input.matrix.time)
		output.all <- rbind(SS.data[,gene.names],output.matrix.time)
	} else {
		input.all <- input.matrix.time
		output.all <- output.matrix.time
	}
	
	colnames(input.all) <- gene.names
	colnames(output.all) <- gene.names
	
	num.samples = dim(input.all)[1]
    
    # compute importances for every target gene
   
	if (ncores==1) {
		# serial computing
		if (verbose) {
		    cat("Using 1 core.\n\n")
		    flush.console()
		}
		
	    for (target.gene.idx in seq(from=1, to=num.genes)) {

            if (verbose) {	
                cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
                flush.console()
			 }

	        target.gene.name <- gene.names[target.gene.idx]
	        # Add target gene to input genes
			these.input.gene.names <- c(target.gene.name, setdiff(input.gene.names,target.gene.name))
			num.input.genes <- length(these.input.gene.names)
		
	        x <- input.all[,these.input.gene.names]
			y <- output.all[,target.gene.name]

		    # set mtry
		    if (class(K) == "numeric") {
		        mtry <- K
		    } else if (K == "sqrt") {
		        mtry <- round(sqrt(num.input.genes))
		    } else {
		        mtry <- num.input.genes
		    } 
		
			# some default parameters 
			nmin <- 1
			permutation_importance <- 0
		
	        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
			          as.single(c(x)),as.single(c(y)),as.integer(nmin),
					  as.integer(ET_randomisation),as.integer(RF_randomisation),
					  as.integer(mtry),as.integer(ntrees),
					  as.integer(bootstrap_sampling),as.integer(permutation_importance),
					  as.double(vector("double",num.input.genes)))[[12]]
					  
			im <- setNames(im,these.input.gene.names)
			im[target.gene.name] <- 0
			
			# normalize variable importances
			im <- im / sum(im)
			
	        weight.matrix[these.input.gene.names, target.gene.name] <- im
	    }
	} else {
		# parallel computing
	    registerDoParallel(); options(cores=ncores)
		
		if (verbose) {
		    message(paste("\nUsing", getDoParWorkers(), "cores."))
		}
		
	    weight.matrix.reg <- foreach(target.gene.name=gene.names, .combine=cbind) %dorng% 
	    {
	        # Add target gene to input genes
			these.input.gene.names <- c(target.gene.name, setdiff(input.gene.names,target.gene.name))
			num.input.genes <- length(these.input.gene.names)
		
	        x <- input.all[,these.input.gene.names]
			y <- output.all[,target.gene.name]

		    # set mtry
		    if (class(K) == "numeric") {
		        mtry <- K
		    } else if (K == "sqrt") {
		        mtry <- round(sqrt(num.input.genes))
		    } else {
		        mtry <- num.input.genes
		    } 
			
			# some default parameters 
			nmin <- 1
			permutation_importance <- 0
		
	        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
			          as.single(c(x)),as.single(c(y)),as.integer(nmin),
					  as.integer(ET_randomisation),as.integer(RF_randomisation),
					  as.integer(mtry),as.integer(ntrees),
					  as.integer(bootstrap_sampling),as.integer(permutation_importance),
					  as.double(vector("double",num.input.genes)))[[12]]
					  
		  	im <- setNames(im,these.input.gene.names)
		  	im[target.gene.name] <- 0
			
		  	# normalize variable importances
		  	im <- im / sum(im)
					  
			im[input.gene.names]
	    }
	    attr(weight.matrix.reg, "rng") <- NULL
	    weight.matrix[input.gene.names,] <- weight.matrix.reg
	}
    return(weight.matrix)
}       
	
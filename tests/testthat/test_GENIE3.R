test_that("GENIE3 tests", {

	exprMatrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
	rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
	colnames(exprMatrix) <- paste("Sample", 1:5, sep="")

	# Default parameters
	weightMatrix <- GENIE3(exprMatrix)

	expect_equal(dim(weightMatrix)[1], 20)
	expect_equal(dim(weightMatrix)[2], 20)
	expect_true(is.numeric(weightMatrix))
	expect_equal(sum(is.na(weightMatrix)), 0)
	expect_equal(sum(weightMatrix < 0), 0)
	expect_equal(sum(diag(weightMatrix)), 0)
	expect_equal(sum(rownames(weightMatrix) == colnames(weightMatrix)), 20)
	expect_equal(sum(rownames(weightMatrix) == rownames(exprMatrix)), 20)

    #### Regulators as number
	regulators <- c(2,4,6)
	weightMatrix <- GENIE3(exprMatrix, nTree=100, regulators=regulators)
	expect_equal(sum(weightMatrix < 0), 0)
	expect_equal(nrow(weightMatrix), 3)
	expect_equal(ncol(weightMatrix), 20)
	expect_true(is.numeric(weightMatrix))
	expect_equal(sum(is.na(weightMatrix)), 0)
	expect_equal(sum(diag(weightMatrix[,rownames(weightMatrix)])), 0)
	# expect_equal(sum(rownames(weightMatrix) == colnames(weightMatrix)), 20)
	expect_equal(sum(colnames(weightMatrix) == rownames(exprMatrix)), 20)

	#### Regulators as name
	regulators <- rownames(exprMatrix)[c(5,9)]
	weightMatrix <- GENIE3(exprMatrix, treeMethod="ET", K=1, nTree=100, regulators=regulators)
	expect_equal(sum(weightMatrix < 0), 0)
	expect_equal(nrow(weightMatrix), 2)
	expect_equal(ncol(weightMatrix), 20)
	expect_true(is.numeric(weightMatrix))
	expect_equal(sum(is.na(weightMatrix)), 0)
	expect_equal(sum(diag(weightMatrix[,rownames(weightMatrix)])), 0)
	# expect_equal(sum(rownames(weightMatrix) == colnames(weightMatrix)), 20)
	expect_equal(sum(colnames(weightMatrix) == rownames(exprMatrix)), 20)
	
	#### Targets as number
	regulators <- c(2,4,6)
	targets <- c(1,5)
	weightMatrix <- GENIE3(exprMatrix, regulators=regulators, targets=targets)
	expect_equal(nrow(weightMatrix), 3)
	expect_equal(ncol(weightMatrix), 2)
	
	#### Targets as name
	targets <- rownames(exprMatrix)[c(5,9)]
	weightMatrix <- GENIE3(exprMatrix, targets=targets)
	expect_equal(nrow(weightMatrix), 20)
	expect_equal(ncol(weightMatrix), 2)
	
	#### Only one target (Note that the regulator is included, and will take all the weight...)
	weightMatrix <- GENIE3(exprMatrix, regulators=2:3, targets=3)
	expect_equal(ncol(weightMatrix), 1)
	
	
	### Other input classes:
	eset <- Biobase::ExpressionSet(assayData=exprMatrix)
	expect_equal(class(GENIE3(eset)), "matrix")
	
	sexp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exprMatrix))
	expect_warning(wm <- GENIE3(sexp))
	expect_equal(class(wm), "matrix")
	
	# sce <- scater::newSCESet(countData=exprMatrix)
	# expect_equal(class(GENIE3(sce)), "matrix")
})


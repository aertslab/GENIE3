test_that("GENIE3 tests", {

	exprMatrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
	rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
	colnames(exprMatrix) <- paste("Sample", 1:5, sep="")

	# Default parameters
	set.seed(123)
	weightMatrix <- GENIE3(exprMatrix, nCores=1)

	testthat::expect_equal(dim(weightMatrix)[1], 20)
	testthat::expect_equal(dim(weightMatrix)[2], 20)
	testthat::expect_true(is.numeric(weightMatrix))
	testthat::expect_equal(sum(is.na(weightMatrix)), 0)
	testthat::expect_equal(sum(weightMatrix < 0), 0)
	testthat::expect_equal(sum(diag(weightMatrix)), 0)
	testthat::expect_equal(sum(rownames(weightMatrix) == colnames(weightMatrix)), 20)
	testthat::expect_equal(sum(rownames(weightMatrix) == rownames(exprMatrix)), 20)

  #### Regulators as number
	regulators <- c(2,4,6)
	weightMatrix <- GENIE3(exprMatrix, nTree=100, regulators=regulators)
	testthat::expect_equal(sum(weightMatrix < 0), 0)
	testthat::expect_equal(nrow(weightMatrix), 3)
	testthat::expect_equal(ncol(weightMatrix), 20)
	testthat::expect_true(is.numeric(weightMatrix))
	testthat::expect_equal(sum(is.na(weightMatrix)), 0)
	testthat::expect_equal(sum(diag(weightMatrix[,rownames(weightMatrix)])), 0)
	# testthat::expect_equal(sum(rownames(weightMatrix) == colnames(weightMatrix)), 20)
	testthat::expect_equal(sum(colnames(weightMatrix) == rownames(exprMatrix)), 20)

	#### Regulators as name
	regulators <- rownames(exprMatrix)[c(5,9)]
	weightMatrix <- GENIE3(exprMatrix, treeMethod="ET", K=1, nTree=100, regulators=regulators)
	testthat::expect_equal(sum(weightMatrix < 0), 0)
	testthat::expect_equal(nrow(weightMatrix), 2)
	testthat::expect_equal(ncol(weightMatrix), 20)
	testthat::expect_true(is.numeric(weightMatrix))
	testthat::expect_equal(sum(is.na(weightMatrix)), 0)
	testthat::expect_equal(sum(diag(weightMatrix[,rownames(weightMatrix)])), 0)
	# testthat::expect_equal(sum(rownames(weightMatrix) == colnames(weightMatrix)), 20)
	testthat::expect_equal(sum(colnames(weightMatrix) == rownames(exprMatrix)), 20)

	#### Targets as number
	regulators <- c(2,4,6)
	targets <- c(1,5)
	weightMatrix <- GENIE3(exprMatrix, regulators=regulators, targets=targets)
	testthat::expect_equal(nrow(weightMatrix), 3)
	testthat::expect_equal(ncol(weightMatrix), 2)

	#### Targets as name
	targets <- rownames(exprMatrix)[c(5,9)]
	weightMatrix <- GENIE3(exprMatrix, targets=targets)
	testthat::expect_equal(nrow(weightMatrix), 20)
	testthat::expect_equal(ncol(weightMatrix), 2)

	#### Only one target (Note that the regulator is included, and will take all the weight...)
	weightMatrix <- GENIE3(exprMatrix, regulators=2:3, targets=3)
	testthat::expect_equal(ncol(weightMatrix), 1)


	### Other input classes:
	eset <- Biobase::ExpressionSet(assayData=exprMatrix)
	testthat::expect_equal(class(GENIE3(eset)), "matrix")

	sexp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exprMatrix))
	wm <- GENIE3(sexp)
	testthat::expect_equal(class(wm), "matrix")

	# sce <- scater::newSCESet(countData=exprMatrix)
	# testthat::expect_equal(class(GENIE3(sce)), "matrix")


	### Multicore
	set.seed(123)
	weightMatrix_multicore1 <- GENIE3(exprMatrix, nCores=4)
	set.seed(123)
	weightMatrix_multicore2 <- GENIE3(exprMatrix, nCores=4)
	testthat::expect_equal(weightMatrix_multicore1, weightMatrix_multicore2)
})


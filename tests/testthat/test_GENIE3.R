test_that("GENIE3 tests", {

	expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
	rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
	colnames(expr.matrix) <- paste("Sample", 1:5, sep="")

	weight.matrix <- GENIE3(expr.matrix)

	expect_equal(dim(weight.matrix)[1], 20)
	expect_equal(dim(weight.matrix)[2], 20)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(sum(diag(weight.matrix)), 0)
	expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 20)
	expect_equal(sum(rownames(weight.matrix) == rownames(expr.matrix)), 20)


	regulators <- c(2,4,6)
	weight.matrix <- GENIE3(expr.matrix, tree.method="ET", K="all", ntrees=100, regulators=regulators)
	zidx <- setdiff(1:20, regulators)
	# expect_equal(sum(weight.matrix[zidx,]), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(nrow(weight.matrix), 3)
	expect_equal(ncol(weight.matrix), 20)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(diag(weight.matrix[,rownames(weight.matrix)])), 0)
	# expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 20)
	expect_equal(sum(colnames(weight.matrix) == rownames(expr.matrix)), 20)


	regulators <- rownames(expr.matrix)[c(5,9)]
	weight.matrix <- GENIE3(expr.matrix, tree.method="ET", K=1, ntrees=100, regulators=regulators)
	zidx <- setdiff(colnames(weight.matrix), regulators)
	# expect_equal(sum(weight.matrix[zidx,]), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(nrow(weight.matrix), 2)
	expect_equal(ncol(weight.matrix), 20)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(diag(weight.matrix[,rownames(weight.matrix)])), 0)
	# expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 20)
	expect_equal(sum(colnames(weight.matrix) == rownames(expr.matrix)), 20)
})


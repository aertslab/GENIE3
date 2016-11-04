test_that("tlGENIE3 works properly", {
	
	data1 <- matrix(sample(1:10, 70, replace=TRUE), nrow=10)
	rownames(data1) <- paste("Gene", 1:10, sep="")
	colnames(data1) <- paste("Time point", 1:7, sep="")
	
	data2 <- matrix(sample(1:10, 50, replace=TRUE), nrow=10)
	rownames(data2) <- paste("Gene", 1:10, sep="")
	colnames(data2) <- paste("Time point", 1:5, sep="")
	
	TS.data <- list(data1,data2)
	
	SS.data <- matrix(sample(1:10, 100, replace=TRUE), nrow=10)
	rownames(SS.data) <- paste("Gene", 1:10, sep="")
	colnames(SS.data) <- paste("Sample", 1:10, sep="")
	
	
	weight.matrix <- tlGENIE3(TS.data)
	
	expect_equal(dim(weight.matrix)[1], 10)
	expect_equal(dim(weight.matrix)[2], 10)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(sum(diag(weight.matrix)), 0)
	expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 10)
	expect_equal(sum(rownames(weight.matrix) == rownames(TS.data[[1]])), 10)
	
	
	weight.matrix <- tlGENIE3(TS.data,SS.data=NULL)
	
	expect_equal(dim(weight.matrix)[1], 10)
	expect_equal(dim(weight.matrix)[2], 10)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(sum(diag(weight.matrix)), 0)
	expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 10)
	expect_equal(sum(rownames(weight.matrix) == rownames(TS.data[[1]])), 10)
	
	
	weight.matrix <- tlGENIE3(TS.data,SS.data=SS.data)
	
	expect_equal(dim(weight.matrix)[1], 10)
	expect_equal(dim(weight.matrix)[2], 10)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(sum(diag(weight.matrix)), 0)
	expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 10)
	expect_equal(sum(rownames(weight.matrix) == rownames(TS.data[[1]])), 10)
	
	
	regulators <- c(2,4,6)
	weight.matrix <- tlGENIE3(TS.data, h=3, tree.method="ET", K="all", ntrees=100, regulators=regulators)
	zidx <- setdiff(1:10, regulators)
	expect_equal(sum(weight.matrix[zidx,]), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(dim(weight.matrix)[1], 10)
	expect_equal(dim(weight.matrix)[2], 10)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(diag(weight.matrix)), 0)
	expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 10)
	expect_equal(sum(rownames(weight.matrix) == rownames(TS.data[[1]])), 10)
	
	
	regulators <- c("Gene5","Gene9")
	weight.matrix <- tlGENIE3(TS.data, tree.method="ET", K=1, ntrees=100, regulators=regulators)
	zidx <- setdiff(rownames(weight.matrix), regulators)
	expect_equal(sum(weight.matrix[zidx,]), 0)
	expect_equal(sum(weight.matrix < 0), 0)
	expect_equal(dim(weight.matrix)[1], 10)
	expect_equal(dim(weight.matrix)[2], 10)
	expect_true(is.numeric(weight.matrix))
	expect_equal(sum(is.na(weight.matrix)), 0)
	expect_equal(sum(diag(weight.matrix)), 0)
	expect_equal(sum(rownames(weight.matrix) == colnames(weight.matrix)), 10)
	expect_equal(sum(rownames(weight.matrix) == rownames(TS.data[[1]])), 10)

})

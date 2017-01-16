test_that("get.link.list works properly", {
	weight.matrix <- matrix(runif(100), nrow=10)
	rownames(weight.matrix) <- paste("Gene", 1:10, sep="")
	colnames(weight.matrix) <- paste("Gene", 1:10, sep="")
	weight.matrix <- weight.matrix - diag(diag(weight.matrix))

	link.list <- get.link.list(weight.matrix)
	expect_equal(nrow(link.list), 90)
	expect_equal(ncol(link.list), 3)
	expect_true(is.numeric(link.list[,"weight"]))
	expect_equal(sum(is.na(link.list[,"weight"])), 0)
	i <- sample(1:90, size=1)
	expect_equal(weight.matrix[link.list[i,"regulatory.gene"],link.list[i,"target.gene"]], link.list[i,"weight"])
	expect_equal(link.list[1,3], max(weight.matrix))
	diag(weight.matrix) <- NA
	expect_equal(link.list[90,3], min(weight.matrix, na.rm=TRUE))

	diag(weight.matrix) <- 0

	nlinks <- sample(1:90, size=1)
	link.list <- get.link.list(weight.matrix, report.max=nlinks)
	expect_equal(nrow(link.list), nlinks)
	expect_equal(ncol(link.list), 3)
	expect_true(is.numeric(link.list[,"weight"]))
	expect_equal(sum(is.na(link.list[,"weight"])), 0)
	i <- sample(1:nlinks, size=1)
	expect_equal(weight.matrix[link.list[i,"regulatory.gene"],link.list[i,"target.gene"]], link.list[i,"weight"])
	expect_equal(link.list[1,3], max(weight.matrix))

	link.list <- get.link.list(weight.matrix, threshold=0.3)
	nlinks <- nrow(link.list)
	expect_true(link.list[nlinks, 3] >= 0.3)
	expect_equal(ncol(link.list), 3)
	expect_true(is.numeric(link.list[,"weight"]))
	expect_equal(sum(is.na(link.list[,"weight"])), 0)
	i <- sample(1:nlinks, size=1)
	expect_equal(weight.matrix[link.list[i,"regulatory.gene"],link.list[i,"target.gene"]], link.list[i,"weight"])
	expect_equal(link.list[1,3], max(weight.matrix))
})


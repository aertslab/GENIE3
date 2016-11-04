test_that("get.link.list works properly", {
	
	weight.matrix <- matrix(runif(100), nrow=10)
	rownames(weight.matrix) <- paste("Gene", 1:10, sep="")
	colnames(weight.matrix) <- paste("Gene", 1:10, sep="")
	weight.matrix <- weight.matrix - diag(diag(weight.matrix))
	
	link.list <- get.link.list(weight.matrix)
	expect_equal(dim(link.list)[1], 90)
	expect_equal(dim(link.list)[2], 3)
	expect_true(is.numeric(link.list[,3]))
	expect_equal(sum(is.na(link.list[,3])), 0)
	i <- sample(1:90, size=1)
	expect_equal(weight.matrix[link.list[i,1],link.list[i,2]], link.list[i,3])
	expect_equal(link.list[1,3], max(weight.matrix))
	diag(weight.matrix) <- NA
	expect_equal(link.list[90,3], min(weight.matrix, na.rm=TRUE))
	
	diag(weight.matrix) <- 0
	
	nlinks <- sample(1:90, size=1)
	link.list <- get.link.list(weight.matrix, report.max=nlinks)
	expect_equal(dim(link.list)[1], nlinks)
	expect_equal(dim(link.list)[2], 3)
	expect_true(is.numeric(link.list[,3]))
	expect_equal(sum(is.na(link.list[,3])), 0)
	i <- sample(1:nlinks, size=1)
	expect_equal(weight.matrix[link.list[i,1],link.list[i,2]], link.list[i,3])
	expect_equal(link.list[1,3], max(weight.matrix))
	
	link.list <- get.link.list(weight.matrix, threshold=0.3)
	nlinks <- dim(link.list)[1]
	expect_true(link.list[nlinks, 3] >= 0.3)
	expect_equal(dim(link.list)[2], 3)
	expect_true(is.numeric(link.list[,3]))
	expect_equal(sum(is.na(link.list[,3])), 0)
	i <- sample(1:nlinks, size=1)
	expect_equal(weight.matrix[link.list[i,1],link.list[i,2]], link.list[i,3])
	expect_equal(link.list[1,3], max(weight.matrix))
})


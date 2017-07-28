test_that("getLinkList works properly", {
  
	weightMatrix <- matrix(runif(100), nrow=10)
	rownames(weightMatrix) <- paste("Gene", 1:10, sep="")
	colnames(weightMatrix) <- paste("Gene", 1:10, sep="")
	weightMatrix <- weightMatrix - diag(diag(weightMatrix))

	linkList <- getLinkList(weightMatrix)
	testthat::expect_equal(nrow(linkList), 90)
	testthat::expect_equal(ncol(linkList), 3)
	testthat::expect_true(is.numeric(linkList[,"weight"]))
	testthat::expect_equal(sum(is.na(linkList[,"weight"])), 0)
	i <- sample(1:90, size=1)
	testthat::expect_equal(weightMatrix[linkList[i,"regulatoryGene"],linkList[i,"targetGene"]], linkList[i,"weight"])
	testthat::expect_equal(linkList[1,3], max(weightMatrix))
	diag(weightMatrix) <- NA
	testthat::expect_equal(linkList[90,3], min(weightMatrix, na.rm=TRUE))

	diag(weightMatrix) <- 0

	nlinks <- sample(1:90, size=1)
	linkList <- getLinkList(weightMatrix, reportMax=nlinks)
	testthat::expect_equal(nrow(linkList), nlinks)
	testthat::expect_equal(ncol(linkList), 3)
	testthat::expect_true(is.numeric(linkList[,"weight"]))
	testthat::expect_equal(sum(is.na(linkList[,"weight"])), 0)
	i <- sample(1:nlinks, size=1)
	testthat::expect_equal(weightMatrix[linkList[i,"regulatoryGene"],linkList[i,"targetGene"]], linkList[i,"weight"])
	testthat::expect_equal(linkList[1,3], max(weightMatrix))

	linkList <- getLinkList(weightMatrix, threshold=0.3)
	nlinks <- nrow(linkList)
	testthat::expect_true(linkList[nlinks, 3] >= 0.3)
	testthat::expect_equal(ncol(linkList), 3)
	testthat::expect_true(is.numeric(linkList[,"weight"]))
	testthat::expect_equal(sum(is.na(linkList[,"weight"])), 0)
	i <- sample(1:nlinks, size=1)
	testthat::expect_equal(weightMatrix[linkList[i,"regulatoryGene"],linkList[i,"targetGene"]], linkList[i,"weight"])
	testthat::expect_equal(linkList[1,3], max(weightMatrix))
})


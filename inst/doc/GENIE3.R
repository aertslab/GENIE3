## ------------------------------------------------------------------------
expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
colnames(expr.matrix) <- paste("Sample", 1:5, sep="")
head(expr.matrix)

## ------------------------------------------------------------------------
library(GENIE3)
weight.matrix <- GENIE3(expr.matrix)

## ------------------------------------------------------------------------
dim(weight.matrix)
weight.matrix[1:5,1:5]

## ------------------------------------------------------------------------
# Genes that are used as candidate regulators
regulators <- c(2, 4, 7)
# Or alternatively:
regulators <- c("Gene2", "Gene4", "Gene7")
weight.matrix <- GENIE3(expr.matrix, regulators=regulators)

## ----eval=FALSE----------------------------------------------------------
#  # Use Extra-Trees method
#  tree.method = "ET"
#  
#  # Number of randomly chosen candidate regulators at each node of a tree
#  K = 7
#  
#  # Number of trees per ensemble
#  ntrees = 50
#  
#  # Run the method with these settings
#  weight.matrix = GENIE3(expr.matrix, tree.method=tree.method, K=K, ntrees=ntrees)

## ----eval=FALSE----------------------------------------------------------
#  weight.matrix <- GENIE3(expr.matrix, ncores=4)

## ------------------------------------------------------------------------
?GENIE3

## ------------------------------------------------------------------------
link.list <- get.link.list(weight.matrix)
dim(link.list)
head(link.list)

## ----eval=FALSE----------------------------------------------------------
#  link.list <- get.link.list(weight.matrix, report.max=5)

## ----eval=FALSE----------------------------------------------------------
#  link.list <- get.link.list(weight.matrix, threshold=0.1)


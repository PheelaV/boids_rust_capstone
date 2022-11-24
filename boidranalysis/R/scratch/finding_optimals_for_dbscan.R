## Example 1: use dbscan on the iris data set
data(iris)
iris <- as.matrix(iris[, 1:4])
## Find suitable DBSCAN parameters:
## 1. We use minPts = dim + 1 = 5 for iris. A larger value can also be used.
## 2. We inspect the k-NN distance plot for k = minPts - 1 = 4
kNNdistplot(x, minPts = 3)
## Noise seems to start around a 4-NN distance of .7
abline(h=100, col = "red", lty = 2)
## Cluster with the chosen parameters
res <- dbscan(x, eps = 100, minPts = 3)
res
pairs(x, col = res$cluster + 1L)
## Use a precomputed frNN object
fr <- frNN(x, eps = 100)
dbscan(fr, minPts = 3)

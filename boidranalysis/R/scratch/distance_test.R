# define a probability density function P
P <- 1:10/sum(1:10)
# define a probability density function Q
Q <- 20:29/sum(20:29)

# combine P and Q as matrix object
x <- rbind(P,Q)


# compute the Euclidean Distance with default parameters
distance(x, method = "euclidean")


r



newTest <- test %>% mutate(
  x = x / sum(x),
  y = y / sum(y)
)

distance(x)
get_dist_mat()

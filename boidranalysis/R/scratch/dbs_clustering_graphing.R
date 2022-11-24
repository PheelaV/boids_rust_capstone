# install.packages("dbscan")
library(dbscan)
# set.seed(665544)
# n <- 100
# x <- cbind(
#   x=runif(10, 0, 5) + rnorm(n, sd = 0.4),
#   y=runif(10, 0, 5) + rnorm(n, sd = 0.4)
# )

tic("general")
no_iter <- 128000
sample_rate <- 32

time_steps <- no_iter / sample_rate

tic("data generation")
boid_data <- flock_return(no_iter = no_iter, init_boids = 515, save_locations_path = data_folder, sample_rate = sample_rate, init_width = 20000, init_height = 20000)
toc()

tic("processing")
boid_data <- boid_data %>% relocate(id, .after = y)
boid_data <- boid_data %>% relocate(x, .after = y)
ordered_boid_data <- order_by_time(boid_data)
toc()

# nt_boid_data <- nth_time_point(ordered_boid_data, time_steps)
# nt_boid_data_2 <- nt_boid_data
# nt_boid_data <- nth_time_point(ordered_boid_data, 1600)
# x <- cbind(
#   x=nt_boid_data$x,
#   y=nt_boid_data$y
# )


tic("flock clustering")
all <- tibble()
for (i in 1:time_steps) {
  nt_boid_data <- nth_time_point(ordered_boid_data, i)
  x <- cbind(
    x=nt_boid_data$x,
    y=nt_boid_data$y
  )
  dbs_res <- dbscan(x, eps = 200, minPts = 3)
  nt_boid_data$cluster = dbs_res$cluster
  all <- all %>% bind_rows(nt_boid_data)
}
toc()

toc()
nth_time_point(all, 1)

library(gganimate)


dbs_res <- dbscan(x, eps = 200, minPts = 2)
dbs_res <- dbscan(nt_boid_data, eps = 200, minPts = 2)
nt_boid_data$cluster = dbs_res$cluster
hullplot(x, dbs_res)
hullplot(tibble(nt_boid_data$x, nt_boid_data$y), dbs_res)
ggplot(nt_boid_data_2, aes(x = x, y = y, col = cluster), scale) + geom_point() + lims(x = c(-10000, 10000), y = c(-10000, 10000))
ggplot(nt_boid_data, aes(x = x, y = y, col = cluster)) + geom_point() + lims(x = c(-10000, 10000), y = c(-10000, 10000))

all <- nt_boid_data %>% bind_rows(nt_boid_data_2)

ggplot(all, aes(x = x, y = y, col = cluster)) + geom_point() + lims(x = c(-10000, 10000), y = c(-10000, 10000))

plot(x, col = dbs_res$cluster)
points(x[dbs_res$cluster == 0, ], pch = 3, col = "grey")

plot(x, pch = 16)
# Connected components on a graph where each pair of points
# with a distance less or equal to eps are connected
d <- dist(x)
components <- comps(d, eps = 220)
plot(x, col = components, pch = 16)
# Connected components in a fixed radius nearest neighbor graph
# Gives the same result as the threshold on the distances above
frnn <- frNN(x, eps = 220)
components <- comps(frnn)
plot(frnn, data = x, col = components)
# Connected components on a k nearest neighbors graph
knn <- kNN(x, 3)
components <- comps(knn, mutual = FALSE)
plot(knn, data = x, col = components)
components <- comps(knn, mutual = TRUE)
plot(knn, data = x, col = components)
# Connected components in a shared nearest neighbor graph
snn <- sNN(x, k = 10, kt = 5)
components <- comps(snn)
plot(snn, data = x, col = components)







# https://cran.r-project.org/web/packages/mstknnclust/vignettes/guide.html
# https://stackoverflow.com/questions/30647537/using-geo-coordinates-as-vertex-coordinates-in-the-igraph-r-package

# library("mstknnclust")
# library(tidyverse)
# boid_data <- flock_return(no_iter = 1000, init_boids = 256, save_locations_path = data_folder, sample_rate = 8, init_width = 4000, init_height = 4000)
ordered_boid_data <- order_by_time(boid_data)
nth_time_point(ordered_boid_data, 100)
test <- tibble(
  id = c(1, 2, 3, 4, 5, 6),
  x = c(1, 2, 3, 18, 19, 20),
  y = c(1, 2, 3, 6, 5, 6)
)

# distinct_boid_no <- nrow(boid_data %>% dplyr::select(id) %>% unique())
#

order_by_time <- function(data) {
  distinct_no <- nrow(data %>% dplyr::select(id) %>% unique())
  ordered_boid_data <- data %>%
    group_by(id) %>%
    mutate(
      order = id + (row_number(id) - 1) * distinct_no,
      timestep = row_number(id)
      ) %>%
    ungroup() %>%
    arrange(order)
}

nth_time_point <- function(data, time_n) {
  distinct_no <- nrow(data %>% dplyr::select(id) %>% unique())
  l <- nrow(data)

  tail <- l - (distinct_no * (time_n - 1))
  head <- distinct_no

  data %>%
    slice_tail(n = tail) %>%
    slice_head(n = head)
}

nth_time_point(ordered_boid_data, 1)

ordered_boid_data <- boid_data %>%
  group_by(id) %>%
  mutate(order = id + (row_number(id) - 1) * distinct_boid_no) %>%
  ungroup() %>%
  arrange(order)
#
# test <- ordered_boid_data %>% slice_tail(n = 100)
test <- ordered_boid_data %>% slice_head(n = 100)

euclidean <- function(a, b) sqrt(sum((a - b)^2))

# take a look at https://cran.r-project.org/web/packages/philentropy/index.html


get_dist_mat <- function(df) {
  l <- nrow(df)
  r <- matrix(0, nrow = l, ncol = l)
  colnames(r) <- df$id
  rownames(r) <- df$id
  for (i in seq_len(l)) {
    for (j in seq_len(l)) {
      r[i, j] <- euclidean(
        c(df$x[i], df$y[i]),
        c(df$x[j], df$y[j])
      )
    }
  }

  r
}

r <- get_dist_mat(test)
results <- mstknnclust::mst.knn(r)

library(philentropy)
# results

# library("igraph")


# igraph::V(results$network)$label.cex <- seq(0.6,0.6,length.out=2)

plot(results$network, vertex.size=8,
     vertex.color=igraph::clusters(results$network)$membership,
     layout=igraph::layout.fruchterman.reingold(results$network, niter=10000),
     # layout=igraph::layout.norm(location),
     layout=igraph::layout.norm(test[c("x", "y")] %>% as.matrix()),
     main=paste("MST-kNN \n Clustering solution \n Number of clusters=",results$cnumber,sep="" ))



# location <- test[c("x", "y")] %>% as.matrix()


location = tibble(id = as.numeric(V(results$network)$name)) %>%
  right_join(test) %>%
  select(x, y) %>%
  as.matrix()


plot(results$network,
     vertex.color=c('red', 'green', 'blue'),
     vertex.label = NA,
     vertex.size = , edge.width = 0.5,
     layout = location,
     rescale=FALSE,
     # asp=0,
     xlim = range(test$x), ylim = range(test$y))



#
# ggplot(test, aes(x, y, label = id)) +
#   geom_point() +
#   geom_text(hjust = 0, vjust = 0)

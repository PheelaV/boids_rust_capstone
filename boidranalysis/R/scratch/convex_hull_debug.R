data <- boid_data


no_time_points <- n_distinct(data$time)

t_data <- nth_time_point(data, 1)
no_flocks <- n_distinct(t_data$cluster_id) - 1

clusters <- data %>% filter(time == 1) %>%
  distinct(cluster_id) %>%
  filter(cluster_id != 0)

volume_t_clusters <- numeric(nrow(clusters))


for (i in 1:length(volume_t_clusters)) {
  c <- clusters[i,]$cluster_id
  cluster_hull <- cxhull(
    as.matrix(
      t_data %>%
        filter(cluster_id == c) %>%
        select(x, y)
    )
  )

  volume_t_clusters[i] <- cluster_hull$volume
}
rm(c)

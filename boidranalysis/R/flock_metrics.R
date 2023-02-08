# gets the average size of flocks for each time point
get_convex_hull_data <- function(data) {
  data %>%
    filter(cluster_id != 0) %>%
    group_by(time, cluster_id) %>%
    filter(n() > 2) %>%
    summarise(hull_area = get_convex_hull(x, y)) %>%
    summarise(average_cluster_volumes = mean(hull_area))
}

# get the number of flocks for each time point
# takes in raw or pre-processed boid data with time and cluster_id
# returns a vector with number of distinct time points as its length
get_no_flocks <- function(data) {
  # for each time point, count the number of flocks, we are treating 0 (set of
  # entities not belonging to any flock
  data %>%
    filter(cluster_id != 0) %>%
    group_by(time) %>%
    distinct(cluster_id) %>%
    summarise(no_flocks = n())
}

# relies on the presence of dx, dy, time, cluster_id
get_average_norm_vel <- function(data) {
  data %>%
    filter(cluster_id != 0) %>%
    select(dx, dy, time, cluster_id) %>%
    group_by(time) %>%
    # group_by(time, cluster_id) %>%
    summarise(
      v0x = sum(abs(dx)) / n(),
      v0y = sum(abs(dy)) / n(),
      vix = abs(sum(dx)),
      viy = abs(sum(dy)),
      count = n(),
      no_flocks = n_distinct(cluster_id)
    ) %>%
    mutate(avg_norm_vel = sqrt(vix^2 + viy^2) / (count * sqrt(v0x^2 + v0y^2))) %>%
    select(time, avg_norm_vel)
}

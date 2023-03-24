library(cxhull)
# area of convex hull
get_convex_hull_data_clean <- function(data) {
  # the 0 cluster is "noise", i.e. all entities not in a distinct cluster
  data <- data %>% filter(cluster_id != 0)

  t_points <- data %>%
    select(time) %>%
    arrange(time) %>%
    distinct()

  # get the number of time points
  no_time_points <- nrow(t_points)

  average_cluster_volume_t <- numeric(no_time_points)

  for (j in 1:no_time_points){
    # t_data <- nth_time_point(data, t) %>% filter(cluster_id != 0)

    t_data <- data %>%
      filter(time == t_points$time[j]) %>%
      select(x, y, cluster_id)
    no_flocks <- n_distinct(t_data$cluster_id) - 1

    clusters <- t_data %>%
      group_by(cluster_id) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      select(cluster_id)

    volume_t_clusters <- numeric(nrow(clusters))

    for (i in 1:length(volume_t_clusters)) {
      c_i <- clusters[i,]$cluster_id
      cluster_hull <- cxhull(
        as.matrix(
          t_data %>%
            filter(cluster_id == c_i) %>%
            select(x, y)
        )
      )

      volume_t_clusters[i] <- cluster_hull$volume
      # volume_t_clusters[i] <- det(var(as.matrix(
      #   t_data %>%
      #     filter(cluster_id == c_i) %>%
      #     select(x, y)
      # )))
    }

    average_cluster_volume_t[j] <- mean(volume_t_clusters)
  }
  return(tibble(average_cluster_volumes = average_cluster_volume_t))
}

# select <- dplyr::select

get_convex_hull_data_noisy <- function(data) {
  t_points <- data %>%
    select(time) %>%
    arrange(time) %>%
    distinct()

  # get the number of time points
  no_time_points <- nrow(t_points)

  t_average_cluster_volume <- numeric(no_time_points)

  for (j in 1:no_time_points) {
    t_data <- data %>%
      filter(time == t_points$time[j], cluster_id != 0) %>%
      select(x, y, cluster_id)

    cluster_ids <- t_data %>%
      group_by(cluster_id) %>%
      summarise(count = n()) %>%
      filter(count > 2) %>%
      select(cluster_id)

      t_cluster_volumes <- numeric(nrow(cluster_ids))

    for (i in 1:nrow(cluster_ids)) {
      # print(cluster_ids$cluster_id[i])
      c_data <- t_data %>% filter(cluster_id == cluster_ids$cluster_id[i])

      # if (nrow(c_data) < 2) {
      #   t_cluster_volumes[i] <- 0
      # } else {
        t_cluster_volumes[i] <- get_convex_hull(c_data$x, c_data$y)
      # }
    }

    t_average_cluster_volume[j] <- mean(t_cluster_volumes)
  }

  return(tibble(average_cluster_volumes = t_average_cluster_volume))

}

get_convex_hull_data_3 <- function(data) {
  data %>%
    filter(cluster_id != 0) %>%
    group_by(time, cluster_id) %>%
    filter(n() > 2) %>%
    summarise(hull_area = get_convex_hull(x, y)) %>%
    summarise(average_cluster_volumes = mean(hull_area))
}

get_convex_hull_data_4 <- function(data) {
  data %>%
    filter(cluster_id != 0) %>%
    group_by(time, cluster_id) %>%
    filter(n() > 2) %>%
    summarise(hull_area = cxhull(
      cbind(x, y)
    )$volume) %>%
    summarise(average_cluster_volumes = mean(hull_area))
}

test <- data %>%
  filter(cluster_id != 0) %>%
  select(time, cluster_id, x, y) %>%
  group_by(time, cluster_id) %>%
  filter(n() > 2) %>%
  group_map(
    ~cxhull(cbind(.x$x, .x$y))$volume
  ) %>%
  rapply(function(row)row[[1]], how="unlist")

cluster_hull <- cxhull(
  as.matrix(
    x, y
  )
)

cxhull(
  as.matrix(
    x, y
  )
)$volume

volume_t_clusters[i] <- cluster_hull$volume

data <- direction_data_by_boid

filtered <- data %>%
  filter(cluster_id != 0)
  group_by(time, cluster_id) %>%
  filter(n() > 2) %>%
  ungroup()


tic("method 1")
test1 <- get_convex_hull_data_clean(data)
toc()

tic("method 2")
test2 <- get_convex_hull_data_noisy(data)
toc()

tic("method 3")
test3 <- get_convex_hull_data_3(data)
test <- toc()
#
tic("method 4")
test4 <- get_convex_hull_data_4(data)
toc()

ggplot(
  tibble(
  t1 = test1$average_cluster_volumes,
  t2 = test2$average_cluster_volumes,
  t3 = test3$average_cluster_volumes,
  t4 = test4$average_cluster_volumes
  ) %>% mutate(t = row_number())
) +
  geom_line(aes(x = t, y = t3, color = "new1"), linewidth = 1) +
  geom_line(aes(x = t, y = t1, color = "old1"), linetype = "dashed", linewidth = 1) +
  geom_line(aes(x = t, y = t2, color = "old2"), linewidth = 1.5, linetype = "dotted") +
  geom_line(aes(x = t, y = t4, color = "new2"), linewidth = .8, linetype = "dashed")

if (FALSE) {
  boid_data_512 %>%
    filter(time == 900, cluster_id == 2) %>%
    select(x, y)


  library(cxhull)

  cluster_hull <- cxhull(
    boid_data_512 %>%
      filter(time == 900, cluster_id == 2) %>%
      select(x, y) %>%
      # mutate(point = c(x, y)) %>%
      mutate(point = map2(x, y, ~c(.x, .y))) %>%
      select(point)
  )

  boid_data %>%
    filter(time == 900, cluster_id == 2) %>%
    select(x, y) %>%
    # mutate(point = c(x, y)) %>%
    mutate(point = map2(x, y, ~c(.x, .y))) %>%
    select(point)

  hull <- cxhull(as.matrix(flock_at_t))

  boid_data %>% filter(time == 50) %>% distinct(cluster_id)

  run_experiment <- function() {
    return(
      flock_detailed(no_iter = 2^14,
                     init_boids = 2^9,
                     save_locations_path = data_folder, # ,
                     sample_rate = 2^5,
                     init_width = 4000,
                     init_height = 4000,
                     sensory_distance = 50,
                     alignment_coef = .02,
                     alignment_trs_coef = 1.15,
                     cohesion_coef = 0.002,
                     cohesion_trs_coef = .95,
                     separation_coef = 4.1,
                     separation_trs_coef = .3,
                     min_speed = .5,
                     max_speed = 2.,
                     max_steering = .65,
                     dbscan_clustering = T,
                     boundary_config = "{\"type\": \"Toroidal\"}") %>%
        tibble()
    )
  }

  tic("run experiment")
  run_experiment()
  toc()

  run_experiment_data_collection_convex <- function(i) {
    return(
      run_experiment() %>% get_convex_hull_data_3() %>%
        mutate(label =  sprintf("t_%d", time)) %>%
        select(label, average_cluster_volumes) %>%
        pivot_wider(values_from = average_cluster_volumes, names_from = label)
    )
  }

  library(tictoc)
  library(parallel)
  library(tidyverse)

  no_cores <- 8
  no_runs <- 1 * no_cores
  tic("total")
  tic("multi core long") # this took 5851.272 sec, 1h 37m
  # now running with half the iterration to get the results in time for capstone supervision
  tic("no flocks")
  results_no_flocks <- mclapply(1:no_runs, run_experiment_data_collection, mc.cores = no_cores)
  toc()

  tic("convex hull")
  results_convex_hull <- mclapply(1:no_runs, run_experiment_data_collection_convex, mc.cores = no_cores)
  toc()

  toc()

  results_no_flocks <- map_dfr(results_no_flocks, bind_rows)
  results_convex_hull <- map_dfr(results_convex_hull, bind_rows)
  ## nooo, there are some errors in the convex hull runs
  # filtered_results_convex_hull <- Filter(function(x) class(x)[1] == "tbl_df", results_convex_hull)


  ncol(results_no_flocks)
  ncol(results_convex_hull)

  results_stats <- tibble(
    t = 1:ncol(results_no_flocks),
    mean_no_flocks = apply(results_no_flocks, 2, mean),
    var_no_flocks = apply(results_no_flocks, 2, var),
    mean_convex_hull = apply(results_convex_hull, 2, mean),
    var_convex_hull = apply(results_convex_hull, 2, var),
  )

  s <- 17
  e <- 32
  results_stats <- tibble(
    t = 1:ncol(results_no_flocks),
    mean_no_flocks = apply(results_no_flocks[s:e,], 2, mean),
    var_no_flocks = apply(results_no_flocks[s:e,], 2, var),
    mean_convex_hull = apply(results_convex_hull[s:e,], 2, mean),
    var_convex_hull = apply(results_convex_hull[s:e,], 2, var),
  )

  # in summary, this is run over 32K iterations, sampling each
  # 32nd step = 1000 data points per run both number of flocks
  # and area of a flock have been run 128 times in total, area
  # of flocks yielded only  # 48 results as there is still some
  # debugging to be done

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_no_flocks)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time, t_x", y = "mean # flocks", title = "mean number of flocks")

  ggplot(
    data = results_stats,
    aes(x = t, y = var_no_flocks)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time, t_x", y = "var # flocks", title = "var number of flocks")

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_convex_hull)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time, t_x", y = "mean flock area", title = "mean flock area (filtered)")

  ggplot(
    data = results_stats,
    aes(x = t, y = var_convex_hull)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time, t_x", y = "var flock area", title = "var flock area (filtered)")
  toc()
}

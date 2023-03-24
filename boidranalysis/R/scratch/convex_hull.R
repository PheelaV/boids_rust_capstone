library(cxhull)
# area of convex hull
get_convex_hull_desc <- function(data) {
    # get the number of time points
    no_time_points <- n_distinct(data$time)

    average_cluster_volume_t <- numeric(no_time_points)

    for (t in 1:no_time_points){
      t_data <- nth_time_point(data, t)
      no_flocks <- n_distinct(t_data$cluster_id) - 1

      clusters <- data %>% filter(time == t) %>% distinct(cluster_id)

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

      average_cluster_volume_t[t] <- mean(volume_t_clusters)
    }
    return(tibble(average_cluster_volumes = average_cluster_volume_t))
}
test <- get_convex_hull_desc(boid_data)
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
      flock_detailed(no_iter = 16000,
                     init_boids = 512,
                     save_locations_path = data_folder, # ,
                     sample_rate = 32,
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
                     dbscan_clustering = T) %>%
        tibble()
    )
  }

  run_experiment_data_collection_convex <- function(i) {
    return(
      run_experiment() %>% get_convex_hull_desc() %>%
        mutate(label =  sprintf("t_%d", row_number())) %>%
        pivot_wider(values_from = average_cluster_volumes, names_from = label)
    )
  }

  library(tictoc)
  library(parallel)
  library(tidyverse)

  tic("total")
  tic("multi core long") # this took 5851.272 sec, 1h 37m
  # now running with half the iterration to get the results in time for capstone supervision
  results_no_flocks <- mclapply(1:128, run_experiment_data_collection, mc.cores = 8)
  results_convex_hull <- mclapply(1:128, run_experiment_data_collection_convex, mc.cores = 8)
  toc()

  results_no_flocks <- map_dfr(results_no_flocks, bind_rows)
  results_convex_hull <- map_dfr(results_convex_hull, bind_rows)
  ## nooo, there are some errors in the convex hull runs
  filtered_results_convex_hull <- Filter(function(x) class(x)[1] == "tbl_df", results_convex_hull)


  results_stats <- tibble(
    t = 1:ncol(results_no_flocks),
    mean_no_flocks = apply(results_no_flocks, 2, mean),
    var_no_flocks = apply(results_no_flocks, 2, var),
    mean_convex_hull = apply(filtered_results_convex_hull, 2, mean),
    var_convex_hull = apply(filtered_results_convex_hull, 2, var),
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

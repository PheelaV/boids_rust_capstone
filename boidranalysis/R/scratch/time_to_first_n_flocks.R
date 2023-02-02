# get number of flocks at time t
.t_no_flocks <- function(data, t) {
  t_data <- nth_time_point(data, t)
  no_flocks <- n_distinct(t_data$cluster_id) - 1
  return(no_flocks)
  # get the number of time points
  no_time_points <- n_distinct(data$time)
}

# get the time to first number of flocks <= n
first_n_flocks_time <- function(data, n) {
  # get the number of time points
  no_time_points <- n_distinct(data$time)

  # for each time point, count the number of flocks, we are treating 0 (set of
  # entities not belonging to any flock, as a distinctive flock by itself and
  # thus subtracting 1 from the final result)

  for (t in 1:no_time_points) {
    t_data <- nth_time_point(data, t)
    no_flocks <- n_distinct(t_data$cluster_id) - 1
    if (no_flocks <= n) {
      return(t)
    }
  }
}

# get the number of flocks for each time point
# takes in raw boid data
# returns a vector wwith number of distinct timepoints as its length
get_no_flocks <- function(data) {
  # get the number of time points
  no_time_points <- n_distinct(data$time)

  # for each time point, count the number of flocks, we are treating 0 (set of
  # entities not belonging to any flock, as a distinctive flock by itself and
  # thus subtracting 1 from the final result)

  results <- numeric(no_time_points)
  for (t in 1:no_time_points) {
    no_flocks <- .t_no_flocks(data, t)
    results[t] <- no_flocks
  }

  return(tibble(no_flocks = results))
}

# usages

if(FALSE) {
  run_experiment <- function() {
    return(
      flock_detailed(no_iter = 16000,
                     init_boids = 512,
                     save_locations_path = data_folder, # ,
                     sample_rate = 16,
                     init_width = 4000,
                     init_height = 4000,
                     sensory_distance = 50,
                     allignment_coef = .02,
                     allignment_trs_coef = 1.15,
                     cohesion_coef = 0.002,
                     cohesion_trs_coef = .95,
                     separation_coef = 4.1,
                     separation_trs_coef = .3,
                     min_speed = .5,
                     max_speed = 2.,
                     max_steering = .65,
                     dbscan_clustering = T)
    )
  }

  run_experiment_data_collection <- function(i) {
    return(
      run_experiment() %>% get_no_flocks() %>%
        mutate(label =  sprintf("t_%d", row_number())) %>%
        pivot_wider(values_from = no_flocks, names_from = label)
    )
  }

  boid_data_512 <- flock_detailed(no_iter = 128000,
                              init_boids = 512,
                              save_locations_path = "", # data_folder,
                              sample_rate = 32,
                              init_width = 4000,
                              init_height = 4000,
                              sensory_distance = 50,
                              allignment_coef = .02,
                              allignment_trs_coef = 1.15,
                              cohesion_coef = 0.002,
                              cohesion_trs_coef = .95,
                              separation_coef = 4.1,
                              separation_trs_coef = .3,
                              min_speed = .5,
                              max_speed = 2.,
                              max_steering = .65,
                              dbscan_clustering = T)

  boid_data_2048 <- flock_detailed(no_iter = 32000,
                              init_boids = 2048,
                              save_locations_path = "", # data_folder,
                              sample_rate = 32,
                              init_width = 4000,
                              init_height = 4000,
                              sensory_distance = 50,
                              allignment_coef = .02,
                              allignment_trs_coef = 1.15,
                              cohesion_coef = 0.002,
                              cohesion_trs_coef = .95,
                              separation_coef = 4.1,
                              separation_trs_coef = .3,
                              min_speed = .5,
                              max_speed = 2.,
                              max_steering = .65,
                              dbscan_clustering = T)


  boid_data <- boid_data_512

  n_distinct(boid_data$time)

  range(boid_data$time)

  boid_data_no_flocks <- boid_data %>% get_no_flocks()

  ggplot(data = boid_data_no_flocks %>% mutate(time = row_number()), aes(x = time, y = no_flocks)) +
    geom_point() +
    theme_bw()

  # calculate on the fly
  boid_data %>% first_n_flocks_time(1)
  first_n_flocks_time(boid_data, 1)

  # or get results from a cached version
  as.numeric((boid_data_no_flocks %>%
    mutate(time = row_number()) %>%
    filter(no_flocks == 3) %>%
    slice_head(n = 1))[1,'time'])
  # when ran with 512 boids, the tt1f = 939
  # for 2048 it was 444

  sprintf("t%d", seq(
    from = 1,
    to = 5,
    by = 1
  ))

  sprintf("t_%d", 1)

  boid_data_no_flocks_512 <- boid_data_512 %>% get_no_flocks()
  boid_data_no_flocks_2048 <- boid_data_2048 %>% get_no_flocks()

  test <- tibble()
  test <- boid_data_no_flocks_512 %>%
    mutate(label =  sprintf("t_%d", row_number())) %>%
    pivot_wider(values_from = no_flocks, names_from = label)

  test <- test %>% bind_rows(
    boid_data_no_flocks_2048 %>%
        mutate(label =  sprintf("t_%d", row_number())) %>%
        pivot_wider(values_from = no_flocks, names_from = label)
  )

  library(tictoc)
  library(parallel)
  library(tidyverse)

  tic("single core") # 83 seconds, most of it is the R bit :<
  results <- tibble()
  for (i in 1:8) {
    results <- results %>% bind_rows(
      run_experiment_data_collection(0)
    )
  }
  # results
  toc()

  detectCores()
  tic("multi core")
  results2 <- mclapply(1:8, run_experiment_data_collection, mc.cores = 4)
  toc()

  results2 <- map_dfr(results2, bind_rows)

  tic("multi core long") # this took 5851.272 sec, 1h 37m
  # now running with half the iterration to get the results in time for capstone supervision
  results3 <- mclapply(1:500, run_experiment_data_collection, mc.cores = 6)
  toc()

  results3 <- map_dfr(results3, bind_rows)
saveRDS(results3, "./results3.rds")

ggplot(data = boid_data_no_flocks %>% mutate(time = row_number()), aes(x = time, y = no_flocks)) +
  scale_y_log10() +
  geom_line() +
  theme_bw()

print("max mean")
print(max(apply(results3, 2, mean)))

########### plotting ###########

results_stats <- tibble(
  t = 1:ncol(results3),
  mean = apply(results3, 2, mean),
  var = apply(results3, 2, var),
)

ggplot(
  data = results_stats,
  aes(x = t, y = mean)
) +
  geom_point() +
  theme_bw() +
  # geom_function(fun = dt, color = "red", args = list(df = 20))
  labs(x = "time, t_x", y = "mean # flocks", title = "1000 simulations, 512 boids, 16K iter, 16sr")

ggplot(
  data = results_stats,
  aes(x = t, y = var)
) +
  geom_point() +
  theme_bw() + labs(x = "time, t_x", y = "var # flocks", title = "1000 simulations, 512 boids, 16K iter, 16sr")

library(circular)

?tapply
par}



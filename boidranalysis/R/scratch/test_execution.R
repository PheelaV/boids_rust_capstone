#!/usr/local/bin/Rscript

source("R/helpers.R")
library(RcppTOML)
library(rextendr)
library(tidyverse)
library(dplyr)
library(tictoc)

{
  wd <- getwd()
  setwd("~/Source/Repos/boids_rust/boidr")
  rextendr::document()
  devtools::load_all(".")
  setwd(wd)
}

{
  tic("10K iter from R:") # this took over 700 seconds
  test_config_2 <- get_config(
    "3noisy_c.toml",
    overwrite = list(
      no_iter = 10000,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}"
    )
  )

  test_boid_data_2 <- rlang::exec(flock_detailed, !!!test_config_2) %>%
    tibble()
  toc()
}

{
  tic("10K iter from R without dbscan:") # this took over 130 seconds
  test_config <- get_config(
    "3noisy_c.toml",
    overwrite = list(
      no_iter = 10000,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}",
      dbscan_clustering = FALSE
    )
  )

  test_boid_data <- rlang::exec(flock_detailed, !!!test_config) %>%
    tibble()
  toc()

}

{
  tic("test execution, 10K without robj return:") # this took 42 seconds
  print(test_execution())
  toc()

}

{
  test_config <- get_config(
    "3noisy_c.toml",
    overwrite = list(
      no_iter = 10000,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}",
      dbscan_clustering = FALSE
    )
  )
  tic("test execution, 10K without robj return:") # 89.263 sec elapsed
  rlang::exec(flock_detailed_no_return, !!!test_config)
  toc()

}


{
  test_config <- get_config(
    "3noisy_c.toml",
    overwrite = list(
      no_iter = 10000,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}",
      dbscan_clustering = TRUE
    )
  )
  tic("test execution, 10K without robj return:") # this took 618.523 sec elapsed
  rlang::exec(flock_detailed_no_return, !!!test_config)
  toc()
}



# 35.727 sec elapsed with steering
# 11.712 sec elapsed without steering

{
  good_config <- get_config(
    "1normal.toml",
    overwrite = list(
      no_iter = 10000,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}",
      dbscan_clustering = FALSE
    )
  )
  # test_config <- get_config
  test_config <- get_config(
    "3noisy_c.toml",
    overwrite = list(
      no_iter = 10000,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}",
      dbscan_clustering = FALSE
    )
  )
  tic("comparing two configs") #
  flock_detailed_no_return(test_config$no_iter,
                           test_config$init_boids,
                           test_config$save_locations_path,
                           test_config$sample_rate,
                           test_config$init_width,
                           test_config$init_height,
                           test_config$sensory_distance,
                           good_config$alignment_coef,
                           good_config$cohesion_coef,
                           good_config$separation_coef,
                           good_config$alignment_trs_coef,
                           good_config$cohesion_trs_coef,
                           good_config$separation_trs_coef,
                           good_config$min_speed,
                           good_config$agent_steering,
                           good_config$max_speed,
                           good_config$max_steering,
                           F, #good_config$dbscan_clustering,
                           test_config$boundary_config,
                           test_config$distance_config,
                           test_config$field_of_vision,
                           test_config$rules_impl,
                           test_config$wander_on,
                           test_config$wander_coef,
                           test_config$wander_rate,
                           test_config$wander_radius,
                           test_config$wander_distance,
                           good_config$baseline_speed,
                           test_config$wander_random)
  toc()
}

{
  test_config <- get_config(
    # "22normal_s.toml",
    "3noisy_c.toml",
    overwrite = list(
      no_iter = 2^13,
      init_width = 1000,
      init_height = 1000,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}"
    )
  )
  test_config$dbscan_clustering = FALSE;
  tic("test execution, experiment_m timing:") # 41.954 sec elapsed
  rlang::exec(flock_detailed_no_return, !!!test_config)
  toc()
}


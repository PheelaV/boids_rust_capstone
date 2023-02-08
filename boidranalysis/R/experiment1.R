source("R/helpers.R")
source("R/flock_metrics.R")

library(tidyr)
library(readr)
library(purrr)
library(parallel)
library(patchwork)

# if(FALSE) {

run_experiment <- function(config, experiment_no_simulations, experiment_name = "experiment1", no_cores = 1) {
  # set up the experiment
  experiment_data_folder <- paste("Data/", experiment_name, "/")
  config$save_locations_path = experiment_data_folder

  # ensure data folder exists
  if (!dir.exists(experiment_data_folder)) {
    dir.create(experiment_data_folder)
  }

  existing_simulations <- length(list.files(experiment_data_folder, pattern="boids-data.*"))
  sample_rate <- config$sample_rate

  # collect data
  run_simulation <- function(i){
    rlang::exec(flock_detailed, !!!config)
    return(i)
  }

  runs <- experiment_no_simulations - existing_simulations
  tic("data generation")
  if(runs > 0) {
    mclapply(1:runs,
             run_simulation,
             mc.cores = no_cores)
  }
  toc()

  simulation_files <- list.files(experiment_data_folder, pattern="boids-data.*")

  # run_experiment() %>% get_convex_hull_data() %>%
  #   mutate(label =  sprintf("t_%d", time)) %>%
  #   select(label, average_cluster_volumes) %>%
  #   pivot_wider(values_from = average_cluster_volumes, names_from = label)

  # collect results

  # file_no is the ith file in a data folder
  # metric is a function that expects boid_data and returns a scalar for each time point
  tic("evaluation")
  collect_results <- function(file_no, metric, preprocess = function(x)x){
    read_csv(paste0(experiment_data_folder, simulation_files[[file_no]])) %>%
      preprocess() %>%
      metric() %>%
      into_experiment_record()
  }
  # TODO: make them all read a file once, instead of for every metric
  # collect_results2 <- function(data, metrics){
  #
  # }

  tic("results convex hull")
  results_convex_hull <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, get_convex_hull_data, preprocess = get_directional_boid_data)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("results no flocks")
  results_no_flocks <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, get_no_flocks, preprocess = get_directional_boid_data)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("results average norm vel")
  results_average_norm_vel <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, metric = get_average_norm_vel, preprocess = get_directional_boid_data)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("visualisation")
  results_stats <- tibble(
    t = 1:ncol(results_no_flocks) * sample_rate,
    mean_no_flocks = apply(results_no_flocks, 2, mean),
    var_no_flocks = apply(results_no_flocks, 2, var),
    mean_convex_hull = apply(results_convex_hull, 2, mean),
    var_convex_hull = apply(results_convex_hull, 2, var),
    mean_average_norm_vel = apply(results_average_norm_vel, 2, mean),
    # var_average_norm_vel = apply(results_average_norm_vel, 2, var),
  )

  experiment_plot_folder <- paste0(experiment_data_folder, "plots/")

  if(!dir.exists(experiment_plot_folder)) {
    dir.create(experiment_plot_folder)
  }

  figure_postfix <- paste0("_", experiment_name, ".jpg")
  figure_height <- 6
  figure_width <- 10
  figure_units <- "in"
  figure_dpi <- "retina" # screen/print/retina or a number

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_no_flocks)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time t", y = "mean_t mean # flocks",
         title = paste(
           "mean number of flocks;",
            "agents:", config$init_boids, ";",
            "iterations:", config$no_iter, ";",
            "sample:", config$sample_rate, ";"
           )
     )
  ggsave(paste0("mean_no_flocks", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(
    data = results_stats,
    aes(x = t, y = var_no_flocks)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time t", y = "var_t mean # flocks",
         title = paste("var number of flocks;",
           "agents:", config$init_boids, ";",
           "iterations:", config$no_iter, ";",
           "sample:", config$sample_rate, ";"
          )
    )

  ggsave(paste0("var_no_flocks", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_convex_hull)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time t", y = "mean_t mean flock area",
         title = paste(
           "mean flock area;",
           "agents:", config$init_boids, ";",
           "iterations:", config$no_iter, ";",
           "sample:", config$sample_rate, ";"
         )
    )
  ggsave(paste0("mean_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(
    data = results_stats,
    aes(x = t, y = var_convex_hull)
  ) +
    geom_point() +
    geom_smooth(method = "gam", level = .9) +
    theme_bw() +
    labs(x = "time t", y = "var_t mean flock area",
         title = paste(
           "var flock area;",
           "agents:", config$init_boids, ";",
           "iterations:", config$no_iter, ";",
           "sample:", config$sample_rate, ";"
         )
    )
  ggsave(paste0("var_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_average_norm_vel)
  ) +
    geom_line() +
    theme_bw() +
    labs(x = "time t", y = "mean average norm vel",
         title = paste(
           "average norm vel;",
           "agents:", config$init_boids, ";",
           "iterations:", config$no_iter, ";",
           "sample:", config$sample_rate, ";"
         )
    )
  ggsave(paste0("mean_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)


  first_n_average_norm_vel <- results_average_norm_vel %>%
    slice_head(n = 4) %>%
    mutate(sim_no = row_number()) %>%
    pivot_longer(cols = starts_with("t_"), values_to = "vals", names_pattern = "t_(.+)", names_transform = as.integer) %>%
    select(sim_no, time = name, vals) %>%
    mutate(label = sprintf("s_%d", sim_no), time = time * sample_rate) %>%
    select(time, vals, label) %>%
    pivot_wider(values_from = vals, names_from = label)

  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

  # f_1 <- ggplot(first_n_average_norm_vel, aes(x = time, y = s_1)) +
  #   geom_line() + labs(x = NULL, y = NULL)
  # f_2 <- ggplot(first_n_average_norm_vel, aes(x = time, y = s_2)) +
  #   geom_line() + labs(x = NULL, y = NULL)
  # f_3 <- ggplot(first_n_average_norm_vel, aes(x = time, y = s_3)) +
  #   geom_line() + labs(x = NULL, y = NULL)
  # f_4 <- ggplot(first_n_average_norm_vel, aes(x = time, y = s_4)) +
  #   geom_line() + labs(x = NULL, y = NULL)
  #
  #   # labs(x = "time t", y = "mean average norm vel")
  #
  # (f_1 / f_2) & (f_3 / f_4)  + plot_annotation(tag_levels = "1") # Uppercase roman numerics


  results_average_norm_vel %>%
    slice_head(n = 5) %>%
    mutate(sim_no = row_number()) %>%
    pivot_longer(cols = starts_with("t_"), values_to = "vals", names_pattern = "t_(.+)", names_transform = as.integer) %>%
    select(sim_no, time = name, vals) %>%
    mutate(label = sprintf("s_%d", sim_no), time = time * sample_rate) %>%
    ggplot(aes(x = time, y = vals)) +
    geom_line() +
    facet_grid(rows = vars(label)) +
    theme_bw() +
    labs(x = "time t", y = "mean average norm vel", title = paste(
      "average norm vel of first n simulations;",
      "agents:", config$init_boids, ";",
      "iterations:", config$no_iter, ";",
      "sample:", config$sample_rate, ";"
    ))
  ggsave(paste0("first_n_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  toc()

  toc()
}

# testing different configs -
# a -"normal"
# b - "stringy balls"
# c - "stringy"
# d - "normal" but >60K iterations
tic("experiment1_a start")

config <- get_config(
  "basic.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Thoroidal\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "experiment1_a")
toc()

tic("experiment1_a2 start")

config <- get_config(
  "basic.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^16,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Thoroidal\"}"
  )
)

no_cores <- 8
experiment_no_simulations <- 120

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "experiment1_a2")
toc()

tic("experiment1_b start")

config <- get_config(
  "string_ball_hybrid.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Thoroidal\"}"
  )
)

# no_cores <- 6
experiment_no_simulations <- 120

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "experiment1_b")
toc()

tic("experiment1_c start")

config <- get_config(
  "string.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Thoroidal\"}"
  )
)

# no_cores <- 6
experiment_no_simulations <- 120

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "experiment1_c")
toc()

# data generation: 891.305 sec elapsed
# `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'
# `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'
# `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'
# `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'
# evaluation: 29.292 sec elapsed
# Warning messages:
#   1: Removed 1 rows containing non-finite values (`stat_smooth()`).
# 2: Removed 1 rows containing missing values (`geom_point()`).
# 3: Removed 1 rows containing non-finite values (`stat_smooth()`).
# 4: Removed 1 rows containing missing values (`geom_point()`).
# 5: Removed 1 rows containing non-finite values (`stat_smooth()`).
# 6: Removed 1 rows containing missing values (`geom_point()`).
# 7: Removed 1 rows containing non-finite values (`stat_smooth()`).
# 8: Removed 1 rows containing missing values (`geom_point()`).
# > toc()
# experiment1_d start: 920.619 sec elapsed
tic("experiment1_d start")

config <- get_config(
  "basic.toml",
  overwrite = list(
    init_boids = 2^10,
    no_iter = 2^16,
    init_width = 8000,
    init_height = 8000,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Thoroidal\"}"
  )
)

# no_cores <- 6
experiment_no_simulations <- 30

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "experiment1_d")
toc()

# toJSON(config)
# }

source("R/helpers.R")
source("R/flock_metrics.R")
source("./R/directions_angles.R")

library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(parallel)
library(patchwork)
library(tictoc)
library(ggplot2)
library(RcppTOML)
# if(FALSE) {

run_experiment <- function(config, experiment_no_simulations, experiment_name = "experiment1", no_cores = 1) {
  # set up the experiment
  experiment_data_folder <- paste0("Data/", experiment_name, "/")
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

  # read_csv(paste0(experiment_data_folder, simulation_files[[file_no]])) %>%
  #   preprocess()

  # file_no is the ith file in a data folder
  # metric is a function that expects boid_data and returns a scalar for each time point
  tic("evaluation")
  collect_results <- function(file_no, metric, preprocess = function(x)x){
    read_csv(paste0(experiment_data_folder, simulation_files[[file_no]])) %>%
      preprocess() %>%
      metric() %>%
      into_experiment_record()
  }
  preprocess = function(x){ get_directional_boid_data2(x, if_else(config$boundary_config == "{\"type\": \"Toroidal\"}", F, T))}

  tic("results average norm vel")
  results_average_norm_vel <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, metric = get_average_norm_vel, preprocess = preprocess)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()
  # tic("results convex hull")
  # results_convex_hull <- mclapply(
  #   1:length(simulation_files),
  #   function(file_no){
  #     collect_results(file_no, get_convex_hull_data, preprocess =  preprocess)
  #   },
  #   mc.cores = no_cores
  # ) %>%
  #   map_dfr(bind_rows)
  # toc()
  #
  # tic("results no flocks")
  # results_no_flocks <- mclapply(
  #   1:length(simulation_files),
  #   function(file_no){
  #     collect_results(file_no, get_no_flocks, preprocess = preprocess)
  #   },
  #   mc.cores = no_cores
  # ) %>%
  #   map_dfr(bind_rows)
  # toc()
  #
  # tic("results average norm vel")
  # results_average_norm_vel <- mclapply(
  #   1:length(simulation_files),
  #   function(file_no){
  #     collect_results(file_no, metric = get_average_norm_vel, preprocess = preprocess)
  #   },
  #   mc.cores = no_cores
  # ) %>%
  #   map_dfr(bind_rows)
  # toc()
  # results voronoi counts: 730.34 sec elapsed
  # tic("results voronoi counts")
  # results_voronoi_counts <- mclapply(
  #   1:length(simulation_files),
  #   function(file_no){
  #     collect_results(file_no, metric = function(x) get_voronoi_area(x, config), preprocess = preprocess)
  #   },
  #   mc.cores = no_cores
  # ) %>%
  #   map_dfr(bind_rows)
  # toc()

  tic("visualisation")
  results_stats <- tibble(
    t = 1:ncol(results_average_norm_vel) * sample_rate,
    # mean_no_flocks = apply(results_no_flocks, 2, mean),
    # var_no_flocks = apply(results_no_flocks, 2, var),
    # mean_convex_hull = apply(results_convex_hull, 2, mean),
    # var_convex_hull = apply(results_convex_hull, 2, var),
    # mean_average_norm_vel = apply(results_average_norm_vel, 2, mean),
    var_average_norm_vel = apply(results_average_norm_vel, 2, var),
    # var_voronoi_counts = apply(results_voronoi_counts, 2, var),
    # mean_voronoi_counts = apply(results_voronoi_counts, 2, mean),
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
  figure_subtitle_config <- paste(
    "agents:", config$init_boids, ";",
    "iter:", config$no_iter, ";",
    "smpl:", config$sample_rate, ";",
    "env:", paste0(config$init_width, "×", config$init_width), ";"
  )

  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = mean_no_flocks)
  # ) +
  #   geom_point() +
  #   geom_smooth(method = "gam", level = .9) +
  #   theme_bw() +
  #   labs(x = "time t", y = "mean_t mean # flocks",
  #        title = paste(
  #          "mean number of flocks;",
  #          "agents:", config$init_boids, ";",
  #          "iterations:", config$no_iter, ";",
  #          "sample:", config$sample_rate, ";"
  #        ),
  #        subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("mean_no_flocks", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = var_no_flocks)
  # ) +
  #   geom_point() +
  #   geom_smooth(method = "gam", level = .9) +
  #   theme_bw() +
  #   labs(x = "time t", y = "var_t mean # flocks",
  #        title = "var number of flocks",
  #        subtitle = figure_subtitle_config
  #   )
  #
  # ggsave(paste0("var_no_flocks", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = mean_convex_hull)
  # ) +
  #   geom_point() +
  #   geom_smooth(method = "gam", level = .9) +
  #   theme_bw() +
  #   labs(x = "time t", y = "mean_t mean flock area",
  #        title = "mean flock area",
  #        subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("mean_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = var_convex_hull)
  # ) +
  #   geom_point() +
  #   geom_smooth(method = "gam", level = .9) +
  #   theme_bw() +
  #   labs(x = "time t", y = "var_t mean flock area",
  #        title = "var flock area",
  #        subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("var_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = mean_average_norm_vel)
  # ) +
  #   geom_line() +
  #   theme_bw() +
  #   labs(x = "time t", y = "mean average norm vel",
  #        title = "average norm vel",
  #        subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("mean_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  #
  ggplot(
    data = results_stats,
    aes(x = t, y = var_average_norm_vel)
  ) +
    geom_line() +
    theme_bw() +
    labs(x = "time t", y = "var average norm vel",
         title = "var average norm vel",
         subtitle = figure_subtitle_config
    )
  ggsave(paste0("var_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = mean_voronoi_counts)
  # ) +
  #   geom_line() +
  #   theme_bw() +
  #   labs(x = "time t", y = "mean Count(voronoi cells area < pi*R^2)",
  #        title = "average voronoi cell area",
  #        subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("mean_count_voronoi_cell_area", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # ggplot(
  #   data = results_stats,
  #   aes(x = t, y = var_voronoi_counts)
  # ) +
  #   geom_line() +
  #   theme_bw() +
  #   labs(x = "time t", y = "var Count(voronoi cells area < pi*R^2)",
  #        title = "var voronoi cell area",
  #        subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("var_count_voronoi_cell_area", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  #
  # first_n_average_norm_vel <- results_average_norm_vel %>%
  #   slice_head(n = 4) %>%
  #   mutate(sim_no = row_number()) %>%
  #   pivot_longer(cols = starts_with("t_"), values_to = "vals", names_pattern = "t_(.+)", names_transform = as.integer) %>%
  #   select(sim_no, time = name, vals) %>%
  #   mutate(label = sprintf("s_%d", sim_no), time = time * sample_rate) %>%
  #   select(time, vals, label) %>%
  #   pivot_wider(values_from = vals, names_from = label)

  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

  # results_average_norm_vel %>%
  #   slice_head(n = 5) %>%
  #   mutate(sim_no = row_number()) %>%
  #   pivot_longer(cols = starts_with("t_"), values_to = "vals", names_pattern = "t_(.+)", names_transform = as.integer) %>%
  #   select(sim_no, time = name, vals) %>%
  #   mutate(label = sprintf("s_%d", sim_no), time = time * sample_rate) %>%
  #   ggplot(aes(x = time, y = vals)) +
  #   geom_line() +
  #   facet_grid(rows = vars(label)) +
  #   theme_bw() +
  #   labs(x = "time t", y = "mean average norm vel",
  #     title = "average norm vel of first n simulations",
  #     subtitle = figure_subtitle_config
  #   )
  # ggsave(paste0("first_n_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  # toc()

  # 4 by 4 faceted plots
  tic("4 by 4 faceted plots")
  no_facets <- 5
  batch_factor <-10
  direction_df <- mclapply(
    1:no_facets,
    function(file_no){
      # somehow there is a few (units) records missing here
      df <- read_csv(paste0(experiment_data_folder, simulation_files[[file_no]])) %>%
        preprocess() %>%
        mutate(result_no = file_no)

      n <- nrow(df)
      section_size <- n %/% batch_factor
      breaks = c(seq(1, n, section_size + if_else(n %% batch_factor != 0, 1, 0)), n + 1)
      labels = 1:batch_factor
      return(
        df %>%
          mutate(timeline_facet = cut(1:n+1, breaks = breaks, labels = labels))
      )
    },
    mc.cores = no_cores
  ) %>%
    purrr::map_dfr(dplyr::bind_rows)


  no_bins <- 180
  count_em_up <- function(vec, range) {
    # define all possible values
    all_values <- min(range):max(range)
    # get value counts
    value_counts <- table(vec)

    # add missing values with a count of 0
    value_counts <- value_counts[match(all_values, names(value_counts))]
    value_counts[is.na(value_counts)] <- 0
    names(value_counts) <- all_values

    return(value_counts)
  }

  direction_counts_df <- direction_df %>%
    group_by(timeline_facet, result_no) %>%
    mutate(heading_bin = cut(headings, breaks = seq(0, 2*pi, length.out = no_bins + 1), labels = F)) %>%
    mutate(bearing_bin = cut(bearings, breaks = seq(-pi, pi, length.out = no_bins + 1), labels = F)) %>%
    reframe(bin = 1:no_bins, heading_count = c(count_em_up(heading_bin, c(1, no_bins))), bearing_count = c(count_em_up(bearing_bin, c(1, no_bins))))

  timeline_labs <- sapply(1:batch_factor, function(x)paste0("k = ",x))
  names(timeline_labs) <- 1:batch_factor

  result_labs <- sapply(1:no_facets, function(x)paste0("run ",x,"."))
  names(result_labs) <- 1:no_facets

  ggplot(direction_counts_df, aes(x = bin, y = heading_count)) +
    coord_polar(theta = "x", start = pi) +
    geom_bar(stat = "identity", fill = "orange", width = .9) +
    labs(x =  paste0("k-th iterration section"), subtitle = figure_subtitle_config) +
    theme_bw() +
    facet_grid(
      result_no~timeline_facet,
      labeller = labeller(
        timeline_facet = timeline_labs,
        result_no = result_labs)
    )
  ggsave(paste0("heading_top_k_overview", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(direction_counts_df, aes(x = bin, y = log(heading_count))) +
    coord_polar(theta = "x", start = pi) +
    geom_bar(stat = "identity", fill = "orange", width = .9) +
    labs(x =  paste0("k-th iterration section"), subtitle = figure_subtitle_config) +
    theme_bw() +
    facet_grid(
      result_no~timeline_facet,
      labeller = labeller(
        timeline_facet = timeline_labs,
        result_no = result_labs)
    )
  ggsave(paste0("heading_log_top_k_overview", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(direction_counts_df, aes(x = bin, y = bearing_count)) +
    coord_polar(theta = "x", start = pi) +
    geom_bar(stat = "identity", fill = "deeppink", width = .9) +
    labs(x = paste0("k-th iterration section "), subtitle = figure_subtitle_config) +
    theme_bw() +
    facet_grid(
      result_no~timeline_facet,
      labeller = labeller(
        timeline_facet = timeline_labs,
        result_no = result_labs)
    )
  ggsave(paste0("bearing_top_k_overview", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)


  ggplot(direction_counts_df, aes(x = bin, y = log(bearing_count))) +
    coord_polar(theta = "x", start = pi) +
    geom_bar(stat = "identity", fill = "deeppink", width = .9) +
    labs(x = paste0("k-th iterration section "), subtitle = figure_subtitle_config) +
    theme_bw() +
    facet_grid(
      result_no~timeline_facet,
      labeller = labeller(
        timeline_facet = timeline_labs,
        result_no = result_labs)
    )
  ggsave(paste0("bearing_log_top_k_overview", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  tau <- 5
  labels = direction_df$result_no

  tic("results voronoi counts")
  curvature_result <- mclapply(
    (direction_df %>% select(result_no) %>% distinct())$result_no,
    function(rn){
      direction_df %>% filter(result_no == rn) %>%
        select(-timeline_facet) %>%
        get_curvature_order_data(config, tau) %>%
        mutate(result_no = rep(rn, n()), label_time = 1:n() * config$sample_rate, label = sprintf("s_%d", result_no))
      # bind_cols(tibble(result_no = rep(rn, length(time)), label_time = 1:length(time) * config$sample_rate))
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  ggplot(curvature_result, aes(x = label_time, y = rop)) +
    theme_bw() +
    labs(
      title = paste0("mean flock path curvature with tau = ", tau),
      subtitle = figure_subtitle_config,
      x = "rate"
    ) +
    geom_line() +
    facet_grid(rows = vars(label))
  ggsave(paste0("curvature_top_k_overview", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  toc()

  tic("direction data graphs")
  # testing these
  direction_df %>%
    arrange(time) %>%
    group_by(result_no, id) %>%
    mutate(bearings_ma =rollapply(bearings, tau ,mean,align='right',fill=NA)) %>%
    ungroup() %>%
    mutate(noise = cluster_id == 0) %>%
    group_by(result_no, time, noise) %>%
    summarise(bearings_ma = mean(bearings_ma)) %>%
    ungroup() %>%
    filter(!is.na(bearings_ma)) %>%
    group_by(result_no) %>%
    mutate(time_label = 1:n() * config$sample_rate, label = sprintf("s_%d", result_no)) %>%
    ggplot(aes(x = time_label, y = bearings_ma)) +
    geom_point(size = 1.2, alpha = .3, aes(color = noise))  +
    stat_summary(aes(y = bearings_ma, color = noise), alpha = .8, fun=mean, geom="line") +
    theme_bw() +
    facet_grid(rows = vars(label)) +
    labs(
      title = paste0("mean average flock curvature tau = ", tau),
      subtitle = figure_subtitle_config,
      x = "time"
    )
  #   mutate(label = sprintf("s_%d", sim_no), time = time * sample_rate) %>%
  #   ggplot(aes(x = time, y = vals)) +
  #   geom_line() +
  #   facet_grid(rows = vars(label))
  ggsave(paste0("dir_data_graph_1", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  direction_df %>%
    arrange(time) %>%
    group_by(result_no, id) %>%
    mutate(breaings_ma=rollapply(bearings, tau,mean,align='right',fill=NA)) %>%
    ungroup() %>%
    group_by(time, cluster_id) %>%
    summarise(breaings_ma = mean(breaings_ma)) %>%
    mutate(noise = cluster_id == 0) %>%
    ungroup() %>%
    filter(!is.na(breaings_ma)) %>%
    mutate(time_label = 1:n() * config$sample_rate) %>%
    ggplot(aes(x = time_label, y = breaings_ma)) +
    geom_point(size = .5, alpha = .1, aes(color = noise))  +
    stat_summary(aes(y = breaings_ma, color = noise), fun=mean, geom="line") +
    theme_bw() +
    labs(
      title = paste0("mean average flock curvature tau = ", tau, "; coerced runs = ", length(table(direction_df$result_no))),
      subtitle = figure_subtitle_config,
      x = "time"
    )
  ggsave(paste0("dir_data_graph_2", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  toc()

  toc()

  toc()
  # preserve config at last
  config$experiment_name = experiment_name
  jsonlite::write_json(config, paste0(experiment_data_folder, "config.json"))
}

# testing different configs -
# a -"normal"
# b - "stringy balls"
# c - "stringy"
# d - "normal" but >60K iterations

# s - variant is ran with BoidTracker instead of SpatHashTracker
# t - variant has been rerun with the toroidial distance instead of the original euclidea

tic("experiment1_a start")

config <- get_config(
  "basic.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
    # distance_config = "{\"type\": \"EucToroidal\"}"
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
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
  )
)

# no_cores <- 6
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
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
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
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
  )
)

# no_cores <- 8
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

# tic("experiment1_d start")
#
# config <- get_config(
#   "basic.toml",
#   overwrite = list(
#     init_boids = 2^10,
#     no_iter = 2^16,
#     init_width = 8000,
#     init_height = 8000,
#     sample_rate = 32,
#     boundary_config = "{\"type\": \"Toroidal\"}",
#     distance_config = "{\"type\": \"EucEnclosed\"}"
#   )
# )
#
# # no_cores <- 6
# experiment_no_simulations <- 30
#
# run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "experiment1_d")
# toc()

# toJSON(config)
# }

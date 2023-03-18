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

  existing_simulations <- length(list.files(experiment_data_folder, pattern="(prepro_)?boids-data.*"))
  sample_rate <- config$sample_rate

  # collect data
  run_simulation <- function(i){
    rlang::exec(flock_detailed, !!!config)
    return(i)
  }

  # preserve config
  config$experiment_name = experiment_name
  jsonlite::write_json(config, paste0(experiment_data_folder, "config.json"))

  runs <- experiment_no_simulations - existing_simulations
  tic("data generation")
  if(runs > 0) {
    mclapply(1:runs,
             run_simulation,
             mc.cores = no_cores)
  }
  toc()

  simulation_files <- list.files(experiment_data_folder, pattern="(prepro_)?boids-data.*")[1:experiment_no_simulations]

  fetch_file <- function(file, init_width, init_height, no_boids, version = 3){
    if (version != 3) {
      return(read_csv(file))
    }

   # this inflates the files by about 1.5x, but is 7.5x times quicker and most importantly is able to chew through >300MB files
    prep_file <-  file
    if(!startsWith(basename(file), "prepro_")) {
      prep_file <- file.path(dirname(file), paste0("prepro_", basename(file)))
    }

    if (!file.exists(prep_file)) {
      # prepro file does not exist
      if (!file.exists(file)) {
        # nor does the original
        stop(paste0("Did not find the raw or preprocessed data file, file: ", file))
      } else {
        # preprocess the file
        preprocess_file(file, init_width = init_width, init_height = init_height, no_boids = no_boids)
        # remove the original
        if (file.exists(prep_file)) {
          file.remove(file)
        } else {
          stop(paste0("Something went wrong with file preprocessing, file: ", file))
        }
      }
    }

    return(read_csv(prep_file))
  }

  # file_no is the ith file in a data folder
  # metric is a function that expects boid_data and returns a scalar for each time point
  tic("evaluation")
  collect_results <- function(file_no, metric, process = function(x)x, fetch = fetch_file){
    fetch_file(
      paste0(experiment_data_folder, simulation_files[[file_no]]),
      init_width = config$init_width,
      init_height = config$init_height,
      no_boids = config$init_boids
      ) |>
      process() |>
      metric() |>
      into_experiment_record()
  }

  # process <-  function(x){ get_directional_boid_data2(x, if_else(config$boundary_config == "{\"type\": \"Toroidal\"}", F, T))}
  process = function(x, version = 3){
    remove_boundary = if_else(config$boundary_config == "{\"type\": \"Toroidal\"}", F, T);

    if (version != 3) {
      get_directional_boid_data2(x, remove_boundary)
    } else {
      get_directional_boid_data3(x, remove_boundary)
    }
  }

  tic("results average norm vel")
  results_average_norm_vel <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, metric = get_average_norm_vel, process = process)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("results convex hull")
  results_convex_hull <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, get_convex_hull_data, process =  process)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("results no flocks")
  results_no_flocks <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, get_no_flocks, process = process)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("results average norm vel")
  results_average_norm_vel <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, metric = get_average_norm_vel, process = process)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()
  # results voronoi counts: 730.34 sec elapsed
  tic("results voronoi counts")
  results_voronoi_counts <- mclapply(
    1:length(simulation_files),
    function(file_no){
      collect_results(file_no, metric = function(x) get_voronoi_area(x, config), process = process)
    },
    mc.cores = no_cores
  ) %>%
    map_dfr(bind_rows)
  toc()

  tic("visualisation")
  results_stats <- tibble(
    t = 1:ncol(results_average_norm_vel) * sample_rate,
    mean_no_flocks = apply(results_no_flocks, 2, mean),
    var_no_flocks = apply(results_no_flocks, 2, var),
    mean_convex_hull = apply(results_convex_hull, 2, mean),
    var_convex_hull = apply(results_convex_hull, 2, var),
    mean_average_norm_vel = apply(results_average_norm_vel, 2, mean),
    var_average_norm_vel = apply(results_average_norm_vel, 2, var),
    var_voronoi_counts = apply(results_voronoi_counts, 2, var),
    mean_voronoi_counts = apply(results_voronoi_counts, 2, mean),
  )

  write_csv(results_stats, paste0(experiment_data_folder, "results_stats.csv"))

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
         ),
         subtitle = figure_subtitle_config
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
         title = "var number of flocks",
         subtitle = figure_subtitle_config
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
         title = "mean flock area",
         subtitle = figure_subtitle_config
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
         title = "var flock area",
         subtitle = figure_subtitle_config
    )
  ggsave(paste0("var_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_average_norm_vel)
  ) +
    geom_line() +
    theme_bw() +
    labs(x = "time t", y = "mean average norm vel",
         title = "average norm vel",
         subtitle = figure_subtitle_config
    )
  ggsave(paste0("mean_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)


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

  ggplot(
    data = results_stats,
    aes(x = t, y = mean_voronoi_counts)
  ) +
    geom_line() +
    theme_bw() +
    labs(x = "time t", y = "mean Count(voronoi cells area < pi*R^2)",
         title = "average voronoi cell area",
         subtitle = figure_subtitle_config
    )
  ggsave(paste0("mean_count_voronoi_cell_area", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  ggplot(
    data = results_stats,
    aes(x = t, y = var_voronoi_counts)
  ) +
    geom_line() +
    theme_bw() +
    labs(x = "time t", y = "var Count(voronoi cells area < pi*R^2)",
         title = "var voronoi cell area",
         subtitle = figure_subtitle_config
    )
  ggsave(paste0("var_count_voronoi_cell_area", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

  first_n_average_norm_vel <- results_average_norm_vel %>%
    slice_head(n = 4) %>%
    mutate(sim_no = row_number()) %>%
    pivot_longer(cols = starts_with("t_"), values_to = "vals", names_pattern = "t_(.+)", names_transform = as.integer) %>%
    select(sim_no, time = name, vals) %>%
    mutate(label = sprintf("s_%d", sim_no), time = time * sample_rate) %>%
    select(time, vals, label) %>%
    pivot_wider(values_from = vals, names_from = label)

  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

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
    labs(x = "time t", y = "mean average norm vel",
      title = "average norm vel of first n simulations",
      subtitle = figure_subtitle_config
    )
  ggsave(paste0("first_n_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)
  toc()

  # 4 by 4 faceted plots
  tic("4 by 4 faceted plots")
  no_facets <- 5
  batch_factor <-10
  direction_df <- mclapply(
    1:no_facets,
    function(file_no){
      # somehow there is a few (units) records missing here
     df <- fetch_file(paste0(experiment_data_folder, simulation_files[[file_no]])) %>%
       process() %>%
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


  # no_bins <- 180
  # count_em_up <- function(vec, range) {
  #   # define all possible values
  #   all_values <- min(range):max(range)
  #   # get value counts
  #   value_counts <- table(vec)
  #
  #   # add missing values with a count of 0
  #   value_counts <- value_counts[match(all_values, names(value_counts))]
  #   value_counts[is.na(value_counts)] <- 0
  #   names(value_counts) <- all_values
  #
  #   return(value_counts)
  # }

  direction_counts_df <- direction_df %>%
  group_by(timeline_facet, result_no) %>%
  get_directional_counts(no_bins = 180)
  # mutate(heading_bin = cut(headings, breaks = seq(0, 2*pi, length.out = no_bins + 1), labels = F)) %>%
  # mutate(bearing_bin = cut(bearings, breaks = seq(-pi, pi, length.out = no_bins + 1), labels = F)) %>%
  # reframe(bin = 1:no_bins, heading_count = c(count_em_up(heading_bin, c(1, no_bins))), bearing_count = c(count_em_up(bearing_bin, c(1, no_bins))))

  write_csv(direction_counts_df, paste0(experiment_data_folder, "direction_counts"))


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

  # quit early for large no_iter numbers as the code bellow is really memory inefficient

  if (config$no_iter > 2^16) {
    toc()
    toc()
    toc()
    return()
  }

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
}

# testing different configs -
# a - toroidal b/d vs reflective/euclidean b/d
# b - "stringy ball hybrid, mars = modified diamonds"
# c - "large stringy"
# e - ?
# f - one run of g with a longer sequence
# g - testing out different noise levels
# h - do longer runs trully stabilise?
tic("experiment2_a1 start")

config <- get_config(
  "basic.toml",
  overwrite = list(
    init_boids = 2^11,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 64,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_a1")
)
toc()

tic("experiment2_a12 start")

config <- get_config(
  "basic2.toml",
  overwrite = list(
    init_boids = 2^11,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 64,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_a12")
)
toc()

tic("experiment2_a2 start")

config <- get_config(
  "basic2.toml",
  overwrite = list(
    init_boids = 2^11,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 64,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_a2")
toc()

tic("experiment2_a3 start")

config <- get_config(
  "basic2.toml",
  overwrite = list(
    init_boids = 2^11,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 64,
    boundary_config = "{\"type\": \"Reflective\"}",
    distance_config = "{\"type\": \"EucEnclosed\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120

run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_a3")
toc()




tic("experiment2_b1 start")

config <- get_config(
  "string_ball_hybrid.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 2000,
    init_height = 2000,
    sample_rate = 64,
    field_of_vision = 200,
    min_speed = 2,
    max_speed = 2,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_b1")
)
toc()


tic("experiment2_b2 start")

config <- get_config(
  "string_ball_hybrid.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 2000,
    init_height = 2000,
    sample_rate = 64,
    min_speed = 2,
    max_speed = 2,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 120
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_b2")
)
toc()

# # this took 13.214 hours on 4 M1 Pro P cores
tic("experiment1_c1 start")

config <- get_config(
  "mars.toml",
  overwrite = list(
    init_boids = 2^13,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 64,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)

no_cores <- 4
experiment_no_simulations <- 120

tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_c1")
)
toc()

tic("experiment1_c2 start")

config <- get_config(
  "large_string.toml",
  overwrite = list(
    init_boids = 2^13,
    no_iter = 2^15,
    init_width = 8000,
    init_height = 8000,
    sample_rate = 64,
    boundary_config = "{\"type\": \"Toroidal\"}"
  )
)

no_cores <- 3
experiment_no_simulations <- 120

tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_c2")
)
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
    boundary_config = "{\"type\": \"Toroidal\"}"
  )
)

no_cores <- 6
experiment_no_simulations <- 30

tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0301_experiment2_d"),
)
toc()

tic("0307_experiment2_e1_1 start")
config <- get_config(
  "density.toml",
  overwrite = list(
    init_boids = 2^8,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)
no_cores <- 8
experiment_no_simulations <- 50
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0307_experiment2_e1_1")
)
toc()

tic("0307_experiment2_e1_2 start")
config <- get_config(
  "density.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)
no_cores <- 8
experiment_no_simulations <- 50
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0307_experiment2_e1_2")
)
toc()

tic("0307_experiment2_e1_3 start")
config <- get_config(
  "density.toml",
  overwrite = list(
    init_boids = 2^10,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)
no_cores <- 8
experiment_no_simulations <- 50
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0307_experiment2_e1_3")
)
toc()

tic("0307_experiment2_e1_4 start")
config <- get_config(
  "density.toml",
  overwrite = list(
    init_boids = 2^11,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)
no_cores <- 8
experiment_no_simulations <- 50
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0307_experiment2_e1_4")
)
toc()

tic("0307_experiment2_e1_5 start")
config <- get_config(
  "density.toml",
  overwrite = list(
    init_boids = 2^11 + 2^10,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)
no_cores <- 4
experiment_no_simulations <- 50
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0307_experiment2_e1_5")
)
toc()

tic("0307_experiment2_e1_6 start")
config <- get_config(
  "density.toml",
  overwrite = list(
    init_boids = 2^12,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
  )
)
no_cores <- 4
experiment_no_simulations <- 50
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0307_experiment2_e1_6")
)
toc()

tic("0314_experiment2_f1 start")
config <- get_config(
  "smaller_string.toml",
  overwrite = list(
    init_boids = 2^10,
    no_iter = 2^15,
    init_width = 1600,
    init_height = 1600,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}",
    sep_bias = TRUE
  )
)
no_cores <- 6
experiment_no_simulations <- 48
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0314_experiment2_f1")
)
toc()

for (eta in 0:30) {
  name <- paste0("0314_experiment2_g_", sprintf("%02d", eta))
  tic(name)

  print(name)
  print(eta / 100)
  config <- get_config(
    "string_m.toml",
    overwrite = list(
      init_boids = 2^10,
      no_iter = 2^15,
      init_width = 1600,
      init_height = 1600,
      sample_rate = 32,
      boundary_config = "{\"type\": \"Toroidal\"}",
      distance_config = "{\"type\": \"EucToroidal\"}",
      sep_bias = TRUE,
      wander_coef = eta/100
    )
  )
  no_cores <- 8
  experiment_no_simulations <- 8
  tryCatch(
    expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = name)
  )
  toc()
}

eta_results_stats = tibble()
for (eta in 0:30) {
  name <- paste0("Data/0314_experiment2_g_", sprintf("%02d", eta))

  eta_results_stats <- eta_results_stats %>% bind_rows(
    read_csv(paste0(name, "/results_stats.csv")) %>%
      mutate(eta = eta / 100)
  )
}

# h and i are the same, except for decreased no_iter
tic("0314_experiment2_i1 start")
config <- get_config(
  "string_s.toml",
  overwrite = list(
    init_boids = 1417,
    no_iter = 2^17,
    init_width = 500,
    init_height = 500,
    sample_rate = 32,
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}",
    sep_bias = TRUE
  )
)

no_cores <- 2
experiment_no_simulations <- 8
tryCatch(
  expr = run_experiment(config, experiment_no_simulations, no_cores = no_cores, experiment_name = "0314_experiment2_i1")
)



deg2rad <- function(deg) deg*pi/180
rad2deg <- function(rad) rad*180/pi

toroidal_vec_pc <- function(pc1, pc2, size, max) {
  if (is.na(pc1) || is.na(pc2)) {
    return(NA)
  }
  d_pc <- pc2 - pc1
  if (abs(d_pc) > max) {
    return(d_pc + if_else(d_pc < 0, size, -size))
  } else {
    return(d_pc)
  }
}

# assumes data is written groups of rows, each representing a single unit, naturally ordered by time
# takes the group and interlaces it with the rest, keeping the time order and thus creating a relay of the simulation as
# the data points have been collected in the first place
# OBSOLETE
order_by_time <- function(data) {
  distinct_no <- nrow(data %>% dplyr::select(id) %>% unique())
  ordered_boid_data <- data %>%
    group_by(id) %>%
    mutate(
      order = id + (row_number(id) - 1) * distinct_no,
      timestep = row_number(id)
    ) %>%potif
    ungroup() %>%
    arrange(order)
}

# given a series of data points ordered by time, assuming the number of units in
# a flock is constant, returns a single n-th "page", which represents the
# location of all entities at the n-th time point
# OBSOLETE
nth_time_point <- function(data, time_n) {
  distinct_no <- nrow(data %>% dplyr::select(id) %>% unique())
  l <- nrow(data)

  tail <- l - (distinct_no * (time_n - 1))
  head <- distinct_no

  data %>%
    slice_tail(n = tail) %>%
    slice_head(n = head)
}

# TODO: everything that assumes a continuous time-series will have to be reworked ;(, like the two functions above
remove_boundary_data <- function(data, sensory_distance, init_width, init_height, time_dependency) {
  # get the row numbers of boundary data
  boundary_rns <-  data%>%
    rowwise() %>%
    mutate(edge = abs(x) > init_width / 2 - sensory_distance || abs(y) > init_height / 2 - sensory_distance) %>%
    ungroup() %>%
    mutate(rn = row_number()) %>%
    filter(edge == T) %>%
    select(rn)

  # because we have two operations which do a diff on the whole data series,
  # if _t is a boundary, up to _t+2 are affected

  # get row numbers for _t+1
  affected_rns <- tibble(rn = boundary_rns$rn + 1)

  # get row numbers for _t+i
  for(i in 2:time_dependency){
    affected_rns <- bind_rows(affected_rns, tibble(rn = boundary_rns$rn + i))
  }

  # combine boundary and affected row numbers
  affected_rns <- bind_rows(affected_rns, boundary_rns) %>%
    # de-duplicate
    distinct(rn) %>%
    # ensure we have not crossed over original row numbers
    filter(rn <= nrow(data))
  # cleanup
  remove(boundary_rns)

  return(affected_rns$rn)
}

library(RcppTOML)

if (FALSE) {
  overwrite = list(sample_rate = 32, no_iter = 2^6)
  A <- get_config("string_ball_hybrid.toml")
  B <- get_config("string_ball_hybrid.toml", overwrite = overwrite)
  B["no_iter"] <- 2^6
  rlang::exec(flock_detailed,!!!A)
  rlang::exec(flock_detailed,!!!B)

  config_name <- "string_ball_hybrid.toml"
}

# given a config name and a list of properties to overwrite, will fetch the
# config settings in a list for simulation execution and do the overwrites
get_config <- function(config_name, overwrite = list()) {
  config_path <- paste0("../boids_app/configs/", config_name)
  if (!file.exists(config_path)) {
    stop("The config file does not exist")
  }

  config <- parseTOML(config_path)

  # must contain all values that are used to call flock_detailed
  converted_config <- list()

  converted_config$no_iter = 2^3
  converted_config$init_boids = config$no_boids
  converted_config$save_locations_path = "" #config$save # config$save_timestamp
  converted_config$sample_rate = config$sample_rate
  converted_config$init_width = config$init_width
  converted_config$init_height = config$init_height
  converted_config$sensory_distance = config$sensory_distance
  converted_config$allignment_coef = config$allignment_coefficient
  converted_config$allignment_trs_coef = config$allignment_treshold_coefficient
  converted_config$cohesion_coef = config$cohesion_coefficient
  converted_config$cohesion_trs_coef = config$cohesion_treshold_coefficient
  converted_config$separation_coef = config$separation_coefficient
  converted_config$separation_trs_coef = config$separation_treshold_coefficient
  converted_config$min_speed = config$min_speed
  converted_config$max_speed = config$max_speed
  converted_config$max_steering = config$max_steering
  converted_config$field_of_vision = config$field_of_vision
  converted_config$dbscan_clustering = TRUE
  converted_config$boundary_config = "{\"type\": \"Toroidal\"}"
  converted_config$distance_config = "{\"type\": \"EucToroidal\"}"
  converted_config$rules_impl = FALSE
  converted_config$wander_on = config$wander_on
  converted_config$wander_rate = config$wander_rate
  converted_config$wander_radius = config$wander_radius
  converted_config$wander_coef = config$wander_coefficient
  converted_config$wander_distance = config$wander_distance
  converted_config$baseline_speed = config$baseline_speed

  if (length(overwrite) == 0) {
    return(converted_config)
  }
  # takes everything on the left side that is different or does not exist in the
  # right side and overwrites it with what is on the right
  m <- tibble(merge(converted_config, overwrite, all = T, sort = F))
  mask <- (m[1,] != m[2,] | is.na(m[1,]))[1,]
  mask[is.na(mask)] <- FALSE
  m[1,][mask] <- m[2,][mask]
  merged_config <- m[1,]

  return(c(merged_config[1, names(converted_config)]))
}

into_experiment_record <- function(x){
  mutate(x, label =  sprintf("t_%d", time)) %>%
    select(3, 2) %>%
    pivot_wider(values_from = 2, names_from = label)
}

# takes in a dataset (any form of DF) and returns a subsample of every n-th record
sub_sample_data <- function(data, start = 1, every_nth = 1) {
  if (start > every_nth) {
    stop("start is greate than every_nth, you will try to subset beyond the range of records")
  }
  data %>%
    slice(seq(start, nrow(data), every_nth))
}

# data <- boid_data
#
# boid_data %>%
#   group_by(id) %>%
#   mutate(lag_x = dplyr::lag(x, order_by = id)) %>%
#   rowwise() %>%
#   mutate(dx = toroidal_vec_pc(x, lag_x, config$init_width, config$init_width / 2))
#
# boid_data %>%
#   group_by(id) %>%
#   do(mutate(., lag_x = lag(x), lag_y = lag(y))) %>%
#   rowwise() %>%
#   mutate( # the toroidial trick will work on eucledian as well
#     dx = toroidal_vec_pc(lag_x, x, config$init_width, config$init_width / 2),
#     dy = toroidal_vec_pc(lag_y, y, config$init_height, config$init_height / 2)
#     ) %>%
#   ungroup() %>%
#   group_by(id) %>%
#   mutate(headings = get_headings2(dx, dy)) %>%
#   slice_tail(n = -1)

get_directional_boid_data1 <- function(data, remove_boundary = F) {
  res <- data %>%
    group_by(id) %>%
    reframe(
      dx = diff(x, 1),
      dy = diff(y, 1),
      id = id[2:n()],
      x = x[2:n()],
      y = y[2:n()],
      headings = get_headings2(dx, dy),
      cluster_id = cluster_id[2:n()],
      n_neighbours = n_neighbours[2:n()],
      time = time[2:n()]
    ) %>%
    group_by(id) %>%
    reframe(
      id = id[2:n()],
      bearings = get_bearings(headings), # returns n-1 records as it does a diff of 1
      headings = headings[2:length(headings)],
      x = x[2:n()],
      y = y[2:n()],
      dx = dx[2:n()],
      dy = dy[2:n()],
      cluster_id = cluster_id[2:n()],
      time = time[2:n()]
    )
  if (!remove_boundary) {
    return(res)
  }


  res %>%
    slice(-tail(
      remove_boundary_data(data,
                           sensory_distance = config$sensory_distance * max(config$allignment_trs_coef, config$cohesion_trs_coef, config$separation_trs_coef),
                           init_width = config$init_width,
                           init_height = config$init_height,
                           time_dependency = 2),
      -2
    ))
}

# version 2 solves continuous boundary condition
get_directional_boid_data2 <- function(data, remove_boundary = F) {
  res <- data %>%
    group_by(id) %>%
    do(mutate(., lag_x = lag(x), lag_y = lag(y))) %>%
    rowwise() %>%
    mutate( # the toroidial trick will work on eucledian as well
      dx = toroidal_vec_pc(lag_x, x, config$init_width, config$init_width / 2),
      dy = toroidal_vec_pc(lag_y, y, config$init_height, config$init_height / 2)
    ) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(headings = get_headings2(dx, dy)) %>%
    slice_tail(n = -1) %>%
    reframe(
      id = id[2:n()],
      bearings = get_bearings(headings), # returns n-1 records as it does a diff of 1
      headings = headings[2:length(headings)],
      x = x[2:n()],
      y = y[2:n()],
      dx = dx[2:n()],
      dy = dy[2:n()],
      cluster_id = cluster_id[2:n()],
      time = time[2:n()]
    )
  if (!remove_boundary) {
    return(res)
  }


  res %>%
    slice(-tail(
      remove_boundary_data(data,
                           sensory_distance = config$sensory_distance * max(config$allignment_trs_coef, config$cohesion_trs_coef, config$separation_trs_coef),
                           init_width = config$init_width,
                           init_height = config$init_height,
                           time_dependency = 2),
      -2
    ))
}

# version 3 relies on boidr to get the dx/dy values (files >300MB were a problem)
get_directional_boid_data3 <- function(data, remove_boundary = F) {
  res <- data %>%
    # group_by(id) %>%
    # do(mutate(., lag_x = lag(x), lag_y = lag(y))) %>%
    # rowwise() %>%
    # mutate( # the toroidial trick will work on eucledian as well
    #   dx = toroidal_vec_pc(lag_x, x, config$init_width, config$init_width / 2),
    #   dy = toroidal_vec_pc(lag_y, y, config$init_height, config$init_height / 2)
    # ) %>%
    # ungroup() %>%
    group_by(id) %>%
    mutate(headings = get_headings2(dx, dy)) %>%
    slice_tail(n = -1) %>%
    reframe(
      id = id[2:n()],
      bearings = get_bearings(headings), # returns n-1 records as it does a diff of 1
      headings = headings[2:length(headings)],
      x = x[2:n()],
      y = y[2:n()],
      dx = dx[2:n()],
      dy = dy[2:n()],
      cluster_id = cluster_id[2:n()],
      time = time[2:n()]
    )
  if (!remove_boundary) {
    return(res)
  }


  res %>%
    slice(-tail(
      remove_boundary_data(data,
                           sensory_distance = config$sensory_distance * max(config$allignment_trs_coef, config$cohesion_trs_coef, config$separation_trs_coef),
                           init_width = config$init_width,
                           init_height = config$init_height,
                           time_dependency = 2),
      -2
    ))
}

my_wedge <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  a[1] * b[2] - a[2] * b[1]
}
library(zoo)

get_curvature_order_data <- function(data, config, tau){
  wedges <- data %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(lead_dx = lead(dx), lead_dy = lead(dy)) %>%
    rowwise() %>%
    mutate(wedge_i = my_wedge(c(dx, dy), c(lead_dx, lead_dy))) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(wedge_tau = rollapply(wedge_i, width = tau, FUN = function(x)sum(x / config$max_speed ^ 2) / tau, fill = NA, align = "left")) %>%
    select(-lead_dx, -lead_dy) %>%
    filter(!is.na(wedge_tau)) %>%
    select(time, cluster_id, id, wedge_tau) %>%
    group_by(time) %>%
    reframe(rop = mean(wedge_tau))

  return(wedges)
}

# no_bins <- 180
#
#
# direction_counts_df <- direction_df %>%
#   group_by(timeline_facet, result_no) %>%
#   mutate(heading_bin = cut(headings, breaks = seq(0, 2*pi, length.out = no_bins + 1), labels = F)) %>%
#   mutate(bearing_bin = cut(bearings, breaks = seq(-pi, pi, length.out = no_bins + 1), labels = F)) %>%
#   reframe(bin = 1:no_bins, heading_count = c(count_em_up(heading_bin, c(1, no_bins))), bearing_count = c(count_em_up(bearing_bin, c(1, no_bins))))
#

get_directional_counts <- function (data, no_bins) {
  .count_em_up <- function(vec, range) {
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

  data %>%
    mutate(heading_bin = cut(headings, breaks = seq(0, 2*pi, length.out = no_bins + 1), labels = F)) %>%
    mutate(bearing_bin = cut(bearings, breaks = seq(-pi, pi, length.out = no_bins + 1), labels = F)) %>%
    reframe(bin = 1:no_bins, heading_count = c(.count_em_up(heading_bin, c(1, no_bins))), bearing_count = c(.count_em_up(bearing_bin, c(1, no_bins)))) %>%
    return()
}

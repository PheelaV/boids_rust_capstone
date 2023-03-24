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

# get_var_bearing <- function(data) {
#   # for each time point, count the number of flocks, we are treating 0 (set of
#   # entities not belonging to any flock
#   data %>%
#     filter(cluster_id != 0) %>%
#     group_by(time) %>%
#     distinct(cluster_id) %>%
#     summarise(no_flocks = n())
# }

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

# relies on boid data to retrieve the number of voronoi areas that are bellow a treshold
# evaluated by the runtime configuration
get_voronoi_area <- function(data, config) {
  sensory_distance = max(
    max(config$cohesion_trs_coef * config$sensory_distance,
        config$separation_trs_coef * config$sensory_distance)
    ,config$alignment_trs_coef * config$sensory_distance)

  treshold <- pi * sensory_distance^2

   data%>%
    group_by(time) %>%
    reframe(voronoi_areas = list(get_voronoi_areas(x, y, config$init_width, config$init_height))) %>%
    rowwise() %>%
    mutate(time, C_voronoi = sum(voronoi_areas < treshold) / config$init_boids) %>%
    select(time, C_voronoi)
}


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
}

if (FALSE) {
  get_voronoi_area(boid_data, config)

  data %>% select(time) %>% distinct() %>% summarise(no_iter = n())
}

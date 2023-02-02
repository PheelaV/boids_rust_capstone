boid_data <- flock_detailed(no_iter = 32000,
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
                            dbscan_clustering = T) %>%
  tibble()

print(boid_data %>% filter(id == 0), n = 100)

# roughly 5% of data is from boundaries
nrow(boid_data %>%
       rowwise() %>%
       filter(abs(x) > init_width / 2 - sensory_distance || abs(y) > init_height / 2 - sensory_distance))/ nrow(boid_data)

# gives us original no_iter
nrow(boid_data) / 512 * 32


sensory_distance <- 50
init_width <- 4000
init_height <- 4000
time_dependency <- 2

remove_boundary_data <- function(data, sensory_distance, init_width, init_height, time_dependency) {
  # get the row numbers of boundary data
  boundary_rns <- data %>%
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

  # return(data %>% slice(-affected_rns$rn))
  return(affected_rns$rn)
}

# ~9.42% data affected
1 - nrow(
  boid_data %>% slice(-(boid_data %>% remove_boundary_data(sensory_distance, init_width, init_height, time_dependency)))
  ) /
nrow(boid_data)


deg2rad <- function(deg) deg*pi/180
rad2deg <- function(rad) rad*180/pi


# assumes data is written groups of rows, each representing a single unit, naturally ordered by time
# takes the group and interlaces it with the rest, keeping the time order and thus creating a relay of the simulation as
# the data points have been collected in the first place
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
  boundary_rns <- boid_data %>%
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
    filter(rn <= nrow(boid_data))
  # cleanup
  remove(boundary_rns)

  return(affected_rns$rn)
}

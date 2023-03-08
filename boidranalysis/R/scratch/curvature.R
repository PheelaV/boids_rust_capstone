
install.packages("stokes")  # uncomment this to install the package
library("stokes")
library(zoo)

wedge(c(1,2), c(3,4))

c(1,2) %*% c(3,4)
c(1,2) %^% c(3,4)

N <- config$init_boids

N

# source: garett brown thesis

1/N * wed

D <-  c(-0.3512260083116, 1.7366174855406)
E <-  c(-0.924598070496, 1.5113961783845)
wedge(as.kform(cbind(D)), as.kform(cbind(E)))
wedge(E, D)

my_wedge(D, E)
my_wedge(E, D)


# the formula is



direction_data_by_boid %>%
  group_by(time, cluster_id) %>%
  rowwise() %>%
  mutate(wedge_i = c(
    direction_data_by_boid %>%
      filter(id == .x$id, time >= .x$time && time <= .x$time + tau) %>%
      rowwise() %>%
      select(v = c(.x$dx, .x$dy))
    ))

direction_data_by_boid %>%
  rowwise() %>%
  filter(time %in% time:(time+tau), id == 1) %>%
  mutate(time, ad = lag(time, n = 2, default = 0))

  c <- lag(c(1:8), n = 2, default = 0)

  direction_data_by_boid %>%
    group_by(cluster_id, time) %>%
    mutate(dx_diff = lead(dx) - dx, dy_diff = lead(dy) - dy) %>%
    mutate(num = map2(dx_diff, dy_diff, ~ (.x ^ 2 + .y ^ 2) / 2 ^ 2))




  df <- data.frame(row_num = 1:6, X = c(2, 3, 2, -1, -2, 1))

  df %>%
    mutate(sum_X = rollapply(wedge_i, width = tau, FUN = sum, fill = NA, align = "left"))

  wedges %>%
    rowwise() %>%
    mutate(wtau = sum(wedges %>% filter(time %in% time:(time+tau), id == id) %>% select(wedge_i)))
  test

  my_wedge(c(-125., -7.67 ), c(-27.7, -74.4))
  # %>%
    mutate(num = map2(dx_diff, dy_diff, ~ (.x ^ 2 + .y ^ 2) / 2 ^ 2))

sort(runif(10, 0, 2 * pi))
# // 1/N * SUM_{i = 1}^{N} 1/Tau * (SUM_{t = 1}^{Tau} (v_i(t) %^% v_i(t+1))/v_0^2 )
# where v_0 is described as the constant velocity agents travel at
# or in average_normalized_velocity it is the average velocity of a particular group


# we have a problem as the number of boids in a particular cluster is not a constant value
# if we want to calculate lag values, we will have to choose a specific number across each lag frame to scale the lot by

# or that would be the clean solution, but we can also just add the lag


config <- get_config(
  # "string_ball_hybrid.toml",
  "basic2.toml",
  overwrite = list(
    init_boids = 2^11,
    no_iter = 2^15,
    init_width = 4000,
    init_height = 4000,
    min_speed = 2,
    max_speed = 2,
    sample_rate = 64,
    # boundary_config = "{\"type\": \"Repulsive\", \"distance\": 100, \"force\": 0.05}"
    boundary_config = "{\"type\": \"Toroidal\"}",
    distance_config = "{\"type\": \"EucToroidal\"}"
    # field_of_vision = 360.
  )
)

my_wedge <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  a[1] * b[2] - a[2] * b[1]
}
boid_data <- read_csv("../boids-data_1677829114108.csv")

direction_data_by_boid <- boid_data %>% get_directional_boid_data2()

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

tau <- 3
get_curvature_order_data(direction_data_by_boid, config, tau) %>%
  ggplot(aes(x = time * tau * 8, y = rop)) +
  geom_line()
summary(wedges$wedge_tau)

###############################################
# generating some mock data here
# then measure on that

circle_point <- function(angle, radius) {
  x <- radius * cos(angle)
  y <- radius * sin(angle)
  return(c(x, y))
}

circle_point_wlk <- function(position, walk_distance, radius, turning_point = c(0, 0)) {
  print(position)
  print(class(position))
  local_pos = position - turning_point
  rad_change = walk_distance / radius
  rot_mat = matrix(c(cos(rad_change), sin(rad_change), -sin(rad_change), cos(rad_change)), byrow = T, nrow = 2)
  new_local_pos = rot_mat %*% local_pos
  new_pos = new_local_pos + turning_point
  return(new_pos)
}

circle_point_wlk(c(1, 0), 1, 1)

velocity = .8
radius = 2
steps2 = floor(2 * pi * radius / velocity * (3/4) )
steps = 10
start = c(radius, 0)

test_data <-
  tibble(
    step = seq(0, steps),
    vel = rep(velocity, steps + 1),
    vel2 = seq(velocity, velocity / 8, length.out = steps + 1),
    # vel = rep(velocity, steps + 1),
    x = rep(start[1], steps + 1),
    y = rep(start[2], steps + 1)
  )

for(i in 2:nrow(test_data)) {
  new_pos = circle_point_wlk(c(test_data[i-1,]$x, test_data[i-1,]$y), test_data[i-1,]$vel, radius * (test_data[i-1,]$vel / test_data[1,]$vel))
  test_data[i,]$x = new_pos[1]
  test_data[i,]$y = new_pos[2]
}

diagonal_steps <- 20
for (i in 1:diagonal_steps) {
  test_data <- test_data %>% add_row(x = test_data[nrow(test_data),]$x + velocity/sqrt(2),  y = test_data[nrow(test_data),]$y + velocity/sqrt(2), step = test_data[nrow(test_data),]$step + 1, vel = velocity)
}

data <- test_data %>%
  mutate(id = 1, cluster_id = 1, step = step + 1) %>%
  rename(time = step) %>%
  get_directional_boid_data2

data %>% ggplot(aes(x = time, y = bearings)) +
  geom_point() +
  lims(y = c(-3, 3))

size = 4
ggplot(test_data, aes(x = x, y = y)) +
  lims(x = c(-size, size), y = c(-size, size)) +
  geom_point()

tau = 8
data_outcome <- data %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(lead_dx = lead(dx), lead_dy = lead(dy)) %>%
  rowwise() %>%
  mutate(wedge_i = my_wedge(c(dx, dy), c(lead_dx, lead_dy))) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(wedge_tau = rollapply(wedge_i, width = tau, FUN = function(x)sum(x) / tau, fill = NA, align = "left"))
print(data_outcome, n = 28)

data_outcome %>% ggplot(aes(x = time, y = bearings)) +
  geom_line()
# %>%
#   select(-lead_dx, -lead_dy) %>%
#   filter(!is.na(wedge_tau)) %>%
#   select(time, cluster_id, id, wedge_tau) %>%
#   group_by(time) %>%
#   reframe(rop = mean(wedge_tau))

summary(direction_data_by_boid$bearings)


# this would be the -pi/+pi turns of an individual

direction_data_by_boid %>%
  filter(id == 1) %>%
  mutate(ma2=rollapply(bearings,5,mean,align='right',fill=NA)) %>%
  ggplot(aes(x = time, y = ma2)) +
  geom_line(aes(color = id))

direction_data_by_boid

boid_data <- read_csv("../boids-data_1677829114108.csv")

direction_data_by_boid <- boid_data %>% get_directional_boid_data2()

# mean overall
direction_data_by_boid %>%
  filter(id == 1) %>%
  arrange(time) %>%
  group_by(id) %>%
  mutate(breaings_ma=rollapply(bearings, 5,mean,align='right',fill=NA)) %>%
  ungroup() %>%
  mutate(noise = cluster_id == 0) %>%
  # filter(cluster_id != 0) %>%
  group_by(time, noise) %>%
  summarise(breaings_ma = mean(breaings_ma)) %>%
  ungroup() %>%
  filter(!is.na(breaings_ma)) %>%
  ggplot(aes(x = time, y = breaings_ma)) +
  geom_point(size = 1.2, alpha = .3, aes(color = noise))  +
  stat_summary(aes(y = breaings_ma, color = noise), alpha = .8, fun=mean, geom="line") +
  theme_bw()


# noise vs no noise
direction_data_by_boid %>%
  filter(id == 1) %>%
  arrange(time) %>%
  group_by(id) %>%
  mutate(breaings_ma=rollapply(bearings, 5,mean,align='right',fill=NA)) %>%
  ungroup() %>%
  group_by(time, cluster_id) %>%
  summarise(breaings_ma = mean(breaings_ma)) %>%
  ungroup() %>%
  filter(!is.na(breaings_ma)) %>%
  mutate(noise = cluster_id == 0) %>%
  ggplot(aes(x = time, y = breaings_ma)) +
  geom_point(size = .5, alpha = .5, aes(color = noise))  +
  stat_summary(aes(y = breaings_ma, color = noise), fun=mean, geom="line") +
  theme_bw()

# abs noise vs no noise
direction_data_by_boid %>%
  # filter(id == 1) %>%
  arrange(time) %>%
  group_by(id) %>%
  mutate(bearings = abs(bearings)) %>%
  mutate(ma2=rollapply(bearings,5,mean,align='right',fill=NA)) %>%
  ungroup() %>%
  group_by(time, cluster_id) %>%
  summarise(ma2 = mean(ma2)) %>%
  ungroup() %>%
  filter(!is.na(ma2)) %>%
  # filter(cluster_id != 0) %>%
  # filter(cluster_id == 1) %>%
  mutate(noise = cluster_id == 0) %>%
  ggplot(aes(x = time, y = ma2)) +
  geom_point(size = .5, alpha = .5, aes(color = noise))  +
  stat_summary(aes(y = ma2, color = noise), fun=mean, geom="line")



# noise vs no noise
direction_data_by_boid %>%
  filter(id == 1) %>%
  arrange(time) %>%
  group_by(id) %>%
  mutate(ma2=rollapply(bearings,10,mean,align='right',fill=NA)) %>%
  ungroup() %>%
  group_by(time, cluster_id) %>%
  summarise(ma2 = mean(ma2)) %>%
  ungroup() %>%
  filter(!is.na(ma2)) %>%
  # filter(cluster_id != 0) %>%
  # filter(cluster_id == 1) %>%
  mutate(noise = cluster_id == 0) %>%
  ggplot(aes(x = time, y = ma2)) +
  geom_point(size = .5, alpha = .5, aes(color = cluster_id))  +
  stat_summary(aes(y = ma2, color = cluster_id), fun=mean, geom="line")

table(direction_data_by_boid$id)



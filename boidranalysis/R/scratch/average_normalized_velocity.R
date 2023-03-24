config <- get_config(
  # "string_ball_hybrid.toml",
  "basic.toml",
  overwrite = list(
    init_boids = 2^9,
    no_iter = 2^12,
    init_width = 4000,
    init_height = 4000,
    sample_rate = 32,
    # boundary_config = "{\"type\": \"Repulsive\", \"distance\": 100, \"force\": 0.05}"
    boundary_config = "{\"type\": \"Toroidal\"}"
  )
)
boid_data <- rlang::exec(flock_detailed, !!!config) %>%
  tibble()

boid_data %>%
  group_by(id) %>%
  reframe(
    dx = diff(x, 1),
    dy = diff(y, 1),
    id = id[2:n()],
    x = x[2:n()],
    y = y[2:n()],
    cluster_id = cluster_id[2:n()],
    n_neighbours = n_neighbours[2:n()],
    time = time[2:n()]
    )

direction_data_by_boid <- boid_data %>%
  get_directional_boid_data()

summary(direction_data_by_boid)


# sample calculation
tibble(
  dx = c(1, 2, -3, -4),
  dy = c(1, -2, -3, 4)
) %>%
  summarise(
    v0x = sum(abs(dx)) / n(),
    v0y = sum(abs(dy)) / n(),
    vix = abs(sum(dx)),
    viy = abs(sum(dy)),
    count = n()
    ) %>%
  mutate(result = sqrt(vix^2 + viy^2) / (count * sqrt(v0x^2 + v0y^2)))

direction_data_by_boid %>%
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
  mutate(anv = sqrt(vix^2 + viy^2) / (count * sqrt(v0x^2 + v0y^2))) %>%
  mutate(no_flocks_norm = (no_flocks - min(no_flocks)) / (max(no_flocks) - min(no_flocks))) %>%
  # mutate(no_flocks_norm = scale(no_flocks, center = F)) %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = anv, color = "anv")) +
  geom_line(aes(y = no_flocks_norm, color = "flocks"))
  # geom_line(aes(x = time, group = cluster_id, color = cluster_id))

direction_data_by_boid %>%
  filter(cluster_id != 0) %>%
  select(dx, dy, time, cluster_id) %>%
  group_by(time) %>%
  group_by(time, cluster_id) %>%
  summarise(
    v0x = sum(abs(dx)) / n(),
    v0y = sum(abs(dy)) / n(),
    vix = abs(sum(dx)),
    viy = abs(sum(dy)),
    count = n()
  ) %>%
  mutate(result = sqrt(vix^2 + viy^2) / (count * sqrt(v0x^2 + v0y^2))) %>%
  ggplot(aes(y = result)) +
geom_point(aes(x = time, group = cluster_id, color = time))




get_average_norm_vel(direction_data_by_boid)


boid_data <- read_csv("./Data/0301_experiment2_a12/boids-data_1678117747179.csv")

direction_df <- boid_data |>
  get_directional_boid_data2()

direction_counts <- direction_df |>
  get_directional_counts(no_bins = 1440)



direction_counts %>% ggplot(aes(x = bin, y = heading_count)) +
  theme_bw() +
  geom_col( color = "salmon") +
  geom_line(color = "blue")


direction_df %>% ggplot(aes(headings)) +
  theme_bw() +
  geom_density()

direction_df %>% ggplot(aes(bearings)) +
  theme_bw() +
  geom_density()


library(zoo)

ma <- 40
direction_df |>
  group_by(id) |>
  arrange(time) |>
  reframe(
    bearing_ma = rollmean(bearings, ma),
    heading_ma = rollmean(headings, ma)
    ) |>
  # filter(id == 0) |>
  ggplot(
    # aes(group = id)
  ) +
  geom_density(aes(scale(bearing_ma)), color = "green") +
  geom_density(aes(scale(heading_ma)), color = "red")
  # coord_polar(start = 2*pi/3.3)

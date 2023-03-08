max_speed <- 10
direction_data_by_boid <- preprocess(boid_data)

no_individuals <- max(direction_data_by_boid$id) + 1

direction_data_by_boid %>%
  # filter(id %in% c(160:179)) %>%
  filter(time >= 100 & time < 200) %>%
  mutate(breaking = sqrt((x - lag(x, default = 0))^2 + (y - lag(y, default = 0))^2) > max_speed) %>%
  mutate(breaking = if_else(time == min(time), F, breaking)) %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(breakage_agent = cumsum(breaking)) %>%
  ungroup() %>%
  mutate(break_global = if_else(breaking, cumsum(breaking), NA)) %>%
  group_by(id) %>%
  fill(break_global) %>%
  ungroup() %>%
  mutate(group_id = if_else(is.na(break_global), id, break_global + no_individuals)) %>%
  ggplot() +
  geom_point(aes(x = x, y = y, color = factor(id %% 10), group = group_id)) +
  scale_color_brewer(palette = "Paired")


tibble(
  x = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2),
  y = c(0, 1, 0, 0, 0, 0, 4, 3, 0, 0)
  ) %>% tidyr::fill(y)

table(test$breaking)
table(boid_data$time)


library(RColorBrewer)

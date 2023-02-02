
direction_data_by_boid <- boid_data %>%
  group_by(id) %>%
  summarise(
    headings = get_headings(x, y)
  ) %>%
  summarise(
    bearings = get_bearings(headings),
    # as bearings are calculated from pairs, to map them one to one, throw away the headings of t_0 to match lengths
    headings = headings[2:length(headings)]
  ) %>%
  ungroup(id)
postfix <- "_before.jpg"

sensory_distance <- 50
init_width <- 4000
init_height <- 4000
time_dependency <- 2

direction_data_by_boid <- boid_data %>%
  group_by(id) %>%
  summarise(
    headings = get_headings(x, y)
  ) %>%
  summarise(
    bearings = get_bearings(headings),
    # as bearings are calculated from pairs, to map them one to one, throw away the headings of t_0 to match lengths
    headings = headings[2:length(headings)]
  ) %>%
  ungroup(id) %>%
  slice(-tail(
    remove_boundary_data(boid_data, sensory_distance, init_width, init_height, time_dependency),
    -2
  ))
postfix <- "_after.jpg"

# we have a time dependency of two - 3rd records becomes the first as it
# merges information from the previous two
# they are already removed and need to be ignored to reach parity


print("range of bearings:")
range(direction_data_by_boid$bearings)

print("range of headings:")
range(direction_data_by_boid$headings)

headings_binned <- direction_data_by_boid %>%
  mutate(
    headings = floor(rad2deg(headings))
  ) %>%
  group_by(headings) %>%
  summarise(
    count = length(id)
  ) %>%
  select(heading = headings, count = count)

bearings_binned <- direction_data_by_boid %>%
  mutate(
    bearings = floor(rad2deg(bearings))
  ) %>%
  group_by(bearings) %>%
  summarise(
    count = length(id)
  ) %>%
  select(bearing = bearings, count = count)

ggplot(headings_binned, aes(x = heading, y = log(count))) +
  coord_polar(theta = "x", start = pi) +
  geom_bar(stat = "identity", fill = "deeppink4", width = .9) +
  # geom_hline(yintercept = seq(0, 500, by = 100), color = "grey80", size = 0.3) +
  # scale_x_continuous(breaks = 0:24, expand = c(.002,0)) +
  labs(x = "Heading w.r.t N") +
  theme_bw()
ggsave(paste0("headings_binned", postfix), path = "./Plots", width = 6, height = 6, units = "in")

summary(direction_data_by_boid %>% select(-id))

ggplot(bearings_binned, aes(x = bearing, y = log10(count))) +
  coord_polar(theta = "x", start = pi) +
  geom_bar(stat = "identity", fill = "deeppink4", width = .9) +
  # geom_hline(yintercept = seq(0, 500, by = 100), color = "grey80", size = 0.3) +
  # scale_x_continuous(breaks = 0:24, expand = c(.002,0)) +
  labs(x = "Bearing w.r.t N") +
  theme_bw()
ggsave(paste0("bearings_binned", postfix), path = "./Plots", width = 6, height = 6, units = "in")


ggplot(direction_data_by_boid, aes(x=headings)) +
  # coord_polar(theta = "x", start = pi) +
  geom_density(color="aquamarine4", fill="aquamarine3") +
  theme_bw()
ggsave(paste0("headings_flat_binned", postfix), path = "./Plots", width = 6, height = 6, units = "in")


ggplot(direction_data_by_boid, aes(x=bearings)) +
  geom_density(color="aquamarine4", fill="aquamarine3") +
  theme_bw()
ggsave(paste0("bearings_flat_binned", postfix), path = "./Plots", width = 6, height = 6, units = "in")




x <- c(0, 2, 3, 5, 6)
y <- c(0, 2, 5, 6, 3)

# this
x[2:length(x)] - x[1:length(x) - 1]
# can be done more easily as this
x
y
diff(x, 1)
diff(y, 1)
#

get_headings <- function(x, y, radians = F) {
  if(radians) {
    mapply(function(x, y) atan2(x, y), diff(x, 1), diff(y,1))
  } else {
    mapply(function(x, y) atan2(x, y) * 180/pi, diff(x, 1), diff(y,1))

  }
}

get_headings(x, y)
get_headings(x, y, Trues) / pi

atan2(1 / 2, sqrt(3) / 2) * 180 / pi
atan2(0.8, -1.4) * 180 / pi


scalar1 <- function(x) {x / sqrt(sum(x^2))}

test <- scalar1(c(2,2))

atan2(test[2], test[1]) * 180/pi
90 - atan2(2, 2) * 180/pi
90 - atan2(5 - 2, 3 - 2) * 180/pi
90 - atan2(6 - 5, 5 - 3) * 180/pi

atan2(3, 1) * 180/pi
atan2(1, 3) * 180/pi

library(ggplot2)
hist(get_headings(x, y))

hist(get_headings(x, y, T) / pi)



setwd("/Users/filipvlcek/Source/Repos/boids_rust")
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

heading_data <- readr::read_csv("test.csv")
headings_by_boid <- heading_data %>%
  group_by(id) %>%
  summarise(
    headings = get_headings(x, y)
  ) %>%
  ungroup(id)

headings_binned <- headings_by_boid %>%
  mutate(
    headings = round(headings, 0)
  ) %>%
  group_by(headings) %>%
  summarise(
    count = length(id)
  ) %>%
  select(heading = headings, count = count)

ggplot(headings_by_boid, aes(x=headings)) +
  geom_density(color="aquamarine4", fill="aquamarine3") +
  theme_bw()

ggplot(headings_binned, aes(x = heading, y = log10(count))) +
  coord_polar(theta = "x", start = pi) +
  geom_bar(stat = "identity", fill = "deeppink4", width = .9) +
  # geom_hline(yintercept = seq(0, 500, by = 100), color = "grey80", size = 0.3) +
  # scale_x_continuous(breaks = 0:24, expand = c(.002,0)) +
  labs(x = "Heading w.r.t N") +
  theme_bw()

install.packages("gganimate")

headings_binned_by_time <-
  headings_by_boid %>%
    mutate(
      headings = round(headings, 0)
    ) %>%
    group_by(id) %>%
    mutate(t = 1:n()) %>%
    # select(t, headings) %>%
    group_by(t) %>%
    group_by(headings) %>%
    mutate(
      count = n()
    ) %>%
    ungroup() %>%
    select(-id)

library(gganimate)

headings_binned_by_time %>% filter(t < 5)

# how many time points do I have
headings_binned_by_time %>% select(t) %>% summarise(unique = n_distinct(t))
# how many seconds of recording do I have
headings_binned_by_time %>% select(t) %>% summarise(unique = n_distinct(t) / 60)


headings_binned_by_time_t <- headings_binned_by_time %>% mutate(count = log10(count))

p <- ggplot(headings_binned_by_time_t %>% filter(t < 562), aes(x = headings, y = count)) +
  coord_polar(theta = "x", start = pi) +
  geom_bar(stat = "identity", fill = "maroon4", width = .9) +
  labs(x = "Heading w.r.t N", y="log10(count)", title = "Boids headings 10s") +
  theme_bw() +
  transition_states(
    t,
    transition_length = 1/120,
    state_length = 1/60
  )
  # enter_fade() +
  # exit_shrink() +
  # ease_aes('sine-in-out')

animate(p, renderer = ffmpeg_renderer(), fps=30, duration = 10)


p  %>% animate()


rendered <- file_renderer(dir = ".", prefix = "gganim_plot", overwrite = TRUE)


animate(p, renderer = rendered)

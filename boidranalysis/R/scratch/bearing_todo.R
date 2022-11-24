source("helpers.R")
library(tidyr)
library(dplyr)
library(ggplot2)



df <- tibble(
    deg =  seq(0, 359, 1),
    x = cos(deg2rad(deg)),
    y = sin(deg2rad(deg)),
  ) %>%
  mutate(
    atan2_vals = atan2(y, x),
    atan2_vals_roudned = round(atan2(y, x), 13)
  )

ggplot(df) +
  geom_col(aes(x = deg, y = abs(atan2_vals_roudned))) +
  scale_x_continuous(
    breaks = seq(0, 359, 45),
    minor_breaks = seq(0, 359, 15)
  ) +
  coord_polar() +
  theme_bw()

df <- df %>%
  mutate(
    pathX = cumsum(x),
    pathY = cumsum(y)
  )

ggplot(df, aes(pathX, pathY)) +
  geom_point()

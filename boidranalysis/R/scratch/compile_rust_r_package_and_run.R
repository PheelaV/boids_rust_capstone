usethis::create_package("/Users/filipvlcek/Source/scratch_code/boidr")

getwd() # /Users/filipvlcek/Source/scratch_code/boidr
# rextendr::use_extendr()

rextendr::document()
devtools::load_all(".")
hello_world()

get_headings <- function(x, y, radians = F) {
  if(radians) {
    mapply(function(x, y) atan2(x, y), diff(x, 1), diff(y,1))
  } else {
    mapply(function(x, y) atan2(x, y) * 180/pi, diff(x, 1), diff(y,1))

  }
}

plot_heading_data <- function(heading_data) {
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

  # ggplot(headings_binned, aes(x = heading, y = log10(count))) +
  #   coord_polar(theta = "x", start = pi) +
  #   geom_bar(stat = "identity", fill = "deeppink4", width = .9) +
  #   labs(x = "Heading w.r.t N") +
  #   theme_bw()
}


library(ggplot2)
library(tidyr)
library(dplyr)

data<- hello_world2(no_iter = 100)
plot_heading_data(as_tibble(data))

data<- hello_world2(no_iter = 1000)
plot_heading_data(as_tibble(data))

data<- hello_world2(no_iter = 10000)
plot_heading_data(as_tibble(data))

data<- hello_world2(no_iter = 10000)
plot_heading_data(as_tibble(data))

data<- hello_world2(no_iter = 100000)
plot_heading_data(as_tibble(data))

data<- hello_world2(no_iter = 1000000)
plot_heading_data(as_tibble(data))

# data<- hello_world2(no_iter = 10000000)
# data<- hello_world2(no_iter = 100000000)

hist(get_headings(data$x, data$y))

heading_data <- as_tibble(data)



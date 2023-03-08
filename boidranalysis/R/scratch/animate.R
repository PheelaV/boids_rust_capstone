library(gganimate)
#> Loading required package: ggplot2


# We'll start with a static plot
p <- ggplot(iris, aes(x = Petal.Width, y = Petal.Length)) +
  geom_point()

plot(p)

anim <- p +
  transition_states(Species,
                    transition_length = 2,
                    state_length = 1)

anim

data_update_rate <- 120
clip_duration = max(boid_data$time) / data_update_rate

max(boid_data$time)
min(boid_data$time)

p <- ggplot(
  boid_data,
  aes(x = x, y=y, colour = time)
) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  scale_colour_gradientn(colours=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))


anim <- p +
  transition_states(time,
                    transition_length = 2,
                    state_length = 1)
# key_frames <- (boid_data %>%
#   select(time) %>%
#   distinct() %>%
#   summarise(no = n()))$no

animate(anim, renderer = av_renderer('animation.mp4'), width = 1920, height = 1080, fps = 30, duration = clip_duration)


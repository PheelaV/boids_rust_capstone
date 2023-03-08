p <- ggplot(
  boid_data,
  aes(x = x, y=y, colour = time)
) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  scale_colour_gradientn(colours=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))


p <- p + transition_time(3) +
  labs(title = "timestep: {frame_time}")

rendered <- file_renderer(dir = ".", prefix = "gganim_plot", overwrite = TRUE)

animate(p)

animate(p, renderer = rendered)

boid_data$timetime_steps <- 3000
time_steps / 120

ffmpeg_renderer(height = 1080, width = 1080)
animate(p, renderer = ffmpeg_renderer(), fps=30, duration = 33)
round(time_steps / 120)
q <- 2

animate(p, renderer = av_renderer('animation.mp4'), width = 720*q, height = 480*q, res = 72*q, fps = 30, duration = 5)
utils::browseURL('animation.mp4')
ffmpeg_renderer(
  gg = "auto",
  ffmpeg = NULL,
  options = list(pix_fmt = "yuv420p")
)


ggplot(
  nth_time_point(all, 100),
  aes(x = x, y=y, colour = cluster)) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  scale_colour_gradientn(colours=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

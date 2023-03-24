library(RColorBrewer)
ggplot(
  data = eta_results_stats,
  aes(x = t, y = mean_no_flocks, color = factor(eta))
) +
  geom_line() +
  # geom_smooth(method = "gam", level = .9) +
  theme_bw() +
  labs(x = "time t", y = "mean_t mean # flocks",
       title = paste(
         "mean number of flocks;",
         "agents:", config$init_boids, ";",
         "iterations:", config$no_iter, ";",
         "sample:", config$sample_rate, ";"
       )
       # subtitle = figure_subtitle_config
  )
# ggsave(paste0("mean_no_flocks", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

ggplot(
  data = eta_results_stats |> arrange(desc(eta)),
  aes(x = t, y = var_no_flocks, color = factor(eta))
) +
  geom_line() +
  # geom_smooth(method = "gam", level = .9) +
  theme_bw() +
  labs(x = "time t", y = "var_t mean # flocks",
       title = "var number of flocks"
       # subtitle = figure_subtitle_config
  )

# ggsave(paste0("var_no_flocks", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

ggplot(
  data = eta_results_stats,
  aes(x = t, y = mean_convex_hull, color = factor(eta))
) +
  geom_line() +
  # geom_smooth(method = "gam", level = .9) +
  theme_bw() +
  labs(x = "time t", y = "mean_t mean flock area",
       title = "mean flock area",
       # subtitle = figure_subtitle_config
  )
# ggsave(paste0("mean_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

ggplot(
  data = eta_results_stats,
    # filter ((eta * 100) %% 5 == 0),
  aes(x = t, y = var_convex_hull, color = eta)
) +
  geom_line(alpha = 0.8) +
  # geom_smooth(method = "gam", level = .9) +
  theme_bw() +
  labs(x = "time t", y = "var_t mean flock area",
       title = "var flock area",
       # subtitle = figure_subtitle_config
  ) +
  scale_colour_gradientn(colours = (brewer.pal(9, 'YlGn')))
# ggsave(paste0("var_area_flock", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)


ggplot(
  data = eta_results_stats |> arrange(desc(eta)),
  aes(x = t, y = mean_average_norm_vel, color = eta, alpha = .8)
) +
  geom_line() +
  theme_bw() +
  labs(x = "time t", y = "mean average norm vel",
       title = "average norm vel",
       # subtitle = figure_subtitle_config
  ) + scale_colour_gradientn(colours = (brewer.pal(9, 'YlGn')))
# ggsave(paste0("mean_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)


ggplot(
  data = eta_results_stats,
  aes(x = t, y = var_average_norm_vel, color = eta)
) +
  geom_line() +
  theme_bw() +
  labs(x = "time t", y = "var average norm vel",
       title = "var average norm vel",
       # subtitle = figure_subtitle_config
  ) + scale_colour_gradientn(colours = (brewer.pal(9, 'YlGn')))
# ggsave(paste0("var_avg_norm_vel", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

ggplot(
  data = eta_results_stats |> arrange(desc(eta)),
  aes(x = t, y = mean_voronoi_counts, color = eta)
) +
  geom_line() +
  theme_bw() +
  labs(x = "time t", y = "mean Count(voronoi cells area < pi*R^2)",
       title = "average voronoi cell area",
       # subtitle = figure_subtitle_config
  ) + scale_colour_gradientn(colours = (brewer.pal(9, 'YlGn')))
# ggsave(paste0("mean_count_voronoi_cell_area", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

ggplot(
  data = eta_results_stats |> arrange(desc(eta)),
  aes(x = t, y = var_voronoi_counts, color = eta)
) +
  geom_line() +
  theme_bw() +
  labs(x = "time t", y = "var Count(voronoi cells area < pi*R^2)",
       title = "var voronoi cell area",
       # subtitle = figure_subtitle_config
  ) + scale_colour_gradientn(colours = (brewer.pal(9, 'YlGn')))
# ggsave(paste0("var_count_voronoi_cell_area", figure_postfix), path = experiment_plot_folder, width = figure_width, height = figure_height, units = figure_units, dpi = figure_dpi)

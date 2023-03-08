library(ggplot2)

# Create a sample dataset
df <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  config_param1 = rep(c("A", "B"), each = 50),
  config_param2 = rep(c("C", "D"), times = 50)
)

# Create a facetted plot with configuration parameters as the categorical variable
ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  facet_grid(config_param1 ~ config_param2) +
  theme_gray()

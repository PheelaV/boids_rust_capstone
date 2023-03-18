my_wedge(c(0.5, -0.5), c(0.5, 0.5))


my_wedge(c(0.5, -0.5), c(0.5, 0.5)) / 0.5^2


my_wedge <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  a[1] * b[2] - a[2] * b[1]
}

expand_grid(x = seq(-3, 3, 0.5), y = seq(-3, 3, 0.5)) %>%
  rowwise() %>%
  mutate(
    sq = if_else(x == 0 | y == 0, max(abs(x), abs(y)) ^ 2, x * y),
    max_sq = max(abs(x), abs(y)) ^ 2
  ) %>%
  mutate(
    w = my_wedge(c(x, y), c(-y, x))
  ) %>%
  ungroup() %>%
  mutate(
    w_sq = w / sq,
    w_max_sq = w / max_sq
  ) %>%
  mutate()
%>% print(n = 120)


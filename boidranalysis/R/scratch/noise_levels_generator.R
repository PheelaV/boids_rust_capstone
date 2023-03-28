tibble(noises = c(seq(from = 0, to = 0.25 - 0.01, by = 0.25 / 2),
                  seq(from = 0.25, to = 0.375 - 0.01, by = 0.5 / 17),
                  seq(from = 0.375, to = 0.625 - 0.01, by = 0.25 / 15),
                  seq(from = 0.625, to = 0.75 - 0.01, by = 0.5 / 17),
                  seq(from = 0.75, to = 1, by = 0.25 / 2))) |>
  mutate(row = row_number()) |>
  ggplot(aes(x = row, y = noises)) + geom_point()


seq(from = 0.375, to = 0.625, by = 0.25 / 15)

##true, true, true) => &
tibble(p = list(
  c(1, -1), c(1, 0), c(1, 1), c(2, 0), c(2, 1),
  c(0, 1), c(2, 2), c(0, -1), c(2, -1), c(1, 2)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##true, true, false) => &
tibble(p = list(
  c(1, 1), c(0, 1), c(-1, 1), c(0, 2), c(1, 2),
  c(1, 0), c(2, 2), c(-1, 0), c(2, 1), c(-1, 2)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##false, true, false) => &
tibble(p = list(
  c(1, 1), c(0, 1), c(-1, 1), c(0, 2), c(-1, 2),
  c(-1, 0), c(-2, 2), c(1, 0), c(1, 2), c(-2, 1)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##false, true, true) => &
tibble(p = list(
  c(-1, 1), c(-1, 0), c(-1, -1), c(-2, 0), c(-2, 1),
  c(0, 1), c(-2, 2), c(0, -1), c(-1, 2), c(-2, -1)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##false, false, true) => &
tibble(p = list(
  c(-1, 1), c(-1, 0), c(-1, -1), c(-2, 0), c(-2, -1),
  c(0, -1), c(-2, -2), c(0, 1), c(-2, 1), c(-1, -2)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##false, false, false) => &
tibble(p = list(
  c(-1, -1), c(0, -1), c(1, -1), c(0, -2), c(-1, -2),
  c(-1, 0), c(-2, -2), c(1, 0), c(-2, -1), c(1, -2)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##true, false, false) => &
tibble(p = list(
  c(-1, -1), c(0, -1), c(1, -1), c(0, -2), c(1, -2),
  c(1, 0), c(2, -2), c(-1, 0), c(-1, -2), c(2, -1)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


##true, false, true) => &
tibble(p = list(
  c(1, -1), c(1, 0), c(1, 1), c(2, 0), c(2, -1),
  c(0, -1), c(2, -2), c(0, 1), c(1, -2), c(2, 1)
))  %>% rowwise() %>% mutate (x = p[[1]], y = p[[2]]) %>% ggplot(aes(x = x, y = y)) + geom_point() +  theme(axis.line = element_line()) +   coord_axes_inside(labels_inside = TRUE) + lims(x = c(-2, 2), y = c(-2, 2))


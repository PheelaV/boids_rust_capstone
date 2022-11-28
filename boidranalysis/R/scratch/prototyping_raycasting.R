scalar1 <- function(x) {x / sqrt(sum(x^2))}
sight_sample = c(1, 2, 3, 4, 5, 6)
n = seq(0, 2)

angle_change_data <- tibble(n0 = numeric(), n1= numeric(), n2 = numeric())
for(ss in sight_sample) {
  angle_change = n / ss - 1 / ss
  angle_change_data <- angle_change_data %>%
    add_row(
      n0 = angle_change[1],
      n1 = angle_change[2],
      n2 = angle_change[3],
    )
}

# theta = deg2rad(90)
theta = 2.3 / 2
rot_mat = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
boid_heading <- c(0,1)

test_data <- tibble(angle_change = c(t(angle_change_data[3,]))) %>%
  rowwise() %>%
  mutate(v1 = list( scalar1(boid_heading + c((rot_mat * angle_change) %*% boid_heading)))) %>%
  select(
    angle_change, v1
  ) %>%
  rowwise() %>%
  mutate(v1_x = v1[1], v1_y = v1[2]) %>%
  mutate(colour = factor(angle_change))

ggplot(test_data, aes(x = v1_x, y = v1_y)) +
  geom_point(colour = "blue") +
  # lims(x = c(-6, 6), y = c(-6, 6)) +
  geom_segment(aes(x = 0, y = 0, xend = v1_x, yend = v1_y, colour = colour))
  # geom_segment(aes(x = 0, y = 0, xend = boid_heading[1], yend = boid_heading[2], colour = "boid heading"))

tibble(x = seq(1, 6) / 3) %>%
  rowwise() %>%
  mutate(y = list((rot_mat * x) %*% c(1,1)))

test_data <- tibble(v1_x = c(3, 6, -5), v1_y = c(3, -4, 2), v1 = list(c(3, 3), c(6, -4), c(-5, 2)))
test_data

test_data <-
  test_data %>%
  rowwise() %>%
  mutate(v1_x = v1[1], v1_y = v1[2]) %>%
  # as_tibble() %>%
  mutate(v2 = list(rot_mat %*% v1)) %>%
  mutate(v2_x = v2[1], v2_y = v2[2]) %>%
  mutate(v3 = list((rot_mat * 1/3) %*% v1)) %>%
  mutate(v3_x = v3[1], v3_y = v3[2])

ggplot(test_data, aes(x = v1_x, y = v1_y)) +
  geom_point(colour = "blue") +
  lims(x = c(-6, 6), y = c(-6, 6)) +
  geom_segment(aes(x = 0, y = 0, xend = v1_x, yend = v1_y, colour = "point"))+
  geom_segment(aes(x = 0, y = 0, xend = v2_x, yend = v2_y, colour = "rot")) +
  geom_segment(aes(x = 0, y = 0, xend = v3_x, yend = v3_y, colour = "rot scaled"))


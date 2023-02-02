direction_data_by_boid$headings

range(direction_data_by_boid$headings)
range(direction_data_by_boid$bearings)
library(circular)
mle.wrappedcauchy(dire)



boid_data %>%
  group_by(id) %>%
  summarise(n = n())

direction_data_by_boid %>%
  group_by(id) %>%
  summarise(n = n())


test <- 1:8
test
test[2:8]


boid_data %>%
  group_by(id) %>%
  summarise(
    headings = get_headings(x, y)
  ) %>%
  summarise(
    bearings = get_bearings(headings),
    # as bearings are calculated from pairs, to map them one to one, throw away the headings of t_0 to match lengths
    headings = headings[2:length(headings)]
  ) %>%
  ungroup(id)

boid_data %>%
  group_by(id) %>%
  summarise(
    headings = get_headings(x, y)
  ) %>%
  group_by(id) %>%
  summarise(n = n())


nrow(boid_data)
nrow(direction_data_by_boid)

nrow(boid_data) - nrow(direction_data_by_boid)
boid_data %>%
  group_by(id) %>%
  slice(3:n()) %>% bind_cols(
    direction_data_by_boid %>% select(headings, bearings)
  )

test_boid_data <- test %>%
  group_by(id) %>%
  summarise(
    headings = get_headings(x, y)
  ) %>%
  summarise(
    bearings = get_bearings(headings),
    # as bearings are calculated from pairs, to map them one to one, throw away the headings of t_0 to match lengths
    headings = headings[2:length(headings)]
  ) %>%
  ungroup(id)

range(test_boid_data$bearings)
range(test_boid_data$headings)

max(test_boid_data$headings)

temp_dat <- boid_data
boid_data <- test_boid_data

########## actual fitting ########## problem while fitting normal, also is it always so shouty?
fit_cauchy <- mle.wrappedcauchy(direction_data_by_boid$bearings, max.iter =  2000)

fit_normal <- mle.wrappednormal(direction_data_by_boid$bearings, max.iter =  2000)

fit_data <- tibble(x = seq(from= 0, to = 2*pi, by = .001)) %>%
  mutate(
    cauchy = dwrappedcauchy(x, mu = fit_cauchy$mu, rho = fit_cauchy$rho),
    norm = dwrappednormal(x, mu = fit_normal$mu, rho = fit_normal$rho, sd = fit_normal$sd)) %>%
ggplot(aes(x = x)) +
  geom_point(aes(y = cauchy, color = "cauchy")) +
  geom_point(aes(y = norm, color = "norm"))

tibble(
  x = seq(from= 0, to = 2*pi, by = .001),
  test_cauchy = dwrappedcauchy(circular(seq(from= 0, to = 2*pi, by = .001)), mu = fit_cauchy$mu, rho = fit_cauchy$rho),
  test_normal = dwrappednormal(circular(seq(from= 0, to = 2*pi, by = .001)), mu = fit_normal$mu, rho = fit_normal$rho, sd = fit_normal$sd)
) %>%
  ggplot(aes(x = x)) +
  geom_point(aes(y = test_cauchy, color = "cauchy")) +
  geom_point(aes(y = test_normal, color = "norm"))
# then the error from pdf can be calculated as L2

sum(bearings_binned$count)
bearings_binned %>%
  mutate(bearing = deg2rad(bearing + 180)) %>%
  mutate(theo_cauchy = dwrappedcauchy(bearing, mu = fit_cauchy$mu, rho = fit_cauchy$rho)) %>%
  mutate(count = count / 510976) %>%
  ggplot(aes(x = bearing)) +
  geom_line(aes(y = count, color = "real")) +
  geom_line(aes(y = theo_cauchy, color = "cauchy fit"))

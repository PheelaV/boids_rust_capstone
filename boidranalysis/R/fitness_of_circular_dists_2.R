library(circular)
library(CircStats)

# generate random wrapped-normal data
random_data <- (rwrpnorm(10000,0,0.8)+pi)%%(2*pi)-pi

# get the density
rd_density <- density(random_data)
# get the function approximation
rd_approx<- approxfun(rd_density$x, rd_density$y)

# fit a wrapped normal
mle <-mle.wrappednormal(random_data)

# define value-range
theta <- seq(-pi,pi,0.01)

ddd <- rd_density(theta)
ddd[which(is.na(rd_density(theta)))] <- 0
mu <- mle$mu
rho <- mle$rho
density <- theta
for(i in 1:length(density)) density[i] <-  dwrpnorm(theta[i], mu, rho)

library(ggplot2)
library(tibble)
ggplot(tibble(theta, density)) +
  geom_line(aes(x = theta, y = density))

ggplot(tibble(theta, density, fitted = ddd)) +
  geom_line(aes(x = theta, y = density, colour = "a")) +
  geom_line(aes(x = theta, y = fitted, colour = "b"))

dwrpnorm(seq(-pi,pi,0.01),mle.wrappednormal(random_data)$mu,mle.wrappednormal(random_data)$rho)
lines(seq(-pi,pi,0.01),dwrpnorm(seq(-pi,pi,0.01),mle.wrappednormal(random_data)$mu,mle.wrappednormal(random_data)$rho),col="red")
lines(seq(-pi,pi,0.01),dwrpcauchy(seq(-pi,pi,0.01),mle.wrappedcauchy(random_data)$mu,mle.wrappedcauchy(random_data)$rho),col="green")



sum((ddd-dwrpnorm(seq(-pi,pi,0.01),mle.wrappednormal(TA)$mu,mle.wrappednormal(TA)$rho))^2)
sum((ddd-dwrpnorm(seq(-pi,pi,0.01),mle.wrappedcauchy(TA)$mu,mle.wrappedcauchy(TA)$rho))^2)

sum((ddd - density)^2)

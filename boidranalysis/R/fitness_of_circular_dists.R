library(circular)
library(CircStats)

# generate random wrapped-normal data
TA <- (rwrpnorm(10000,0,0.8)+pi)%%(2*pi)-pi

# get the density
d <- density(TA)
# get the function approximation
dd <- approxfun(d$x, d$y)

# fit a wrapped normal
mle <-mle.wrappednormal(TA)

# define value-range
theta <- seq(-pi,pi,0.01)

ddd <- dd(theta)
ddd[which(is.na(dd(theta)))] <- 0
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

dwrpnorm(seq(-pi,pi,0.01),mle.wrappednormal(TA)$mu,mle.wrappednormal(TA)$rho)
lines(seq(-pi,pi,0.01),dwrpnorm(seq(-pi,pi,0.01),mle.wrappednormal(TA)$mu,mle.wrappednormal(TA)$rho),col="red")
lines(seq(-pi,pi,0.01),dwrpcauchy(seq(-pi,pi,0.01),mle.wrappedcauchy(TA)$mu,mle.wrappedcauchy(TA)$rho),col="green")



sum((ddd-dwrpnorm(seq(-pi,pi,0.01),mle.wrappednormal(TA)$mu,mle.wrappednormal(TA)$rho))^2)
sum((ddd-dwrpnorm(seq(-pi,pi,0.01),mle.wrappedcauchy(TA)$mu,mle.wrappedcauchy(TA)$rho))^2)

sum((ddd - density)^2)

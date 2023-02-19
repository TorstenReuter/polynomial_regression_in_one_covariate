#######################################################################
#### Finding and plotting the Term c(alpha) in Theorem 5.2.        ####
#### Quadratic Regression with standard normal covariate           ####
#### working under the (false) assumption that b = 0               ####
#######################################################################

library(ggplot2)

# m2 in the information matrix of the design without interior interval
m2 <- function(ga) {
  y <- ga + sqrt(2/pi)*qnorm(1-ga/2,0,1)*exp(-(qnorm(1-ga/2,0,1)^2)/2) 
  return(y)
}

# m4 in the information matrix of the design without interior interval
m4 <- function(ga) {
  y <- sqrt(2/pi)*(qnorm(1-ga/2,0,1)^3)*exp(-(qnorm(1-ga/2,0,1)^2)/2) + 3*m2(ga)
  return(y)
}

# Vector of values of alpha
ga <- seq(0.02,0.999,0.001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qnorm(1 - ga/2,0,1)
  c <- ga*t^2 * ((2*m2(ga)- (ga*t^2))/(ga*m4(ga) - m2(ga)^2) - 1/m2(ga))
}

c <- cfunc(ga)
dfc <- data.frame("alpha" = ga, "c" = c)


plot1 <- ggplot() +
  geom_line(data = dfc, aes(x = alpha, y = c), size = 1, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("c(",alpha,")")))

plot1
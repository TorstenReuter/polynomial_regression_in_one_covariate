###################################################################
####### Linear Regression with half-t_3-distributed         #######
####### covariate. Checking if b>0 for all alpha.           #######
###################################################################

library(ggplot2)
library(extraDistr)

# m1 in Information Matrix of design without left interval
m1 <- function(ga) {
  t <- qht(1 - ga, nu = 3)
  y <- 6*sqrt(3)/((3 + t^2)*pi)
  return(y)
}


# m2 in Information Matrix of design without left interval
m2 <- function(ga) {
  t <- qht(1 -ga, nu = 3)
  y <- 3 + (6*sqrt(3)*t/((3 + t^2)*pi)) - 6*atan(t/sqrt(3))/pi
  return(y)
}

# Vector of values for alpha
ga <- seq(0.001,0.999,0.001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qht(1 -ga, nu = 3)
  c <- ga*t*((2*m1(ga) - ga*t) / (m1(ga)^2 - ga*m2(ga)))
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
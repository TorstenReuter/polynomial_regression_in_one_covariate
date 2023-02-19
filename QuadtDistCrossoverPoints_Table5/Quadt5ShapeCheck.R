######################################################################
#### Table 5                                                      ####
#### Quadratic Regression, t_5 distributed covariate              ####
#### Finding the crossover point where the D-optimal subsampling  ####
#### design switches from 3 disjoint intervals in its support     ####
#### to only two intervals.                                       ####
#### This is analogous to Theorem 5.6                             ####
######################################################################

library(ggplot2)

# m2 in Information Matrix of design without interior interval
m2 <- function(x) {
  y <- (10/(6*pi)) * (-((2*sqrt(5)*x*(x^2 - 5))/(x^2 + 5)^2) - 2*atan(x/sqrt(5)) + pi)
  return(y)
}


# m4 in Information Matrix of design without interior interval
m4 <- function(x) {
  y <- (50/(6*pi)) * (((10*sqrt(5)*x*(x^2 + 3))/(x^2 + 5)^2) - 6*atan(x/sqrt(5)) + 3*pi)  
  return(y)
}

# Vector of values for alpha
ga <- seq(0.00001,0.99999,0.00001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qt(1 - ga/2, df = 5)
  c <- ga*t^2 * ((2*m2(t)- (ga*t^2))/(ga*m4(t) - m2(t)^2) - 1/m2(t))
}

c <- cfunc(ga)
dfc <- data.frame("alpha" = ga, "c" = c)

plot1 <- ggplot() +
  geom_line(data = dfc, aes(x = alpha, y = c), size = 1, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0.082065, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("c(",alpha,")")))

plot1


######################################################################
#### Table 5                                                      ####
#### Quadratic Regression, t_6 distributed covariate              ####
#### Finding the crossover point where the D-optimal subsampling  ####
#### design switches from 3 disjoint intervals in its support     ####
#### to only two intervals.                                       ####
#### This is analogous to Theorem 5.6                             ####
######################################################################

library(ggplot2)

# m2 in Information Matrix of design without interior interval
m2 <- function(x) {
  y <- (3/2) * (1 - (((15 + x^2)*x^3)/(6 + x^2)^(2.5)))
  return(y)
}


# m4 in Information Matrix of design without interior interval
m4 <- function(x) {
  y <- (27/2) * (1 - ((x^5)/(6 + x^2)^(2.5)))
  return(y)
}

# Vector of values for alpha
ga <- seq(0.01,0.99999,0.00001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qt(1 - ga/2, df = 6)
  c <- ga*t^2 * ((2*m2(t)- (ga*t^2))/(ga*m4(t) - m2(t)^2) - 1/m2(t))
}

c <- cfunc(ga)
dfc <- data.frame("alpha" = ga, "c" = c)

plot1 <- ggplot() +
  geom_line(data = dfc, aes(x = alpha, y = c), size = 1, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0.346695, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("c_6")))

plot1


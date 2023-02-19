######################################################################
#### Table 5                                                      ####
#### Quadratic Regression, t_7 distributed covariate              ####
#### Finding the crossover point where the D-optimal subsampling  ####
#### design switches from 3 disjoint intervals in its support     ####
#### to only two intervals.                                       ####
#### This is analogous to Theorem 5.6                             ####
######################################################################

library(ggplot2)

# m2 in Information Matrix of design without interior interval
m2 <- function(x) {
  y <- (7/15) * (3 - ((2*sqrt(7)*x*(x^2 + 21)*(3*x^2 - 7))/ 
                        (pi*(7 + x^2)^3))
                 - (6*atan(x/sqrt(7)))/pi)
  return(y)
}


# m4 in Information Matrix of design without interior interval
m4 <- function(x) {
  y <- (49/15) * (3 - ((2*sqrt(7)*x*(x^2 - 21)*(3*x^2 + 7))/ 
                          (pi*(7 + x^2)^3))
                        - (6*atan(x/sqrt(7)))/pi)  
  return(y)
}

# Vector of values for alpha
ga <- seq(0.01,0.99999,0.00001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qt(1 - ga/2, df = 7)
  c <- ga*t^2 * ((2*m2(t)- (ga*t^2))/(ga*m4(t) - m2(t)^2) - 1/m2(t))
}

c <- cfunc(ga)
dfc <- data.frame("alpha" = ga, "c" = c)

plot1 <- ggplot() +
  geom_line(data = dfc, aes(x = alpha, y = c), size = 1, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0.503735, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("c_7")))

plot1


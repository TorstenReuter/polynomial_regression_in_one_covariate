######################################################################
#### Table 5                                                      ####
#### Quadratic Regression, t_9 distributed covariate              ####
#### Finding the crossover point where the D-optimal subsampling  ####
#### design switches from 3 disjoint intervals in its support     ####
#### to only two intervals.                                       ####
#### This is analogous to Theorem 5.6                             ####
######################################################################

library(ggplot2)

# m2 in Information Matrix of design without interior interval
m2 <- function(x) {
  y <- (9/(35*pi)) * ((-6*x*(-3645 + 1971*x^2 + 165*x^4 + 5*x^6))/ 
                      ((9 + x^2)^4)
                      - 10*atan(x/3) + 5*pi)
  return(y)
}


# m4 in Information Matrix of design without interior interval
m4 <- function(x) {
  y <- (243/(35*pi)) * ((-6*(-3 + x)*x*(3 + x)*(x^4 + 42*x^2 + 81))/ 
                        ((9 + x^2)^4)
                      - 2*atan(x/3) + pi)  
  return(y)
}

# Vector of values for alpha
ga <- seq(0.01,0.99999,0.00001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qt(1 - ga/2, df = 9)
  c <- ga*t^2 * ((2*m2(t)- (ga*t^2))/(ga*m4(t) - m2(t)^2) - 1/m2(t))
}

c <- cfunc(ga)
dfc <- data.frame("alpha" = ga, "c" = c)

plot1 <- ggplot() +
  geom_line(data = dfc, aes(x = alpha, y = c), size = 1, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0.666995, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("c_9")))

plot1


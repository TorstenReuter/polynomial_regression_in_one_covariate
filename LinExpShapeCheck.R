###################################################################
####### Linear Regression with Exp(1) distributed covariate #######
####### Checking if b>0 for all alpha                       #######
###################################################################

library(ggplot2)

# m1 in Information Matrix of design without left interval
m1 <- function(ga) {
  y <- ((log(1/ga) + 1)*ga)
  return(y)
}


# m2 in Information Matrix of design without left interval
m2 <- function(ga) {
  y <- ((2 + 2*log(1/ga) + log(1/ga)^2)*ga)
  return(y)
}

# Vector of values for alpha
ga <- seq(0.001,0.999,0.001)

# The term c(alpha)
cfunc <- function(ga) {
  c <- ga*log(1/ga)*((2*m1(ga) - ga*log(1/ga)) / (m1(ga)^2 - ga*m2(ga)))
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
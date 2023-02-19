#####################################################################
#### Plotting Figure 5.                                          ####
#### Script for finding and plotting the percentag of mass on    ####
#### [-b,b] and [a,1} of the D-optimal subsampling design for    #### 
#### quadratic regression given a uniform covariate on [-1,1].   ####
#####################################################################

library(ggplot2)

# Vector of values for alpha
ga <- seq(0.0001,0.9999,0.0001)

# Calculating a and b
a <- (1 - ((1/2)*(1 - ga + sqrt((45 - 15*ga + 15*ga^2 - 45*ga^3 + 20*ga^4 -
                             4*ga*sqrt(5)*sqrt(45 - 90*ga + 90*ga^2 -
                                                 75*ga^3 + 57*ga^4 - 27*ga^5 + 5*ga^6))/(45*(1-ga)))))) / (2*ga)
b <- (((1/2)*(1 - ga + sqrt((45 - 15*ga + 15*ga^2 - 45*ga^3 + 20*ga^4 -
                               4*ga*sqrt(5)*sqrt(45 - 90*ga + 90*ga^2 -
                                                   75*ga^3 + 57*ga^4 - 27*ga^5 + 5*ga^6))/(45*(1-ga))))) - (1 - ga))/(ga)






dfa <- data.frame("alpha" = ga, "a" = a)
dfb <- data.frame("alpha" = ga, "b" = b)

plot1 <- ggplot() +
  geom_line(data = dfa, aes(x = alpha, y = a), size = 1, linetype = "solid") +
  geom_hline(yintercept=1/3, linetype="dotted", size = 1) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("(1-a)/2",alpha)))

plot1

plot2 <- ggplot() +
  geom_line(data = dfb, aes(x = alpha, y = b), size = 1, linetype = "solid") +
  geom_hline(yintercept=1/3, linetype="dotted", size = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("b/",alpha)))

plot2

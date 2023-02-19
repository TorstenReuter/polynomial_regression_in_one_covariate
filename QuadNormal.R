###################################################################
#### Script for finding and plotting the D-optimal subsampling ####
#### design in quadratic regression, given a standard normal   ####
#### covariate                                                 ####
###################################################################

library(nleqslv)
library(ggplot2)
library(gridExtra)


# set alpha
ga <- 0.5

# set the range in which the sensitivity function should be displayed
# this still needs to be set individually
DDrange <- 1.25  # ga 0.5
#DDrange <- 1.55  # ga 0.3
#DDrange <- 2.15  # ga 0.1
#DDrange <- 2.85  # ga 0.01 

######### Calculating a and b and finding the values in Table 2 ############

# Initial guess of boundary points
# this still needs to be set individually 
# but xstart <- c(3,0.1) works for a wide range of alpha
xstart <- c(3,0.1) # ga 0.01 or 0.1 or 0.3 or 0.5



# Regression function
f <- function(x) {
  y <- numeric(3)
  y[1] <- 1
  y[2] <- x
  y[3] <- x^2
  y
}

# m2 in Information Matrix x[1] = a, x[2] = b
Term1 <- function(x) {
  y <- sqrt(2/pi)*(x[1]*exp(-x[1]^2/2) - x[2]*exp(-x[2]^2/2)) + 
    1 + 2*(pnorm(x[2]) - pnorm(x[1])) 
  y
}


# Term 2 in Information Matrix where m4 = Term2 + 3*Term1; x[1] = a, x[2] = b
Term2 <- function(x) {
  y <- sqrt(2/pi)*((x[1]^3)*exp(-x[1]^2/2) - (x[2]^3)*exp(-x[2]^2/2))  
  y
}


# System of non-linear Equation to solve
dslnex <- function(x) {
  M <- cbind(c(ga,0,Term1(x)),c(0,Term1(x),0),c(Term1(x),0,Term2(x)+3*Term1(x)))
  y <- numeric(2)
  y[1] <- pnorm(x[2]) - pnorm(x[1]) + 0.5 - ga/2
  y[2] <- t(f(x[1]))%*%solve(M)%*%f(x[1]) - t(f(x[2]))%*%solve(M)%*%f(x[2])
  y
}
# Solver
para <- nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))$x

# Finding the values in Table 2
a <- para[1]
b <- para[2]

a
b

1 - pnorm(a)

2*pnorm(b) - 1

(2*pnorm(b) - 1)/ga

########################## Plotting #################################


n <- 100000
points <- seq(-4, 4, length.out = n)


CDFphi <- dnorm(points, mean = 0, sd = 1)
CDFg <- seq(-4, 4, length.out = n)
for(i in 1:n) {
  if ((CDFg[i] > -a && CDFg[i] < -b) || (CDFg[i] < a && CDFg[i] > b) ) {
    CDFg[i] <- 0
  } else {
    CDFg[i] <- dnorm(CDFg[i], mean = 0, sd = 1)
  }
}


dfphi <- data.frame("x" = points, "CDF" = CDFphi)
dfg <- data.frame("x" = points, "CDF" = CDFg)

plot1 <- ggplot() +
  geom_line(data = dfphi, aes(x = x, y = CDF, color = "black"), size = 1.4, linetype = "dashed") +
  geom_line(data = dfg, aes(x = x, y = CDF, color = "black"), size = 1, linetype = "solid") +
  scale_color_identity(name = element_blank(),
                       breaks = c("blue", "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(-3.5, 3.5)) +
  xlab("x") + 
  ylab("Density")



plot1



# Now plot the directional derivative (DD)

# define basic functions that appear in DD
f <- function(x) {
  y <- numeric(3)
  y[1] <- 1
  y[2] <- x
  y[3] <- x^2
  y
}

Term1 <- function(x) {
  y <- sqrt(2/pi)*(x[1]*exp(-x[1]^2/2) - x[2]*exp(-x[2]^2/2)) + 
    1 + 2*(pnorm(x[2]) - pnorm(x[1])) 
  y
}


Term2 <- function(x) {
  y <- sqrt(2/pi)*((x[1]^3)*exp(-x[1]^2/2) - (x[2]^3)*exp(-x[2]^2/2))  
  y
}

M <- cbind(c(ga,0,Term1(para)),c(0,Term1(para),0),c(Term1(para),0,Term2(para)+3*Term1(para)))


# DD
FPhi <- function(x){
  y <- ga*(t(f(x)) %*% solve(M) %*% f(x))
  y
}

n2 <- 10000
points2 <- seq(-DDrange, DDrange, length.out = n2)
DD <- seq(-DDrange, DDrange, length.out = n2)

for(i in 1:n2) {
  DD[i] <- FPhi(DD[i])
}

dfFphi <- data.frame("x" = points2, "FPhi" = DD)

# set the lower end of the vertical bars
lowerylimit <- FPhi(-DDrange)

plot2 <- ggplot() +
  geom_line(data = dfFphi, aes(x = x, y = FPhi, color = "black"), size = 1) +
  geom_hline(yintercept=FPhi(a), linetype="dotted", color = "black") +
  geom_segment(aes(x = -a, y = min(dfFphi$FPhi), xend = -a, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = a, y = min(dfFphi$FPhi), xend = a, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = -b, y = min(dfFphi$FPhi), xend = -b, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = b, y = min(dfFphi$FPhi), xend = b, yend = max(dfFphi$FPhi)), linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(-3.5, 3.5)) +
  scale_color_identity(name = element_blank(),
                       breaks = c("red")) +
  xlab("x") + 
  ylab("Sensitivity function")


plot2


plot1 <- ggplot_gtable(ggplot_build(plot1))
plot2 <- ggplot_gtable(ggplot_build(plot2))

plot2$widths <- plot1$widths

grid.arrange(plot1, plot2, ncol=1)


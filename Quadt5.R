######################################################################
#### Table 4                                                      ####
#### Quadratic Regression, t_5 distributed covariate              ####
#### Script for finding and plotting the D-optimal subsampling    ####
#### design and finding the values in Table 4                     ####
######################################################################

library(nleqslv)
library(ggplot2)
library(gridExtra)


# set alpha in (0, 0.82). 
# As the D-optimal subsampling design has only 2 disjoint intervals 
# in its support for larger alpha.
ga <- 0.05

# Range in which sensitivity function will be displayed
DDrange <- 3  


################## Calculating a and b #############################

# Initial guess of a and b 
if (ga < 0.003){
  xstart <- c(10,0.05)
} else if (ga < 0.008) {
  xstart <- c(8,0.1)
} else if (ga < 0.02) {
  xstart <- c(4,0.1)
} else {
  xstart <- c(3,0.1)
}


# Regression function
f <- function(x) {
  y <- numeric(3)
  y[1] <- 1
  y[2] <- x
  y[3] <- x^2
  y
}

# m2 in Information Matrix
m2 <- function(x) {
  y <- (10/(3*pi)) * (((sqrt(5)*x[2]*(x[2]^2 - 5))/(x[2]^2 + 5)^2) + atan(x[2]/sqrt(5))) + 
    (10/(6*pi)) * (-((2*sqrt(5)*x[1]*(x[1]^2 - 5))/(x[1]^2 + 5)^2) - 2*atan(x[1]/sqrt(5)) + pi)
  y
}


# m4 in Information Matrix
m4 <- function(x) {
  y <- (50/(3*pi)) * (((-5*sqrt(5)*x[2]*(x[2]^2 + 3))/(x[2]^2 + 5)^2) + 3*atan(x[2]/sqrt(5))) + 
    (50/(6*pi)) * (((10*sqrt(5)*x[1]*(x[1]^2 + 3))/(x[1]^2 + 5)^2) - 6*atan(x[1]/sqrt(5)) + 3*pi)  
  y
}



# System of non-linear equations to solve
dslnex <- function(x) {
  M <- cbind(c(ga,0,m2(x)),
             c(0,m2(x),0),
             c(m2(x),0,m4(x)))
  Minv <- solve(M)
  y <- numeric(2)
  y[1] <- pt(x[2], df = 5) - pt(x[1], df = 5) + 0.5 - ga/2
  y[2] <- t(f(x[1]))%*%Minv%*%f(x[1]) - t(f(x[2]))%*%Minv%*%f(x[2])
  y
}

# Solveing for a and b and finding the values in Table 4
para <- nleqslv(xstart, dslnex,method = "Newton", control=list(trace=1,btol=.01))$x
a <- para[1]
b <- para[2]

para
1 - pt(a, df = 5)
2*(pt(b, df = 5) - (1/2))
check <- 2*(1 - pt(a, df = 5) + pt(b, df = 5) - (1/2))
check

percentageonmiddle <- (pt(b, df = 5) - (1/2)) / (pt(b, df = 5) - (1/2) + 1 - pt(a, df = 5)) 
percentageonmiddle
########################## Plotting #################################


n <- 100000
points <- seq(-5, 5, length.out = n)


CDFphi <- dt(points, df = 5)
CDFg <- seq(-5, 5, length.out = n)
for(i in 1:n) {
  if ((CDFg[i] > -a && CDFg[i] < -b) || (CDFg[i] < a && CDFg[i] > b) ) {
    CDFg[i] <- 0
  } else {
    CDFg[i] <- dt(CDFg[i], df = 5)
  }
}


dfphi <- data.frame("x" = points, "CDF" = CDFphi)
dfg <- data.frame("x" = points, "CDF" = CDFg)

plot1 <- ggplot() +
  geom_line(data = dfphi, aes(x = x, y = CDF), size = 1.4, linetype = "dashed") +
  geom_line(data = dfg, aes(x = x, y = CDF), size = 1, linetype = "solid") +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
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

M <- cbind(c(ga,0,m2(para)),c(0,m2(para),0),c(m2(para),0,m4(para)))


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
  geom_line(data = dfFphi, aes(x = x, y = FPhi), size = 1) +
  geom_hline(yintercept=FPhi(a), linetype="dotted", color = "black") +
  geom_segment(aes(x = -a, y = min(dfFphi$FPhi), xend = -a, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = a, y = min(dfFphi$FPhi), xend = a, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = -b, y = min(dfFphi$FPhi), xend = -b, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = b, y = min(dfFphi$FPhi), xend = b, yend = max(dfFphi$FPhi)), linetype="dotted") +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("x") + 
  ylab("Sensitivity function")


plot2


plot1 <- ggplot_gtable(ggplot_build(plot1))
plot2 <- ggplot_gtable(ggplot_build(plot2))

plot2$widths <- plot1$widths

grid.arrange(plot1, plot2, ncol=1)


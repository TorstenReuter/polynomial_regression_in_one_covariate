####################################################################
#### Finding and Plotting the D-optimal subsampling design      ####
#### for linear regression given a Exp(lamda) distributed       ####
#### covariate                                                  ####
####################################################################


library(nleqslv)
library(ggplot2)
library(gridExtra)


# set 0 < alpha < 1 and lambda > 0 
ga <- 0.5
gl <- 1

# set the range in which the sensitivity function should be displayed
# this still needs to be set individually
DDrange <- 4      # ga = 0.5, 0.3
#DDrange <- 2.15   # ga = 0.1
#DDrange <- 20     # ga = 0.01 


############### Calculating a and b and values of Table 1 ##################

# Initial guess of a and b 
# this still needs to be set individually
xstart <- c(0.5, 2)    # ga 0.5-0.3
#xstart <- c(0.1, 3.3)  # ga 0.1
#xstart <- c(0.1, 5)    # ga 0.01



# Regression function 
f <- function(x) {
  y <- numeric(2)
  y[1] <- 1
  y[2] <- x
  y
}

# m_1 in Information Matrix x[1] = b; x[2] = a 
Term1 <- function(x) {
  y <-  (1/gl) + (x[2] + 1/gl)*exp(-gl*x[2]) - (x[1] + 1/gl)*exp(-gl*x[1]) 
  y
}


# m_2 in Information Matrix x[1] = b; x[2] = a
Term2 <- function(x) {
  y <-  (2/gl^2) + (2/gl^2 + x[2]^2 + 2*x[2]/gl)*exp(-gl*x[2]) + (-x[1]^2 - 2*x[1]/gl - 2/gl^2)*exp(-gl * x[1])
  y
}


# System of non-linear equations to solve
dslnex <- function(x) {
  M <- cbind(c(ga,Term1(x)),c(Term1(x),Term2(x)))
  y <- numeric(2)
  y[1] <- 1 - exp(-gl * x[1]) + exp(-gl * x[2]) - ga
  y[2] <- t(f(x[1]))%*%solve(M)%*%f(x[1]) - t(f(x[2]))%*%solve(M)%*%f(x[2])
  y
}
# Solver
para <- nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))$x

# Finding values for Table 1
b <- para[1]
a <- para[2]

a
b

check <- 1 + exp(-gl * a) - exp(-gl * b)
check

Integral_0tob <- 1 - exp(-gl * b)
Integral_0tob

Integral_atoInfty <- exp(-gl * a)
Integral_atoInfty

percentage_0tob <- Integral_0tob/ga
percentage_0tob

######################## Plotting for Figure 1 #############################
PDF <- function(x) {
  y <- gl * exp(-gl*x)
  y
}

n <- 10000
points <- seq(0, 5, length.out = n)


CDFphi <- PDF(points)
CDFg <- seq(0, 5, length.out = n)
for(i in 1:n) {
  if ((CDFg[i] < b) || (CDFg[i] >  a)) {
    CDFg[i] <- PDF(CDFg[i])
  } else {
    CDFg[i] <- 0
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
  scale_x_continuous(limits = c(0, 5)) +
  xlab("x") + 
  ylab("Density")



plot1


# Define basic functions that appear in sensitivity function
f <- function(x) {
  y <- numeric(2)
  y[1] <- 1
  y[2] <- x
  y
}


M <- cbind(c(ga,Term1(para)),c(Term1(para),Term2(para)))


# Sensitivity function
FPhi <- function(x){
  y <- ga*(t(f(x)) %*% solve(M) %*% f(x))
  y
}

n2 <- 100000
points2 <- seq(0, DDrange, length.out = n2)
DD <- seq(0, DDrange, length.out = n2)

for(i in 1:n2) {
  DD[i] <- FPhi(DD[i])
}

dfFphi <- data.frame("x" = points2, "FPhi" = DD)


plot2 <- ggplot() +
  geom_line(data = dfFphi, aes(x = x, y = FPhi, color = "black"), size = 1) +
  geom_hline(yintercept=FPhi(a), linetype="dotted", color = "black") +
  geom_segment(aes(x = -a, y = min(dfFphi$FPhi), xend = -a, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = a, y = min(dfFphi$FPhi), xend = a, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = -b, y = min(dfFphi$FPhi), xend = -b, yend = max(dfFphi$FPhi)), linetype="dotted")+
  geom_segment(aes(x = b, y = min(dfFphi$FPhi), xend = b, yend = max(dfFphi$FPhi)), linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(0, 5)) +
  scale_color_identity(name = element_blank(),
                       breaks = c("red")) +
  xlab("x") + 
  ylab("Sensitivity function")


plot2


plot1 <- ggplot_gtable(ggplot_build(plot1))
plot2 <- ggplot_gtable(ggplot_build(plot2))

plot2$widths <- plot1$widths

grid.arrange(plot1, plot2, ncol=1)



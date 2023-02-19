#####################################################################
#### Script for finding and plotting the D-optimal subsampling   ####
#### design for quadratic regression given a uniform covariate   ####
#### on [-1,1].                                                  ####
#####################################################################

library(nleqslv)
library(ggplot2)
library(gridExtra)


# set alpha
ga <- 0.5

################## Calculating a and b #############################



# set the range in which the sensitivity funcition should be displayed
# this still needs to be set individually
DDrange <- 0.82  # ga 0.5
#DDrange <- 0.9  # ga ga 0.3
#DDrange <- 1.05  # ga 0.01 or 0.1


# calculating a and b
a <- (1/2)*(1 - ga + sqrt((45 - 15*ga + 15*ga^2 - 45*ga^3 + 20*ga^4 -
                           4*ga*sqrt(5)*sqrt(45 - 90*ga + 90*ga^2 -
                           75*ga^3 + 57*ga^4 - 27*ga^5 + 5*ga^6))/(45*(1-ga))))
b <- a - (1 - ga)
para <- cbind(a,b)


(1/2) - (1/2) * a
(1/2)*b

########################## Plotting #################################

n <- 100000
points <- seq(-1.2, 1.2, length.out = n)

CDFphi <- seq(-1.2, 1.2, length.out = n)
for(i in 1:n) {
  if ((CDFphi[i] < -1) || (CDFphi[i]  > 1) ) {
    CDFphi[i] <- 0
  } else {
    CDFphi[i] <- 1/2
  }
}


CDFg <- seq(-1.2, 1.2, length.out = n)
for(i in 1:n) {
  if ((CDFg[i] > -para[1] && CDFg[i] < -para[2]) || 
        (CDFg[i] < para[1] && CDFg[i] > para[2]) || 
         CDFg[i] < -1 || 
         CDFg[i] > 1) {
    CDFg[i] <- 0
  } else {
    CDFg[i] <- 1/2
  }
}


dfphi <- data.frame("x" = points, "CDF" = CDFphi)
dfg <- data.frame("x" = points, "CDF" = CDFg)

plot1 <- ggplot() +
  geom_line(data = dfphi, aes(x = x, y = CDF), size = 1.4, linetype = "dashed") +
  geom_line(data = dfg, aes(x = x, y = CDF), size = 1, linetype = "solid") +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("x") + 
  ylab("Density")



plot1


# define basic functions that appear in sensitivity function
f <- function(x) {
  y <- numeric(3)
  y[1] <- 1
  y[2] <- x
  y[3] <- x^2
  y
}

# m2 in Information Matrix
Term1 <- function(x) {
  y <- 1/3 - (1/3)*x[1]^3 + (1/3)*x[2]^3
  y
}


# m4 in Information Matrix
Term2 <- function(x) {
  y <- 1/5 - (1/5)*x[1]^5 + (1/5)*x[2]^5  
  y
}

M <- cbind(c(ga,0,Term1(para)),c(0,Term1(para),0),c(Term1(para),0,Term2(para)))


# sensitivity function
FPhi <- function(x){
  y <- ga * (t(f(x)) %*% solve(M) %*% f(x))
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
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("x") + 
  ylab("Sensitivity function")



plot2


plot1 <- ggplot_gtable(ggplot_build(plot1))
plot2 <- ggplot_gtable(ggplot_build(plot2))

plot2$widths <- plot1$widths

grid.arrange(plot1, plot2, ncol=1)


####################################################################
#### Table 6 / Figure 6                                         ####
#### Quadratic Regression, standard normal covariate            ####
#### Calculating the Efficiencies of uniform random subsampling ####
#### and an IBOSS-like subsampling that allocates a third of    ####
#### its mass on either of the three intervals                  ####
####################################################################

library(ggplot2)
library(nleqslv)


# Term 1 in Information Matrix
Term1 <- function(x) {
  y <- sqrt(2/pi)*(x[1]*exp(-x[1]^2/2) - x[2]*exp(-x[2]^2/2)) + 
    1 + 2*(pnorm(x[2]) - pnorm(x[1])) 
  y
}


# Term 2 in Information Matrix
Term2 <- function(x) {
  y <- sqrt(2/pi)*((x[1]^3)*exp(-x[1]^2/2) - (x[2]^3)*exp(-x[2]^2/2))  
  y
}

########## efficiency ######################
efficiency <- function(ga, x){
  Mstar <- cbind(c(ga,0,Term1(x)),c(0,Term1(x),0),c(Term1(x),0,Term2(x)+3*Term1(x)))
  dstar <- det(Mstar)
  d <- 2*(ga^3)
  eff <- (d / dstar)^(1/dim(Mstar)[1])
  return(eff)
}

####### efficiency of IBOSS like design #############
IBOSSeff <- function(ga, x){
  Mstar <- cbind(c(ga,0,Term1(x)),c(0,Term1(x),0),c(Term1(x),0,Term2(x)+3*Term1(x)))
  dstar <- det(Mstar)
  a <- qnorm(1 - ga/3, mean = 0, sd = 1)
  b <- qnorm(0.5 + ga/6) 
  anb <- c(a, b)
  MIBOSS <- cbind(c(ga,0,Term1(anb)),c(0,Term1(anb),0),c(Term1(anb),0,Term2(anb)+3*Term1(anb)))
  dIBOSS <- det(MIBOSS)
  eff <- (dIBOSS / dstar)^(1/dim(Mstar)[1])
  return(eff)
}



# Regression function
f <- function(x) {
  y <- numeric(3)
  y[1] <- 1
  y[2] <- x
  y[3] <- x^2
  y
}


############ calculating efficiency for all ga ###############
ga <- seq(0.0001,0.9999,0.0001)


eff <- vector()
effIBOSS <- vector()
for (alpha in ga) {
  if (alpha < 0.003){
    xstart <- c(4,0.01)
  } else if (alpha < 0.1){
    xstart <- c(3,0.1)
  } else if (alpha >= 0.1 & alpha < 0.5) {
    xstart <- c(2.2,0.1)
  } else if (alpha >= 0.5 & alpha < 0.7) {
    xstart <- c(1.5,0.1)
  } else {
    xstart <- c(1,0.1)
  }
  dslnex <- function(x) {
    M <- cbind(c(alpha,0,Term1(x)),c(0,Term1(x),0),c(Term1(x),0,Term2(x)+3*Term1(x)))
    y <- numeric(2)
    y[1] <- pnorm(x[2]) - pnorm(x[1]) + 0.5 - alpha/2
    y[2] <- t(f(x[1]))%*%solve(M)%*%f(x[1]) - t(f(x[2]))%*%solve(M)%*%f(x[2])
    y
  }
  para <- nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))$x
  eff <- c(eff, efficiency(alpha, para))
  effIBOSS <- c(effIBOSS, IBOSSeff(alpha, para))
}

dfeff <- data.frame("alpha" = ga, "Efficiency" = eff)
dfIBOSSeff <- data.frame("alpha" = ga, "Efficiency" = effIBOSS)

plot1 <- ggplot() +
  geom_line(data = dfeff, aes(x = alpha, y = eff), size = 1, linetype = "solid") +
  geom_line(data = dfIBOSSeff, aes(x = alpha, y = effIBOSS), size = 1, linetype = "dashed") +
  geom_hline(yintercept=1, linetype="dotted", size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("Efficiency"))) + 
  ylim(0,1)

plot1


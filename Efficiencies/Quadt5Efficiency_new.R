####################################################################
#### Table 6 / Figure 6                                         ####
#### Quadratic Regression, t_5 distributed covariate            ####
#### Calculating the Efficiencies of uniform random subsampling ####
#### and an IBOSS-like subsampling that allocates a third of    ####
#### its mass on either of the three intervals                  ####
####################################################################

library(ggplot2)
library(nleqslv)

TwoIntervalm2 <- function(ga) {
  y <- (10/(6*pi)) * 
    (-((2*sqrt(5)*qt(1 - ga/2, df = 5)*(qt(1 - ga/2, df = 5)^2 - 5))/
         (qt(1 - ga/2, df = 5)^2 + 5)^2) - 
       2*atan(qt(1 - ga/2, df = 5)/sqrt(5)) + pi)
  return(y)
}


# m4 in Information Matrix of design without interior interval
TwoIntervalm4 <- function(ga) {
  y <- (50/(6*pi)) * 
    (((10*sqrt(5)*qt(1 - ga/2, df = 5)*(qt(1 - ga/2, df = 5)^2 + 3))/
        (qt(1 - ga/2, df = 5)^2 + 5)^2) - 
       6*atan(qt(1 - ga/2, df = 5)/sqrt(5)) + 3*pi)  
  return(y)
}



# m2 in Information Matrix with interior interval
ThreeIntervalm2 <- function(x) {
  y <- (10/(3*pi)) * 
    (((sqrt(5)*x[2]*(x[2]^2 - 5))/(x[2]^2 + 5)^2) + atan(x[2]/sqrt(5))) + 
    (10/(6*pi)) * 
    (-((2*sqrt(5)*x[1]*(x[1]^2 - 5))/(x[1]^2 + 5)^2) - 2*atan(x[1]/sqrt(5)) + pi)
  return(y)
}


# m4 in Information Matrix with interior interval
ThreeIntervalm4 <- function(x) {
  y <- (50/(3*pi)) * 
    (((-5*sqrt(5)*x[2]*(x[2]^2 + 3))/(x[2]^2 + 5)^2) + 3*atan(x[2]/sqrt(5))) + 
    (50/(6*pi)) * 
    (((10*sqrt(5)*x[1]*(x[1]^2 + 3))/(x[1]^2 + 5)^2) - 6*atan(x[1]/sqrt(5)) + 3*pi)  
  return(y)
}

##### Finding det(M(\xi^*)) ################
dstarfunc <- function(ga){
  if (ga < 0.082065){
    if (ga < 0.003){
      xstart <- c(10,0.05)
    } else if (ga < 0.008) {
      xstart <- c(8,0.1)
    } else if (ga < 0.02) {
      xstart <- c(4,0.1)
    } else {
      xstart <- c(3,0.1)
    }
    dslnex <- function(x) {
      M <- cbind(c(ga,0,ThreeIntervalm2(x)),
                 c(0,ThreeIntervalm2(x),0),
                 c(ThreeIntervalm2(x),0,ThreeIntervalm4(x)))
      y <- numeric(2)
      y[1] <- pt(x[2], df = 5) - pt(x[1], df = 5) + 0.5 - alpha/2
      y[2] <- t(f(x[1]))%*%solve(M)%*%f(x[1]) - t(f(x[2]))%*%solve(M)%*%f(x[2])
      y
    }
    para <- nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))$x
    Mstar <- cbind(c(ga,0,ThreeIntervalm2(para)),
                      c(0,ThreeIntervalm2(para),0),
                      c(ThreeIntervalm2(para),0,ThreeIntervalm4(para)))
    dstar <- det(Mstar)
  } else {
    Mstar <- cbind(c(ga,0,TwoIntervalm2(ga)),
                   c(0,TwoIntervalm2(ga),0),
                   c(TwoIntervalm2(ga),0,TwoIntervalm4(ga)))
    dstar <- det(Mstar)
  }
  return(dstar)
}


#### efficiency of uniform random subsampling ####
efficiency <- function(ga){
  dstar <- dstarfunc(ga)
  d <- (125*ga^3 / 3) - (125*ga^3/27)
  eff <- (d / dstar)^(1/3)
  return(eff)
}


####### efficiency of IBOSS like design #############
IBOSSeff <- function(ga){
  dstar <- dstarfunc(ga)
  a <- qt(1 - ga/3, df = 5)
  b <- qt(0.5 + ga/6, df = 5) 
  anb <- c(a, b)
  MIBOSS <- cbind(c(ga,0,ThreeIntervalm2(anb)),
                  c(0,ThreeIntervalm2(anb),0),
                  c(ThreeIntervalm2(anb),0,ThreeIntervalm4(anb)))
  dIBOSS <- det(MIBOSS)
  eff <- (dIBOSS / dstar)^(1/3)
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
  eff <- c(eff, efficiency(alpha))
  effIBOSS <- c(effIBOSS, IBOSSeff(alpha))
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


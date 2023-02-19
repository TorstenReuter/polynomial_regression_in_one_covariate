####################################################################
#### Table 6 / Figure 6                                         ####
#### Linear Regression, Exp(1) distributed covariate            ####
#### Calculating the Efficiencies of uniform random subsampling ####
#### and an IBOSS-like subsampling that allocates half of its   ####
#### mass on either interval                                    ####
####################################################################

library(nleqslv)
library(ggplot2)

# setting lambda. If changed, the IBOSS-like design is wrong. 
gl <- 1

# m_1 in Information Matrix x[1] = b; x[2] = a 
Term1 <- function(x) {
  y <-  (1/gl) + (x[2] + 1/gl)*exp(-gl*x[2]) - (x[1] + 1/gl)*exp(-gl*x[1]) 
  return(y)
}

# m_2 in Information Matrix x[1] = b; x[2] = a
Term2 <- function(x) {
  y <-  (2/gl^2) + (2/(gl^2) + x[2]^2 + ((2*x[2])/gl))*exp(-gl*x[2]) + (-x[1]^2 - 2*x[1]/gl - 2/gl^2)*exp(-gl * x[1])
  return(y)
}

########## efficiency of \xi_0 ######################
efficiency <- function(ga, gl, x){
  Mstar <- cbind(c(ga,Term1(x)),c(Term1(x),Term2(x)))
  dstar <- det(Mstar)
  d <- (ga^2) / (gl^2)
  eff <- (d / dstar)^(1/dim(Mstar)[1])
  return(eff)
}

####### efficiency of IBOSS like design only for gl = 1 #############
IBOSSeff <- function(ga, gl, x){
  Mstar <- cbind(c(ga,Term1(x)),c(Term1(x),Term2(x)))
  dstar <- det(Mstar)
  left <- -log(1 - ga/2)
  right <- -log(ga/2) 
  leftright <- c(left, right)
  MIBOSS <- cbind(c(ga,Term1(leftright)),c(Term1(leftright),Term2(leftright)))
  dIBOSS <- det(MIBOSS)
  eff <- (dIBOSS / dstar)^(1/dim(Mstar)[1])
  return(eff)
}

# Regression function
f <- function(x) {
  y <- numeric(2)
  y[1] <- 1
  y[2] <- x
  return(y)
}

############ calculating efficiency for all ga ###############
ga <- seq(0.0001,0.9999,0.0001)

eff <- vector()
effIBOSS <- vector()
for (alpha in ga) {
  if (alpha < 0.05){
    xstart <- c(0,10)
  } else if (alpha < 0.1) {
    xstart <- c(0.1,5)
  } else if (alpha < 0.5) {
    xstart <- c(0.3,4)
  } else if (alpha < 0.7) {
    xstart <- c(1,3)
  } else {
    xstart <- c(1,2)
  }
  dslnex <- function(x) {
    M <- cbind(c(alpha,Term1(x)),c(Term1(x),Term2(x)))
    y <- numeric(2)
    y[1] <- 1 - exp(-gl * x[1]) + exp(-gl * x[2]) - alpha
    y[2] <- t(f(x[1]))%*%solve(M)%*%f(x[1]) - t(f(x[2]))%*%solve(M)%*%f(x[2])
    y
  }
  para <- nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))$x
  eff <- c(eff, efficiency(alpha, gl, para))
  effIBOSS <- c(effIBOSS, IBOSSeff(alpha, gl, para))
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



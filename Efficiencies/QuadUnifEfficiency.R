####################################################################
#### Table 6 / Figure 6                                         ####
#### Quadratic Regression, uniform on [-1,1] covariate          ####
#### Calculating the Efficiencies of uniform random subsampling ####
#### and an IBOSS-like subsampling that allocates a third of    ####
#### its mass on either of the three intervals                  ####
####################################################################

library(ggplot2)
library(gridExtra)


################## Functions for a and b #############################

aFunc <- function(ga){
  a <- (1/2)*(1 - ga + sqrt((45 - 15*ga + 15*ga^2 - 45*ga^3 + 20*ga^4 -
                               4*ga*sqrt(5)*sqrt(45 - 90*ga + 90*ga^2 -
                                                   75*ga^3 + 57*ga^4 - 27*ga^5 + 5*ga^6))/(45*(1-ga)))) 
  a
}

bFunc <- function(ga){
  b <- ((1/2)*(1 - ga + sqrt((45 - 15*ga + 15*ga^2 - 45*ga^3 + 20*ga^4 -
                               4*ga*sqrt(5)*sqrt(45 - 90*ga + 90*ga^2 -
                               75*ga^3 + 57*ga^4 - 27*ga^5 + 5*ga^6))/(45*(1-ga))))) -
                              (1 - ga)
  b
}


# m2 in Information Matrix of optimal design
m2star <- function(ga) {
  y <- 1/3 - (1/3)*aFunc(ga)^3 + (1/3)*bFunc(ga)^3
  y
}


# m4 in Information Matrix of optimal design
m4star <- function(ga) {
  y <- 1/5 - (1/5)*aFunc(ga)^5 + (1/5)*bFunc(ga)^5  
  y
}

# m2 in Information Matrix of design
m2 <- function(x) {
  y <- 1/3 - (1/3)*x[1]^3 + (1/3)*x[2]^3
  y
}

# m4 in Information Matrix of design
m4 <- function(x) {
  y <- 1/5 - (1/5)*x[1]^5 + (1/5)*x[2]^5
  y
}

# determinant of the information matrix of the optimal design
detMstar <- function(ga){
  Mstar <- cbind(c(ga,0,m2star(ga)),c(0,m2star(ga),0),c(m2star(ga),0,m4star(ga)))
  return(det(Mstar))
} 
  


# determinat of  the information matrix of uniform random subsampling ##########

detMnull <- function(ga){
  Mnull <- cbind(c(ga,0,ga/3),c(0,ga/3,0),c(ga/3,0,ga/5))
  return(det(Mnull))
}

######## efficiency of uniform random subsampling ############

efficiency <- function(ga){
  eff <- (detMnull(ga) / detMstar(ga))^(1/3)
  return(eff)
} 

####### efficiency of IBOSS like design #############
IBOSSeff <- function(ga){
  dstar <- detMstar(ga)
  a <- 1 - 2*ga/3
  b <- ga/3 
  anb <- c(a, b)
  MIBOSS <- cbind(c(ga,0,m2(anb)),c(0,m2(anb),0),c(m2(anb),0,m4(anb)))
  dIBOSS <- det(MIBOSS)
  eff <- (dIBOSS / dstar)^(1/dim(MIBOSS)[1])
  return(eff)
}

######## plotting ###########

ga <- seq(0.0001,0.9999,0.0001)


eff <- vector()
for (alpha in ga) {
  eff <- c(eff, efficiency(alpha))
}

effIBOSS <- vector()
for (alpha in ga) {
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












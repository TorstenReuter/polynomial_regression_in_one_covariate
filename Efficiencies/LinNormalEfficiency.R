####################################################################
#### Table 6                                                    ####
#### Linear Regression, standar normal covariate                ####
#### Calculating the Efficiencies of uniform random subsampling ####
####################################################################

library(nleqslv)
library(ggplot2)




########## efficiency ######################
efficiency <- function(ga){
  Mstar <- cbind(c(ga,0),c(0,ga + sqrt(2/pi)*qnorm(1-ga/2)*exp(-(qnorm(1-ga/2)^2)/2)))
  dstar <- det(Mstar)
  d <- (ga^2)
  eff <- (d / dstar)^(1/dim(Mstar)[1])
  return(eff)
}


############ calculating efficiency for all ga ###############
ga <- seq(0.00001,0.9999,0.00001)

eff <- vector()
for (alpha in ga) {
  eff <- c(eff, efficiency(alpha))
}

dfeff <- data.frame("alpha" = ga, "Efficiency" = eff)

plot1 <- ggplot() +
  geom_line(data = dfeff, aes(x = alpha, y = eff), size = 1, linetype = "solid") +
  geom_hline(yintercept=1, linetype="dotted", size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("Efficiency"))) +
  ylim(0,1)

plot1


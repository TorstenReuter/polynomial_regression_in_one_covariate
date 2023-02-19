######################################################################
#### Table 5                                                      ####
#### Quadratic Regression, t_30 distributed covariate             ####
#### Finding the crossover point where the D-optimal subsampling  ####
#### design switches from 3 disjoint intervals in its support     ####
#### to only two intervals.                                       ####
#### This is analogous to Theorem 5.6                             ####
######################################################################

library(ggplot2)

# m2 in Information Matrix of design without interior interval
m2 <- function(x) {
  y <- (30/(28672*(30 + x^2)^(29/2)))*
    (-660382083462799072265625 * x^3 - 
       114466227800218505859375 * x^5 - 
       13081854605739257812500 * x^7 - 
       1065928893800976562500 * x^9 - 
       64601751139453125000 * x^11 - 
       2981619283359375000 * x^13 - 
       106013130075000000 * x^15 - 
       2910164355000000 * x^17 - 
       61266618000000 * x^19 - 
       972486000000 * x^21 - 
       11275200000 * x^23 - 
       90201600 * x^25 - 
       445440 * x^27 - 
       1024 * x^29 + 
       489776025600000000000000 * sqrt(30 + x^2) + 
       228562145280000000000000 * x^2 * sqrt(30 + x^2) + 
       49521798144000000000000 * x^4 * sqrt(30 + x^2) + 
       6602906419200000000000 * x^6 * sqrt(30 + x^2) + 
       605266421760000000000 * x^8 * sqrt(30 + x^2) + 
       40351094784000000000 * x^10 * sqrt(30 + x^2) + 
       2017554739200000000 * x^12 * sqrt(30 + x^2) + 
       76859228160000000 * x^14 * sqrt(30 + x^2) + 
       2241727488000000 * x^16 * sqrt(30 + x^2) + 
       49816166400000 * x^18 * sqrt(30 + x^2) + 
       830269440000 * x^20 * sqrt(30 + x^2) + 
       10063872000 * x^22 * sqrt(30 + x^2) + 
       83865600 * x^24 * sqrt(30 + x^2) + 
       430080 * x^26 * sqrt(30 + x^2) + 
       1024 * x^28 * sqrt(30 + x^2))
  return(y)
}


# m4 in Information Matrix of design without interior interval
m4 <- function(x) {
  y <- (1350/(372736*(30 + x^2)^(29/2)))*
    (-114466227800218505859375*x^5 - 
       13081854605739257812500*x^7 - 
       1065928893800976562500*x^9 - 
       64601751139453125000*x^11 - 
       2981619283359375000*x^13 - 
       106013130075000000*x^15 - 
       2910164355000000*x^17 - 
       61266618000000*x^19 - 
       972486000000*x^21 - 
       11275200000*x^23 - 
       90201600*x^25 - 
       445440*x^27 - 
       1024*x^29 + 
       489776025600000000000000*sqrt(30 + x^2) + 
       228562145280000000000000*x^2 * sqrt(30 + x^2) + 
       49521798144000000000000*x^4 * sqrt(30 + x^2) + 
       6602906419200000000000*x^6 * sqrt(30 + x^2) + 
       605266421760000000000*x^8 * sqrt(30 + x^2) + 
       40351094784000000000*x^10 * sqrt(30 + x^2) + 
       2017554739200000000*x^12 * sqrt(30 + x^2) + 
       76859228160000000*x^14 * sqrt(30 + x^2) + 
       2241727488000000*x^16 * sqrt(30 + x^2) + 
       49816166400000*x^18 * sqrt(30 + x^2) + 
       830269440000*x^20 * sqrt(30 + x^2) + 
       10063872000*x^22 * sqrt(30 + x^2) + 
       83865600*x^24 * sqrt(30 + x^2) +
       430080*x^26 * sqrt(30 + x^2) + 
       1024*x^28 * sqrt(30 + x^2))
  return(y)
}

# Vector of values for alpha
ga <- seq(0.01,0.99999,0.00001)

# The term c(alpha)
cfunc <- function(ga) {
  t <- qt(1 - ga/2, df = 30)
  c <- ga*t^2 * ((2*m2(t)- (ga*t^2))/(ga*m4(t) - m2(t)^2) - 1/m2(t))
}

c <- cfunc(ga)
dfc <- data.frame("alpha" = ga, "c" = c)


plot1 <- ggplot() +
  geom_line(data = dfc, aes(x = alpha, y = c), size = 1, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0.925835, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab(expression(paste(alpha))) + 
  ylab(expression(paste("c_30")))

plot1



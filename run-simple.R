# simple HyperHMM example
# NOTE! ensure that the HyperHMM executable has been compiled and is in the current working directory. if not, this code will try to compile it, but in a pretty dumb way!

library(ggplot2)

# source to wrap external HyperHMM call
source("hyperhmm-wrap.R")
# source for plots
source("hypercube-plots.R")

# create simple observation set
obs = matrix(data=c(0,0,1, 0,1,1), byrow=T, nrow=2)
# fit hypercubic model
fitted = hyperhmm(obs)
# summary plots
grid.arrange( plot.bubbles2(fitted[[1]], formatted=T),
              plot.hypercube2(fitted[[4]]),
              plot.pfg(fitted[[4]]),
              nrow = 1)

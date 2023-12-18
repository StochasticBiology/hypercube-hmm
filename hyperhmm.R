library(Rcpp)
library(RcppArmadillo)

sourceCpp("hyperhmm-r.cpp")

# source code for plots
library(ggpubr)
source("hypercube-plots.R")

# construct toy example
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
# do inference
fitted = HyperHMM(m)

# convert output transition set to previous format
str.trans = apply(fitted[[4]], 1, paste, collapse=" ")

# set of plots
plot.bubs = plot.bubbles2(fitted[[1]], formatted=TRUE)
plot.cube = plot.hypercube2(str.trans, use.width = T, node.labels=F, seg.labels = T, threshold=0, rotate.phi=F)
plot.diag = plot.pfg(str.trans, pfg.layout="matrix")

# arrange plots together
ggarrange(plot.bubs, plot.cube, plot.diag, nrow=1)

# some less trivial examples
m.1 = matrix(c(0,0,0,0,0,
                   1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   0,0,0,0,0,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1), byrow=TRUE, ncol = 5)

m.2 = matrix(c(1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   1,1,1,1,1,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1,
                   1,1,1,1,1), byrow=TRUE, ncol = 5)
             
HyperHMM(m.2)
HyperHMM(m.2, initialstates = m.1, outputinput = 1)

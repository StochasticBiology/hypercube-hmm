### simple demos of R embedding of HyperHMM

# source code for inference
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

# source code for plots
library(ggpubr)
source("hypercube-plots.R")

### first -- just do inference on a very simple example. we model two observations, 001 and 011
# construct toy example
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
# do inference
fitted = HyperHMM(m)

m = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0), byrow=TRUE, ncol=3)
fitted = HyperHMM(m)

# produce a set of plots. here we use a syntax to demonstrate back-compatibility
plot.bubs = plot.bubbles2(fitted[[1]], formatted=TRUE)
plot.cube = plot.hypercube2(fitted[[4]], use.width = T, node.labels=F, seg.labels = T, threshold=0, rotate.phi=F)
plot.diag = plot.pfg(fitted[[4]], pfg.layout="matrix")
# and here slightly more intuitively referencing the named elements of the returned list
plot.bubs = plot.bubbles2(fitted$stats, formatted=TRUE)
plot.cube = plot.hypercube2(fitted$viz, use.width = T, node.labels=F, seg.labels = T, threshold=0, rotate.phi=F)
plot.diag = plot.pfg(fitted$viz, pfg.layout="matrix")

# and arrange plots together
ggarrange(plot.bubs, plot.cube, plot.diag, nrow=1)

### some less trivial examples using the two-pathway synthetic model
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
             
fit.1 = HyperHMM(m.2)
plot.standard(fit.1)
fit.2 = HyperHMM(m.2, initialstates = m.1, outputinput = 1)
plot.standard(fit.2)

### now use some data -- from the ovarian cancer case study
cgh.raw = readLines("Data/ovarian.txt")
cgh.mat = do.call(rbind, lapply(strsplit(cgh.raw, ""), as.numeric))
fit.cgh = HyperHMM(cgh.mat)
plot.standard(fit.cgh)
# rather nice hypercube visualisation
plot.hypercube.flux(fit.cgh)

# this script does a no-frills plot set for precomputed HyperHMM output
# you'll need to provide the "label" of the experiment, and in the working directory the code will need the HyperHMM output:
# mean_[label].txt
# sd_[label].txt
# graph_viz_[label].txt
# transitions_[label].txt

source("hyperhmm-wrap.R")
source("hypercube-plots.R")

# replace this label with the label of the experiment you're plotting
label = "simple_case1_L5"

# read file output. calling this "wrapper" function with simulate = F just asks the code to retrieve precomputed output
fitted = hyperhmm(NULL, label=label, simulate = F)

# no-frills plot commands. these plots have some customisable stylings -- see hypercube-plots.R for the various options
plot.bubs = plot.bubbles2(fitted[[1]], formatted=T)
plot.cube = plot.hypercube2(fitted[[4]], use.width = T, node.labels=F, seg.labels = T, threshold=0, rotate.phi=F)
plot.diag = plot.pfg(fitted[[4]], pfg.layout="matrix")

# arrange plots together
grid.arrange(plot.bubs, plot.cube, plot.diag, nrow=1)

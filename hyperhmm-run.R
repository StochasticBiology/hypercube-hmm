### simple demos of R embedding of HyperHMM

# source code for inference
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

# source code for plots
library(ggpubr)
source("hypercube-plots.R")

# read in cross-sectional data and return a matrix
cube.read.crosssectional = function(fname) {
  data.raw = readLines(fname)
  data.mat = do.call(rbind, lapply(strsplit(data.raw, ""), as.numeric))
  return(data.mat)
}

# read in longitudinal data and return a list of two matrices
cube.read.longitudinal = function(fname) {
  data.list = list()
  data.raw = read.table(fname, header=FALSE, colClasses = "character")
  data.list$from = do.call(rbind, lapply(strsplit(data.raw[,1], ""), as.numeric))
  data.list$to = do.call(rbind, lapply(strsplit(data.raw[,2], ""), as.numeric))
  return(data.list)
}

### synthetic case studies
plot.out = fit.data = list()
files = c("simple_case1_L5", 
          "simple_case2_L5", "simple_case2_L7", "simple_case2_L9", 
          "simple_case4_L5", 
          "double_case2_L5", "double_case2_L7", "double_case2_L9")
for(expt in 1:length(files)) {
  print(paste0("Expt ", expt, ": ", files[expt]))
  data.mat = cube.read.crosssectional(paste0("Data/", files[expt], ".txt"))
  fit.data[[expt]] = HyperHMM(data.mat)
  plot.out[[expt]] = plot.standard(fit.data[[expt]])
}
# put all these plots together
ggarrange(plotlist = plot.out)

### synthetic dataset with high-order interactions
hiorder.mat = cube.read.crosssectional("Data/hi-order.txt")
fit.hiorder = HyperHMM(hiorder.mat)
plot.hiorder = plot.standard(fit.hiorder)

### ovarian cancer dataset
cgh.mat = cube.read.crosssectional("Data/ovarian.txt")
fit.cgh = HyperHMM(cgh.mat)
plot.cgh = plot.standard(fit.cgh)

### TB dataset
tb.list = cube.read.longitudinal("Data/tb_drug.txt")
fit.tb = HyperHMM(tb.list$to, initialstates=tb.list$from, nboot=1)
plot.tb = plot.standard(fit.tb)

# put these more specialised plots together
ggarrange(plot.hiorder, plot.cgh, plot.tb)

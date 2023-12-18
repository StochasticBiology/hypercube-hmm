# this script is part of a pair comparing HyperHMM and HyperTraPS for different datasets
# this one launches external processes doing the inference for each case
# when these have completed, run double-approach-retrieve.R to pull the data and plot

# pull the source needed to wrap external calls to the two approaches
source("hypertraps-wrap.R")
source("hyperhmm-wrap.R")

ht.length = 1       # length parameter for the HyperTraPS runs (10^(n+2))
start.hhmm = T      # actually set HyperHMM runs going?
start.ht = T        # actually set HyperTraPS runs going?
nboot = 10          # HyperHMM bootstraps
random.walkers = 0  # random walk for each?
param.kernel = 5    # HyperTraPS perturbation kernel

# the various experiments
filenames = c("double_case2_L5", "double_case2_L7", "double_case2_L9", "hi-order", "ovarian", "simple_case10_L5", "simple_case1_L5", "simple_case2_L5", "simple_case2_L7", "simple_case2_L9", "simple_case4_L5")

for(fname in filenames) {
  # read data from file into observations matrix
  rl = readLines(paste(c(fname, ".txt"), collapse=""))
  for(i in 1:length(rl)) {
    s = strsplit(rl[i], split=" ")
    a = strsplit(s[[1]][1], split="")[[1]]
    ad = unlist(lapply(a, as.numeric))
    if(i == 1) {
      obs = matrix(ad, byrow=T, nrow=1)
    } else {
      obs = rbind(obs, ad)
    }
  }
  mydata = obs
  rownames(mydata) = NULL
  precursors = NULL         # we'll treat longitudinal cases separately, after this loop
  
  label = fname
  ### HyperHMM
  if(start.hhmm == T) {
    # fork external process
    fitted = hyperhmm(as.matrix(mydata), precursors=NULL, label=label, nboot=100, random.walkers=1, fork=T)
  }
  
  ### HyperTraPS
  if(start.ht == T) {
    ht.simulate = "fork"
    # fork external processes
    post.sample = hypertraps(as.matrix(mydata), param.length=ht.length, param.kernel=param.kernel, label=label, simulate = ht.simulate)
  }
}

# now the longitudinal cases
fname = "tb_drug"
rl = readLines(paste(c(fname, ".txt"), collapse=""))
for(i in 1:length(rl)) {
  # read observations and precursors into matrices
  s = strsplit(rl[i], split=" ")
  a = strsplit(s[[1]][1], split="")[[1]]
  d = strsplit(s[[1]][2], split="")[[1]]
  ad = unlist(lapply(a, as.numeric))
  dd = unlist(lapply(d, as.numeric))
  if(i == 1) {
    prec = matrix(ad, byrow=T, nrow=1)
    obs = matrix(dd, byrow=T, nrow=1)
  } else {
    prec = rbind(prec, ad)
    obs = rbind(obs, dd)
  }
}
label = fname
# fork the external code
if(start.hhmm == T) {
  fitted = hyperhmm(as.matrix(obs), precursors=as.matrix(prec), label=label, nboot=3, random.walkers=1, fork=T)
}
if(start.ht == T) {
  param.kernel = 5
  ht.simulate = "fork"
  
  post.sample = hypertraps(as.matrix(obs), precursors=as.matrix(prec), param.length=ht.length, param.kernel=param.kernel, label=label, simulate = ht.simulate)
}

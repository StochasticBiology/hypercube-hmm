# this script compares different pathway inference approaches for an ovarian cancer dataset
# the script should work straightforwardly *if* the required packages are installed. "Oncotree" is a standard one and should be easy
# "TRONCO" may be a bit harder. it is available via BioConductor; see https://bioconductor.org/packages/release/bioc/html/TRONCO.html. the installation idea is
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("TRONCO")
# it may also require the installation of some libraries external to R. a dependency that caused us some trouble was the R library "cgdsr", which requires R's "devtools"
# the R "wrapper" for HyperHMM may also be a bit awkward. it needs the executable HyperHMM file (compiled from C++) to be present in the working directory -- see the comments therein

hypertraps.precomputed = F
source("hypertraps-wrap.R")
source("hyperhmm-wrap.R")
ht.length = 1
start.hhmm = F
start.ht = T

filenames = c("double_case2_L5", "double_case2_L7", "double_case2_L9", "hi-order", "ovarian", "simple_case10_L5", "simple_case1_L5", "simple_case2_L5", "simple_case2_L7", "simple_case2_L9", "simple_case4_L5")

for(fname in filenames) {
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

label = fname
### HyperHMM
if(start.hhmm == T) {
nboot = 10
random.walkers = 0
precursors = NULL

fitted = hyperhmm(as.matrix(mydata), precursors=NULL, label=label, nboot=100, random.walkers=1, fork=T)
}

### HyperTraPS

if(start.ht == T) {


param.kernel = 5
precursors = NULL
ht.simulate = "fork"

post.sample = hypertraps(as.matrix(mydata), param.length=ht.length, param.kernel=param.kernel, label=label, simulate = ht.simulate)
}
}

fname = "tb_drug"
rl = readLines(paste(c(fname, ".txt"), collapse=""))
for(i in 1:length(rl)) {
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
if(start.hhmm == T) {
fitted = hyperhmm(as.matrix(obs), precursors=as.matrix(prec), label=label, nboot=100, random.walkers=1, fork=T)
}
if(start.ht == T) {
param.kernel = 5
ht.simulate = "fork"

post.sample = hypertraps(as.matrix(obs), precursors=as.matrix(prec), param.length=ht.length, param.kernel=param.kernel, label=label, simulate = ht.simulate)
}
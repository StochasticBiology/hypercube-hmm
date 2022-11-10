# this script compares different pathway inference approaches for an ovarian cancer dataset
# the script should work straightforwardly *if* the required packages are installed. "Oncotree" is a standard one and should be easy
# "TRONCO" may be a bit harder. it is available via BioConductor; see https://bioconductor.org/packages/release/bioc/html/TRONCO.html. the installation idea is
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("TRONCO")
# it may also require the installation of some libraries external to R. a dependency that caused us some trouble was the R library "cgdsr", which requires R's "devtools"
# the R "wrapper" for HyperHMM may also be a bit awkward. it needs the executable HyperHMM file (compiled from C++) to be present in the working directory -- see the comments therein

library(Oncotree)
library(gridExtra)
library(tictoc)

hypertraps.precomputed = T
mhn.working = F
nwalk = 5000
ht.length = 1

filenames = c("double_case2_L5", "hi-order", "ovarian")

for(fname in filenames) {
  expt = fname
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
L = ncol(mydata)
rownames(mydata) = NULL

  ### TRONCO approaches 

library(TRONCO)

# format for TRONCO approaches
dataset_plain = import.genotypes(mydata, event.type='point')
dataset_plain = change.color(dataset_plain, 'point', "Lightblue")
alterations = events.selection(as.alterations(dataset_plain), filter.freq = .05)
gene.hypotheses = dataset_plain$types
dataset_plain.clean = events.selection(dataset_plain, filter.in.names=c(as.genes(alterations), gene.hypotheses))
print(dataset_plain.clean)

# TRONCO fits
tic(paste(c("Capri ", expt), collapse=""))
model.capri = tronco.capri(dataset_plain.clean)
toc(log=TRUE)
model.caprese = tronco.caprese(dataset_plain.clean)

# TRONCO plots
tronco.plot(model.capri)
tronco.plot(model.caprese)

### Oncotree

# fit and bootstrap
tic(paste(c("Oncotree ", expt), collapse=""))
ov.tree <- oncotree.fit(mydata)
toc(log=TRUE)
#ov.b1 <- bootstrap.oncotree(ov.tree, R=100, type="parametric")
#ov.b1

# outputs
#opar <- par(mfrow=c(3,2), mar=c(2,0,0,0))
#plot(ov.b1, nboots=4)
#plot(ov.b1, nboots=4, fix.nodes=TRUE)
#par(opar)

### Mutual Hazard Networks

# required source files

if(mhn.working == T){
setwd("MHN/")
source("UtilityFunctions.R")
source("ModelConstruction.R")
source("Likelihood.R")
source("RegularizedOptimization.R")
setwd("../")

# convert data to required format
pD <- Data.to.pD(mydata)

# fit
Theta.BC <- Learn.MHN(pD, lambda=0.01)

# output
colnames(Theta.BC) <- colnames(mydata)
rownames(Theta.BC) <- colnames(Theta.BC)
View(exp(Theta.BC))
} else {
  tmp = read.csv(paste(c("out-mhn-", fname, ".csv"), collapse=""), header=T)
  Theta.BC = log(as.matrix(tmp[,2:ncol(tmp)]))
}

### HyperHMM

source("hyperhmm-wrap.R")

label = fname
nboot = 10
random.walkers = 0
precursors = NULL

tic(paste(c("HyperHMM ", expt), collapse=""))
fitted = hyperhmm(as.matrix(mydata), label=label, simulate = F)
toc(log=TRUE)

### HyperTraPS

source("hypertraps-wrap.R")

label = fname
param.kernel = 5
precursors = NULL
ht.simulate = F


tic(paste(c("HyperTraPS ", expt), collapse=""))
post.sample = hypertraps(as.matrix(mydata), param.length=ht.length, param.kernel=param.kernel, label=label, simulate = ht.simulate)
toc(log=TRUE)

# simulate MHN
trans.list.mhn = c()
bubble.probs.mhn = matrix(0, nrow=L, ncol=L)
for(i in 1:nwalk) {
  state = rep(0,L)
  for(this.step in 1:L) {
    probs = diag(Theta.BC)
    for(j in 1:L) {
      for(k in 1:L) {
        probs[j] = probs[j]+state[k]*Theta.BC[j,k]
      }
    }
    if(this.step == L) { next.step = which(state == 0)}
    else {
    next.step = sample(which(state == 0), size=1, prob=exp(probs[which(state==0)]))
    }
    old.state = state
    state[next.step] = 1
    bubble.probs.mhn[next.step,this.step] = bubble.probs.mhn[next.step,this.step]+1
    trans.list.mhn = c(trans.list.mhn, paste(paste(old.state, collapse=""), paste(state, collapse="")))
  }
}

# simulate OncoTree
parents = ov.tree$parent$parent.num
trans.list.ot = c()
bubble.probs.ot = matrix(0, nrow=L, ncol=L)
for(i in 1:nwalk) {
  state = rep(0,L)
  for(this.step in 1:L) {
    possibles = c()
    for(j in 1:L) {
      required = parents[j+1]-1
      if(required == 0) {
        if(state[j] == 0) {
          possibles = c(possibles, j)
        }
      } else {
        if(state[j] == 0 & state[required] == 1) {
          possibles = c(possibles, j)
        }
      }
    }
    next.step = sample(possibles, size=1)
    old.state = state
    state[next.step] = 1
    bubble.probs.ot[next.step,this.step] = bubble.probs.ot[next.step,this.step]+1
    trans.list.ot = c(trans.list.ot, paste(paste(old.state, collapse=""), paste(state, collapse="")))
  }
}

# simulate Capri
parents = model.capri$model$capri_aic$parents.pos
trans.list.capri = c()
bubble.probs.capri = matrix(0, nrow=L, ncol=L)
for(i in 1:nwalk) {
  state = rep(0,L)
  for(this.step in 1:L) {
    possibles = c()
    for(j in 1:L) {
      required = unlist(parents[j])
      if(required[1] == -1) {
        if(state[j] == 0) {
          possibles = c(possibles, j)
        }
      } else {
        if(state[j] == 0 & !any(state[required]) == 0) {
          possibles = c(possibles, j)
        }
      }
    }
    next.step = sample(possibles, size=1)
    old.state = state
    state[next.step] = 1
    bubble.probs.capri[next.step,this.step] = bubble.probs.capri[next.step,this.step]+1
    trans.list.capri = c(trans.list.capri, paste(paste(old.state, collapse=""), paste(state, collapse="")))
   # print(trans.list.capri[length(trans.list.capri)])
  }
}

# simulate HyperTraPS
trans.list.ht = c()
bubble.probs.ht = matrix(0, nrow=L, ncol=L)
for(post in 1:nrow(post.sample)) {
for(i in 1:(max(ceiling(nwalk/nrow(post.sample)), 1))) {
  post.matrix = t(matrix(unlist(post.sample[post,1:(L*L)]), nrow=L))
  
  state = rep(0,L)
  for(this.step in 1:L) {
    probs = diag(post.matrix)
    for(j in 1:L) {
      for(k in 1:L) {
        probs[j] = probs[j]+state[k]*post.matrix[j,k]
      }
    }
    if(this.step == L) { next.step = which(state == 0)}
    else {
    next.step = sample(which(state == 0), size=1, prob=exp(probs[which(state==0)]))
    }
    old.state = state
    state[next.step] = 1
    bubble.probs.ht[next.step,this.step] = bubble.probs.ht[next.step,this.step]+1
    trans.list.ht = c(trans.list.ht, paste(paste(old.state, collapse=""), paste(state, collapse="")))
  }
}
}


source("hypercube-plots.R")
if(expt == "ovarian") {sl = F} else {sl = T}
o.1 = plot.bubbles2(bubble.probs.capri) + scale_size(range = c(0,20))
o.2 = plot.bubbles2(bubble.probs.ot) + scale_size(range = c(0,20))
o.3 = plot.bubbles2(bubble.probs.mhn) + scale_size(range = c(0,20))
o.4 = plot.bubbles2(bubble.probs.ht) + scale_size(range = c(0,20))
o.5 = plot.bubbles2(fitted[[1]], formatted=T) + scale_size(range = c(0,20))
p.1 = plot.hypercube2(trans.list.capri, use.width= T, node.labels=F, seg.labels = sl, threshold=0)
p.2 = plot.hypercube2(trans.list.ot, use.width= T, node.labels=F, seg.labels = sl, threshold=0)
p.3 = plot.hypercube2(trans.list.mhn, use.width= T, node.labels=F, seg.labels = sl, threshold=0)
p.4 = plot.hypercube2(trans.list.ht, use.width= T, node.labels=F, seg.labels = sl, threshold=0)
p.5 = plot.hypercube2(fitted[[4]], use.width = T, node.labels=F, seg.labels = sl, threshold=0)
q.1 = plot.pfg(trans.list.capri, pfg.layout="matrix")
q.2 = plot.pfg(trans.list.ot, pfg.layout="matrix")
q.3 = plot.pfg(trans.list.mhn, pfg.layout="matrix")
q.4 = plot.pfg(trans.list.ht, pfg.layout="matrix")
q.5 = plot.pfg(fitted[[4]], pfg.layout="matrix")

png(paste(c("comp-", fname, ".png"), collapse=""), width=1200*sf, height=800*sf, res=72*sf)
grid.arrange(o.1, o.2, o.3, o.4, o.5, p.1, p.2, p.3, p.4, p.5, q.1, q.2, q.3, q.4, q.5, nrow=3)
dev.off()
}
tic.log(format=TRUE)
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

# source for plots
source("hypercube-plots.R")

hypertraps.precomputed = T   # have we precomputed the HyperTraPS output, or should we run it now? recommended to precompute
ht.length = 3                # MCMC length parameter for HyperTraPS output
mhn.working = F              # is the MHN code working, or should we try and retrieve pre-saved outputs?
nwalk = 5000                 # number of walkers to simulate in summarising outputs
bubscale = 12                # scale for bubble plots
sf = 3                       # resolution factor for plots

# experiments to compare
filenames = c("double_case2_L5", "hi-order", "ovarian")

# loop through experiments
for(fname in filenames) {
  # read in data
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
  # required source files if the code is working
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
    # if the code isn't working, try and retrieve precomputed values
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
  
  #### there follows code to simulate walkers on the output of each approach, to compare them. 
  
  # simulate MHN
  # initialise records
  trans.list.mhn = c()
  bubble.probs.mhn = matrix(0, nrow=L, ncol=L)
  # loop through walkers
  for(i in 1:nwalk) {
    # initialise state
    state = rep(0,L)
    # loop through steps on the hypercube
    for(this.step in 1:L) {
      # work out transition probabilities given this state
      probs = diag(Theta.BC)
      for(j in 1:L) {
        for(k in 1:L) {
          probs[j] = probs[j]+state[k]*Theta.BC[j,k]
        }
      }
      # choose next step
      if(this.step == L) { next.step = which(state == 0)}
      else {
        next.step = sample(which(state == 0), size=1, prob=exp(probs[which(state==0)]))
      }
      # record this step
      old.state = state
      state[next.step] = 1
      bubble.probs.mhn[next.step,this.step] = bubble.probs.mhn[next.step,this.step]+1
      trans.list.mhn = c(trans.list.mhn, paste(paste(old.state, collapse=""), paste(state, collapse="")))
    }
  }
  
  # simulate OncoTree
  parents = ov.tree$parent$parent.num
  # initialise records
  trans.list.ot = c()
  bubble.probs.ot = matrix(0, nrow=L, ncol=L)
  # loop through walkers
  for(i in 1:nwalk) {
    # initialise state
    state = rep(0,L)
    # loop through steps on the hypercube
    for(this.step in 1:L) {
      possibles = c()
      # work out which steps are possible given this state, from parent structure of output
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
      # choose next step
      next.step = sample(possibles, size=1)
      # record this step
      old.state = state
      state[next.step] = 1
      bubble.probs.ot[next.step,this.step] = bubble.probs.ot[next.step,this.step]+1
      trans.list.ot = c(trans.list.ot, paste(paste(old.state, collapse=""), paste(state, collapse="")))
    }
  }
  
  # simulate Capri
  parents = model.capri$model$capri_aic$parents.pos
  # initialise records
  trans.list.capri = c()
  bubble.probs.capri = matrix(0, nrow=L, ncol=L)
  # loop through walkers
  for(i in 1:nwalk) {
    # initialise state
    state = rep(0,L)
    # loop through steps on the hypercube
    for(this.step in 1:L) {
      # work out which steps are possible given this state (from parent structure of output)
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
      # choose next step
      next.step = sample(possibles, size=1)
      # record this step
      old.state = state
      state[next.step] = 1
      bubble.probs.capri[next.step,this.step] = bubble.probs.capri[next.step,this.step]+1
      trans.list.capri = c(trans.list.capri, paste(paste(old.state, collapse=""), paste(state, collapse="")))
      # print(trans.list.capri[length(trans.list.capri)])
    }
  }
  
  # simulate HyperTraPS
  # initialise records
  trans.list.ht = c()
  bubble.probs.ht = matrix(0, nrow=L, ncol=L)
  # loop through posterior samples
  for(post in 1:nrow(post.sample)) {
    for(i in 1:(max(ceiling(nwalk/nrow(post.sample)), 1))) {
      # interpret output as matrix
      post.matrix = t(matrix(unlist(post.sample[post,1:(L*L)]), nrow=L))
      # initialise state
      state = rep(0,L)
      # loop through steps on the hypercube
      for(this.step in 1:L) {
        # work out probabilities given this state
        probs = diag(post.matrix)
        for(j in 1:L) {
          for(k in 1:L) {
            probs[j] = probs[j]+state[k]*post.matrix[j,k]
          }
        }
        # choose next step
        if(this.step == L) { next.step = which(state == 0)}
        else {
          next.step = sample(which(state == 0), size=1, prob=exp(probs[which(state==0)]))
        }
        # record this step
        old.state = state
        state[next.step] = 1
        bubble.probs.ht[next.step,this.step] = bubble.probs.ht[next.step,this.step]+1
        trans.list.ht = c(trans.list.ht, paste(paste(old.state, collapse=""), paste(state, collapse="")))
      }
    }
  }
  
  # experiments-specific plot styling
  if(expt == "ovarian") {sl = F} else {sl = T}
  if(expt == "hi-order") { rotate.phi = T } else {rotate.phi = F}
  
  # enforce integer ticks for bubble plots
  int.y = scale_y_continuous(breaks = function(x) unique(floor(seq(0, (max(x) + 1) * 1.1))))
  int.x = scale_x_continuous(breaks = function(x) unique(floor(seq(0, (max(x) + 1) * 1.1))))
  
  # construct various summary plots
  o.1 = plot.bubbles2(bubble.probs.capri) + scale_size(range = c(0,bubscale))  + int.x + int.y
  o.2 = plot.bubbles2(bubble.probs.ot) + scale_size(range = c(0,bubscale))  + int.x + int.y
  o.3 = plot.bubbles2(bubble.probs.mhn) + scale_size(range = c(0,bubscale))  + int.x + int.y
  o.4 = plot.bubbles2(bubble.probs.ht) + scale_size(range = c(0,bubscale))  + int.x + int.y
  o.5 = plot.bubbles2(fitted[[1]], formatted=T) + scale_size(range = c(0,bubscale))   + int.x + int.y
  p.1 = plot.hypercube2(trans.list.capri, use.width= T, node.labels=F, seg.labels = sl, threshold=0, rotate.phi=rotate.phi)
  p.2 = plot.hypercube2(trans.list.ot, use.width= T, node.labels=F, seg.labels = sl, threshold=0, rotate.phi=rotate.phi)
  p.3 = plot.hypercube2(trans.list.mhn, use.width= T, node.labels=F, seg.labels = sl, threshold=0, rotate.phi=rotate.phi)
  p.4 = plot.hypercube2(trans.list.ht, use.width= T, node.labels=F, seg.labels = sl, threshold=0, rotate.phi=rotate.phi)
  p.5 = plot.hypercube2(fitted[[4]], use.width = T, node.labels=F, seg.labels = sl, threshold=0, rotate.phi=rotate.phi)
  q.1 = plot.pfg(trans.list.capri, pfg.layout="matrix")+ scale_edge_width(range = c(0,2))
  q.2 = plot.pfg(trans.list.ot, pfg.layout="matrix")+ scale_edge_width(range = c(0,2))
  q.3 = plot.pfg(trans.list.mhn, pfg.layout="matrix")+ scale_edge_width(range = c(0,2))
  q.4 = plot.pfg(trans.list.ht, pfg.layout="matrix")+ scale_edge_width(range = c(0,2))
  q.5 = plot.pfg(fitted[[4]], pfg.layout="matrix")+ scale_edge_width(range = c(0,2))
  
  # panel labelling layer
  labmod= theme(plot.tag=element_text(size=7*sf))
  
  # write trellis plot to file
  png(paste(c("comp-", fname, ".png"), collapse=""), width=800*sf, height=1200*sf, res=72*sf)
  grid.arrange(o.1+labs(tag = "A.i")+labmod, p.1+labs(tag = "ii")+labmod, q.1+labs(tag = "iii")+labmod,
               o.2+labs(tag = "B.i")+labmod, p.2+labs(tag = "ii")+labmod, q.2+labs(tag = "iii")+labmod,
               o.3+labs(tag = "C.i")+labmod, p.3+labs(tag = "ii")+labmod, q.3+labs(tag = "iii")+labmod,
               o.4+labs(tag = "D.i")+labmod, p.4+labs(tag = "ii")+labmod, q.4+labs(tag = "iii")+labmod,
               o.5+labs(tag = "E.i")+labmod, p.5+labs(tag = "ii")+labmod, q.5+labs(tag = "iii")+labmod,
               nrow=5)
  dev.off()
  
  
}
tic.log(format=TRUE)

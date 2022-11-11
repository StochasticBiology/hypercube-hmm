# this script is part of a pair comparing HyperHMM and HyperTraPS for different datasets
# this one pulls the output from previously-run processes and plots them
# first run double-approach-start.R to set these going

# code for wrapping approaches (here needed for retrieval)
source("hypertraps-wrap.R")
source("hyperhmm-wrap.R")

# source for plots
source("hypercube-plots.R")

# HyperTraPS length parameter -- make sure this matches the one in double-approach-start.R
ht.length = 3

nwalk = 1000   # number of walkers to simulate for summarised output
sf = 3         # resolution factor for plots
bubscale = 10  # scale for bubble plots

filenames = c("double_case2_L5", "double_case2_L7", "double_case2_L9", "hi-order", "ovarian", "simple_case10_L5", "simple_case1_L5", "simple_case2_L5", "simple_case2_L7", "simple_case2_L9", "simple_case4_L5", "tb_drug")

# initialise lists for all plots (to help with construction big figures)
bub.ht = cube.ht = pfg.ht = bub.hm = cube.hm = pfg.hm = list()
index = 1

# loop through experiments
for(fname in filenames) {
  # read in data. don't need this except to determine L
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
  
  ### HyperHMM pull (no simulation)
  label = fname
  fitted = hyperhmm(as.matrix(mydata), simulate=F, label=label)
  ### HyperTraPS pull (no simulation)
  post.sample = hypertraps(as.matrix(mydata), param.length=ht.length, param.kernel=5, label=label, simulate = "no")
  
  # simulate walkers on HyperTraPS posteriors to get transition list
  # initialist transition list and bubble matrix
  trans.list.ht = c()
  bubble.probs.ht = matrix(0, nrow=L, ncol=L)
  # loop through posterior samples
  for(post in 1:nrow(post.sample)) {
    for(i in 1:(max(ceiling(nwalk/nrow(post.sample)), 1))) {
      # interpret output sample as matrix
      post.matrix = t(matrix(unlist(post.sample[post,1:(L*L)]), nrow=L))
      #initialise simulated state
      state = rep(0,L)
      # loop through steps
      for(this.step in 1:L) {
        probs = diag(post.matrix)
        # compute transition probabilities given current state
        for(j in 1:L) {
          for(k in 1:L) {
            probs[j] = probs[j]+state[k]*post.matrix[j,k]
          }
        }
        # pick next step
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
  
  # experiment-specific plot characteristics
  if(fname == "tb_drug") { sl = F } else { sl = T }
  if(L <= 5) { nl = T; ps = 1} else {nl = F; ps=0.2}
  if(fname == "hi-order") { rotate.phi = T } else {rotate.phi = F}
  threshold = 20
  # enforce integer axis labels for bubble plots
  int.y = scale_y_continuous(breaks = function(x) unique(floor(seq(0, (max(x) + 1) * 1.1))))
  int.x = scale_x_continuous(breaks = function(x) unique(floor(seq(0, (max(x) + 1) * 1.1))))
  
  # construct the various plots for HyperTraPS and HyperHMM
  bub.ht[[index]] = o.4 = plot.bubbles2(bubble.probs.ht) + scale_size(range = c(0,bubscale)) + int.x + int.y
  cube.ht[[index]] = p.4 = plot.hypercube2(trans.list.ht, use.width= T, node.labels=nl, seg.labels = sl, p.size = ps, threshold=threshold, rotate.phi = rotate.phi)
  pfg.ht[[index]] = q.4 = plot.pfg(trans.list.ht, pfg.layout="matrix") + scale_edge_width(range = c(0,2))
  bub.hm[[index]] = o.5 = plot.bubbles(fitted[[1]]) + scale_size(range = c(0,bubscale))+ theme(legend.position="none") + int.x + int.y
  cube.hm[[index]] = p.5 = plot.hypercube2(fitted[[4]], use.width = T, node.labels=nl, seg.labels = sl, p.size = ps, threshold=threshold, rotate.phi = rotate.phi)
  pfg.hm[[index]] = q.5 = plot.pfg(fitted[[4]], pfg.layout="matrix") + scale_edge_width(range = c(0,2))
  
  index = index+1
  # labelling layer
  labmod= theme(plot.tag=element_text(size=7*sf))
  
  # write plot to file
  png(paste(c(fname, "-2.png"), collapse=""), width=800*sf, height=400*sf, res=72*sf)
  grid.arrange(o.4+labs(tag = "A.i")+labmod, p.4+labs(tag = "ii")+labmod, q.4+labs(tag = "iii")+labmod, 
               o.5+labs(tag = "B.i")+labmod, p.5+labs(tag = "ii")+labmod, q.5+labs(tag = "iii")+labmod, nrow=2)
  dev.off()
}

# manuscript plots
png("old-first-fig.png", width=600*sf, height=480*sf, res=72*sf)
grid.arrange(bub.hm[[8]]+labs(tag = "A.i")+labmod, cube.hm[[8]]+labs(tag = "ii")+labmod, pfg.hm[[8]]+labs(tag = "iii")+labmod, 
             bub.hm[[1]]+labs(tag = "B.i")+labmod, cube.hm[[1]]+labs(tag = "ii")+labmod, pfg.hm[[1]]+labs(tag = "iii")+labmod, 
             bub.hm[[4]]+labs(tag = "C.i")+labmod, cube.hm[[4]]+labs(tag = "ii")+labmod, pfg.hm[[4]]+labs(tag = "iii")+labmod, nrow = 3)
dev.off()

png("old-second-fig.png", width=600*sf, height=800*sf, res=72*sf)
grid.arrange(bub.ht[[6]]+labs(tag = "A.i")+labmod, bub.ht[[7]]+labs(tag = "ii")+labmod, bub.ht[[11]]+labs(tag = "iii")+labmod, 
             bub.hm[[6]]+labs(tag = "B.i")+labmod, bub.hm[[7]]+labs(tag = "ii")+labmod, bub.hm[[11]]+labs(tag = "iii")+labmod, 
             bub.ht[[1]]+labs(tag = "C.i")+labmod, bub.ht[[2]]+labs(tag = "ii")+labmod, bub.ht[[3]]+labs(tag = "iii")+labmod, 
             bub.hm[[1]]+labs(tag = "D.i")+labmod, bub.hm[[2]]+labs(tag = "ii")+labmod, bub.hm[[3]]+labs(tag = "iii")+labmod,nrow=4)              
dev.off()


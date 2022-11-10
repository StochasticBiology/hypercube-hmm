#  single approach cases

library(tictoc)

source("hypertraps-wrap.R")
source("hyperhmm-wrap.R")
nwalk = 1000
sf = 3
ht.length = 1

filenames = c("double_case2_L5", "double_case2_L7", "double_case2_L9", "hi-order", "ovarian", "simple_case10_L5", "simple_case1_L5", "simple_case2_L5", "simple_case2_L7", "simple_case2_L9", "simple_case4_L5")
#ilenames = filenames[1]

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
  L = ncol(mydata)
  rownames(mydata) = NULL
  
  ### HyperHMM
  
  source("hyperhmm-wrap.R")
  
  label = fname
  nboot = 10
  random.walkers = 0
  precursors = NULL
  
  fitted = hyperhmm(as.matrix(mydata), simulate=F, label=label)
  
  post.sample = hypertraps(as.matrix(mydata), param.length=ht.length, param.kernel=5, label=label, simulate = "no")
  
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
  sl = T
  o.4 = plot.bubbles2(bubble.probs.ht) + scale_size(range = c(0,20))
  p.4 = plot.hypercube2(trans.list.ht, use.width= T, node.labels=F, seg.labels = sl, threshold=0)
  q.4 = plot.pfg(trans.list.ht, pfg.layout="matrix")
  o.5 = plot.bubbles(fitted[[1]]) + scale_size(range = c(0,20))+ theme(legend.position="none")
  p.5 = plot.hypercube2(fitted[[4]], use.width = T, node.labels=F, seg.labels = sl, threshold=0)
  q.5 = plot.pfg(fitted[[4]], pfg.layout="matrix")
  
  labmod= theme(plot.tag=element_text(size=7*sf))
  
  png(paste(c(fname, "-2.png"), collapse=""), width=800*sf, height=400*sf, res=72*sf)
  #grid.arrange(o.4, o.5, p.4, p.5, q.4, q.5, nrow=3)
  grid.arrange(o.4+labs(tag = "A.i")+labmod, p.4+labs(tag = "ii")+labmod, q.4+labs(tag = "iii")+labmod, 
               o.5+labs(tag = "B.i")+labmod, p.5+labs(tag = "ii")+labmod, q.5+labs(tag = "iii")+labmod, nrow=2)
  dev.off()
}


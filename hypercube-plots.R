# various plots for the output of hypercubic inference (or other approaches coerced to hypercubic format)

library(stringr)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(gridExtra)
library(igraph)
library(gtools)

# redundant? old, less flexible bubble plot function
plot.bubbles = function(stats.df, labels=NULL) {
  message("Building bubble plot")
  ## bubbles
  if(is.null(labels)) {
    g.1 = ggplot(stats.df, aes(x=order, y=feature)) + geom_point(aes(size=mean), colour="#CCCCCC") + geom_point(aes(size=sd), shape = 1, colour="#444444") + theme_classic() 
  } else {
    g.1 = ggplot(stats.df, aes(x=order, y=feature)) + geom_point(aes(size=mean), colour="#CCCCCC") + geom_point(aes(size=sd), shape = 1, colour="#444444") + scale_y_continuous(breaks=length(fitted[[3]]):1, labels=fitted[[3]]) + theme_classic() 
  }
  g.1
}

plot.bubbles2 = function(bp,             # output data structure. either just a matrix of probabilities (formatted == F) or a dataframe output from the HyperHMM wrapper with means and sds (formatted == T)
                         labels=NULL,    # labels for feature names
                         formatted=F     # dataframe formatted or not? see above
) {
  message("Building bubble plot")
  
  # if we've just got a matrix of probabilities, pull it into long form
  if(formatted == F) {
    bp.df = data.frame()
    for(i in 1:nrow(bp)) {
      for(j in 1:ncol(bp)) {
        bp.df = rbind(bp.df, data.frame(order=j, feature=i, prob=bp[i,j]))
      }
    }
  } else { bp.df = bp; bp.df$prob = bp.df$mean  } 
  ## plot bubbles
  if(is.null(labels)) {
    g.1 = ggplot(bp.df, aes(x=order, y=feature)) + geom_point(aes(size=prob), colour="#CCCCCC") + theme_classic() + theme(legend.position = "none")
  } else {
    g.1 = ggplot(bp.df, aes(x=order, y=feature)) + geom_point(aes(size=prob), colour="#CCCCCC") + scale_y_continuous(breaks=length(labels):1, labels=labels) + theme_classic() + theme(legend.position = "none")
  }
  return(g.1)
}

# redundant? older less flexible version -- will be sparsely commented
plot.hypercube = function(translist, use.width = T,           # use line width to display edge weights?
                          duplicate.offset = 0.,   # vertical offset for nodes in identical positions
                          lab.size = 3,            # size for edge labels
                          p.size = 1,              # point size
                          node.labels = T,         # node labels, yes or no?
                          threshold = 0,           # ignore edges under a threshold in the hypercube plot
                          break.redundancy = 0) {
  
  message("Building hypercube plot")  
  ## hypercube
  # get unique set of transitions and associated counts
  l = unique(translist)
  counts = rep(0, length(l))
  for(i in 1:length(l)) {
    set = which(translist == l[i])
    counts[i] = length(set)
  }
  #  l = l[counts > threshold]
  
  # split into lists of source and destination nodes
  srcs = dests = list()
  n = 1
  for(line in l) {
    s = strsplit(line, " ")
    srcs[[n]] = s[[1]][1]
    dests[[n]] = s[[1]][2]
    n = n + 1
  }
  # set string length and 0^L string
  len = nchar(srcs[[1]]) 
  zero = paste(rep("0", len), collapse="")
  
  # produce useful vectors
  srcs = unlist(srcs)
  dests = unlist(dests)
  nodes = unique(union(srcs, dests))
  nnodes = length(nodes)
  
  # produce list storing where incoming edges to each node come from
  ins = list()
  for(node in nodes) {
    refs = which(dests == node)
    refcodes = srcs[refs]
    ins[[length(ins)+1]] = which(nodes %in% refcodes)
  }
  
  ########### first produce hypercube visualisation
  
  message("Calculating embedding")
  # spherical polars: r, theta, phi
  # r = 1 everywhere
  rs = rep(1, nnodes)
  # theta is just set by number of 1s in a string
  thetas = unlist(lapply(nodes, function(s) { return(str_count(s, "0")*3.14159/len) }) )
  
  # initialise phis
  phis = rep(-1, nnodes)
  
  # phi for 0^L is zero; phis for the first set of nodes are evenly distributed over cos(0, pi)
  phis[which(nodes == zero)] = 0
  first.nodes = sort(dests[which(srcs == zero)])
  refs = c()
  for(i in 1:length(first.nodes)) {
    refs = c(refs, which(nodes == first.nodes[i]))
  }
  for(i in 1:length(refs)) {
    if(i > 1) {
      phis[refs[i]] = (i-1)*3.14159/(length(refs)-1)
    } else {
      phis[refs[i]] = 0
    }
  }
  
  message("Iterating")
  # iterate while we still have phi values to find
  change = T
  while(change == T) {
    change = F
    for(i in 1:nnodes) {
      if(phis[i] == -1) {
        # we haven't got a phi for this node yet
        if(!(-1 %in% unlist(phis[ins[[i]]]))) {
          # we have got phis for all its ancestors
          # this node's phi is set to the mean of its ancestors'
          phis[i] = mean(unlist(phis[ins[[i]]]))
          # remember we made a change in this iteration, so keep going
          change = T
        }
      }
    }
  }
  
  message("More coords")
  # dataframes for spherical and cartesian coordinates
  spcoords = data.frame(r = rs, theta = thetas, phi = phis, label = nodes)
  
  transcoords = data.frame()
  for(this.theta in unique(spcoords$theta)) {
    subset = spcoords[spcoords$theta == this.theta,]
    this.rank = rank(subset$phi, ties.method="average")
    subset$phi = 3.14159*(this.rank-1)/(max(this.rank)-1)
    if(nrow(subset) == 1) { subset$phi = 0 }
    transcoords = rbind(transcoords, subset)
  }
  
  if(break.redundancy) {
    spcoords = transcoords
  }
  
  coords = data.frame(x = spcoords$r*cos(spcoords$phi)*sin(spcoords$theta), y = spcoords$r*sin(spcoords$phi)*sin(spcoords$theta), z = spcoords$r*cos(spcoords$theta), label = spcoords$label)
  #coords = coords[nodes.count > threshold,]
  
  # dataframe for line segments in cartesian coords
  segments = data.frame()
  seglabels = data.frame()
  safe.nodes = rep(F, nrow(coords))
  for(i in 1:length(srcs)) {
    if(counts[i] > threshold) {
      src = which(coords$label == srcs[i])
      dest = which(coords$label == dests[i])
      safe.nodes[src] = safe.nodes[dest] = T
      label = paste(c("+", which(unlist(strsplit(srcs[i], split="")) !=unlist(strsplit(dests[i], split="")))), collapse="")
      segment = data.frame(src.x = coords[src,]$x, src.y = coords[src,]$y, src.z = coords[src,]$z, dest.x = coords[dest,]$x, dest.y = coords[dest,]$y, dest.z = coords[dest,]$z, count=counts[i])
      segments = rbind(segments, segment)
      seglabels = rbind(seglabels, data.frame(x = (segment$src.x + segment$dest.x)/2, y = (segment$src.y + segment$dest.y)/2,  z = (segment$src.z + segment$dest.z)/2, label=label))
    }
  }
  
  base = data.frame(src.x=0,src.z=0,dest.x=0,dest.z=0,count=0)
  
  message("Make plot")
  # plot
  if(use.width) {
    cube.plot = ggplot() +
      geom_segment(data=segments, aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x, size=count/2000), color="#AAAAAA") +
      scale_size_identity() +
      geom_point(data = coords[safe.nodes == T,], aes(x=z, y=x), size=p.size, color="#AAAAAA") +
      geom_text(data=seglabels, aes(x=z,y=x,label=label), color="#888888", size=lab.size) +
      xlim(-1,1) + ylim(-1,1) + theme_void() + theme(legend.position="none")
    if(node.labels) { cube.plot = cube.plot + geom_text_repel(data = coords[safe.nodes == T,], aes(x=z,y=x,label=label)) }
  } else {
    cube.plot = ggplot() +
      geom_segment(data=segments, aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x), color="#AAAAAA") +
      geom_point(data = coords[safe.nodes == T,], aes(x=z, y=x), size=p.size, color="#AAAAAA") +
      geom_text(data=seglabels, aes(x=z,y=x,label=label), color="#888888", size=lab.size) +
      xlim(-1,1) + ylim(-1,1) + theme_void() + theme(legend.position="none")
    if(node.labels) { cube.plot = cube.plot + geom_text_repel(data = coords[safe.nodes == T,], aes(x=z,y=x,label=label)) }
  }
  cube.plot
}

# plot graph of ordered pair acquisitions
plot.pfg = function(translist,              # list of transitions between states
                    pfg.layout = "matrix",  # graph loyout
                    curvature = 1           # geometric parameter for edge curviness
                    ) {
  
  message("Building PFG")
  ## PFG
  # get the set of pairwise first-next feature labels throughout each pathway
  edges=data.frame()
  message("Reading")
  zeroes = strsplit(translist[1], split=" ")[[1]][1]
  for(i in 1:1000) {
    #print(i)
    # get source and destination states
    src = strsplit(translist[i], " ")[[1]][1]
    dest = strsplit(translist[i], " ")[[1]][2]
    srcn = as.numeric(strsplit(src, "")[[1]])
    destn = as.numeric(strsplit(dest, "")[[1]])
    # identify changed feature
    change = which(srcn-destn != 0)
    # add this and last change to list of ordered pairs
    if(length(change) == 1) {
      if(src != zeroes) {
        edges = rbind(edges, data.frame(src=lastchange, dest=change))
      } else {
        edges = rbind(edges, data.frame(src=0, dest=change))
      }
      lastchange = change
    }
  }
  # get unique edges in this adjacency list and counts
  uedges = unique(edges)
  ucount = 0*uedges$src
  for(i in 1:nrow(uedges)) {
    ucount[i] = length(which(edges$src == uedges$src[i] & edges$dest == uedges$dest[i]))
  }
  
  message("Building adj mat")
  # construct graph from these edges
  uedges = rbind(uedges, uedges[1,])
  ucount[length(ucount)+1] = 1
  g = graph.data.frame(uedges)
  E(g)$weight = ucount
  maxw = max(ucount)
  sumw = sum(ucount)
  
  # sort nodes to a canonical order -- helps when we're comparing different cases, especially with the matrix layout
  s <- mixedsort(names(V(g)))
  new.g = igraph::permute(g, match(V(g)$name, s))
  g = new.g
  V(g)$name[1] = "-"
  
  message("Building plot")
  # plot PFG
  if(pfg.layout == "tree") {
    g.3 = ggraph(g, layout="tree") + geom_edge_bend(aes(edge_width=exp(weight/sumw), edge_alpha = weight/sumw), strength=curvature,  arrow=arrow()) + geom_node_point() + geom_node_label(aes(label=name), nudge_x = 0.05, nudge_y=-0.05) + theme_void() + theme(legend.position = "none")
  } else if(pfg.layout == "matrix") {
    g.3 = ggraph(g, layout="matrix") + geom_edge_bend(aes(edge_width=exp(weight/sumw), edge_alpha = weight/sumw), strength=curvature,  arrow=arrow()) + geom_node_point() + geom_node_label(aes(label=name), nudge_x = 0.05, nudge_y=-0.05) + theme_void() + theme(legend.position = "none")
  } else {
    g.3 = ggraph(g) + geom_edge_bend(aes(edge_width=exp(weight/sumw), edge_alpha = weight/sumw), strength=curvature,  arrow=arrow()) + geom_node_point() + geom_node_label(aes(label=name), nudge_x = 0.05, nudge_y=-0.05) + theme_void() + theme(legend.position = "none")
  }
  return(g.3)
}

# binary to decimal function
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

plot.hypercube2 = function(translist,               # set of transitions
                           use.width = T,           # use line width to display edge weights?
                           duplicate.offset = 0.,   # vertical offset for nodes in identical positions
                           lab.size = 3,            # size for edge labels
                           p.size = 1,              # point size
                           node.labels = T,         # node labels, yes or no?
                           seg.labels = T,          # line segment labels?
                           threshold = 0,           # ignore edges under a threshold in the hypercube plot
                           break.redundancy = 0,    # itself redundant now?
                           rotate.phi = F           # rotate states out of the page (in case of trajectories bunched up near the top/bottom)
                           ) {
  
  message("Building hypercube plot")  
  
  ## hypercube
  # get unique set of transitions and associated counts
  l = unique(translist)
  counts = rep(0, length(l))
  for(i in 1:length(l)) {
    set = which(translist == l[i])
    counts[i] = length(set)
  }
  #  l = l[counts > threshold]
  
  # split into lists of source and destination nodes
  srcs = dests = list()
  n = 1
  for(line in l) {
    s = strsplit(line, " ")
    srcs[[n]] = s[[1]][1]
    dests[[n]] = s[[1]][2]
    n = n + 1
  }
  # set string length and 0^L string
  len = nchar(srcs[[1]]) 
  zero = paste(rep("0", len), collapse="")
  
  # produce useful vectors
  srcs = unlist(srcs)
  dests = unlist(dests)
  
  all.nodes = apply(expand.grid(rep(list(0:1),len)), 1, paste, collapse="")
  nodes = all.nodes
  nnodes = length(nodes)
  
  # produce list storing where incoming edges to each node come from
  ins = list()
  for(node in nodes) {
    refs = which(dests == node)
    refcodes = srcs[refs]
    ins[[length(ins)+1]] = which(nodes %in% refcodes)
  }
  
  # produce hypercube visualisation
  
  message("Calculating embedding")
  # spherical polars: r, theta, phi
  # r = 1 everywhere
  rs = rep(1, nnodes)
  # theta is just set by number of 1s in a string
  
  thetas = unlist(lapply(nodes, function(s) { return(str_count(s, "0")*3.14159/len) }) )
  
  # initialise phis
  phis = rep(-1, nnodes)
  zero.counts = unlist(lapply(nodes, function(s) { return(str_count(s, "0")) }) )
  for(zeroes in 0:len) {
    refs = which(zero.counts == zeroes)
    these.phis = (0:(length(refs)-1))/length(refs)*3.14159
    phis[refs] = these.phis
  }
  
  # dataframes for spherical and cartesian coordinates
  spcoords = data.frame(r = rs, theta = thetas, phi = phis, label = nodes)
  # rotate phi values if required
  if(rotate.phi == T) {
    spcoords$phi = spcoords$phi + 3.14159/2
  }
  
  coords = data.frame(x = spcoords$r*cos(spcoords$phi)*sin(spcoords$theta), y = spcoords$r*sin(spcoords$phi)*sin(spcoords$theta), z = spcoords$r*cos(spcoords$theta), label = spcoords$label)
  
  # dataframe for line segments in cartesian coords
  segments = data.frame()
  seglabels = data.frame()
  safe.nodes = rep(F, nrow(coords))
  for(i in 1:length(srcs)) {
    if(counts[i] > threshold) {
      src = which(coords$label == srcs[i])
      dest = which(coords$label == dests[i])
      safe.nodes[src] = safe.nodes[dest] = T
      label = paste(c("+", which(unlist(strsplit(srcs[i], split="")) !=unlist(strsplit(dests[i], split="")))), collapse="")
      segment = data.frame(src.x = coords[src,]$x, src.y = coords[src,]$y, src.z = coords[src,]$z, dest.x = coords[dest,]$x, dest.y = coords[dest,]$y, dest.z = coords[dest,]$z, count=counts[i])
      segments = rbind(segments, segment)
      seglabels = rbind(seglabels, data.frame(x = (segment$src.x + segment$dest.x)/2, y = (segment$src.y + segment$dest.y)/2,  z = (segment$src.z + segment$dest.z)/2, label=label))
    }
  }
  
  base = data.frame(src.x=0,src.z=0,dest.x=0,dest.z=0,count=0)
  
  message("Make plot")
  # plot
  if(use.width) {
    cube.plot = ggplot() +
      geom_segment(data=segments, aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x, size=count/max(count)), color="#CCCCCC") +
      scale_size_identity() +
      geom_point(data = coords[safe.nodes == T,], aes(x=z, y=x), size=p.size, color="#444444") +
      xlim(-1,1) + ylim(-1,1) + theme_void() + theme(legend.position="none")
    if(node.labels == T) { cube.plot = cube.plot + geom_text_repel(data = coords[safe.nodes == T,], aes(x=z,y=x,label=label)) }
    if(seg.labels == T) {cube.plot = cube.plot + geom_text(data=seglabels, aes(x=z,y=x,label=label), color="#888888", size=lab.size) }
  } else {
    cube.plot = ggplot() +
      geom_segment(data=segments, aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x), color="#CCCCCC") +
      geom_point(data = coords[safe.nodes == T,], aes(x=z, y=x), size=p.size, color="#444444") +
      xlim(-1,1) + ylim(-1,1) + theme_void() + theme(legend.position="none")
    if(node.labels == T) { cube.plot = cube.plot + geom_text_repel(data = coords[safe.nodes == T,], aes(x=z,y=x,label=label)) }
    if(seg.labels == T) { cube.plot = cube.plot + geom_text(data=seglabels, aes(x=z,y=x,label=label), color="#888888", size=lab.size) }
    
  }
  return(cube.plot)
}

plot.standard = function(fitted) {
  plot.bubs = plot.bubbles2(fitted[[1]], formatted=TRUE)
  plot.cube = plot.hypercube2(fitted[[4]], use.width = T, node.labels=F, seg.labels = T, threshold=0, rotate.phi=F)
  plot.diag = plot.pfg(fitted[[4]], pfg.layout="matrix")
  # and arrange plots together
  return(ggarrange(plot.bubs, plot.cube, plot.diag, nrow=1))
}

DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

BinToDec <- function(state) {
  this.ref = 0
  for(j in 1:length(state)) {
    this.ref = this.ref + state[j]*(2**(length(state)-j))
  }
  return(this.ref)
}

# adapted HyperTraPS plot for HyperHMM outputs
plot.hypercube.flux = function(my.post, thresh = 0.05, node.labels = TRUE, use.probability = FALSE) {
  ### produce hypercube subgraph
  bigL = my.post$L
  if(use.probability == TRUE) {
    trans.p = my.post$transitions[my.post$transitions$Probability > thresh & my.post$transitions$Bootstrap == 0,]
  } else {
    trans.p = my.post$transitions[my.post$transitions$Flux > thresh & my.post$transitions$Bootstrap == 0,]
  }
  trans.g = graph_from_data_frame(trans.p[,2:ncol(trans.p)])
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = V(trans.g)$name
  layers = str_count(bs, "1")
  if(use.probability == TRUE) {
    this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
      geom_edge_link(aes(edge_width=Probability, edge_alpha=Probability)) +
      scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,0.5)) + theme_graph() #+
  } else {
  this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
    geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux)) +
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,0.5)) + theme_graph() #+
   }
  if(node.labels == TRUE) {
    this.plot = this.plot + geom_node_text(aes(label=binname),angle=45,size=2) 
  }
  return(this.plot)
}

# load required libraries
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(gridExtra)
library(igraph)

for(expt in 1:1) {
  use.width = T           # use line width to display edge weights?
  duplicate.offset = 0.   # vertical offset for nodes in identical positions
  lab.size = 3            # size for edge labels
  p.size = 1              # point size
  node.labels = T         # node labels, yes or no?
  threshold = 0           # ignore edges under a threshold in the hypercube plot
  curvature = 1           # curvature of edges in the PFG plot

  # read data from different experiments
  if(expt == 1) {
    tl = readLines("graph_viz_simple1.txt.txt"); p.size = 2
  }
  
  # get unique set of transitions and associated counts
  l = unique(tl)
  counts = rep(0, length(l))
  for(i in 1:length(l)) {
    set = which(tl == l[i])
    counts[i] = length(set)
  }
    
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
  
  
  # dataframes for spherical and cartesian coordinates
  spcoords = data.frame(r = rs, theta = thetas, phi = phis, label = nodes)
  
  transcoords = data.frame()
  for(this.theta in unique(spcoords$theta)) {
    subset = spcoords[spcoords$theta == this.theta,]
    this.rank = rank(subset$phi, ties.method="random")
    subset$phi = 3.14159*(this.rank-1)/(max(this.rank)-1)
    if(nrow(subset) == 1) { subset$phi = 0 }
    transcoords = rbind(transcoords, subset)
  }
  
  spcoords = transcoords
  
  coords = data.frame(x = spcoords$r*cos(spcoords$phi)*sin(spcoords$theta), y = spcoords$r*sin(spcoords$phi)*sin(spcoords$theta), z = spcoords$r*cos(spcoords$theta), label = spcoords$label)
  
  
  # dataframe for line segments in cartesian coords
  segments = data.frame()
  seglabels = data.frame()
  for(i in 1:length(srcs)) {
    if(counts[i] > threshold) {
      src = which(coords$label == srcs[i])
      dest = which(coords$label == dests[i])
      label = paste(c("+", which(unlist(strsplit(srcs[i], split="")) !=unlist(strsplit(dests[i], split="")))), collapse="")
      segment = data.frame(src.x = coords[src,]$x, src.y = coords[src,]$y, src.z = coords[src,]$z, dest.x = coords[dest,]$x, dest.y = coords[dest,]$y, dest.z = coords[dest,]$z, count=counts[i])
      segments = rbind(segments, segment)
      seglabels = rbind(seglabels, data.frame(x = (segment$src.x + segment$dest.x)/2, y = (segment$src.y + segment$dest.y)/2,  z = (segment$src.z + segment$dest.z)/2, label=label))
    }
  }
  
  base = data.frame(src.x=0,src.z=0,dest.x=0,dest.z=0,count=0)
  
  # plot
  if(use.width) {
    cube.plot = ggplot() +
      geom_segment(data=segments, aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x, size=count/5000), color="#AAAAAA") +
      scale_size_identity() +
      geom_point(data = coords, aes(x=z, y=x), size=p.size, color="#AAAAAA") +
      geom_text(data=seglabels, aes(x=z,y=x,label=label), color="#888888", size=lab.size) +
      xlim(-1,1) + ylim(-1,1) + theme_void() + theme(legend.position="none")
    if(node.labels) { cube.plot = cube.plot + geom_text_repel(data = coords, aes(x=z,y=x,label=label)) }
  } else {
    cube.plot = ggplot() +
      geom_segment(data=segments, aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x), color="#AAAAAA") +
      geom_point(data = coords, aes(x=z, y=x), size=p.size, color="#AAAAAA") +
      geom_text(data=seglabels, aes(x=z,y=x,label=label), color="#888888", size=lab.size) +
      xlim(-1,1) + ylim(-1,1) + theme_void() + theme(legend.position="none")
    if(node.labels) { cube.plot = cube.plot + geom_text_repel(data = coords, aes(x=z,y=x,label=label)) }
  }
  g1 = cube.plot
  
  
  ########### now produce probabilistic feature graph

  # get the set of pairwise first-next feature labels throughout each pathway
  edges=data.frame()
  zeroes = strsplit(tl[1], split=" ")[[1]][1]
  for(i in 1:1000) {
    src = strsplit(tl[i], " ")[[1]][1]
    dest = strsplit(tl[i], " ")[[1]][2]
    srcn = as.numeric(strsplit(src, "")[[1]])
    destn = as.numeric(strsplit(dest, "")[[1]])
    change = which(srcn-destn != 0)
    if(src != zeroes) {
      edges = rbind(edges, data.frame(src=lastchange, dest=change))
    } else {
      edges = rbind(edges, data.frame(src=0, dest=change))
    }
    lastchange = change
  }
  # get unique edges and counts
  uedges = unique(edges)
  ucount = 0*uedges$src
  for(i in 1:nrow(uedges)) {
    ucount[i] = length(which(edges$src == uedges$src[i] & edges$dest == uedges$dest[i]))
  }

  # construct graph from these edges
  uedges = rbind(uedges, uedges[1,])
  ucount[length(ucount)+1] = 1
  g = graph.data.frame(uedges)
  E(g)$weight = ucount
  maxw = max(ucount)

  # plot PFG
  g2 = ggraph(g) + geom_edge_bend(aes(edge_width=exp(weight/maxw), edge_alpha = weight/maxw), strength=curvature,  arrow=arrow()) + geom_node_point() + geom_node_label(aes(label=name)) + theme_void() + theme(legend.position = "none")

  if(expt >= 4) {
    png(paste(c("output-", expt, ".png"), collapse=""), width=2400, height=800, res=150)
  } else {
    png(paste(c("output-", expt, ".png"), collapse=""), width=1600, height=500, res=150)
  }
  grid.arrange(g1, g2, nrow=1)
  dev.off()
}


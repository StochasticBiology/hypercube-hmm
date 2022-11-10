### Wrapper script for external calls to HyperHMM

library(stringr)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(gridExtra)
library(igraph)

## try to find (or compile) the HyperHMM executable in this directory for external calling

# possible names of executable
cmds = c("hyperhmm.ce", "hyperhmm.exe", "hyperhmm")

# see if any are here
found = F
for(cmd in cmds) {
  if(cmd %in% list.files()) { commandname = cmd; found = T }
}
if(found == F) {
  message("Couldn't find HyperHMM executable in this folder. Attempting an automatic compilation.")
  # detect platform -- try g++ on unix and cl (Visual Studio) on Windows
  if(.Platform$OS.type == "unix") {
#    system("g++ hyperhmm.cpp -lm -larmadillo -o hyperhmm.ce")
    system("g++ hyperhmm.cpp -I/opt/homebrew/Cellar/armadillo/11.4.2/include -L/opt/homebrew/Cellar/armadillo/11.4.2/lib -larmadillo -std=c++11 -o hyperhmm.ce")
  } else {
    system("cl /EHsc hyperhmm.cpp")
  }
  # after attempted compilation, check for executable again
  for(cmd in cmds) {
    if(cmd %in% list.files()) { commandname = cmd; found = T }
  }
  if(found == F) {}
  message("Couldn't find HyperHMM executable in this folder. Smart automatic compilation isn't yet supported; please compile the source code (using, for example, gcc on Mac/Linux or Visual Studio on Windows) and run this script from the folder containing the executable.")
  stop()
}

## wrapper function for a HyperHMM call
hyperhmm = function(obs, precursors=NA, nboot=1, random.walkers=0, label="label", simulate = T, fork=F) {
  
  # catch various issues
  if(!is.matrix(obs)) { message("Couldn't interpret observations as a matrix!"); return() }
  if(ncol(obs) < 2) { message("Inference with < 2 features is meaningless!"); return() }
  if(ncol(obs) > 15) { message("You're welcome to try > 15 features but this may use a lot of (too much?) memory."); }
  if(nrow(obs) == 0) { message("No data in observation matrix!"); return() }

  # get number and names of features
  L = ncol(obs)
  features = colnames(obs)
  
  # deal with non-binary values
  if(any(obs != 0 & obs != 1)) {
    message("Non-binary observations found. Interpreting all nonzero entries as presence.")
    obs[obs != 0 & obs != 1] = 1
  }
  
  # compile into strings
  obs.sums = rowSums(obs)
  obs.rows = apply(obs, 1, paste, collapse="")
  final.rows = obs.rows
  
  # check to see if precursor states have been provided; set cross-sectional flag accordingly
  if(!is.matrix(precursors)) {
    message("No precursor observations; assuming cross-sectional setup.")
    cross.sectional = 1
  } else {
    cross.sectional = 0
    
    # deal with issues in precursor set
    if(!is.matrix(precursors)) { message("Couldn't interpret precursor observations as a matrix!"); return() }
    if(nrow(precursors) != nrow(obs)) { message("Number of observations and precursors didn't match!"); return() }
    if(ncol(precursors) != ncol(obs)) { message("Number of features in observations and precursors didn't match!"); return() }
    
    if(any(precursors != 0 & precursors != 1)) {
      message("Non-binary precursors found. Interpreting all nonzero entries as presence.")
      precursors[precursors != 0 & precursors != 1] = 1
    }

    # check to make sure all precursor states come before their descendants
    precursor.sums = rowSums(precursors)
    if(any(precursor.sums > obs.sums)) {
      message("Precursors found more advanced than observations!")
      return() 
    }
    
    # compile precursor-observation pairs into strings
    precursor.rows = apply(precursors, 1, paste, collapse="")
    for(i in 1:length(obs.rows)) {
      final.rows[i] = paste(c(precursor.rows[i], obs.rows[i]), collapse=" ")
    }
  }
  if(any(obs.sums == 0)) {
    message("Dropping 0^L observations")
  }
  final.rows = final.rows[obs.sums != 0]
  
  # create input file for HyperHMM
  filename = paste(c("hhmm-in-", label, ".txt"), collapse="")
  write(final.rows, filename)
  
  # create call to HyperHMM
  if(simulate == T) {
    if(fork == T) {
  syscommand = paste(c("./", commandname, " ", filename, " ", L, " ", nboot, " ", label, " ", cross.sectional, " ", random.walkers, " &"), collapse="")
  message(paste(c("Attempting to externally execute ", syscommand), collapse=""))
  system(syscommand)
  message("Forked: not retrieving output")
  return(NULL)
    } else {
      syscommand = paste(c("./", commandname, " ", filename, " ", L, " ", nboot, " ", label, " ", cross.sectional, " ", random.walkers), collapse="")
      message(paste(c("Attempting to externally execute ", syscommand), collapse=""))
      system(syscommand)
    }
  }

  # attempt to read the output of the run
  
  # BUG FIXED NEEDS PUSH: seems there's an underscore in the name for longitudinal but not cross-sectional?
  mean.filename = paste(c("mean_", label, ".txt"), collapse="")
  sd.filename = paste(c("sd_", label, ".txt"), collapse="")
  transitions.filename = paste(c("transitions_", label, ".txt"), collapse="")
  viz.filename = paste(c("graph_viz_", label, ".txt"), collapse="")
  
  message(paste(c("Attempting to read output... e.g. ", mean.filename), collapse=""))
  if(!(mean.filename %in% list.files()) | !(sd.filename %in% list.files()) | !(transitions.filename %in% list.files()) | !(viz.filename %in% list.files())) {
    message(paste(c("Couldn't find file output of externally executed HyperHMM!"), collapse=""))
    return()
  }
  mean.table = read.table(mean.filename)
  sd.table = read.table(sd.filename)
  transitions = read.csv(transitions.filename, sep=" ")
  # BIG FIXED NEEDS PUSH -- double .txt on viz name
  viz.tl = readLines(viz.filename);
  
  message("All expected files found.")
  L = nrow(mean.table)
  # pull the wide output data into long format
  stats.df = data.frame(feature=rep(0,L*L), order=rep(0,L*L), mean=rep(0, L*L), sd=rep(0, L*L))
  index = 1
  for(i in 1:L) {
    for(j in 1:L) {
      stats.df$feature[index] = L-i+1
      stats.df$order[index] = j
      stats.df$mean[index] = mean.table[i,j]
      stats.df$sd[index] = sd.table[i,j]
      index = index+1 
    }
  }
  
  # create a list of useful objects and return it
  fitted = list(stats.df, transitions, features, viz.tl)
  return(fitted)
}

plot.hyperhmm = function(fitted, use.width = T,           # use line width to display edge weights?
                         duplicate.offset = 0.,   # vertical offset for nodes in identical positions
                         lab.size = 3,            # size for edge labels
                         p.size = 1,              # point size
                         node.labels = T,         # node labels, yes or no?
                         threshold = 0,           # ignore edges under a threshold in the hypercube plot
                         break.redundancy = 0,
                         tree.layout = T,
                         curvature = 1) {

  message("Building bubble plot")
  ## bubbles
  if(is.null(fitted[[3]])) {
    g.1 = ggplot(fitted[[1]], aes(x=order, y=feature)) + geom_point(aes(size=mean), colour="#AAAAAA") + geom_point(aes(size=sd), shape = 1, colour="#888888") + theme_classic() 
  } else {
    g.1 = ggplot(fitted[[1]], aes(x=order, y=feature)) + geom_point(aes(size=mean), colour="#AAAAAA") + geom_point(aes(size=sd), shape = 1, colour="#888888") + scale_y_continuous(breaks=length(fitted[[3]]):1, labels=fitted[[3]]) + theme_classic() 
  }
  message("Building hypercube plot")  
  ## hypercube
  # get unique set of transitions and associated counts
  l = unique(fitted[[4]])
  counts = rep(0, length(l))
  for(i in 1:length(l)) {
    set = which(fitted[[4]] == l[i])
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
  g.2 = cube.plot
  
  message("Building PFG")
  ## PFG
  # get the set of pairwise first-next feature labels throughout each pathway
  edges=data.frame()
  message("Reading")
  zeroes = strsplit(fitted[[4]][1], split=" ")[[1]][1]
  for(i in 1:1000) {
    src = strsplit(fitted[[4]][i], " ")[[1]][1]
    dest = strsplit(fitted[[4]][i], " ")[[1]][2]
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
  
  message("Building adj mat")
  # construct graph from these edges
  uedges = rbind(uedges, uedges[1,])
  ucount[length(ucount)+1] = 1
  g = graph.data.frame(uedges)
  E(g)$weight = ucount
  maxw = max(ucount)
  
  message("Building plot")
  # plot PFG
  if(tree.layout) {
    g.3 = ggraph(g, layout="tree") + geom_edge_bend(aes(edge_width=exp(weight/maxw), edge_alpha = weight/maxw), strength=curvature,  arrow=arrow()) + geom_node_point() + geom_node_label(aes(label=name), nudge_x = 0.05, nudge_y=-0.05) + theme_void() + theme(legend.position = "none")
  } else {
    g.3 = ggraph(g) + geom_edge_bend(aes(edge_width=exp(weight/maxw), edge_alpha = weight/maxw), strength=curvature,  arrow=arrow()) + geom_node_point() + geom_node_label(aes(label=name), nudge_x = 0.05, nudge_y=-0.05) + theme_void() + theme(legend.position = "none")
  }
  
  message("Laying out")
  grid.arrange(g.1, g.2, g.3, nrow=1)
}


                                                                                  
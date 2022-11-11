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
    system("g++ hyperhmm.cpp -lm -larmadillo -o hyperhmm.ce")
    # ^ that's a sort of default compile command on a Unix box, but others will differ. My Mac for example needs
    # system("g++ hyperhmm.cpp -I/opt/homebrew/Cellar/armadillo/11.4.2/include -L/opt/homebrew/Cellar/armadillo/11.4.2/lib -larmadillo -std=c++11 -o hyperhmm.ce")
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
hyperhmm = function(obs,                 # matrix of observations
                    precursors=NA,       # precursors (e.g. ancestors) -- blank for cross-sectional
                    nboot=1,             # number of bootstrap resamples
                    random.walkers=0,    # run random walkers for each resample? 0 no 1 yes
                    label="label",       # label for file I/O 
                    simulate = T,        # actually run HyperHMM? If not, try and pull precomputed output using "label"
                    fork=F               # if we're running it, fork to a new process, or keep linear?
                    ) {
  
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


## wrapper function for a HyperTraPS call
hypertraps = function(obs,                   # matrix of observations
                      precursors=NA,         # matrix of precursors (e.g. for ancestor-descendant samples)
                      param.length = 2,      # MCMC chain length: 10^(n+2)
                      param.kernel = 5,      # MCMC perturbation kernel -- see source code
                      label="label",         # label for file I/O
                      simulate="fork"        # "linear" to call in code; "fork" to fork; other to not simulate and just retrieve info
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
    for(i in 1:length(obs.rows)) {
      final.rows[2*(i-1)+1] = paste(rep(0, L), collapse="")
      final.rows[2*(i-1)+2] = paste(obs.rows[i], collapse="")
    }
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
      final.rows[2*(i-1)+1] = paste(precursor.rows[i], collapse="")
      final.rows[2*(i-1)+2] = paste(obs.rows[i], collapse="")
    }
  }
  
  # create input file for HyperTraPS
  filename = paste(c("ht-in-", label, ".txt"), collapse="")
  write(final.rows, filename)
  
  # create call to HyperTraPS
  if(simulate == "linear") {
    syscommand = paste(c("./hypertraps-dt.ce ", filename, " 1 ", param.length, " ", param.kernel, " 0"), collapse="")
    message(paste(c("Attempting to externally execute ", syscommand), collapse=""))
    system(syscommand)
  } else if(simulate == "fork") {
    syscommand = paste(c("./hypertraps-dt.ce ", filename, " 1 ", param.length, " ", param.kernel, " 0 &"), collapse="")
    message(paste(c("Attempting to externally execute ", syscommand), collapse=""))
    system(syscommand)
    message("Not retrieving info -- process forked")
    return(NULL)
  }
  
  # attempt to read the output of the run
  message("Attempting to read output...")
  
  # construct expected posterior filename
  post.filename = paste(c(filename, "-posterior-0-1-", param.length, "-", param.kernel, ".txt"), collapse="")
  
  message(post.filename)
  if(!(post.filename %in% list.files())) {
    message("Couldn't find file output of externally executed HyperTraPS!")
    return()
  }
  post.table = read.table(post.filename)
  
  message("All expected files found.")
  return(post.table)
}

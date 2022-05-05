# hypercube-hmm
Hypercubic inference using Hidden Markov Models

Code to infer the transition probabilities on a hypercubic transition network given some observations of emitted signals from a hidden Markov model on that network. Visualisations of the inferred parameterised model and its summary dynamics are also performed. 

The inference code uses the C++ Armadillo library \cite{sanderson2016armadillo}. Bubble plot visualisations use Python libraries matplotlib \cite{Hunter:2007}, pandas \cite{mckinney-proc-scipy-2010}, seaborn \cite{Waskom2021}, and numpy \cite{harris2020array}. Other visualisations use R libraries stringr \cite{stringr}, ggplot2 \cite{ggplot2}, ggrepel \cite{ggrepel}, gridExtra \cite{gridExtra}, and igraph \cite{igraph}.

`run-example.sh` is a Bash script wrapping a simple, fast example case. `run.sh` wraps the full set of HyperHMM experiments.

`hyperhmm.cpp` is the inference code. Compile with, for example,

`g++ hyperhmm.cpp -o hyperhmm.ce`

The code takes several command line arguments:

`./hyperhmm.ce [datafile] [number of features] [number of bootstrap resamples] [output file label] [cross-sectional data (0 or 1)] [simulate random walkers for each sample (0 or 1)]`

- `datafile` -- a datafile containing cross-sectional or longitudinal observations as bit strings (see below)
- `number of features` -- the number of features/traits/characters in the system (i.e. the length of a bit string in the datafile)
- `number of bootstrap resamples` -- 0 for a single point estimate; >0 runs the given number of resamples in addition
- `output file label` -- a string to label output files
- `cross-sectional data` -- if 1, data are assumed to be independent snapshots and the datafile should have one column. If 0, data are assumed to be observed transitions between two states and the datafile should have two columns (ancestor and descendant state)
- `simulate random walkers for each sample` -- random walkers are simulated on the point estimate hypercube to summarise the dynamics. If 1, they are simulated on every bootstrapped hypercube too and the summary is over all resamples. This can make the code take much longer for simple systems.

For example,
`./hyperhmm.ce Data/simple_case1_L5.txt 5 100 simple1 1 0`

`bubble_plots.py` summarises the outputs of the various experiments using `hyperhmm.cpp` as bubble plots.

`cube-pfg.R` produces hypercube visualisations and probabilistic feature graphs from the outputs.

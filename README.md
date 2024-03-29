# hypercube-hmm
Hypercubic inference using Hidden Markov Models, from https://academic.oup.com/bioinformatics/article/39/1/btac803/6895098 

![image](https://github.com/StochasticBiology/hypercube-hmm/assets/50171196/893b96c8-12ec-4511-b314-567955870ab2)

Code for HyperHMM [1] to infer the transition probabilities on a hypercubic transition network given some observations of emitted signals from a hidden Markov model on that network. Visualisations of the inferred parameterised model and its summary dynamics are also performed. This is the expectation-maximisation cousin of HyperTraPS https://github.com/StochasticBiology/HyperTraPS (simple implementation here https://github.com/StochasticBiology/hypertraps-simple).

We also compare the behaviour of HyperHMM with HyperTraPS [2,3], Oncotrees [4], conjunctive Bayes networks (CBN) via TRONCO [5], and mutual hazard networks (MHN) [6].

Requirements
-------

The inference code uses the C++ `armadillo` library [7]. Other visualisations use R libraries `stringr` [8], `ggplot2` [9], `ggrepel` [10], `gridExtra` [11], and `igraph` [12].

The R implementation requires `Rcpp` [13] and `RcppArmadillo` [14].

Contents
=======

Examples
-------

_In R:_

Make sure you have the `Rcpp` and `RcppArmadillo` libraries installed. Then you can access HyperHMM functionality with

```
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")
```

The function `HyperHMM`, described further below, performs inference given (at least) a matrix of observations. For example,

```
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
HyperHMM(m)
```

Have a look at `hyperhmm-demos.R` for some simple examples, and `hyperhmm-run.R` for more examples from the paper.

_From the command line:_

You'll need to compile the C++ code. Install the `armadillo` library on your machine, then the command to compile from the Terminal may look something like
```
g++ hyperhmm.cpp -larmadillo -o hyperhmm.ce
```

Then a simple example could be run with

```
./hyperhmm.ce --obs Data/simple_case1_L5.txt
```

Inference
-------

The code takes several command-line arguments:

| Argument | R | Command-line | Default |
|----------|---|--------------|---------|
| Input data | obs=*matrix* | --obs *filename* | None (required) |
| Label for output files | [none] | --label *labelstring* | *filename*-out |
| Random number seed (for bootstrapping) | seed=*N* | --seed *N* | 1 |
| Number of bootstrap resamples | nboot=*N*| --nboot *N* | 100 |
| Simulate random walkers for each resample? | fullsample=*N* | --fullsample | (off) |
| Initial states for longitudinal data | initialstates=*matrix* | [in input data, see below] | (none)
| Longitudinal data? | [implied by initialstates] | --longitudinal | (off) |

The code outputs several objects. In R, these are elements of a named list that is returned by `HyperHMM`; from the command line, they are datafiles. These objects describe the maximum likelihood inferred transition probabilities between states (`$transitions` in R; `transitions_labelstring.txt` from the command-line), the mean and bootstrap standard deviation (`$stats`; wide-form `mean_labelstring.txt` and `sd_labelstring.txt`), and a set of transitions from sampled trajectories of random walkers simulated on the inferred network (`$viz`; `graph_viz_labelstring.txt`).

Data
------
Synthetic and published data is in `Data`. The ovarian cancer dataset is from [15]; the tuberculosis dataset is from [16].

_In R:_

Data is provided as a matrix via the required first argument "obs", where each row gives an independent snapshot observation of length L, for example

`0011`<br>
`0111`<br>
`1101`

For longitudinal data, a matrix with the same number of rows and columns can be provided describing the initial states for each observation via the argument "initialstates", for example

`0001`<br>
`0011`<br>
`1001`

Have a look at `hyperhmm-demos.R` for some examples, and `hyperhmm-run.R` for more examples from the paper.

_For the command line:_

For cross-sectional observations, data should be provided as a single-column file where each row gives an independent snapshot observation of length L, for example

`0011`<br>
`0111`<br>
`1101`

For longitudinal observations, including those derived from estimated phylogenies, data should be provided as a two-column file where each row gives an independent transition observation between two states of length L, for example

`0001 0011`<br>
`0100 0111`<br>
`1001 1101`

The initial state is assumed to be all `0`s, and the system evolves by acquiring `1`s.

Examples
------
_In R:_

| Command | Description |
|---------|---------|
| `HyperHMM(m)` | Cross-sectional observations m |
| `HyperHMM(m.2, initialstates=m.1)` | Longitudinal observations m.1 -> m.2 |
| `HyperHMM(m, nboot=1000, seed=2)` | As top row with 1000 resamples and random number seed 2 |

Output and plotting
------

*In R*, `HyperHMM` returns a named list with the following elements

| Element | Description |
|---------|---------|
| *stats* | "Bubble plot" summary statistics: the mean and bootstrap s.d. of the probability that a feature is acquired at a given ordering step in the accumulation process |
| *transitions* | The weights and probability fluxes for each transition between states, for each bootstrap resample and for the original dataset (corresponding to *Bootstrap* == 0 in the data frame) |
| *features* | Empty by default, to allow for feature names |
| *viz* | Collection of observed transitions from simulated random walkers, represented as string containing before and after states, for use in sampled visualisation |
| *L* | The number of features in the system 

Plotting is done in R, using functions in `hypercube-plots.R`. These include "bubble" plots for mean feature orderings, hypercube visualisations illustrating the complete transition network, and ordering graphs showing pairwise acquisition orderings.

These include:

| Element | Description | Optional arguments
|---------|---------|------|
| `plot.bubbles2` | "Bubble plot" summary statistics: the mean and bootstrap s.d. of the probability that a feature is acquired at a given ordering step in the accumulation process | |
| `plot.hypercube2` | Visualisation of the hypercubic transition graph | |
| `plot.pfg` | Visualisation of which feature acquisitions follow which others | |
| `plot.standard` | Collection of the three above plot types | |
| `plot.hypercube.flux` | (Rather improved) visualisation of the hypercubic transition graph -- recommended over the above |

Other approaches and manuscript links
-------

From earlier development, the code base includes several R scripts that "wrap" external calls to HyperHMM, including preparing data and pulling output for plots. This setup was less streamlined than fully integrating HyperHMM into the R environment via Rcpp. There is also code that "wraps" HyperTraPS, and explores other inference approaches. These are now all in the `wrappers` directory and will currently have some path issues relative to the other files.

The code `hyperhmm-wrap.R` provides a function to call HyperHMM from R; `hypertraps-wrap.R` does the same for HyperTraPS.

The analyses and figures in the associated manuscript are reproduced with the various `...-start.R` and `...-retrieve.R` scripts. If you want to run these, pull the `Data/` contents into the working directory first. *As of a 2023 update, the HyperHMM results in the paper are also produced with a single R script `hyperhmm-run.R` in the root directory.*

Specifically, `double-approach-....R` compares HyperHMM and HyperTraPS, and `other-approaches-....R` compares HyperHMM, HyperTraPS, Oncotrees, MHN, and CBN (TRONCO).

Old command-line formatting
------

In its first published form, the following was true. This input format is preserved for backwards compatibility but the setup above is recommended for flexibility.

The code takes several command line arguments:

`./hyperhmm.ce [datafile] [number of features] [number of bootstrap resamples] [output file label] [cross-sectional data (0 or 1)] [simulate random walkers for each sample (0 or 1)]`

- `datafile` -- a datafile containing cross-sectional or longitudinal observations as bit strings (see below)
- `number of features` -- the number of features/traits/characters in the system (i.e. the length of a bit string in the datafile)
- `number of bootstrap resamples` -- 0 for a single point estimate; >0 runs the given number of resamples in addition
- `output file label` -- a string to label output files
- `cross-sectional data` -- if 1, data are assumed to be independent snapshots and the datafile should have one column. If 0, data are assumed to be observed transitions between two states and the datafile should have two columns (ancestor and descendant state)
- `simulate random walkers for each sample` -- random walkers are simulated on the point estimate hypercube to summarise the dynamics. If 1, they are simulated on every bootstrapped hypercube too and the summary is over all resamples. This can make the code take much longer for simple systems.

For example,
`./hyperhmm.ce Data/simple_case1_L5.txt 5 100 simple1 1 1`

References
=====

[1] Moen, M.T. and Johnston, I.G., 2023. HyperHMM: efficient inference of evolutionary and progressive dynamics on hypercubic transition graphs. Bioinformatics, 39(1), p.btac803.

[2] Johnston, I.G. and Williams, B.P., 2016. Evolutionary inference across eukaryotes identifies specific pressures favoring mitochondrial gene retention. Cell systems, 2(2), pp.101-111.

[3] Greenbury, S.F., Barahona, M. and Johnston, I.G., 2020. HyperTraPS: inferring probabilistic patterns of trait acquisition in evolutionary and disease progression pathways. Cell systems, 10(1), pp.39-51.

[4] Szabo,A. and Pappas,L. (2022) Oncotree: Estimating Oncogenetic Trees. R package version 0.3.4.

[5] De Sano,L. et al. (2016) TRONCO: an R package for the inference of cancer progression models from heterogeneous genomic data. Bioinformatics, 32, 1911–1913.

[6] Schill,R. et al. (2020) Modelling cancer progression using mutual hazard networks. Bioinformatics, 36, 241–249.

[7] Conrad Sanderson and Ryan Curtin. Armadillo: a template-based c++ library for linear algebra. Journal of Open Source Software, 1(2):26, 2016.

[8] Hadley Wickham. stringr: Simple, Consistent Wrappers for Common String Operations, 2019. R package version 1.4.0.

[9] Hadley Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

[10] Kamil Slowikowski. ggrepel: Automatically Position Non-Overlapping Text Labels with ’ggplot2’, 2021. R package version 0.9.1.

[11] Baptiste Auguie. gridExtra: Miscellaneous Functions for ”Grid” Graphics, 2017. R package version 2.3.

[12] Gabor Csardi and Tamas Nepusz. The igraph software package for complex network research. InterJournal, Complex Systems:1695, 2006.

[13] Dirk Eddelbuettel and Romain Francois (2011). Rcpp:Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1-18, <doi:10.18637/jss.v040.i08>.

[14] Eddelbuettel D, Sanderson C (2014). “RcppArmadillo: Accelerating R with high-performance C++ linear algebra.” _Computational Statistics and Data Analysis_, *71*, 1054-1063. doi:10.1016/j.csda.2013.02.005 <https://doi.org/10.1016/j.csda.2013.02.005>.

[15] Turid Knutsen et al. The interactive online sky/m-fish & cgh database and the entrez cancer chromosomes search database: linkage of chromosomal aberrations with the genome sequence. Genes, Chromosomes and Cancer, 44(1):52–64, 2005.

[16] Nicola Casali et al. Evolution and transmission of drug-resistant tuberculosis in a russian population. Nature genetics, 46(3):279–286, 2014.

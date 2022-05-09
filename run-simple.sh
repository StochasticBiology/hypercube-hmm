# compile the code
g++ hyperhmm.cpp -o hyperhmm.ce

# do inference on simple synthetic case; 10 bootstrap resamples with random walker sampling for each
./hyperhmm.ce Data/simple_case1_L5.txt 5 10 simple1 1 1

# plot "bubble plot" summaries of this output
python3 bubble_plots_simple.py

# plot hypercubic and probabilistic feature graphs
Rscript cube-pfg-simple.R

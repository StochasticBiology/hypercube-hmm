# wrapper function doing inference and producing summary plots
processcase() {
  ./hyperhmm.ce $1 $2 $3 $4 $5 $6 
  python3 bubble_plots_arg.py $4
  Rscript cube-pfg-arg.R $4
}

# individual experiments from the manuscript
# third argument: 10 bootstrap resamples run for speed. use 100 (or more) for precision
processcase Data/simple_case1_L5.txt 5 10 simple1-5 1 1 > tmp1 &
processcase Data/simple_case2_L5.txt 5 10 simple2-5 1 1 > tmp2 &
processcase Data/simple_case2_L7.txt 7 10 simple2-7 1 1 > tmp3 &
processcase Data/simple_case2_L9.txt 9 10 simple2-9 1 1 > tmp4 &
processcase Data/simple_case4_L5.txt 5 10 simple4-5 1 1 > tmp5 &
processcase Data/simple_case10_L5.txt 5 10 simple10-5 1 1 > tmp6 &
processcase Data/double_case2_L5.txt 5 10 double2-5 1 1 > tmp7 &
processcase Data/double_case2_L7.txt 7 10 double2-7 1 1 > tmp8 &
processcase Data/double_case2_L9.txt 9 10 double2-9 1 1 > tmp9 &
processcase Data/ovarian.txt 7 10 ovarian 1 1 > tmp10 &
processcase Data/tb_drug.txt 10 10 tb_drug 0 1 > tmp11 &



library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

m.1 = matrix(c(0,0,0,0,0,
                   1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   0,0,0,0,0,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1), byrow=TRUE, ncol = 5)

m.2 = matrix(c(1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   1,1,1,1,1,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1,
                   1,1,1,1,1), byrow=TRUE, ncol = 5)
             
HyperHMM(m.2)
HyperHMM(m.2, initialstates = m.1, outputinput = 1)

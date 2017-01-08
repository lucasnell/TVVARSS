TVVARSS
========

Time-varying vector autoregressive state-space model
--------

Initial code by Anthony Ives

Implementation in C++ by Lucas Nell

The main files in this repository are the following:

- `TVVARSS.cpp`: C++ version of the `TVVARSS.ml` function
- `TVVARSS_cpp.R`: the `cpp_TVVARSS` function that uses the C++ version of `TVVARSS.ml`
- `TVVARSS_sim_test_22Dec16.R`: very preliminary tests of `cpp_TVVARSS`

The C++ implementation used the "Armadillo" C++ library and interfaced with R using
the packages `Rcpp` and `RcppArmadillo`.


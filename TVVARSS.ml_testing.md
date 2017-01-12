Testing `TVVARSS.ml`
================
Lucas Nell
2017-01-12

This document outlines the results of using optimization and simulated annealing on R and C++ versions of `TVVARSS.ml`. The code to conduct these tests is in [this file](./initial_testing_files/TVVARSS.ml_testing.R).

Optimization
============

### Using `optim(..., method = 'Nelder-Mead')`

    R version took 1 minutes
    Rcpp version took 0.01 minutes

    Are they equal?
    Component "par": Mean relative difference: 0.244417
    Component "value": Mean relative difference: 0.04055309
    Component "counts": Mean relative difference: 0.886023

### Using `optim(..., method = 'BFGS')`

    R version took 2.92 minutes
    Rcpp version took 0.02 minutes

    Are they equal?
    Component "par": Mean relative difference: 2.608607
    Component "value": Mean relative difference: 2.524436
    Component "counts": Mean relative difference: 0.872134

### Using `bobyqa(...)`

    R version took 0.45 minutes
    Rcpp version took 0.01 minutes

    Are they equal?
    Component "par": Mean relative difference: 0.0001922985
    Component "fval": Mean relative difference: 0.002921203
    Component "feval": Mean relative difference: 0.6392727

Simulated annealing
===================

    R version took 8.98 minutes
    Rcpp version took 0.58 minutes

    Are they equal?
    Component "trace.mat": Mean relative difference: 2.694283e-08

Testing `TVVARSS`
================
Lucas Nell
2017-02-07

This document outlines the results of the following:

1.  Using optimization and simulated annealing on `TVVARSS.ml` and `cpp_TVVARSS_ml`. The code to conduct these tests is in [this file](./test_files/TVVARSS.ml_testing.R).
2.  Overall output from `TVVARSS` and `cpp_TVVARSS`. The code to conduct these tests is in [this file](./test_files/TVVARSS_testing.R).

Comparisons between R and C++ versions used the `all.equal` function, which tests for 'near equality'. The tolerance for near equality is 1.49 × 10<sup>−8</sup>.

Optimization and simulated annealing
====================================

Optimization
------------

#### Using `optim(..., method = 'Nelder-Mead')`

    R version took 0.99 minutes
    Rcpp version took 0.04 minutes

    Are they equal?
    TRUE

#### Using `optim(..., method = 'BFGS')`

    R version took 3.1 minutes
    Rcpp version took 0.03 minutes

    Are they equal?
    Component "par": Mean relative difference: 2.466744
    Component "value": Mean relative difference: 1.639933
    Component "counts": Mean relative difference: 0.8042328

#### Using `bobyqa(...)`

    R version took 0.47 minutes
    Rcpp version took 0.02 minutes

    Are they equal?
    Component "par": Mean relative difference: 0.0003638024
    Component "fval": Mean relative difference: 0.001088931
    Component "feval": Mean relative difference: 0.2065455

Simulated annealing
-------------------

    R version took 8.17 minutes
    Rcpp version took 0.35 minutes

    Are they equal?
    Component "trace.mat": Mean relative difference: 2.740943e-08

Overall output from `TVVARSS` and `cpp_TVVARSS`
===============================================

Methods `'Nelder-Mead'`, `'BFGS'`, and `'bobyqa'` were run, with and without simulated annealing, to see how output and performance compared. Times below are for the 6 total runs for each version of `TVVARSS`:

    R version took 49.35 minutes
    cpp version took 2.21 minutes

Using `method = 'Nelder-Mead'`
------------------------------

#### With annealing

    Are they equal?
    TRUE

#### No annealing

    Are they equal?
    TRUE

Using `method = 'BFGS'`
-----------------------

#### With annealing

    Are they equal?
    TRUE

#### No annealing

    Are they equal?
    Modes: list, NULL
    Lengths: 43, 0
    names for target but not for current
    Attributes: < Modes: list, NULL >
    Attributes: < Lengths: 1, 0 >
    Attributes: < names for target but not for current >
    Attributes: < current is not list-like >
    current is not list-like

This looks weird because the R version here had an error, so returned `NULL`. Below is the error returned and the output from `traceback()` run afterward:

    Error in solve.default(FF) : 
      system is computationally singular: reciprocal condition number = 3.97576e-17

    6: solve.default(FF)
    5: solve(FF) at TVVARSS_forC_17Dec16.R#82
    4: fn(par, ...)
    3: (function (par) 
       fn(par, ...))(c(2442101.41194779, -3417820.0737009, -1299053.23694122, 
       3731843.14805551, -17889292.6552748, -3290188.25617733, 4547907.20822766, 
       -1474431.03066228, 270.996097653055, 1307468.69164327, 95126.5692356095, 
       -3088894.64420066, -3443284.04369214, 4977648.28422679, 150664.747666296, 
       1821492.68435145))
    2: optim(fn = TVVARSS.ml, par = par, X = tX, U = tU, par.fixed = par.fixed, 
           method = "BFGS", control = optim.control) at TVVARSS_forC_17Dec16.R#420
    1: TVVARSS(X, Tsamplefract = 0.9, method = "BFGS", show.fig = FALSE, 
           annealing = FALSE)

Using `method = 'bobyqa'`
-------------------------

#### With annealing

    Are they equal?
    Component "se": Mean relative difference: 0.06499206
    Component "su": Mean relative difference: 0.02594734
    Component "Sb0": Mean relative difference: 0.01249422
    Component "Sb": Mean relative difference: 0.1126508
    Component "B0": Mean relative difference: 0.02988222
    Component "B": Mean relative difference: 0.08901546
    Component "logLik": Mean relative difference: 0.01867835
    Component "AIC": Mean relative difference: 0.01768283
    Component "B0.fitted": Mean relative difference: 0.4543676
    Component "B.fitted": Mean relative difference: 0.4304262
    Component "X.fitted": Mean relative difference: 0.04126096
    Component "PP.fitted": Mean relative difference: 0.4624978
    Component "eigen.fitted": Modes: complex, numeric
    Component "eigen.fitted": Mean relative Mod difference: 0.3038433
    Component "opt.par": Mean relative difference: 0.04722132

#### No annealing

    Are they equal?
    Component "se": Mean relative difference: 0.01180563
    Component "su": Mean relative difference: 1.253447
    Component "Sb0": Mean relative difference: 7.997683
    Component "Sb": Mean relative difference: 0.02654061
    Component "B0": Mean relative difference: 1.861446
    Component "B": Mean relative difference: 0.02131687
    Component "logLik": Mean relative difference: 0.001173521
    Component "AIC": Mean relative difference: 0.001081522
    Component "B0.fitted": Mean relative difference: 1.790884
    Component "B.fitted": Mean relative difference: 0.1685675
    Component "X.fitted": Mean relative difference: 0.01180828
    Component "PP.fitted": Mean relative difference: 0.2763867
    Component "eigen.fitted": Mean relative Mod difference: 0.2197468
    Component "opt.par": Mean relative difference: 0.07720379

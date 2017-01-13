Testing `TVVARSS`
================
Lucas Nell
2017-01-13

This document outlines the results of the following:

1.  Using optimization and simulated annealing on `TVVARSS.ml` and `cpp_TVVARSS_ml`. The code to conduct these tests is in [this file](./initial_testing_files/TVVARSS.ml_testing.R).
2.  Overall output from `TVVARSS` and `cpp_TVVARSS`. The code to conduct these tests is in [this file](./initial_testing_files/TVVARSS_testing.R).

Comparisons between R and C++ versions used the `all.equal` function, which tests for 'near equality'. The tolerance for near equality is 1.49 × 10<sup>−8</sup>.

Optimization and simulated annealing
====================================

Optimization
------------

#### Using `optim(..., method = 'Nelder-Mead')`

    R version took 1 minutes
    Rcpp version took 0.01 minutes

    Are they equal?
    Component "par": Mean relative difference: 0.244417
    Component "value": Mean relative difference: 0.04055309
    Component "counts": Mean relative difference: 0.886023

#### Using `optim(..., method = 'BFGS')`

    R version took 2.92 minutes
    Rcpp version took 0.02 minutes

    Are they equal?
    Component "par": Mean relative difference: 2.608607
    Component "value": Mean relative difference: 2.524436
    Component "counts": Mean relative difference: 0.872134

#### Using `bobyqa(...)`

    R version took 0.45 minutes
    Rcpp version took 0.01 minutes

    Are they equal?
    Component "par": Mean relative difference: 0.0001922985
    Component "fval": Mean relative difference: 0.002921203
    Component "feval": Mean relative difference: 0.6392727

Simulated annealing
-------------------

    R version took 8.98 minutes
    Rcpp version took 0.58 minutes

    Are they equal?
    Component "trace.mat": Mean relative difference: 2.694283e-08

Overall output from `TVVARSS` and `cpp_TVVARSS`
===============================================

Methods `'Nelder-Mead'`, `'BFGS'`, and `'bobyqa'` were run, with and without simulated annealing, to see how output and performance compared. Times below are for the 6 total runs for each version of `TVVARSS`:

    R version took 96.69 minutes
    cpp version took 4.05 minutes

Using `method = 'Nelder-Mead'`
------------------------------

#### With annealing

    Are they equal?
    Component "se": Mean relative difference: 7.369287e-06
    Component "su": Mean relative difference: 3.267878e-05
    Component "Sb0": Mean relative difference: 3.174063e-06
    Component "Sb": Mean relative difference: 7.188568e-06
    Component "B0": Mean relative difference: 1.604839e-05
    Component "B": Mean relative difference: 7.586208e-06
    Component "logLik": Mean relative difference: 2.343666e-06
    Component "AIC": Mean relative difference: 2.22845e-06
    Component "B0.fitted": Mean relative difference: 0.00202612
    Component "B.fitted": Mean relative difference: 0.00227323
    Component "X.fitted": Mean relative difference: 0.0001535365
    Component "PP.fitted": Mean relative difference: 0.001932475
    Component "eigen.fitted": Mean relative Mod difference: 0.002329896
    Component "opt.par": Mean relative difference: 9.404633e-06

#### No annealing

    Are they equal?
    Component "se": Mean relative difference: 0.07513969
    Component "su": Mean relative difference: 1.289825
    Component "Sb0": Mean relative difference: 0.6557525
    Component "Sb": Mean relative difference: 0.8638878
    Component "B0": Mean relative difference: 0.2299503
    Component "B": Mean relative difference: 0.1497434
    Component "logLik": Mean relative difference: 0.005663374
    Component "AIC": Mean relative difference: 0.005226255
    Component "B0.fitted": Mean relative difference: 1.122417
    Component "B.fitted": Mean relative difference: 0.2494728
    Component "X.fitted": Mean relative difference: 0.05635341
    Component "PP.fitted": Mean relative difference: 0.8868769
    Component "eigen.fitted": Mean relative Mod difference: 0.1633977
    Component "opt.par": Mean relative difference: 0.2691596

Using `method = 'BFGS'`
-----------------------

#### With annealing

    Are they equal?
    TRUE

R and C++ versions had output all within 1.49 × 10<sup>−8</sup> of each other.

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

This looks weird because the R version here had an error, so returned `NULL`. Below is the error returned and the output from `traceback()` run afterward::

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
    Component "se": Mean relative difference: 0.2280353
    Component "su": Mean relative difference: 0.03555844
    Component "Sb0": Mean relative difference: 0.06097886
    Component "Sb": Mean relative difference: 0.3197025
    Component "B0": Mean relative difference: 0.1698673
    Component "B": Mean relative difference: 0.57555
    Component "logLik": Mean relative difference: 0.06561574
    Component "AIC": Mean relative difference: 0.06239371
    Component "B0.fitted": Mean relative difference: 0.5248492
    Component "B.fitted": Mean relative difference: 0.9133515
    Component "X.fitted": Mean relative difference: 0.05170573
    Component "PP.fitted": Mean relative difference: 0.6607509
    Component "eigen.fitted": Modes: complex, numeric
    Component "eigen.fitted": Mean relative Mod difference: 0.8285955
    Component "opt.par": Mean relative difference: 0.2463107

#### No annealing

    Are they equal?
    Component "se": Mean relative difference: 0.01664413
    Component "su": Mean relative difference: 0.3265023
    Component "Sb0": Mean relative difference: 0.9173934
    Component "Sb": Mean relative difference: 0.2940105
    Component "B0": Mean relative difference: 2.451531
    Component "B": Mean relative difference: 0.08640339
    Component "logLik": Mean relative difference: 0.01676588
    Component "AIC": Mean relative difference: 0.01547324
    Component "B0.fitted": Mean relative difference: 0.9571879
    Component "B.fitted": Mean relative difference: 0.4245081
    Component "X.fitted": Mean relative difference: 0.007991256
    Component "PP.fitted": Mean relative difference: 0.8627313
    Component "eigen.fitted": Mean relative Mod difference: 0.502622
    Component "opt.par": Mean relative difference: 0.2287242

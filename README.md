TVVARSS
========

##### Time-varying vector autoregressive state-space model


### This folder

The files in this folder include the following:

- [`TVVARSS.cpp`](./TVVARSS.cpp): 
  C++ version of the `TVVARSS.ml` function (`cpp_TVVARSS_ml`)
- [`TVVARSS_cpp.R`](./TVVARSS_cpp.R): 
  the C++ version of `TVVARSS` (`cpp_TVVARSS`) that uses `cpp_TVVARSS_ml`
- [`TVVARSS_testing.Rmd`](./TVVARSS_testing.Rmd): 
  comparisons of C++ to R versions of `TVVARSS` and `TVVARSS.ml` functions
- [`TVVARSS_testing.md`](./TVVARSS_testing.md): 
  output from the above .Rmd file (*this is the one to view on GitHub*)




### `test_files` subfolder


The `test_files` folder contains the following files:

- [`TVVARSS.ml.R`](./test_files/TVVARSS.ml.R):
  current version of the inner function `TVVARSS.ml`, used for testing
- [`TVVARSS.ml_pars.R`](./test_files/TVVARSS.ml_pars.R):
  initial setup of parameters used to test `TVVARSS.ml`
- [`TVVARSS.ml_testing.R`](./test_files/TVVARSS.ml_testing.R):
  testing of R and C++ versions of `TVVARSS.ml` function
- [`TVVARSS.ml_testing.RData`](./test_files/TVVARSS.ml_testing.RData):
  resulting objects from tests of R and C++ versions of `TVVARSS.ml` function
- [`TVVARSS_forC_17Dec16.R`](./test_files/TVVARSS_forC_17Dec16.R):
  original R version of `TVVARSS` used for comparison to C++ version
- [`TVVARSS_sim_test_17Dec16.R`](./test_files/TVVARSS_sim_test_17Dec16.R):
  initial simulation test conducted on the original R version by ARI
- [`TVVARSS_testing.R`](./test_files/TVVARSS_testing.R):
  tests comparing output from R and C++ versions of `TVVARSS`
- [`TVVARSS_testing.RData`](./test_files/TVVARSS_testing.RData):
  resulting objects from tests comparing output from R and C++ versions of `TVVARSS`




### Libraries and packages

The C++ implementation uses the "Armadillo" C++ library and interfaces with R using
the packages `Rcpp` and `RcppArmadillo`.



# Compiling

You may get the following error message when running `Rcpp::sourceCpp('TVVARSS.cpp')`:

```r
ld: warning: directory not found for option '-L/usr/local/lib/gcc/x86_64-apple-darwin13.0.0/4.8.2'
```

### If homebrew is not installed

You can install gfortran the following way in bash (i.e., in the Terminal app on macOS)
([link to solution source](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)):

```bash
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

### If homebrew is installed

If you have homebrew installed, you'll have to do this instead
([link to solution on StackOverflow](http://stackoverflow.com/a/35947039/5016095)):

```bash
brew reinstall gcc48 --with-fortran
mkdir -p ~/.R
printf "FLIBS=\"\"\nF77=\"gfortran-4.8\"\nFC=\"gfortran-4.8\"\n" >> ~/.R/Makevars
```

The ~/.R/Makevars file should look like this:

```
FLIBS=""
F77="gfortran-4.8"
FC="gfortran-4.8"
```

Now restart R.


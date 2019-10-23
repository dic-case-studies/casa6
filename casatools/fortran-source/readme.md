
## Basic Fortran Libraries

[casacode](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa/browse) and [casacore](https://github.com/casacore/casacore) each contain a small amount of Fortran code used in compute intensive areas of the code base. It is believed that using Fortran in these locations provides performance improvements.

At this point, the only real option for open-source Fortran compilers is [gfortran](https://gcc.gnu.org/wiki/GFortran). The unfortunate side effect of using gfortran is that most applications built with gfortran are required to link against a shared library (libgfortran) that comes with gfortran and supports some of the capabilities underlying the language. When building and running on a single host, this does not pose a problem. However, when one needs to build on one host and run on many Linux hosts, then the problems begin. These problems are primarily due to the gradual changes to libgfortran.

Our problems are compounded by the fact that we rely on [numpy](https://www.numpy.org), and numpy also depends upon libgfortran. In the general case, loading libgfortran from multiple versions of gfortran does not work. For this reason, casatools setup.py has an option to link casatools (as well as the other Fortran based libraries) against the libgfortran that comes with numpy.

Currently the libraries that are optionally built against the numpy libgfortran:

* [rpfits](http://www.atnf.csiro.au/computing/software/rpfits.html)
* [lapack](http://www.netlib.org/lapack/) and [blas](http://www.netlib.org/blas/)

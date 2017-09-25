
## CASAtools

CASAtools is a self-contained python module that provides the tools from the [CASA](http://casa.nrao.edu/) project. The module only contains the non-GUI tools which are the [SWIG](http://swig.org) bound C++ functionality from CASA.

## Maturity

This is alpha software in every sense of the word. Currently the goal is just to make this usable on the Linux platforms which the [CASA](http://casa.nrao.edu) project supports for development, [RedHat Enterprise Linux](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) release 6 and release 7.

## Check Out

* git clone https://open-bitbucket.nrao.edu/scm/casa/CASAtools.git

## Build

The CASA development RPMs must be installed. This can be accomplished by adding the CASA [YUM](https://en.wikipedia.org/wiki/Yum_(.rpm)) repository to /etc/yum.repos.d and installing:

* -bash-4.2# yum install casa-toolset-2

After the CASA build environment is installed, the CASAtools module can be built like:

* scripts/gcw-pick
* aclocal
* autoconf
* ./configure
* ./setup.py build


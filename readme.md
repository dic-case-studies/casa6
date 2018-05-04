
## CASAtasks

If you are looking for the place to check out [CASA](http://casa.nrao.edu/) with [Git](https://en.wikipedia.org/wiki/Git), this is not the right place. The repository to check out [CASA](http://casa.nrao.edu/) can be found [here](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa/browse). At some point in the future, this may form the underlying kernel for a future version of [CASA](http://casa.nrao.edu/), but at this point, this package is still in initial development.

CASAtasks is a self-contained python module that provides the tasks from the [CASA](http://casa.nrao.edu/) project. This package depends on the [CASAtools](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse) python module being found in your **PYTHONPATH** at build time and runtime. The CASAtasks are stateless routines and recipes built on CASAtools.

## Maturity

This is is still in intial development (read *pre-alpha* software), and at this point, it is not intended for use outside of the CASA project. Currently it builds wherever [CASAtools](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse) can be built, and the only extra Python dependency is [matplotlib](https://matplotlib.org). Most of the development of this module has been done with Python 3.6, but we plan to make an effort to get this to work with Python 2.7.

## Building CASAtasks

#### Install Dependencies

First make sure that the version of Python you want to use is available. [MacPorts](https://www.macports.org) has been used to successfully install Python 3.6 on [OSX](http://en.wikipedia.org/wiki/MacOS). On [RedHat](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux), the [Software Collections](https://developers.redhat.com/products/softwarecollections/overview/) has been used to install Python 3.6.:
```
-bash-4.2# yum install rh-python36
-bash-4.2# yum install rh-python36-python-numpy
-bash-4.2# yum install rh-python36-python-scipy
-bash-4.2# yum install rh-python36-python-tkinter
```
Unfortunately, Software Collections does not supply an RPM for [matplotlib](https://matplotlib.org). For installing matplotlib, [pip](https://pypi.org/project/pip) can be used:
```
-bash-4.2# pip install matplotlib
```

Second build and install [CASAtools](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse) using your choice of Python and the instructions found [here](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse).

#### Checkout

Checkout the CASAtasks source code:

```
-bash-4.2$ git clone -q --recursive https://open-bitbucket.nrao.edu/scm/casa/CASAtasks.git
```

#### Build

After all of the dependencies have been installed and the source code for CASAtasks is available we can build CASAtasks. Make sure that =which python= returns the version of python that was used to build [CASAtools](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse). Then build the tasks with:
```
-bash-4.2$ cd CASAtasks
-bash-4.2$ PYTHONPATH=../CASAtools/build/lib.macosx-10.12-x86_64-3.6 ./setup.py build
```
**Substitute** the path to your build of [CASAtools](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse) in the build line above.

#### Run Available Tests

Tests are actively being added, but you can run the existing tests with:
```
-bash-4.2$ PYTHONPATH=../CASAtools/build/lib.macosx-10.12-x86_64-3.6 ./setup.py test
```

Again, **substitute** the path to your build of [CASAtools](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtools/browse) in the line above.

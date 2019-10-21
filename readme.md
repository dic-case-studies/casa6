## CASA 6

CASA 6 marks a major new chapter in CASA's development. With CASA 6, tools and
tasks will be available as Python wheels from a [PyPI repository](https://casa-pip.nrao.edu/).
These wheels can then be installed in any version of Python 3.6. In time CASA will
provide wheels for other versions of Python.

CASA 6 divides CASA's functionality into two Python packages:

  * [casatools](casatools/readme.md) -- the C++ tools with a minimal Python layer
  * [casatasks](casatasks/readme.md) -- a pure Python layer which provides a higher level of abstraction

These packages can be installed in Python from NRAO's [PyPI repository](https://casa-pip.nrao.edu/).
The recommended way to do this is using a virtual environment:
```
-bash-4.2$ python3 -m venv casa6
-bash-4.2$ source casa6/bin/activate
-bash-4.2$ pip install --extra-index-url https://casa-pip.nrao.edu:443/repository/pypi-group/simple casatasks
-bash-4.2$ #sanity check
           #              -bash-4.2$ python
           #              >>> import casatasks as ct
           #              >>> help(ct)
           #              >>> exit()
           #              -bash-4.2$ 
-bash-4.2$ deactivate
```
However, users are free to install CASA's wheels as they like.

In general, [casatasks](casatasks/readme.md) is identical to the CASA 5 version of the tasks.
The primary difference is that the ```inp```, ```go```, ```default```, and *subparameters*
are no longer part of the tasks. This functionality has been moved up into another module
called [casashell](https://open-bitbucket.nrao.edu/projects/CASA/repos/casashell/). This
module provides an environment that is very similar to the standard CASA prompt. In this
layer parameters can be specified by setting global variables, checked with ``inp``, and
the task can be invoked with ```go```. This change was made because this functionality
depends upon [IPython](https://ipython.org/) and the fetching, setting, and clearing of
global variables by access to the task call stack. The former may not be available and
the latter is problematic for packages loaded into the user's Python installation. For
users who want this environment, it is accessible with the casashell module either from
the bash command line if the casashell module has been installed:
```
-bash-4.2$ python3 -m casashell
```
or from Python 3:
```
-bash-4.2$ python3
Python 3.6.7 (default, Aug 16 2019, 18:25:33) 
[GCC 5.3.1 20160406 (Red Hat 5.3.1-6)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from casashell import start_casa
>>> start_casa([])
IPython 7.8.0 -- An enhanced Interactive Python.

Using matplotlib backend: TkAgg
in startup
CASA 5.9.9-922 -- Common Astronomy Software Applications [2019.214]

CASA <1>:
```

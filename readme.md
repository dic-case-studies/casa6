## CASA 6

CASA 6 marks a major new chapter in CASA's development. With CASA 6, tools and
tasks will be available as Python wheels from a [PyPI repository](https://casa-pip.nrao.edu/).
These wheels can then be installed in any version of Python 3.6. In time CASA will
provide wheels for other versions of Python.

### Organization

CASA 6 divides CASA's functionality into two Python packages:

  * [casatools](casatools/readme.md) -- the C++ tools with a minimal Python layer
  * [casatasks](casatasks/readme.md) -- a pure Python layer which provides a higher level of abstraction

[casatasks](casatasks/readme.md) is very similar to the CASA 5 tasks (but see discussion of inp/go below).
However, with [casatools](casatools/readme.md) there was some redundancy in CASA 5. This redundancy has
been removed. For example, the image analysis tool has always been called ```image``` but it was exposed
in CASA 5 as ```ia``` (and public instance of the ```image``` tool) and ```iatool``` a synonym for
```image```. The instance and synonym are no longer available for CASA 6 tools, but it is easy to import
the CASA 6 tool in the CASA 5 style:
```
>>> from casatools import image as iatool
>>> ia = iatool( )
```
These changes are sumarized in [this document](https://casa.nrao.edu/download/devel/docs/casa6/CASA-ToolNames.pdf).

### Installation

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

### Inp/Go and Subparameters

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

### CASA's GUIs

The primary GUIs [casaviewer](https://open-bitbucket.nrao.edu/projects/CASA/repos/casaviewer/browse)
and [casaplotms](https://open-bitbucket.nrao.edu/projects/CASA/repos/casaplotms/browse) are available
as separate modules. These allow the GUIs to be loaded and used.

Some of the GUI tools available as part of CASA 5 are only available in the packaged tar-file based
disribution of CASA 6. The GUIs that are currently only available in the monolithic vesion of
CASA 6 currently includes:

  1. casabrowser
  2. casafeather
  3. casalogger
  4. casaplotserver


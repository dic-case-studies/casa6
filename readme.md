
## CASAtools

CASAtools is a self-contained python module that provides the tools from the [CASA](http://casa.nrao.edu/) project. The module only contains the non-GUI tools which are the [SWIG](http://swig.org) bound C++ functionality from CASA.

## Maturity

This is alpha software in every sense of the word. Currently the goal is just to make this usable on the Linux platforms which the [CASA](http://casa.nrao.edu) project supports for development, [RedHat Enterprise Linux](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) release 6 and release 7.

## Check Out

* git clone https://open-bitbucket.nrao.edu/scm/casa/CASAtools.git

## Build

The CASA development RPMs must be installed. This can be accomplished by adding the CASA [YUM](https://en.wikipedia.org/wiki/Yum_(.rpm)) repository to /etc/yum.repos.d and installing:

* -bash-4.2# yum install ccache
* -bash-4.2# yum install casa-toolset-2

After the CASA build environment is installed, the CASAtools module can be built like:

* scripts/gcw-pick
* aclocal
* autoconf
* ./configure
* ./setup.py build

## Available Tools

| Tool Name           | Description                                         |
| ------------------- | --------------------------------------------------- |
| atmosphere          | Atmosphere model                                    |
| componentlist       | A tool for the manipulation of groups of components |
| coordsys            | Operations on CoordinateSystems                     |
| functional          | Functionals handling                                |
| imagepol            | Polarimetric analysis of images                     |
| image               | Operations on images                                |
| measures            | measures tool                                       |
| quanta              | quanta tool handles units and quantities            |
| regionmanager       | Create and manipulate regions of interest           |
| synthesisimager     | tool for synthesis imaging                          |
| synthesisimstore    | tool for synthesis imaging                          |
| synthesisnormalizer | tool for synthesis imaging                          |
| synthesisutils      | tool for synthesis imaging                          |
| table               | Access tables from casapy                           |

## Changes from Standard CASA

While the goal was to simply reconstitute the [CASA tools](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa/browse/gcwrap/tools) within a unencumbered python module, deviations were required as work progressed. These deviations are divided into categories based upon whether the change relates to the way the tools behave or the way the XML files are structured.

### Changes to Behavior

1. __Parameter Order__ -- with standard CASA, the [SWIG](http://swig.org) binding is exposed directly to the user. While this is often fine, for CASAtools more control was desired over how parameter checking happens and how coersion happens. Because of this, a wrapper is created over the SWIG binding. This wrapper explicitly specifies the order of the parameters. The input parameters come first, followed by the output parameters. Within the input and output parameters, the order is determined by the XML file. With standard CASA, all parameters are largely unordered with binding controlled by the way [python argument lists](https://stackoverflow.com/questions/3394835/args-and-kwargs) are handled.

### XML Changes

1. __String Constants__ -- (**developer**) default values for strings, should **not** include quotes. The standard CASA XML processing would strip out opening and closing quotes, but CASAtools XML processing does not.

2. __Path Type Added__ -- (**developer**/**automatic**) a new parameter type called *path* was added to CASAtools. This *path* type evolved from the standard CASA convention of adding the *mustexist* attribute to string typed parameters. The *automatic* XML translation converts string parameters which include the *mustexist* attribute to path parameters. C++ receives the *path* parameters as a string, but when the *mustexist* attribut of the *path* parameter is set to *true*, the CASAtools binding layer will expand the required path using the **CASADATA** (colon separated) path list to locate the required file or directory. In the cases where the *mustexist* attribute was not used, developers must make any changes to convert a *string* parameter to a *path* parameter.

3. __Type Element Added__ -- (**developer**) CASAtools adds a new *<type>* XML element which has precedence over the *type* attribute, e.g. *<param type="...">*. This addition allows types that can be passed to a *any*/*variant* param to be enumerated. This allows *path* strings to be expanded by the CASAtools binding layer when the *path* is part of an *any*/*variant* parameter. When more than one *<type>* is specified for a parameter, the parameter is assumed to be an *any*/*variant* parameter.

4. __Bool__ -- (**automatic**) The majority of standard CASA XML files indicate boolean type parameters with *bool*. However, a small number of them use *boolean*. With CASAtools, only *bool* is accepted.

5. __Array Values__ -- (**automatic**) With CASAtools, the behavior *<value>* elements for vector initialization has been rationalized:
  * **<value/>** -- empty vector (zero elements)
  * **<value><value/></value>** -- vector with one element initalized to the default initialization for the vector element type
  * **```xml <value><value>0.0</value></value>```** -- vector with one element initialized as specified, more internal *<value>* elements can be used to increase the default size of the vector


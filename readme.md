
## CASAtools

If you are looking for the place to check out [CASA](http://casa.nrao.edu/) with [Git](https://en.wikipedia.org/wiki/Git), this is not the right place. The repository to check out [CASA](http://casa.nrao.edu/) can be found [here](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa/browse). At some point in the future, this may form the underlying kernel for a future version of [CASA](http://casa.nrao.edu/), but at this point, this package is still in gestation.

CASAtools is a self-contained python module that provides the tools from the [CASA](http://casa.nrao.edu/) project. This module only contains the non-GUI tools which are the [SWIG](http://swig.org) bound C++ functionality from CASA.

CASAtools can be built and used with either **Python** **2.7**, **3.4** or **3.6**. However, *[CASAtasks](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtasks/browse) can only be built with Python 3*. So if you plan to build [CASAtasks](https://open-bitbucket.nrao.edu/projects/CASA/repos/CASAtasks/browse), please ensure that you have a Python 3 interpreter available.

## Maturity

This is alpha software in every sense of the word. Currently it builds on [RedHat Enterprise Linux](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) and OSX (RHEL6, RHEL7 and OSX 10.12) using the standard [CASA](http://casa.nrao.edu) development environment.

## Building CASAtools

#### Install Dependencies

Those already working on CASA 5 can install all of the CASAtools dependencies with [YUM](https://en.wikipedia.org/wiki/Yum_(.rpm)):

```
-bash-4.2$ yum install casa-toolset-3
```
After doing this, you may need to add RedHat's Python 3 to your $PATH with:
```
-bash-4.2$ source /opt/rh/rh-python36/enable
```
In addition, [Java 8](https://java.com/en/download/) or greater is required. You can check this by looking for 1.8 in the Java version output:
```
-bash-4.2$ /usr/bin/java -version
java version "1.8.0_121"
Java(TM) SE Runtime Environment (build 1.8.0_121-b13)
Java HotSpot(TM) 64-Bit Server VM (build 25.121-b13, mixed mode)
-bash-4.2$
```

##### Notable Dependencies
1. [GNU Scientific Library](https://www.gnu.org/software/gsl/) version **2.2** or greater
1. [SWIG](http://www.swig.org) version **3.0** or greater
1. [OpenMPI](https://www.open-mpi.org) version 1.10 *required when parallelization is enabled*


#### Checkout

Checkout the CASAtools source code:

```
-bash-4.2$ git clone -q --recursive https://open-bitbucket.nrao.edu/scm/casa/casatools.git
```

#### Build

While CASAtools can be built with either Python 2.7 or Python 3, [CASAtasks](https://open-bitbucket.nrao.edu/projects/CASA/repos/casatasks/browse) requires Python 3 so it is probably best to build with Python 3 from the start. Make sure that ```-bash-4.2$ which python``` returns a path to the Python 3 executable you want to use because [CASAtasks](https://open-bitbucket.nrao.edu/projects/CASA/repos/casatasks/browse) builds with this executable.

With the CASA build environment set up, the CASAtools module can be built like:
```
-bash-4.2$ cd casatools
-bash-4.2$ scripts/gcw-pick
-bash-4.2$ autoconf
-bash-4.2$ ./configure
-bash-4.2$ ./setup.py build
```
The `gcw-pick` script adjusts the standard CASA source tree for building with `setup.py`, and once CASAtools is integrated with CASA this step will not be necessary. `gcw-pick` may run for quite a while...

A particular version of Python can be selected at configure time like:
```
-bash-4.2$ PYTHON_VERSION=3.4 ./configure
```
but as noted above, [CASAtasks](https://open-bitbucket.nrao.edu/projects/CASA/repos/casatasks/browse) will use whatever version of python is in your path so it is best to configure your ```PATH``` to get a particular version of Python.

#### Alternate Developer Build For Linux

On Linux, one can alternately use **make**, which supports incremental and parallel builds that are useful for development. In this case, the steps are
```
-bash-4.2$ cd casatools
-bash-4.2$ scripts/gcw-pick
-bash-4.2$ autoconf
-bash-4.2$ ./configure
-bash-4.2$ ./setup.py genmake [--debug]
-bash-4.2$ make
```
The **--debug** parameter of the **./setup.py genmake** is optional. If provided, the resulting makefile is configured to use the -g option when compiling C++ code. Otherwise, the makefile is configured to produce an optimized build using the -O2 option when compiling C++ files. The **make** command optionally takes the **-j** option with the number of parallel threads to use while building which can significantly improve build times on multi-core machines, eg
```
-bash-4.2$ make -j8
```
To execute a parallel build with eight threads.

The **./setup.py genmake** command creates a file named *makefile* with the make commands. In general, this file should only be modified if one knows what they are doing. Once the *makefile* has been generated, one can make changes to most source files (eg, .cc and .h files) and simply rerun the **make** command to generate an updated build. However, there are cases when earlier commands will have to be rerun. For example, if *ac/templates/setup.py.in* is modified, commands starting with **autoconf** will have to be rerun to generate a new *setup.py* file.

A non-exhaustive list of when **./setup.py genmake** followed by **make** will need to be rerun includes:
* modification of cerberus files. These are copied by **./setup.py genmake**; **make** does not deal with them
* changes in tool dependencies when modifying *&lt;tool&gt;.xml* files. If one does not change tool dependencies when modifying *&lt;tool&gt;.xml* files, they need only rerun **make**. However, if a new tool dependency is introduced or an existing dependency removed (eg via a function return value), one will have to start by rerunning **./setup.py genmake** so the new dependency can propagate to the the *makefile*.

#### See If It Works

If the build completes successfully, try loading the CASAtools module:

```
-bash-4.2$ PYTHONPATH=build/lib.linux-x86_64-3.6 python3
Python 3.6.3 (default, Jan  9 2018, 10:19:07)
[GCC 4.8.5 20150623 (Red Hat 4.8.5-11)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from casatools import image
>>> ia = image( )
>>> ia.fromshape("mytest.im",[20,20])
2019-03-25 16:15:08     INFO    ImageFactory::createImage       Created Paged image 'mytest.im' of shape [20, 20] with float valued pixels.
True
>>>
-bash-4.2$
```
CASAtools can be configured using *~/.casa/toolrc.py*. For example, extra data search paths can be added like:
```
-bash-4.2$ cat ~/.casa/toolrc.py
datapath=[ "~/develop/casa/data/unittests" ]
-bash-4.2$
```

#### Run Available Tests

A number of tests have been brought into CASAtools from CASA. Like the rest of CASAtools, the testing infrastructure is still in the process of being refined (and in need of some refactoring). However, with the proper setup it works for *some of the people, some of the time*. To run the tests, you need to check out the subset of the CASA data repository needed for CASAtools unit tests. This data repository will use about **14GB** of space. You can check out the necessary data repository like:
```
-bash-4.2$ git clone --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casa-data.git unittests
-bash-4.2$ cd unittests
-bash-4.2$ git show HEAD:casatools-unittests | bash
```
After the checkout is complete, you must update your CASAtools RC file to indicate where the unit test data can be found. In my case it looks like:
```
-bash-4.2$ cat ~/.casa/toolrc.py
datapath=[ "~/develop/casa/data/unittests" ]
-bash-4.2$
```
The `datapath` list specifies directories where CASAtools should look for data files, and should include the directory we just checked out. After all of this is set, the full suite of tests can be run with:
```
-bash-4.2$ ./setup.py test
```
There are bout 48 tests in total. Running all of them takes about 12 minutes on my laptop.

*To run a single set of test*, you can just execute the test file:
```
-bash-4.2$ PYTHONPATH=build/lib.macosx-10.12-x86_64-3.6 python tests/tools/coordsys/test_coordsys.py
```
If the ```casatools``` module is already available in your ```PYTHONPATH```, you would not need to
explicitly specify it on the command line. This would be the case, for example, if you installed
CASAtools via a binary *PyPI* wheel. If you want to run a single test from the set, you can specify
it on the command line like:
```
-bash-4.2$ PYTHONPATH=build/lib.linux-x86_64-3.6 python tests/tools/coordsys/test_coordsys.py coordsys_test.test_constructor
```
Here, ```test_constructor``` is one test within the ```coordsys_test``` *TestCase* specified in ```test_coordsys.py```.

#### Notes

If some time has passed since the last build, you should (sometimes) remove *xml-casa-assembly-1.7.jar*, e.g.
```
-bash-4.2$ rm ./scripts/java/xml-casa-assembly-1.7.jar
-bash-4.2$ scripts/gcw-pick
```
before rebuilding because this [JAR](https://en.wikipedia.org/wiki/JAR_(file_format)) file is automatically fetched from a download site by *gcw-pick*. However, it is not fetched if it already exists. Deleting the current copy will result in a new copy being fetched which *may* be newer.

## Available Tools

| Tool Name            | Description                                                     |
| -------------------- | --------------------------------------------------------------- |
| agentflagger         | Tool for manual and automated flagging                          |
| atcafiller           | Filler for ATNF/ATCA RPFITS data                                |
| atmosphere           | Atmosphere model                                                |
| calanalysis          | Get and fit data from a calibration table (CASA 3.4 and later). |
| calibrater           | Synthesis calibration (self- and cross-)                        |
| componentlist        | A tool for the manipulation of groups of components             |
| coordsys             | Operations on CoordinateSystems                                 |
| functional           | Functionals handling                                            |
| imagemetadata        | Operations on imagemetadata                                     |
| imagepol             | Polarimetric analysis of images                                 |
| imager               | tool for synthesis imaging                                      |
| image                | Operations on images                                            |
| iterbotsink          | tool for synthesis imaging                                      |
| logsink              | tool for logsink                                                |
| measures             | measures tool                                                   |
| miriadfiller         | Tool for the importmiriad task                                  |
| msmetadata           | Operations to retrieve metadata from a measurment set           |
| mstransformer        | Tool to apply spw and frame transformations in MS               |
| ms                   | Operations on measurement sets                                  |
| quanta               | quanta tool handles units and quantities                        |
| regionmanager        | Create and manipulate regions of interest                       |
| sakura               | New single dish tool interface using sakura                     |
| sdm                  | Manipulate or examine SDM datasets                              |
| simulator            | Tool for simulation                                             |
| singledishms         | New single dish tool interface to process an MS                 |
| spectralline         | spectral line tool                                              |
| synthesisdeconvolver | tool for synthesis imaging                                      |
| synthesisimager      | tool for synthesis imaging                                      |
| synthesisimstore     | tool for synthesis imaging                                      |
| synthesisnormalizer  | tool for synthesis imaging                                      |
| synthesisutils       | tool for synthesis imaging                                      |
| table                | Access tables from casapy                                       |
| vlafiller            | tool for VLA filler tasks                                       |
| vpmanager            | Tool for specifying voltage patterns and primary beams          |


## Tool Initialization

The user initalization and customization file for the CASAtools module is `~/.casa/toolrc.py`. If this file does not exist, then an attempt is mad to import the RC values from casatoolrc (i.e. `from casatoolrc import *`). If the `casatoolrc` module does not exist in `PYTHONPATH`, then the default values for all initialization state is used.

## Changes from Standard CASA

While the goal was to simply reconstitute the [CASA tools](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa/browse/gcwrap/tools) within a unencumbered python module, deviations were required as work progressed. These deviations are divided into categories based upon whether the change relates to the way the tools behave or the way the XML files are structured. The *gcw-pick* script includes an automatic XML update function which is able to make some of the more basic changes to the XML code automatically as the XML files are pulled from the *casa-source/gcwrap* tree into CASAtools.

### Changes to Behavior

1. __Parameter Order__ --- with standard CASA, the [SWIG](http://swig.org) binding is exposed directly to the user. While this is often fine, for CASAtools more control was desired over how parameter checking happens and how coersion happens. Because of this, a wrapper is created over the SWIG binding. This wrapper explicitly specifies the order of the parameters. The input parameters come first, followed by the output parameters. Within the input and output parameters, the order is determined by the XML file. With standard CASA, all parameters are largely unordered with binding controlled by the way [python argument lists](https://stackoverflow.com/questions/3394835/args-and-kwargs) are handled.

2. __Utils Tool__ --- the standard CASA *uttool* (AKA *utils* or just *ut*), has been converted from a toolbox/grab bag **tool** to a **singleton object** called *ctsys* in CASAtools, e.g.:
    ```
    In [1]: from casatools import ctsys
    ```

3. __ctsys.resolve( )__ --- a new member function was added to resolve the path to an data file based upon **CASADATA** path (as is done for `<type mustexist="true">path</type>`)

3. __rc file__ --- the rc file, which is evaluated at startup to configure CASAtools, is `~/.casa/toolrc.py`; if this file is not found, then an attempt is mad to import the RC values from casatoolrc (i.e. `from casatoolrc import *`)

### Xml Changes

1. __String Constants__ --- (**developer**) default values for strings, should **not** include quotes. The standard CASA XML processing would strip out opening and closing quotes, but CASAtools XML processing does not.

2. __Path Type Added__ --- (**developer**/**automatic**) a new parameter type called *path* was added to CASAtools. This *path* type evolved from the standard CASA convention of adding the *mustexist* attribute to string typed parameters. The *automatic* XML translation converts string parameters which include the *mustexist* attribute to path parameters. C++ receives the *path* parameters as a string, but when the *mustexist* attribut of the *path* parameter is set to *true*, the CASAtools binding layer will expand the required path using the **CASADATA** (colon separated) path list to locate the required file or directory. In the cases where the *mustexist* attribute was not used, developers must make any changes to convert a *string* parameter to a *path* parameter.

3. __Type Element Added__ --- (**developer**) CASAtools adds a new *`<type>`* XML element which has precedence over the *type* attribute, e.g. *`<param type="...">`*. This addition allows types that can be passed to a *any*/*variant* param to be enumerated. This allows *path* strings to be expanded by the CASAtools binding layer when the *path* is part of an *any*/*variant* parameter. When more than one *`<type>`* is specified for a parameter, the parameter is assumed to be an *any*/*variant* parameter.

4. __Bool__ --- (**automatic**) The majority of standard CASA XML files indicate boolean type parameters with *bool*. However, a small number of them use *boolean*. With CASAtools, only *bool* is accepted.

5. __StringArray__ --- (**automatic**) The majority of standard CASA XML files indicate the vector of strings type parameters with *stringArray*. However, a small number of them use *stringarray*. With CASAtools, only *stringArray* is accepted.

5. __Array Values__ --- (**developer**) With CASAtools, the behavior *`<value>`* elements for vector initialization has been rationalized:
    * **`<value/>`** --- empty vector (zero elements)
    * **`<value><value/></value>`** --- vector with one element initalized to the default initialization for the vector element type
    * **`<value><value>0.0</value></value>`** --- vector with one element initialized as specified, more internal *`<value>`* elements can be used to increase the default size of the vector

   It is not possible to automatically adjust these because with the old XML translation `<value/>` will sometimes result in a one element vector and other times a zero element vector.

## XML

### Tool Specification

<pre>
    &lt;tool name="<i><font color="red">tool-name</font></i>"&gt;
        <font color="blue">&lt;shortdescription&gt;</font><i>one-line description</i><font color="blue">&lt;/shortdescription&gt;</font>
        <font color="blue">&lt;description&gt;</font><i>paragraph description</i><font color="blue">&lt;/description&gt;</font>

        <font color="#7CFC00">&lt;needs&gt;</font><i><font color="red">tool-name</font></i><font color="#7CFC00">&lt;/needs&gt;</font>

        <font color="blue">&lt;code&gt;</font>
            <font color="blue">&lt;include&gt;</font>path to headerfile<font color="blue">&lt;/include&gt;</font>
            <font color="blue">&lt;private&gt;</font>
                <font color="blue">&lt;include&gt;</font>path to headerfile<font color="blue">&lt;/include&gt;</font>
            <font color="blue">&lt;/private&gt;</font>
        <font color="blue">&lt;/code&gt;</font>

        <font color="#7CFC00">&lt;method</font> name="<i><font color="red">method-name</font></i>"<font color="#7CFC00">&gt;</font>
            <font color="blue">&lt;shortdescription&gt;</font><i>one-line description</i><font color="blue">&lt;/shortdescription&gt;</font>
            <font color="blue">&lt;description&gt;</font><i>paragraph description</i><font color="blue">&lt;/description&gt;</font>
            <font color="blue">&lt;input&gt;</font>
                <font color="#7CFC00">&lt;param</font> name="<font color="red">param-name</font>"<font color="#7CFC00">&gt;</font>
                    <font color="blue">&lt;description&gt;</font><i>short description</i><font color="blue">&lt;/description&gt;</font>
                    <font color="#7CFC00">&lt;type</font> <font color="blue">units="<i><font color="red">unit-name</font></i>" mustexist="true"</font><font color="#7CFC00">&gt;</font><i><font color="red">type-name</font></i><font color="#7CFC00">&lt;/type&gt;</font>
                    <b><font color="blue">&lt;value&gt;value-content&lt;/value&gt;</font></b>
                <font color="#7CFC00">&lt;/param&gt;</font>
            <font color="blue">&lt;/input&gt;</font>
            <font color="blue">&lt;output&gt;</font>
                <font color="#7CFC00">&lt;param</font> name="<font color="red">param-name</font>"<font color="#7CFC00">&gt;</font>
                    <font color="blue">&lt;description&gt;</font><i>short description</i><font color="blue">&lt;/description&gt;</font>
                    <font color="#7CFC00">&lt;type</font> <font color="blue">units="<i><font color="red">unit-name</font></i>" mustexist="true"</font><font color="#7CFC00">&gt;</font><i><font color="red">type-name</font></i><font color="#7CFC00">&lt;/type&gt;</font>
                    <b><font color="blue">&lt;value&gt;value-content&lt;/value&gt;</font></b>
                <font color="#7CFC00">&lt;/param&gt;</font>
            <font color="blue">&lt;/output&gt;</font>

            <font color="blue">&lt;returns</font> type="<font color="red">type-name</font>"<font color="blue">/&gt;</font>
        <font color="#7CFC00">&lt;/method&gt;</font>
    &lt;/tool&gt;
</pre>

* <font color="red">Red</font> text indicate unique names. They may be developer specified (e.g. tool-name), they may be predefined enumerations (e.g. unit-name), or they may be both (e.g. type-name).
* <font color="blue">Blue</font> XML elements or attributes indicates optional items (zero or one).
* <font color="#7CFC00">Green</font> XML elements can be repeated.
* Black XML elements or attributes must be supplied (once)

## testing links

[more information for testing](test.md)

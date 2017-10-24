
## CASAtools

CASAtools is a self-contained python module that provides the tools from the [CASA](http://casa.nrao.edu/) project. The module only contains the non-GUI tools which are the [SWIG](http://swig.org) bound C++ functionality from CASA.

## Maturity

This is alpha software in every sense of the word. Currently the goal is just to make this usable on the Linux platforms which the [CASA](http://casa.nrao.edu) project supports for development, [RedHat Enterprise Linux](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) release 6 and release 7.

## Building CASAtools

#### Install Dependencies

The CASA development RPMs must be installed. This can be accomplished by adding the CASA [YUM](https://en.wikipedia.org/wiki/Yum_(.rpm)) repository to /etc/yum.repos.d and installing:

```
-bash-4.2# yum install ccache
-bash-4.2# yum install casa-toolset-2
```

#### Checkout

Checkout the CASAtools source code:

```
-bash-4.2$ git clone --recursive https://open-bitbucket.nrao.edu/scm/casa/CASAtools.git
```

#### Build

After the CASA build environment is installed, the CASAtools module can be built like:

```
-bash-4.2$ cd CASAtools
-bash-4.2$ scripts/gcw-pick
-bash-4.2$ aclocal
-bash-4.2$ autoconf
-bash-4.2$ ./configure
-bash-4.2$ ./setup.py build
```

#### See If It Works

If the build completes successfully, try loading the CASAtools module:

```
-bash-4.2$ type ipython
ipython is hashed (/opt/casa/02/bin/ipython)
-bash-4.2$ PYTHONPATH=build/lib.linux-x86_64-2.7 ipython 
Python 2.7.12 (default, Apr  4 2017, 16:53:53) 
Type "copyright", "credits" or "license" for more information.

IPython 5.1.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: from CASAtools import image

In [2]: ia = image( )

In [3]: ia.fromshape("mytest.im",[20,20])
2017-10-11 20:49:13	INFO		Created Paged image 'mytest.im' of shape [20, 20] with float valued pixels.
Out[3]: True

In [4]: exit
-bash-4.2$ 
```

#### Notes

If some time has passed since the last build, you should (sometimes) remove *xml-casa-assembly-1.0.jar*, e.g.
```
-bash-4.2$ rm ./scripts/java/xml-casa-assembly-1.0.jar
-bash-4.2$ scripts/gcw-pick
```
before rebuilding because this [JAR](https://en.wikipedia.org/wiki/JAR_(file_format)) file is automatically fetched from a download site by *gcw-pick*. However, it is not fetched if it already exists. Deleting the current copy will result in a new copy being fetched which *may* be newer.

## Available Tools

| Tool Name           | Description                                                     |
| ------------------- | --------------------------------------------------------------- |
| atmosphere          | Atmosphere model                                                |
| calanalysis         | Get and fit data from a calibration table (CASA 3.4 and later). |
| calibrater          | Synthesis calibration (self- and cross-)                        |
| componentlist       | A tool for the manipulation of groups of components             |
| coordsys            | Operations on CoordinateSystems                                 |
| functional          | Functionals handling                                            |
| imagepol            | Polarimetric analysis of images                                 |
| image               | Operations on images                                            |
| measures            | measures tool                                                   |
| msmetadata          | Operations to retrieve metadata from a measurment set           |
| ms                  | Operations on measurement sets                                  |
| quanta              | quanta tool handles units and quantities                        |
| regionmanager       | Create and manipulate regions of interest                       |
| synthesisimager     | tool for synthesis imaging                                      |
| synthesisimstore    | tool for synthesis imaging                                      |
| synthesisnormalizer | tool for synthesis imaging                                      |
| synthesisutils      | tool for synthesis imaging                                      |
| table               | Access tables from casapy                                       |


## Changes from Standard CASA

While the goal was to simply reconstitute the [CASA tools](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa/browse/gcwrap/tools) within a unencumbered python module, deviations were required as work progressed. These deviations are divided into categories based upon whether the change relates to the way the tools behave or the way the XML files are structured. The *gcw-pick* script includes an automatic XML update function which is able to make some of the more basic changes to the XML code automatically as the XML files are pulled from the *casa-source/gcwrap* tree into CASAtools.

### Changes to Behavior

1. __Parameter Order__ --- with standard CASA, the [SWIG](http://swig.org) binding is exposed directly to the user. While this is often fine, for CASAtools more control was desired over how parameter checking happens and how coersion happens. Because of this, a wrapper is created over the SWIG binding. This wrapper explicitly specifies the order of the parameters. The input parameters come first, followed by the output parameters. Within the input and output parameters, the order is determined by the XML file. With standard CASA, all parameters are largely unordered with binding controlled by the way [python argument lists](https://stackoverflow.com/questions/3394835/args-and-kwargs) are handled.

2. __Utils Tool__ --- the standard CASA *uttool* (AKA *utils* or just *ut*), has been converted from a toolbox/grab bag **tool** to a **singleton object** called *ctsys* in CASAtools, e.g.:
    ```
    In [1]: from CASAtools import ctsys
    ```

3. __ctsys.resolve( )__ --- a new member function was added to resolve the path to an data file based upon **CASADATA** path (as is done for `<type mustexist="true">path</type>`)

### XML Changes

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
    &lt;tool&gt;
        &lt;method name="<i><font color="red">tool-name</font></i>"&gt;
            <font color="blue">&lt;needs&gt;<i>tool-name-1</i>&lt;needs/&gt;</font>
            <font color="LightBlue">&lt;needs&gt;<i>tool-name-2</i>&lt;needs/&gt;</font>
        &lt;method/&gt;
    &lt;tool/&gt;
</pre>

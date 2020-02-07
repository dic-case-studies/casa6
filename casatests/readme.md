## CASAtests

CASAtests is a self-contained python module that provides special test scripts for CASA. Many of the tests in this package depend on the [casatestutils](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse) python module being pre-installed in your system.

### Install dependencies

First make sure you have installed the CASAtools, the CASAtasks and the CASAtestutils before you run the tests. Below is a quick guide on how to install
them with pip.
```
-bash-4.2# python3 -m venv casa6_0
-bash-4.2# source casa6_0/bin/activate
-bash-4.2# pip3 install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools
-bash-4.2# pip3 install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks
-bash-4.2# pip3 install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatestutils
```

### Data sets for the tests

Download the datasets from https://open-bitbucket.nrao.edu/scm/casa/casa-data-req.git

Add to $HOME/.casa/toolrc.py the following line, pointing to the locations of your datasets
```
-bash-4.2# cat $HOME/.casa/toolrc.py
datapath=[ "/your-data-dir/casa-data/" , "/your-data-dir/casa-data-req/","/your-data-dir/mydata/"]
```
### Run tests

All tests added to casatests must be able to run with Python.

#### Run a single test using Python

The following performance test will look for a dataset located in <your_data_dir>/visibilities/alma/.
Make sure your toolrc.py datapath has an entry pointing to <your_data_dir>. 

```
-bash-4.2# python3 ./casatests/performance/test_perf_tclean_mem_setweighting.py
```

#### Run a test using pytest (it needs to be pre-installed with pip3)
```
-bash-4.2# pytest ./casatests/performance/test_perf_tclean_mem_setweighting.py 
```

#### Run a single test using a casalith tarball
First, install the casatestutils wheel using the pip3 from the tarball and run the test with python3 or the casa binary. Example:
```
-bash-4.2# ./casa-6.1.0/bin/pip3 install https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatestutils
-bash-4.2# ./casa-6.1.0/bin/python3 ./casatests/performance/test_perf_tclean_mem_setweighting.py
-bash-4.2# ./casa-6.1.0/bin/casa -c ./casatests/performance/test_perf_tclean_mem_setweighting.py
```

#### Run all stakeholders test
......

### Add a new test to casatests

The following is the naming convention for tests added to the directories inside casatests. 

#### performance
   * test\_perf\_[taskname]\_[description]\_[...].py
   * e.g. test\_perf\_tclean\_mem\_setweighting.py
   * e.g. test\_perf\_tclean\_runtime\_cube\_model\_write.py   

#### stakeholders
   * test\_stk\_[stakeholder]\_[taskname]\_[description].py
   * e.g. test\_stk\_alma\_pipeline\_tclean\_cycle8.py

#### e2e: End-to-End
   * test\_e2e\_[description]\_[...].py
   * e.g. test\_e2e\_alma\_m100\_if.py
   * e.g. test\_e2e\_alma\_m100\_sd.py

#### pipeline
   * test\_pipe\_[telescope]\_[use-case]\_[description]\_[...].py
   * e.g. test\_pipe\_vlass\_calib\_[...].py
   * e.g. test\_pipe\_alma\_full\_ephem.py

#### benchmarks
   * test\_bench\_[description]\_[...].py
   




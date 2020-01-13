## CASAtests

CASAtests is a self-contained python module that provides special test scripts for CASA. Many of the tests in this package depend on the [casatestutils](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse) python module being pre-installed in your system.

#### Install Dependencies

First make sure you have installed the CASAtasks and/or the CASAtools before you run the tests. Below is a quick guide on how to install
them from the Python wheels.
```
-bash-4.2# python3 -m venv casa6_0
-bash-4.2# source casa6_0/bin/activate
-bash-4.2# pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools
-bash-4.2# pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks
-bash-4.2# pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatestutils
```

#### Data sets for the tests

Download from ...

Add to $HOME/.casa/toolrc.py the following line, pointing to the locations of your datasets
```
-bash-4.2# cat $HOME/.casa/toolrc.py
datapath=[ "/casadata/user/casa-data/" , "/casadata/user/casa-data-req/","/casadata/user/mydata/"]
```

#### Run a single test using Python

This performance test will need a dataset located in <your_data_dir>/performance/.
Make sure your toolrc.py datapath has an entry pointing to <your_data_dir>. The script expects to find
a directory called casa-perf inside your_data_dir, where it will look for the files needed
to run the tests.


```
-bash-4.2# python3 casatests/performance/test_perf_tclean_mem_setweighting.py
```

#### Run a test using pytest
```
-bash-4.2# pytest casatests/performance/test_perf_tclean_mem_setweighting.py 
```

#### Run all stakeholders tests

......







## CASAtests

CASAtests is a self-contained python module that provides special test scripts for CASA. Many of the tests in this package depend on the [casatestutils](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse) python module being pre-installed in your system.

### Install Dependencies

First make sure you have installed the CASAtasks and/or the CASAtools before you run the tests. Below is a quick guide on how to install
them from the Python wheels.
```
-bash-4.2# python3 -m venv casa6_0
-bash-4.2# source casa6_0/bin/activate
-bash-4.2# pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools
-bash-4.2# pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks
-bash-4.2# pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatestutils
```

### Data sets for the tests

Download from ...

Add to $HOME/.casa/toolrc.py the following line, pointing to the locations of your datasets
```
-bash-4.2# cat $HOME/.casa/toolrc.py
datapath=[ "/casadata/user/casa-data/" , "/casadata/user/casa-data-req/","/casadata/user/mydata/"]
```
### Run tests

All tests added to casatests must be able to run with Python.

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

### Add a new test to casatests

The following is the naming convention for tests added to the directories inside casatests. 

##### performance
   * test_perf_[taskname]_[description]_[...].py
   * e.g. test_perf_tclean_memory_setweighting.py
   * e.g. test_perf_tclean_runtime_cube_model_write.py   

##### stakeholders
   * test_stk_[taskname]_[description]_[...].py
   * test_stk_tclean_alma_pipeline.py

##### e2e: End-to-End
   * test_e2e_[description]_[...].py
   * test_e2e_alma_m100_if.py
   * test_e2e_alma_m100_sd.py

#### pipeline
   * test_pipe_[telescope]_[use-case]_[description]_[...].py
   * test_pipe_vlass_calib_<...>.py
   * test_pipe_alma_full_ephem.py

#### benchmarks
   * test_bench_[description]_[...].py
   




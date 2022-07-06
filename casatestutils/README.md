# casatestutils Package

##casatestutils 
A generic testhelper module for use with CASA testing.

Documentation: https://open-confluence.nrao.edu/display/CASA/casatestutils%3A+A+generic+test+helper+module

## runtest.py
[runtest.py](casatestutils/runtest.py) is a single test wrapper to run CASA Python tests. The script can run one or
more scripts and will create a directory called ```nosedir``` in the working directory and inside nosedir it will create
a separated directory with the test name to hold any files created by each test script. 

With the latest version, runtest.py will also create a file called ```short_summary.log``` inside each 
test directory with a single line summary of the exit status of each test of the script. If more than one test 
script is run at the same time, runtest.py will create a file called ```summary_of_failed.log``` inside the top 
directory nosedir. This file will give the error messages of a ny test cases that failed or was skipped. In case all
test cases are successful, the file will not be created.

The script uses unittest and pytest and has the following command line options.

### Get help (`runtest.py -h or --help`)
```
usage: runtest.py [-h] [-i] [-v] [-x] [-s test [test ...]] [-f [FILE]]
                  [-e [MAPFILE]] [-b BRANCH] [-p PKG] [-w WORK_DIR]
                  [-n NCORES] [-t TEST_PATHS] [-l TEST_LIST] [-c TEST_CONFIG]
                  [-j TEST_GROUP] [-m PMODE] [--bamboo] [-r RCDIR]
                  [--ignore_list IGNORE_LIST]
```
Execute it with a casalith tarball or python
```
./casa.6.4.0.16/bin/casa -c ./runtest.py <path-to>/test_task_tclean.py

./python3 ./runtest.py <path-to>/test_task_tclean.py
```

### Find test scripts
**runtest.py** can run a test script from any location. If given a test script name with ending ".py",
it will get it from the given location. If ".py" is ommited, it will pull it from the CASA git
repository using the trunk branch.

### Display the list of existing tests (`runtest.py -i or --list`)
Run the -i or --list option to see the available test scripts along with the CASA JIRA
components associated with them. The JIRA components are defined in the file ```component_to_test_map.json```.
stored in casatestutils. See next section.
```
python3 runtest.py -i
python3 runtest.py --list
...
test_agentflagger ['agentflagger', 'casatools', 'default']
test_asdmsummary ['asdmsummary']
test_bandpass ['bandpass', 'Calibration']
...
```
### Map test scripts to CASA JIRA components (`runtest.py -j or --TEST_GROUP`)

Execute the tests associated to any CASA JIRA component. The mapping of test scripts
to JIRA components is defined in a JSON file available in [component_to_test_map.json](casatestutils/component_to_test_map.json).
For example, the test scripts associated with the JIRA component Flagging will be 
checked out from the casa6 git trunk and executed when using this feature.
```
python3 runtest.py -j Flagging
python3 runtest.py --TEST_GROUP Flagging
...

Namespace(bamboo=False, branch=None, classes=None, dry_run=True, file=None, list=False, mapfile=None, ncores=2, pkg=None, pmode=None, rcdir=None, test_config=None, test_group='Flagging', test_list=None, test_paths=None, verbose=False, work_dir=None)
Testing Components['Flagging']

Testnames: ['test_flagcmd', 'test_flagdata', 'test_flagmanager']
Cleaning: /opt/casa/Tests/verification/CAS-13640/py38/nosedir/
Tests: ['test_flagcmd', 'test_flagdata', 'test_flagmanager']
Setting Working Directory: /opt/casa/Tests/verification/CAS-13640/py38/nosedir/test_flagcmd/
Fetching Tests From Git Main Since No Local Test is Given
```

### Dry-run the tests (`runtest.py -x or --dry-run`)
```
python3 runtest.py -x <path-to>/test_req_task_listobs.py
python3 runtest.py --dry-run <path-to>/test_req_task_listobs.py
...
collected 88 items                                                                                                                                                                           
<Module test_req_task_listobs.py>
================================================== 7 warnings in 10.65s ==========================================
```
### Run a test script
Run from a local test script, from git trunk or from a JIRA branch. When **runtest.py**
checks out the test script from git trunk or a branch, the test script should not be given
with the ".py" extension.

#### from a local test script
```
python3 runtest.py /path-to/test_mytask.py
```
#### from the CASA git trunk
```
python3 runtest.py test_task_sdpolaverage
...
Testnames: ['test_task_sdpolaverage']
Fetching Tests From Git Main Since No Local Test is Given
...
```
#### from a CASA JIRA branch (`runtest.py -b or --branch`)
```
python3 runtest.py -b CAS-13640 test_tclean[test_onefield_clark]
...
Testnames: ['test_task_tclean[test_onefield_clark]']

CHECKING OUT BRANCH: CAS-13640
...
test_tclean.py::test_onefield::test_onefield_clark 
[onefield] Test_Onefield_clark : mfs with clark minor cycle  
...
```
#### Run a test using the CASA executable
```
./casa.6.4.0.16/bin/casa -c ./runtest.py <path-to>/test_task_tclean.py
```
#### Run specific test cases from a local script
This is useful when knowing a priori the name of the test case or class to run.
```
./python3 -c ./runtest.py <path-to>/test_task_tclean.py[test_onefield_clark,test_onefield_mem]
```

#### Run a test inside a casashell
```
python3 -m casashell
> from casatestutils import runtest
> runtest.run(['test_task_tclean'])                       # pull test script from git trunk
> runtest.run(['/path-to/test_task_flagdata.py'])         # run a local test script
> runtest.run(['test_ttask_clean[test_onefield_clark]'])  # pull test script from git trunk
```

#### Run tests from a JIRA component and ignore some tests (`runtest.py -j <component> --ignore_list`)
This option is used to ignore a test from a test suite when using the -j option.
The parameter can be a list of comma separated tests or a JSON file in the same structure as component_to_test_map.py.
In this example, the component Flagging will pull 3 test scripts for flagcmd, flagdata and flagmanager.
and the parameters are asking to ignore test_flagdata and test_flagcmd. Only test_flagmanager.py will run in this example.
```
./python3 ./runtest.py -j Flagging --ignore_list [test_task_flagdata,test_task_flagcmd]
```

#### Run a test similar to Bamboo setup
Experimental for Developers, mainly for Test Infrastructure Team use.

```
Required Flags
-m, --pmode     : Parallelization mode: serial, parallel, both
-w, --work_dir  : Path to Working Directory to Unpack Tar/Dmg and run tests
-p, --pkg       : Tarball / Dmg to be tested
--bamboo        : Flag to Tell runtest similar to Bamboo
-j / --test_group  or -l / --test_list : A test group ( Component) or test list must be passed 

Optional Flags
-n, --ncores    : Number of Cores to Use for MPI Tests ( Default to 2)
-r, --rcdir     : Casa rcdir 
##### Examples
python3 runtest.py --bamboo -n 4 -p casa-6.4.3-3-py3.6.tar.xz -m serial -j asdmsummary -w /path/to/working/directory
python3 runtest.py --bamboo -p casa-6.4.3-3-py3.6.tar.xz -m serial  -w /path/to/working/directory --test_list test_coordsys,test_tclean,

```

#### Additional Notes

--bamboo option requires [testrunner module] (https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse/casatestutils/testrunner). 

#### Verify runtest.py
Run the unit tests of runtest.py in the command-line or in a Jupyter notebook

Run locally using a casa tarball:
```
<tarball>/bin/pip3 install jupyter

<tarball>/bin/python3 -m jupyter notebook --generate-config

<tarball>/bin/python3 -m jupyter notebook --browser=firefox --ip='*' --NotebookApp.toke='' --NotebookApp.password='' tests/nb_test_runtest.ipynb
```
Run in [Google Colab](https://colab.research.google.com/drive/1lunhY-8iLot2H0UwFJ98_IWsmBM2TIMd?usp=sharing)

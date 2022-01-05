# casatestutils Package

##casatestutils 
A generic testhelper module for use with CASA testing.

Documentation: https://open-confluence.nrao.edu/display/CASA/casatestutils%3A+A+generic+test+helper+module

## runtest.py
[runtest.py](casatestutils/runtest.py) is a single test wrapper to run CASA Python tests. 
The script uses unittest and pytest and has the following command line options.

### Get help (`runtest.py -h or --help`)
```
usage: runtest.py [-h] [-i] [-v] [-x] [-s test [test ...]] [-f [FILE]]
                  [-e [MAPFILE]] [-b BRANCH] [-p PKG] [-w WORK_DIR]
                  [-n NCORES] [-t TEST_PATHS] [-l TEST_LIST] [-c TEST_CONFIG]
                  [-j TEST_GROUP] [-m PMODE] [--bamboo]
```
Execute it with a casalith tarball or python
```
./casa.6.4.0.16/bin/casa -c ./runtest.py <path-to>/test_tclean.py

./python3 ./runtest.py <path-to>/test_tclean.py
```

### Find test scripts
runtest.py can run a test script from any location. If given a test script name with ending ".py",
it will get it from the given location. If ".py" is ommited, it will pull it from the CASA git
repository using the trunk branch.

### Display the list of existing tests (`runtest.py -i or --list`)
Run the -i or --list option to see the available test scripts along with the CASA JIRA
components associated with them.
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
Run from a local test script, from git trunk or from a JIRA branch

#### from a local test script
```
python3 runtest.py /path-to/test_mytask.py
```
#### from the CASA git trunk
```
python3 runtest.py test_sdpolaverage
...
Testnames: ['test_sdpolaverage']
Fetching Tests From Git Main Since No Local Test is Given
...
```
#### from a CASA JIRA branch (`runtest.py -b or --branch`)
```
python3 runtest.py -b CAS-13640 test_tclean[test_onefield_clark]
...
Testnames: ['test_tclean[test_onefield_clark]']

CHECKING OUT BRANCH: CAS-13640
...
test_tclean.py::test_onefield::test_onefield_clark 
[onefield] Test_Onefield_clark : mfs with clark minor cycle  
...
```
#### run a test using the CASA executable
```
./casa.6.4.0.16/bin/casa -c ./runtest.py <path-to>/test_tclean.py
```
#### run a test inside a casashell
```
python3 -m casashell
> from casatestutils import runtest
> runtest.run['test_tclean']  # pull test script from git trunk
> runtest.run['/path-to/test_flagdata.py']  # run a local test script
> runtest.run['test_tclean[test_onefield_clark]']  # pull test script from git trunk
```
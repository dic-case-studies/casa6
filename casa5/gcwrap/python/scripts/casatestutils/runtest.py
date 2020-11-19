from __future__ import print_function
import argparse
import os
import shutil
import logging
import sys
import traceback
import subprocess
import unittest
import json
sys.path.insert(0,'')


####################            Constants

# Python Version
PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

# CASA Directory
CASA_DIR = os.environ["CASAPATH"].split()[0]

# mem mode variables
HAVE_MEMTEST=True
MEM = 0
try:
    import memTest
except ImportError:
    HAVE_MEMTEST = False

# cov mode variables
HAVE_COVTEST=True
COV = 0
try:
    import coverage
except ImportError:
    HAVE_COVTEST = False

# pybot mode variables
HAVE_ROBOT = True
USE_PYBOT = 0
try:
    import robot
except ImportError:
    HAVE_ROBOT = False

#### PYTEST IMPORT
HAVE_PYTEST = True
try:
    import pytest
except ImportError:
    HAVE_PYTEST = False

#### NOSE IMPORT
HAVE_NOSE = True
try:
    import nose
except ImportError:
    HAVE_NOSE = False

IS_CASA6 = False
CASA6 = False
HAVE_CASA6 = False

# JIRA BRANCH TO CHECKOUT
JIRA_BRANCH = None

try:
    import casatools
    import casatasks
    CASA6 = True
    IS_CASA6 = True
    HAVE_CASA6 = True
    
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

if not IS_CASA6:
    import testwrapper
    from testwrapper import UnitTest
    

# Use Nose attribute Functionality
RUN_SUBTEST = False

# Dry run of Tests
DRY_RUN = False

# Define which tests to run
whichtests = 0

######
def main(testnames, verbose = True, DRY_RUN = False):
    if IS_CASA6:
        #sys.path.append(os.path.abspath(os.path.basename(__file__)))
        if HAVE_PYTEST:
            cwd = os.getcwd() + "/"
            workpath = os.getcwd() +"/nosedir/"
            workdir = os.getcwd() +"/nosedir/"

            clean_working_directory(workpath)
            # Copy Tests to Working Directory
            os.makedirs(workdir)
            print("Tests: {}".format(testnames))
            gittest = True

            
            for testname in testnames:
                nonlocalTest = False
                cmd = []
                if testname.startswith("../") or testname.startswith("/"):
                    try:
                        real_path = os.path.realpath(testname)
                        testname = testname.split("/")[-1]
                        nonlocalTest = True
                    except:
                        traceback.print_exc()
                        
                # Copy Test To nosedir Directory if in cwd
                if testname.startswith("test"):
                    test = testname
                    # Check if specific tests are requested
                    if "[" and "]" in test:
                        testname = getname(test)
                        tests = gettests(test)
                        #print(testname)
                        #print(tests)

                        teststring = ""
                        #print(len(tests))
                        if len(tests) == 1:
                            teststring = tests[0]
                        elif len(tests) > 1:
                            print(tests)
                            teststring = " or ".join(tests)
                        #workdir = os.getcwd() +"/workdir/nosedir/{}/".format(testname)
                        #os.makedirs(workdir)
                        #cmd = [ workdir ]
                        #cmd = [ ".".join([testname,'py']), "-k {}".format(teststring)]
                        cmd = ["-k {}".format(teststring)] + cmd
                        test = testname
                        
                    # Set up Test Working Directory
                    if not os.path.exists(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])):
                        print("Setting Working Directory: {}".format(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                        os.makedirs(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                        cmd = [ workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]) ] + cmd
                        

                        
                    # Check to see if tests need to be pulled from git. Only needs to be done once
                    if not test.endswith(".py") and gittest == True:
                        print("Fetching Tests From Git Since No Local Test is Given")
                        git_fetch_casa_tests( workpath + 'casa6')
                        os.makedirs(workdir + "tests/")
                        gather_all_tests(workpath +'casa6/', workdir + "tests/")
                        gittest = False
                        
                    if test.endswith(".py"):
                        
                        try:
                            if nonlocalTest:
                                if "[" in real_path:
                                    real_path = real_path.split("[")[0]
                                print("Copying: {} to {}".format(real_path, workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                                shutil.copy2(real_path, workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                            else:
                                print("Copying: {} to {}".format(test, workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                                shutil.copy2(test, workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                        except:
                            traceback.print_exc()
                    else:
                        try:
                            print("Copying: {} to {}".format(workdir + "tests/",test), workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                            shutil.copy2("{}{}.py".format(workdir + "tests/",test), workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                        except:
                            traceback.print_exc()
                            
                    if verbose:
                        cmd = ["--verbose"] + cmd


                    if DRY_RUN:
                        cmd = ["--collect-only"] + cmd
                    
                    if not os.path.isdir(workpath + '/xml/{}/'.format(test if not test.endswith(".py") else test[:-3])):
                        os.makedirs(workpath + '/xml/{}/'.format(test if not test.endswith(".py") else test[:-3]))
                    xmlfile = workpath + 'xml/{}/nose.xml'.format(test if not test.endswith(".py") else test[:-3])

                    #########
                    #============================================================ warnings summary =====================================================
                    # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436
                    # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436:
                    # PytestDeprecationWarning: The 'junit_family' # default value will change to 'xunit2' in pytest 6.0.
                    # Add 'junit_family=xunit1' to your pytest.ini file to keep the current format in future versions of pytest and silence thiswarning.
                    #  _issue_warning_captured(deprecated.JUNIT_XML_DEFAULT_FAMILY, config.hook, 2)
                    #########
                    cmd = ["--junitxml={}".format(xmlfile)] + ["-s"] + cmd
                    #print("Running Command: pytest {}".format(cmd))
                    #print("Work Path: {}".format(workpath))
                    if len(os.listdir(workpath)) < 1: # If only the XML dir was created
                        print("No Tests to Run")
                        sys.exit()
                    else:
                        print("Running Command: pytest {}".format(cmd))
                        pytest.main(cmd)
                        
def readfile(FILE):
    # It will skip lines that contain '#' and
    # it will only read words starting with test
    if(not os.path.exists(FILE)):
        print('ERROR: List of tests does not exist')
        return []

    List = []
    infile = open(FILE, "r")
    for newline in infile:
        if newline.__contains__('#'):
            continue
        
        if newline.startswith('test'):
            words = newline.split()
            List.append(words[0])

    infile.close()
    return List

def list_tests():
    print('Full list of unit tests')
    print('-----------------------')
    if IS_CASA6:
        git_fetch_casa_tests(os.getcwd() +"/testlist/casa6")
        gather_all_tests(os.getcwd() +"/testlist/casa6",os.getcwd() +"/testlist/")
        tests = sorted(os.listdir(os.getcwd() +"/testlist/"))
        for test in tests:
            if test.startswith("test"):
                print(test)
    else:
        for t in readfile(LISTofTESTS):
            print(t)

def git_fetch_casa_tests_branch(path, branch):

    cwd = os.getcwd()
    if os.path.exists(path):
        try:
            os.rmdir(path)
        except:
            shutil.rmtree(path)
    os.makedirs(path)

    os.chdir(path)

    FNULL = open(os.devnull, 'w')
    print("CHECKING OUT BRANCH: {}".format(branch))
    subprocess.call(["git","clone", "https://open-bitbucket.nrao.edu/scm/casa/casa6.git"], stdout=FNULL, stderr=subprocess.STDOUT)
    os.chdir(path + 'casa6')
    subprocess.call(["git","checkout","{}".format(branch)], stdout=FNULL, stderr=subprocess.STDOUT)
    shutil.move(path + 'casa6/casatools', path + 'casatools')
    shutil.move(path + 'casa6/casatasks', path + 'casatasks')
    shutil.move(path + 'casa6/casatestutils', path + 'casatestutils')
    shutil.move(path + 'casa6/casatests', path + 'casatests')
    shutil.move(path + 'casa6/casa5', path + 'casa5')
    shutil.rmtree(path + 'casa6')

    os.chdir(cwd)


def gather_all_tests(path, workpath):

    if sys.version_info[0] > 2:
        import pathlib
        for filename in pathlib.Path(path).rglob('test_*.py'):
            #print(filename)
            shutil.copy2(filename, workpath)

def gather_single_tests(path, workpath,test):

    foundTest = False
    print("Searching for Test: {}.py".format(test))
    if sys.version_info[0] > 2:
        import pathlib
        for filename in pathlib.Path(path).rglob('{}.py'.format(test)):
            foundTest = True
            shutil.copy2(filename, workpath)
    
    if not foundTest:
        print("Could Not Find Test: {}.py".format(test))

def print_test_classes(testnames):
    for test in testnames:
        if verbose:
            pytest.main(["--collect-only", "-q", "-v", test])
        else:
            pytest.main(["--collect-only", test])

def clean_working_directory(workdir):

    print("Cleaning: {}".format(workdir))
    if os.path.exists(workdir):
        try:
            os.rmdir(workdir)
        except:
            shutil.rmtree(workdir)
    

def getname(testfile):
    '''Get the test name from the command-line
       Ex: from test_clean[test1], returns test_clean'''
    n0 = testfile.rfind('[')
    n1 = testfile.rfind(']')
    if n0 != -1:
        return testfile[:n0]

def gettests(testfile):
    '''Get the list of specific tests from the command-line
       Ex: from test_clean[test1,test3] returns [test1,test3]'''
    n0 = testfile.rfind('[')
    n1 = testfile.rfind(']')
    if n0 != -1:
        temp = testfile[n0+1:n1]
        tests = temp.split(',')
        return tests

def getclasses(tests):
    if HAVE_PYTEST:
        for test in tests:
            if test.endswith(".py"):
                pytest.main(["--fixtures",test])
    else:
        logger.debug("Start def getclasses()")
        '''Get the classes of a test script It will copy the test script to /tmp and remove it afterwards'''
        here = os.getcwd()
        tmpdir = '/tmp'
        for filename in testnames:
            if not filename.startswith("test_"):
                #print "Cannot Get Classes for Regression Test"
                logger.error("Cannot Get Classes for Regression Test: %s",filename)
                return
        try:
            os.chdir(tmpdir)
            
            for filename in testnames:
                tt = UnitTest(filename)
                tt.copyTest(copyto=tmpdir)
                classes = tt.getTestClasses(filename)
                for c in classes:
                    pprint.pprint('Class '+c.__name__)
                    for attr, value in c.__dict__.iteritems():
                        if len(attr) >= len("test") and attr[:len("test")] == "test":
                            print('\t{}'.format(c(attr)))
                os.remove(filename+'.py')
                os.remove(filename+'.pyc')
            os.chdir(here)
        except:
            print('--> ERROR: Cannot copy script to {}'.format(tmpdir))
            logger.error('Failed to open file', exc_info=True)
            return

def haslist(name):
    '''Check if specific list of tests have been requested'''
    n0 = name.rfind('[')
    n1 = name.rfind(']')
    if n0 == -1:
        return False
    return True

def run(testnames):

    if IS_CASA6:
        #sys.path.append(os.path.abspath(os.path.basename(__file__)))
        if HAVE_PYTEST:
            cwd = os.getcwd() + "/"
            workpath = os.getcwd() +"/nosedir/"
            workdir = os.getcwd() +"/nosedir/"

            clean_working_directory(workpath)
            # Copy Tests to Working Directory
            os.makedirs(workdir)

            if RUN_ALL:
                print("Gathering All Tests")
                git_fetch_casa_tests( workpath + 'casa6')
                os.makedirs(workdir + "all/")
                gather_all_tests(workpath +'casa6/', workdir + "all/")
                cmd = [ workdir ]
                cmd = ["--continue-on-collection-errors"] + cmd
                if verbose:
                    cmd = ["--verbose"] + cmd


                if DRY_RUN:
                    cmd = ["--collect-only"] + cmd
                    
                if not os.path.isdir(workpath + '/xml/all/'):
                    os.makedirs(workpath +  '/xml/all/')
                xmlfile = workpath + '/xml/all/nose.xml'

                #########
                #============================================================ warnings summary =====================================================
                # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436
                # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436:
                # PytestDeprecationWarning: The 'junit_family' # default value will change to 'xunit2' in pytest 6.0.
                # Add 'junit_family=xunit1' to your pytest.ini file to keep the current format in future versions of pytest and silence this warning.
                #  _issue_warning_captured(deprecated.JUNIT_XML_DEFAULT_FAMILY, config.hook, 2)
                #########
                cmd = ["--junitxml={}".format(xmlfile)] + ["-s"] + cmd

                if len(os.listdir(workdir + "all/")) == 0:
                    print("No Tests to Run")
                else:
                    print("Running Command: pytest {}".format(cmd))
                    pytest.main(cmd)
                
            else:
                print("Tests: {}".format(testnames))
                gittest = True

                for testname in testnames:
                    cmd = []
                    
                    # Copy Test To nosedir Directory if in cwd
                    if testname.startswith("test"):
                        test = testname
                        # Check if specific tests are requested
                        if "[" and "]" in test:
                            testname = getname(test)
                            tests = gettests(test)
                            #print(testname)
                            #print(tests)

                            teststring = ""
                            #print(len(tests))
                            if len(tests) == 1:
                                teststring = tests[0]
                            elif len(tests) > 1:
                                print(tests)
                                teststring = " or ".join(tests)
                            #workdir = os.getcwd() +"/workdir/nosedir/{}/".format(testname)
                            #os.makedirs(workdir)
                            #cmd = [ workdir ]
                            #cmd = [ ".".join([testname,'py']), "-k {}".format(teststring)]
                            cmd = ["-k {}".format(teststring)] + cmd
                            test = testname
                            
                        # Set up Test Working Directory
                        if not os.path.exists(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])):
                            print("Setting Working Directory: {}".format(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                            os.makedirs(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                            cmd = [ workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]) ] + cmd
                            

                            
                        # Check to see if tests need to be pulled from git. Only needs to be done once
                        if not test.endswith(".py") and gittest == True:
                            if JIRA_BRANCH is not None:
                                git_fetch_casa_tests_branch(workpath + 'casa6/', JIRA_BRANCH)
                                os.makedirs(workdir + "tests/")
                                gather_all_tests(workpath +'casa6/', workdir + "tests/")
                                gittest = False

                            else:
                                print("Fetching Tests From Git Main Since No Local Test is Given")
                                git_fetch_casa_tests( workpath + 'casa6/')
                                os.makedirs(workdir + "tests/")
                                gather_all_tests(workpath +'casa6/', workdir + "tests/")
                                gittest = False
                            
                        if test.endswith(".py"):
                            try:
                                print("Copying: {} to {}".format(test, workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                                shutil.copy2(test, workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                            except:
                                traceback.print_exc()
                        else:
                            try:
                                print("Copying: {} to {}".format(workdir + "tests/",test), workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                                shutil.copy2("{}{}.py".format(workdir + "tests/",test), workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                            except:
                                traceback.print_exc()
                                
                        # https://docs.pytest.org/en/stable/usage.html
                        if verbose:
                            cmd = ["--verbose"] + ["--tb=short"] + cmd
                        elif not verbose:
                            cmd = ["-ra"] + ["--tb=line"] + cmd

                        if DRY_RUN:
                            cmd = ["--collect-only"] + cmd
                        
                        if not os.path.isdir(workpath + '/xml/{}/'.format(test if not test.endswith(".py") else test[:-3])):
                            os.makedirs(workpath + '/xml/{}/'.format(test if not test.endswith(".py") else test[:-3]))
                        xmlfile = workpath + 'xml/{}/nose.xml'.format(test if not test.endswith(".py") else test[:-3])

                        #########
                        #============================================================ warnings summary =====================================================
                        # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436
                        # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436:
                        # PytestDeprecationWarning: The 'junit_family' # default value will change to 'xunit2' in pytest 6.0.
                        # Add 'junit_family=xunit1' to your pytest.ini file to keep the current format in future versions of pytest and silence thiswarning.
                        #  _issue_warning_captured(deprecated.JUNIT_XML_DEFAULT_FAMILY, config.hook, 2)
                        #########
                        cmd = ["--junitxml={}".format(xmlfile)] + ["-s"] + ["--disable-pytest-warnings"] + cmd
                        #print("Running Command: pytest {}".format(cmd))
                        #print("Work Path: {}".format(workpath))
                        if len(os.listdir(workpath)) < 1: # If only the XML dir was created
                            print("No Tests to Run")
                            sys.exit()
                        else:

                            myworkdir = os.getcwd()
                            os.chdir("{}".format(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                            print("Test Directory: {}".format(os.getcwd()))
                            print("Running Command: pytest {}".format(cmd))
                            pytest.main(cmd)
                            os.chdir(myworkdir)
                        
                    ##################################################
                    ########## Real Path ##########
                    ##################################################
                     # Copy Test To nosedir Directory assuming it's in another location
                    elif testname.startswith("/"):
                        testpath = testname.split("[")[0]
                        cmd = []
                        dirname = testname.split("/")[-1]
                        test = dirname
                        if "[" and "]" in test:
                            testname = getname(test)
                            tests = gettests(test)
                            #print(testname)
                            #print(tests)

                            teststring = ""
                            #print(len(tests))
                            if len(tests) == 1:
                                teststring = tests[0]
                            elif len(tests) > 1:
                                print(tests)
                                teststring = " or ".join(tests)
                            #workdir = os.getcwd() +"/workdir/nosedir/{}/".format(testname)
                            #os.makedirs(workdir)
                            #cmd = [ workdir ]
                            #cmd = [ ".".join([testname,'py']), "-k {}".format(teststring)]
                            cmd = ["-k {}".format(teststring)] + cmd
                            dirname = testname
                        
                        dirname = "{}".format(dirname if not dirname.endswith(".py") else dirname[:-3])

                        # Set up Test Working Directory
                        if not os.path.exists(workdir + "{}/".format(dirname)):
                            print("Setting Working Directory: {}".format(workdir + "{}/".format(dirname)))
                            os.makedirs(workdir + "{}/".format(dirname))
                            cmd = [ workdir + "{}/".format(dirname) ] + cmd
                        try:
                            shutil.copy2(testpath, workdir + "{}/".format(dirname))
                        except:
                            traceback.print_exc()
     
                        if verbose:
                            cmd = ["--verbose"] + ["--tb=short"] + cmd
                        elif not verbose:
                            cmd = ["-ra"] + ["--tb=line"] + cmd


                        if DRY_RUN:
                            cmd = ["--collect-only"] + cmd
                        
                        if not os.path.isdir(workpath + '/xml/{}/'.format(dirname)):
                            os.makedirs(workpath + '/xml/{}/'.format(dirname))
                        xmlfile = workpath + 'xml/{}/nose.xml'.format(dirname)

                        #########
                        #============================================================ warnings summary =====================================================
                        # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436
                        # ./lib/py/lib/python3.6/site-packages/_pytest/junitxml.py:436:
                        # PytestDeprecationWarning: The 'junit_family' # default value will change to 'xunit2' in pytest 6.0.
                        # Add 'junit_family=xunit1' to your pytest.ini file to keep the current format in future versions of pytest and silence thiswarning.
                        #  _issue_warning_captured(deprecated.JUNIT_XML_DEFAULT_FAMILY, config.hook, 2)
                        #########
                        cmd = ["--junitxml={}".format(xmlfile)] + ["-s"] + ["--disable-pytest-warnings"] + cmd
                        #print("Running Command: pytest {}".format(cmd))
                        #print("Work Path: {}".format(workpath))
                        if len(os.listdir(workpath)) < 1: # If only the XML dir was created
                            print("No Tests to Run")
                            sys.exit()
                        else:
                            myworkdir = os.getcwd()
                            os.chdir(workdir + "{}/".format(dirname))
                            print("Test Directory: {}".format(os.getcwd()))
                            print("Running Command: pytest {}".format(cmd))
                            pytest.main(cmd)
                            os.chdir(myworkdir)
                            

            os.chdir(cwd)

    if not IS_CASA6: # If in CASA5
        if HAVE_NOSE:
            # Global variable used by regression framework to determine pass/failure status
            global regstate
            regstate = False
                
            listtests = testnames
            if listtests == '--Help' or listtests == []:
                usage()
                sys.exit()
                
            if listtests == '--list':
                list_tests()
                sys.exit()
                        
            if listtests == 'all':
                whichtests = 0
                # Get the full list of tests from file
                listtests = readfile(LISTofTESTS)
                if listtests == []:
                    raise Exception('List of tests "{}" is empty or does not exist'.format(LISTofTESTS))

            elif (type(testnames) != type([])):
                if (os.path.isfile(testnames)):
                    # How to prevent it from opening a real test???
                    whichtests = 1
                    listtests = readfile(testnames)
                    if listtests == []:
                        raise Exception('List of tests is empty')
                else:
                    raise Exception('List of tests does not exist')
                    
            else:
                # run specific tests
                whichtests = 1


            # Directories
            PWD = os.getcwd()
            WDIR = PWD+'/nosedir/'
            
            # Create a working directory
            workdir = WDIR
            print("Creating work directory: {}".format(workdir))
            if os.access(workdir, os.F_OK) is False:
                os.makedirs(workdir)
            else:
                shutil.rmtree(workdir)
                os.makedirs(workdir)
            
            # Move to working dir
            os.chdir(workdir)
            
            # Create a directory for nose's xml files
            xmldir = WDIR+'xml/'
            if os.access(xmldir, os.F_OK) is False:
                os.makedirs(xmldir)
            else:
                shutil.rmtree(xmldir)
                os.makedirs(xmldir)
            
            print("Starting unit tests for {}: ".format(listtests))

            # ASSEMBLE and RUN the TESTS
            if not whichtests:
                '''Run all tests'''
                list = []
                testlist_to_execute= []
                # awells CAS-10844 Fix
                suiteList = []
                for f in listtests:
                    suite = unittest.TestSuite()
                    try:
                        tests = UnitTest(f).getUnitTest()
                        if RUN_SUBTEST:
                            testlist_to_execute = testlist_to_execute + getsubtests(f,tests)
                        for test in tests:
                            suite.addTest(test)
                        suiteList.append(suite)
                    except:
                        traceback.print_exc()
                list = suiteList

            elif (whichtests == 1):
                '''Run specific tests'''
                list = []
                testlist_to_execute= []
                for f in listtests:
                    if not haslist(f):
                        testcases = UnitTest(f).getUnitTest()
                        list = list+testcases

                        if RUN_SUBTEST:
                            testlist_to_execute = testlist_to_execute + getsubtests(f,list)
                    else:
                        ff = getname(f)
                        tests = gettests(f)
                        # allow splitting of large test groups into smaller chunks
                        # large long running test groups make parallel test scheduling
                        # a bit more complicated so splitting them to smaller groups
                        # helps
                        # syntax: [testsplit:chunk_index-number_of_chunks]
                        if len(tests) == 1 and tests[0].startswith('testsplit:'):
                            import math
                            testcases = UnitTest(ff).getUnitTest()
                            chk, nchk = map(int, tests[0].split(':')[1].split('-'))
                            if chk > nchk or chk < 1:
                                raise ValueError('testsplit chunk must be 1 <= nchunks')
                            nchk = min(len(testcases), nchk)
                            chksz = int(math.ceil(len(testcases) / float(nchk)))
                            offset = (chk - 1) * chksz
                            print("running tests {} to {}".format(offset, min(offset + chksz, len(testcases))))
                            testcases = testcases[offset:offset + chksz]
                        else:
                            testcases = UnitTest(ff).getUnitTest(tests)
                        list = list+testcases

            if RUN_SUBTEST:


                if len(testlist_to_execute) == 0:
                    raise ValueError("Cannot Find Tests with Attribute:'{}'".format(ATTR_VAL))
                if not whichtests:
                    for i in range(0,len(list)):
                        tmp = []
                        for item in list[i]:
                            if [item._testMethodName,item.__module__] in testlist_to_execute:
                                tmp.append(item)
                        list[i] =  unittest.TestSuite(tmp)
                else:
                    tmp = []
                    for item in list:
                        if [item._testMethodName,item.__module__] in testlist_to_execute:
                            tmp.append(item)
                    list = tmp

            if (len(list) == 0):
                os.chdir(PWD)
                raise Exception('ERROR: There are no valid tests to run')


            # Run all tests and create a XML report
            xmlfile = xmldir+'nose.xml'

            try:
                if (HAVE_MEMTEST and MEM):
                    regstate = nose.run(argv=[sys.argv[0],"-d","-s","--with-memtest","--verbosity=2",
                                    "--memtest-file="+xmlfile], suite=list, addplugins=[memTest.MemTest()])
                else:
                    regstate = nose.run(argv=[sys.argv[0],"-d","-s","--with-xunit","--verbosity=2",
                                    "--xunit-file="+xmlfile], suite=list)

                os.chdir(PWD)
            except:
                print("Failed to run one or more tests")
                traceback.print_exc()
            else:
                os.chdir(PWD)


if __name__ == "__main__":

    # mem mode variables
    print("HAVE_MEMTEST: {}".format(HAVE_MEMTEST))
    print("HAVE_COVTEST: {}".format(HAVE_COVTEST))
    print("HAVE_ROBOT: {}".format(HAVE_ROBOT))
    print("HAVE_PYTEST: {}".format(HAVE_PYTEST))
    print("HAVE_NOSE: {}".format(HAVE_NOSE))
    print("HAVE_CASA6: {}".format(HAVE_CASA6))

    verbose = False

    #original_datapath = casa['dirs']['data']
    # List of tests to run
    testnames = []
        
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--all", action='store_true',help='run all tests defined in trunk/gcwrap/python/scripts/tests/uTest_list.json.')
    parser.add_argument("-l", "--list",action='store_true',help='print the list of tests & tags defined in ')
    parser.add_argument("-v", "--verbose",action='store_true',help="Verbose Test Execution")
    parser.add_argument("-x", "--dry-run",action='store_true',help="dry run Test Execution")
    parser.add_argument("-s", "--classes",nargs='+',metavar='test',help='print the classes from a test script')
    parser.add_argument("-f", "--file",nargs='?', type=argparse.FileType('r'),help='run the tests defined in an ASCII file <list>; one test per line')
    parser.add_argument("-m", "--mem",action='store_true',help='show the memory used by the tests and the number of files left open.')
    
    # Component Arguments
    parser.add_argument("-p","--components", help='Comma separated list of components', required=False)
    parser.add_argument("-e","--mapfile", nargs='?', type=argparse.FileType('r'), help='Component to test map file', required=False)

    parser.add_argument("-b","--branch", help='JIRA Branch to checkout', required=False)

    if not IS_CASA6:
        if "-c" in sys.argv:
            i = sys.argv.index("-c")
        args, unknownArgs = parser.parse_known_args(sys.argv[i+2:])
    else:
        args, unknownArgs = parser.parse_known_args()
    
    print(args)
    # print(unknownArgs)

    if args.components is not None:
        components = args.components
        components = components.split(",")
        print("Testing Components" + str(components))
        
        #
        if args.mapfile is not None:
            component_to_test_map = json.load(args.mapfile)
        else:
            import casatestutils as _;
            with open("{}/{}".format(_.__path__[0], "component_to_test_map.json")) as ctt:
                component_to_test_map = json.load(ctt)
        
        for c in components:
            _isComponent = False
            component = c.strip()
            for myDict in component_to_test_map["testlist"]:
                if component in myDict["testGroup"]:
                    _isComponent = True
                    testnames.append(myDict["testScript"])
            if not _isComponent:
                print("No Tests for Component: {}".format(component))

    if args.verbose:
        verbose = True
    if args.list:
        list_tests()
        sys.exit()
    ## Dry Run
    DRY_RUN = False
    if args.dry_run:
        DRY_RUN = True

    ## RUN ALL
    RUN_ALL = False
    if args.all:
        RUN_ALL = True
        testnames = ["all"]

    if args.classes is not None:
        print(args.classes)
        getclasses(args.classes)
        os._exit(0)
        

    if args.file is not None:
        logger.info('Reading Test List from %s: ', args.file)
        for line in args.file:
            try:
                logger.debug("Adding Test %s from file %s",re.sub(r'[\n\r]+', '',line),args.file)
                testnames.append(re.sub(r'[\n\r]+', '',line))
            except:
                raise Exception(" The list should contain one test per line.")
                
    if args.mem: # run specific tests in mem mode
        logger.info('Setting Mem Mode')
        MEM = 1


    if args.branch is not None:
        JIRA_BRANCH = args.branch

    for arg in unknownArgs:
        if arg.startswith(("-", "--")):
            raise ValueError('unrecognized argument: %s'%(arg))
            sys.exit()
        #
        else:
            if arg.startswith("test") and not RUN_ALL:
                testnames.append(arg)
            # local Path file
            #elif arg.startswith("./") and arg.endswith(".py") and not RUN_ALL:
            elif arg.startswith("./") and ".py" in arg and not RUN_ALL:

                testnames.append(arg[2:])

            elif (arg.startswith("../") or arg.startswith("/"))  and ".py" in arg and not RUN_ALL:
                try:
                    real_path = os.path.realpath(arg)
                    #print("real_path: {}".format(real_path))
                    testnames.append(real_path)
                except:
                    traceback.print_exc()
                
    #If no tests are given, no subet tag or --all option
    if testnames == [] or len(testnames) == 0:
        #TODO
        print("List of tests is empty")
        parser.print_help(sys.stderr)
        sys.exit(1)
        #sys.exit()
        #raise Exception("List of tests is empty")

#    print(testnames)
#    sys.exit()
    try:
        run(testnames)

    except:
        traceback.print_exc()


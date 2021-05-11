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
import datetime
import platform
from testrunner.shell_runner import ShellRunner
from testrunner import xvfb_helper
from testrunner.xunit import Xunit
default_timeout = 1800
sys.path.insert(0,'')

########################################################################################################################
###########################################            Functions            ############################################
########################################################################################################################

class casa_test:
    def __init__(self,
                 name,
                 path,
                 test_group=None,
                 test_type=None,
                 maintainer=None,
                 email=None,
                 options=None,
                 comment=None,
                 timeout = default_timeout):
        self.name = name
        self.path = path
        self.test_type = test_type
        self.maintainer = maintainer
        self.email = email
        self.options = options
        self.test_group = test_group
        self.comment = comment
        self.timeout = timeout

    def __eq__(self, other):
        return other is not None and \
               self.name == other.name and \
               self.path == other.path and \
               self.options == other.options

    def __hash__(self):
        return hash(('name', self.name,'path', self.path, 'options', self.options))

def fetch_tests(work_dir, branch):

    repo_path = "https://open-bitbucket.nrao.edu/scm/casa/"
    source_dir=work_dir + "/casasources"
    if not os.path.exists(source_dir):
        os.makedirs(source_dir)
    repositories = ["casampi", "casaplotms", "almatasks", "casa6","casaviewer"]
    # All of the repositositories have their tests in different directories
    # so we need a mapping
    def get_repo_test_paths(x):
        return {
            "casa6": ["/casa6/casatests/regression/","/casa6/casatests/stakeholder/","/casa6/casatasks/tests/","/casa6/casatools/tests/"],
            "casampi": ["/casampi/src/casampi/tests"],
            "casaplotms": ["/casaplotms/tests/plotms"],
            "almatasks": ["/almatasks/tests/tasks"],
            "casaviewer": ["/casaviewer/tests/tasks"]
        }[x]

    # Clone the repository and checkout branch
    for repo in repositories:
        cmd = ("git clone " + repo_path + repo).split()
        print(cmd)
        r = ShellRunner()
        r.runshell(cmd, default_timeout, source_dir)
        cmd = ("git checkout " + branch).split()
        print(cmd)
        r = ShellRunner()
        r.runshell(cmd, default_timeout, source_dir + "/" + repo)
        for x in get_repo_test_paths(repo):
            test_paths.append(source_dir + "/" + x)
    return test_paths

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

def git_fetch_casa_tests(path):

    cwd = os.getcwd()
    if os.path.exists(path):
        try:
            os.rmdir(path)
        except:
            shutil.rmtree(path)
    os.makedirs(path)

    os.chdir(path)
    subprocess.call(["git","init", "--quiet"])
    FNULL = open(os.devnull, 'w')
    subprocess.call(["git","remote","add","-f","origin", "https://open-bitbucket.nrao.edu/scm/casa/casa6.git"], stdout=FNULL, stderr=subprocess.STDOUT)
    subprocess.call(["git","config","core.sparseCheckout","true"])

    print("casatasks/tests/tasks", file=open(".git/info/sparse-checkout", "a"))
    
    print("casatools/tests/tools/agentflagger", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/calanalysis", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/componentlist", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/coordsys", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/image", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/imagepol", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/measures", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/ms", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/msmetadata", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/regionmanager", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/sdm", file=open(".git/info/sparse-checkout", "a"))
    print("casatools/tests/tools/vpmmanager", file=open(".git/info/sparse-checkout", "a"))

    print("casatests/benchmarks", file=open(".git/info/sparse-checkout", "a"))
    print("casatests/e2e", file=open(".git/info/sparse-checkout", "a"))
    print("casatests/performance", file=open(".git/info/sparse-checkout", "a"))
    print("casatests/pipeline", file=open(".git/info/sparse-checkout", "a"))
    print("casatests/stakeholders", file=open(".git/info/sparse-checkout", "a"))

    subprocess.call(["git","pull","--quiet","origin","master"])
    os.chdir(cwd)

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
        #logger.debug("Start def getclasses()")
        '''Get the classes of a test script It will copy the test script to /tmp and remove it afterwards'''
        here = os.getcwd()
        tmpdir = '/tmp'
        for filename in testnames:
            if not filename.startswith("test_"):
                #print "Cannot Get Classes for Regression Test"
                #logger.error("Cannot Get Classes for Regression Test: %s",filename)
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
            #logger.error('Failed to open file', exc_info=True)
            return

def haslist(name):
    '''Check if specific list of tests have been requested'''
    n0 = name.rfind('[')
    n1 = name.rfind(']')
    if n0 == -1:
        return False
    return True

def unpack_dmg(pkg, work_dir, outputdir):
    mountpoint = work_dir + "/mnt"
    if not os.path.exists(mountpoint):
        os.makedirs(mountpoint)

    print ("Unpacking dmg: " + pkg + " to " +  outputdir)
    cmd = ("hdiutil attach " + pkg + " -mountpoint " + mountpoint).split()
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    installpath = outputdir + "/CASA.app"
    cmd = ("ditto " + mountpoint + "/CASA.app " + outputdir + "/CASA.app").split()
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    cmd = ("hdiutil detach " + mountpoint).split()
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    return installpath
    
def unpack_tarball(pkg, outputdir):
    print ("Unpacking tarball: " + pkg + " to " +  outputdir)
    cmd = ("tar -xf " + pkg + " -C " + outputdir).split()
    print(cmd)
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    installpath = outputdir + "/" + os.path.basename(pkg).replace(".tar.xz","")
    return installpath

def get_casatestutils_exec_path(pkg_dir):
    for currentpath, folders, files in os.walk(pkg_dir):
        for file in files:
            #print(">>>" + os.path.join(currentpath, file))
            if currentpath.endswith('casatestutils') and file == 'runtest.py':
                return(os.path.join(currentpath, file))

def unpack_pkg(pkg, work_dir, outputdir):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    if platform.system() == "Linux":
        installpath = unpack_tarball(pkg, outputdir)
        print ("Package root: " + installpath)
        exec_path = installpath + "/bin"
    elif platform.system() == "Darwin":
        installpath = unpack_dmg(pkg,work_dir, outputdir)
        print("Package root: " + installpath)
        exec_path = installpath + "/Contents/MacOS"
    else:
        raise Exception("Unknown operating system")
    if exec_path is None:
        raise Exception ("Couldn't find casa executable path")
    casatestutils_exec_path = get_casatestutils_exec_path(installpath)
    if casatestutils_exec_path == None:
        raise Exception("Couldn't find casatestutils")
    return exec_path, casatestutils_exec_path

def run_bamboo(pkg, work_dir, branch = None, test_group = None, test_list= None, test_paths = [], test_config_path=None, ncores=2, verbosity=False):

    if test_list is not None:
        test_list = [x.strip() for x in test_list.split(',')]
    if args.test_group is not None:
        test_group = [x.strip() for x in test_group.split(',')]
    # Unpack the distribution
    print ("Test list: " + str (test_list))
    print ("Test group: " + str (test_group))
    if pkg is None:
        raise Exception("Missing pkg")
    if work_dir is None:
        raise Exception("Missing work_dir")
    exec_path, casatestutils_exec_path = unpack_pkg(pkg, work_dir, work_dir + "/pkg")
    #exec_path = "/tmp/work/pkg/CASA.app/Contents/MacOS"
    print("Executable path: " + exec_path)
    print("casatestutils path: " + casatestutils_exec_path)

    # Start Xvfb on Linux
    xvfb = xvfb_helper.XvfbHelper()
    if sys.platform != "darwin":
        xvfb.start_virtual_frame_buffer()

    test_config = None
    if test_config_path == None:
       test_config_path = work_dir + "/casasources/casa6/casatestutils/casatestutils/component_to_test_map.json"
    # Read the JSON configuration
    print ("Reading config from" + test_config_path)
    with open(test_config_path ) as f:
      test_config = json.load(f)

    # Get the actual tests as list
    test_config_elems = test_config['testlist']

    if args.branch == None:
        branch = "master"

    # Clone a default set of repositories to if test paths are not provided from command line
    if len(test_paths) == 0 :
        test_paths = fetch_tests(str(work_dir), branch)
    # Match the test names provided in the JSON file to the actual test locations.
    tests_to_run = []
    for x in test_config_elems:
        for dir in test_paths:
            for currentpath, folders, files in os.walk(dir):
                for file in files:
                    if file == (x['testScript']+".py"):
                        test_location = os.path.join(currentpath, file)
                        if verbosity:
                            print("Found: " + test_location)
                            print("Script:", x['testScript'], test_location, x['testGroup'])
                        if x['timeout'] < 1:
                            timeout = default_timeout
                        else:
                            timeout = x['timeout']

                        opts = x['testOptions']
                        opts.sort()

                        tests_to_run.append(casa_test(x['testScript'],
                                                      test_location,
                                                      x['testGroup'],
                                                      x['testType'],
                                                      x['maintainer'],
                                                      x['maintainerEmail'],
                                                      tuple(opts),
                                                      x['comment'],
                                                      timeout))

    # Filter tests by test list
    if test_list is not None and len(test_list)>0:
        print ("Test list provided. Filtering tests.")
        tmp_tests_to_run = []
        for test in test_list:
            found = False
            for t1 in tests_to_run:
                if test == t1.name:
                    tmp_tests_to_run.append(t1)
                    found = True
            if not found:
                # If there is a test in the list but no configuration, add it without any options.
                # This can be useful for testing new scripts.
                print ("Test " + test + " configuration not found. Searching for the test...")
                print ("dir: " + dir )
                # Potential test location
                for dir in test_paths:
                    # Search for files (needed for casatools and casatasks)
                    for currentpath, folders, files in os.walk(dir):
                        for file in files:
                            if file == test + ".py":
                                test_location = os.path.join(currentpath, file)
                                print("Found: " + test_location)
                                print("No JSON configuration found. Test will be added to execution list without options.")
                                tmp_tests_to_run.append(casa_test(test, test_location, "","","","","",""))
                                found = True
            if not found:
                raise Exception("Couldn't locate test: " + test)
        tests_to_run = tmp_tests_to_run

    # Filter by Jira components
    if test_group is not None and len(test_group)>0:
        print("Jira component list provided. Filtering tests.")
        tmp_tests_to_run = []
        for jira_component in test_group:
            found = False
            for t1 in tests_to_run:
                if jira_component.lower() in [x.lower() for x in t1.test_group]:
                    tmp_tests_to_run.append(t1)
                    found = True
            # Throw an exception is user provides a component that doesn't exist
            if not found:
                print ("WARNING: No tests found for jira_component " + jira_component + ". Check the contents of " + test_config_path)
        # Remove duplicates
        tests_to_run = set(tmp_tests_to_run)

    print("Running these tests:")
    for t in tests_to_run:
        if verbose:
            print(t.name + " : " +
                  t.path + " : " +
                  t.test_type + " : " +
                  t.maintainer + " : " +
                  t.email + " : " +
                  str(t.options) + " : " +
                  str(t.test_group) + " : " +
                  str(t.timeout) + " : " +
                  t.comment)
        else:
            print(t.name)

    # Run tests
    for test in tests_to_run:
        r = ShellRunner()
        xunit = Xunit()
        casaopts = " --nogui --nologger --log2term"

        # Skip MPI on Darwin for now
        if "mpi" in test.options and sys.platform != "darwin":
            print("Running tests in MPI mode")
            casa_exe = exec_path + "/mpicasa"
            casaopts = "-n " + str(ncores) + " " + exec_path + "/casa" + " --nogui --nologger --log2term --agg "
        else:
            casa_exe = exec_path + "/casa"
      
        assert (test != None)
        #cmd = (casa_exe + " " + casaopts + " -c " + casatestutils_exec_path + " " + test.path).split()
        cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
        cwd = work_dir + "/" + test.name
        if not os.path.exists(cwd):
            os.makedirs(cwd)
        print("Running cmd " + str(cmd))
        print("in "  + cwd)
        starttime = datetime.datetime.now()

        output = r.runshell(cmd, test.timeout,cwd)

        endtime = datetime.datetime.now()
        runtime = endtime - starttime
        xunit.append_result(test.name, str(runtime), len(output), output)
    xunit.generateXml("suite")

    # Close Xvfb on exit
    import atexit
    @atexit.register
    def goodbye():
        print("Stopping Xvfb.")
        if sys.platform != "darwin":
            xvfb.signal_stop_virtual_frame_buffer()
        print("Xvfb stopped.")
########################################################################################################################
######################################            Imports / Constants            #######################################
########################################################################################################################
# Python Version
PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

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
    from __main__ import default
    from tasks import *
    from taskinit import *
except ImportError:
    CASA6 = True
    IS_CASA6 = True
    HAVE_CASA6 = True

# Use Nose attribute Functionality
RUN_SUBTEST = False

# Dry run of Tests
DRY_RUN = False

# Define which tests to run
whichtests = 0

########################################################################################################################
#############################################            Main            ###############################################
########################################################################################################################

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

########################################################################################################################
#############################################            Main            ###############################################
########################################################################################################################
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
        
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("-a", "--all", action='store_true',help='run all tests defined in trunk/gcwrap/python/scripts/tests/uTest_list.json.')
    parser.add_argument("-i", "--list",action='store_true',help='print the list of tests & tags defined in ')
    parser.add_argument("-v", "--verbose",action='store_true',help="Verbose Test Execution")
    parser.add_argument("-x", "--dry-run",action='store_true',help="dry run Test Execution")
    parser.add_argument("-s", "--classes",nargs='+',metavar='test',help='print the classes from a test script')
    parser.add_argument("-f", "--file",nargs='?', type=argparse.FileType('r'),help='run the tests defined in an ASCII file <list>; one test per line')
    parser.add_argument("-m", "--mem",action='store_true',help='show the memory used by the tests and the number of files left open.')
    
    # Component Arguments

    parser.add_argument("-e","--mapfile", nargs='?', type=argparse.FileType('r'), help='Component to test map file', required=False)
    parser.add_argument("-b","--branch", help='JIRA Branch for test repository checkouts', required=False)

    # casa-build-utils Arguments
    parser.add_argument('-p','--pkg', help='Tarball or dmg', required=False)
    parser.add_argument('-w','--work_dir', help='Working directory.', required=False)
    parser.add_argument('-n','--ncores', help='Number of cores for MPI tests', default=2)
    parser.add_argument('-t','--test_paths', help='A comma separated list of paths containing tests.', required=False)
    parser.add_argument('-l','--test_list', help='Filter tests by a comma separated list of tests', required=False)
    parser.add_argument('-c','--test_config',  help='Test configuration file', required=False)
    parser.add_argument('-j','--test_group',  help='Filter tests by a comma separated list of components', required=False)
    # parser.add_argument('-b','--branch',  help='Branch for test repository checkouts', required=False) ## Duplicate

    parser.add_argument('--bamboo', help='Set Bamboo Flag to True',default=False,action='store_true', required=False)

    if not IS_CASA6:
        if "-c" in sys.argv:
            i = sys.argv.index("-c")
        args, unknownArgs = parser.parse_known_args(sys.argv[i+2:])
    else:
        args, unknownArgs = parser.parse_known_args()
    
    print(args)
    #args=parser.parse_args()
    if args.pkg:
        print("Package: " + args.pkg)
    print("Test configuration file: " + str(args.test_config))
    print("Number of cores: " + str(args.ncores))
    print("Workdir: " + str(args.work_dir))

    print("Operating system: " +  platform.system())
    # print(unknownArgs)


    if args.test_group is not None:
        components = args.test_group
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
        #logger.info('Reading Test List from %s: ', args.file)
        for line in args.file:
            try:
                #logger.debug("Adding Test %s from file %s",re.sub(r'[\n\r]+', '',line),args.file)
                testnames.append(re.sub(r'[\n\r]+', '',line))
            except:
                raise Exception(" The list should contain one test per line.")
                
    if args.mem: # run specific tests in mem mode
        #logger.info('Setting Mem Mode')
        MEM = 1


    if args.branch is not None:
        JIRA_BRANCH = args.branch

    test_paths = []
    if args.test_paths is not None:
        test_paths = [x.strip() for x in args.test_paths.split(',')]

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

                #testnames.append(arg[2:])
                try:
                    real_path = os.path.realpath(arg[2:])
                    #print("real_path: {}".format(real_path))
                    testnames.append(real_path)
                except:
                    traceback.print_exc()

            elif (arg.startswith("../") or arg.startswith("/"))  and ".py" in arg and not RUN_ALL:
                try:
                    real_path = os.path.realpath(arg)
                    #print("real_path: {}".format(real_path))
                    testnames.append(real_path)
                except:
                    traceback.print_exc()

    try:

        if args.bamboo:
            run_bamboo(args.pkg, args.work_dir, args.branch, args.test_group, args.test_list, test_paths, args.test_config, args.ncores, args.verbose)

        else:
            #If no tests are given, no subet tag or --all option
            #print(testnames)
            if testnames == [] or len(testnames) == 0:
                #TODO
                print("List of tests is empty")
                parser.print_help(sys.stderr)
                sys.exit(1)
                #sys.exit()
                #raise Exception("List of tests is empty")
            run(testnames)
    except:
        traceback.print_exc()



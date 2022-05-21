########################################################################################################################
############################################            Imports            #############################################
########################################################################################################################
import argparse
import os
import shutil
import sys
import traceback
import subprocess
import unittest
import json
import datetime
import platform

########################################################################################################################
######################################            Imports / Constants            #######################################
########################################################################################################################

default_timeout = 1800
sys.path.insert(0,'')

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
verbose = False

# JIRA BRANCH TO CHECKOUT
JIRA_BRANCH = None

try:
    from __main__ import default
    from tasks import *
    from taskinit import *
except ImportError:
    CASA6 = True
    IS_CASA6 = True

# Dry run of Tests
DRY_RUN = False

########################################################################################################################
###########################################            Functions            ############################################
########################################################################################################################
# At the moment, this needs to be a sep function due to repr and escape characters, try/ except for osx
def write_conftest_linux(filepath):
    string = """
import pytest
import inspect
import os

@pytest.mark.trylast
def pytest_configure(config):
    terminal_reporter = config.pluginmanager.getplugin('terminalreporter')
    config.pluginmanager.register(TestDescriptionPlugin(terminal_reporter), 'testdescription')

class TestDescriptionPlugin:

    def __init__(self, terminal_reporter):
        self.terminal_reporter = terminal_reporter
        self.desc = None
        self.funcn = None

    def pytest_runtest_protocol(self, item):
        #from pprint import pprint
        #d = item.__dict__
        #pprint(d, indent=2)
        self.desc = inspect.getdoc(item.obj)
        #print(item._nodeid)
        self.funcn = item._nodeid

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_logstart(self, nodeid, location):
        #print("Verbosity Level: {}".format(self.terminal_reporter.verbosity))
        if self.terminal_reporter.verbosity == 0:
            yield
            self.terminal_reporter.write(f'\\n{self.funcn} \\n')
        else:
            self.terminal_reporter.write('\\n')
            yield
            if self.desc:
                    self.terminal_reporter.write(f'\\n{self.desc} \\n')
            else:
                    self.terminal_reporter.write(f'\\n')

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(item, call):
        outcome = yield
        report = outcome.get_result()
        #print(dir(report))
        report.start = call.start
        report.stop = call.stop
        if report.when=='teardown':
            filepath = os.path.join(os.getcwd(),'short_summary.log')

            file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
            file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid,))
            file_obj.close()

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(item, call):
        outcome = yield
        report = outcome.get_result()
        if report.when=='call':
            filepath = os.path.join(os.getcwd(),'short_summary.log')
            # write short summary to file
            file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
            file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid))
            file_obj.close()

            # Write not pass to Textfile
            if report.outcome != 'passed':
                file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
                file_obj.write("\\tDuration: {}s\\n".format(round(report.duration,5)))
                if report.outcome == 'failed':
                    file_obj.write("\\tMessage : {}\\n".format(report.longrepr.reprcrash.message))
                file_obj.close()
                filepath = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')),'summary_of_failed.log')
                file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
                file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid))
                file_obj.write("\\tDuration: {}s\\n".format(round(report.duration,5)))
                if report.outcome == 'failed':
                    file_obj.write("\\tMessage : {}\\n".format(report.longrepr.reprcrash.message))
                file_obj.close()
    """
    file_obj = open(filepath,'w')
    file_obj.write(string)
    file_obj.close()

# At the moment, this needs to be a sep function due to repr and escape characters, try/ except for osx
def write_conftest_osx(filepath):
    string = """
import pytest
import inspect
import os

@pytest.mark.trylast
def pytest_configure(config):
    terminal_reporter = config.pluginmanager.getplugin('terminalreporter')
    try:
        config.pluginmanager.unregister(TestDescriptionPlugin(terminal_reporter), 'testdescription')
    except:
        pass
    config.pluginmanager.register(TestDescriptionPlugin(terminal_reporter), 'testdescription')

class TestDescriptionPlugin:

    def __init__(self, terminal_reporter):
        self.terminal_reporter = terminal_reporter
        self.desc = None
        self.funcn = None

    def pytest_runtest_protocol(self, item):
        #from pprint import pprint
        #d = item.__dict__
        #pprint(d, indent=2)
        self.desc = inspect.getdoc(item.obj)
        #print(item._nodeid)
        self.funcn = item._nodeid

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_logstart(self, nodeid, location):
        #print("Verbosity Level: {}".format(self.terminal_reporter.verbosity))
        if self.terminal_reporter.verbosity == 0:
            yield
            self.terminal_reporter.write(f'\\n{self.funcn} \\n')
        else:
            self.terminal_reporter.write('\\n')
            yield
            if self.desc:
                    self.terminal_reporter.write(f'\\n{self.desc} \\n')
            else:
                    self.terminal_reporter.write(f'\\n')

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(item, call):
        outcome = yield
        report = outcome.get_result()
        if report.when=='call':
            filepath = os.path.join(os.getcwd(),'short_summary.log')
            # write short summary to file
            file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
            file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid))
            file_obj.close()

            # Write not pass to Textfile
            if report.outcome != 'passed':
                file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
                file_obj.write("\\tDuration: {}s\\n".format(round(report.duration,5)))
                if report.outcome == 'failed':
                    file_obj.write("\\tMessage : {}\\n".format(report.longrepr.reprcrash.message))
                file_obj.close()
                filepath = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')),'summary_of_failed.log')
                file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
                file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid))
                file_obj.write("\\tDuration: {}s\\n".format(round(report.duration,5)))
                if report.outcome == 'failed':
                    file_obj.write("\\tMessage : {}\\n".format(report.longrepr.reprcrash.message))
                file_obj.close()
    """
    file_obj = open(filepath,'w')
    file_obj.write(string)
    file_obj.close()

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

def read_conf(conf):
    with open(conf) as f:
        lines = [line.rstrip() for line in f]
    return dict(x.split('==') for x in lines)

def fetch_tests(work_dir, branch):

    repo_path = "https://open-bitbucket.nrao.edu/scm/casa/"
    source_dir=work_dir + "/casasources"
    if not os.path.exists(source_dir):
        os.makedirs(source_dir)
    repositories = ["casa6", "casampi", "casaplotms", "almatasks","casaviewer"]
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
        try:
            r = ShellRunner()
            r.runshell(cmd, default_timeout, source_dir)
        except:
            cwd = os.getcwd()
            os.chdir(source_dir)
            subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
            os.chdir(cwd)
        if branch != 'master' and repo != 'casa6':
            cmd = 'git ls-remote --heads {}{} {} | wc -l'.format(repo_path, repo, branch )
            #print(cmd)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
            out = proc.stdout.read()
            if int(out)== 0: 
                if repo in ['casaplotserver', 'casaplotms','casaviewer','casampi','almatasks','casatelemetry']:
                    print("build.conf location: " + source_dir + "/casa6/build.conf" )
                    branchtag = "tags/{}".format(read_conf(source_dir+"/casa6/build.conf")[repo])
                    print("branchtag" + branchtag)
                cmd = ("git checkout " + branchtag).split()
            else:
                cmd = ("git checkout " + branch).split()
        elif branch == 'master' and repo != 'casa6' and os.path.isfile(source_dir+"/casa6/build.conf"):
            if repo in ['casaplotserver', 'casaplotms','casaviewer','casampi','almatasks','casatelemetry']:
                branchtag = "tags/{}".format(read_conf(source_dir+"/casa6/build.conf")[repo])
            cmd = ("git checkout " + branchtag).split()
        else:
            cmd = ("git checkout " + branch).split()
        print(cmd)
        try:
            r = ShellRunner()
            r.runshell(cmd, default_timeout, source_dir + "/" + repo)
        except:
            cwd = os.getcwd()
            os.chdir(source_dir + "/" + repo)
            subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
            os.chdir(cwd)
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
        if os.path.isdir(os.getcwd() +"/testlist/"):
            shutil.rmtree(os.getcwd() +"/testlist/")
        os.makedirs(os.getcwd() +"/testlist/")
        testpaths = fetch_tests(os.getcwd() +"/testlist/", 'master')
        for path in testpaths:
            gather_all_tests(path, os.getcwd() +"/testlist/")
        tests = sorted(os.listdir(os.getcwd() +"/testlist/"))
        for test in tests:
            if test.startswith("test_"):
                print(test)
    else:
        for t in readfile(LISTofTESTS):
            print(t)

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
        pytest.main(["--collect-only", "-q", "-v", test])

def clean_working_directory(workdir):

    print("Cleaning: {}".format(workdir))
    if os.path.exists(workdir):
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
   
    installpath = None
    
    print("outputdir contents:" + outputdir)
    for root, dirs, files in os.walk(outputdir):
        for d in dirs:
            print(" " + d)
            if d.startswith("casa-"):
                installpath = d
                print("installpath: " + installpath)
        break

    if installpath is None:
        raise  RuntimeError("Couldn't find a directory that looks like a Casa distribution. Expected directory name to start with 'casa-'")  
    return outputdir + "/" + installpath

def get_casatestutils_exec_path(pkg_dir):
    # Since runtest is no longer part of casatestutils, this may be removed.
    for currentpath, folders, files in os.walk(pkg_dir):
        for file in files:
            #print(">>>" + os.path.join(currentpath, file))
            if currentpath.endswith('casatestutils') and file == 'runtest.py':
                return(os.path.join(currentpath, file))
    return "/dev/null/"

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

########################################################################################################################
##############################################            Run            ###############################################
########################################################################################################################

def run(testnames, branch=None, DRY_RUN=False):

    if IS_CASA6:

        if HAVE_PYTEST:
            cwd = os.getcwd() + "/"
            workpath = os.getcwd() +"/nosedir/"
            workdir = os.getcwd() +"/nosedir/"

            clean_working_directory(workpath)
            # Copy Tests to Working Directory
            os.makedirs(workdir)

            # Remove Duplicates
            # Since working directory is based on script name
            # We need to remove multiple calls to the same script
            setlist = []
            for duplicate in list(set([ x.split("/")[-1] for x in testnames])):
                inlist = True
                for test in testnames:
                    if duplicate in test:
                        if inlist:
                            setlist.append(test)
                            inlist = False
            testnames = setlist
            print("Tests: {}".format(sorted(testnames)))
            gittest = True
            if branch ==None:
                branch = 'master'
            # Only Checkout When Needed
            if any([False if ".py" in x else True for x in testnames ]):
                testpaths = fetch_tests(workdir, branch)
                os.makedirs(workdir + "tests/")
                for path in testpaths:
                    gather_all_tests(path, workdir + "tests/")
                print(workdir + "tests/")

            for testname in testnames:
                cmd = []

                # Copy Test To nosedir Directory if in cwd
                if testname.startswith("test"):
                    test = testname
                    # Check if specific tests are requested
                    if "[" and "]" in test:
                        testname = getname(test)
                        tests = gettests(test)

                        teststring = ""
                        if len(tests) == 1:
                            teststring = tests[0]
                        elif len(tests) > 1:
                            print(tests)
                            teststring = " or ".join(tests)

                        cmd = ["-k {}".format(teststring)] + cmd
                        test = testname

                    # Set up Test Working Directory
                    if not os.path.exists(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])):
                        print("Setting Working Directory: {}".format(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3])))
                        os.makedirs(workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]))
                        cmd = [ workdir + "{}/".format(test if not test.endswith(".py") else test[:-3]) ] + cmd

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
                    
                    cmd = ["--verbose"] + ["-ra"] + ["--tb=short"] + cmd
                    if DRY_RUN:
                        cmd = ["--collect-only"] + cmd
                    
                    if not os.path.isdir(workpath + '/xml/{}/'.format(test if not test.endswith(".py") else test[:-3])):
                        os.makedirs(workpath + '/xml/{}/'.format(test if not test.endswith(".py") else test[:-3]))
                    xmlfile = workpath + 'xml/{}/nose.xml'.format(test if not test.endswith(".py") else test[:-3])

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
                        conf_name = os.path.join(os.getcwd(),"conftest.py")
                        if platform.system() == 'Darwin':
                            write_conftest_osx(conf_name)
                        else:
                            write_conftest_linux(conf_name)
                        try: 
                            from casampi.MPIEnvironment import MPIEnvironment
                            if MPIEnvironment.is_mpi_enabled:
                                pytest.main(cmd)
                            else:
                                subprocess.run([sys.executable,"-m","pytest"] + cmd, env={**os.environ})
                        except:
                            subprocess.run([sys.executable,"-m","pytest"] + cmd, env={**os.environ})
                        os.remove(conf_name)
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
                        teststring = ""
                        if len(tests) == 1:
                            teststring = tests[0]
                        elif len(tests) > 1:
                            print(tests)
                            teststring = " or ".join(tests)
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

                    cmd = ["--verbose"] + ["-ra"] + ["--tb=short"] + cmd

                    if DRY_RUN:
                        cmd = ["--collect-only"] + cmd
                    
                    if not os.path.isdir(workpath + '/xml/{}/'.format(dirname)):
                        os.makedirs(workpath + '/xml/{}/'.format(dirname))
                    xmlfile = workpath + 'xml/{}/nose.xml'.format(dirname)
                    cmd = ["--junitxml={}".format(xmlfile)] + ["-s"] + ["--disable-pytest-warnings"] + cmd
                    if len(os.listdir(workpath)) < 1: # If only the XML dir was created
                        print("No Tests to Run")
                        sys.exit()
                    else:
                        myworkdir = os.getcwd()
                        os.chdir(workdir + "{}/".format(dirname))
                        print("Test Directory: {}".format(os.getcwd()))
                        print("Running Command: pytest {}".format(cmd))
                        conf_name = os.path.join(os.getcwd(),"conftest.py")
                        if platform.system() == 'Darwin':
                            write_conftest_osx(conf_name)
                        else:
                            write_conftest_linux(conf_name)
                        try: 
                            from casampi.MPIEnvironment import MPIEnvironment
                            if MPIEnvironment.is_mpi_enabled:
                                pytest.main(cmd)
                            else:
                                subprocess.run([sys.executable,"-m","pytest"] + cmd, env={**os.environ})
                        except:
                            subprocess.run([sys.executable,"-m","pytest"] + cmd, env={**os.environ})
                        os.remove(conf_name)
                        os.chdir(myworkdir)
            os.chdir(cwd)

        if not HAVE_PYTEST:
            raise ImportError('No Module Named Pytest. Pytest is Required for runtest.py')


########################################################################################################################
#######################################            Run Bamboo Option            ########################################
########################################################################################################################

def run_bamboo(pkg, work_dir, branch = None, test_group = None, test_list= None, test_paths = [], test_config_path=None, ncores=2, verbosity=False, pmode=None, tests_to_ignore=None):

    if test_list is not None:
        test_list = [x.strip() for x in test_list.split(',')]
    if args.test_group is not None:
        test_group = [x.strip() for x in test_group.split(',')]

    # Unpack the distribution
    print ("run_bamboo")
    print ("Test list: " + str (test_list))
    print ("Test group: " + str (test_group))

    if pkg is None:
        raise Exception("Missing pkg")
    if work_dir is None:
        raise Exception("Missing work_dir")
    exec_path, casatestutils_exec_path = unpack_pkg(pkg, work_dir, work_dir + "/pkg")

    print("Executable path: " + exec_path)
    print("casatestutils path: " + casatestutils_exec_path)

    # Start Xvfb on Linux
    xvfb = xvfb_helper.XvfbHelper()
    if sys.platform != "darwin":
        xvfb.start_virtual_frame_buffer()

    if args.branch == None:
        branch = "master"

    print ("run_bamboo fetch_tests branch" + branch)
    
    # Clone a default set of repositories to if test paths are not provided from command line
    if len(test_paths) == 0 :
        test_paths = fetch_tests(str(work_dir), branch)

    test_config = None
    if test_config_path == None:
       test_config_path = work_dir + "/casasources/casa6/casatestutils/casatestutils/component_to_test_map.json"
    # Read the JSON configuration
    print ("Reading config from: " + test_config_path)
    with open(test_config_path ) as f:
      test_config = json.load(f)

    # Get the actual tests as list
    test_config_elems = test_config['testlist']

    print("Test Paths: ", test_paths)

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
        #print(test_group)
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

    print("Subset tests:")
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
    print("")

    if tests_to_ignore is not None:
        print("\nTests to Ignore: ",tests_to_ignore )
        indices = []
        for i, t in enumerate(tests_to_run):
            if t.name in tests_to_ignore:
                indices.append(i)
        tests_to_run = [v for i,v in enumerate(tests_to_run) if i not in indices]

    for test in tests_to_run:
        r = ShellRunner()
        xunit = Xunit()

        #pmodes: pmodes = ['serial','parallel','both']
        # Skip MPI on Darwin for now
        if "mpi" in test.options and sys.platform != "darwin" and ( pmode == 'parallel' or pmode == 'both'):
            print("Running test: {} in MPI mode".format(test.name))
            casa_exe = exec_path + "/mpicasa"
            casaopts = "-n " + str(ncores) + " " + exec_path + "/casa" + " --nogui --nologger --log2term --agg " + rcdir + " "
            assert (test != None)
            cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
            cwd = work_dir + "/" + test.name
            if pmode == 'both': 
                cwd = work_dir + "/" + test.name + '_mpi'

            print("Running cmd " + str(cmd) + "in " + cwd)
            if not os.path.exists(cwd):
                os.makedirs(cwd)
            starttime = datetime.datetime.now()
            output = r.runshell(cmd, test.timeout,cwd)
            endtime = datetime.datetime.now()
            runtime = endtime - starttime
            xunit.append_result(test.name, str(runtime), len(output), output)
            print("")

            if pmode == 'both':
                if "casampi" in test.name:
                    continue
                print("Running test: {} in Serial mode".format(test.name))
                casaopts = " --nogui --nologger --log2term"
                casa_exe = exec_path + "/casa"
                assert (test != None)
                cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
                cwd = work_dir + "/" + test.name

                print("Running cmd " + str(cmd) + "in " + cwd)
                if not os.path.exists(cwd):
                    os.makedirs(cwd)
                starttime = datetime.datetime.now()
                output = r.runshell(cmd, test.timeout,cwd)
                endtime = datetime.datetime.now()
                runtime = endtime - starttime
                xunit.append_result(test.name, str(runtime), len(output), output)
                print("")

        elif "mpi" not in test.options and sys.platform != "darwin" and ( pmode == 'parallel' or pmode == 'both'):
            # Special Case when you have tests in the HPC / parallel list but no mpi test option
            print("Running test: {} in Serial mode".format(test.name))
            casaopts = " --nogui --nologger --log2term"
            casa_exe = exec_path + "/casa"
            assert (test != None)
            cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
            cwd = work_dir + "/" + test.name

            print("Running cmd " + str(cmd) + "in " + cwd)
            if not os.path.exists(cwd):
                os.makedirs(cwd)
            starttime = datetime.datetime.now()
            output = r.runshell(cmd, test.timeout,cwd)
            endtime = datetime.datetime.now()
            runtime = endtime - starttime
            xunit.append_result(test.name, str(runtime), len(output), output)
            print("")

        elif pmode == 'serial':
            if "casampi" in test.name:
                continue
            print("Running test: {} in Serial mode".format(test.name))
            casaopts = " --nogui --nologger --log2term"
            casa_exe = exec_path + "/casa"
            assert (test != None)
            cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
            cwd = work_dir + "/" + test.name

            print("Running cmd " + str(cmd) + "in " + cwd)
            if not os.path.exists(cwd):
                os.makedirs(cwd)
            starttime = datetime.datetime.now()
            output = r.runshell(cmd, test.timeout,cwd)
            endtime = datetime.datetime.now()
            runtime = endtime - starttime
            xunit.append_result(test.name, str(runtime), len(output), output)
            print("")

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
########################################            Main-Start-Up            ###########################################
########################################################################################################################

if __name__ == "__main__":

    print("HAVE_MEMTEST: {}".format(HAVE_MEMTEST))
    print("HAVE_COVTEST: {}".format(HAVE_COVTEST))
    print("HAVE_ROBOT: {}".format(HAVE_ROBOT))
    print("HAVE_PYTEST: {}".format(HAVE_PYTEST))
    print("HAVE_NOSE: {}".format(HAVE_NOSE))
    print("")

    # List of tests to run
    testnames = []
    test_paths = []

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument("-i", "--list",action='store_true',help='print the list of tests & tags defined in component_to_test_map.json')
    parser.add_argument("-v", "--verbose",action='store_true',help="Verbose Test Execution")
    parser.add_argument("-x", "--dry-run",action='store_true',help="dry run Test Execution")
    parser.add_argument("-s", "--classes",nargs='+',metavar='test',help='print the classes from a test script') # copy of Dry-Run
    parser.add_argument("-f", "--file",nargs='?', type=argparse.FileType('r'),help='run the tests defined in an ASCII file <list>; one test per line')

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
    parser.add_argument('-m','--pmode',  help='Parallelization mode: serial, parallel, both', required=False)
    parser.add_argument('--bamboo', help='Set Bamboo Flag to True',default=False,action='store_true', required=False)
    parser.add_argument('-r','--rcdir',  help='Casa rcdir', required=False)

    parser.add_argument('--ignore_list',  help='map file of tests to ignore', required=False)

    if not IS_CASA6:
        if "-c" in sys.argv:
            i = sys.argv.index("-c")
        args, unknownArgs = parser.parse_known_args(sys.argv[i+2:])
    else:
        args, unknownArgs = parser.parse_known_args()
    
    print(args)
    print("")

    tests_to_ignore = None
    if args.ignore_list is not None:
        if args.ignore_list.endswith(".json"):
            ignore_test_map = json.load(open(args.ignore_list))
            tests_to_ignore = [x["testScript"].strip() for x in ignore_test_map["testlist"]]
        else:
            tests_to_ignore = [x.strip() for x in args.ignore_list.split(",")]

    print("Operating system: " +  platform.system())
    print("")
    rcdir=""

    if args.rcdir is not None:
        rcdir="--rcdir=" + args.rcdir
    print("rcdir: " + rcdir)

    if args.test_group is not None:
        components = args.test_group
        components = [x.strip() for x in components.split(",")]
        print("Testing Components" + str(components))
        print("")
        
        if not args.bamboo:
            if args.mapfile is not None:
                component_to_test_map = json.load(args.mapfile)
            else:
                try:
                    import casatestutils as _;
                    with open("{}/{}".format(_.__path__[0], "component_to_test_map.json")) as ctt:
                        component_to_test_map = json.load(ctt)
                except:
                    print("No JSON file to Map")
            no_test_components = []
            for c in components:
                _isComponent = False
                component = c.strip()
                for myDict in component_to_test_map["testlist"]:
                    #print(component, myDict["testGroup"])
                    if component in myDict["testGroup"] or component in myDict["testType"]:
                        _isComponent = True
                        if myDict["testScript"] not in testnames: testnames.append(myDict["testScript"])
                if not _isComponent:
                    print("No Tests for Component: {}".format(component))
                    no_test_components.append(component)

            if (len(no_test_components) > 0) and (len(testnames)==0):
                print("No Test Suite for Component(s): {} Using Component 'default'".format(no_test_components))
                component = 'default'
                for myDict in component_to_test_map["testlist"]:
                    if component in myDict["testGroup"]:
                        _isComponent = True
                        testnames.append(myDict["testScript"])

    if args.verbose:
        verbose = True
    if args.list:
        try:
            tmp = {}
            import casatestutils as _;
            with open("{}/{}".format(_.__path__[0], "component_to_test_map.json")) as ctt:
                component_to_test_map = json.load(ctt)
            for myDict in component_to_test_map["testlist"]:
                tmp[myDict["testScript"]] = myDict["testGroup"]
            for key, value in tmp.items():
                print(key,value)
        except:
            list_tests()
        sys.exit(1)

    ## Dry Run
    if args.dry_run or (args.classes is not None):
        DRY_RUN = True

    if args.file is not None:
        #logger.info('Reading Test List from %s: ', args.file)
        for line in args.file:
            try:
                #logger.debug("Adding Test %s from file %s",re.sub(r'[\n\r]+', '',line),args.file)
                testnames.append(re.sub(r'[\n\r]+', '',line))
            except:
                raise Exception(" The list should contain one test per line.")

    if args.branch is not None:
        JIRA_BRANCH = args.branch


    if args.test_paths is not None:
        test_paths = [x.strip() for x in args.test_paths.split(',')]

    temp_storage = []
    for arg in unknownArgs:
        if arg.startswith(("-", "--")):
            raise ValueError('unrecognized argument: %s'%(arg))
            sys.exit()
        else:
            if '[' in arg:
                tests = [x.strip() for x in arg.split("],")]
                for i in range(len(tests)):
                    test = tests[i]
                    if '[' in test and not test.endswith("]"):
                        tests[i] = tests[i] + "]"
                for i in range(len(tests)):
                    test = tests[i]
                    #print(tests)
                    if test.find(",") < test.find('['):
                        temp_storage = temp_storage + test.split(',',1)
                    else:
                        temp_storage.append(test)
                tests = temp_storage
            else:
                tests = [x.strip() for x in arg.split(",")]
            for test in tests:
                try:
                    testcases = None
                    # Check if testcases are provided
                    if "[" in test:
                        listarray = test.split("[")
                        if not listarray[0].endswith(".py"):
                            testnames.append(test)
                        else:
                            test = listarray[0]
                            testcases = listarray[1]
                            real_path = os.path.realpath(test)
                            if ("test_" not in real_path) or  ("test_" not in real_path) or ( os.path.exists(real_path) ==False):
                                print("{} is not a Test File".format(test))
                                continue
                            else:
                                if testcases is not None: real_path = os.path.realpath(test) + "[" + testcases
                                testnames.append(real_path)

                    # Check if test is real path are provided
                    elif test.endswith(".py"):
                        real_path = os.path.realpath(test)
                        if ("test_" not in real_path) or  ("test_" not in real_path) or ( os.path.exists(real_path) ==False):
                            print("{} is not a Test File".format(test))
                            continue
                        else:
                            testnames.append(real_path)

                    # else Assume test exists in bitbucket
                    else:
                        testnames.append(test)
                except:
                    traceback.print_exc()
    try:
        if args.bamboo:
            from testrunner.shell_runner import ShellRunner
            from testrunner import xvfb_helper
            from testrunner.xunit import Xunit
            if args.pkg:
                print("Package: " + args.pkg)
            print("Test configuration file: " + str(args.test_config))
            print("Number of cores: " + str(args.ncores))
            print("Workdir: " + str(args.work_dir))
            print("branch: " + str(args.branch))
            pmodes = ['serial','parallel','both']
            if args.pmode not in pmodes:
                raise Exception("Invalid pmode: '{}'. Valid modes: '{}'".format(args.pmode ,str(pmodes)))


            run_bamboo(args.pkg, args.work_dir, args.branch, args.test_group, args.test_list, test_paths, args.test_config, args.ncores, args.verbose, args.pmode, tests_to_ignore)

        else:
            #If no tests are given, no subet tag or --all option
            #print("Testnames: {}".format(testnames))
            if args.test_paths is not None:
                tests = []
                test_paths = [x.strip() for x in args.test_paths.split(',')]
                if len(testnames) == 0:
                    for test_path in test_paths:
                        for root, dirs, files in os.walk(test_path):
                            for file in files:
                                if file.endswith(".py") and file.startswith("test_"):
                                     tests.append(os.path.join(root, file))
                else:
                    for test_path in test_paths:
                        #print(test_path)
                        for test in testnames:
                            if not test.endswith(".py"): 
                                test = test + ".py"
                            #print(test)
                            for root, dirs, files in os.walk(test_path):
                                for file in files:
                                    if file == test:
                                        tests.append(os.path.join(root, file))
                testnames = tests

            if tests_to_ignore is not None:
                print("\nTests to Ignore: ",tests_to_ignore )
                indices = []
                for i, t in enumerate(testnames):
                     if t.split("/")[-1].replace(".py","") in tests_to_ignore:
                        indices.append(i)
                testnames = [v for i,v in enumerate(testnames) if i not in indices]

            if testnames == [] or len(testnames) == 0:
                print("List of tests is empty")
                parser.print_help(sys.stderr)
                sys.exit(1)
            #print(testnames)
            #sys.exit()
            run(testnames, args.branch, DRY_RUN)
    except:
        traceback.print_exc()


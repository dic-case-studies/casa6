
""" Script to run unit tests from the command line as: 
    casapy [casa-options] -c runUnitTest.py testname1 testname2 ...
    casapy [casa-options] -c runUnitTest.py testname1[test_r,test23] testname2...
    casapy [casa-options] -c runUnitTest.py --Help
    casapy [casa-options] -c runUnitTest.py --list
    casapy [casa-options] -c runUnitTest.py --list
    casapy [casa-options] -c runUnitTest.py --file Tests.txt
    casapy [casa-options] -c runUnitTest.py --classes test_listobs
    
    or from inside casapy:
    runUnitTest.main(['testname']) 
    runUnitTest.main()
    
    NOTE: It will search for scripts in the casapy installation directory, which usually is in:
           <casa_install_dir>/python/2.6/tests"""

# The main class in testwrapper.py is:
# class UnitTest, methods: getUnitTest(), getFuncTest()
#
# UnitTest.getUnitTest() --> for new scripts as described in ....doc....
# UnitTest.getFuncTest() --> for old tests, which contain functions data() and run()


import os
import sys
import getopt
import traceback
import unittest
import string
import re
import shutil
import pprint
import nose
from taskinit import casalog
import casa_stack_manip

##
## testwrapper.py depends upon the current directory being in the path because
## it changes to the directory where the test is located and then imports it.
## CASA no longer leaves empty strings in sys.path to avoid confusion when
## stray files are in the current directory.
##
sys.path.insert(0,'')

PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

CASA_DIR = os.environ["CASAPATH"].split()[0]
TESTS_DIR = CASA_DIR + "/" + os.environ["CASAPATH"].split()[1] + '/lib/python' + PYVER + '/tests/'
#DATA_DIR = CASA_DIR+'/data/'
#print 'HELLOR DATA_DIR'
#print DATA_DIR
if not os.access(TESTS_DIR, os.F_OK):
    if os.access(CASA_DIR+'/lib64', os.F_OK):
        TESTS_DIR = CASA_DIR+'/lib64/python' + PYVER + '/tests/'
    elif os.access(CASA_DIR+'/lib', os.F_OK):
        TESTS_DIR = CASA_DIR+'/lib/python'+ PYVER +'/tests/'
    else:            #Mac release
        TESTS_DIR = CASA_DIR+'/Resources/python/tests/'

HAVE_MEMTEST=True
try:
    import memTest
except:
    HAVE_MEMTEST = False

import testwrapper
from testwrapper import *


RUN_SUBTEST = False

# Tests included in the following file are run automatically by
# Hudson. This is also the list of tests run when no options are given
# to this program
LISTofTESTS = TESTS_DIR+'unittests_list.txt'


# memory mode variable
MEM = 0

def usage():
    print '========================================================================='
    print '\nRunUnitTest will execute Python unit test(s) of CASA tasks.'
    print 'Usage:\n'
    print 'casapy [casapy-options] -c runUnitTest.py [options] test_name\n'
    print 'Options:'
    print '  no option              print this message and exit.'
    print '  -a or --all            run all tests defined in '
    print '                         trunk/gcwrap/python/scripts/tests/unittests_list.txt.'
    print '  <test_name>            run only <test_name> (more tests are separated by spaces).'
    print '  -f or --file <list>    run the tests defined in an ASCII file <list>; one test per line.'
    print '  -d or --datadir <dir>  set an env. variable to a directory, TEST_DATADIR=<dir> '
    print '                         that can be used inside the tests.'
    print '  -m or --mem            show the memory used by the tests and the number of files left open.'
    print '  -g or --debug          set casalog.filter to DEBUG.'
    print '  -l or --list           print the list of tests from '
    print '                         trunk/gcwrap/python/scripts/tests/unittests_list.txt.'
    print '  -s or --classes        print the classes from a test script (those returned by suite()).'
    print '  -H or --Help           print this message and exit.\n'
    print 'NOTE: it will look for tests in the install directory, which usually is \r'
    print '      <casa_install_dir>/python/2.7/tests'
    print 'See documentation in: http://www.eso.org/~scastro/ALMA/CASAUnitTests.htm\n'
    print '=========================================================================='

def list_tests():
    print 'Full list of unit tests'
    print '-----------------------'
    for t in readfile(LISTofTESTS):
        print t
    

def haslist(name):
    '''Check if specific list of tests have been requested'''
    n0 = name.rfind('[')
    n1 = name.rfind(']')
    if n0 == -1:
        return False
    return True

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

def readfile(FILE):
    # It will skip lines that contain '#' and
    # it will only read words starting with test
    if(not os.path.exists(FILE)):
        print 'ERROR: List of tests does not exist'
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

def settestdir(datadir):
    '''Set an environmental variable for the data directory'''
    absdatadir = os.path.abspath(datadir)
    os.environ.__setitem__('TEST_DATADIR',absdatadir)
    return

def getclasses(testnames):
    '''Get the classes of a test script
       It will copy the test script to /tmp
       and remove it afterwards'''

    here = os.getcwd()
    tmpdir = '/tmp'
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
                        print '\t%s'%c(attr)
                  
            os.remove(filename+'.py')       
            os.remove(filename+'.pyc')       
        
        os.chdir(here)
    except:
        print '--> ERROR: Cannot copy script to %s'%tmpdir
        return

def getsubtests(filename,list=[]):
    f = filename
    testlist_to_execute = []
    here = os.getcwd()
    tmpdir = '/tmp'
    os.chdir(tmpdir)
    tt = UnitTest(f)
    tt.copyTest(copyto=tmpdir)
    classes = tt.getTestClasses(f)
    for c in classes:
        # Check if class has @attr(tag=''). Note: class @attr takes priority over func attr
        if 'tag' in c.__dict__:
            if c.tag == ATTR_VAL:
                for attr, value in c.__dict__.iteritems():
                    if len(attr) >= len("test") and attr[:len("test")] == "test":
                        testlist_to_execute.append([attr,value.__module__])
        else:
            # Check if functions within each class has @attr(tag = '') or func.tag = ''
            for attr, value in c.__dict__.iteritems():
                if len(attr) >= len("test") and attr[:len("test")] == "test":
                    if hasattr(value,'tag'):
                        if value.tag == ATTR_VAL:
                          testlist_to_execute.append([attr,value.__module__])
    os.remove(f+'.py')
    os.remove(f+'.pyc')
    os.chdir(here)

    return testlist_to_execute

# Define which tests to run    
whichtests = 0

def main(testnames):

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
            raise Exception, 'List of tests \"%s\" is empty or does not exist'%LISTofTESTS

    elif (type(testnames) != type([])):
        if (os.path.isfile(testnames)):
            # How to prevent it from opening a real test???
            whichtests = 1
            listtests = readfile(testnames)
            if listtests == []:
                raise Exception, 'List of tests is empty'
        else:
            raise Exception, 'List of tests does not exist'
            
    else:
        # run specific tests
        whichtests = 1


    # Directories
    PWD = os.getcwd()
    WDIR = PWD+'/nosedir/'
    
    # Create a working directory
    workdir = WDIR
    print 'Creating work directory '+ workdir
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
    
    print "Starting unit tests for %s: " %(listtests)

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
                    print 'running tests %d to %d' % \
                        (offset, min(offset + chksz, len(testcases)))
                    testcases = testcases[offset:offset + chksz]
                else:
                    testcases = UnitTest(ff).getUnitTest(tests)
                list = list+testcases                

    if RUN_SUBTEST:


        if len(testlist_to_execute) == 0:
            raise ValueError, "Cannot Find Tests with Attribute:'%s'"%(ATTR_VAL)
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
        raise Exception, 'ERROR: There are no valid tests to run'
                                                                     
                
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
        print "Failed to run one or more tests"
        traceback.print_exc()
    else:
        os.chdir(PWD)

def main_with_rethrow(testnames):
    """
    Sets the __rethrow_casa_exceptions switch, and restores its value at the end.
    This can be used to run tests with behavior similar to the casatasks version of tassk
    in CASA6. That is, exceptions are raised normally by all tasks.
    """
    frame = casa_stack_manip.stack_frame_find()
    orig_rethrow = frame.get('__rethrow_casa_exceptions', False)
    try:
        frame['__rethrow_casa_exceptions'] = True
        main(testnames)
    finally:
        frame['__rethrow_casa_exceptions'] = orig_rethrow


# ------------------ NOTE ---------------------------------------------
# Once CASA moves to Python 2.7, the getpopt module should be replaced
# by argparse. The next section will need to be updated accordingly
# ---------------------------------------------------------------------
if __name__ == "__main__":
    ## flush output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)
    
    # Get command line arguments
    if "-c" in sys.argv:
        # If called with ... -c runUnitTest.py from the command line,
        # then parse the command line parameters
        i = sys.argv.index("-c")
        if len(sys.argv) >= i + 2 and \
               re.compile("runUnitTest\.py$").search(sys.argv[i + 1]):
            
        
            try:
                # Get only this script options
                opts,args=getopt.getopt(sys.argv[i+2:], "Halmgs:f:d:r:", ["Help","all","list","mem",
                                                                     "debug","classes=","file=",
                                                                     "datadir=","attr="])
                
            except getopt.GetoptError, err:
                # Print help information and exit:
                print str(err) # will print something like "option -a not recognized"
                usage()
                os._exit(2)
                
            # List of tests to run
            testnames = []
            
            # Boolean for file with tests.
            # One could allow the use of --file with specific tests given in
            # the command line by removing this option and appending to the
            # testnames list in the args handling
            hasfile = False
            alltests = False
            
            #If no option is given, show the Help page
            if opts == [] and args == []:
                usage()
                os._exit(0)
                
            
            # All other options       
            for o, a in opts:
                if o in ("-H", "--Help"):
                    usage()
                    os._exit(0) 
                if o in ("-l", "--list"):
                    list_tests()
                    os._exit(0)
                if o in ("-s", "--classes"): 
                    testnames.append(a)
                    getclasses(testnames)
                    os._exit(0)
                if o in ("-m", "--mem"):
                    # run specific tests in mem mode            
                    MEM = 1
                elif o in ("-g", "--debug"):
                    #Set the casalog to DEBUG
                    casalog.filter('DEBUG')
                elif o in ("-d", "--datadir"):
                    # This will create an environmental variable called
                    # TEST_DATADIR that can be read by the tests to use
                    # an alternative location for the data. This is used 
                    # to test tasks with MMS data
                    # directory with test data
                    datadir = a
                    if not os.path.isdir(datadir):                            
                        raise Exception, 'Value of --datadir is not a directory -> '+datadir  
                    
                    # Set an environmental variable for the data directory
                    settestdir(datadir)
                    if not os.environ.has_key('TEST_DATADIR'):    
                        raise Exception, 'Could not create environmental variable TEST_DATADIR'                        
                        
                elif o in ("-a", "--all"):
                    alltests = True
                    whichtests = 0
                    testnames = 'all'
                    break
                elif o in ("-f", "--file"):
                    hasfile = True
                    testnames = a

                elif o in ("-r", "--attr"):
                    RUN_SUBTEST = True
                    ATTR_VAL = a
                    
                else:
                    assert False, "unhandled option"


            # Deal with other arguments
            if args != [] and not hasfile and not alltests:
                testnames = args
                                        
        else:
            testnames = []
        
    else:
        # Not called with -c (but possibly execfile() from iPython)
        testnames = []

                    
    frame = casa_stack_manip.stack_frame_find()
    orig_rethrow = frame.get('__rethrow_casa_exceptions', False)
    try:
        frame['__rethrow_casa_exceptions'] = True
        main_with_rethrow(testnames)
    except:
        traceback.print_exc()
    finally:
        frame['__rethrow_casa_exceptions'] = orig_rethrow

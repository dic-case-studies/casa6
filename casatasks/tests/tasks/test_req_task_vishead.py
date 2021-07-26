########################################################################
# test_req_task_vishead.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_vishead/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import vishead, casalog
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy
# get is_python3 and is_CASA6
from casatasks.private.casa_transition import *

def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

def chmod_recursive(path, mode):

    os.chmod(path, mode)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root, d), mode)
        for f in files:
            os.chmod(os.path.join(root, f), mode)

if CASA6:
    datapath = casatools.ctsys.resolve('unittest/vishead/')
    casaimagepath = casatools.ctsys.resolve('unittest/vishead/ngc5921.ms')
else:
    datapath = os.path.join(os.environ.get('CASAPATH').split()[0], 'casatestdata/unittest/vishead/')
    casaimagepath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/vishead/ngc5921.ms'

# used in a type comparison
if is_python3:
    numpy_str_ = numpy.str_
else:
    numpy_str_ = numpy.string_

logpath = casalog.logfile()

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/vishead/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print('WARN: directory '+DATADIR+' does not exist')

print('vishead tests will use data from '+datapath)

input_file = 'n4826_16apr98.ms'  # 128 channels
#if testmms:
#    input_file = 'n4826_16apr98.mms'
stars = "*************"
stop_on_first_error = False

# Keeps track of number of passes, failures
# Private class
class tester:
    def __init__(self):
        self.total = 0
        self.fail = 0
        self.current_test = ""

    def start(self, msg):
        self.total += 1
        print()
        print(stars + " Test %s (" % self.total + msg + ") start " + stars)
        self.current_test = msg

    def end(self, condition, error_msg):
        status = "OK"
        if not is_true(condition):
            if is_python3:
                print(error_msg, file=sys.stderr)
            else:
                print >> sys.stderr, error_msg
            self.fail += 1
            status = "FAIL"
            if stop_on_first_error:
                raise Exception("Halt!")
        print(stars + " Test " + self.current_test + " " + status + " " + stars)

    def done(self):
        print("%s/%s tests passed" % (self.total-self.fail, self.total))
        if self.fail > 0:
            raise Exception("%s/%s failures" % (self.fail, self.total))
        else:
            print("All tests passed, congratulations!")


def is_true(x):
    # numpy array comparison yields a list of booleans
    # (not a single boolean). If necessary, convert from
    # list of booleans to single boolean (all elements must be True)
    if type(x) != bool and type(x) != numpy.bool_:
        return False not in x
    else:
        return x

class vishead_test(unittest.TestCase):

    def setUp(self):
        if (os.path.exists(input_file)):
            os.system('rm -rf ' + input_file)

        os.system('cp -RH ' + os.path.join(datapath, input_file) + ' ' + input_file)
        if not is_CASA6:
            default('vishead')

    def tearDown(self):
        casalog.setlogfile(logpath)

        if os.path.exists(input_file):
            os.system('rm -rf ' + input_file)

        if os.path.exists('testlog.log'):
            os.remove('testlog.log')

        if os.path.exists('ngc5921.ms'):
            shutil.rmtree('ngc5921.ms', ignore_errors=True)

    def test_takesMeasurementSet(self):
        ''' 1. test_takesMeasurementSet: Check that vishead takes a MeasurementSet (*.ms)'''
        casalog.setlogfile('testlog.log')
        vishead(vis=casaimagepath, mode='list')
        with open('testlog.log') as logf:
            success = ('SEVERE  vishead::image::open' not in logf.read())
        self.assertTrue(success)

    def test_listMode(self):
        ''' 2. test_listMode: Check that vishead in list mode produces a list in the logger with origin=vishead and as output as a Python array'''
        casalog.setlogfile('testlog.log')
        output = vishead(vis=casaimagepath, mode='list')
        keywords = ["cal_grp", "field", "fld_code", "freq_group_name", "log", "observer", "project", "ptcs", "release_date", \
        "schedule", "schedule_type", "spw_name", "source_name", "telescope"]
        anyWordsPresent = False
        for word in keywords:
            if word in output:
                anyWordsPresent = True
        self.assertTrue(anyWordsPresent)

    def test_summaryMode(self):
        ''' 3. test_summaryMode: Check that vishead in summary mode produces output in the logger with origin=summary'''
        casalog.setlogfile('testlog.log')
        vishead(vis=casaimagepath, mode='summary')
        with open('testlog.log') as logf:
            success = ('Antennas:' in logf.read())
        self.assertTrue(success)

    def test_getModeDefaultIndex(self):
        ''' 4. test_getModeDefaultIndex: Check that vishead in get mode returns the entire array of values for a keyword  if no hdindex is specified '''
        keywordValue = vishead(vis=casaimagepath, mode='get', hdkey='field')
        success = len(keywordValue[0]) == 3
        self.assertTrue(success)

    def test_getModeSpecifiedIndex(self):
        ''' 5. test_getModeSpecifiedIndex: Check that vishead in get mode returns a single value for a keyword if an hdindex is specified '''
        keywordValue = vishead(vis=casaimagepath, mode='get', hdkey='field', hdindex='1')
        success = ('1445+09900002_0' in keywordValue)
        self.assertTrue(success)

    def test_putMode(self):
        ''' 6. test_putMode: Check that vishead in put mode changes the value of a keyword '''
        shutil.copytree(casaimagepath, os.path.join(os.getcwd(), 'ngc5921.ms'))
        # Python 2 vs python 2
        if sys.version_info[0] < 3:
            import stat
            chmod_recursive('ngc5921.ms', stat.S_IRWXU )
        else:
            chmod_recursive('ngc5921.ms', 0o777)
        vishead(vis='ngc5921.ms', mode='put', hdkey='field', hdindex='1', hdvalue='TEST')
        testKeywordValue = vishead(vis='ngc5921.ms', mode='get', hdkey='field', hdindex='1')
        success = ('TEST' in testKeywordValue)
        self.assertTrue(success)

    # Merged test cases

    def test_list(self):
        '''Vishead: List mode'''
        t = tester()

        # os.system('pwd')
        # os.system('find ./pointingtest.ms -type f | xargs cksum | grep OBSERVATION | grep -v svn')

        t.start("vishead( '" + input_file + "', 'list', [])")
        orig_hdr = vishead(input_file, 'list', [])  # default listitems seems to work when
        # run manually, but not from here.
        print("Original header =")
        print(orig_hdr)
        t.end(type(orig_hdr) == type({'key1': 'val1', 'key2': 'val2'})
              and orig_hdr['source_name'][0][2] == 'NGC4826',
              "... is a bad header")

    def test_summary(self):
        '''Vishead: Summary mode'''
        t = tester()
        t.start("summary")
        vishead(input_file, 'summary')
        t.end(True, "summary failed")

    def test_accessors(self):
        '''Vishead: Test put/get modes'''
        t = tester()
        orig_hdr = vishead(input_file, 'list', [])

        # Test the set/get value routines.  All of them
        for keyword in orig_hdr:
            print("List value of %s:" % (keyword), orig_hdr[keyword])

            # Test getting.
            valref = vishead(input_file, mode='get', hdkey=keyword, hdindex='')
            print("Read value:     ", valref)
            # sys.stdout.flush()
            for j in range(2):
                t.start("get " + keyword + "[%d]" % j)
                val = valref[j]

                # numpy array comparison returns a list of booleans
                # therefore we have to manually traverse (sigh...) larger data
                # structures containing numpy arrays, in order to check
                # for equality
                if type(orig_hdr[keyword][j]) is dict:
                    are_equal = (orig_hdr[keyword][j].keys() == val.keys())
                    for k in val.keys():
                        are_equal = (are_equal and is_true(orig_hdr[keyword][j][k] == val[k]))
                        if not is_true(orig_hdr[keyword][j][k] == val[k]):
                            print(orig_hdr[keyword][j][k] == val[k])
                else:
                    are_equal = (val == orig_hdr[keyword][j])
                if hasattr(are_equal, 'all'):
                    are_equal = are_equal.all()
                t.end(are_equal, \
                      "Got " + str(val) + ", expected " + str(orig_hdr[keyword][j]))

            # Test putting.
            # Put does not yet use the ref part of valref.
            val = valref[0]
            if type(val) is dict:
                print(str(keyword) + ' is probably a column ' + \
                      'with variable length arrays, don\'t try to write that')
                # because the task doesn't support it
                continue

            if len(val) == 1:
                if type(val[0]) == numpy_str_:
                    myval = 'the_coolest_' + val[0]
                else:
                    myval = 42.0 + val[0]

                t.start("put/get " + keyword)
                print("New value:      ", myval)
                vishead(input_file, mode='put', hdkey=keyword, hdindex='', hdvalue=myval)

                newval = vishead(input_file, mode='get', hdkey=keyword, hdindex='')[0]
                print("Read new value: ", newval)

                t.end(newval == myval, "Got " + str(newval) + ", expected " + str(myval))
            else:
                # read/write full column
                all_values = vishead(input_file, mode='get', hdkey=keyword)[0]
                vishead(input_file, mode='put', hdkey=keyword, hdindex='', hdvalue=all_values)

                i = 0
                for e in val:
                    if type(e) == numpy_str_:
                        myval = 'the_coolest_' + e
                    else:
                        myval = 42.0 + e

                    t.start("put/get " + keyword + '[' + str(i) + ']')

                    print("New value:      ", myval)
                    vishead(input_file, mode='put', hdkey=keyword, hdindex=str(i),
                            hdvalue=myval)

                    newval = vishead(input_file, mode='get', hdkey=keyword, hdindex=str(i))[0]
                    print("Read new value: ", newval)

                    t.end(newval == myval, "Got " + str(newval) + ", expected " + str(myval))

                    i += 1

            # imhead( input_file, 'put', 'object', val['value'] )

        t.done()

def suite():
    return[vishead_test]

# Main #
if __name__ == '__main__':
    unittest.main()


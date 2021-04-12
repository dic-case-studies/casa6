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
    casaimagepath = casatools.ctsys.resolve('unittest/vishead/ngc5921.ms')
else:
    casaimagepath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/vishead/ngc5921.ms'

logpath = casalog.logfile()

class vishead_test(unittest.TestCase):

    def setUp(self):
        if not CASA6:
            default(vishead)

    def tearDown(self):
        casalog.setlogfile(logpath)
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

def suite():
    return[vishead_test]

# Main #
if __name__ == '__main__':
    unittest.main()


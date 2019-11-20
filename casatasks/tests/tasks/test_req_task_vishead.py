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

if CASA6:
    casaimagepath = casatools.ctsys.resolve('visibilities/vla/ngc5921.ms')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        casaimagepath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/ngc5921.ms'
    else:
        casaimagepath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/ngc5921.ms'

logpath = casalog.logfile()

class vishead_test(unittest.TestCase):

    def setUp(self):
        if not CASA6:
            default(vishead)

    def tearDown(self):
        casalog.setlogfile(logpath)
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')

    def test_takesMeasurementSet(self):
        ''' 1. test_takesMeasurementSet: Check that vishead takes a MeasurementSet (*.ms)'''
        casalog.setlogfile('testlog.log')
        vishead(vis=casaimagepath, mode='list')
        success = ('SEVERE  vishead::image::open' not in open('testlog.log').read())
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
        success = ('Antennas:' in open('testlog.log').read())
        self.assertTrue(success)

    def test_getModeDefaultIndex(self):
        ''' 4. test_getModeDefaultIndex: Check that vishead in get mode returns the entire array of values for a keyword 
        if no hdindex is specified '''
        keywordValue = vishead(vis=casaimagepath, mode='get', hdkey='field')
        success = len(keywordValue[0]) == 3
        self.assertTrue(success)

    def test_getModeSpecifiedIndex(self):
        ''' 5. test_getModeSpecifiedIndex: Check that vishead in get mode returns a single value for a keyword if
        an hdindex is specified '''
        keywordValue = vishead(vis=casaimagepath, mode='get', hdkey='field', hdindex='1')
        success = ('1445+09900002_0' in keywordValue)
        self.assertTrue(success)

    def test_putMode(self):
        ''' 6. test_putMode: Check that vishead in put mode changes the value of a keyword '''
        vishead(vis=casaimagepath, mode='put', hdkey='field', hdindex='1', hdvalue='TEST')
        testKeywordValue = vishead(vis=casaimagepath, mode='get', hdkey='field', hdindex='1')
        success = ('TEST' in testKeywordValue)
        vishead(vis=casaimagepath, mode='put', hdkey='field', hdindex='1', hdvalue='1445+09900002_0')
        self.assertTrue(success)

def suite():
    return[vishead_test]

# Main #
if __name__ == '__main__':
    unittest.main()


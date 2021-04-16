##########################################################################
# test_req_task_listhistory.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_listhistory/about
#
#
##########################################################################
CASA6 = False
try:
    print("Importing CASAtools")
    import casatools
    from casatools import ctsys
    from casatools.platform import bytes2str
    from casatasks import listhistory, casalog
    CASA6 = True
except ImportError:
    print ("Cannot import CASAtools using taskinit")
    import commands
    from __main__ import default
    from tasks import *
    from taskinit import *
    import casa_stack_manip

import sys
import os
import subprocess
import unittest
import shutil

logpath = casalog.logfile()

if CASA6:
    datapath = casatools.ctsys.resolve('unittest/listhistory/Itziar.ms')
    fakepath = casatools.ctsys.resolve('unittest/listhistory/')
    #filepath = casatools.ctsys.resolve('testlog.log')
else:
    dataroot = os.environ.get('CASAPATH').split()[0] + '/'
    datapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/listhistory/Itziar.ms'
    fakepath = dataroot + '/casatestdata/unittest/listhistory/'
    #filepath = 'testlog.log'

# This is for tests that check what the parameter validator does when parameters are
# given wrong types - these don't exercise the task but the parameter validator!
if CASA6:
    validator_exc_type = AssertionError
else:
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
    validator_exc_type = RuntimeError

class listhistory_test(unittest.TestCase):
    
    def setUp(self):
        if not CASA6:
            default(listhistory)
        else:
            pass
    
    def tearDown(self):
        casalog.setlogfile(logpath)
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
    
    def test_takesMS(self):
        '''test takesMS: Check that list history takes a valid MS and refuses incorrect inputs'''
        casalog.setlogfile('testlog.log')
        listhistory(datapath)
        with open('testlog.log') as tlog:
            self.assertFalse('SEVERE' in tlog.read())

        with self.assertRaises(RuntimeError):
            listhistory(fakepath)
        with open('testlog.log') as tlog:
            self.assertTrue('SEVERE' in tlog.read())

        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(validator_exc_type):
                listhistory('fake')
        else:
            listhistory('fake')
            self.assertTrue('Argument vis failed to verify' in open('testlog.log').read())
        
    
    def test_logfile(self):
        '''test logfile: Checks to see that a log file is written and populated'''
        casalog.setlogfile('testlog.log')
        listhistory(datapath)
        with open ('testlog.log') as tlf:
            self.assertTrue('History table entries' in tlf.read())
        
    # Here is the start of the merged test cases
    # -----------------------------------------------
    
    def test_emptyInput(self):
        '''Test 1: Empty input should return False'''
        # CASA5 tasks return False, casatasks throw exceptions
        myms = ''
        if CASA6 or\
           casa_stack_manip.stack_frame_find().get('__rethrow_casa_exceptions', False):
            self.assertRaises(Exception,listhistory,myms)
        else:
            res = listhistory(myms)
            self.assertFalse(res)
            
    def test_returnNone(self):
        '''Test 2: Good input should return None'''
        res = listhistory(datapath)
        self.assertEqual(res,None)
        
    def test_listLen(self):
        '''Test 3: Compare length of reference and new lists'''
        logfile= "testlog.log"

        open(logfile,"w").close( )
        casalog.setlogfile(logfile)
        res = listhistory(datapath)

        refnum = 13
        # In CASA6, this +1 accounts for the following log line (which is not in CASA5):
        # Task listhistory complete. Start time: 2020-10-19 11:33:40.195569 End time: ...
        if CASA6:
            refnum += 1

        # Get only the relevant lines in the logfile, between 'Begin/End Task'
        newfile= "newlisth.log"
        cmd="sed -n \"/Begin Task/,/End Task/p\" %s > %s " %(logfile,newfile)
        print(cmd)
        os.system(cmd)
        logfile = newfile

        if CASA6:
            cmd=['wc', '-l', logfile]
            print(cmd)
            output = bytes2str(subprocess.check_output(cmd))
        else:
            cmd="wc -l %s" % logfile
            print(cmd)
            output=commands.getoutput(cmd)

        num = int(output.split()[0])
        self.assertEqual(refnum,num)
    
def suite():
    return[listhistory_test]

if __name__ == '__main__':
    unittest.main()


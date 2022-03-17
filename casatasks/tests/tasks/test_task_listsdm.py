##########################################################################
# test_task_listsdm.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.information.listsdm.html
#
# test_readsdm: Makes sure a sdm can be opened as well and that non-sdms won't be
# test_logOut: Makes sure that the log file contains the proper information
# test_dictOut: Makes sure that the python dict is returned and contains the proper keys
#
##########################################################################
import os
import unittest
import shutil

import casatools
from casatasks import listsdm, casalog

# Real and fake datapaths
datapath = casatools.ctsys.resolve('unittest/listsdm/TOSR0001_sb1308595_1.55294.83601028935')
falsepath = casatools.ctsys.resolve('unittest/listsdm/')
filepath = casatools.ctsys.resolve('testlog.log')
        
logpath = casalog.logfile()
contained = ['baseband', 'chanwidth', 'end', 'field', 'intent', 'nchan', 'nsubs', 'reffreq', 'source', 'spws', 'start', 'timerange']
# Need to find out what should be displayed in the logger
containedLog = ['Baseband', 'ChWidth', 'FieldName', 'Intent', 'Chans', 'Scan', 'Timerange', 'SpwID', 'Observed from']

class listsdm_test(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def tearDown(self):
        casalog.setlogfile(logpath)
    
    @classmethod
    def tearDownClass(cls):
        if os.path.exists(filepath):
            os.remove(filepath)
    
    def test_readSDM(self):
        '''test readsdm: Makes sure the sdm can be opened without error and fails when looking at fake paths'''
        self.assertTrue(listsdm(sdm=datapath))

        exp_exc = FileNotFoundError
        with self.assertRaises(exp_exc):
            listsdm(sdm=falsepath)

        exp_exc = AssertionError
        with self.assertRaises(exp_exc):
            listsdm(sdm='fake')

    def test_logOut(self):
        '''test logout: Tests to make sure the logfile contains the proper infromation (Requires further documentation)'''
        casalog.setlogfile('testlog.log')
        listsdm(sdm=datapath)
        for item in containedLog:
            self.assertTrue(item in open('testlog.log').read())

    def test_dictOut(self):
        '''test dictout: Tests to make sure that the python dictionary contains all the proper keys'''
        sdmDict = listsdm(sdm=datapath)
        if all (k in sdmDict[1] for k in contained):
            pass
        else:
            self.fail()

if __name__ == '__main__':
    unittest.main()

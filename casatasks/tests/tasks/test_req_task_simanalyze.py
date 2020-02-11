##########################################################################
# test_req_task_simanalyze.py
#
# Copyright (C) 2020
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
# https://open-jira.nrao.edu/browse/CAS-3669
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_simanalyze/about
#
# Test case: requirement
#
##########################################################################
 
 
####    Imports     ####
import os
import unittest
 
CASA6 = False
 
try:
    import casatools # not a good idea inside the casashell...
    import casatasks
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
 
# DATA #
if CASA6:
    dataroot = casatools.ctsys.resolve()
    configpath = casatools.ctsys.resolve((os.path.join
                                          (dataroot, 
                                           'alma/simmos/sma.extended.cfg'))
    imagepath = casatools.ctsys.resolve((os.path.join
                                         (dataroot,
                                          'nrao/VLA/CalModels/3C286_Q.im/')))
else:
    dataroot = os.environ.get('CASAPATH').split()[0]
    configpath = (dataroot + 
                  os.path.join(dataroot, 'alma/simmos/sma.extended.cfg'))
    imagepath = (dataroot + 
                 os.path.join(dataroot,'nrao/VLA/CalModels/3C286_Q.im/'))

logpath = casalog.logfile()

####    Tests     ####
class simanalyze_test(unittest.TestCase):

    ### Set Up
    @classmethod
    def setUpClass(cls):
        '''A class method called before tests in an individual class run'''
        pass
  
    def setUp(self):
        '''Method called to prepare the test fixture.  This is called immediately before calling the test method'''
        if not CASA6:
            default(taskname)
        pass
 
    ### Teardown
    def tearDown(self):
        '''Method called immediately after the test method has been called and the result recorded'''
        pass
 
    @classmethod
    def tearDownClass(cls):
        '''A class method called after tests in an individual class have run'''
        pass
 
    ### Test Cases
    def test_imaging_False_analysis_False(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_False_analysis_True(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_interferometric_only_analysis_False(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_total_power_only_analysis_False(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_interferometric_and_total_power_analysis_False(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 
 
    def test_imaging_True_interferometric_only_analysis_True(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_total_power_only_analysis_True(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_interferometric_and_total_power_analysis_True(self):
        '''test_someFunctionality: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

####    Suite: Required for CASA5     ####
def suite():
    return[simanalyze_test]
  
####    Imports     ####
if __name__ == '__main__':
    unittest.main()

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
    import casatools
    import casatasks
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
 
# DATA #
if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms/')
    datacopy = casatools.ctsys.resolve('uid___X02_X3d737_X1_01_small.ms/')
    calpath = casatools.ctsys.resolve(os.path.join(os.path.dirname(os.path.abspath(casatools.__file__)),'__data__/nrao/VLA/CalModels/3C138_K.im'))
    filepath = casatools.ctsys.resolve('testlog.log')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
    datacopy = 'uid___X02_X3d737_X1_01_small.ms/'
    calpath = os.environ.get('CASAPATH').split()[0] + '/data/nrao/VLA/CalModels/3C138_K.im'
    filepath = 'testlog.log'

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
 
 
 
##########################################################################################################################################################################
##########################################################                  EXAMPLE                    ###################################################################
##########################################################################################################################################################################
##########################################################################
# test_req_task_delmod.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_delmod/about
#
# Test-removefield: Check that only the selected fields are removed (This part is broken)
#
##########################################################################
#
#CASA6 = False
#try:
#    import casatools
#    from casatasks import delmod, rmtables, clearcal, casalog, ft
#    CASA6 = True
#    tb = casatools.table()
#except ImportError:
#    from __main__ import default
#    from tasks import delmod, rmtables, clearcal, casalog, ft
#    from taskinit import tbtool as table
#    tb = table()
#
#import os
#import unittest
#import shutil
#import glob
#from filecmp import dircmp
#
# DATA #
#if CASA6:
#    datapath = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms/')
#    datacopy = casatools.ctsys.resolve('uid___X02_X3d737_X1_01_small.ms/')
#    calpath = casatools.ctsys.resolve(os.path.join(os.path.dirname(os.path.abspath(casatools.__file__)),'__data__/nrao/VLA/CalModels/3C138_K.im'))
#    filepath = casatools.ctsys.resolve('testlog.log')
#else:
#    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'):
#        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
#    else:
#        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
#    datacopy = 'uid___X02_X3d737_X1_01_small.ms/'
#    calpath = os.environ.get('CASAPATH').split()[0] + '/data/nrao/VLA/CalModels/3C138_K.im'
#    filepath = 'testlog.log'
#
#logpath = casalog.logfile()
#
#class delmod_test(unittest.TestCase):
#
#    def setUp(self):
#        shutil.copytree(datapath, datacopy)
#        os.chmod(datacopy, 493)
#        for root, dirs, files in os.walk(datacopy):
#            for d in dirs:
#                os.chmod(os.path.join(root, d), 493)
#            for f in files:
#                os.chmod(os.path.join(root, f), 493)
#        clearcal(datacopy, addmodel=True)
#        if not CASA6:
#            default(delmod)
#
#    def tearDown(self):
#        print('TABLE IS BEING REMOVED')
#        casalog.setlogfile(logpath)
#        rmtables(datacopy)
#        if os.path.exists(filepath):
#            os.remove(filepath)
#
#    def test_removefield(self):
#        ''' test_removefield:Check that the field selection paramter works as intended'''
#        ft(datacopy, model=calpath)
#        delmod(datacopy, field='0')
#        dcmp = dircmp(datacopy, datapath)
#        self.assertTrue(len(dcmp.diff_files) > 0)
#
#def suite():
#    return[delmod_test]
#
#if __name__ == '__main__':
#    unittest.main()

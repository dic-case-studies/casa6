##########################################################################
# test_req_tool_su_getoptinumsize.py
#
# Copyright (C) 2021
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
# [https://open-jira.nrao.edu/browse/CAS-13590]
#
# Based on the requirements listed in CASADocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatools.synthesisutils.html
#
# Test case: requirement
# test single term psf fitting 
# test a wrong nterms (>1) for single term psf fitting
# test a wrong psfcutoff (<0, 1) 
# test mulit-term psf fitting (nterms=2)
# 
##########################################################################
 
 
####    Imports     ####
import os
import sys
import unittest
import copy
import shutil
 
# Example of importing helper functions
from casatestutils import testhelper as th
#from casatestutils.imagerhelpers import TestHelpers
 
is_CASA6 = False
 
try:
    from casatools import synthesisutils
    #from casatasks import ...
    is_CASA6 = True
    su = synthesisutils() 
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    su = casac.synthesisutils() 
 
####    Tests     ####
class su_getoptimumsize_test(unittest.TestCase):
    ### su.getoptimumsize takes a single parameter, size.
    
    ### Set Up
    def setUp(self):
        pass            
    ### Teardown
    def tearDown(self):
        pass 
    ### Test Cases
    def test_default(self):
        '''Test default size'''
        ret = su.getOptimumSize()
        self.assertEqual(ret,100)

    def test_2(self):
        '''Test odd non-optimal number'''
        ret = su.getOptimumSize(size=501)
        self.assertEqual(ret,512)

    def test_3(self):
        '''Test even non-optimal number '''
        ret = su.getOptimumSize(size=510)
        self.assertEqual(ret,512)

####    Suite: Required for CASA5     ####
def suite():
    return[su_getoptimumsize_test]
  
####    Imports     ####
if __name__ == '__main__':
    unittest.main()
 
 

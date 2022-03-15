##########################################################################
# test_task_polfromgain.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.calibration.polfromgain.html
#
#
##########################################################################
import os
import unittest
import shutil
import time
from casatestutils import generate_weblog
from casatestutils import stats_dict
from casatestutils import add_to_dict
from casatestutils.compare import compare_caltables
from casatestutils.compare import compare_dictionaries

import casatools
from casatasks import polfromgain
tb = casatools.table()
import sys
import os

import numpy

# DATA #
datapath = casatools.ctsys.resolve('unittest/polfromgain/ngc5921.ms/')
caldata = casatools.ctsys.resolve('unittest/polfromgain/ngcgain.G0/')
refcal = casatools.ctsys.resolve('unittest/polfromgain/polfromgainCalCompare.cal')
        
datacopy = 'gaincopy.ms'
calpath = 'polfromgainout.cal'
refcopy = 'refcopy.cal'
test_dict = {}
        
class polfromgain_test(unittest.TestCase):
    
    def setUp(self):
        pass
            
        shutil.copytree(datapath, datacopy)
        shutil.copytree(refcal, refcopy)
    
    def tearDown(self):
        shutil.rmtree(datacopy)
        
        if os.path.exists(refcopy):
            shutil.rmtree(refcopy)
        
        if os.path.exists(calpath):
            shutil.rmtree(calpath)
    
    @classmethod
    def tearDownClass(cls):
        generate_weblog("polfromgain", test_dict)
    
    @stats_dict(test_dict)
    def test_dictOutput(self):
        ''' Check that the task takes a MS and outputs a dictionary '''
        result = polfromgain(datacopy, tablein=caldata)
        add_to_dict(self, output=test_dict, dataset=datacopy, resulting_dict = result)
        
        self.assertTrue(type(result) == type({}))
        
    @stats_dict(test_dict)
    def test_calTable(self):
        ''' Check that using the caltable param creates a new cal file '''
        polfromgain(datacopy, tablein=caldata, caltable=calpath)
        cal_exists = os.path.exists(calpath)
        add_to_dict(self, output=test_dict, dataset=datacopy)
        
        self.assertTrue(cal_exists)
        
    @stats_dict(test_dict)
    def test_dictContains(self):
        ''' Check that the dict contains all the correct infromation '''
        add_to_dict(self, output=test_dict, dataset=datacopy)
        refdict = {'1331+30500002_0': {}, '1445+09900002_0': {'Spw0': [1.0, 0.0085134637703004248, 0.025951269611231682, 0.0], 'SpwAve': [1.0, 0.0085134637703004248, 0.025951269611231682, 0.0]}, 'N5921_2': {}}
        result = polfromgain(datacopy, tablein=caldata)
        print(result)
        print(refdict)
        print(compare_dictionaries(result, refdict))
        
        keyCheck = result.keys() == refdict.keys()
        Spw0Check = numpy.all(numpy.isclose(refdict['1445+09900002_0']['Spw0'],result['1445+09900002_0']['Spw0'], rtol = 8e-7, atol=1e-8))
        SpwAveCheck = numpy.all(numpy.isclose(refdict['1445+09900002_0']['SpwAve'],result['1445+09900002_0']['SpwAve'], rtol = 8e-7, atol=1e-8))
        
        if Spw0Check and SpwAveCheck and keyCheck:
            self.assertTrue(True)
        else:
            self.fail(msg='{} != {}'.format(result,refdict))
            
        
        #self.assertTrue(casaTestHelper.compare_dictionaries(result, refdict, rtol = 8e-7))
        
    @stats_dict(test_dict)
    def test_calTableContains(self):
        ''' Check that the contents of the calibration table match the reference file '''
        result = polfromgain(datacopy, tablein=caldata, caltable=calpath)
        add_to_dict(self, output=test_dict, dataset=datacopy, resulting_dict=result)
        
        
        self.assertTrue(compare_caltables(calpath, refcopy))

if __name__ == '__main__':
    unittest.main()

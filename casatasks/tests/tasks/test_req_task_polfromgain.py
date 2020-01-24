##########################################################################
# test_req_task_polfromgain.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_polfromgain/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import polfromgain
    CASA6 = True
    tb = casatools.table()
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

import os
import unittest
import shutil
import time
import casaTestHelper
import numpy

# DATA #
if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/vla/ngc5921.ms/')
    caldata = casatools.ctsys.resolve('caltables/ngcgain.G0/')
    refcal = casatools.ctsys.resolve('caltables/polfromgainCalCompare.cal')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/ngc5921.ms/'
        caldata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ngcgain.G0/'
        refcal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/polfromgainCalCompare.cal'

    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/ngc5921.ms/'
        caldata = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ngcgain.G0/'
        refcal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/polfromgainCalCompare.cal'
        
datacopy = 'gaincopy.ms'
calpath = 'polfromgainout.cal'
refcopy = 'refcopy.cal'
test_dict = {}
        
class polfromgain_test(unittest.TestCase):
    
    
    def setUp(self):
        if not CASA6:
            default(polfromgain)
            
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
        casaTestHelper.generate_weblog("polfromgain", test_dict)
    
    @casaTestHelper.stats_dict(test_dict)
    def test_dictOutput(self):
        ''' Check that the task takes a MS and outputs a dictionary '''
        result = polfromgain(datacopy, tablein=caldata)
        casaTestHelper.add_to_dict(self, output=test_dict, dataset=datacopy, resulting_dict = result)
        
        self.assertTrue(type(result) == type({}))
        
    @casaTestHelper.stats_dict(test_dict)
    def test_calTable(self):
        ''' Check that using the caltable param creates a new cal file '''
        polfromgain(datacopy, tablein=caldata, caltable=calpath)
        cal_exists = os.path.exists(calpath)
        casaTestHelper.add_to_dict(self, output=test_dict, dataset=datacopy)
        
        self.assertTrue(cal_exists)
        
    @casaTestHelper.stats_dict(test_dict)
    def test_dictContains(self):
        ''' Check that the dict contains all the correct infromation '''
        casaTestHelper.add_to_dict(self, output=test_dict, dataset=datacopy)
        refdict = {'1331+30500002_0': {}, '1445+09900002_0': {'Spw0': [1.0, 0.0085134637703004248, 0.025951269611231682, 0.0], 'SpwAve': [1.0, 0.0085134637703004248, 0.025951269611231682, 0.0]}, 'N5921_2': {}}
        result = polfromgain(datacopy, tablein=caldata)
        print(result)
        print(refdict)
        print(casaTestHelper.compare_dictionaries(result, refdict))
        
        keyCheck = result.keys() == refdict.keys()
        Spw0Check = numpy.all(numpy.isclose(refdict['1445+09900002_0']['Spw0'],result['1445+09900002_0']['Spw0'], rtol = 8e-7, atol=1e-8))
        SpwAveCheck = numpy.all(numpy.isclose(refdict['1445+09900002_0']['SpwAve'],result['1445+09900002_0']['SpwAve'], rtol = 8e-7, atol=1e-8))
        
        if Spw0Check and SpwAveCheck and keyCheck:
            self.assertTrue(True)
        else:
            self.fail(msg='{} != {}'.format(result,refdict))
            
        
        #self.assertTrue(casaTestHelper.compare_dictionaries(result, refdict, rtol = 8e-7))
        
    @casaTestHelper.stats_dict(test_dict)
    def test_calTableContains(self):
        ''' Check that the contents of the calibration table match the reference file '''
        result = polfromgain(datacopy, tablein=caldata, caltable=calpath)
        casaTestHelper.add_to_dict(self, output=test_dict, dataset=datacopy, resulting_dict=result)
        
        
        self.assertTrue(casaTestHelper.compare_caltables(calpath, refcopy))
    
    
def suite():
    return[polfromgain_test]

if __name__ == '__main__':
    unittest.main()

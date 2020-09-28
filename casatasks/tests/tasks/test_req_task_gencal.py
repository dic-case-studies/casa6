##########################################################################
# test_req_task_gencal.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_accor/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import gencal, rmtables
    tb = casatools.table()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import filecmp
import testhelper as th

### DATA ###

if CASA6:
    evndata = casatools.ctsys.resolve('visibilities/other/n08c1.ms/')
    vlbadata = casatools.ctsys.resolve('visibilities/vlba/ba123a.ms/')
    vlbacal = casatools.ctsys.resolve('caltables/ba123a.gc/')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        evndata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/other/n08c1.ms/'
        vlbadata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vlba/ba123a.ms/'
        vlbacal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ba123a.gc/'
    else:
        evndata = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/other/n08c1.ms/'
        vlbadata = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vlba/ba123a.ms/'
        vlbacal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ba123a.gc/'


caltab = 'cal.A'
evncopy = 'evn_copy.ms'
vlbacopy = 'vlba_copy.ms'

class gencal_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        shutil.copytree(evndata, evncopy)
        shutil.copytree(vlbadata, vlbacopy)
    
    def setUp(self):
        if not CASA6:
            default(gencal)
            
    def tearDown(self):
        rmtables(caltab)
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(evncopy)
        shutil.rmtree(vlbacopy)
    
    
    def test_gainCurve(self):
        ''' Test calibration table produced when gencal is run on an MS with a GAIN_CURVE table '''
        
        gencal(vis=vlbacopy, caltable=caltab, caltype='gc')
        
        self.assertTrue(os.path.exists(caltab))
        self.assertTrue(th.compTables(caltab, vlbacal, ['WEIGHT']))
        

    def test_noGainCurve(self):
        ''' Test that when gencal is run on an MS with no GAIN_CURVE table it creates no calibration table '''

        try:
            gencal(vis=evncopy, caltable=caltab, caltype='gc')
        except:
            pass
        
        self.assertFalse(os.path.exists(caltab))
        

def suite():
    return[gencal_test]

if __name__ == '__main__':
    unittest.main()

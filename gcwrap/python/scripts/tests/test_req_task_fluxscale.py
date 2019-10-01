##########################################################################
# test_req_task_fluxscale.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_fluxscale/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import fluxscale
    CASA6 = True
    tb = casatools.table()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy as np

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/vla/CalMSwithModel.ms')
    datapath2 = casatools.ctsys.resolve('visibilities/vla/nepModel.ms')
    gCal = casatools.ctsys.resolve('caltables/ModelGcal.G0')
    nepCal = casatools.ctsys.resolve('caltables/nepModel.G0')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/CalMSwithModel.ms/'
        datapath2 = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/nepModel.ms/'
        gCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ModelGcal.G0/'
        nepCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/nepModel.G0/'
        
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/CalMSwithModel.ms/'
        datapath2 = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/nepModel.ms/'
        gCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ModelGcal.G0/'
        nepCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/nepModel.G0/'
        
def getParam(data):
    tb.open(data)
    datamean = np.mean(tb.getcol('CPARAM'))
    tb.close()
    
    return datamean
    
fluxout = 'fluxout.cal'
fluxout2 = 'fluxout2.cal'

listtext = 'listtext.txt'

class fluxscale_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass
        
    def setUp(self):
        
        if not CASA6:
            default(fluxscale)
            
        
    def tearDown(self):
        
        if os.path.exists(fluxout):
            shutil.rmtree(fluxout)
            
        if os.path.exists(fluxout2):
            shutil.rmtree(fluxout2)
            
        if os.path.exists(listtext):
            os.remove(listtext)
        
    @classmethod
    def tearDownClass(cls):
        pass
        
    
    def test_vis(self):
        '''
            test_vis
            ------------------------
            
            Check that the task takes a MS and creates a file
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'])
        
        self.assertTrue(os.path.exists(fluxout))
        
        
    def test_reference(self):
        '''
            test_reference
            ------------------------
            
            Check that changing the referenced field will generate a different output
        '''
        
        res1 = fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'])
        res2 = fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['1'])
        
        self.assertFalse(res1 == res2)
        
        
    def test_transfer(self):
        '''
            test_transfer
            ------------------------
            
            Check that specifying transfer gives the field name to transfer the flux scale to
        '''
        
        res1 = fluxscale(vis=datapath2, caltable=nepCal, fluxtable=fluxout, reference=['0'])
        res2 = fluxscale(vis=datapath2, caltable=nepCal, fluxtable=fluxout2, reference=['0'], transfer=['4'])
        
        self.assertFalse(res1 == res2)
        
        
    def test_listfile(self):
        '''
            test_listfile
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], listfile=listtext)
        
        self.assertTrue(os.path.exists(listtext))
        
        
    def test_append(self):
        '''
            test_append
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], append=False)
        tb.open(fluxout)
        length1 = len(tb.getcol('TIME'))
        tb.close()
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], append=True)
        tb.open(fluxout)
        length2 = len(tb.getcol('TIME'))
        tb.close()
        
        self.assertTrue(length1 < length2)
        
        
    def test_refspwmap(self):
        '''
            test_refspwmap
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], refspwmap=[1,1,1,1])
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'], refspwmap=[2,2,2,2])
        
        res1 = getParam(fluxout)
        res2 = getParam(fluxout2)
        
        self.assertFalse(np.all(res1 == res2))
        
        
    def test_gainthreshold(self):
        '''
            test_gainthreshold
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], gainthreshold=0.01)
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'])
        
        res1 = getParam(fluxout)
        res2 = getParam(fluxout2)
        
        self.assertFalse(np.all(res1 == res2))
        
        
    def test_antenna(self):
        '''
            test_antenna
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], antenna='0~7')
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'])
        
        res1 = getParam(fluxout)
        res2 = getParam(fluxout2)
        
        self.assertFalse(np.all(res1 == res2))
        
        
    def test_timerange(self):
        '''
            test_timerange
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], timerange='04:33:23~05:35:17')
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'])
        
        res1 = getParam(fluxout)
        res2 = getParam(fluxout2)
        
        self.assertFalse(np.all(res1 == res2))
        
        
    def test_scan(self):
        '''
            test_scan
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], scan='0~4')
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'])
        
        res1 = getParam(fluxout)
        res2 = getParam(fluxout2)
        
        self.assertFalse(np.all(res1 == res2))
        
        
    def test_incremental(self):
        '''
            test_incremental
            ------------------------
        '''
        
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], incremental=True)
        fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'])
        
        res1 = getParam(fluxout)
        res2 = getParam(fluxout2)
        
        self.assertFalse(np.all(res1 == res2))
        
        
    def test_fitorder(self):
        '''
            test_fitorder
            ------------------------
        '''
        
        res1 = fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout, reference=['0'], fitorder=40)
        res2 = fluxscale(vis=datapath, caltable=gCal, fluxtable=fluxout2, reference=['0'])
        
        self.assertFalse(np.all(res1['1']['spidx'] == res2['1']['spidx']))
        
    """
    def test_display(self):
        '''
            test_display
            ------------------------
        '''
        
        fluxscale()
        
        tb.open()
        datamean = np.mean(tb.getcol())
        tb.close()
        
        self.assertTrue()"""
        
        
def suite():
    return[fluxscale_test]

if __name__ == "__main__":
    unittest.main()

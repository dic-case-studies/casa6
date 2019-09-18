##########################################################################
# test_req_task_fixplanets.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_fixplanets/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import fixplanets
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
    datapath = casatools.ctsys.resolve('visibilities/vla/gaincaltest2.ms/')
    refpath = casatools.ctsys.resolve('visibilities/alma/nep2-shrunk.ms/')
    ephem = casatools.ctsys.resolve('caltables/Titan_55197-59214dUTC_J2000.tab/')
    fakefile = casatools.ctsys.resolve('text/fake.txt/')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/gaincaltest2.ms/'
        refpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/nep2-shrunk.ms/'
        ephem = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/Titan_55197-59214dUTC_J2000.tab/'
        fakefile = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/text/fake.txt/'
        
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/gaincaltest2.ms/'
        refpath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/nep2-shrunk.ms/'
        ephem = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/Titan_55197-59214dUTC_J2000.tab/'
        fakefile = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/text/fake.txt/'
    
dataSOURCE = datapath + 'SOURCE'
dataFIELD = datapath + 'FIELD'
    
copypath = 'copypath.ms'
    
copySOURCE1 = 'copypath.ms/SOURCE'
copyFIELD1 = 'copypath.ms/FIELD'

copypath2 = 'copypath2.ms'
    
copySOURCE2 = 'copypath2.ms/SOURCE'
copyFIELD2 = 'copypath2.ms/FIELD'

refcopy1 = 'refcopy1.ms'
refcopy2 = 'refcopy2.ms'



def compRows(vis1, vis2, subtable):
    
    tb.open(vis1)
    res1 = tb.getcol(subtable)
    tb.close()
    
    tb.open(vis2)
    res2 = tb.getcol(subtable)
    tb.close()
    
    return np.all(res1 == res2)
    
    
class fixplanets_test(unittest.TestCase):
    
    
    def setUp(self):
        if not CASA6:
            default(fixplanets)
            
        shutil.copytree(datapath, copypath)
        shutil.copytree(datapath, copypath2)
        shutil.copytree(refpath, refcopy1)
        shutil.copytree(refpath, refcopy2)
            
    def tearDown(self):
        if os.path.exists(copypath):
            shutil.rmtree(copypath)
            
        if os.path.exists(copypath2):
            shutil.rmtree(copypath2)
            
        if os.path.exists(refcopy1):
            shutil.rmtree(refcopy1)
            
        if os.path.exists(refcopy2):
            shutil.rmtree(refcopy2)
    
    @classmethod
    def tearDownClass(cls):
        
        if os.path.exists('fixplanetstemp2-copypath.ms'):
            shutil.rmtree('fixplanetstemp2-copypath.ms')
    
    def test_direction(self):
        '''
            test_direction
            ----------------
            
            Check that the direction parameter will change the direction in the souce and field tables
        '''
        
        fixplanets(copypath, direction='J2000 19h30m00 -40d00m00', field='1')
        
        sourcediff = compRows(copySOURCE1, dataSOURCE, 'DIRECTION')
        fielddiff = compRows(copyFIELD1, dataFIELD, 'DELAY_DIR')
        
        self.assertFalse(np.all([sourcediff, fielddiff]) == True)
        self.assertFalse(fixplanets(copypath, field='1', direction='wrongdir'))
        self.assertFalse(fixplanets(copypath, field='1', direction=fakefile))
        
    def test_takesvis(self):
        '''
            test_takesvis
            ---------------
            
            Check that the task takes a MS
        '''
        
        self.assertTrue(fixplanets(copypath))
        
    def test_fixuvw(self):
        '''
            test_fixuvw
            -------------
            
            Check that the uvw coordinates are recalculated
        '''
        
        fixplanets(copypath, field='1', direction=ephem, fixuvw=True)
        fixplanets(copypath2, field='1', direction=ephem, fixuvw=False)
        
        self.assertFalse(compRows(copypath, copypath2, 'UVW') == True)
        
        
    def test_attachephemeris(self):
        '''
            test_attachephemeris
            ----------------------
            
            Check that an ephemeris gets attached to the MS
        '''
        
        
        fixplanets(copypath, field='1', direction=ephem)
        tb.open(copypath+'/FIELD')
        cols = tb.colnames()
        tb.close()
        
        self.assertTrue('EPHEMERIS_ID' in cols)
        
    def test_refant(self):
        '''
            test_refant
            -------------
            
            Check that the refant is set to the provided value if the pointing table is being used
        '''
        
        fixplanets(refcopy1, field='2', refant='2')
        fixplanets(refcopy2, field='2')
        
        source1 = refcopy1 + '/SOURCE'
        source2 = refcopy2 + '/SOURCE'
        
        self.assertFalse(compRows(source1, source2, 'DIRECTION') == True)
        
        
    def test_reftime(self):
        '''
            test_reftime
            --------------
            
            Check that the timestamp provided by reftime is used if the pointing table infromation is being used
        '''
        
        fixplanets(refcopy1, field='2', reftime='4968126851')
        fixplanets(refcopy2, field='2')
        
        source1 = refcopy1 + '/SOURCE'
        source2 = refcopy2 + '/SOURCE'
        
        self.assertFalse(compRows(source1, source2, 'DIRECTION') == True)
        
    def test_median(self):
        '''
            test_median
            -------------
            
            Check that the median value for reftime behaves properly
        '''
        
        fixplanets(refcopy1, field='2', reftime='median')
        fixplanets(refcopy2, field='2')
        
        source1 = refcopy1 + '/SOURCE'
        source2 = refcopy2 + '/SOURCE'
        
        self.assertFalse(compRows(source1, source2, 'DIRECTION') == True)
        
    def test_misspointing(self):
        '''
            test_misspointing
            -------------------
            
            Check failure to find the POINTING table rows depending on antenna and time
        '''
        
        self.assertFalse(fixplanets(copypath, field='1'))
        
    
    
def suite():
    return[fixplanets_test]

if __name__ == "__main__":
    unittest.main()

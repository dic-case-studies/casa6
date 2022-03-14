##########################################################################
# test_task_fixplanets.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.fixplanets.html#
#
##########################################################################
import sys
import os
import unittest
import shutil
import numpy as np

import casatools
from casatasks import fixplanets
tb = casatools.table()

### DATA ###
datapath = casatools.ctsys.resolve('unittest/fixplanets/')

# Input data
msfile = 'gaincaltest2.ms'
reffile = 'nep2-shrunk.ms'
ephem = 'Titan_55197-59214dUTC_J2000.tab'
fakefile = 'fake.txt'

dataSOURCE = datapath + msfile + '/SOURCE'
dataFIELD = datapath + msfile + '/FIELD'
    
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
        if not os.path.exists(copypath):           
            shutil.copytree(os.path.join(datapath, msfile), copypath)
        if not os.path.exists(copypath2):           
            shutil.copytree(os.path.join(datapath, msfile), copypath2)
        if not os.path.exists(refcopy1):           
            shutil.copytree(os.path.join(datapath,reffile), refcopy1)
        if not os.path.exists(refcopy2):           
            shutil.copytree(os.path.join(datapath,reffile), refcopy2)
        if not os.path.exists(ephem):           
            shutil.copytree(os.path.join(datapath,ephem), ephem)
        if not os.path.exists(fakefile):           
            shutil.copyfile(os.path.join(datapath,fakefile), fakefile)
            
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
        shutil.rmtree(ephem, ignore_errors=True)
        os.system('rm -rf ' + fakefile)
    
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
        with self.assertRaises(RuntimeError):
            fixplanets(copypath, field='1', direction='wrongdir')
        with self.assertRaises(RuntimeError):
            fixplanets(copypath, field='1', direction=fakefile)
        
    def test_takesvis(self):
        '''
            test_takesvis
            ---------------
            
            Check that the task takes a MS
        '''
        
        try:
            fixplanets(copypath)
        except Exception as exc:
            self.fail('Unexpected exception: {}'.format(exc))

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
        
        with self.assertRaises(RuntimeError):
            fixplanets(copypath, field='1')

if __name__ == "__main__":
    unittest.main()

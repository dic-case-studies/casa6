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

import casatools
from casatools import ms as mstool
from casatools import msmetadata as msmdtool
from casatasks import fixplanets
CASA6 = True
tb = casatools.table()

import sys
import os
import unittest
import shutil
import numpy as np

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

# data from test_fixplanets
outms = 'uid___A002_X1c6e54_X223-thinned.ms'
inpms = os.path.join(datapath, outms)
outms2 = 'uid___A002_X1c6e54_X223-thinned.mms/'
inpms2 = os.path.join(datapath,outms2)

mymst = mstool()
mymsmdt = msmdtool()



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
        # test_fixplanets setup
        res = None
        shutil.rmtree(outms, ignore_errors=True)
        shutil.copytree(inpms, outms)
        shutil.rmtree(outms2, ignore_errors=True)
        os.system('cp -R ' + inpms2 + ' ' + outms2)

    def tearDown(self):
            if os.path.exists(copypath):
                shutil.rmtree(copypath)

            if os.path.exists(copypath2):
                shutil.rmtree(copypath2)

            if os.path.exists(refcopy1):
                shutil.rmtree(refcopy1)

            if os.path.exists(refcopy2):
                shutil.rmtree(refcopy2)

            if os.path.exists('titan.eml'):
                os.remove('titan.eml')

            if os.path.exists('titan.eml.tab'):
                shutil.rmtree('titan.eml.tab')

            # test_fixplanets teardown
            shutil.rmtree(outms, ignore_errors=True)
            shutil.rmtree(outms2, ignore_errors=True)
    
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

    # Merged test cases from test_fixplanets

    def verify(self, thems, thefield, theref):
        therval = True
        mymsmdt.open(thems)
        thefieldids = mymsmdt.fieldsforname(thefield)
        mymsmdt.close()
        mymst.open(thems)
        thedir = mymst.getfielddirmeas(fieldid=thefieldids[0])
        mymst.close()
        print("Read direction result %s" % thedir)
        if not (thedir['refer'] == theref):
            print("ERROR: reference not as expected: expected %s, got %s" % (theref, thedir['refer']))
            therval = False
        return therval

    def test_standardASDM2011(self):
        '''test1: Does a standard fixplanets work on an MS imported from an ASDM from April 2011'''
        for myms in [outms, outms2]:
            fixplanets(myms, 'Titan', True)

    def test_setDirectionASDM2011(self):
        '''test2: Does the setting of a given direction work on an MS imported from an ASDM from April 2011'''
        for myms in [outms, outms2]:
            fixplanets(myms, 'Titan', False, 'J2000 0h0m0s 0d0m0s')
            self.assertTrue(self.verify(myms, 'Titan', 'J2000'))

    def test_missingRefExpectedError(self):
        '''test3: Does the setting of a given direction with ref !=J2000 and != sol.sys. object give the expected error?'''
        for myms in [outms, outms2]:
            with self.assertRaises(RuntimeError):
                fixplanets(myms, 'Titan', False, 'B1950 0h0m0s 0d0m0s')

    def test_setDirectionWithSolRef(self):
        '''test4: Does the setting of a given direction work with a sol system ref frame?'''
        for myms in [outms, outms2]:
            fixplanets(myms, 'Titan', False, 'SATURN 0h0m0s 0d0m0s')
            self.assertTrue(self.verify(myms, 'Titan', 'SATURN'))

    def test_reftimeASDM2011(self):
        '''test5: Does a standard fixplanets work on an MS imported from an ASDM from April 2011 with parameter reftime'''
        for myms in [outms, outms2]:
            fixplanets(vis=myms, field='Titan', fixuvw=True, reftime='median')

    def test_outOfBoundsReftime(self):
        '''test6: Does a standard fixplanets with out of bounds parameter reftime give the expected error'''
        for myms in [outms, outms2]:
            with self.assertRaises(TypeError):
                fixplanets(vis=myms, field='Titan', fixuvw=True, reftime='2012/07/11/08:41:32')

    def test_incorrectRefTime(self):
        '''test7: Does a standard fixplanets with wrong parameter reftime give the expected error'''
        for myms in [outms, outms2]:
            with self.assertRaises(TypeError):
                fixplanets(vis=myms, field='Titan', fixuvw=True, reftime='MUDIAN')

    def test_withEphemeris(self):
        '''test8: Does a fixplanets with an ephemeris work'''
        for myms in [outms, outms2]:
            fixplanets(vis=myms, field='Titan', fixuvw=True,
                       direction=os.path.join(datapath, 'Titan_55437-56293dUTC.tab'))

            self.assertTrue(os.path.exists(myms + '/FIELD/EPHEM0_Titan.tab'))
            self.assertTrue(self.verify(myms, 'Titan', 'APP'))

    def test_ephemerisMimeFormat(self):
        '''test9: Does a fixplanets with an ephemeris in mime format work'''
        os.system('cp ' + os.path.join(datapath, 'titan.eml') + ' .')
        for myms in [outms, outms2]:
            os.system("rm -rf titan.eml.tab")
            fixplanets(vis=myms, field='Titan', fixuvw=True, direction='titan.eml')

            self.assertTrue(os.path.exists(myms + '/FIELD/EPHEM0_Titan.tab'))
            self.assertTrue(self.verify(myms, 'Titan', 'J2000'))


if __name__ == "__main__":
    unittest.main()

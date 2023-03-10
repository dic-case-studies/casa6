##########################################################################
# test_req_task_importuvfits.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_importuvfits/about
#
# test_make tests to see that a valid ms is created
# test_vlaOldName tests the output for the antnamescheme parameter set to old on VLA data
# test_vlaNewName tests the output for the antnamescheme parameter set to new on VLA data
# test_evlaOldName tests the output for the antnamescheme parameter set to old on EVLA data
# test_evlaNewName tests the output for the antnamescheme parameter set to new on EVLA data
# test_carmaOldName tests the output for the antnamescheme parameter set to old on CARMA data
# test_carmaNewName tests the output for the antnamescheme parameter set to new on CARMA data
# test_weightspectexist tests that the WEIGHT_SPECTRUM column has been created in the MS
# test_weightspectpos tests that all values in the WEIGHT_SPECTRUM column are positive
# test_negToFlag tests that negitive values for WEIGHT in the uvfits file have flipped the FLAG val to TRUE
# test_valicWeight tests that the WEIGHTS values are the sum of the corrispoding WEIGHTS_SPECTRUM values
# test_invalidinput tests that non existing inputs will not be accepted
# test_overwrite tests that existing ms files are not overwritten
#
##########################################################################
CASA6 = False
try:
    import casatools
    from casatools import ctsys
    from casatasks import casalog, importuvfits, exportuvfits, rmtables
    ms = casatools.ms()
    tb = casatools.table()
    qa = casatools.quanta()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    import casa_stack_manip

import gc
import math
import os
import unittest
import traceback
import shutil
import sys
import numpy as np


if CASA6:
    mergedDataRoot = ctsys.resolve('unittest/importuvfits')
    vlapath = ctsys.resolve('unittest/importuvfits/3C219D_CAL.UVFITS')
    path = ctsys.resolve('unittest/importuvfits/refim_Cband.G37line.ms')

    exportuvfits(vis=path, fitsfile='EVLAUV.UVFITS')
    evlapath = ctsys.resolve('EVLAUV.UVFITS')
    carmapath = ctsys.resolve('unittest/importuvfits/mirsplit.UVFITS')
    
    #filepath = ctsys.resolve('EVLAUV.UVFITS')
    
    #testlogpath = ctsys.resolve('testlog.log')
else:
    dataroot = os.environ.get('CASAPATH').split()[0] + '/'
    mergedDataRoot = dataroot + 'casatestdata/unittest/importuvfits'
    vlapath = dataroot + 'casatestdata/unittest/importuvfits/3C219D_CAL.UVFITS'
    carmapath = dataroot + 'casatestdata/unittest/importuvfits/mirsplit.UVFITS'
    exportuvfits(vis=dataroot + 'casatestdata/unittest/importuvfits/refim_Cband.G37line.ms', fitsfile='EVLAUV.UVFITS', overwrite=True)

    evlapath = 'EVLAUV.UVFITS'

logpath = casalog.logfile()

class importuvfits_test(unittest.TestCase):
    # 06/13/2010: This seemed to be the only MS in the regression repo
    # that is a good test of padwithflag.
    inpms = 'cvel/input/ANTEN_sort_hann_for_cvel_reg.ms'

    origms = 'start.ms'  # Just a copy of inpms
    fitsfile = 'hanningsmoothed.UVF'
    msfromfits = 'end.ms'

    records = {}
    need_to_initialize = True  # Do once, at start.
    do_teardown = False  # Do once, after initializing and filling records.

    # Its value here should not really matter.

    def setUp(self):
        if not CASA6:
            default(importuvfits)

    def tearDown(self):
        casalog.setlogfile(logpath)

        rmtables('test_set.ms')
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')

        if os.path.exists('kf.ms'):
            shutil.rmtree('kf.ms')
        if os.path.exists('xyz.uvfits'):
            os.remove('xyz.uvfits')

        if self.do_teardown:
            self.qa.done( )
            shutil.rmtree(self.origms)
            shutil.rmtree(self.msfromfits)
            os.remove(self.fitsfile)
            self.do_teardown = False

    @classmethod
    def tearDownClass(cls):
        os.remove(evlapath)
        gc.collect()

    def test_make(self):
        '''test_make: Tests to make sure a valid ms is created'''
        importuvfits(fitsfile=vlapath, vis='test_set.ms')
        self.assertTrue(os.path.exists('test_set.ms'), msg='the MS is not generated by importuvfits')
        try:
            # check if test_4.ms in a valid ms
            ms.open('test_set.ms')
            ms.close()

        except:
            self.fail('output is not a valid ms')

    # Don't do what these two are doing. We don't want to be testing listobs in this test
    # Try to find a way to do this with as few casa tasks as possible
    def test_vlaOldName(self):
        '''test_vlaOldName: Tests to make sure the antnamescheme old functions on VLA data'''
        importuvfits(fitsfile=vlapath, vis='test_set.ms', antnamescheme='old')
        tb.open('test_set.ms'+'/ANTENNA')
        namelist = tb.getcol("NAME").tolist()
        self.assertFalse('VA01' in namelist)
        self.assertTrue('1' in namelist)
        tb.close()

    def test_vlaNewName(self):
        '''test_vlaNewName: Tests to make sure the antnamescheme new functions on VLA data'''
        importuvfits(fitsfile=vlapath, vis='test_set.ms', antnamescheme='new')
        tb.open('test_set.ms'+'/ANTENNA')
        namelist = tb.getcol("NAME").tolist()
        self.assertTrue('VA01' in namelist)
        tb.close()

    def test_evlaOldName(self):
        '''test_evlaOldName: Tests to make sure the antnamescheme old functions on EVLA/JVLA data'''
        importuvfits(fitsfile=evlapath, vis='test_set.ms', antnamescheme='old')
        tb.open('test_set.ms'+'/ANTENNA')
        namelist = tb.getcol("NAME").tolist()
        self.assertFalse('EA01' in namelist)
        self.assertTrue('1' in namelist)
        tb.close()

    def test_evlaNewName(self):
        '''test_evlaNewName: Tests to make sure the antnamescheme new functions on EVLA/JVLA data'''
        importuvfits(fitsfile=evlapath, vis='test_set.ms', antnamescheme='new')
        tb.open('test_set.ms'+'/ANTENNA')
        namelist = tb.getcol("NAME").tolist()
        self.assertTrue('EA01' in namelist)
        tb.close()

    def test_carmaOldName(self):
        '''test_carmaOldName: Tests to make sure the antnamescheme old functions on CARMA data'''
        importuvfits(fitsfile=carmapath, vis='test_set.ms', antnamescheme='old')
        tb.open('test_set.ms'+'/ANTENNA')
        namelist = tb.getcol("NAME").tolist()
        self.assertFalse('CA1' in namelist)
        self.assertTrue('1' in namelist)
        tb.close()

    def test_carmaNewName(self):
        '''test_carmaNewName: Tests to make sure the antnamescheme new functions on CARMA data'''
        importuvfits(fitsfile=carmapath, vis='test_set.ms', antnamescheme='new')
        tb.open('test_set.ms'+'/ANTENNA')
        namelist = tb.getcol("NAME").tolist()
        self.assertTrue('CA1' in namelist)
        tb.close()

    # Should I come back and test these with all the possible arrays too? VLA/EVLA/CARMA
    def test_weightspectexist(self):
        '''test_weightspectexist: Test for the existence of the WEIGHT_SPECTRUM column in the MS'''
        importuvfits(fitsfile=evlapath, vis='test_set.ms', antnamescheme='old')
        tb.open('test_set.ms')
        try:
            tb.getcol('WEIGHT_SPECTRUM')
        except RuntimeError:
            self.fail('The WEIGHTS_SPECTRUM column is not genereated by importuvfits')
        tb.close()

    def test_weightspectpos(self):
        '''test_weightspecpos: Test that the values in WEIGHT_SPECTRUM are positive'''
        # All values were the same in the WEIGHT_SPECTRUM col so I'm just going to test with the one.
        importuvfits(fitsfile=evlapath, vis='test_set.ms')
        tb.open('test_set.ms')
        valcheck = tb.getcol('WEIGHT_SPECTRUM')[0][0][0]
        self.assertTrue(valcheck >= 0.0, msg='WEIGHT_SPECTRUM is negative. It should always be positive')
        tb.close()

    # Need a way of reading the UV fits file before conversion to ms
    def test_negToFlag(self):
        '''test_negToFlag: Tests that FLAG is set to True when the WEIGHT in the fits file was negative'''
        importuvfits(fitsfile=evlapath, vis='test_set.ms')
        tb.open('test_set.ms')
        testFLAG = tb.getcol('FLAG')
        # Get the index for where FLAG=True. This will be used later when comparing to the values in the fits file.
        flagTrue = np.where(testFLAG == True)
        # For now this just checks that some values were flipped to True.
        self.assertTrue(True in testFLAG, msg='There were no Weight values Flagged (Temp test condition)')
        tb.close()

    # Save this one for later, it's a little more involved than the others
    def test_validweight(self):
        # Still confused on why exactly this is the case or if it is what the documentation was talking about
        '''test_validweight: Tests that WEIGHT vals are the sum of WEIGHT_SPECTRUM values'''
        importuvfits(fitsfile=evlapath, vis='test_set.ms', antnamescheme='new')
        tb.open('test_set.ms')
        weightSpecCol = tb.getcol('WEIGHT_SPECTRUM')
        weightCol = tb.getcol('WEIGHT')

        for i in range(100):
            # Come back and try and increase precision, how many sig figs does it need to be?
            weights = (weightCol[:, i])
            spects = (np.sum(weightSpecCol[:, :, i], axis=1))

            self.assertTrue(np.isclose(weights[0], spects[0], atol=3e-5))
            self.assertTrue(np.isclose(weights[1], spects[1], atol=3e-5))
        tb.close()

    def test_invalidinput(self):
        '''test_invalidinput: Tests to see if the given fits file is valid, or if given an invalid filename'''
        # Try for existing non uvfits file types
        casalog.setlogfile('testlog.log')
        if CASA6 or\
           casa_stack_manip.stack_frame_find().get('__rethrow_casa_exceptions', False):
            if CASA6:
                exc_type = AssertionError
            else:
                exc_type = RuntimeError

            with self.assertRaises(exc_type):
                importuvfits(fitsfile='fake.uvfits', vis='test_set.ms')
        else:
            importuvfits(fitsfile='fake.uvfits', vis='test_set.ms')
            self.assertTrue('failed to verify' in open('testlog.log').read(), msg='Verified a non-existing uvfits')

    # not talked about in the documentation
    def test_overwrite(self):
        '''test_overwrite: Tests to make sure files aren't overwritten'''
        casalog.setlogfile('testlog.log')
        importuvfits(fitsfile=vlapath, vis='test_set.ms')
        with self.assertRaises(RuntimeError):
            importuvfits(fitsfile=vlapath, vis='test_set.ms')
        self.assertTrue('user does not want to remove it.' in open('testlog.log').read(), msg='No warning saying that the file will not overwrite was displayed')

    # Merged test cases from test_importuvfits

    def test_receptor_angle(self):
        """CAS-7081: Test receptor angle is preserved"""
        msname = os.path.join(mergedDataRoot, "uvfits_test.ms")
        self.assertTrue(ms.open(msname), "Input dataset not found")
        uvfits = "xyz.uvfits"
        self.assertTrue(ms.tofits(uvfits), "Failed to write uvfits")
        ms.done()
        feed = "/FEED"
        tb.open(msname + feed)
        rec_ang = "RECEPTOR_ANGLE"
        expec = tb.getcol(rec_ang)
        tb.done()
        importname = "kf.ms"
        importuvfits(fitsfile=uvfits, vis=importname)
        tb.open(importname + feed)
        got = tb.getcol(rec_ang)
        tb.done()
        self.assertTrue(np.max(np.abs(got-expec)) < 1e-7, "Receptor angles not preserved")

def suite():
    return[importuvfits_test]

if __name__ == '__main__':
    unittest.main()

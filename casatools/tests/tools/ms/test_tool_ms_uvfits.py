##########################################################################
# test_tool_ms_uvfits.py
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
# Tests for the UVFITS I/O using the ms tool
#
#
##########################################################################
import os
import numpy
import sys
import shutil
import unittest

import numpy as np

from casatools import ms as mstool
from casatools import msmetadata as msmdtool
from casatools import table as tbtool
from casatools import ctsys, quanta

'''
Unit tests for UVFITS I/O tasks.

Features tested:
  0. Can multiple spws with the same # of channels be exported to UVFITS
     using padwithflags?
  1. When that UVFITS file is read back in, is its data still correct?
'''

datapath = 'unittest/uvfits/'

def check_eq(val, expval, tol=None):
    """Checks that val matches expval within tol."""
    try:
        if tol:
            are_eq = abs(val - expval) < tol
        else:
            are_eq = val == expval
        if hasattr(are_eq, 'all'):
            are_eq = are_eq.all()
        if not are_eq:
            raise ValueError('!=')
    except ValueError:
        raise ValueError("%r != %r" % (val, expval))


class uvfits_test(unittest.TestCase):
    # 06/13/2010: This seemed to be the only MS in the regression repo
    # that is a good test of padwithflag.
#    inpms = 'cvel/input/ANTEN_sort_hann_for_cvel_reg.ms'

#    origms = 'start.ms'               # Just a copy of inpms
    fitsfile = 'hanningsmoothed.UVF'
    msfromfits = 'end.ms'
    
    records = {}
    need_to_initialize = True    # Do once, at start.
    do_teardown        = False   # Do once, after initializing and filling records.
                                 # Its value here should not really matter.

    def setUp(self):
        self.qa = quanta( )
        self.msname = ''
        self.fitsname = ''
        #pass
        #if self.need_to_initialize:
        #    self.initialize()

    #def initialize(self):
        # The realization that need_to_initialize needs to be
        # a class variable more or less came from
        # http://www.gossamer-threads.com/lists/python/dev/776699
    #    self.__class__.need_to_initialize = False

    #    if not os.path.exists(self.origms):
            # Copying is technically unnecessary for split,
            # but self.self.origms is shared by other tests, so making
            # it readonly might break them.
    #        shutil.copytree(datapath + self.inpms, self.origms)

    #    if os.path.exists(self.fitsfile):
    #        os.remove(self.fitsfile)

     #   try:
     #       exportuvfits(self.origms, self.fitsfile, padwithflags=True)
     #       self.records['exported'] = os.path.exists(self.fitsfile)

      #      if self.records['exported']:
      #          importuvfits(self.fitsfile, self.msfromfits)
     #   except Exception, e:
     #       print "Error exporting or importing uv data"
     #       raise e


    def tearDown(self):
        shutil.rmtree(self.msname, ignore_errors=True)
        if os.path.exists(self.fitsname):
            os.system('rm -rf '+self.fitsname)
        if self.do_teardown:
            self.qa.done( )
            shutil.rmtree(self.origms)
            shutil.rmtree(self.msfromfits)
            os.remove(self.fitsfile)
            self.do_teardown = False

    #def test_sts(self):
    #    """Subtables, time avg. without correlation selection"""
    #    self.check_subtables('', [(4, 1)])
        
    #def test_data(self):
    #    """DATA[2],   time avg. without correlation selection"""
    #    check_eq(self.records['']['data'],
    #             numpy.array([[ 0.14428490-0.03145669j],
    #                          [-0.00379944+0.00710297j],
    #                          [-0.00381106-0.00066403j],
    #                          [ 0.14404297-0.04763794j]]),
    #             0.0001)
        
    #def test_wt(self):
    #    """WEIGHT[5], time avg. without correlation selection"""
    #    check_eq(self.records['']['weight'],
    #             numpy.array([143596.34375, 410221.34375,
    #                          122627.1640625, 349320.625]),
    #             1.0)
    
    def test_stokes(self):
        """Verify fix to CAS_4283, uvfits files containing actual Stokes parameters will not be imported"""
        myms = mstool()
        self.msname = "my.ms"
        inpfits = ctsys.resolve(datapath + "1331+305_I.UVFITS")
        self.assertRaises(Exception, myms.fromfits, self.msname, inpfits)

    def test_receptor_angle(self):
        """CAS-7081: Test receptor angle is preserved"""
        myms = mstool()
        inpms = ctsys.resolve(datapath + "uvfits_test.ms")
        self.assertTrue(myms.open(inpms), "Input dataset not found")
        self.fitsname = "xyz.uvfits"
        self.assertTrue(myms.tofits(self.fitsname), "Failed to write uvfits")
        myms.done()
        feed = "/FEED"
        mytb = tbtool()
        mytb.open(inpms + feed)
        rec_ang = "RECEPTOR_ANGLE"
        expec = mytb.getcol(rec_ang)
        mytb.done()
        self.msname = "ke.ms"
        self.assertTrue(myms.fromfits(self.msname, self.fitsname), "Failed uvfits import")
        myms.done()
        mytb.open(self.msname + feed)
        got = mytb.getcol(rec_ang)
        mytb.done()
        self.assertTrue(np.max(np.abs(got-expec)) < 1e-7, "Receptor angles not preserved")

    def test_diameters(self):
        """CAS-5818: Verify bogus dish diameters in AN table are not used but normal algorithm is used instead"""
        myms = mstool()
        inpfits = ctsys.resolve(datapath + "CTR_CHI_TR2.RWYCP-10rows-ANT-DIAMTER-0")
        self.msname = "CAS-5818.ms"
        self.assertTrue(myms.fromfits(self.msname, inpfits), "Failed to import uvfits file")
        myms.done()
        mymd = msmdtool()
        mymd.open(self.msname)
        diam = mymd.antennadiameter(-1)
        mymd.done()
        expec = "25m"
        for i in diam.keys():
            self.assertTrue(self.qa.eq(diam[i], expec), "Unexpected diameter for antenna " + i)

    def test_filename_extensions(self):
        """CAS-7696: Verify we turn off fits filename extension support when necessary"""
        myms = mstool()
        inpfits = ctsys.resolve(datapath + "name10rows+000")
        self.msname = "CAS-7696.ms"
        self.assertTrue(myms.fromfits(self.msname, inpfits), "Failed to import uvfits file")
        myms.done()

    def test_export_overwrite(self):
        """CAS-5492: test the overwrite parameter when exporting MSes to uvfits"""
        myms = mstool()
        inpms = ctsys.resolve(datapath + "uvfits_test.ms")
        myms.open(inpms)
        self.fitsname = "CAS-5492.uvfits"
        self.assertTrue(myms.tofits(self.fitsname))
        # fail because overwrite=False
        self.assertRaises(RuntimeError, myms.tofits, fitsfile=self.fitsname,overwrite=False)
        # succeed because overwrite=True
        self.assertTrue(myms.tofits(self.fitsname, overwrite=True))
        myms.done()
            
    def test_badscan(self):
        """CAS-10054: Tests intermittent incorrect scan number in last row of single-scan dataset"""
        myms = mstool()
        inpfits = ctsys.resolve(datapath + "3c273.fits7")
        self.msname = "ngc4826.tutorial.3c273.7.ms"
        self.assertTrue(myms.fromfits(self.msname, inpfits), "Failed to import uvfits file")
        myms.done()

        mytb = tbtool()
        mytb.open(self.msname)
        scans=mytb.getcol('SCAN_NUMBER')
        mytb.close()

        nrows=len(scans)

        print('Last row has scan='+str(scans[nrows-1])+' ; (should be 1).')
        self.assertFalse(scans[nrows-1]==2, "Last row has wrong scan number: "+str(scans[nrows-1]) )
        # the following verifies that _all_ scan numbers are correct (and lists unique values)
        self.assertTrue(sum(scans==1)==nrows, "Unexpected scan number found: "+str(np.unique(scans)) )

if __name__ == '__main__':
    unittest.main()

import os
import numpy
import sys
import shutil
import unittest

import numpy as np

from casatasks import importuvfits

from casatools import ctsys
from casatools import ms as mstool
from casatools import msmetadata as msmdtool
from casatools import table as tbtool
from casatools import quanta

'''
Unit tests for UVFITS I/O tasks.

Features tested:
  0. Can multiple spws with the same # of channels be exported to UVFITS
     using padwithflags?
  1. When that UVFITS file is read back in, is its data still correct?
'''

datapath = 'regression/unittest/uvfits/'
datapath2 = 'regression/'

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


class importuvfits_test(unittest.TestCase):
    # 06/13/2010: This seemed to be the only MS in the regression repo
    # that is a good test of padwithflag.
    inpms = 'cvel/input/ANTEN_sort_hann_for_cvel_reg.ms'

    origms = 'start.ms'               # Just a copy of inpms
    fitsfile = 'hanningsmoothed.UVF'
    msfromfits = 'end.ms'
    
    records = {}
    need_to_initialize = True    # Do once, at start.
    do_teardown        = False   # Do once, after initializing and filling records.
                                 # Its value here should not really matter.
    
    def setUp(self):
        self.qa = quanta( )
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
    
    def test_receptor_angle(self):
        """CAS-7081: Test receptor angle is preserved"""
        myms = mstool()
        msname = ctsys.resolve(datapath + "uvfits_test.ms")
        self.assertTrue(myms.open(msname), "Input dataset not found")
        uvfits = "xyz.uvfits"
        self.assertTrue(myms.tofits(uvfits), "Failed to write uvfits")
        myms.done()
        feed = "/FEED"
        mytb = tbtool()
        mytb.open(msname + feed)
        rec_ang = "RECEPTOR_ANGLE"
        expec = mytb.getcol(rec_ang)
        mytb.done()
        importname = "kf.ms"
        importuvfits(fitsfile=uvfits, vis=importname)
        mytb.open(importname + feed)
        got = mytb.getcol(rec_ang)
        mytb.done()
        self.assertTrue(np.max(np.abs(got-expec)) < 1e-7, "Receptor angles not preserved")

def suite():
    return [importuvfits_test]        

if __name__ == '__main__':
    unittest.main()

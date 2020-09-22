import os
import numpy
import sys
import shutil
import unittest

import numpy as np

from casatasks import exportuvfits

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


class exportuvfits_test(unittest.TestCase):
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
    
    def test_export_overwrite(self):
        """CAS-5492: test the overwrite parameter when exporting MSes to uvfits"""
        myms = mstool()
        msname = ctsys.resolve(datapath + "uvfits_test.ms")
        myms.open(msname)
        fitsname = "CAS-5492.uvfits"
        self.assertTrue(myms.tofits(fitsname))
        # fail because overwrite=False
        self.assertRaises(Exception, myms.tofits, fitsfile=fitsname, overwrite=False)
        # succeed because overwrite=True
        self.assertTrue(myms.tofits(fitsname, overwrite=True))
        myms.done()
        self.assertRaises(Exception, exportuvfits, msname, fitsname, overwrite=False)
        try:
            exportuvfits(msname, fitsname, overwrite=True)
        except Exception as exc:
            self.fail('Unexpected exception: {}'.format(exc))
            
def suite():
    return [exportuvfits_test]        

if __name__ == '__main__':
    unittest.main()

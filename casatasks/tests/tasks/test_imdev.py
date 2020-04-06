import os
import sys
import shutil
import unittest
import math
import numpy
import numbers

from casatools import image as iatool
from casatools import regionmanager as rgtool
from casatools import ctsys
from casatools import table
from casatasks import imdev

_rg = rgtool( )

#run using
# `which casa` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py --mem test_ia_deviation
#
'''
Unit tests for task ia.deviation().
'''

datapath = ctsys.resolve('regression/unittest/ia_deviation/')

input0 = datapath + "100x100x2.im"
ref0 = datapath + "ref0.im"
ref1 = datapath + "ref1.im"
ref2 = datapath + "ref2.im"
ref3 = datapath + "ref3.im"
ref4 = datapath + "ref4.im"
ref5 = datapath + "ref5.im"
ref6 = datapath + "ref6.im"
ref7 = datapath + "ref7.im"

class imdev_test(unittest.TestCase):

    def _compare(self, resold, resnew, helpstr):
        mytype = type(resold)
        self.assertTrue(mytype == type(resnew), helpstr + ": types differ")
        if mytype == dict:
            for k in resold.keys():
                self._compare(resold[k], resnew[k], helpstr)
        elif mytype == numpy.ndarray:
            oldarray = resold.ravel()
            newarray = resnew.ravel()
            self.assertTrue(
                len(oldarray) == len(newarray),
                helpstr + ": array lengths not equal"
            )
            for i in range(len(oldarray)):
                self._compare(oldarray[i], newarray[i], helpstr)
        elif mytype == str:
            self.assertTrue(
                resold == resnew,
                helpstr + ": string inequality, old = " + resold + ", new = " + resnew
            )
        elif isinstance(resold, numbers.Integral) or mytype == numpy.int32:
            self.assertTrue(
                resold == resnew,
                helpstr + ": integral inequality, old = " + str(resold) + ", new = " + str(resnew)
            )
        elif isinstance(resold, numbers.Real):
            self.assertTrue(
                resold == resnew
                or abs(resnew/resold - 1) < 1e-6,
                helpstr + "float inequality: old = " + str(resold)
                + ", new = " + str(resnew)
            )
        else:
            self.assertTrue(False, "Unhandled type " + str(mytype))

    def setUp(self):
        self.res = None
        self._myia = iatool()
    
    def tearDown(self):
        self._myia.done()
        tb = table( )
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done( )

    def test001(self):
        """Every pixel is a grid point"""
        outname = "out0.im"
        imdev(
            input0, outname, grid=[1,1], xlength="4pix", ylength="4pix",
            stattype="npts", interp="cub",anchor=[0,0], statalg="cl"
        )
        self._myia.open(ref0)
        expec = self._myia.getchunk()
        self._myia.open(outname)
        got = self._myia.getchunk()
        self._myia.done()
        self._compare(got, expec, "imstatimage test 1")
        
    def test_circle(self):
        """test circles work correctly CAS-10296"""
        myia = self._myia
        imagename = "mycirc.im"
        myia.fromshape(imagename, [100, 100])
        bb = myia.getchunk()
        bb[:] = 1
        myia.putchunk(bb)
        myia.done()
        outfile = "mycirc_out.im"
        imdev(
            imagename=imagename,outfile=outfile, xlength="40pix",
            ylength="", stattype="sum", grid=[20,20]
        )
        myia.open(outfile)
        self.assertTrue(
            numpy.isclose(myia.getchunk()[50,50], 1257.0, 1e-7),
            "incorrect grid pixel value"
        )
        myia.done()

def suite():
    return [imdev_test]

if __name__ == '__main__':
    unittest.main()

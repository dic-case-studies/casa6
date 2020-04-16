########################################################################
# test_req_task_imtrans.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
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
# CAS-12700
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_imtrans/about
#
#
##########################################################################

from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import unittest

try:
    from casatools import ctsys, image, table
    from casatasks import imtrans
    _tb = table( )
    ctsys_resolve = ctsys.resolve
    is_CASA6 = True
except ImportError:
    from tasks import *
    from taskinit import *
    import casac
    from __main__ import *
    image = iatool
    # not a local tool
    _tb = tb
    is_CASA6 = False
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        data_root = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'
    else:
        data_root = os.environ.get('CASAPATH').split()[0] + '/casa-data-req'
    def ctsys_resolve(apath):
        return os.path.join(data_root, apath)

datapath = 'image'

good_image = "reorder_in.fits"
cas_2364im = "CAS-2364.im"

def run_imtrans(imagename, outfile, order):
    return imtrans(imagename=imagename, outfile=outfile, order=order)

class imtrans_test(unittest.TestCase):
    
    def setUp(self):
        shutil.copy(ctsys_resolve(os.path.join(datapath,good_image)), good_image)
    
    def tearDown(self):
        os.remove(good_image)
        self.assertTrue(len(_tb.showcache()) == 0)

    def test_exceptions(self):
        """imtrans: Test various exception cases"""
        def testit(imagename, outfile, order):
            # CASA6 tasks always throw exceptions, CASA5 tasks return False
            if is_CASA6:
                self.assertRaises(Exception, run_imtrans, imagename, outfile, order)
            else:
                self.assertFalse(run_imtrans(imagename, outfile, order))

        # blank imagename
        testit("", "blah", "012")
        
        # not enough specified axes
        testit(good_image, "blah", "01")
        testit(good_image, "blah", 10)
        testit(good_image, "blah", ["r", "d"])
        
        # too many specified axes
        testit(good_image, "blah", "0123")
        testit(good_image, "blah", 1230)
        testit(good_image, "blah", ["r", "d", "f", "s"])

        # Bogus axes specification
        testit(good_image, "blah", "123")
        testit(good_image, "blah", ["r", "d", "s"])
        testit(good_image, "blah", ["r", "d", "r"])
        testit(good_image, "blah", 103)
        
    def test_straight_copy(self):
        """No actual transposing"""
        imagename = good_image
        myia = image()
        myia.open(imagename)
        expecteddata = myia.getchunk()
        expectednames = myia.coordsys().names()
        myia.close()
        count = 0
        for order in ["012", 12, ['r', 'd', 'f'], ["righ", "declin", "freq"]]:
            outfile = "straight_copy_" + str(count)
            self.assertTrue(run_imtrans(imagename, outfile, order))
            myia.open(outfile)
            gotdata = myia.getchunk()
            gotnames = myia.coordsys().names()
            myia.done(remove=True)
            self.assertTrue((expecteddata == gotdata).all())
            self.assertTrue(expectednames == gotnames)
            count += 1

    def test_transpose(self):
        """Test transposing"""
        imagename = good_image
        myia = image()
        myia.open(imagename)
        expecteddata = myia.getchunk()
        expectednames = myia.coordsys().names()
        myia.done()
        count = 0
        for order in ["120", 120, ['d', 'f', 'r'], ["declin", "freq", "righ"]]:
            outname = "transpose_" + str(count)
            self.assertTrue(run_imtrans(imagename, outname, order))
            myia.open(outname)
            gotdata = myia.getchunk()
            inshape = expecteddata.shape
            for i in range(inshape[0]):
                for j in range(inshape[1]):
                    for k in range(inshape[2]):
                        self.assertTrue(expecteddata[i][j][k] == gotdata[j][k][i])
            gotnames = myia.coordsys().names()
            myia.done(remove=True)
            self.assertTrue(expectednames[0] == gotnames[2])
            self.assertTrue(expectednames[1] == gotnames[0])
            self.assertTrue(expectednames[2] == gotnames[1])
            count += 1

    def test_cas_2364(self):
        "test CAS-2364 fix"
        shutil.copytree(ctsys_resolve(os.path.join(datapath, cas_2364im)), cas_2364im)
        order = "0132"
        out1 = "blah2.im"
        self.assertTrue(imtrans(imagename=cas_2364im, outfile=out1, order=order))
        myia = image()
        # to verify fix, just open the image. bug was that exception was thrown when opening output from reorder
        myia.open(out1)
        self.assertTrue(myia)
        myia.done(remove=True)
        shutil.rmtree(cas_2364im)

    def test_history(self):
        """Test history records are written"""
        myia = image()
        imagename = "zz.im"
        myia.fromshape(imagename, [10,10,4,10])
        order = "3210"
        outfile = "zz_out.im"
        self.assertTrue(
            imtrans(imagename=imagename, outfile=outfile, order=order)
        )
        myia.open(outfile)
        msgs = myia.history()
        myia.done(remove=True)
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        teststr = "imtrans"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
        myia.open(imagename)
        myia.done(remove=True)

    def test_imageinfo(self):
        """Verify image info is copied"""
        myia = image()
        imname = "imageinfo_test.im"
        myia.fromshape(imname, [10,20,30])
        myia.setbrightnessunit("Jy/beam")
        myia.setrestoringbeam("4arcmin", "3arcmin", "0deg")
        outfile = "imageinfo_test_out.im"
        self.assertTrue(imtrans(imname, outfile, "201"))
        myia.open(outfile)
        self.assertEqual(
            myia.brightnessunit(), "Jy/beam",
            "brightness unit not copied"
        )
        self.assertTrue(
            myia.restoringbeam(), "restoring beam not copied"
        )
        myia.done(remove=True)
        myia.open(imname)
        myia.done(remove=True)

def suite():
    return [imtrans_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

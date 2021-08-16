########################################################################
# test_req_tool_ia_transpose.py
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
    data_root = os.environ.get('CASAPATH').split()[0] + '/casatestdata/'
    def ctsys_resolve(apath):
        return os.path.join(data_root, apath)

datapath = 'unittest/ia_transpose/'
good_image = "reorder_in.fits"
cas_2364im = "CAS-2364.im"

def run_transpose(imagename, outfile, order):
    myia = image()
    myia.open(imagename)
    res = myia.transpose(outfile=outfile, order=order)
    myia.done()
    return res

def run_imtrans(imagename, outfile, order):
    return imtrans(imagename=imagename, outfile=outfile, order=order)

class ia_transpose_test(unittest.TestCase):
    
    def setUp(self):
        shutil.copy(ctsys_resolve(os.path.join(datapath, good_image)), good_image)
    
    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)
        # make sure directory is clean as per verification test requirement
        cwd = os.getcwd()
        for filename in os.listdir(cwd):
            file_path = os.path.join(cwd, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    # CASA 5 tests need this directory
                    if filename != 'xml':
                        shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

    def test_exceptions(self):
        """imtrans: Test various exception cases"""
        #def testit(imagename, outfile, order):
        #    self.assertRaises(Exception, run_transpose, imagename, outfile, order)

        # blank imagename
        self.assertRaises(Exception, run_transpose, "", "blah", "012")
        
        # not enough specified axes
        self.assertRaises(Exception, run_transpose, good_image, "blah", "01")
        self.assertRaises(Exception, run_transpose, good_image, "blah", 10)
        self.assertRaises(Exception, run_transpose, good_image, "blah", ["r", "d"])
        
        # too many specified axes
        self.assertRaises(Exception, run_transpose, good_image, "blah", "0123")
        self.assertRaises(Exception, run_transpose, good_image, "blah", 1230)
        self.assertRaises(Exception, run_transpose, good_image, "blah", ["r", "d", "f", "s"])

        # Bogus axes specification
        self.assertRaises(Exception, run_transpose, good_image, "blah", "123")
        self.assertRaises(Exception, run_transpose, good_image, "blah", ["r", "d", "s"])
        self.assertRaises(Exception, run_transpose, good_image, "blah", ["r", "d", "r"])
        self.assertRaises(Exception, run_transpose, good_image, "blah", 103)
        
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
            newim = run_transpose(imagename, outfile, order)
            gotdata = newim.getchunk()
            gotnames = newim.coordsys().names()
            newim.done()
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
            for outname in ["transpose_" + str(count), ""]:
                newim = run_transpose(imagename, outname, order)
                gotdata = newim.getchunk()
                inshape = expecteddata.shape
                for i in range(inshape[0]):
                    for j in range(inshape[1]):
                        for k in range(inshape[2]):
                            self.assertTrue(expecteddata[i][j][k] == gotdata[j][k][i])
                gotnames = newim.coordsys().names()
                newim.done()
                self.assertTrue(expectednames[0] == gotnames[2])
                self.assertTrue(expectednames[1] == gotnames[0])
                self.assertTrue(expectednames[2] == gotnames[1])
            count += 1

    def test_cas_2364(self):
        "test CAS-2364 fix"
        shutil.copytree(ctsys_resolve(os.path.join(datapath, cas_2364im)), cas_2364im)
        order = "0132"
        out1 = "blahxx.im"
        myia = image()
        myia.open(cas_2364im)
        trans = myia.transpose(out1, order)
        myia.done()
        trans.done()
        self.assertTrue(len(_tb.showcache()) == 0)
        # to verify fix, just open the image. bug was that exception was thrown when opening output from reorder
        myia.open(out1)
        self.assertTrue(myia)
        myia.done()

    def test_history(self):
        """Test history records are written"""
        myia = image()
        myia.fromshape("", [10,10,4,10])
        order = "3210"
        kk = myia.transpose("",order=order)
        myia.done()
        msgs = kk.history()
        kk.done()
        teststr = "ia.transpose"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")    
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
        
    def test_imageinfo(self):
        """Verify image info is copied"""
        myia = image()
        myia.fromshape("",[10,20,30])
        myia.setbrightnessunit("Jy/beam")
        myia.setrestoringbeam("4arcmin", "3arcmin", "0deg")
        kk = myia.transpose("", "201")
        self.assertEqual(
            kk.brightnessunit(), "Jy/beam",
            "brightness unit not copied"
        )
        self.assertTrue(
            kk.restoringbeam(), "restoring beam not copied"
        )
        kk.done()

def suite():
    return [ia_transpose_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

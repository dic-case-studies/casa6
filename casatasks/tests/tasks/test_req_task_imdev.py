##########################################################################
# test_req_task_imdev.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_imdev/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import imdev
    tb = casatools.table()
    ia = casatools.image()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy
from filecmp import dircmp
import numbers

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('image/orion_tfeather.im')
    datapath2 = casatools.ctsys.resolve('image/ngc5921.clean.image')
    stokespath = casatools.ctsys.resolve('image/image_input_processor.im/')
    interppath = casatools.ctsys.resolve('image/f2h_quantile.im')
    #old test data
    oldPath = casatools.ctsys.resolve('image/')
    input0 = oldPath + "100x100x2.im"
    ref0 = oldPath + "ref0.im"


else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/orion_tfeather.im'
        datapath2 = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image'
        stokespath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/image_input_processor.im/'
        interppath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/f2h_quantile.im/'
        # test_imdev data path
        oldPath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/'
        input0 = oldPath + "100x100x2.im"
        ref0 = oldPath + "ref0.im"

        
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/orion_tfeather.im'
        datapath2 = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image'
        stokespath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/image_input_processor.im/'
        interppath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/f2h_quantile.im/'
        # test_imdev data path
        oldPath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/'
        input0 = oldPath + "100x100x2.im"
        ref0 = oldPath + "ref0.im"
        

def imArray(image):

    tb.open(image)
    arrayVal = tb.getcol('map')
    tb.close()

    return arrayVal

output = 'testimage.im'
output2 = 'testimage2.im'
output3 = 'testimage3.im'



class imdev_test(unittest.TestCase):

    
    def setUp(self):
        if not CASA6:
            default(imdev)

        self.res = None
        self._myia = casatools.image()
            
    def tearDown(self):
        if os.path.exists(output):
            shutil.rmtree(output)
            
        if os.path.exists(output2):
            shutil.rmtree(output2)
            
        if os.path.exists(output3):
            shutil.rmtree(output3)
            
        if os.path.exists('testcopy.im'):
            shutil.rmtree('testcopy.im')

        if os.path.exists("mycirc_out.im"):
            shutil.rmtree("mycirc_out.im")

        if os.path.exists("mycirc.im"):
            shutil.rmtree("mycirc.im")

        self._myia.done()
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done()


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
    
    def test_outfile(self):
        ''' Check that the outfile parameter passes the name out the output image to be produced '''

        imdev(imagename=datapath, outfile=output)
        self.assertTrue(os.path.exists(output))
        
    def test_region(self):
        ''' Check that the region parameter selects a different section than the default '''

        imdev(imagename=datapath, outfile=output)
        imdev(imagename=datapath, outfile=output2, region='circle[[5h35m21s, -5d24m12s], 10.0arcsec]')

        origRes = imArray(output)
        finRes = imArray(output2)

        print(numpy.array_equal(origRes, finRes))

        self.assertFalse(numpy.array_equal(origRes, finRes))
        
    def test_box(self):
        ''' Check that the box parameter properly selects a rectangular region '''

        imdev(imagename=datapath, outfile=output)
        imdev(imagename=datapath, outfile=output2, box='0,0,50,50')

        origRes = imArray(output)
        finRes = imArray(output2)

        print(numpy.array_equal(origRes, finRes))

        self.assertFalse(numpy.array_equal(origRes, finRes))
    
    def test_chans(self):
        ''' Check that the chans parameter selects a different channel '''

        imdev(imagename=datapath2, outfile=output)
        imdev(imagename=datapath2, outfile=output2, chans='1')

        origRes = imArray(output)
        finRes = imArray(output2)

        print(origRes.shape, finRes.shape)
        print(numpy.array_equal(origRes, finRes))

        print(datapath2)

        print("OrigRes mean: ", numpy.mean(origRes))
        print("finRes mean: ", numpy.mean(finRes))

        self.assertFalse(numpy.array_equal(origRes, finRes))
        
    def test_stokes(self):
        '''
            test_stokes
            -------------
            
            NOTE: Need to find another data set with stokes options
            Come back to this one
        '''

        imdev(imagename=stokespath, outfile=output)
        imdev(imagename=stokespath, outfile=output2, stokes='I')

        origRes = imArray(output)
        finRes = imArray(output2)

        print(numpy.array_equal(origRes, finRes))

        self.assertFalse(numpy.array_equal(origRes, finRes))
        
        
    def test_mask(self):
        ''' Check that mask selection masks a portion of the original image '''

        datacopy = 'testcopy.im'
        shutil.copytree(datapath2, datacopy)

        imdev(imagename=datapath2, outfile=output)
        imdev(imagename=datapath2, outfile=output2, mask='"testcopy.im">0.1')

        origRes = imArray(output)
        finRes = imArray(output2)

        dcmp = dircmp(output, output2)

        self.assertTrue(os.path.exists(output2))

        self.assertTrue(len(dcmp.diff_files) > 0)
        
        
    def test_overwrite(self):
        ''' Check that the overwrite parameter = True overwrites a file of the existing  and raises no error '''

        imdev(imagename=datapath, outfile=output)
        imdev(imagename=datapath, outfile=output, overwrite=True)
        
        
    def test_grid(self):
        ''' Check that the grid parameter changes the grid spacing  '''

        imdev(imagename=datapath, outfile=output, xlength=10)
        imdev(imagename=datapath, outfile=output2, xlength=10, grid=[2,2])

        origRes = imArray(output)
        finRes = imArray(output2)

        print(numpy.array_equal(origRes, finRes))

        self.assertFalse(numpy.array_equal(origRes, finRes))
        
        
    def test_anchor(self):
        ''' Check that this selects the anchor pixel position '''

        imdev(imagename=datapath, outfile=output, xlength=10, grid=[4,5])
        imdev(imagename=datapath, outfile=output2, xlength=10, grid=[4,5], anchor=[0,0])

        origRes = imArray(output)
        finRes = imArray(output2)

        print(numpy.array_equal(origRes, finRes))

        self.assertFalse(numpy.array_equal(origRes, finRes))
        
        
        
    def test_xlength(self):
        ''' Check that this parameter sets the x coordinate length of the bos, or the diameter of the circle. Cirle is used if ylength is an empty string '''

        imdev(imagename=datapath, outfile=output)
        imdev(imagename=datapath, outfile=output2, xlength=10)
        imdev(imagename=datapath, outfile=output3, xlength=10, ylength=10)

        origRes = imArray(output)
        xlenRes = imArray(output2)
        xylenRes = imArray(output3)

        print(numpy.array_equal(origRes, xlenRes))
        print(numpy.array_equal(origRes, xylenRes))
        
        self.assertFalse(numpy.array_equal(origRes, xlenRes))
        self.assertFalse(numpy.array_equal(origRes, xylenRes))
        
        
    def test_ylength(self):
        ''' Check that this gives the y coordinate length of a box. This returns a different image than the default '''

        imdev(imagename=datapath, outfile=output)
        imdev(imagename=datapath, outfile=output2, ylength=10)
        
        origRes = imArray(output)
        ylenRes = imArray(output2)
        
        self.assertFalse(numpy.array_equal(origRes, ylenRes))
        
        
    def test_interp(self):
        ''' Check that the use of different interpolation algorithms creates different image files '''

        # TODO: Needs work, how to force differing interpolations
        datacopy = 'testcopy.im'
        shutil.copytree(interppath, datacopy)

        imdev(imagename=input0, outfile=output, interp="cubic", xlength='4pix', ylength='4pix', stattype='sigma',
              grid=[3, 3], anchor=[0, 0], statalg='cl')
        imdev(imagename=input0, outfile=output2, interp="linear", xlength='4pix', ylength='4pix', stattype='sigma',
              grid=[3, 3], anchor=[0, 0], statalg='cl')



        # dcmp = dircmp(output, output2)
        print(imArray(output).shape)
        print(imArray(output2).shape)
        print("is equal: ", numpy.array_equal(imArray(output), imArray(output2)))
        print(datapath)

        res1 = imArray(output)
        res2 = imArray(output2)

        self.assertFalse(numpy.array_equal(res1, res2))


    def test_stattype(self):
        '''
            test_stattype
            ----------------
        '''

        imdev(imagename=datapath, outfile=output)
        imdev(imagename=datapath, outfile=output2, stattype='median')
        
        dcmp = dircmp(output, output2)
        
        self.assertTrue(len(dcmp.diff_files) > 0)
        
        
    def test_statalg(self):
        ''' Check that changing the stat alg from classic to chauenet produces a different image '''

        imdev(imagename=datapath, outfile=output, xlength=10)
        imdev(imagename=datapath, outfile=output2, xlength=10, statalg='chauvenet')
        
        dcmp = dircmp(output, output2)
        
        self.assertTrue(len(dcmp.diff_files) > 0)
    
    
    def test_zscore(self):
        ''' Check that using the zscore parameter generates a different image '''

        imdev(imagename=datapath, outfile=output, xlength=10, statalg='chauvenet')
        imdev(imagename=datapath, outfile=output2, xlength=10, statalg='chauvenet', zscore=2)
        
        dcmp = dircmp(output, output2)
        
        self.assertTrue(len(dcmp.diff_files) > 0)
        
        
    def test_maxiter(self):
        ''' Check that using the maxiter parameter generates a different image '''

        imdev(imagename=datapath, outfile=output, xlength=10, statalg='chavenet')
        imdev(imagename=datapath, outfile=output2, xlength=10, statalg='chauvenet', maxiter=2)
        
        dcmp = dircmp(output, output2)
        
        self.assertTrue(len(dcmp.diff_files) > 0)

    # test cases from test_imdev

    def test_allGridPoints(self):
        """Every pixel is a grid point"""
        imdev(
            input0, output, grid=[1, 1], xlength="4pix", ylength="4pix",
            stattype="npts", interp="cub", anchor=[0, 0], statalg="cl"
        )
        self._myia.open(ref0)
        expec = self._myia.getchunk()
        self._myia.open(output)
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
            imagename=imagename, outfile=outfile, xlength="40pix",
            ylength="", stattype="sum", grid=[20, 20]
        )
        myia.open(outfile)
        self.assertTrue(
            numpy.isclose(myia.getchunk()[50, 50], 1257.0, 1e-7),
            "incorrect grid pixel value"
        )
        myia.done()


        
        
        
def suite():
    return[imdev_test]

if __name__ == '__main__':
    unittest.main()

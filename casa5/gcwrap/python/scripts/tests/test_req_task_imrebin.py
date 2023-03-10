##########################################################################
# test_req_task_imrebin.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_imrebin/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import imrebin, casalog
    CASA6 = True
    _tb = casatools.table()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
    ia = iatool()
import sys
import os
import unittest
import shutil
import numpy as np

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('unittest/imrebin/orion_tfeather.im/')
    stokespath = casatools.ctsys.resolve('unittest/imrebin/image_input_processor.im/')
    tb = casatools.table()
    ia = casatools.image()

else:
    datapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/imrebin/orion_tfeather.im/'
    stokespath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/imrebin/image_input_processor.im/'
        
def makeImage():
    
    imagename = "gen.im"
    ia.fromshape(imagename, [10, 10, 10])
    bb = ia.getchunk()
    for i in range(10):
        bb[i,5,:] = i
        bb[i,0:5,:] = i+1
        bb[i,6:10,:] = i+2
    ia.putchunk(bb)
    
    ia.done()
        
    return imagename

def makeDegImage():
    
    imagename = "gendeg.im"
    ia.fromshape(imagename, [10, 10, 1, 1])
    bb = ia.getchunk()
    for i in range(10):
        bb[i,5,:] = i
        bb[i,0:5,:] = i+1
        bb[i,6:10,:] = i+2
    ia.putchunk(bb)
    
    ia.done()

def makeFloatImage():
    
    imagename = "genfloat.im"
    ia.fromshape(imagename, [10, 10, 10])
    bb = ia.getchunk()
    for i in range(10):
        bb[i,5,:] = i+.3
        bb[i,0:5,:] = i+1.3
        bb[i,6:10,:] = i+2.3
    ia.putchunk(bb)
    
    ia.done()
        
    return imagename

def makeCompImage():
    
    imagename = "gencomp.im"
    putArr = np.array([[complex(j,2) for i in range(10)] for j in range(10)])
    ia.fromarray(outfile=imagename, pixels=putArr)
    
    ia.done()
        
    return imagename

# function from merged test case

def alleqnum(x,num,tolerance=0):
    if len(x.shape)==1:
        for i in range(x.shape[0]):
            if not (abs(x[i]-num) < tolerance):
                print("x[",i,"]=", x[i])
                return False
    if len(x.shape)==2:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                if not (abs(x[i][j]-num) < tolerance):
                    print("x[",i,"][",j,"]=", x[i][j])
                    return False
    if len(x.shape)==3:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                for k in range(x.shape[2]):
                    if not (abs(x[i][j][k]-num) < tolerance):
                        print("x[",i,"][",j,"][",k,"]=", x[i][j][k])
                        return False
    if len(x.shape)==4:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                for k in range(x.shape[2]):
                    for l in range(x.shape[3]):
                        if not (abs(x[i][j][k][l]-num) < tolerance):
                            print("x[",i,"][",j,"][",k,"][",l,"]=", x[i][j][k])
                            return False
    if len(x.shape)>4:
        stop('unhandled array shape in alleq')
    return True
        
useImage = 'gen.im'
useFloat = 'genfloat.im'
useComp = 'gencomp.im'
useDeg = 'gendeg.im'

rebinned = 'rebinned.im'
rebinned2 = 'rebinned2.im'

logpath = casalog.logfile()
testlog = 'testlog.log'

teardownList = ['dkfajfas.im', 'dx.im', 'erzvd.im', 'gd.im', 'gk.im', 'kbesd.im',
                'kjfasd.im', 'kyzb5.im', 'maskim', 'outxdkd.im', 'st.im', 'vcsfea.im',
                'xxyy.im', 'aa.im', 'ab.im', 'zz.im', 'zz_out.im']
        
class imrebin_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        
        makeImage()
        makeFloatImage()
        makeCompImage()
        makeDegImage()
        
    def setUp(self):
        if not CASA6:
            self._myia = iatool()
            default(imrebin)
        else:
            self._myia = casatools.image()
    
    def tearDown(self):
        self._myia.done()
        casalog.setlogfile(logpath)
        
        if os.path.exists(testlog):
            os.remove(testlog)
        
        if os.path.exists(rebinned):
            shutil.rmtree(rebinned)
            
        if os.path.exists(rebinned2):
            shutil.rmtree(rebinned2)
    
    @classmethod
    def tearDownClass(cls):
        
        shutil.rmtree(useImage)
        shutil.rmtree(useFloat)
        shutil.rmtree(useComp)
        shutil.rmtree(useDeg)
        for item in teardownList:
            if os.path.exists(item):
                shutil.rmtree(item)
        
    
    def test_newsize(self):
        '''
            test_newsize
            --------------
            
            Check that the new image is downscaled in size by the appropriate factor
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2])
        
        tb.open(useImage)
        regImage = tb.getcol('map')
        tb.close()
        
        tb.open(rebinned)
        binImage = tb.getcol('map')
        tb.close()
        
        sizeReg = np.shape(regImage)
        sizeBin = np.shape(binImage)
        
        print(sizeReg[0], sizeBin[0])
        self.assertTrue(sizeReg[0] / 2 == sizeBin[0])
    
    def test_floatValue(self):
        '''
            test_floatValue
            -----------------
            
            Check that the task supports images with float values
        '''
        
        imrebin(imagename=useFloat, outfile=rebinned, factor=[2,2])
        self.assertTrue(rebinned)
        
    def test_compValue(self):
        '''
            test_compValue
            ----------------
            
            Check that the task support images with complex values
            TODO come back to this one to make sure the complex component isn't being discarded
        '''
        
        imrebin(imagename=useComp, outfile=rebinned, factor=[2,2])
        self.assertTrue(os.path.exists(rebinned))

    def test_outAverage(self):
        '''
            test_outAverage
            -----------------
            
            Check that the Output pixel values are the average of the input pixel values
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2])
        
        tb.open(useImage)
        oneChan = tb.getcol('map')[:,:,0,0]
        tb.close()
        
        tb.open(rebinned)
        oneChanFinal = tb.getcol('map')[:,:,0,0]
        tb.close()
        
        testArr = [[0 for i in range(5)] for j in range(5)]
        
        for i in range(0,np.shape(oneChan)[0],2):
            for j in range(0,np.shape(oneChan)[0],2):
                binned = float(oneChan[i][j] + oneChan[i+1][j] + oneChan[i][j+1] + oneChan[i+1][j+1]) / float(4)
                testArr[int(i/2)][int(j/2)] = binned
                
        self.assertTrue(np.all(testArr == oneChanFinal))
        
        
    def test_polNoRebin(self):
        '''
            test_polNoRebin
            -----------------
            
            Check that the polarization axis cannot be rebinned
        '''
        
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(RuntimeError):
                imrebin(imagename=datapath, outfile=rebinned, factor=[2,2,2])
        else:
            casalog.setlogfile(testlog)
            imrebin(imagename=datapath, outfile=rebinned, factor=[2,2,2])
            self.assertTrue('SEVERE' in open(testlog).read())
            
    def test_factor(self):
        '''
            test_factor
            -------------
            
            Check that the factors array must contain at least one element, and fewer than (or equal to) the number of input image axes.
            All these values must be positive.
        '''
        
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(RuntimeError):
                imrebin(imagename=datapath, outfile=rebinned, factor=[])
            with self.assertRaises(RuntimeError):
                imrebin(imagename=datapath, outfile=rebinned, factor=[2,2,2,2,2])
            with self.assertRaises(RuntimeError):
                imrebin(imagename=datapath, outfile=rebinned, factor=[2,-2])
            if CASA6:
                with self.assertRaises(AssertionError):
                    imrebin(imagename=datapath, outfile=rebinned, factor=[2.2,2.2])
            with self.assertRaises(RuntimeError):
                imrebin(imagename=datapath, outfile=rebinned, factor=[1,1,1])
            
        else:
            casalog.setlogfile(testlog)
            imrebin(imagename=datapath, outfile=rebinned, factor=[])
            self.assertTrue('SEVERE' in open(testlog).read())
            casalog.setlogfile(logpath)
            os.remove(testlog)
            
            casalog.setlogfile(testlog)
            imrebin(imagename=datapath, outfile=rebinned, factor=[2,2,2,2,2])
            self.assertTrue('SEVERE' in open(testlog).read())
            casalog.setlogfile(logpath)
            os.remove(testlog)
            
            casalog.setlogfile(testlog)
            imrebin(imagename=datapath, outfile=rebinned, factor=[2,-2])
            self.assertTrue('SEVERE' in open(testlog).read())
            
            casalog.setlogfile(testlog)
            imrebin(imagename=datapath, outfile=rebinned, factor=[1,1,1])
            self.assertTrue('SEVERE' in open(testlog).read())
            
            imrebin(imagename=datapath, outfile=rebinned, factor=[2,2.2])
            imrebin(imagename=datapath, outfile=rebinned2, factor=[2,2])
        
            tb.open(rebinned)
            imFloat = tb.getcol('map')
            tb.close()
            
            tb.open(rebinned2)
            imInt = tb.getcol('map')
            tb.close()
            
            self.assertTrue(np.all(imFloat == imInt))
        
        
        
    def test_axisRemain(self):
        '''
            test_axisRemain
            -----------------
            
            Check that if the number of elements in the factors array is fewer than the number of axes then the remaining axes are not rebinned
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2])
        imrebin(imagename=useImage, outfile=rebinned2, factor=[2,2,2])
        
        tb.open(rebinned)
        sizeFewer = np.shape(tb.getcol('map'))
        tb.close()
        
        tb.open(rebinned2)
        sizeEqual = np.shape(tb.getcol('map'))
        tb.close()
        
        self.assertTrue(sizeFewer[2] > sizeEqual[2])
        
    def test_crop(self):
        '''
            test_crop
            -----------
            
            Check that crop = True crops off the extra pixels off the end of the axis
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[3,3], crop=True)
        imrebin(imagename=useImage, outfile=rebinned2, factor=[3,3], crop=False)
        
        tb.open(rebinned)
        cropBin = np.shape(tb.getcol('map'))
        tb.close()
        
        tb.open(rebinned2)
        noCropBin = np.shape(tb.getcol('map'))
        tb.close()
        
        self.assertFalse(np.all(cropBin == noCropBin))
        
    def test_dropDeg(self):
        '''
            test_dropDeg
            --------------
            
            Check that degenerate axis are dropped
        '''
        
        imrebin(imagename=useDeg, outfile=rebinned, factor=[2,2], dropdeg=False)
        imrebin(imagename=useDeg, outfile=rebinned2, factor=[2,2], dropdeg=True)
        
        tb.open(rebinned)
        noDrop = np.shape(tb.getcol('map'))
        tb.close()
        
        tb.open(rebinned2)
        withDrop = np.shape(tb.getcol('map'))
        tb.close()
        
        self.assertTrue(noDrop != withDrop)
        
        
        
    def test_region(self):
        '''
            test_region
            -------------
            
            Check that the region parameter selects the region to be rebinned
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[3,3], region='centerbox[[5pix,5pix],[3pix,3pix]]')
        
        tb.open(rebinned)
        regSelect = tb.getcol('map')
        tb.close()
        
        self.assertTrue(np.all(regSelect == 6))
        
    def test_box(self):
        '''
            test_box
            ----------
            
            Check that the box parameter properly selects a subset of data
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2], box='0,0,2,2')
        
        tb.open(rebinned)
        boxSelect = tb.getcol('map')
        tb.close()
        
        self.assertTrue(np.all(boxSelect == 1.5))
        
    def test_chans(self):
        '''
            test_chans
            ------------
            
            Check that the channel selection paramter properly selects a subset of the data
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2], chans='0')
        
        tb.open(rebinned)
        selectChan = np.shape(tb.getcol('map'))
        tb.close()
        
        self.assertTrue(selectChan[2] == 1)
        
        
    def test_stokes(self):
        '''
            test_stokes
            -------------
            
            Check that the stokes selection parameter properly selects a subset of the data
        '''
        
        imrebin(imagename=stokespath, outfile=rebinned, factor=[2,2], stokes='i')
        
        tb.open(rebinned)
        outSize = np.shape(tb.getcol('map'))
        tb.close()
    
        self.assertTrue(outSize[3] == 1)
        
    def test_mask(self):
        '''
            test_mask
            -----------
            
            Check that the mask parameter masks the correct areas based on the selection made
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2], mask='gen.im>2')
        
        tb.open(rebinned)
        maskCheck = tb.getcol('map')[:,:,0,0]
        tb.close()
        
        print(maskCheck[0][0],)
        self.assertTrue(maskCheck[0][0] == maskCheck[0][1] == maskCheck[0][2] == 0)
        
    def test_overwrite(self):
        '''
            test_overwrite
            ----------------
            
            Check that overwrite = True is required to overwrite the output file
        '''
        
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2])
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(RuntimeError):
                imrebin(imagename=useImage, outfile=rebinned, factor=[2,2])
        else:
            casalog.setlogfile(testlog)
            imrebin(imagename=useImage, outfile=rebinned, factor=[2,2])
            self.assertTrue('SEVERE' in open(testlog).read())
            
        imrebin(imagename=useImage, outfile=rebinned, factor=[2,2], overwrite=True)
        self.assertTrue(os.path.exists(rebinned))

    # ------ merged test cases ------

    def test_stretch(self):
        """ ia.rebin(): Test stretch parameter"""
        yy = self._myia
        mymask = "maskim"
        yy.fromshape(mymask, [200, 200, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [200, 200, 1, 10]
        imagename = "aa.im"
        yy.fromshape(imagename, shape)
        yy.addnoise()
        yy.done()
        outfile = "ab.im"
        self.assertRaises(
            Exception, imrebin, imagename=imagename, outfile=outfile, factor=[2, 2, 1, 1],
            mask=mymask + ">0", stretch=False, overwrite=True
        )
        imrebin(
            imagename=imagename, outfile=outfile, factor=[2, 2, 1, 1],
            mask=mymask + ">0", stretch=True, overwrite=True
        )
        yy.open(outfile)
        self.assertTrue((yy.shape() == [100, 100, 1, 10]).all())
        yy.done()

    def test_general(self):
        """ ia.rebin(): General tests"""
        # tests moved from imagetest_regression.py and modified

        myia = self._myia
        shp2 = [20, 40]
        d2 = myia.makearray(1.0, [shp2[0], shp2[1]])
        #
        imagename = "st.im"
        myim2 = myia.newimagefromarray(outfile=imagename, pixels=d2)
        myim2.done()
        outfile = "gk.im"
        self.assertRaises(
            Exception, imrebin, imagename=imagename, outfile=outfile,
            factor=[-100, 2], overwrite=True
        )
        imrebin(
            imagename=imagename, outfile=outfile, overwrite=True,
            factor=[2, 2]
        )
        myia.open(outfile)
        p = myia.getchunk()
        self.assertTrue(alleqnum(p, 1.0, tolerance=0.0001))
        myia.done()

    def test_multibeam(self):
        """Test multiple beams"""
        myia = self._myia
        imagename = "gd.im"
        myia.fromshape(imagename, [10, 10, 10])
        myia.setrestoringbeam(
            major="4arcsec", minor="2arcsec", pa="0deg",
            channel=0, polarization=0
        )
        outfile = "dx.im"
        imrebin(
            imagename=imagename, outfile=outfile,
            factor=[2, 2, 1]
        )

        self.assertRaises(
            Exception, imrebin, imagename=imagename, outfile=outfile,
            factor=[2, 2, 2]
        )

    def test_crop(self):
        """Test crop parameter"""
        myia = self._myia
        imagename = "xxyy.im"
        myia.fromshape(imagename, [20, 20, 20])
        factor = [3, 3, 3]
        myia.done()
        outfile = "outxdkd.im"
        imrebin(imagename=imagename, outfile=outfile, factor=factor, crop=True)
        myia.open(outfile)
        self.assertTrue((myia.shape() == [6, 6, 6]).all())
        myia.done()
        imrebin(imagename=imagename, outfile=outfile, factor=factor, crop=False, overwrite=True)
        myia.open(outfile)
        self.assertTrue((myia.shape() == [7, 7, 7]).all())
        myia.done()

    def test_dropdeg(self):
        """Test dropdeg parameter"""
        myia = self._myia
        imagename = "kjfasd.im"
        myia.fromshape(imagename, [20, 20, 1])
        factor = [5, 5]
        myia.done()
        outfile = "dkfajfas.im"
        imrebin(imagename=imagename, outfile=outfile, factor=factor, dropdeg=True)
        myia.open(outfile)
        self.assertTrue((myia.shape() == [4, 4]).all())
        myia.done()

    def test_box(self):
        """Test use of box"""
        myia = self._myia
        imagename = "erzvd.im"
        myia.fromshape(imagename, [30, 30, 1])
        factor = [5, 5]
        myia.done()
        outfile = "vcsfea.im"
        imrebin(imagename=imagename, outfile=outfile, factor=factor, box="5,5,25,25", crop=True)
        myia.open(outfile)
        self.assertTrue((myia.shape() == [4, 4, 1]).all())
        myia.done()

    def test_dropdeg2(self):
        """ axes that become degenerate when regridded are dropped if dropdeg=True: CAS-5836"""
        myia = self._myia
        imagename = "kbesd.im"
        myia.fromshape(imagename, [20, 20, 20])
        factor = [1, 1, 20]
        myia.done()
        outfile = "kyzb5.im"
        imrebin(
            imagename=imagename, outfile=outfile,
            factor=factor, dropdeg=True
        )
        myia.open(outfile)
        self.assertTrue((myia.shape() == [20, 20]).all())
        myia.done()

    def test_history(self):
        """Test history writing"""
        myia = self._myia
        imagename = "zz.im"
        factor = [1, 1, 20]
        myia.fromshape(imagename, [20, 20, 20])
        myia.done()
        outfile = "zz_out.im"
        imrebin(imagename=imagename, outfile=outfile, factor=factor)
        myia.open(outfile)
        msgs = myia.history()
        myia.done()
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        teststr = "imrebin"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
    
    
    
def suite():
    return[imrebin_test]

if __name__ == '__main__':
    unittest.main()

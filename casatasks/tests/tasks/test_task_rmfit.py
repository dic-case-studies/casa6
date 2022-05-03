##########################################################################
# test_task_rmfit.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.rmfit.html
#
#
##########################################################################
import sys
import os
import numpy
import unittest
import shutil
from filecmp import dircmp
import math
import casatools
from casatasks import rmfit
myia = casatools.image()
tb = casatools.table()
mypo = casatools.imagepol()
myia = casatools.image()
ctsys_resolve = casatools.ctsys.resolve


## DATA ## 
casaim = ctsys_resolve('unittest/rmfit/ngc5921.clean.image')
eq_beams = ctsys_resolve('unittest/rmfit/pol_eq_beams.fits')
neq_beams = ctsys_resolve('unittest/rmfit/pol_neq_beams.fits')
outfile = 'out.im'

def table_comp(im1, im2):
    tb.open(im1)
    table1 = tb.getcol('map')
    tb.close()
    
    tb.open(im2)
    table2= tb.getcol('map')
    tb.close()
    
    return((table1 == table2).all())
    
    
class rmfit_test(unittest.TestCase):
    
    def setUp(self):
        myia.fromshape(outfile, [20, 20, 4, 20])
        myia.addnoise()
        myia.done()
    
    def tearDown(self):
        mypo.done()
        for f in (
            outfile, 'rm.im', 'out2.im', 'rm2.im', 'rm1.im',
            'rm_input.im', 'xx.im', 'yy.im'
        ):
            if os.path.exists(f):
                shutil.rmtree(f)
    
    def test_makesImage(self):
        '''
            test_makesImage
            ------------------
            
            This test checks that a rotation measure image is generated when the task is run.
        '''
        
        rmfit(imagename=outfile, rm='rm.im')
        self.assertTrue(os.path.exists('rm.im'))
        
    def test_needsQUV(self):
        '''
            test_needsQUV
            ----------------
            
            This test checks that if the image provided doesn't have Stokes Q, U, or V the task will fail to execute
        '''
        
        with self.assertRaises(RuntimeError):
            rmfit(imagename=casaim, rm='rm.im')

    def test_multiImage(self):
        '''
            test_multiImage
            ------------------
            
            This test checks that rmfit can take multiple image names
        '''
        
        outfile2 = 'out2.im'
        myia.fromshape(outfile2, [20, 20, 4, 20])
        myia.addnoise()
        #csys = myia.coordsys()
        #refval = csys.referencevalue()['numeric']
        #refval[3] = 1.5e9
        #csys.setreferencevalue(refval)
        #myia.setcoordsys(csys.torecord())
        myia.done()
        rmfit(imagename=[outfile, outfile2], rm='rm.im')
        
        self.assertTrue(os.path.exists('rm.im'))
        
    
    def test_sigmaParam(self):
        '''
            test_sigmaParam
            -----------------
            
            Test that the sigma parameter gives the thermal noise in Stokes U and Q
        '''
        
        try:
            rmfit(imagename=outfile, rm='rm.im', sigma=3)
        except Exception:
            self.fail()
        rmfit(imagename=outfile, rm='rm2.im')
        
        self.assertFalse(table_comp('rm.im', 'rm2.im'))
        
    def test_maxpaerr(self):
        '''
            test_maxpaerr
            ---------------
            
            Test that the maxpaerr parameter changes the max allowed position angle
        '''

        try:
            rmfit(imagename=outfile, rm='rm.im', maxpaerr=1)
        except Exception:
            self.fail()
        rmfit(imagename=outfile, rm='rm2.im')
        
        self.assertFalse(table_comp('rm.im', 'rm2.im'))
        
    def test_rmfgparam(self):
        '''
            test_rmfgparam
            -----------------
            
            Test that rmfg is used to specify a foreground rmvalue which is subtracted
        '''
        
        try:
            rmfit(imagename=outfile, rm='rm.im', rmfg=1)
        except Exception:
            self.fail()
        rmfit(imagename=outfile, rm='rm2.im')
        
        self.assertFalse(table_comp('rm.im', 'rm2.im'))
        
    def test_freqAmount(self):
        '''
            test_freqAmount
            -----------------
            
            This test checks that rmfit will only work if the provided image has more than two frequency channels
        '''
        
        outfile2 = 'out2.im'
        myia.fromshape(outfile2, [20,20,4,1])
        myia.addnoise()
        myia.done()
        
        try:
            rmfit(imagename=outfile, rm='rm.im')
        except Exception:
            self.fail()

        with self.assertRaises(RuntimeError):
            rmfit(imagename=outfile2, rm ='rm2.im')
            
    def test_axisDiff(self):
        '''
            test_axisDiff
            ---------------
            
            This tests checks that the task will still run with multiple images with differing spectal axis dimensions
        '''
        
        outfile2 = 'out2.im'
        myia.fromshape(outfile2, [20,20,4,5])
        myia.addnoise()
        myia.done()

        try:
            rmfit(imagename=[outfile, outfile2], rm='rm.im')
        except Exception:
            self.fail()
 
    def test_rmmax(self):
        '''
            test_rmmax
            ------------
            
            This test checks that the rmmax param gives the max absolute RM value that should be solved for
        '''
        
        #TODO I'm not sure exactly how to do this one. I'm unsure what effects that this task is supposed to have on the image
        # come back to this one
        
        rmfit(imagename=outfile, rm='rm.im', rmmax=100000000, sigma=0.001)
        rmfit(imagename=outfile, rm='rm2.im', rmmax=10.0, sigma=0.001)
        
        tb.open('rm.im')
        rmCol = tb.getcol('map')
        tb.close()
        
        tb.open('rm2.im')
        rm2Col = tb.getcol('map')
        tb.close()
        
        self.assertFalse(numpy.all(rmCol == rm2Col))
        
    # Merged in test cases
    def test_rmfit_basics(self):
        """Sanity tests for task rmfit"""
        #myia = iatool()
        outfile = "xx.im"
        myia.fromshape(outfile, [20, 20, 4, 20])
        myia.addnoise()
        myia.done()
        myrm = "rm1.im"
        try:
            rmfit(imagename=outfile, rm=myrm)
        except Exception:
            self.fail()
        myia.open(myrm)
        self.assertTrue((myia.shape() == [20, 20]).all())
        got1 = myia.statistics(list=True, verbose=True)['sumsq']
        myia.done()
     
        # test concatenation of images
        outfile = "yy.im"
        myia.fromshape(outfile, [20, 20, 4, 20])
        myia.addnoise()
        csys = myia.coordsys()
        refval = csys.referencevalue()['numeric']
        refval[3] = 1.5e9
        csys.setreferencevalue(refval)
        myia.setcoordsys(csys.torecord())
        myia.done()
        images = ["xx.im", "yy.im"]
        myrm = "rm2.im"
        try:
            rmfit(imagename=images, rm=myrm)
        except Exception:
            self.fail()
        myia.open(myrm)
        self.assertTrue((myia.shape() == [20, 20]).all())
        got2 = myia.statistics(list=True, verbose=True)['sumsq']
        myia.done()
        self.assertTrue(abs(got1 - got2) > 0.1)
        tb.done()
        self.assertTrue(len(tb.showcache()) == 0)

    def test_algorithm(self):
        """Test rotation measure computation algorithm"""
        #myia = iatool()
        imagename = "rm_input.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        csys = myia.coordsys()
        incr = csys.increment()['numeric']
        incr[3] = 1000*incr[3]
        csys.setincrement(incr)
        myia.setcoordsys(csys.torecord())
        pixvals = myia.getchunk()
        # U values all 1
        U = 1
        pixvals[:,:,2,:] = U
        c = 29979245800.0/100
        RM = 9.6
        pa0deg = 22.5
        pa0 = pa0deg/180*math.pi
        for chan in range(myia.shape()[3]):
            freq = myia.toworld([0,0,0,chan])['numeric'][3]
            lam = c/freq
            Q = U/math.tan(2*(pa0 + RM*lam*lam))
            pixvals[:,:,1,chan] = Q
        myia.putchunk(pixvals)
        myia.done()
        rmim = "rm.im"
        pa0im = "pa0.im"
        sigma = 10e-8
        rmfit(imagename=imagename, rm=rmim, pa0=pa0im, sigma=sigma)
        myia.open(rmim)
        stats = myia.statistics(list=True, verbose=True)
        self.assertTrue((abs(stats['min'][0] - RM)) < 1e-4)
        self.assertTrue((abs(stats['max'][0] - RM)) < 1e-4)
        myia.done(remove=True)
        myia.open(pa0im)
        stats = myia.statistics(list=True, verbose=True)
        self.assertTrue((abs(stats['min'][0] - pa0deg)) < 1e-4)
        self.assertTrue((abs(stats['max'][0] - pa0deg)) < 1e-4)
        myia.done(remove=True)
        tb.done()
        self.assertTrue(len(tb.showcache()) == 0)

if __name__ == '__main__':
    unittest.main()
        
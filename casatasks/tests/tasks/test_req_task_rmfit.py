##########################################################################
# test_req_task_rmfit.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_rmfit/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import rmfit
    CASA6 = True
    myia = casatools.image()
    tb = casatools.table()
    # TODO: Check if the toolname is the same in casa5 and casa6
    mypo = casatools.imagepol()
    myia = casatools.image()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    # not sure what the tool is in casa5 (need to check)
    mypo = imagepol()
    myia = iatool()
import sys
import os
import numpy
import unittest
import shutil
from filecmp import dircmp
import math


## DATA ## 

if CASA6:
    casaim = casatools.ctsys.resolve('unittest/rmfit/ngc5921.clean.image/')
    eq_beams = casatools.ctsys.resolve('unittest/rmfit/pol_eq_beams.fits')
    neq_beams = casatools.ctsys.resolve('unittest/rmfit/pol_neq_beams.fits')
else:
    casaim = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/rmfit/ngc5921.clean.image/'
    eq_beams = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/rmfit/pol_eq_beams.fits'
    neq_beams = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/rmfit/pol_eq_beams.fits'
    
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
        
        if not CASA6:
            default(rmfit)
        
        myia.fromshape(outfile, [20, 20, 4, 20])
        myia.addnoise()
        myia.done()
    
    def tearDown(self):
        
        mypo.done()
        
        shutil.rmtree(outfile)
        
        if os.path.exists('rm.im'):
            shutil.rmtree('rm.im')
        
        if os.path.exists('out2.im'):
            shutil.rmtree('out2.im')
            
        if os.path.exists('rm2.im'):
            shutil.rmtree('rm2.im')
            
        if os.path.exists('rm1.im'):
            shutil.rmtree('rm1.im')
            
        if os.path.exists('rm_input.im'):
            shutil.rmtree('rm_input.im')
            
        if os.path.exists('xx.im'):
            shutil.rmtree('xx.im')
            
        if os.path.exists('yy.im'):
            shutil.rmtree('yy.im')
    
    
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
        
        if CASA6:
            with self.assertRaises(RuntimeError):
                rmfit(imagename=casaim, rm='rm.im')
        else:
            self.assertTrue(True != rmfit(imagename=casaim, rm='rm.im'))
        
        
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
            
        if CASA6:
            with self.assertRaises(RuntimeError):
                rmfit(imagename=outfile2, rm ='rm2.im')
        else:
            self.assertTrue(True != rmfit(imagename=outfile2, rm='rm2.im'))
            
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
        
        
        
def suite():
    return[rmfit_test]

if __name__ == '__main__':
    unittest.main()
        

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
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    myia = iatool()
import sys
import os
import numpy
import unittest
import shutil
from filecmp import dircmp


## DATA ## 

if CASA6:
    casaim = casatools.ctsys.resolve('image/ngc5921.clean.image/')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        casaim = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image/'
    else:
        casaim = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image/'
    
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
        
        shutil.rmtree(outfile)
        
        if os.path.exists('rm.im'):
            shutil.rmtree('rm.im')
        
        if os.path.exists('out2.im'):
            shutil.rmtree('out2.im')
            
        if os.path.exists('rm2.im'):
            shutil.rmtree('rm2.im')
    
    
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
        
        self.assertTrue(rmfit(imagename=outfile, rm='rm.im', sigma=3))
        rmfit(imagename=outfile, rm='rm2.im')
        
        self.assertFalse(table_comp('rm.im', 'rm2.im'))
        
    def test_maxpaerr(self):
        '''
            test_maxpaerr
            ---------------
            
            Test that the maxpaerr parameter changes the max allowed position angle
        '''
        
        self.assertTrue(rmfit(imagename=outfile, rm='rm.im', maxpaerr=1))
        rmfit(imagename=outfile, rm='rm2.im')
        
        self.assertFalse(table_comp('rm.im', 'rm2.im'))
        
    def test_rmfgparam(self):
        '''
            test_rmfgparam
            -----------------
            
            Test that rmfg is used to specify a foreground rmvalue which is subtracted
        '''
        
        self.assertTrue(rmfit(imagename=outfile, rm='rm.im', rmfg=1))
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
        
        self.assertTrue(rmfit(imagename=outfile, rm='rm.im'))
        
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
    
        self.assertTrue(rmfit(imagename=[outfile, outfile2], rm='rm.im'))
        
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
        
        
        
def suite():
    return[rmfit_test]

if __name__ == '__main__':
    unittest.main()
        

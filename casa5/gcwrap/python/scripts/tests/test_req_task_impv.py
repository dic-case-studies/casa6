##########################################################################
# test_req_task_impv.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_impv/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import impv, casalog
    CASA6 = True
except ImportError:
    from __main__ import *
    from tasks import *
    from taskinit import *
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
    
import sys
import os
import unittest
import shutil
import numpy as np
from filecmp import dircmp

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('image/ngc5921.clean.image')
    qa = casatools.quanta()
    mytb = casatools.table()
    myia = casatools.image()

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image'
        myia = ia
        mytb = tb
        
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image'
        myia = ia
        mytb = tb
        
testfile = 'testing.im'
testfile2 = 'testing2.im'
testfile3 = 'testing3.im'

logpath = casalog.logfile()
logname = 'testlog.log'

def makeImage():
    
    imagename = "gen.im"
    myia.fromshape(imagename, [10, 10, 10])
    bb = myia.getchunk()
    # basic sanity test, no rotation involved
    for i in range(10):
        bb[i,5,:] = i
        bb[i,0:5,:] = i+1
        bb[i,6:10,:] = i+2
    myia.putchunk(bb)
    expeccoord = myia.toworld([1,5,0])['numeric'][2]
    mycsys = myia.coordsys()
    units = mycsys.units()
    expinc = mycsys.increment()["numeric"]
    expinc = [
        abs(
            qa.convert(
                qa.quantity(expinc[0], units[0]), "arcsec"
            )["value"]
        ),
        expinc[2]
    ]
    myia.done()
        
    return imagename

class impv_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def setUp(self):
        if not CASA6:
            default(impv)
            
    def tearDown(self):
        if os.path.exists(testfile):
            shutil.rmtree(testfile)
            
        if os.path.exists(testfile2):
            shutil.rmtree(testfile2)
            
        if os.path.exists(testfile3):
            shutil.rmtree(testfile3)
            
        if os.path.exists('gen.im'):
            shutil.rmtree('gen.im')
            
        casalog.setlogfile(logpath)
        
        if os.path.exists(logname):
            os.remove(logname)
    
    @classmethod
    def tearDownClass(cls):
        pass
    
    def test_createsImage(self):
        '''
            test_createsImage
            -------------------
            
            Check that an output image name is generated if an outfile is given
        '''
        
        impv(imagename=datapath, outfile=testfile, start=[10,15], end=[110,120])
        self.assertTrue(os.path.exists(testfile))
        
        # Check that the outfile needs to be given
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(UnboundLocalError):
                impv(imagename=datapath, outfile='', start=[10,15], end=[110,120])
        else:
            casalog.setlogfile(logname)
            impv(imagename=datapath, outfile='', start=[10,15], end=[110,120])
            
            self.assertTrue('SEVERE' in open(logname).read())
        
        
    def test_modelength(self):
        '''
            test_modelength
            -----------
            
            If mode='coords' use start and end values, if 'length' use center, length, and pa
            The use of parameters for the other mode is not allowed
        '''
        # Catch if using the wrong combo of length and mode are allowed
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(UnboundLocalError):
                impv(imagename=datapath, outfile=testfile, start=[10,15], end=[110,120], mode='length')
        else:
            casalog.setlogfile(logname)
            impv(imagename=datapath, outfile=testfile, start=[10,15], end=[110,120], mode='length')
            self.assertTrue('SEVERE' in open(logname).read())
            
        self.assertFalse(os.path.exists(testfile))
        
    def test_modeStartEnd(self):
        ''' Check that start and end is required '''
        # Catch if there is no start or end given
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(UnboundLocalError):
                impv(imagename=datapath, outfile=testfile, mode='coords')
        else:
            casalog.setlogfile(logname)
            impv(imagename=datapath, outfile=testfile, mode='coords')
            self.assertTrue('SEVERE' in open(logname).read())
            
    def test_modeUnsupported(self):
        ''' Check unsupported mode values '''
        impv(imagename=datapath, outfile=testfile, start=[10,15], end=[110,120], mode='c')

        self.assertTrue(os.path.exists(testfile))

        
    def test_region(self):
        '''
            test_region
            -------------
            
            Check that the region parameter sets the spectral extent of the final image
            If this is left blank then the entire spectral range is used.
        '''
        
        impv(imagename=datapath, outfile=testfile, start=[10,15], end=[110,120], region='box[[0pix,0pix],[255pix,255pix]]')
        self.assertTrue(os.path.exists(testfile))
        
        
    def test_startEnd(self):
        '''
            test_startEnd
            ---------------
            
            Check that the start and end parameters can be represented multiple ways
            TODO come back to this for the alternate coord systems selections
        '''
        
        imagename = makeImage()
        
        impv(imagename=imagename, outfile=testfile, start=[2,5], end=[7,5])
        impv(imagename=imagename, outfile=testfile2, start=["3.00000038arcmin", "0'"], end=["2.15980000e+04'", "0arcmin"])
        impv(imagename=imagename, outfile=testfile3, start=["0h0m12s", "0d0m0s" ], end=["2.15980000e+04'", "0arcmin"])
        
        self.assertTrue(os.path.exists(testfile))
        self.assertTrue(os.path.exists(testfile2))
        self.assertTrue(os.path.exists(testfile3))
        
        
    def test_center(self):
        '''
            test_center
            -------------
            
            Check that the center parameter can be represented in multiple ways
            TODO come back to this for the alternate coord system selecions
        '''
        
        imagename = makeImage()
        
        impv(imagename=imagename, outfile=testfile, center=[4.5, 5], pa='90deg', length=5,mode='length')
        impv(imagename=imagename, outfile=testfile2, center=["0:0:02", "0.0.0"], pa='90deg', length=5,mode='length')
        self.assertTrue(os.path.exists(testfile))
        self.assertTrue(os.path.exists(testfile2))
        
    def test_length(self):
        '''
            test_length
            -------------
            
            Check that the length parameter can be specified by numeric value or valid quantity
        '''
        
        impv(imagename=datapath, outfile=testfile, center=[45,50], pa='45deg', length=5, mode='length')
        self.assertTrue(os.path.exists(testfile))
        
        impv(imagename=datapath, outfile=testfile2, center=[45,50], pa='45deg', length='5arcmin', mode='length')
        self.assertTrue(os.path.exists(testfile2))
        
        
    def test_pa(self):
        '''
            test_pa
            ---------
            
            Check that the pa parameter works when provided with a valid quanitity
            
            NOTE: Does the default case need to be execured again, even though it is run in other parts of the tests?
            This may just be adding time without any real benifit to the test
        '''
        
        impv(imagename=datapath, outfile=testfile, center=[45,50], pa='45deg', length=5, mode='length')
        self.assertTrue(os.path.exists(testfile))
        
        
        
    def test_width(self):
        '''
            test_width
            ------------
            
            Check that the width parameer works when provided with a valid string quantity or a quantity record
        '''
        
        impv(imagename=datapath, outfile=testfile, center=[45,50], pa='45deg', length=5, mode='length', width='4arcsec')
        self.assertTrue(os.path.exists(testfile))
        
        impv(imagename=datapath, outfile=testfile2, center=[45,50], pa='45deg', length=5, mode='length', width=qa.quantity('4arcsec'))
        self.assertTrue(os.path.exists(testfile2))
        
        
    def test_unit(self):
        '''
            test_unit
            -----------
            
            Check that this parameter gives the unit for the offset axis in the resulting image
        '''
        
        impv(imagename=datapath, outfile=testfile, center=[45,50], pa='45deg', length=5, mode='length', unit='deg')
        self.assertTrue(os.path.exists(testfile))
        
        
    def test_tableRecord(self):
        '''
            test_tableRecord
            ------------------
            
            Check that the result is written to the output image as a table record, meaning it can be retrieved with the table tool.
        '''
        
        impv(imagename=datapath, outfile=testfile, center=[45,50], pa='45deg', length=5, mode='length')
        mytb.open(testfile)
        mytb.close()
    
    
def suite():
    return[impv_test]

if __name__ == '__main__':
    unittest.main()

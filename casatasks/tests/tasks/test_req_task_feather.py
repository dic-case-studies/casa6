##########################################################################
# test_req_task_feather.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_feather/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import feather, casalog
    tb = casatools.table()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import os
import unittest
import shutil
import numpy as np

### DATA ###

if CASA6:
#    interpath = casatools.ctsys.resolve('image/orion_tfeather.im/')
#    sdpath = casatools.ctsys.resolve('image/orion_tsdmem.image/')
    datapath = casatools.ctsys.resolve('unittest/feather/')

else:
    datapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/feather/'
#    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
#        interpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/orion_tfeather.im/'
#        sdpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/orion_tsdmem.image/'
        
#    else:
#   sdpath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/orion_tsdmem.image/'


#Input files
intimg = 'orion_tfeather.im'
sdimg = 'orion_tsdmem.image'

# Output files
output = 'feathered.im'
output2 = 'feathered2.im'

logpath = casalog.logfile()
logname = 'testlog.log'

def get_map(infile):

    tb.open(infile)
    res = tb.getcol('map')
    tb.close()
    
    return res

class feather_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def setUp(self):
        if not os.path.exists(intimg):
            os.symlink(os.path.join(datapath, intimg), intimg)
        if not os.path.exists(sdimg):
            os.symlink(os.path.join(datapath, sdimg), sdimg)
        
        if not CASA6:
            default(feather)
            
    def tearDown(self):
        if os.path.exists(output):
            shutil.rmtree(output)
            
        if os.path.exists(output2):
            shutil.rmtree(output2)
            
        if os.path.exists(logname):
            os.remove(logname)
            
        casalog.setlogfile(logpath)
    
    @classmethod
    def tearDownClass(cls):
        os.unlink(intimg)
        os.unlink(sdimg)
    
    def test_combine(self):
        '''
            test_combine
            --------------
            
            Check that interferometric and Single dish images can be combined
        '''
        
        feather(imagename=output, highres=intimg, lowres=sdimg)
        self.assertTrue(os.path.exists(output))
        
    def test_imagename(self):
        '''
            test_imagename
            ----------------
            
            Check that the imagename parameter gives the name of the output image file
        '''
        
        feather(imagename=output, highres=intimg, lowres=sdimg)
        feather(imagename=output2, highres=intimg, lowres=sdimg)
        
        self.assertTrue(os.path.exists(output))
        self.assertTrue(os.path.exists(output2))
        
    def test_highres(self):
        '''
            test_highres
            --------------
            
            Check that the interferometric image is provided with this parameter
            This parameter is nessisary to run the task
        '''
        
        if CASA6:
            with self.assertRaises(AssertionError):
                feather(imagename=output, lowres=sdimg)
        else:
            casalog.setlogfile(logname)
            feather(imagename=output, lowres=sdimg)
            self.assertTrue(('SEVERE' in open(logname).read()))
            
        
    def test_lowres(self):
        '''
            test_lowres
            -------------
            
            Check that the single dish image is provided with this parameter
            This parameter is nessisary to run the task
        '''
        
        if CASA6:
            with self.assertRaises(AssertionError):
                feather(imagename=output, highres=intimg)
        else:
            casalog.setlogfile(logname)
            feather(imagename=output, highres=intimg)
            self.assertTrue('SEVERE' in open(logname).read())
        
        
    def test_sdfactor(self):
        '''
            test_sdfactor
            ---------------
            
            Check that differing sdfactors results in differing image files
        '''
        
        feather(imagename=output, highres=intimg, lowres=sdimg)
        feather(imagename=output2, highres=intimg, lowres=sdimg, sdfactor=0.5)
        
        res1 = get_map(output)
        res2 = get_map(output2)
        
        self.assertFalse(np.all(np.isclose(res1, res2)))
        
    def test_effdishdiam(self):
        '''
            test_effdishdiam
            ------------------
            
            Check that chaging the effective dish diameter results in differing image files
        '''
        
        feather(imagename=output, highres=intimg, lowres=sdimg)
        feather(imagename=output2, highres=intimg, lowres=sdimg, effdishdiam=1)
        
        res1 = get_map(output)
        res2 = get_map(output2)
        
        self.assertFalse(np.all(np.isclose(res1, res2)))
        
        if CASA6:
            with self.assertRaises(RuntimeError):
                feather(imagename=output2, highres=intimg, lowres=sdimg, effdishdiam=1000)
        else:
            casalog.setlogfile(logname)
            feather(imagename=output2, highres=intimg, lowres=sdimg, effdishdiam=1000)
            self.assertTrue('SEVERE' in open(logname).read())
        
    def test_lowpassfiltersd(self):
        '''
            test_lowpassfiltersd
            ----------------------
            
            Check that lowpassfiltersd = True results in a different image than the default
        '''
        
        feather(imagename=output, highres=intimg, lowres=sdimg)
        feather(imagename=output2, highres=intimg, lowres=sdimg, lowpassfiltersd=True)
        
        res1 = get_map(output)
        res2 = get_map(output2)
        
        self.assertFalse(np.all(np.isclose(res1, res2)))
    
    
def suite():
    return[feather_test]

if __name__ == '__main__':
    unittest.main()

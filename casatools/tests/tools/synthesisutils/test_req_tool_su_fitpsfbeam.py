##########################################################################
# test_req_tool_su_fitpsfbeam.py
#
# Copyright (C) 2021
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
# [https://open-jira.nrao.edu/browse/CAS-13590]
#
# Based on the requirements listed in CASADocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatools.synthesisutils.html
#
# Test case: requirement
# test single term psf fitting 
# test a wrong nterms (>1) for single term psf fitting
# test a wrong psfcutoff (<0, 1) 
# test mulit-term psf fitting (nterms=2)
# 
##########################################################################
 
 
####    Imports     ####
import os
import sys
import unittest
import copy
import shutil
 
# Example of importing helper functions
from casatestutils import testhelper as th
from casatestutils.imagerhelpers import TestHelpers
 
is_CASA6 = False
 
try:
    from casatools import ctsys, synthesisutils, image
    #from casatasks import ...
    is_CASA6 = True
    su = synthesisutils() 
    _ia = image()
    # Location of input data
    #datapath = ctsys.resolve('unittest/synthesisutils/')
    datapath = '/export/home/murasame/casasrc/cas13590dev/testdata/'
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    su = casac.synthesisutils() 
    _ia = imtool()
    # Location of input data
    datapath = os.path.join(os.environ['CASAPATH'].split()[0], 'casatestdata/unittest/synthesisutils/')
 
####    Tests     ####
class su_fitpsfbeam_test(unittest.TestCase):
    ### su.fitpsBeam has three parameters, imagename, nterms, and psfcutoff.
    ### nterms and psfcutoff have the default values but other values need to be tested.
    ### Approriate multiterm psfs need to be present for nterms>1
    
    ### Set Up
    def setUp(self):
        # input images
        self.inputims=[]
        psfim1 = 'su_fitpsfbeam_test_mfs.psf'
        self.inputims.append(psfim1)
        psfim1fullpath = os.path.join(datapath, psfim1)
        shutil.copytree(psfim1fullpath,psfim1)
        # base name for multiterm 
        psfim2base = 'su_fitpsfbeam_test_mtmfs.psf.tt'
        for ext in ['0','1','2']:
            psfname = psfim2base+ext
            psf2ttfullpath = os.path.join(datapath, psfname) 
            shutil.copytree(psf2ttfullpath,psfname)
            self.inputims.append(psfname)
            
    ### Teardown
    def tearDown(self):
        for img in self.inputims:
            shutil.rmtree(img)
 
    ### Test Cases
    def test_mfs(self):
        '''Test that fitting of mfs psf works '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]
        _ia.open(psfim)
        origbeam = _ia.restoringbeam()
        # Set different values for the beam
        resetbeam = copy.deepcopy(origbeam)
        resetbeam['major']['value']=1.0
        resetbeam['minor']['value']=1.0
        resetbeam['positionangle']['value']=0.0
        _ia.setrestoringbeam(major=resetbeam['major'], minor=resetbeam['minor'], pa=resetbeam['positionangle'])
        _ia.done()
 
        # expected values
        bmref= {'major': {'unit': 'arcsec', 'value': 52.93011474609375}, 'minor': {'unit': 'arcsec', 'value': 49.147438049316406}, 'positionangle': {'unit': 'deg', 'value': -87.39420318603516}}
 
        ret = su.fitPsfBeam(imagename=psfim[:-4])   
        _ia.open(psfim)
        newbeam = _ia.restoringbeam()
        _ia.done()
        print("psfim=",psfim)
        print("origbeam=",origbeam)
        print("newbeam=",newbeam)     

        self.assertTrue(ret)
        # test against expected values
        # - original beam in psf image is different from fitted values
        self.assertDictContainsSubset(newbeam, bmref)

    def test_mfs_wrong_nterms(self):
        '''Test that it catches if nterms is inconsistent with input psf (nterms=2, for a single term  psf) '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]

        self.assertRaises(RuntimeError, su.fitPsfBeam, imagename=psfim[:-4],nterms=2)   

    def test_mfs_wrong_psfcutoff(self):
        '''Test that psfcutoff is given outside the allowed range (1.0) '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]

        ret1 = su.fitPsfBeam(imagename=psfim[:-4],psfcutoff=1.0)   
        ret2 = su.fitPsfBeam(imagename=psfim[:-4],psfcutoff=-1.0)   
        self.assertFalse(ret1)
        self.assertFalse(ret2)

    def test_mfs_largerpsfcutoff(self):
        '''Test that psfcutoff with a valid (larger) number  works '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]
        _ia.open(psfim)
        origbeam = _ia.restoringbeam()
        # Set different values for the beam
        resetbeam = copy.deepcopy(origbeam)
        resetbeam['major']['value']=1.0
        resetbeam['minor']['value']=1.0
        resetbeam['positionangle']['value']=0.0
        _ia.setrestoringbeam(major=resetbeam['major'], minor=resetbeam['minor'], pa=resetbeam['positionangle'])
        _ia.done()
 
        # expected values
        bmref= {'major': {'unit': 'arcsec', 'value': 49.89558029174805}, 'minor': {'unit': 'arcsec', 'value': 46.80238342285156}, 'positionangle': {'unit': 'deg', 'value': -88.28898620605469}}
 
        ret = su.fitPsfBeam(imagename=psfim[:-4],psfcutoff=0.5)   
        _ia.open(psfim)
        newbeam = _ia.restoringbeam()
        _ia.done()
        print("psfim=",psfim)
        print("origbeam=",origbeam)
        print("newbeam=",newbeam)     

        self.assertTrue(ret)
        # test against expected values
        # - original beam in psf image is different from fitted values
        self.assertDictContainsSubset(newbeam, bmref)

    def test_mtmfs_nterms2(self):
        '''Test that fitting of multiterm  psf works '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfimtt0 = self.inputims[1]
        _ia.open(psfimtt0)
        origbeam = _ia.restoringbeam()
        # Set different values for the beam
        resetbeam = copy.deepcopy(origbeam)
        resetbeam['major']['value']=1.0
        resetbeam['minor']['value']=1.0
        resetbeam['positionangle']['value']=0.0
        _ia.setrestoringbeam(major=resetbeam['major'], minor=resetbeam['minor'], pa=resetbeam['positionangle'])
        _ia.done()
 
        # expected values
        bmref= {'major': {'unit': 'arcsec', 'value': 52.93011474609375}, 'minor': {'unit': 'arcsec', 'value': 49.147438049316406}, 'positionangle': {'unit': 'deg', 'value': -87.39420318603516}}
 
        #print("psfimtt0==",psfimtt0)
        ret = su.fitPsfBeam(imagename=psfimtt0[:-8],nterms=2)   
        _ia.open(psfimtt0)
        newbeam = _ia.restoringbeam()
        _ia.done()
        #print("psfimitt0=",psfimtt0)
        #print("origbeam=",origbeam)
        #print("newbeam=",newbeam)     

        self.assertTrue(ret)
        # test against expected values
        # - original beam in psf image is different from fitted values
        self.assertDictContainsSubset(newbeam, bmref)
 
####    Suite: Required for CASA5     ####
def suite():
    return[su_fitpsfbeam_test]
  
####    Imports     ####
if __name__ == '__main__':
    unittest.main()
 
 

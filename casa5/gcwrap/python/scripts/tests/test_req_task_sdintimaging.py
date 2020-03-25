##########################################################################
# test_req_task_sdintimaging.py
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
# Parent ticket for sdintiamging task implementation
# https://open-jira.nrao.edu/browse/CAS-12613 
# Ticket for the test generation
# https://open-jira.nrao.edu/browse/CAS-12898
#
# Based on the requirements listed in plone found here (TBD):
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_sdintimaging/about
#
#
##########################################################################
# List of tests
#Single Pointing Simulation:
#
#1. Wideband Multi-Term Imaging: single dish - interferometer combination (SD+INT)
#specmode='mfs', deconvolver='mtmfs', gridder='standard', usedata='sdint'
#testname: test_singlepointing_mfs_sdint 
#
#2. Wideband Multi-Term Imaging: inteferometer only (INT-only) 
#specmode='mfs', deconvolver='mtmfs', gridder='standard', usedata='int'
#test name: test_singlepointing_mfs_intonly
#
#3. Wideband Multi-Term Imaging: single dish only (SD-only)
#specmode='mfs', deconvolver='mtmfs', gridder='standard', usedata='sd'
#testname: test_singlepointing_mfs_sdonly
#
#4. Cube imaging: single dish - interferometer combination (SD+INT)
#specmode='cube',  deconvolver='multiscale', gridder='standard', usedasta='sdint'
#testname: test_singlepointing_cube_sdint
#
#5. Cube imaging: interferometer only (INT-only)
#specmode='cube', deconvolver='multiscale', gridder='standard', usedata='int' 
#testname: test_singlepointing_cube_intonly
#
#6. Cube imaging: single dish only (SD-only)
#specmode='cube', deconvolver='multiscale', gridder='standard', usedata='sd'
#testname: test_singlepointing_cube_sdonly
#
#Mosaic Simulation:
#7. Wideband Multi-Term Mosaic Imaging: single dish - interferometer combination (SD+INT)
#specmode='mfs',  deconvolver='mtmfs', gridder='mosaic', usedata='sdint'
#testname: test_mosaic_mfs_sdint
#
#8. Wideband Multi-Term Mosaic Imaging: interforometer only (INT-only)
#specmode='mfs',  deconvolver='mtmfs', gridder='mosaic', usedata='int'
#testname: test_mosaic_mfs_intonly
#
#9. Wideband Multi-Term Mosaic Imaging: single dish only (SD-only)
#specmode='mfs',  deconvolver='mtmfs', gridder='mosaic', usedata='sd'
#testname: test_mosaic_mfs_sdonly
#
#10. Cube Imaging: single dish  - interferometer combination (SD+INT)
#specmode='cube',  deconvolver='multiscale', gridder='mosaic', usedata='sdint'
#testname: test_mosaic_cube_sdint
#
#11. Cube Imaging: interferometer only (INT-only)
#specmode='cube',  deconvolver='multiscale', gridder='mosaic', usedata='int'
#testname: test_mosaic_cube_int
#
#12. Cube Imaging: single dish (SD-only)
# specmode='cube',  deconvolver='multiscale', gridder='mosaic', usedata='sd'
#testname: test_mosaic_cube_sd
###########################################################################

####    Imports     ####
import os
import sys
import unittest
#import casaTestHelper as th
import shutil
import numpy as np
from casatestutils.imagerhelpers import TestHelpers 

CASA6 = False
try: 
    from casatools import ctsys
    from casatasks import casalog, sdintimaging
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
    sys.path.append(os.path.abspath(os.path.basename(__file__)))
    #from testhelper_imager import TestHelpers

     #refdatapath = '/export/home/murasame2/casadev/imagerRefact/sdint/orig_scripts/WidebandSDINT/Data'
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from parallel.parallel_task_helper import ParallelTaskHelper
    from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
    #from imagerhelpers.testhelper_imager import TestHelpers

th = TestHelpers()

if CASA6:
    #refdatapath = ctsys.resolve('regression/unittest/sdintimaging')
    visdatapath = ctsys.resolve('visibilities/evla')
    imdatapath = ctsys.resolve('image') 
    maskdatapath = ctsys.resolve('text') 
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
    #refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdintimaging/'
        visdatapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/evla/'
        imdatapath = os.environ.get('CASAPATH').split()[0] +'/data/casa-data-req/image/' 
        maskdatapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/text/'
    else:
        visdatapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/evla/'
        imdatapath = os.environ.get('CASAPATH').split()[0] +'/image/' 
        maskdatapath = os.environ.get('CASAPATH').split()[0] + '/text/'
    
    #if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
    #    refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/evla/'
    #else:
    #    refdatapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/evla/'

    #For local testing
    #visdatapath = '/home/vega/rurvashi/TestCASA/VerificationTests/WBSDINT/Data/'
    #imdatapath = '/home/vega/rurvashi/TestCASA/VerificationTests/WBSDINT/Data/'
    #maskdatapath = '/home/vega/rurvashi/TestCASA/VerificationTests/WBSDINT/Data/'


class testref_base(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self.epsilon = 0.05
        self.msfile = ""
        self.img = "tst"
        #self.cfcache = 'cfcach'
        # To use subdir in the output image names in some tests (CAS-10937)
        #self.img_subdir = 'refimager_tst_subdir'
        self.parallel = False
        self.nnode = 0
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True
            self.PH = PyParallelImagerHelper()
            self.nnode = len(self.PH.getNodeList())
        #self.th = TestHelpers()
 

    def tearDown(self):
        # Default: delete all (input and output data)
        self.delData()
        # leave for input and output (e.g. for debugging)
        #self.delData(delinput=False, deloutput=False)

    @classmethod
    def tearDownClass(cls):
        pass
    

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self,inputdata={}):
        # clean-up 
        #os.system('rm -rf ' + self.img_subdir)
        os.system('rm -rf ' + self.img+'*')

        if inputdata!={}:
            if 'msname' in inputdata:
                self.msfile=inputdata['msname']
                if (os.path.exists(self.msfile)):
                    os.system('rm -rf ' + self.msfile)
                shutil.copytree(os.path.join(visdatapath,self.msfile), self.msfile)
            if 'sdimage' in inputdata:
                self.sdimage=inputdata['sdimage']
                if (os.path.exists(self.sdimage)):
                    os.system('rm -rf ' + self.sdimage)
                shutil.copytree(os.path.join(imdatapath,self.sdimage), self.sdimage)
            if 'sdpsf' in inputdata:
                self.sdpsf=inputdata['sdpsf']
                if (os.path.exists(self.sdpsf)):
                    os.system('rm -rf ' + self.sdpsf)
                shutil.copytree(os.path.join(imdatapath,self.sdpsf), self.sdpsf)
            if 'mask' in inputdata:
                self.mask=inputdata['mask']
                if (os.path.exists(self.mask)):
                    os.system('rm -rf ' + self.mask)
                origmask = os.path.join(maskdatapath,self.mask)
                if os.path.exists(origmask):
                    if os.path.isfile(origmask):
                        shutil.copyfile(origmask, self.mask)
                    else:
                        shutil.copytree(origmask, self.mask)

    def delData(self,delinput=True, deloutput=True):
        ''' delete input and output data or some '''
        #if msname != "":
        #    self.msfile=msname
        #if (os.path.exists(self.msfile)):
        #    os.system('rm -rf ' + self.msfile)
        #os.system('rm -rf ' + self.img_subdi)
        if delinput:
            if self.msfile!='':
                os.system('rm -rf ' + self.msfile)
            if self.sdimage!='':
                os.system('rm -rf ' + self.sdimage)
            if self.sdpsf!='':
                os.system('rm -rf ' + self.sdpsf)
            if self.mask!='':
                os.system('rm -rf ' + self.mask)
        if deloutput:
            os.system('rm -rf ' + self.img+'*')

### functional tests for sdintimaging start here ####

class test_singlepointing(testref_base):

    def setUp(self):
        super(test_singlepointing, self).setUp()
        # casa6
        # super().setUp() 
        # Set common parameters here
        self.imsize=800
        self.cell = '9.0arcsec'
        self.phasecenter = 'J2000 19:59:28.500 +40.44.01.50'
        self.nchan=3
        self.reffreq='1.5GHz'
        self.scales=[0,12,20,40,60,80,100]
        self.pblimit=-0.1  ## Set to negative value as the SinglePointing simulation has no primary beams
        self.interpolation='nearest'
        ################ Niter = 1000 orginal runtest parameter
        ################  Required to check absolute numerical accuracy (i.e. does it converge correctly)
        #self.niter=1000
        #self.cycleniter= 200
        ################ Niter=100  for quicker execution 
        ################ Good enough for basic checks and to catch numerical changes.
        self.niter=100
        self.cycleniter=50

    # Test 1
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_mfs_sdint(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='mfs', usedata='sdint')
        """ [singlePointing] Test_singlepointing_mfs_sdint """
        ######################################################################################
        # Test single field imaging for sd+int combination - mfs 
        # main parameters to be tested: specmode='mfs', usedata='sdint', gridder='standard'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_standard.ms',
                   'sdimage':'papersky_standard.sdimage',
                   'sdpsf':'papersky_standard.sdpsf',
                   'mask':'papersky_standard.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='mtmfs'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)

        imname=self.img+'.sp_mfs_sdint'

        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0,parallel=self.parallel)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 0.991, [400,400,0,0]),
                                   (outimg+'.image.tt0', 1.1820, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.28, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.956, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.179, [300,400,0,0]) ])      # extended emission with alpha=0
        

    #Test 2
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_mfs_intonly(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='mfs', usedata='int')
        """ [singlePointing] Test_singlepointing_mfs_intonly """
        ######################################################################################
        # Test single field imaging for intonly - mfs (should work without sdimage and sdpsf being set))
        # main parameters to be tested: specmode='mfs', usedata='int', gridder='standard'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_standard.ms',
                     'sdimage':'papersky_standard.sdimage',
                     'sdpsf':'papersky_standard.sdpsf',
                     'mask':'papersky_standard.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        #
        # Other secondary non-default parameters: 
        deconvolver='mtmfs'
        cycleniter=20
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.sp_mfs_intonly'
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0,parallel=self.parallel)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),
                                   (outimg+'.image.tt0', 2.25, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.21, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.819, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', -3.6, [300,400,0,0]) ])      # extended emission with alpha=0 ( will be steep for intonly)
        ## Since this is int_only, the values will be wrong.

    # Test 3
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_mfs_sdonly(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='mfs', usedata='sd')
        """ [singlePointing] Test_singlepointing_mfs_sdonly """
        ######################################################################################
        # Test single field imaging for sdonly - mfs 
        # main parameters to be tested: specmode='mfs', usedata='sd', gridder='standard'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_standard.ms',
                   'sdimage':'papersky_standard.sdimage',
                   'sdpsf':'papersky_standard.sdpsf',
                   'mask':'papersky_standard.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='mtmfs'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.sp_mfs_sdonly'
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0,parallel=self.parallel)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),
                                   (outimg+'.image.tt0', 7.91, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 15.3, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.13, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.018, [300,400,0,0]) ])      # extended emission with alpha=0

    #Test4
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_cube_sdint(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='cube', usedata='sdint')
        """ [singlePointing] Test_singlepointing_cube_sdint """
        ######################################################################################
        # Test single field imaging for sdint - cube 
        # main parameters to be tested: specmode='cube', usedata='sdint', gridder='standard'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_standard.ms',
                   'sdimage':'papersky_standard.sdimage',
                   'sdpsf':'papersky_standard.sdpsf',
                   'mask':'papersky_standard.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='multiscale'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.sp_cube_sdint'
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0,parallel=self.parallel)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 0.99, [400,400,0,0]),
                                   (outimg+'.psf', 0.99, [400,400,0,1]),
                                   (outimg+'.image', 1.66, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.459, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.07, [350,433,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 0.16, [300,400,0,1]) ])      # extended emission with alpha=0
        ## Check multiple channels. point source flux is same, extended emission will be different because of resolution change.

    #Test5
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_cube_intonly(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='cube', usedata='int')
        """ [singlePointing] Test_singlepointing_cube_intonly """
        ######################################################################################
        # Test single field imaging for int - cube 
        # main parameters to be tested: specmode='cube', usedata='int', gridder='standard'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_standard.ms',
                   'sdimage':'papersky_standard.sdimage',
                   'sdpsf':'papersky_standard.sdpsf',
                   'mask':'papersky_standard.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='multiscale'
        cycleniter=20
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.sp_cube_intonly'
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0,parallel=self.parallel)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [400,400,0,0]),
                                   (outimg+'.psf', 1.0, [400,400,0,1]),
                                   (outimg+'.image', 2.13, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.59, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.424, [350,433,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', -0.03, [300,400,0,1]) ])      # extended emission with alpha=0

    #Test6
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_cube_sdonly(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='cube', usedata='sd')
        """ [singlePointing] Test_singlepointing_cube_sdonly """
        ######################################################################################
        # Test single field imaging for sd - cube 
        # main parameters to be tested: specmode='cube', usedata='sd', gridder='standard'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_standard.ms',
                   'sdimage':'papersky_standard.sdimage',
                   'sdpsf':'papersky_standard.sdpsf',
                   'mask':'papersky_standard.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='multiscale'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.sp_cube_sdonly'

        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0,parallel=self.parallel)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [400,400,0,0]),
                                   (outimg+'.psf', 1.0, [400,400,0,1]),
                                   (outimg+'.image', 32.43, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 64.11, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 10.76, [350,433,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 22.64, [300,400,0,1]) ])      # extended emission with alpha=0


class test_mosaic(testref_base):

    def setUp(self):
        # common parameters
        super(test_mosaic, self).setUp()
        self.imsize=1500
        self.cell = '9.0arcsec'
        self.phasecenter = 'J2000 19:59:28.500 +40.44.01.50'
        self.nchan=3
        self.reffreq='1.5GHz'
        self.scales=[0,12,20,40,60,80,100]
        self.pblimit=0.1  # Set to a positive value since this is mosaics with PBs.  NOTE : pblimit=-0.1 gives NaNs (need to fix)
        self.interpolation='nearest'
        ################ Niter = 1000 orginal runtest parameter
        ################  Required to check absolute numerical accuracy (i.e. does it converge correctly)
        #self.niter=1000
        #self.cycleniter= 200
        ################ Niter=100  for quicker execution 
        ################ Good enough for basic checks and to catch numerical changes.
        self.niter=100
        self.cycleniter= 50

     #Test7
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mosaic_mfs_sdint(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='mfs', usedata='sdint')
        """ [Mosaic] Test_mosaic_mfs_sdint """
        ######################################################################################
        # Test mosaic imaging for sdint - mfs 
        # main parameters to be tested: specmode='mfs', usedata='sdint', gridder='mosaic'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_mosaic.ms',
                   'sdimage':'papersky_mosaic.sdimage',
                   'sdpsf':'papersky_mosaic.sdpsf',
                   'mask':'papersky_mosaic.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='mtmfs'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.mos_mfs_sdint'
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2,parallel=self.parallel)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 0.9905, [750,750,0,0]),
                                   (outimg+'.image.tt0', 1.106, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.264, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.949, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.231, [650,720,0,0]) ])      # extended emission with alpha=0

    #Test8
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mosaic_mfs_intonly(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='mfs', usedata='int')
        """ [Mosaic] Test_mosaic_mfs_intonly """
        ######################################################################################
        # Test mosaic imaging for int - mfs 
        # main parameters to be tested: specmode='mfs', usedata='int', gridder='mosaic'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_mosaic.ms',
                   'sdimage':'papersky_mosaic.sdimage',
                   'sdpsf':'papersky_mosaic.sdpsf',
                   'mask':'papersky_mosaic.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='mtmfs'
        cycleniter=20
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.mos_mfs_intonly'
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2,parallel=self.parallel)
        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [750,750,0,0]),
                                   (outimg+'.image.tt0', 2.021, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.63, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.962, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', -1.02, [650,720,0,0]) ])      # extended emission with alpha=0 (steep with intonly)

    #Test9
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mosaic_mfs_sdonly(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='mfs', usedata='sd')
        """ [Mosaic] Test_mosaic_mfs_sdonly """
        ######################################################################################
        # Test mosaic imaging for sd - mfs 
        # main parameters to be tested: specmode='mfs', usedata='sd', gridder='mosaic'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_mosaic.ms',
                   'sdimage':'papersky_mosaic.sdimage',
                   'sdpsf':'papersky_mosaic.sdpsf',
                   'mask':'papersky_mosaic.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='mtmfs'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.mos_mfs_sdonly'
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2,parallel=self.parallel)
        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [750,750,0,0]),
                                   (outimg+'.image.tt0', 7.756, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 15.68, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.12, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.013, [650,720,0,0]) ])      # extended emission with alpha=0

    #Test10
#    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mosaic_cube_sdint(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='cube', usedata='sdint')
        """ [Mosaic] Test_mosaic_cube_sdint """
        ######################################################################################
        # Test mosaic imaging for sdint - cube
        # main parameters to be tested: specmode='cube', usedata='sdint', gridder='mosaic'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_mosaic.ms',
                   'sdimage':'papersky_mosaic.sdimage',
                   'sdpsf':'papersky_mosaic.sdpsf',
                   'mask':'papersky_mosaic.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='multiscale'
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.mos_cube_sdint'
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2,parallel=self.parallel)
        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 0.99, [750,750,0,0]),
                                   (outimg+'.psf', 0.99, [750,750,0,1]),
                                   (outimg+'.image', 1.57, [700,783,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.55, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.027, [700,783,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 0.199, [650,720,0,1]) ])      # extended emission with alpha=0

    #Test11
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mosaic_cube_intonly(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='cube', usedata='int')
        """ [Mosaic] Test_mosaic_cube_intonly """
        ######################################################################################
        # Test mosaic imaging for intonly - cube
        # main parameters to be tested: specmode='cube', usedata='int', gridder='mosaic'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_mosaic.ms',
                   'sdimage':'papersky_mosaic.sdimage',
                   'sdpsf':'papersky_mosaic.sdpsf',
                   'mask':'papersky_mosaic.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='multiscale'
        cycleniter=20
        # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.mos_cube_intonly'
        ret = sdintimaging(usedata='int', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2,parallel=self.parallel)
        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [750,750,0,0]),
                                   (outimg+'.psf', 1.0, [750,750,0,1]),
                                   (outimg+'.image', 2.044, [700,783,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 1.023, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.175, [700,783,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 0.05, [650,720,0,1]) ])      # extended emission with alpha=0

    #Test12
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mosaic_cube_sdonly(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='cube', usedata='sd')
        """ [Mosaic] Test_mosaic_cube_sdonly """
        ######################################################################################
        # Test mosaic imaging for sdonly - cube
        # main parameters to be tested: specmode='cube', usedata='sd', gridder='mosaic'
        # with the default weighting (='natural')
        ######################################################################################
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'msname':'papersky_mosaic.ms',
                   'sdimage':'papersky_mosaic.sdimage',
                   'sdpsf':'papersky_mosaic.sdpsf',
                   'mask':'papersky_mosaic.true.im.masklist'}
        # data specific parameters 
        # imsize, cell, phasecenter, reffreq, nchan, scales 
        # set to the default values for sdgain (1.0) and dishdia (100.0)
        #
        # Other secondary non-default parameters: 
        deconvolver='multiscale'
       # iterations may need to be shorten for the final version of test
        self.prepData(inputdata=inputdata)
        imname=self.img+'.mos_cube_sdonly'
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2,parallel=self.parallel)
        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [750,750,0,0]),
                                   (outimg+'.psf', 1.0, [750,750,0,1]),
                                   (outimg+'.image', 32.70, [700,783,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 66.239, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 10.92, [700,783,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 25.50, [650,720,0,1]) ])      # extended emission with alpha=0


def suite():
    return[test_singlepointing,test_mosaic]

if __name__ == '__main__':
    unittest.main() 

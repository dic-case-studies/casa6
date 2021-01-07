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
#
#Special Cases
#13. Single Pointing Test with SD+INT data, with different channels flagged in SD and INT
#testname: test_singlepointing_mfs_sdint_flagged
#
#14. Single Pointing Test with SD+INT data, with different channels flagged in SD and INT
#testname: test_singlepointing_cube_sdint_flagged
#
#15. Single Pointing Test with SD+INT data, with sdpsf="" and internal auto-calculation.
#testname: test_singlepointing_mfs_sdint_autopsf
#
#16. Single pointing test with INT-only data from refim_point.ms : Compare with tclean cube
#testname: test_intonly_cube_compare_with_tclean
#
#17. Single pointing test with INT-only data from refim_point.ms : Compare with tclean mtmfs
#testname: test_intonly_mfs_compare_with_tclean
#
###########################################################################

####    Imports     ####
import os
import sys
import unittest
#import casaTestHelper as th
import shutil
import numpy as np
from casatestutils.imagerhelpers import TestHelpers 
import inspect

CASA6 = False
try: 
    from casatools import ctsys
    from casatasks import casalog, sdintimaging, flagdata, tclean
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
    refdatapath = ctsys.resolve('regression/unittest/clean/refimager/')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/casa-data-req'):
        visdatapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/evla/'
        imdatapath = os.environ.get('CASAPATH').split()[0] +'/casa-data-req/image/' 
        maskdatapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/text/'
    else:
        visdatapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/evla/'
        imdatapath = os.environ.get('CASAPATH').split()[0] +'/image/' 
        maskdatapath = os.environ.get('CASAPATH').split()[0] + '/text/'

    refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/clean/refimager/'

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

#        self.parallel = False
#        self.nnode = 0
#        if ParallelTaskHelper.isMPIEnabled():
#            self.parallel = True
#            self.PH = PyParallelImagerHelper()
#            self.nnode = len(self.PH.getNodeList())
        #self.th = TestHelpers()
 

    def tearDown(self):
        # Default: delete all (input and output data)
        #self.delData()
        # leave for input and output (e.g. for debugging)
        self.delData(delinput=False, deloutput=False)

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
            if 'refmsname' in inputdata:
                self.refmsfile=inputdata['refmsname']
                if (os.path.exists(self.refmsfile)):
                    os.system('rm -rf ' + self.refmsfile)
                shutil.copytree(os.path.join(refdatapath,self.refmsfile), self.refmsfile)
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
            if hasattr(self,'msfile') and self.msfile!='':
                os.system('rm -rf ' + self.msfile)
            if hasattr(self,'refmsfile') and self.refmsfile!='':
                os.system('rm -rf ' + self.refmsfile)
            if hasattr(self,'sdimage') and self.sdimage!='':
                os.system('rm -rf ' + self.sdimage)
            if hasattr(self,'sdpsf') and self.sdpsf!='':
                os.system('rm -rf ' + self.sdpsf)
            if hasattr(self,'mask') and self.mask!='':
                os.system('rm -rf ' + self.mask)
        if deloutput:
            os.system('rm -rf ' + self.img+'*')


    def checkfinal(self,pstr=""):
          pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  runUnitTest.main(['test_req_task_sdintimaging["+ inspect.stack()[1][3] +"]'])"
          casalog.post(pstr,'INFO')
          if( pstr.count("( Fail") > 0 ):
              print(pstr)
              self.fail("\n"+pstr)

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
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        if self.niter==100:
            incycleniter=20 # overwrite the initial setup for niter=100 to make the test pass for 6.1 (need furhter investigation)

        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=0, cycleniter=incycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 0.991, [400,400,0,0]),
                                   (outimg+'.image.tt0', 1.187, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.262, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.954, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.174, [300,400,0,0]) ])      # extended emission with alpha=0
        
        self.checkfinal(pstr=report)

    #Test 2
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),
                                   (outimg+'.image.tt0', 1.09, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.1, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.996, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', -2.35, [300,400,0,0]) ])      # extended emission with alpha=0 ( will be steep for intonly)
        ## Since this is int_only, the values will be wrong.
        self.checkfinal(pstr=report)


    # Test 3
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),
                                   (outimg+'.image.tt0', 7.91, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 15.3, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.137, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.018, [300,400,0,0]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


    #Test4
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 0.99, [400,400,0,0]),
                                   (outimg+'.psf', 0.99, [400,400,0,1]),
                                   (outimg+'.image', 1.66, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.459, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.07, [350,433,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 0.168, [300,400,0,1]) ])      # extended emission with alpha=0
        ## Check multiple channels. point source flux is same, extended emission will be different because of resolution change.
        self.checkfinal(pstr=report)


    #Test5
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [400,400,0,0]),
                                   (outimg+'.psf', 1.0, [400,400,0,1]),
                                   (outimg+'.image', 1.48, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.347, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.013, [350,433,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', -0.036, [300,400,0,1]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


    #Test6
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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

        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [400,400,0,0]),
                                   (outimg+'.psf', 1.0, [400,400,0,1]),
                                   (outimg+'.image', 18.65, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 33.15, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 7.93, [350,433,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 15.33, [300,400,0,1]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)

    # Test 13
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_mfs_sdint_flagged(self):
        """ [singlePointing] Test_singlepointing_mfs_sdint_flagged """
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

        ## Flag some data
        flagdata(vis=self.msfile, spw='0:2')  ## Last channel of interferometer data

        imname=self.img+'.sp_mfs_sdint'

        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 0.990, [400,400,0,0]),
                                   (outimg+'.image.tt0', 1.144, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.371, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -1.29, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.11, [300,400,0,0]) ])      # extended emission with alpha=0
        
        self.checkfinal(pstr=report)

    #Test 14
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_cube_sdint_flagged(self):
        """ [singlePointing] Test_singlepointing_cube_sdint_flagged """
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

        ## Flag some data
        flagdata(vis=self.msfile, spw='0:2')  ## Last channel of interferometer data


        imname=self.img+'.sp_cube_sdint'
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)

        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 0.99, [400,400,0,0]),
                                   (outimg+'.psf', 0, [400,400,0,2]),
                                   (outimg+'.image', 1.66, [350,433,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.459, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 0, [350,433,0,2]),    # point source of 1 Jy
                                   (outimg+'.image', 0, [300,400,0,2]) ])      # extended emission with alpha=0
        ## Check multiple channels. point source flux is same, extended emission will be different because of resolution change.
        self.checkfinal(pstr=report)

    #Test 15 
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_singlepointing_mfs_sdint_autopsf(self):
        # Equivalent to onetest(runtype='SinglePointing', specmode='mfs', usedata='sdint')
        """ [singlePointing] Test_singlepointing_mfs_sdint_autopsf """
        ########################
        # Same as Test1, but with auto SD PSF calculation.
        ###########################
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

        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf="", vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.0)


        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 0.990, [400,400,0,0]),
                                   (outimg+'.image.tt0', 1.189, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.261, [300,400,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.939, [350,433,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.0736, [300,400,0,0]) ])      # extended emission with alpha=0
        
        self.checkfinal(pstr=report)


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
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2)

        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 0.9905, [750,750,0,0]),
                                   (outimg+'.image.tt0', 1.098, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.268, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.95, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.248, [650,720,0,0]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


    #Test8
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2)
        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [750,750,0,0]),
                                   (outimg+'.image.tt0', 1.05, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 0.14, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -1.016, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', -0.78, [650,720,0,0]) ])      # extended emission with alpha=0 (steep with intonly)
        self.checkfinal(pstr=report)


    #Test9
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2)
        outimg = imname+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', 
                                     outimg+'.residual.tt0', outimg+'.image.tt0', 
                                     outimg+'.image.tt1',outimg+'.alpha'], 
                           imgval=[(outimg+'.psf.tt0', 1.0, [750,750,0,0]),
                                   (outimg+'.image.tt0', 7.756, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.image.tt0', 15.68, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.alpha', -0.12, [700,783,0,0]),    # point source with alpha=-1
                                   (outimg+'.alpha', 0.013, [650,720,0,0]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


    #Test10
#    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2)
        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 0.99, [750,750,0,0]),
                                   (outimg+'.psf', 0.99, [750,750,0,1]),
                                   (outimg+'.image', 1.554, [700,783,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.485, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 1.014, [700,783,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 0.187, [650,720,0,1]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


    #Test11
    #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='int', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2)
        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [750,750,0,0]),
                                   (outimg+'.psf', 1.0, [750,750,0,1]),
                                   (outimg+'.image', 1.452, [700,783,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 0.41, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 0.917, [700,783,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 0.050149, [650,720,0,1]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


    #Test12
  #  @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=imname,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,pbmask=0.2)
        outimg = imname+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf', 
                                     outimg+'.residual', outimg+'.image'], 
                           imgval=[(outimg+'.psf', 1.0, [750,750,0,0]),
                                   (outimg+'.psf', 1.0, [750,750,0,1]),
                                   (outimg+'.image', 18.17, [700,783,0,0]),    # point source of 1 Jy
                                   (outimg+'.image', 33.16, [650,720,0,0]),        # extended emission with alpha=0
                                   (outimg+'.image', 7.932, [700,783,0,1]),    # point source of 1 Jy
                                   (outimg+'.image', 16.319, [650,720,0,1]) ])      # extended emission with alpha=0
        self.checkfinal(pstr=report)


######################################
##### Compare with tclean
######################################
class test_compare_sdint_tclean(testref_base):

    def setUp(self):
        # common parameters
        super(test_compare_sdint_tclean, self).setUp()

    #Test16
    #16. Single pointing test with INT-only data from refim_point.ms : Compare with tclean cube
    #testname: test_intonly_cube_compare_with_tclean
    #
    def test_intonly_cube_compare_with_tclean(self):
        """ [Compare] Test_intonly_cube_compare_with_tclean """
        # inputdata: set of the data to be copied from the data repos or else where during setup. 
        inputdata={'refmsname':'refim_point.ms'}
        self.prepData(inputdata=inputdata)
        imname1=self.img+'.sdint'
        ret1 = sdintimaging(vis=self.refmsfile,imagename=imname1, usedata='int', imsize=200, cell='10.0arcsec', nchan=3, spw='0:0~2', pblimit=0.1, interpolation='nearest',specmode='cube',niter=15, cycleniter=5, gridder='mosaic',mosweight=False, deconvolver='multiscale',scales=[0])

        imname2=self.img+'.tclean'
        ret2 = tclean(vis=self.refmsfile,imagename=imname2, imsize=200, cell='10.0arcsec', nchan=3, spw='0:0~2', pblimit=0.1, interpolation='nearest',specmode='cube',niter=15, cycleniter=5, gridder='mosaic',mosweight=False, deconvolver='multiscale',scales=[0])

        outimname1 = imname1+'.joint.cube'

        report=th.checkall(imgexist=[outimname1+'.psf', outimname1+'.image',
                                     imname2+'.psf', imname2+'.image'], 
                           imgval=[(outimname1+'.psf', 1.0, [100,100,0,0]),
                                   (imname2+'.psf', 1.0, [100,100,0,0]),
                                   (outimname1+'.residual', 0.809179, [100,100,0,0]),  ## End of minor cycle : 0.818269. Changes to 0.809179 after major cycle. 
                                   (imname2+'.residual', 0.809179, [100,100,0,0]),
                                   (outimname1+'.image', 1.3766, [100,100,0,0]),
                                   (imname2+'.image', 1.3766, [100,100,0,0]),
                                   (outimname1+'.image', 1.3561, [100,100,0,1]),
                                   (imname2+'.image', 1.3561, [100,100,0,1]) ])
        self.checkfinal(pstr=report)

    #Test17
    #17. Single pointing test with INT-only data from refim_point.ms : Compare with tclean mtmfs
    #testname: test_intonly_mfs_compare_with_tclean

    def test_intonly_mfs_compare_with_tclean(self):
        """ [Compare] Test_intonly_mfs_compare_with_tclean """
        inputdata={'refmsname':'refim_point.ms'}
        self.prepData(inputdata=inputdata)
        imname1=self.img+'.sdint'
        ret1 = sdintimaging(vis=self.refmsfile,imagename=imname1, usedata='int', imsize=200, cell='10.0arcsec', nchan=5, reffreq='1.5GHz', start='1.0GHz',width='200.0MHz', interpolation='nearest',specmode='mfs',niter=10, cycleniter=5, gridder='standard', deconvolver='mtmfs',scales=[0])

        imname2=self.img+'.tclean'
        ret2 = tclean(vis=self.refmsfile,imagename=imname2, imsize=200, cell='10.0arcsec', nchan=5, reffreq='1.5GHz', specmode='mfs',niter=10, cycleniter=5, gridder='standard', deconvolver='mtmfs',nterms=2, scales=[0])

        outimname1 = imname1+'.joint.multiterm'

        report=th.checkall(imgexist=[outimname1+'.psf.tt0', outimname1+'.image.tt0',
                                     imname2+'.psf.tt0', imname2+'.image.tt0'], 
                           imgval=[(outimname1+'.psf.tt0', 1.0, [100,100,0,0]),
                                   (imname2+'.psf.tt0', 1.0, [100,100,0,0]),
                                   (outimname1+'.image.tt0', 1.04, [100,100,0,0]),
                                   (imname2+'.image.tt0', 1.04, [100,100,0,0]),
                                   (outimname1+'.alpha', -1.06, [100,100,0,0]),
                                   (imname2+'.alpha', -1.06, [100,100,0,0]) ])
        self.checkfinal(pstr=report)



def suite():
    return[test_singlepointing,test_mosaic,test_compare_sdint_tclean]

if __name__ == '__main__':
    unittest.main() 

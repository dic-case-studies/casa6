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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here (TBD):
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_sdintimaging/about
#
#
##########################################################################


####    Imports     ####
import os
import sys
import unittest
import shutil
import numpy as np
from casatestutils.imagerhelpers import TestHelpers

CASA6 = False
try: 
    from casatools import ctsys
    from casatasks import casalog, sdintimaging
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
    sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
    #from testhelper_imager import TestHelpers
    #import casaTestHelper as th
    #refdatapath = '/export/home/murasame2/casadev/imagerRefact/sdint/orig_scripts/WidebandSDINT/Data'
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from parallel.parallel_task_helper import ParallelTaskHelper
    from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
    #from imagerhelpers.testhelper_imager import TestHelpers
    #import casaTestHelper as th

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
        #if ParallelTaskHelper.isMPIEnabled():
        #    self.parallel = True
        #    self.PH = PyParallelImagerHelper()
        #    self.nnode = len(self.PH.getNodeList())
        #self.th = TestHelpers()
 

    def tearDown(self):
        # delete all (input and output data)
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
        self.pblimit=-0.1
        # orginal runtest parameter
        #self.niter=1000
        #self.cycleniter= 200
        # for quicker execution 
        self.niter=100
        self.cycleniter= 20
        self.interpolation='nearest'

    # Test 1
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'], imgval=[(outimg+'.psf.tt0', 0.991, [400,400,0,0]),(outimg+'.image.tt0', 1.1820, [350,433,0,0])])

    #Test 2
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
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'], imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),(outimg+'.image.tt0', 1.0932, [350,433,0,0])])

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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'], imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),(outimg+'.image.tt0', 17.199, [311,356,0,0])])
    #Test 3
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'], imgval=[(outimg+'.psf.tt0', 1.0, [400,400,0,0]),(outimg+'.image.tt0', 17.199, [311,356,0,0])])

    #Test4
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test5
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
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test6
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)

        outimg = self.img+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

class test_mosaic(testref_base):

    def setUp(self):
        # common parameters
        super(test_mosaic, self).setUp()
        self.imsize=1500
        pblimit=0.1
        self.cell = '9.0arcsec'
        self.phasecenter = 'J2000 19:59:28.500 +40.44.01.50'
        self.nchan=3
        self.reffreq='1.5GHz'
        self.scales=[0,12,20,40,60,80,100]
        self.pblimit=-0.1
        # orginal runtest parameter
        #self.niter=1000
        #self.cycleniter= 200
        # for quicker execution 
        self.niter=100
        self.cycleniter= 20
        self.interpolation='nearest'

     #Test7
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)
        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test8
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
        ret = sdintimaging(usedata='int', vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)
        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test9
    def test_mosaic_mfs_sdonly(self):
        # Equivalent to onetest(runtype='Mosaic', specmode='mfs', usedata='sd')
        """ [Mosaic] Test_mosaic_mfs_sd """
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)
        outimg = self.img+'.joint.multiterm'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test10
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
        ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)
        outimg = self.img+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test11
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
        ret = sdintimaging(usedata='int', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)
        outimg = self.img+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])

    #Test12
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
        ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='cube', gridder='mosaic', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, mask=self.mask, interactive=0,parallel=self.parallel)
        outimg = self.img+'.joint.cube'
        report=th.checkall(imgexist=[outimg+'.psf.tt0', outimg+'.residual.tt0', outimg+'.image.tt0', outimg+'.image.tt1',outimg+'.alpha'])


def suite():
    return[test_singlepointing,test_mosaic]

if __name__ == '__main__':
    unittest.main() 

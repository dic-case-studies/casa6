from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import unittest
import inspect
import numpy as np
import operator

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
     from casatools import ctsys, quanta, measures, image, vpmanager, calibrater
     from casatasks import casalog, delmod, imsubimage, tclean, uvsub, imhead, imsmooth, immath, widebandpbcor
     from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
     from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

     sys.path.append(os.path.abspath(os.path.basename(__file__)))
     from testhelper_imager import TestHelpers

     _ia = image( )
     _vp = vpmanager( )
     _cb = calibrater( )
     _qa = quanta( )
     _me = measures( )
     #refdatapath = ctsys.resolve('regression/unittest/clean/refimager')
     refdatapath = '/export/home/murasame2/casadev/imagerRefact/sdint/orig_scripts/WidebandSDINT/Data'
     
else:
     from __main__ import default
     from tasks import *
     from taskinit import *
     from parallel.parallel_task_helper import ParallelTaskHelper
     from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
     from imagerhelpers.testhelper_imager import TestHelpers

     _ia = iatool( )
     _vp = vptool( )
     _cb = cbtool( )
     # not local tools
     _qa = qa
     _me = me
     #refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/clean/refimager'
     refdatapath = '/export/home/murasame2/casadev/imagerRefact/sdint/orig_scripts/WidebandSDINT/Data'

     
## List to be run
def suite():
     return [test_singlepointing]
#     return [test_onefield, test_iterbot, test_multifield,test_stokes,test_cube, test_widefield,test_mask, test_modelvis,test_startmodel,test_widefield_failing]
 
## Base Test class with Utility functions
class testref_base(unittest.TestCase):

     def setUp(self):
          self.epsilon = 0.05
          self.msfile = ""
          self.img = "tst"
          self.sdimage = ""
          self.sdpsf = ""
          self.mask = ""
          self.cfcache = 'cfcach'
          # To use subdir in the output image names in some tests (CAS-10937)
          self.img_subdir = 'refimager_tst_subdir'
          self.parallel = False
          self.nnode = 0
          if ParallelTaskHelper.isMPIEnabled():
              self.parallel = True
              self.PH = PyParallelImagerHelper()
              self.nnode = len(self.PH.getNodeList())

          self.th = TestHelpers()

     def tearDown(self):
          """ don't delete it all """
#          self.delData()

     # Separate functions here, for special-case tests that need their own MS.
     def prepData(self,inputdata={}):
          os.system('rm -rf ' + self.img_subdir)
          os.system('rm -rf ' + self.img+'*')
          if inputdata!={}:
              if inputdata.has_key('msname'):
                  self.msfile=inputdata['msname']
                  if (os.path.exists(self.msfile)):
                      os.system('rm -rf ' + self.msfile)
                  shutil.copytree(os.path.join(refdatapath,self.msfile), self.msfile)
              if inputdata.has_key('sdimage'):
                  self.sdimage=inputdata['sdimage']
                  if (os.path.exists(self.sdimage)):
                      os.system('rm -rf ' + self.sdimage)
                  shutil.copytree(os.path.join(refdatapath,self.sdimage), self.sdimage)
              if inputdata.has_key('sdpsf'):
                  self.sdpsf=inputdata['sdpsf']
                  if (os.path.exists(self.sdpsf)):
                      os.system('rm -rf ' + self.sdpsf)
                  shutil.copytree(os.path.join(refdatapath,self.sdpsf), self.sdpsf)
              if inputdata.has_key('mask'):
                  self.mask=inputdata['mask']
                  if (os.path.exists(self.mask)):
                      os.system('rm -rf ' + self.mask)
                  origmask = os.path.join(refdatapath,self.mask)
                  if os.path.exists(origmask):
                      if os.path.isfile(origmask):
                          shutil.copyfile(origmask, self.mask)
                      else:
                          shutil.copytree(origmask, self.mask)
                  #shutil.copytree(os.path.join(refdatapath,self.mask), self.mask)
          
     def prepCfcache(self,cfcache=""):
         if (os.path.exists(self.cfcache)):
               os.system('rm -rf ' + self.cfcache)
         if cfcache!="":
               self.cfcache=cfcache 
         if (os.path.exists(self.cfcache)):
               os.system('rm -rf ' + self.cfcache)
         shutil.copytree(os.path.join(refdatapath,self.cfcache), self.cfcache)

     def delData(self,msname=""):
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.cfcache)):
               os.system('rm -rf ' + self.cfcache)
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          os.system('rm -rf ' + self.img_subdir)
          os.system('rm -rf ' + self.img+'*')

     def prepInputmask(self,maskname=""):
          if maskname!="":
              self.maskname=maskname
          if (os.path.exists(self.maskname)):
              os.system('rm -rf ' + self.maskname)
          shutil.copytree(os.path.join(refdatapath,self.maskname), self.maskname)

     def checkfinal(self,pstr=""):
          #pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  casa -c `echo $CASAPATH | awk '{print $1}'`/gcwrap/python/scripts/regressions/admin/runUnitTest.py test_refimager["+ inspect.stack()[1][3] +"]"
          pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  runUnitTest.main(['test_tclean["+ inspect.stack()[1][3] +"]'])"
          casalog.post(pstr,'INFO')
          if( pstr.count("(Fail") > 0 ):
               self.fail("\n"+pstr)

##############################################
## Unit test for sdintimaging task
class test_singlepointing(testref_base):

     def setUp(self):
          # common parameters
          super(test_singlepointing, self).setUp()
          self.imsize=800
          self.cell = '9.0arcsec'
          self.phasecenter = 'J2000 19:59:28.500 +40.44.01.50'
          self.nchan=3
          self.reffreq='1.5GHz'
          self.scales=[0,12,20,40,60,80,100]
          self.pblimit=-0.1
          self.niter=1000
          self.cycleniter= 200
          self.interpolation='nearest'
          
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
          ret = sdintimaging(usedata='sdint', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, interactive=0,parallel=self.parallel)

          #report=self.th.checkall(imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'], imval=[(self.img+'.psf', 1.0, [50,50,0,0])])
          #self.checkfinal(pstr=report)


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
          # iterations may need to be shorten for the final version of test
          self.prepData(inputdata=inputdata)
          ret = sdintimaging(usedata='int', vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, interactive=0,parallel=self.parallel)

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
          ret = sdintimaging(usedata='sd', sdimage=self.sdimage, sdpsf=self.sdpsf, vis=self.msfile,imagename=self.img,imsize=self.imsize,cell=self.cell,phasecenter=self.phasecenter, specmode='mfs', gridder='standard', nchan=self.nchan, reffreq=self.reffreq, pblimit=self.pblimit,interpolation=self.interpolation, deconvolver=deconvolver, scales=self.scales, niter=self.niter, cycleniter=self.cycleniter, interactive=0,parallel=self.parallel)

if is_CASA6:
     if __name__ == '__main__':
          unittest.main()

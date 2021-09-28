##########################################################################
##########################################################################
#
# Test programs for the refactored imager :  test_tclean
#
# Each of the following categories (classes) has a set of tests within it.
#
#  test_onefield                 # basic tests, deconvolution algorithms
#  test_iterbot                   # iteration control options for mfs and cube
#  test_multifield               # multiple fields of same type and with different shapes/deconvolvers/gridders
#  test_stokes                    # multiple stokes planes, imaging with flagged correlations..
#  test_cube                      # all things cube. Spectral frame setup, handling empty channels, etc
#  test_widefield                # facets, wprojection, imagemosaic, mosaicft, awproject
#  test_mask                      # input mask options : regridding, mask file, automasking, etc
#  test_modelvis                # saving models (column/otf), using starting models, predict-only (setjy)
#  test_ephemeris                # ephemeris tests for gridder standard and mosaic, mode mfs and cubesource
#
# To run from within casa 5:  
#
#  runUnitTest.main(['test_tclean'])                                              # Run all tests
#  runUnitTest.main(['test_tclean[test_onefield]'])                               # Run tests from test_onefield
#  runUnitTest.main(['test_tclean[test_onefield_mtmfs]'])                         # Run one specific test
#  runUnitTest.main(['test_tclean[test_onefield_mtmfs,test_onefield_hogbom]'])    # Multiple specific tests
#
# To see the full list of tests :   grep "\"\"\" \[" test_tclean.py
#
#  These tests need data stored in casatestdata/unittest/tclean
#  The datasets are symliked to the above directory. If using cp to copy them locally,
#  Use cp -RH
#
#  For a developer build, to get the datasets locally 
#
#  --- Get the test data repo :  svn co https://svn.cv.nrao.edu/svn/casatestdata casatestdata
#  --- Add a link to casatestdata inside $CASAPATH
# ########################################################################
# SKIPPED TESTS 
# More tests were added to skip (as of 2019,04,26)
#
# (as of 2019.02.05 - Seven tests total)
# The following tests are currently skipped as the supports of the particular
# modes are not available in parallel mode yet
# =>   
#     test_multifield_both_cube_diffshape:
#     test_multifield_cube_mfs
#     test_multifield_cube_mtmfs
#     test_cube_21
#
# The following tests in pricipal should be working but curently broken 
# until fixes to test or code are properly made.
# =>  test_multifield_facets_mfs
#     test_multifield_facets_mtmfs
#     test_cube_D1
# 
# Added to skip at least for 5.5
#     test_cube_chanchunks
#     test_cube_chanchunks_savemodel (possible race conditions)
#     test_modelvis_2 (possible race conditions)
#     test_modelvis_3 (possible race conditions)
#     test_modelvis_5 (possible race conditions)
#     test_modelvis_6 (possible race conditions)
#     test_modelvis_7 (possible race conditions)
#     test_modelvis_8 (possible race conditions)
#     test_modelvis_9 (possible race conditions)
#     test_modelvis_10 (possible race conditions)
#     test_modelvis_11 (possible race conditions)
#     test_startmodel_with_mask_mfs(possible race conditions)
#     test_startmodel_with_mask_mtmfs(possible race conditions)

#Ressurected from skipping after some fixes
#     test_mask_5
#     test_iterbot_cube_2
#     test_multifield_both_cube
##########################################################################
#
#  Datasets
#
#  refim_twochan.ms : 2 channels, one 1Jy point source with spectral index of -1.0
#  refim_twopoints_twochan.ms : Two point sources, 1Jy and 5Jy, both with spectral index -1.0. For multifield tests.
#  refim_point.ms : 1-2 GHz, 20 channels, 1 spw, one 1Jy point source with spectral index -1.0.
#  refim_point_withline.ms : refim_point with a 'line' added into 3 channels (just topo)
#  refim_mawproject.ms : Two pointing wideband mosaic with 1 point source in between the two pointings
#  refim_mawproject_offcenter.ms : Two pointing wideband mosaic with 1 point source at center of one pointing
#  refim_point_stokes.ms : RR=1.0, LL=0.8, RL and LR are zero. Stokes I=0.9, V=0.1, U,Q=0.0
#  refim_point_linRL.ms : I=1, Q=2, U=3, V=4  in circular pol basis.
#  venus_ephem_test.ms : 7-point mosaic of Venus (ephemeris), Band 6, 1 spw, averaged to 1 chan
#
##########################################################################

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
     from casatasks import casalog, delmod, imsubimage, tclean, uvsub, imhead, imsmooth, immath, widebandpbcor, impbcor, flagdata, makemask
     from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
     from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
     from casatasks import impbcor

     from casatestutils.imagerhelpers import TestHelpers

     _ia = image( )
     _vp = vpmanager( )
     _cb = calibrater( )
     _qa = quanta( )
     _me = measures( )
     
     refdatapath = ctsys.resolve('unittest/tclean/')
     #refdatapath = "/export/home/riya/rurvashi/Work/ImagerRefactor/Runs/UnitData"
     #refdatapath = "/home/vega/rurvashi/TestCASA/ImagerRefactor/Runs/WFtests"
else:
     from __main__ import default
     from tasks import *
     from taskinit import *
     from parallel.parallel_task_helper import ParallelTaskHelper
     from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
     #from imagerhelpers.testhelper_imager import TestHelpers

     _ia = iatool( )
     _vp = vptool( )
     _cb = cbtool( )
     # not local tools
     _qa = qa
     _me = me

     refdatapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/tclean/'
     #refdatapath = "/export/home/riya/rurvashi/Work/ImagerRefactor/Runs/UnitData"
     #refdatapath = "/home/vega/rurvashi/TestCASA/ImagerRefactor/Runs/WFtests"
 
from casatestutils.imagerhelpers import TestHelpers

## List to be run
def suite():
     return [test_onefield, test_iterbot, test_multifield,test_stokes, test_modelvis, test_cube, test_mask, test_startmodel, test_widefield, test_pbcor, test_mosaic_mtmfs, test_mosaic_cube, test_ephemeris, test_hetarray_imaging, test_wproject, test_errors_failures]
 
## Base Test class with Utility functions
class testref_base(unittest.TestCase):

     def setUp(self):
          self.epsilon = 0.05
          self.msfile = ""
          self.img = "tst"
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
          self.check_final = self.th.check_final

     def tearDown(self):
          """ don't delete it all """
#          self.delData()

     # Separate functions here, for special-case tests that need their own MS.
     def prepData(self,msname=""):
          os.system('rm -rf ' + self.img_subdir)
          os.system('rm -rf ' + self.img+'*')
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          shutil.copytree(os.path.join(refdatapath,self.msfile), self.msfile)
          
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

##############################################
##############################################

##Task level tests : one field, 2chan.
class test_onefield(testref_base):
     
     def test_onefield_defaults(self):
          """ [onefield] Test_Onefield_defaults : Defaults """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',interactive=0,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'], imgval=[(self.img+'.psf', 1.0, [50,50,0,0])])
          self.assertTrue(self.check_final(pstr=report))

     def test_onefield_clark(self):
          """ [onefield] Test_Onefield_clark : mfs with clark minor cycle """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='clark',interactive=0,parallel=self.parallel) #,phasecenter='J2000 19h59m57.5s +40d49m00.077s') # default is clark
          #off center#ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=1000,interactive=0,phasecenter='J2000 19h59m57.5s +40d49m00.077s') # default is clark
          #compare with clean#clean(vis=self.msfile,imagename=self.img+'.old',imsize=200,cell='8.0arcsec',niter=1000,psfmode='clark',phasecenter='J2000 19h59m57.5s +40d49m00.077s') # default is clark
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imgval=[(self.img+'.psf',1.0,[50,50,0,0])])
          self.assertTrue(self.check_final(pstr=report))

     def test_onefield_hogbom(self):
          """ [onefield] Test_Onefield_hogbom : mfs with hogbom minor cycle """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)#,phasecenter='J2000 19h59m57.5s +40d49m00.077s')
          report=self.th.checkall(ret=ret, peakres=0.35, modflux=0.77, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imgval=[(self.img+'.psf',1.0,[50,50,0,0])]  , tfmask=[(self.img+'.image',['mask0']), (self.img+'.pb',['mask0'])] )
          self.assertTrue(self.check_final(pstr=report))

     def test_onefield_mem(self):
          """ [onefield] Test_Onefield_mem : mfs with mem minor cycle """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=10,deconvolver='mem',interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=12.7, modflux=6.98, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imgval=[(self.img+'.psf',1.0,[100,100,0,0])])
          self.assertTrue(self.check_final(pstr=report))

     def test_onefield_multiscale(self):
          """ [onefield] Test_Onefield_multiscale : mfs with multiscale minor cycle """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=10,deconvolver='multiscale',scales=[0,20,40,100],interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.823, modflux=3.816, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imgval=[(self.img+'.psf',1.0,[100,100,0,0])])
          self.assertTrue(self.check_final(pstr=report))

     def test_onefield_mtmfs(self):
          """ [onefield] Test_Onefield_mtmfs : mt-mfs with minor cycle iterations """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10, imgexist=[self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0', self.img+'.model.tt0',self.img+'.model.tt1',self.img+'.alpha'], imgval=[(self.img+'.psf.tt0',1.0,[50,50,0,0]),(self.img+'.psf.tt1',1.039e-05,[50,50,0,0])])
          ## iterdone=11 only because of the return (iterdone_p+1) in MultiTermMatrixCleaner::mtclean() !
          self.assertTrue(self.check_final(pstr=report))

     def test_onefield_autonames(self):
          """ [onefield] Test_Onefield_autonames : Test auto increment of image names """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',parallel=self.parallel)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',restart=False,parallel=self.parallel)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',restart=False,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf',self.img+'_2.psf',self.img+'_3.psf'] )
          self.assertTrue(self.check_final(pstr=report))

     # weighting test
     def test_onefield_weighting(self):
          """ [onefield] Test_Onefield_weighting : mfs with different weighting (natural, uniform, briggs, radial, superuniform)"""
          self.prepData('refim_twochan.ms')
          # default = natural
          ret0 = tclean(vis=self.msfile,imagename=self.img+'0',imsize=100,cell='8.0arcsec',niter=10,weighting='natural', interactive=0,parallel=self.parallel) 
          # uniform
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,weighting='uniform', interactive=0,parallel=self.parallel) 
          report=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image', self.img+'.model'], imgval=[(self.img+'.psf',1.0,[50,50,0,0])])
#          self.assertTrue(self.check_final(pstr=report))

          # briggs r=-2
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,weighting='briggs', robust=-2, interactive=0,parallel=self.parallel)     
          report2=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'2.psf', self.img+'2.residual', self.img+'2.image', self.img+'2.model'], imgval=[(self.img+'2.psf',1.0,[50,50,0,0])])

          # briggs r=0.5(default)
          ret3 = tclean(vis=self.msfile,imagename=self.img+'3',imsize=100,cell='8.0arcsec',niter=10,weighting='briggs', robust=0.5, interactive=0,parallel=self.parallel)     
          report3=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'3.psf', self.img+'3.residual', self.img+'3.image', self.img+'3.model'], imgval=[(self.img+'3.psf',1.0,[50,50,0,0])])

          # briggs r=2
          ret4 = tclean(vis=self.msfile,imagename=self.img+'4',imsize=100,cell='8.0arcsec',niter=10,weighting='briggs', robust=2, interactive=0,parallel=self.parallel)     
          report4=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'4.psf', self.img+'4.residual', self.img+'4.image', self.img+'4.model'], imgval=[(self.img+'4.psf',1.0,[50,50,0,0]),(self.img+'4.sumwt',3430533.5,[0,0,0,0])])

          # radial
          ret5 = tclean(vis=self.msfile,imagename=self.img+'5',imsize=100,cell='8.0arcsec',niter=10,weighting='radial', interactive=0,parallel=self.parallel)     
          report5=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'5.psf', self.img+'5.residual', self.img+'5.image', self.img+'5.model'], imgval=[(self.img+'5.psf',1.0,[50,50,0,0])])

          # superuniform
          ret6 = tclean(vis=self.msfile,imagename=self.img+'6',imsize=100,cell='8.0arcsec',niter=10,weighting='superuniform', interactive=0,parallel=self.parallel)     
          report6=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'6.psf', self.img+'6.residual', self.img+'6.image', self.img+'6.model'], imgval=[(self.img+'6.psf',1.0,[50,50,0,0])])

          # briggs r=0.5(default) with mtmfs (to test SIImageStoreMultiTerm)
          ret7 = tclean(vis=self.msfile,imagename=self.img+'7',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs', weighting='briggs', robust=0.5, interactive=0,parallel=self.parallel)     
          report7=self.th.checkall(ret=ret, peakres=0.263, modflux=0.575, iterdone=10, imgexist=[self.img+'7.psf.tt0', self.img+'7.residual.tt0', self.img+'7.image.tt0', self.img+'7.model.tt0'], imgval=[(self.img+'7.psf.tt0',1.0,[50,50,0,0]),(self.img+'7.psf.tt1',0.0898,[50,50,0,0]),(self.img+'7.sumwt.tt0',1532169.875,[0,0,0,0]),(self.img+'7.sumwt.tt1',137693.875,[0,0,0,0])])


          # beamareas: uniform < briggs-r=-2 < briggs r=0.5 < briggs r=+2 < natural, ...
          # by default, it checks if im1's beam < im2's beam
          print("Test beamarea of tst0.image (natural) is greater than beamarea of tst.image (uniform)")
          self.assertTrue(self.th.check_beam_compare(self.img+'.image', self.img+'0.image'))
          print("Test beamarea of tst2.image (briggs -2) is greater than beamarea of tst.image (uniform)")
          self.assertTrue(self.th.check_beam_compare(self.img+'.image', self.img+'2.image'))
          print("Test beamarea of tst3.image (briggs 0.5) is greater than beamarea of tst2.image (briggs -2))")
          self.assertTrue(self.th.check_beam_compare(self.img+'2.image', self.img+'3.image'))
          print("Test beamarea of tst4.image (briggs 2) is greater than beamarea of tst3.image (briggs 0.5))")
          self.assertTrue(self.th.check_beam_compare(self.img+'3.image', self.img+'4.image'))
     
          self.assertTrue(self.check_final(pstr = report+report2+report3+report4+report5+report6+report7))
          
     def test_onefield_pcwdT_and_pcwdF(self):
     
          self.prepData('refim_twochan.ms')
          
          ret = tclean(self.msfile , imagename=self.img+'1', imsize=20, cell='8.0arcsec', niter=0, nchan=1, spw='0:1', interactive=0, gridder='standard',perchanweightdensity=True,specmode='cube',weighting='briggs',robust=0.5)
	
          ret = tclean(self.msfile , imagename=self.img+'2', imsize=20, cell='8.0arcsec', niter=0, nchan=1, spw='0:1', interactive=0, gridder='standard',perchanweightdensity=False,specmode='cube',weighting='briggs',robust=0.5)
          
          _ia.open(self.img+'1.psf')
          pcwdT_img = _ia.getchunk()
          _ia.close()
        
          _ia.open(self.img+'2.psf')
          pcwdF_img = _ia.getchunk()
          _ia.close()
      
          abs_dif = np.sum(np.abs(pcwdF_img-pcwdT_img))
          _, report1 = self.th.check_val_less_than(abs_dif, 3, valname='abs_dif',testname ="test_onefield_pcwdT_and_pcwdF")

          self.assertTrue(self.check_final(report1))   
     
          
        

     def test_onefield_twoMS(self):
          """ [onefield] Test_Onefield_twoMS : One field, two input MSs, also
          test automatic fallback to 'data' column when no 'corrected' data
          column"""
          ms1 = 'refim_point_onespw0.ms'
          ms2 = 'refim_point_onespw1.ms'
          self.prepData(ms1)
          self.prepData(ms2)
#          try:
#               ## This run should fail with an exception
#               ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10)
#               correct=False
#          except Exception as e:
#              correct=True
#          self.assertTrue(correct)
          ## This run should go smoothly.
          ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf',self.img+'.residual'])
          self.delData(ms1)
          self.delData(ms2)
          self.assertTrue(self.check_final(pstr=report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. onefield with two MSs, briggs weighing. Enable this when CAS011978 is fixed")
     def test_onefield_twoMS_Briggs(self):
          """ [onefield] Test_Onefield_twoMS with Briggs weighting: One field, two input MSs (for verification of CAS-11978)"""
          ms1 = 'refim_point_onespw0.ms'
          ms2 = 'refim_point_onespw1.ms'
          self.prepData(ms1)
          self.prepData(ms2)
#          try:
#               ## This run should fail with an exception
#               ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10)
#               correct=False
#          except Exception as e:
#              correct=True
#          self.assertTrue(correct)
          ## This run should go smoothly. 
          ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10,weighting='briggs', interactive=0, parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.365259, modflux=0.798692, imgexist=[self.img+'.psf',self.img+'.residual'])
          self.delData(ms1)
          self.delData(ms2)
          self.assertTrue(self.check_final(pstr=report))


     def test_onefield_twoMS_diffcolumns(self):
          """ [onefield] Test_Onefield_twoMS_diffcolumns : One field, two input MSs, one with data and one with data and corrected """
          ms1 = 'refim_point_onespw0.ms'
          ms2 = 'refim_point_onespw1.ms'
          self.prepData(ms1)
          self.prepData(ms2)

          ## Make corrected_data column for one of them
          _cb.open(ms2)
          _cb.close()

          ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10,datacolumn='corrected',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf',self.img+'.residual'])
          self.delData(ms1)
          self.delData(ms2)
          self.assertTrue(self.check_final(pstr=report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Erratic in parallel")
     def test_onefield_briggsabs(self):
          """[onefield] test_onefield_briggsabs: """
          self.prepData('refim_point_withline.ms')
          delmod(self.msfile)
          imnat=self.img+"_nat"
          imbriggs0=self.img+"_briggsabs_0"
          imbriggs_2=self.img+"_briggsabs_2"
          imbriggs_2_2=self.img+"_briggsabs_2_2pix"
          retnat = tclean(vis=self.msfile,imagename=imnat,imsize=100,cell='8.0arcsec',specmode='mfs',deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='natural', parallel=self.parallel)
          ret0 = tclean(vis=self.msfile,imagename=imbriggs0,imsize=100,cell='8.0arcsec',specmode='mfs', perchanweightdensity=True,deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='briggsabs', robust=0, noise='1Jy',parallel=self.parallel)
          ret_2=tclean(vis=self.msfile,imagename=imbriggs_2,imsize=100,cell='8.0arcsec',specmode='mfs', perchanweightdensity=True,deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='briggsabs', robust=-2.0, noise='1Jy', parallel=self.parallel)
###          ret_2_1=tclean(vis=self.msfile,imagename=imbriggs_2_2,imsize=100,cell='8.0arcsec',specmode='mfs', perchanweightdensity=True,deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='briggsabs', robust=-2.0, npixels=2, noise='1Jy', parallel=self.parallel)

          self.assertTrue(os.path.exists(imnat+'.image') and os.path.exists(imbriggs0+'.image') and os.path.exists(imbriggs_2+'.image') )
          ###briggsabs 0 should be natural
          self.assertTrue(self.th.check_beam_compare(imbriggs0+'.image', imnat+'.image', operator.eq))
          self.assertTrue(self.th.check_beam_compare(imbriggs_2+'.image', imbriggs0+'.image'))
        ###  self.assertTrue(self.th.check_beam_compare(imbriggs_2_2+'.image', imbriggs_2+'.image', operator.le))

     def test_onefield_restart_mfs(self):
          """ [onefield] : test_onefield_restart_mfs : Check calcpsf,calcres and ability to restart and continue"""
          ## TODO : Need to add and use info in the return record, when only a major cycle is done. Then check nmajorcycle.
          self.prepData('refim_twochan.ms')

          ## Only psf
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=True,calcres=False,deconvolver='clark',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf'], imgexistnot=[self.img+'.residual', self.img+'.image'],nmajordone=1)

          ## Only residual
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=False,calcres=True,deconvolver='clark',restoration=False,parallel=self.parallel)
          report1=self.th.checkall(imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'],nmajordone=1)

          ## Start directly with minor cycle and do only the last major cycle.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='clark',parallel=self.parallel)
          report2=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, imgexist=[self.img+'.psf',self.img+'.residual', self.img+'.image'],nmajordone=1)

          ## Re-start from existing model image and continue on...
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='clark',parallel=self.parallel)
          report3=self.th.checkall(ret=ret, peakres=0.161, modflux=0.991, imgexist=[self.img+'.psf',self.img+'.residual', self.img+'.image'],nmajordone=1)

          self.assertTrue(self.check_final(pstr=report+report1+report2+report3))


     def test_onefield_restart_mtmfs(self):
          """ [onefield] : test_onefield_restart_mtmfs : Check calcpsf,calcres and ability to restart and continue"""
          ## TODO : Need to add and use info in the return record, when only a major cycle is done. Then check nmajorcycle.
          self.prepData('refim_twochan.ms')

          ## Only psf
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=True,calcres=False,deconvolver='mtmfs',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf.tt0', self.img+'.psf.tt1'], imgexistnot=[self.img+'.residual.tt0', self.img+'.image.tt0'],nmajordone=1)

          ## Only residual
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=False,calcres=True,deconvolver='mtmfs',restoration=False,parallel=self.parallel)
          report1=self.th.checkall(imgexist=[self.img+'.psf.tt0',self.img+'.psf.tt1', self.img+'.residual.tt0', self.img+'.residual.tt1'], imgexistnot=[self.img+'.image.tt0'],nmajordone=1)

          ## Start directly with minor cycle and do only the last major cycle.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='mtmfs',parallel=self.parallel)
          report2=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, imgexist=[self.img+'.psf.tt1',self.img+'.residual.tt1', self.img+'.image.tt1', self.img+'.alpha'],nmajordone=1,imgval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          ## Re-start from existing model image and continue on...
          ## ( If restart from modified residuals... the alpha is -1.25xx which is wrong. 
          ##   In this case, need to do calcres=True which will do extra first major cycle (nmajor=2) )
          ## But... current code (as of r33373) makes appropriate restored image but does not mess up residuals.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='mtmfs',parallel=self.parallel)
          report3=self.th.checkall(ret=ret, peakres=0.136, modflux=0.988, imgexist=[self.img+'.psf.tt1',self.img+'.residual.tt1', self.img+'.image.tt1', self.img+'.alpha'],nmajordone=1,imgval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          ### Calcres=True and restart (to test CAS-10337)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=True,deconvolver='mtmfs',parallel=self.parallel)
          report4=self.th.checkall(ret=ret, peakres=0.0477, modflux=1.077, imgexist=[self.img+'.psf.tt1',self.img+'.residual.tt1', self.img+'.image.tt1', self.img+'.alpha'],nmajordone=2,imgval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          self.assertTrue(self.check_final(pstr=report+report1+report2+report3+report4))

     def test_onefield_all_outputs_mfs(self):
          """ [onefield] : test_onefield_all_outputs_mfs : Make all output images even when not needed """
          self.prepData('refim_twochan.ms')

          ## Make only partial outputs
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='hogbom')
#          report1=self.th.checkall(imgexist=[self.img+'.psf', self.img+'.residual'],imgexistnot=[self.img+'.image',self.img+'.model'],nmajordone=1)

          ## Make all outputs
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='hogbom',restoration=True,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'2.psf', self.img+'2.residual',self.img+'2.image',self.img+'2.model'],nmajordone=1)
 
          self.assertTrue(self.check_final(pstr = report2))

     def test_onefield_all_outputs_mtmfs(self):
          """ [onefield] : test_onefield_all_outputs_mtmfs : Make all output images even when not needed """
          self.prepData('refim_twochan.ms')

          ## Make only partial outputs
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs')
#          report1=self.th.checkall(imgexist=[self.img+'.psf.tt0', self.img+'.psf.tt1'],imgexistnot=[self.img+'.image.tt0',self.img+'.model.tt0', self.img+'.alpha'],nmajordone=1)

          ## Make all outputs
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',restoration=True,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'2.psf.tt0', self.img+'2.psf.tt1',self.img+'2.image.tt0',self.img+'2.model.tt0', self.img+'2.alpha'],nmajordone=1)
 
          self.assertTrue(self.check_final(pstr=report2))

     def test_onefield_restore_mtmfs_niter0(self):
          """ [onefield] : test_onefield_restore_mtmfs_niter0 : Niter=0 run followed by restoration without a model"""
          self.prepData('refim_twochan.ms')

          ## This test also checks the principal solution calculation on the dirty images.

          ## niter=0 run 
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',restoration=False,parallel=self.parallel)
          report1=self.th.checkall(imgexist=[self.img+'.psf.tt0', self.img+'.psf.tt1'], imgexistnot=[self.img+'.model.tt0', self.img+'.model.tt0'],nmajordone=1)
          ## restore only 
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',calcres=False,calcpsf=False,restoration=True,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.alpha'],nmajordone=0,
                             imgval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          ## niter=0 and restore ( in one step )
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',restoration=True,parallel=self.parallel)
          report3=self.th.checkall(imgexist=[self.img+'2.image.tt0', self.img+'2.alpha'],nmajordone=1,
                             imgval=[(self.img+'.alpha',-1.0,[50,50,0,0])] )

          self.assertTrue(self.check_final(pstr=report1+report2+report3))

     def test_onefield_rectangular_pixels(self):
          """ [onefield] : test_onefield_rectangular_pixels : Test restoration with rectangular pixels (cas-7171)"""
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell=['10.0arcsec','30.0arcsec'],niter=10,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'], imgval=[(self.img+'.image',0.482,[50,49,0,0])] )

          self.assertTrue(self.check_final(report))


     def test_onefield_mtmfs_2spws_2MSs(self):
          """ [onefield] Test_onefield_mtmfs_2spws_2MSs : MT-MFS on multi-spws in separate MSs, to test default reffreq and coordinate system generation (CAS-9518) """

          ms1 = 'refim_point_onespw0.ms'
          ms2 = 'refim_point_onespw1.ms'
          self.prepData(ms1)
          self.prepData(ms2)
          
          ret = tclean(vis=[ms1,ms2],imagename=self.img,imsize=100,cell='8.0arcsec',
                       interactive=0,niter=10,deconvolver='mtmfs',parallel=self.parallel)

          checkims = [self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0',self.img+'.model.tt0']  
          
          ## For parallel run, check sub workdirectory images also.
          #if self.parallel==True:
          #     checkims = checkims + self.th.getNParts( imprefix=self.img, 
          #                                              imexts=['residual.tt0','residual.tt1',
          #                                                      'psf.tt0','psf.tt1',
          #                                                      'model.tt0','model.tt1']) 
          report = self.th.checkall(ret=ret, 
                                     peakres=0.409, modflux=0.764, iterdone=10, nmajordone=2,
                                     imgexist=checkims, 
                                     imgval=[(self.img+'.alpha',-2.0,[50,50,0,0]),
                                            (self.img+'.sumwt.tt0', 94050.05,[0,0,0,0]) ,
                                            (self.img+'.sumwt.tt1', 0.006198,[0,0,0,0]) ], 
                                     reffreq= [(self.img+'.image.tt0',1489984775.68)] )
          
          self.assertTrue(self.check_final(report))

     def test_onefield_mtmfs_nterms1(self):
          """ [onefield] Test_Onefield_mtmfs_nterms1 : mt-mfs with nterms=1 (CAS-11364, CAS-11367) """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',nterms=1,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.369, modflux=0.689, iterdone=10, imgexist=[self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0', self.img+'.model.tt0'], imgval=[(self.img+'.psf.tt0',1.0,[50,50,0,0]),(self.img+'.image.tt0',1.05,[50,50,0,0])])
          ## iterdone=11 only because of the return (iterdone_p+1) in MultiTermMatrixCleaner::mtclean() !
          self.assertTrue(self.check_final(pstr=report))
          
     def test_onefield_mtmfs_smallscalebias(self):
          """ [onefield] Test_Onefield_mtmfs : mt-mfs with minor cycle iterations and smallscalebias = 0.9 """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=10,deconvolver='mtmfs',nterms=1,interactive=0,parallel=self.parallel,smallscalebias=0.9,scales=[0,20,40,100])
          report=self.th.checkall(ret=ret, peakres=0.73153, modflux=2.9194, iterdone=10, imgexist=[self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0', self.img+'.model.tt0'],imgval=[(self.img+'.image.tt0',0.526,[100,100,0,0])])
          self.assertTrue(self.check_final(pstr=report))   

     def test_onefield_gridders(self):
          """ [onefield] Test_Onefield_gridders : Check all single field gridder equivalent names are accepted """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',gridder='ft', interactive=0,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'], imgval=[(self.img+'.psf', 1.0, [50,50,0,0])])
          ret2 = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',gridder='gridft', interactive=0,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'], imgval=[(self.img+'.psf', 1.0, [50,50,0,0])])
          self.assertTrue(self.check_final(pstr=report+report2))


     def test_onefield_cube_restoringbeam(self):
          """ [onefield] Test explicit restoring beams for cube : Test peak flux with niter=0, compared with smoothing vs restoringbeam"""
          
          self.prepData('refim_point.ms')
          
          ret1 = tclean(vis=self.msfile,imagename=self.img,
                        imsize=100,cell='10.0arcsec',interpolation='nearest',
                        interactive=0,niter=0,specmode='cube',
                        parallel=self.parallel)
          imsmooth(imagename=self.img+'.image', targetres=True, major='120.0arcsec', minor='120.0arcsec', pa='0deg',outfile=self.img+'.smoothed.image',overwrite=True)

          ret2 = tclean(vis=self.msfile,imagename=self.img+'.rest',
                        imsize=100,cell='10.0arcsec',interpolation='nearest',
                        interactive=0,niter=0,specmode='cube',restoringbeam='120.0arcsec',
                        parallel=self.parallel)
          
          header = imhead(self.img+'.rest.image',verbose=False)
          
          estr = "["+inspect.stack()[1][3]+"] Has single restoring beam ? : " + self.th.verdict('restoringbeam' in header) + "\n"

          report = self.th.checkall(imgexist=[self.img+'.rest.image'],
                                     imgval=[(self.img+'.image',1.36,[50,50,0,2]),
                                            (self.img+'.smoothed.image',1.54,[50,50,0,2]),
                                            (self.img+'.rest.image',1.54,[50,50,0,2]),
                                            (self.img+'.image',0.79,[50,50,0,18]),
                                            (self.img+'.smoothed.image',1.21,[50,50,0,18]),
                                            (self.img+'.rest.image',1.21,[50,50,0,18])  ])

          ## Note : In this test, setting niter=2000 will get all the runs to produce the same correct values.
          
          ## Pass or Fail (and why) ?
          self.assertTrue(self.check_final(estr+report))


     def test_onefield_mtmfs_restoringbeam(self):
          """ [onefield] Test explicit restoring beams for mtmfs : Test peak flux with niter=0, compared with smoothing vs   restoringbeam"""
          
          self.prepData('refim_point.ms')
          
          ret1 = tclean(vis=self.msfile,imagename=self.img,
                        imsize=100,cell='10.0arcsec',
                        interactive=0,niter=0,specmode='mfs', deconvolver='mtmfs',
                        parallel=self.parallel)

          imsmooth(imagename=self.img+'.image.tt0', targetres=True,
                   major='120.0arcsec', minor='120.0arcsec', pa='0deg',
                   outfile=self.img+'.smoothed.image.tt0',overwrite=True)
          imsmooth(imagename=self.img+'.image.tt1', targetres=True,
                   major='120.0arcsec', minor='120.0arcsec', pa='0deg',
                   outfile=self.img+'.smoothed.image.tt1',overwrite=True)
          #os.system('rm -rf '+code+'trest_1_smoothed.alpha')
          immath(imagename=[self.img+'.smoothed.image.tt0', self.img+'.smoothed.image.tt1'],
                 mode='evalexpr', expr='IM1/IM0',outfile=self.img+'.smoothed.alpha')

          ret2 = tclean(vis=self.msfile,imagename=self.img+'.rest',
                        imsize=100,cell='10.0arcsec',
                        interactive=0,niter=0,specmode='mfs', deconvolver='mtmfs',restoringbeam='120.0arcsec',
                        parallel=self.parallel)
          
          report = self.th.checkall(imgexist=[self.img+'.rest.image.tt0'],
                                     imgval=[(self.img+'.image.tt0',1.06,[50,50,0,0]),
                                            (self.img+'.alpha',-1.03,[50,50,0,0]),
                                            (self.img+'.smoothed.image.tt0',1.5,[50,50,0,0]),
                                            (self.img+'.smoothed.alpha',-2.19,[50,50,0,0]),
                                            (self.img+'.rest.image.tt0',1.5,[50,50,0,0]),
                                            (self.img+'.rest.alpha',-2.19,[50,50,0,0])  ])

          ## Note : In this test, setting niter=100 will get all the runs to produce the same, correct alpha=-1.0

          ## Pass or Fail (and why) ?
          self.assertTrue(self.check_final(report))

     def test_onefield_projections(self):
          """ [onefield] Test_Onefield_projections : test selected projections  """
          self.prepData('refim_twochan.ms')
          # default projection = SIN
          ret = tclean(vis=self.msfile,imagename=self.img+'SIN',imsize=100,cell='8.0arcsec',interactive=0,parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'NCP',projection='NCP',imsize=100,cell='8.0arcsec',interactive=0,parallel=self.parallel)
          ret3 = tclean(vis=self.msfile,imagename=self.img+'TAN',projection='TAN',imsize=100,cell='8.0arcsec',interactive=0,parallel=self.parallel)
          ret4 = tclean(vis=self.msfile,imagename=self.img+'ARC',projection='ARC',imsize=100,cell='8.0arcsec',interactive=0,parallel=self.parallel)
          # Current fails with "wcs wcsset_error: Invalid parameter value" for HEALPix
          #ret5 = tclean(vis=self.msfile,imagename=self.img+'HPX',projection='HPX',imsize=100,cell='8.0arcsec',interactive=0,parallel=self.parallel)
          testname=inspect.stack()[0][3]
          report=self.th.checkall(ret=ret, imgexist=[self.img+'SIN.image', self.img+'NCP.image', self.img+'TAN.image',self.img+'ARC.image'], imgval=[(self.img+'SIN.psf',1.0,[50,50,0,0])])
          retSIN = imhead(self.img+"SIN.image", mode='list')
          retNCP = imhead(self.img+"NCP.image", mode='list')
          retTAN = imhead(self.img+"TAN.image", mode='list')
          retARC = imhead(self.img+"ARC.image", mode='list')

          checkimage = "["+testname+"] The image in SIN projection : (" + self.th.verdict(retSIN['projection']=='SIN') + ")\n"
          # in serial 'NCP' is added in projection key but in parallel, this seems to be trancated.
          checkimage += "["+testname+"] The image in NCP projection : (" + self.th.verdict(retNCP['projection'].find('SIN ([0, 1.16122]')==0) + ")\n"
          checkimage += "["+testname+"] The image in TAN projection : (" + self.th.verdict(retTAN['projection']=='TAN') + ")\n"
          checkimage += "["+testname+"] The image in ARC projection : (" + self.th.verdict(retARC['projection']=='ARC') + ")\n"
          
          self.assertTrue(self.check_final(pstr=checkimage+report))
          
     def test_onefield_psf_fit(self):
          """ [onefield] test_onefield_psf_fit : test psf fitting algorithm for different pixels per beam """

          self.prepData('refim_point.ms')

          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell=['28.8arcsec'],niter=0,calcres=False,parallel=self.parallel)
          _ia.open(self.img+'.psf')
          hdr = _ia.summary(list=False)
          _ia.close()
          beam_3ppb = [hdr['restoringbeam']['major']['value'], hdr['restoringbeam']['minor']['value'], hdr['restoringbeam']['positionangle']['value']]
          _, report1 = self.th.check_val(beam_3ppb[0], 55.77472686767578, valname='Beam Major Axis', exact=False)

          self.prepData('refim_point.ms')

          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell=['14.4arcsec'],niter=0,calcres=False,parallel=self.parallel)
          _ia.open(self.img+'.psf')
          hdr = _ia.summary(list=False)
          _ia.close()
          beam_5ppb = [hdr['restoringbeam']['major']['value'], hdr['restoringbeam']['minor']['value'], hdr['restoringbeam']['positionangle']['value']]
          _, report2 = self.th.check_val(beam_5ppb[0], 51.443878173828125, valname='Beam Major Axis', exact=False)

          self.prepData('refim_point.ms')

          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell=['3.6arcsec'],niter=0,calcres=False,parallel=self.parallel)
          _ia.open(self.img+'.psf')
          hdr = _ia.summary(list=False)
          _ia.close()
          beam_20ppb = [hdr['restoringbeam']['major']['value'], hdr['restoringbeam']['minor']['value'], hdr['restoringbeam']['positionangle']['value']]
          _, report3 = self.th.check_val(beam_20ppb[0], 51.55984878540039, valname='Beam Major Axis', exact=False)

          self.assertTrue(self.check_final(report1+report2+report3))

##############################################
##############################################

##Task level tests : iteration controls
class test_iterbot(testref_base):

     def test_iterbot_mfs_1(self):
          """ [iterbot] Test_Iterbot_Mfs_1 : Zero Iterations """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=0,interactive=0,restoration=False,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_2(self):
          """ [iterbot] Test_Iterbot_Mfs_2 : Iterations with low gain """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,gain=0.1,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10,nmajordone=2,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_3(self):
          """ [iterbot] Test_Iterbot_Mfs_3 : Cycleniter test """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,cycleniter=3,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10, nmajordone=5,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_4(self):
          """ [iterbot] Test_Iterbot_Mfs_4 : Iterations with high gain """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10, gain=0.5,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.024, modflux=1.274, iterdone=10, nmajordone=3,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_5(self):
          """ [iterbot] Test_Iterbot_Mfs_5 : Threshold test """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,threshold='0.1Jy',gain=0.5,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.0924, modflux=1.129, iterdone=5, nmajordone=3,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_6(self):
          """ [iterbot] Test_Iterbot_Mfs_6 : Cycleniter and threshold """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10, cycleniter=3, threshold='0.1Jy',gain=0.5,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.0924, modflux=1.129, iterdone=5, nmajordone=3,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_7(self):
          """ [iterbot] Test_Iterbot_Mfs_7 : Threshold + cyclefactor to trigger major cycles earlier """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,threshold='0.01Jy', gain=0.5,cyclefactor=10.0,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.024, modflux=1.274, iterdone=10, nmajordone=9,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_8(self):
          """ [iterbot] Test_Iterbot_Mfs_8 : minpsffraction to trigger major cycles earlier. """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=20,threshold='0.01Jy', minpsffraction = 0.5,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.16127, modflux=0.9919, iterdone=20, nmajordone=4,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_9(self):
          """ [iterbot] Test_Iterbot_Mfs_9 : maxpsffraction """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=20,threshold='0.01Jy', minpsffraction=0.8,maxpsffraction=0.5,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, peakres=0.16127, modflux=0.9919, iterdone=20, nmajordone=4,imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_cube_1(self):
          """ [iterbot] Test_Iterbot_cube_1 : iteration counting across channels (>niter) """
          self.prepData('refim_point_withline.ms')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='clark',niter=10,threshold='0.75Jy',interactive=0,parallel=self.parallel)
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone'])
          else:
            ret=retpar 
          report=self.th.checkall(ret=ret, iterdone=90,nmajordone=2,imgexist=[self.img+'.psf', self.img+'.residual'])
          ## Only chans 6 and 7 reach cycleniter, others reach threshold in fewer than 10 iters per chan.

          self.assertTrue(self.check_final(report))

     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily for 5.5")
     # test_imager_helper issue - now fixed and working
     def test_iterbot_cube_2(self):
          """ [iterbot] Test_Iterbot_cube_2 : High threshold, iterate only on line channels. """
          self.prepData('refim_point_withline.ms')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='clark',niter=10,threshold='1.75Jy',interactive=0,parallel=self.parallel)

          ret={}
          if self.parallel:
            # peakres and modflux is determined from node1 
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'peakres', 'modflux'])
          else:
            ret=retpar 
          report=self.th.checkall(ret=ret, peakres=1.73, modflux=0.407,iterdone=12,nmajordone=2,imgexist=[self.img+'.psf', self.img+'.residual'])

          self.assertTrue(self.check_final(report))


     def test_iterbot_cube_3(self): # test for returned summary/plot for no iteration case 
          """ [iterbot] Test_Iterbot_cube_3 : Very high threshold, no iteration (verification of CAS-8576 fix) """
          self.prepData('refim_point_withline.ms')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='clark',niter=10,threshold='3.5Jy',interactive=0,parallel=self.parallel)
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone'])
          else:
            ret=retpar 
          report=self.th.checkall(ret=ret,iterdone=0,nmajordone=1,imgexist=[self.img+'.psf', self.img+'.residual'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_divergence(self): 
          """ [iterbot] Test_Iterbot_divergence : Use negative loop gain to make it diverge (verification of CAS-9244 fix) """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', niter=50,cycleniter=5, gain=-0.2,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret,iterdone=10,nmajordone=3,imgexist=[self.img+'.psf', self.img+'.image'])

          self.assertTrue(self.check_final(report))

     def test_iterbot_mfs_deconvolvers(self):
          """ [iterbot] : test_iterbot_deconvolvers : Do all minor cycle algorithms respond in the same way to iteration controls ? Now they do. """
          # clark and hogbom reach niter first, but multiscale gets to cyclethreshold first. Check peakres and iterdone.
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,threshold='0.1Jy', interactive=0,deconvolver='clark',parallel=self.parallel)
          report1=self.th.checkall(ret=ret1, peakres=0.3922, modflux=0.732, iterdone=10, nmajordone=2,imgexist=[self.img+'1.psf', self.img+'1.residual', self.img+'1.image'])

          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,threshold='0.1Jy', interactive=0,deconvolver='hogbom',parallel=self.parallel)
          report2=self.th.checkall(ret=ret2, peakres=0.3530, modflux=0.7719, iterdone=10, nmajordone=2,imgexist=[self.img+'2.psf', self.img+'2.residual', self.img+'2.image'])

          ret3 = tclean(vis=self.msfile,imagename=self.img+'3',imsize=100,cell='8.0arcsec',niter=10,threshold='0.1Jy', interactive=0,deconvolver='multiscale',parallel=self.parallel,smallscalebias=0.6)
          report3=self.th.checkall(ret=ret3, peakres=0.3922, modflux=0.7327, iterdone=10, nmajordone=2,imgexist=[self.img+'3.psf', self.img+'3.residual', self.img+'3.image'])
     

          self.assertTrue(self.check_final(report1+report2+report3))

     def test_iterbot_cube_deconvolvers(self):
          """ [iterbot] : test_iterbot_deconvolvers : Do all minor cycle algorithms respond in the same way to iteration controls ? Let's see ! """
          # clark, hogbom and multiscale reach niter. Check peakres and iterdone
          self.prepData('refim_eptwochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=200,cell='10.0arcsec',niter=100,threshold='0.1Jy', interactive=0,deconvolver='clark',specmode='cube',interpolation='nearest',parallel=self.parallel)
          report1=self.th.checkall(ret=ret1, iterdone=114, nmajordone=2,imgexist=[self.img+'1.psf', self.img+'1.residual', self.img+'1.image'],
                                   imgval=[(self.img+'1.image',0.935,[100,100,0,0]), (self.img+'1.image',0.282,[100,100,0,2])  ])

          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=200,cell='10.0arcsec',niter=100,threshold='0.1Jy', interactive=0,deconvolver='hogbom',specmode='cube',interpolation='nearest',parallel=self.parallel)
          report2=self.th.checkall(ret=ret2, iterdone=116, nmajordone=2,imgexist=[self.img+'2.psf', self.img+'2.residual', self.img+'2.image'],
                                   imgval=[(self.img+'2.image',0.939,[100,100,0,0]), (self.img+'2.image',0.282,[100,100,0,2])  ])

          ret3 = tclean(vis=self.msfile,imagename=self.img+'3',imsize=200,cell='10.0arcsec',niter=100,threshold='0.1Jy', interactive=0,deconvolver='multiscale',specmode='cube',interpolation='nearest',parallel=self.parallel,smallscalebias=0.6,scales=[0,6,10,20,40])
          report3=self.th.checkall(ret=ret3, iterdone=100, nmajordone=3,imgexist=[self.img+'3.psf', self.img+'3.residual', self.img+'3.image'], 
                                   imgval=[(self.img+'3.image',0.888,[100,100,0,0]), (self.img+'3.image',0.1601,[100,100,0,2])  ])

          self.assertTrue(self.check_final(report1+report2+report3))

          
     def test_iterbot_cube_tol(self): 
          """ [iterbot] Test_Iterbot_cube_tol :threshold test to allow a tolerance (1/100)  (verification of CAS-11278 fix) """
          self.prepData('refim_point_withline.ms')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='hogbom',niter=1000000,threshold='0.50001Jy',gain=0.1,cycleniter=5,interactive=0,parallel=self.parallel)
           
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'stopcode'])
          else:
            ret=retpar 

          report=self.th.checkall(ret=ret,iterdone=151,nmajordone=4,stopcode=2,imgexist=[self.img+'.psf', self.img+'.residual'])

          self.assertTrue(self.check_final(report))


     def test_iterbot_cube_nsigma(self): 
          """ [iterbot] Test_Iterbot_cube_nsigma : nsigma threshold for cube"""
          self.prepData('refim_point_withline.ms')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='hogbom',niter=1000000,threshold='0.000001Jy', nsigma=10.0, gain=0.5,interactive=0, parallel=self.parallel)
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'stopcode'])
          else:
            ret=retpar 

          report=self.th.checkall(ret=ret,iterdone=407,nmajordone=11,stopcode=8,imgexist=[self.img+'.psf', self.img+'.residual'])

          self.assertTrue(self.check_final(report))


     def test_iterbot_mfs_restart_updatedmask(self): 
          """ [iterbot] Test_mfs_restart_updatedmask : restart mfs with an updated mask (CAS-13508 fix verification) """
          self.prepData('refim_twopoints_twochan.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=512,cell='8.0arcsec',specmode='mfs',deconvolver='hogbom',niter=2,parallel=self.parallel)
          os.system('rm -rf '+self.img+'.mask')

          # Create a new mask around the 1Jy source only. 
          mask_around_dim_source='circle[[256pix,256pix],30pix]'
          makemask(inpimage=self.img+'.residual', inpmask=mask_around_dim_source, output=self.img+'.mask.dim.source',overwrite=True,mode='copy')
          mask_around_dim_source=self.img+'.mask.dim.source'

          # 2nd clean with the updated mask
          retpar=tclean(vis=self.msfile, imagename=self.img,imsize=512,cell='8.0arcsec',specmode='mfs',niter=1,calcres=False, calcpsf=False, restart=True, mask=mask_around_dim_source, interactive=0,parallel=self.parallel)

          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'stopcode'])
          else:
            ret=retpar 
          
          # check initial cyclethreshold
          # true value: 0.14600408
          report=self.th.checkall(ret=ret,imgexist=[self.img+'.psf', self.img+'.residual'], imgval=[(self.img+'.model',0.200,[256,256,0,0])], firstcyclethresh=0.14600408)
          #report=self.th.checkall(ret=ret,imgexist=[self.img+'.psf', self.img+'.residual'], imgval=[(self.img+'.model',0.200,[256,256,0,0])], firstcyclethresh=0.55600408)
          self.assertTrue(self.check_final(report))

     def test_iterbot_cube_restart_updatedmask(self): 
          """ [iterbot] Test_cube_restart_updatedmask : restart cube with an updated mask (CAS-13508 fix verification) """
          self.prepData('refim_twopoints_twochan.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=512,cell='8.0arcsec',specmode='cube',deconvolver='hogbom',niter=4,parallel=self.parallel)
          os.system('rm -rf '+self.img+'.mask')

          # Create a new mask around the 1Jy source only. 
          mask_around_dim_source='circle[[256pix,256pix],30pix]'
          makemask(inpimage=self.img+'.residual', inpmask=mask_around_dim_source, output=self.img+'.mask.dim.source',overwrite=True,mode='copy')
          mask_around_dim_source=self.img+'.mask.dim.source'

          # 2nd clean with the updated mask
          retpar=tclean(vis=self.msfile, imagename=self.img,imsize=512,cell='8.0arcsec',specmode='cube',niter=10,calcres=False, calcpsf=False, restart=True, mask=mask_around_dim_source, interactive=0,parallel=self.parallel)

          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'stopcode'])
          else:
            ret=retpar 
          
          # check initial cyclethreshold
          # true value: 
          # 0.68944335
          report=self.th.checkall(ret=ret,imgexist=[self.img+'.psf', self.img+'.residual'], imgval=[(self.img+'.model',0.83,[256,256,0,0]),(self.img+'.model', 0.27, [256,256,0,1])], firstcyclethresh=0.68944335)

          self.assertTrue(self.check_final(report))


     def test_iterbot_mfs_restart_pbmask(self): 
          """ [iterbot] Test_mfs_restart_updatedmask : restart mfs with pbmask (CAS-13508 fix verification) """
          self.prepData('refim_twopoints_twochan.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=512,cell='8.0arcsec',specmode='mfs',deconvolver='hogbom',niter=2,parallel=self.parallel)
          os.system('rm -rf '+self.img+'.mask')

          # 2nd clean with pbmask 
          retpar=tclean(vis=self.msfile, imagename=self.img,imsize=512,cell='8.0arcsec',specmode='mfs',niter=1,calcres=False, calcpsf=False, restart=True, usemask='pb', pbmask=0.8, interactive=0,parallel=self.parallel)

          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'stopcode'])
          else:
            ret=retpar 
          
          # check initial cyclethreshold
          # true value: 0.14600408
          report=self.th.checkall(ret=ret,imgexist=[self.img+'.psf', self.img+'.residual'], imgval=[(self.img+'.model',0.20,[256,256,0,0])], firstcyclethresh=0.1460498)

          self.assertTrue(self.check_final(report))

     def test_iterbot_cube_restart_pbmask(self): 
          """ [iterbot] Test_cube_restart_updatedmask : restart cube with pbmask (CAS-13508 fix verification) """
          self.prepData('refim_twopoints_twochan.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=512,cell='8.0arcsec',specmode='cube',deconvolver='hogbom',niter=4,parallel=self.parallel)
          os.system('rm -rf '+self.img+'.mask')

          # 2nd clean with pbmask 
          retpar=tclean(vis=self.msfile, imagename=self.img,imsize=512,cell='8.0arcsec',specmode='cube',niter=10,calcres=False, calcpsf=False, restart=True, usemask='pb', pbmask=0.8, interactive=0,parallel=self.parallel)

          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone', 'stopcode'])
          else:
            ret=retpar 
          
          # check initial cyclethreshold
          # true value:0.68944335  
          report=self.th.checkall(ret=ret,imgexist=[self.img+'.psf', self.img+'.residual'], imgval=[(self.img+'.model',0.37,[256,256,0,0]),(self.img+'.model', 0.15, [256,256,0,1])],firstcyclethresh=0.68944335)

          self.assertTrue(self.check_final(report))
##############################################
##############################################
##############################################

##Task level tests : multi-field, 2chan.
### For some of these tests, do the same with uvsub and compare ? 
class test_multifield(testref_base):
     
     def test_multifield_both_mfs(self):
          """ [multifield] Test_Multifield_both_mfs : Two fields, both mfs """
          self.prepData("refim_twopoints_twochan.ms")
          if is_CASA6 or not ParallelTaskHelper.isMPIEnabled():
               logstart = self.th.get_log_length()
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, 
                        iterdone=13, 
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image'],
                        imgval=[(self.img+'.image',1.075,[50,50,0,0]),
                               (self.img+'1.image',5.590,[40,40,0,0]),
                               (self.img+'.residual',0.04,[30,18,0,0])])

          if is_CASA6 or not ParallelTaskHelper.isMPIEnabled(): # the Casa 6 Casa5 Subtree Pull Request bamboo plan creates multiple log files -> this test doesn't work there
               report += self.th.check_logs(logstart, expected=[r"::deconvolve[\t ]+\[tst\] iters=", r"::deconvolve[\t ]+\[tst1\] iters="])

          self.assertTrue(self.check_final(report))

     def test_multifield_both_mtmfs(self):
          """ [multifield] Test_Multifield_both_mtmfs : Two fields, both mt-mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\n\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='mtmfs',interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, 
                        iterdone=12,
                        nmajordone=2,
                        imgexist=[self.img+'.image.tt0', self.img+'1.image.tt0',self.img+'.image.tt1', self.img+'1.image.tt1', self.img+'.alpha', self.img+'1.alpha'],
                        imgval=[(self.img+'.image.tt0',1.094,[50,50,0,0]),
                               (self.img+'1.image.tt0',5.577,[40,40,0,0]),
                               (self.img+'.alpha',-0.90,[50,50,0,0]),
                               (self.img+'1.alpha',-1.0,[40,40,0,0])])
          self.assertTrue(self.check_final(report))


     def test_multifield_both_cube(self):
          """ [multifield] Test_Multifield_both_cube : Two fields, both cube"""
          self.prepData("refim_twopoints_twochan.ms")
          if is_CASA6 or not ParallelTaskHelper.isMPIEnabled():
               logstart = self.th.get_log_length()
          #self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\n')
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nimagename='+self.img+'2\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:48.895 +40.55.58.543\n')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest',parallel=self.parallel)
          ret={}
          #if self.parallel:
            #ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone'])
            #if self.nnode < 2:
            #  iterdone_expected=42  # single server case = serial
            #else:
            #  iterdone_expected=46
          #else:
          iterdone_expected=42
          ret=retpar 
          report=self.th.checkall(ret=ret, 
                        #iterdone=42,
                        iterdone=iterdone_expected,
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image'],
                        imgval=[(self.img+'.image',1.434,[50,50,0,0]),
                               (self.img+'1.image',7.44,[40,40,0,0]),
                               (self.img+'2.image',7.42,[51,40,0,0]),
                               (self.img+'.image',0.758,[50,50,0,1]),
                               (self.img+'1.image',3.715,[40,40,0,1]),
                               (self.img+'2.image',3.675,[51,40,0,1]) ])

          prefix = r"::deconvolve::MPIServer-[0-9]+" if self.parallel else "::deconvolve"
          # the channel output is a questionmark, because the iterators could give each MPI process a single channel (which is also why we can't count on there being a ":C1")
          if is_CASA6 or not ParallelTaskHelper.isMPIEnabled(): # the Casa 6 Casa5 Subtree Pull Request bamboo plan creates multiple log files -> this test doesn't work there
               report += self.th.check_logs(logstart, expected=[r"Processing channels in range \[[0-9]+, [0-9]+\]",
                                                                prefix+r"[\t ]+\[tst(:C0)?\] iters=",# prefix+r"[\t ]+\[tst(:C1)?\] iters=",
                                                                prefix+r"[\t ]+\[tst1(:C0)?\] iters=",# prefix+r"[\t ]+\[tst1(:C1)?\] iters=",
                                                                prefix+r"[\t ]+\[tst2(:C0)?\] iters="])# prefix+r"[\t ]+\[tst2(:C1)?\] iters="])
          self.assertTrue(self.check_final(report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Diffirent nchans of cubes in multi-field imaging is  not supported in parallel mode")
     def test_multifield_both_cube_diffshape(self):
          """ [multifield] Test_Multifield_both_cube : Two fields, both cube but different nchans"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nnchan=3\n')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest',parallel=self.parallel)
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone'])
          else:
            ret=retpar 
          report=self.th.checkall(ret=ret, 
                        iterdone=22,
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image'],
                        imgval=[(self.img+'.image',1.434,[50,50,0,0]),
                               (self.img+'1.image',7.44,[40,40,0,0]),
                               (self.img+'.image',0.762,[50,50,0,1]),
                               (self.img+'1.image',3.70,[40,40,0,1]) ])
          self.assertTrue(self.check_final(report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Mixing cube and mfs in multi-field imaging  is not supported in parallel mode")
     def test_multifield_cube_mfs(self):
          """ [multifield] Test_Multifield_cube_mfs : Two fields, one cube and one mfs"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nspecmode=mfs\nnchan=1\n')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest',parallel=self.parallel)
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone'])
            iterdone_expected=22
          else:
            iterdone_expected=15
            ret=retpar 
          report=self.th.checkall(ret=ret, 
                        #iterdone=15,
                        iterdone=iterdone_expected,
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image'],
                        imgval=[(self.img+'.image',1.4,[50,50,0,0]),
                               (self.img+'1.image',5.6,[40,40,0,0]),
                               (self.img+'.image',0.75,[50,50,0,1])])
          self.assertTrue(self.check_final(report))

     def test_multifield_mfs_mtmfs(self):
          """ [multifield] Test_Multifield_mfs_mtmfs : Two fields, one mt-mfs and one mfs (i.e. different deconvolvers)"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nreffreq=1.5GHz\ndeconvolver=mtmfs\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, 
                        iterdone=13,
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image.tt0',self.img+'1.alpha'],
                        imgval=[(self.img+'.image',1.094,[50,50,0,0]),
                               (self.img+'1.image.tt0',5.57,[40,40,0,0]), 
                               (self.img+'1.alpha', -1.0, [40,40,0,0])  ])
          self.assertTrue(self.check_final(report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Mixing cube and mtmfs in multi-field imaging is not supported in parallel mode")
     def test_multifield_cube_mtmfs(self):
          """ [multifield] Test_Multifield_cube_mtmfs : Two fields, one cube and one mtmfs"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nreffreq=1.5GHz\ndeconvolver=mtmfs\nspecmode=mfs\n')
          retpar = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest',parallel=self.parallel)
          ret={}
          if self.parallel:
            ret=self.th.mergeParaCubeResults(retpar, ['iterdone', 'nmajordone'])
          else:
            ret=retpar 
          report=self.th.checkall(ret=ret, 
                        iterdone=15,  # two chans in one field, and one chan in the other
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image.tt0',self.img+'1.alpha'],
                        imgval=[(self.img+'.image',1.427,[50,50,0,0]),
                               (self.img+'1.image.tt0',5.575,[40,40,0,0]),
                               (self.img+'.image',0.762,[50,50,0,1]) , 
                               (self.img+'1.alpha', -1.0, [40,40,0,0])  ])
          self.assertTrue(self.check_final(report))


     def test_multifield_diff_gridders(self):
          """ [multifield] Test_Multifield_diff_gridders : Two fields, both mfs, gridft and wproject """
          self.prepData("refim_twopoints_twochan.ms")
#          ##Outlier uses gridft
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]\ngridder=gridft')
#          ## Outlier uses wproject but with different number of planes as the main field
#          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nmask=circle[[40pix,40pix],10pix]\ngridder=wproject\nwprojplanes=6')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',gridder='wproject',wprojplanes=4,interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, 
                        iterdone=13,
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image'],
                        imgval=[(self.img+'.image',1.075,[50,50,0,0]),
                               (self.img+'1.image',5.58,[40,40,0,0])])
          self.assertTrue(self.check_final(report))


     def test_multifield_autonames(self):
          """ [multifield] Test_Multifield_4 : Test auto increment of image names """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',outlierfile=self.img+'.out.txt',parallel=self.parallel)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',outlierfile=self.img+'.out.txt',restart=False,parallel=self.parallel)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',restart=False) # no outlier...
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',outlierfile=self.img+'.out.txt',restart=False,parallel=self.parallel)

          report=self.th.checkall(imgexist=[self.img+'.psf',self.img+'1.psf',self.img+'_2.psf',self.img+'1_2.psf',self.img+'_3.psf',self.img+'_4.psf',self.img+'1_4.psf'], imgexistnot=[self.img+'1_3.psf'] )
          self.assertTrue(self.check_final(report))


### TODO :  Either put a check so that if any fields overlap, an error is thrown. Or, do sensible model choosing for some modes but detect and complain for other modes where it's harder to pick which model image to use.
     def test_multifield_overlap_mfs(self):
          """ [multifield] Test_Multifield_overlap_mfs : Two overlapping image fields, both mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[200,200]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:59:02.426 +40.51.14.559')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:58:39.580 +40.55.55.931",outlierfile=self.img+'.out.txt',niter=20,deconvolver='hogbom',interactive=0,parallel=self.parallel)

          report=self.th.checkall(ret=ret, 
                        iterdone=40, ## both images see the brightest source. 
                        nmajordone=2,
                        imgexist=[self.img+'.image', self.img+'1.image'],
                        imgval=[(self.img+'.image',5.575,[48,51,0,0]),
                               (self.img+'1.image',5.574,[130,136,0,0])]) ## both images have correct flux (not twice or zero !)
          self.assertTrue(self.check_final(report))


     def test_multifield_overlap_mtmfs(self):
          """ [multifield] Test_Multifield_overlap_mtmfs : Two overlapping image fields, both mt-mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[200,200]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:59:02.426 +40.51.14.559\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:58:39.580 +40.55.55.931",outlierfile=self.img+'.out.txt',niter=20,deconvolver='mtmfs',interactive=0,parallel=self.parallel)
          report=self.th.checkall(ret=ret, 
                        iterdone=39, ## both images see the brightest source.
                        nmajordone=2,
                        imgexist=[self.img+'.image.tt0', self.img+'1.image.tt0'],
                        imgval=[(self.img+'.image.tt0',5.52,[48,51,0,0]),
                                (self.img+'1.image.tt0',5.53,[130,136,0,0]),
                               (self.img+'.alpha',-0.965,[48,51,0,0]),
                               (self.img+'1.alpha',-0.965,[130,136,0,0])]) 
          self.assertTrue(self.check_final(report))


     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Facetted mfs imaging test in parallel is skipped temporarily until a fix is found. ")
     def test_multifield_facets_mfs(self):
          """ [multifield] Test_Multifield_mfs_facets : Facetted imaging (mfs) """
          self.prepData("refim_twopoints_twochan.ms")
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',phasecenter="J2000 19:59:00.2 +40.50.15.50",facets=2,deconvolver='hogbom',niter=30,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf'],imgval=[(self.img+'.psf',1.0,[100,100,0,0]),(self.img+'.image',5.56,[127,143,0,0]) ] )
          self.assertTrue(self.check_final(report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Facetted mtmfs imaging test in parallel is skipped temporarily until a fix is found. ")
     def test_multifield_facets_mtmfs(self):
          """ [multifield] Test_facets_mtmfs : Facetted imaging (mt-mfs) """
          self.prepData("refim_twopoints_twochan.ms")
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',phasecenter="J2000 19:59:00.2 +40.50.15.50",facets=2,deconvolver='mtmfs',niter=30,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.alpha'],imgval=[(self.img+'.psf.tt0',1.0,[100,100,0,0]),(self.img+'.image.tt0',5.56,[127,143,0,0]),(self.img+'.alpha',-1.0,[127,143,0,0]) ] )
          self.assertTrue(self.check_final(report))


#     def test_multifield_cube_chunks(self):
#          """ [multifield] Test_Multifield_cube_chunks : Two fields, two sections of the same cube"""
#          self.prepData("refim_point.ms")
#          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nnchan=5\nstart=5')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='10.0arcsec',specmode='cube',nchan=5,start=0,outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,interpolation='nearest')
#          report=self.th.checkall(ret=ret, 
#                        iterdone=38,
#                        nmajordone=2,
#                        imgexist=[self.img+'.image', self.img+'1.image'],
#                        imgval=[(self.img+'.image',1.434,[50,50,0,0]),
#                               (self.img+'1.image',7.452,[40,40,0,0]),
#                               (self.img+'.image',0.762,[50,50,0,1]),
#                               (self.img+'1.image',3.702,[40,40,0,1]) ])
#          self.assertTrue(self.check_final(report))
##############################################
##############################################

##Task level tests : Stokes imaging options
class test_stokes(testref_base):

     def test_stokes_mfs_I(self):
          """ [stokes] Test_Stokes_I_mfs mfs with stokes I"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='I',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,0])])
          self.assertTrue(self.check_final(report))

     def test_stokes_mfs_IV(self):
          """ [stokes] Test_Stokes_mfs_IV : mfs with stokes IV"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IV',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,0]),(self.img+'.image',4.0,[50,50,1,0])  ])
          self.assertTrue(self.check_final(report))

     def test_stokes_mfs_QU(self):
          """ [stokes] Test_Stokes_mfs_QU : mfs with stokes QU"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='QU',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',2.0,[50,50,0,0]),(self.img+'.image',3.0,[50,50,1,0])  ])
          self.assertTrue(self.check_final(report))

     def test_stokes_mfs_Q(self):
          """ [stokes] Test_Stokes_mfs_Q : mfs with stokes Q"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='Q',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',2.0,[50,50,0,0]) ] )
          self.assertTrue(self.check_final(report))

     def test_stokes_mfs_U(self):
          """ [stokes] Test_Stokes_mfs_U : mfs with stokes U"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='U',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',3.0,[50,50,0,0]) ] )
          self.assertTrue(self.check_final(report))

     def test_stokes_mfs_V(self):
          """ [stokes] Test_Stokes_mfs_V : mfs with stokes V"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='V',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',4.0,[50,50,0,0]) ] )
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_I(self):
          """ [stokes] Test_Stokes_cube_I : cube with stokes I"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='I',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,0]),(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',1.0,[50,50,0,2]) ] )
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_IV(self):
          """ [stokes] Test_Stokes_stokes_IV : cube with stokes V"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IV',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,0]),(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',1.0,[50,50,0,2]),  (self.img+'.image',4.0,[50,50,1,0]),(self.img+'.image',4.0,[50,50,1,1]),(self.img+'.image',4.0,[50,50,1,2])] )
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_QU(self):
          """ [stokes] Test_Stokes_stokes_QU : cube with stokes QU"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='QU',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',2.0,[50,50,0,0]),(self.img+'.image',2.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,0,2]),  (self.img+'.image',3.0,[50,50,1,0]),(self.img+'.image',3.0,[50,50,1,1]),(self.img+'.image',3.0,[50,50,1,2])] )
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_Q(self):
          """ [stokes] Test_Stokes_cube_Q : cube with stokes Q"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='Q',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',2.0,[50,50,0,0]),(self.img+'.image',2.0,[50,50,0,1]) ,(self.img+'.image',2.0,[50,50,0,2]) ])
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_U(self):
          """ [stokes] Test_Stokes_cube_U : cube with stokes U"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='U',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',3.0,[50,50,0,0]),(self.img+'.image',3.0,[50,50,0,1]) ,(self.img+'.image',3.0,[50,50,0,2]) ])
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_V(self):
          """ [stokes] Test_Stokes_cube_V : cube with stokes V"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='V',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',4.0,[50,50,0,0]),(self.img+'.image',4.0,[50,50,0,1]) ,(self.img+'.image',4.0,[50,50,0,2]) ])
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_IQUV_fromRL(self):
          """ [stokes] Test_Stokes_cube_IQUV_fromRL : cube with stokes IQUV"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,3,1]) ])
          self.assertTrue(self.check_final(report))

     def test_stokes_cube_IQUV_fromXY(self):
          """ [stokes] Test_Stokes_cube_IQUV_fromXY : cube with stokes IQUV"""
          self.prepData('refim_point_linXY.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube',interpolation='nearest',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,3,1]) ])
          self.assertTrue(self.check_final(report))

     def test_stokes_mtmfs_Q(self):
          """ [stokes] Test_Stokes_mtmfs_Q : mtmfs with stokes Q"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='Q',deconvolver='mtmfs',nterms=2,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0'], imgexistnot=[self.img+'.image.alpha'], imgval=[(self.img+'.image.tt0',2.0,[50,50,0,0]) ] )
          self.assertTrue(self.check_final(report))

     def test_stokes_mtmfs_IQUV(self):
          """ [stokes] Test_Stokes_mtmfs_IQUV : mtmfs with stokes IQUV"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',deconvolver='mtmfs',nterms=2,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0'],imgexistnot=[self.img+'.image.alpha'], imgval=[(self.img+'.image.tt0',1.0,[50,50,0,0]),(self.img+'.image.tt0',2.0,[50,50,1,0]), (self.img+'.image.tt0',3.0,[50,50,2,0]),(self.img+'.image.tt0',4.0,[50,50,3,0]) ])
          _ia.open(self.img+'.image.tt0')
          if _ia.brightnessunit() == "Jy/beam":
               report = report + "(Pass : Units are Jy/beam in the restored image)\n"
          else:
               report = report + "(Fail : Units are not Jy/beam in the restored image)\n"
          _ia.close()
          self.assertTrue(self.check_final(report))


#     def test_stokes_cube_I_flags(self):
#          """ [onefield] Test_Stokes_cube_I_flags : cube with stokes I and only XY or YX flagged"""
#          self.prepData('refim_point_linXY.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube')
#          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,4,1]) ])

#     def test_stokes_cube_pseudo_I_flags(self):
#          """ [onefield] Test_Stokes_cube_pseudo_I_flags : cube with stokes I and one of XX or YY flagged"""
#          self.prepData('refim_point_linXY.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube')
#          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,4,1]) ])

##############################################
##############################################

##Task level tests : cube.
class test_cube(testref_base):

#     def __init__(self,methodName='runTest'):
#          testref_base.__init__(self,methodName)
#          self.test_cube_0.__func__.__doc__ %="aaaa"

     def setUp(self):
          super(test_cube, self).setUp()

          ## Setup some variables to use in all the tests

          ## chan 4 (TOPO)
          qfstart=_qa.quantity("1.2GHz")
          #qvstart=_qa.quantity("-59958.5km/s")
          # for restf=1.25GHz
          qvstart=_qa.quantity("11991.7km/s")
          # ch10
          #qvstart=_qa.quantity("16788.4km/s")

          #mfstart=_me.frequency('LSRK',_qa.quantity("1.09999GHz"))
          # ch4 (for rest 1.25GHz)
          mfstart=_me.frequency('LSRK',_qa.quantity("1.199989GHz"))
          mvstart=_me.radialvelocity('BARY',_qa.quantity("11977.6km/s"))
          #dop = _me.todoppler('radio',mfstart,_qa.quantity('1.0GHz'))
          mfstart10=_me.frequency('LSRK',_qa.quantity(" 1.17999GHz"))                                                        
          # doppler with ch4 freq
          dop = _me.todoppler('radio',mfstart,_qa.quantity('1.25GHz'))                                              

          #1chan width 
          #qvwidth = _qa.quantity("11991.700km/s")
          #qvwidth = _qa.quantity("4796.7km/s")
          qvwidth = _qa.quantity("11991.7km/s")
          mvwidth = _me.radialvelocity('TOPO',qvwidth)

          # restf = 1.25GHz
          # vel range: 59961.1 -  -31174.7 km/s (lsrk/radio)
          #            74952.3 -  -28238.3 km/s (lsrk/optical)  

          self.testList = {
                      0:{'imagename':'Cubetest_chandefstdefwidth','spw':'0','start':0,'width':1,'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, default start and width, LSRK'},
                      1:{'imagename':'Cubetest_chandefstdefwidthtopo','spw':'0','start':0,'width':1, 'outframe':'TOPO','veltype':'radio',
                        'desc':'channel, default start and width, TOPO'},
                      2:{'imagename':'Cubetest_chandefstwidth2','spw':'0','start':0,'width':2, 'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, default start, width=2, LSRK'},
                      3:{'imagename':'Cubetest_chanst5wd1','spw':'0','start':5,'width':1, 'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, start=5, default width, LSRK'},
                      # this will result in blank channnel images (calcChanFreqs requires start and width in channel       
                      # mode to be given in chan index                                                                 
                      4:{'imagename':'Cubetest_chandefstwd1spwsel','spw':'0:5~19','start':0,'width':1, 'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, spw=0:5~19, LSRK'},
                      #5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'40MHz','outframe':'TOPO',
                      #  'desc':'frequency, default start, width=\'40MHz\', TOPO'},
                      # data set changed!
                      5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'100MHz','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, default start, width=\'100MHz\'(2 x chanwidth), TOPO'},
                      6:{'imagename':'Cubetest_freqst2defwd','spw':'0','start':'1.1GHz','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, start=\'1.1GHz\', default width, TOPO'},
                      7:{'imagename':'Cubetest_freqst2defwdspwsel','spw':'0:4~19','start':'1.1GHz','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, start=\'1.1GHz\', default width, spw=0:4~19, TOPO'},
                      8:{'imagename':'Cubetest_freqst10wdm','spw':'0','start':'1.5GHz','width':'-50MHz','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, start=\'1.5GHz\', width=\'-50MHz\', TOPO'},
                      9:{'imagename':'Cubetest_veldefstwd2','spw':'0','start':'','width':'23983.4km/s','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, default start, width=\'23983.4km/s\', TOPO'},
                     10:{'imagename':'Cubetest_veldefstwd2m','spw':'0','start':'','width':'-23983.4km/s','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, default start, width=\'-23983.4m/s\', TOPO'},
                     11:{'imagename':'Cubetest_velst4defwd','spw':'0','start':'11991.7km/s','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, start=\'11991.7km/s\', default width, TOPO'},
                     12:{'imagename':'Cubetest_velst4defwdbary','spw':'0','start':'11977.6km/s','width':'','outframe':'BARY','veltype':'radio',
                        'desc':'velocity, start=\'11977.6km/s\', default width, BARY'},
                     # currently 13 is not quite properly working, investigating - 2014.08.27 TT 
                     # for refim_point.ms ch9=-41347.8km/s (opt)
                     #13:{'imagename':'Cubetest_optvelst10wdeflsrk','spw':'0','start':'-49962.6km/s','width':'',
                     13:{'imagename':'Cubetest_optvelst19wdlsrk','spw':'0','start':'-41347.8km/s','width':'20000km/s',
                        'veltype':'optical','outframe':'LSRK',
                     ##   'desc':'velocity, start=\'74952.3km/s\', default width, veltype=optical LSRK'},
                     #   'desc':'velocity, start=\'-49962.6km/s\', default width, veltype=optical LSRK'},
                        'desc':'velocity, start=\'-41347.5km/s\', default width , veltype=optical LSRK'},
                     14:{'imagename':'Cubetest_stqfreqdefwd','spw':'0','start':qfstart,'width':'', 'veltype':'radio','outframe':'',
                        'desc':'frequency, start(quanity)=%s, default width, veltype=radio LSRK' % qfstart},
                     15:{'imagename':'Cubetest_stmfreqdefwd','spw':'0','start':mfstart,'width':'', 'veltype':'radio','outframe':'',
                        'desc':'frequency, start=%s, default width, veltype=radio LSRK' % mfstart},
                     16:{'imagename':'Cubetest_stqveldefwd','spw':'0','start':qvstart,'width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity(quantity), start=%s, default width, TOPO ' % qvstart},
                     17:{'imagename':'Cubetest_stmveldefwd','spw':'0','start':mvstart,'width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity(measure), start=%s, default width(outframe=TOPO will be overridden)' % mvstart},
                     18:{'imagename':'Cubetest_veldefstqvwidth','spw':'0','start':'','width':qvwidth,'outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, default start, width(quantity)=%s' % qvwidth},
                     19:{'imagename':'Cubetest_veldefstmvwidth','spw':'0','start':'','width':mvwidth,'outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, default start, width(measure)=%s, TOPO' % mvwidth},
                     20:{'imagename':'Cubetest_stdopdefwd','spw':'0','start':dop,'width':'','outframe':'LSRK','veltype':'radio',
                        'desc':'doppler, start=%s, default width, LSRK' % dop},
                     # with a gap in spw channel sel
                     21:{'imagename':'Cubetest_st4gap','spw':'0:4~9;12~14','start':4,'width':'','outframe':'LSRK','veltype':'radio',
                        'desc':'channel, start=%s, default width, channel gap (10-11) LSRK' % 4},
                     # stride > 1
                     22:{'imagename':'Cubetest_st4stride2','spw':'0:0~10^2','start':0,'width':'','outframe':'LSRK','veltype':'radio', 'interpolation':'nearest',
                        'desc':'channel, start=%s, default width, step=2 LSRK nearest' % 0},
                     23:{'imagename':'Cubetest_defstspwchansel4','spw':'0:4~13','start':'','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'spw with channel selection( 0:4~13 ), default start, LSRK nearest'}
                    }
          
#          self.test_cube_0.__func__.__doc__ %=self.testList[0]['desc']
     

     def run_cubetclean(self, testid):
          """ core function to execute a cube tclean """
          if 'interpolation' in self.testList[testid]:
              interpolation = self.testList[testid]['interpolation']
          else:
              interpolation = 'linear'

          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw=self.testList[testid]['spw'],\
                       imagename=self.img+self.testList[testid]['imagename'],\
                       start=self.testList[testid]['start'],\
                       width=self.testList[testid]['width'],\
                       veltype=self.testList[testid]['veltype'],\
                       outframe=self.testList[testid]['outframe'], \
                       interpolation=interpolation,parallel=self.parallel)
          return ret

     def test_cube_0(self):
          """ [cube] Test_Cube_0 new """
          testid=0
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50002,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',999988750)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_1(self):
          """ [cube] Test_Cube_1  """
          testid=1
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50002,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO', 9.9999999e8)

          self.assertTrue(self.check_final(report+report2))

     def test_cube_2(self):
          """ [cube] Test_Cube_2  """
          testid=2
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.4643,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.02498846e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_3(self):
          """ [cube] Test_Cube_3  """
          # start = 5 (1.25GHZ IN TOPO)
          testid=3
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.2000,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.249985937e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_4(self):
          """ [cube] Test_Cube_4  """
          testid=4
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          # NEED CHECK!!!
          #report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          #imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.5000,
          #[50,50,0,0])])
          #report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.23998593e9)

     def test_cube_5(self):
          """ [cube] Test_Cube_5  """
          # width by freq (2x chanw) result should be the same as #2
          testid=5
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.4643,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.025e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_6(self):
          """ [cube] Test_Cube_6  """ 
          # start in freq=1.1GHz (=chan5)
          testid=6
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.36365354,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.1e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_7(self):
          """ [cube] Test_Cube_7  """
          # start 1.1GHz(TOPO)=chan2 spw=4~19
          # Currently different behaviors between serial and parallel on non-overlapping data and image
          # parameter selections.
          # serial: result in chan 0&1 psf blanked  
          # parallel: spw channel selection will be ignored and tuneselectdata will 
          # select overlapping data and image selections (this seems to me more correct? behavior)
          # as of 2019.01.08, this is no longer true, psf blanked for chan 0 and 1 for parallel
          testid=7
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          # parallel 
          ##if self.parallel:
          #    report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          #    imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.36,
          #    [50,50,0,0]),(self.img+self.testList[testid]['imagename']+'.image',1.2000,
          #    [50,50,0,3])])
          #else: # serial
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',0.0,
          [50,50,0,0]),(self.img+self.testList[testid]['imagename']+'.image',1.2000,
          [50,50,0,3])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.1e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_8(self):
          """ [cube] Test_Cube_8  """
          # start =1.5GHz(chan10)  width=-50MHz TOPO (descending freq)
          testid=8
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.42858946,
          [50,50,0,9])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.5e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_9(self):
          """ [cube] Test_Cube_9  """
          # width in vel (=23983.4km/s=2xChanW) def start (=cube will be ascending order in vel)
          testid=9
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.46184647,
          [50,50,0,9])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.925e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_10(self):
          """ [cube] Test_Cube_10  """
          # width in vel = -23983.4m/s def start (cube will be in descending order in vel)
          testid=10
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.46184647,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.025e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_11(self):
          """ [cube] Test_Cube_11  """
          # start 11991.7km/s (chan4)
          testid=11
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001776,
          [50,50,0,4])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.2e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_12(self):
          """ [cube] Test_Cube_12  """
          # start 11977.6km/s (BARY) = chan4
          testid=12
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
          [50,50,0,4])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','BARY',1.200058783e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_13(self):
          """ [cube] Test_Cube_13  """
          # 
          testid=13
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          # use own tclean command as nchan need to modify
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,specmode='cube',nchan=8,restfreq=['1.25GHz'],
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',spw=self.testList[testid]['spw'],
                       imagename=self.img+self.testList[testid]['imagename'],start=self.testList[testid]['start'], 
                       width=self.testList[testid]['width'],veltype=self.testList[testid]['veltype'],
                       outframe=self.testList[testid]['outframe'],parallel=self.parallel)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          #report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          #imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
          #[50,50,0,4])])
          #report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.2000e9)

     def test_cube_14(self):
          """ [cube] Test_Cube_14  """
          # start = quantity ('1.2GHz') frame default(LSRK)
          testid=14
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.25000215,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.2e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_15(self):
          """ [cube] Test_Cube_15  """
          # measure freq in LSRK ch4
          testid=15
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image', 1.25001216,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199989e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_16(self):
          """ [cube] Test_Cube_16  """
          # start quantity vel=11991.7km/s outframe=topo (ascending vel order)
          testid=16
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001776,
          [50,50,0,4])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.2000e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_17(self):
          """ [cube] Test_Cube_17  """
          # start measure vel=11977.6km/s BARY, outframe=TOPO will be overridedden (ascending vel order)
          testid=17
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
          [50,50,0,4])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','BARY',1.200058783e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_18(self):
          """ [cube] Test_Cube_18  """
          # defaut start, width in vel (quantity) +11991.7km/s (TOPO, radio)=datachan width, will be
          # ascending order in vel so highet DATA channel will be chan 0 in the image (image chan0=1.45GHz)
          testid=18
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001764,
          [50,50,0,9])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.45e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_19(self):
          """ [cube] Test_Cube_19  """
          # default start, width in vel (measure) +11991.7km/s (TOPO, radio)
          testid=19
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001764,
          [50,50,0,9])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.45e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_20(self):
          """ [cube] Test_Cube_20  """
          # doppler (with ch4 LSRK freq, rest freq=1.25GHz)
          testid=20
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.5000546,
          [50,50,0,4])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199989152e9)
          self.assertTrue(self.check_final(report+report2))

     # looks like it passes now....
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Data sel with channel gaps is not supported in parallel mode")
     def test_cube_21(self):
          """ [cube] Test_Cube_21  """
          # data sel with channel gap (10,11 excluded) 4~9, 12~14
          testid=21
          self.testList[testid]['interpolation']='nearest'
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.250001562, [50,50,0,0]),
                 (self.img+self.testList[testid]['imagename']+'.image',0.0, [50,50,0,6]),
                 (self.img+self.testList[testid]['imagename']+'.image',0.0, [50,50,0,7])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199986500e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_22(self):
          """ [cube] Test_Cube_22  """
          # stride (step=2) use nearest interpolation (other interpotion methods
          # may not work well...)
          testid=22
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.5000546,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',0.999988750387e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_23(self):
          """ [cube] Test_Cube_23  """
          testid=23
          print(" : " , self.testList[testid]['desc'])
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.2500156,
          [50,50,0,0])])
          report2 = self.th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.20e9)
          self.assertTrue(self.check_final(report+report2))

     # following tests for cube image spectral channel order for the data with decreasing channel frequecies
     def test_cube_descendF1(self):
           # default start and width
          # first image channel = first data channel, image channel frequecy descreases with increasing channel number
          """ [cube] Test_Cube_DescendF1: specmode cube with descending frequency data, default start and width  """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK',1.949978e9, -0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF2(self):
          # first image channel = data channel 5, image channel frequecy descreases with increasing channel number
          """ [cube] Test_Cube_DescendF2: specmode cube with descending frequency data, start in channel no. with default width  """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start=5, width='', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK',1.699981e9,-0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF3(self):
          # cube image should be identical with test_cube_descendF2
          # first image channel = data channel 5, width=1, ->image channel frequecy descreases with increasing channel number
          """ [cube] Test_Cube_DescendF3: specmode cube with descending frequency data, start in channel no. with  width=1  """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start=5, width=1, imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK',1.699981e9,-0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF4(self):
          # start in channel no., width=-1  
          # first image channel = data channel 9, -> channel frequecy increases with increasing image channel number  
          """ [cube] Test_Cube_DescendF4: specmode cube with descending frequency data, start in freuquency  with  default width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start=9, width=-1, imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK', 1.499983125e9,0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF5(self):
          # start in frequency, default width
          # first image channel = data channel 9, -> channel frequecy increases with increasing channel number
          """ [cube] Test_Cube_DescendF5: specmode cube with descending frequency data, start in freuquency  with  default width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='1.499983125GHz', width='', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK',1.499983125e9,0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF6(self):
          # start in frequency, positive width
          # first image channel = data channel 9, -> channel frequecy increases with increasing channel number
          """ [cube] Test_Cube_DescendF6: specmode cube with descending frequency data, start in freuquency  with  a positive width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='1.499983125GHz', width='0.049999438GHz', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK',1.499983125e9,0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF7(self):
          # start in frequency, negative width
          # first image channel = data channel 5, -> channel frequecy decreases with increasing channel number
          """ [cube] Test_Cube_DescendF7: specmode cube with descending frequency data, start in freuquency  with  a negative width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='1.699981GHz', width='-0.049999438GHz', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image','LSRK',1.699981e9,-0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF8(self):
          # start='',  a positive frequency width
          # the data channel 0 - 9  will be selected, but since width >0 for the ms, image channel order will be reversed w.r.t data channel order
          """ [cube] Test_Cube_DescendF8: specmode cube with descending frequency data, width in  a positive freuquency  with  default start """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='', width='0.049999438GHz', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 0.999989e9, 0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF9(self):
          # start='',  a netative frequency width
          # the data channel 0 - 9 will be selected, but since width <0 for the ms, image channel order will be the same order as data channel
          """ [cube] Test_Cube_DescendF9: specmode cube with descending frequency data, width in  a positive freuquency  with  default start """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='', width='-0.049999438GHz', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.44998369263e9, -0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF10(self):
          # start in velocity , width=''
          # first image channel = data channel 5, -> channel velocity increases  with increasing channel number
          """ [cube] Test_Cube_DescendF10: specmode cube with descendign frequency data, start in velocity with defualt width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='-107920.698km/s', width='', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.699980875e9, -0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF11(self):
          # start in velocity , a positive vel width
          # first image channel = data channel 5, -> channel velocity increases  with increasing channel number
          """ [cube] Test_Cube_DescendF11: specmode cube with descendign frequency data, start in velocity with defualt width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='-107920.698km/s', width='1.1991563418e4km/s', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.699980875e9, -0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF12(self):
          # start in velocity , a negative vel width
          # first image channel = data channel 9, -> channel velocity decreases with increasing channel number
          """ [cube] Test_Cube_DescendF12: specmode cube with descendign frequency data, start in velocity with defualt width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='-59954.444km/s', width='-1.1991563418e4km/s', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.49998312558e9, 0.049999438e9)
          self.assertTrue(self.check_final(report))

     def test_cube_descendF13(self):
          # width  in a positive velocity, default start
          # the data channel 10-19 (lower side in vel) will be selected, since vel width >0, image channel order will be reversed order w.r.t data channel order
          """ [cube] Test_Cube_DescendF13: specmode cube with descendign frequency data, start in velocity with defualt width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='', width='1.1991563418e4km/s', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.449983688e9, -0.0499994375194e9)
          self.assertTrue(self.check_final(report))


     def test_cube_descendF14(self):
          # width  in a negative velocity, default start
          # the data channel 10-19 (lower side in vel) will be selected, since vel width <0, image channel order will the same order w.r.t data channel order
          """ [cube] Test_Cube_DescendF14: specmode cube with descendign frequency data, start in velocity with defualt width """
          self.prepData('refim_point_descendingfreqs.ms')
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw='0', start='', width='-1.1991563418e4km/s', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )

          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 0.999988750387e9, 0.049999438e9)
          self.assertTrue(self.check_final(report))

     # now rans OK with CAS-9386
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Cubedata mode test in parallel is skipped temporarily until a fix is found. ")
     def test_cube_D1(self):
          """ [cube] Test_Cube_D1 : specmode cubedata - No runtime doppler corrections """
          self.prepData('refim_Cband.G37line.ms')
          ret = tclean(vis=self.msfile,field='1',spw='0:105~135',specmode='cubedata',nchan=30,start=105,width=1,veltype='radio',
                       imagename=self.img,imsize=256,cell='0.01arcmin',phasecenter=1,deconvolver='hogbom',niter=10,parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',86.254,[128,128,0,18])])
          ## line is smoother
          self.assertTrue(self.check_final(report))

     def test_cube_D2(self):
          """ [cube] Test_Cube_D2 : specmode cube - WITH doppler corrections """
          self.prepData('refim_Cband.G37line.ms')
          ret = tclean(vis=self.msfile,field='1',spw='0:105~135',specmode='cube',nchan=30,start=105,width=1,veltype='radio',
                       imagename=self.img,imsize=256,cell='0.01arcmin',phasecenter=1,deconvolver='hogbom',niter=10,parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',92.1789,[128,128,0,20])])
          ## line is tighter
          self.assertTrue(self.check_final(report))
     def test_cube_perchanweight_briggs(self):
          """[cube] test_cube_perchanweight_briggs: """
          self.prepData('refim_point_withline.ms')
          delmod(self.msfile)
          imnat=self.img+"_nat"
          imbriggs0=self.img+"_briggs0"
          imbriggs_2=self.img+"_briggs_2"
          imbriggs_3=self.img+"_briggs_3"
          retnat = tclean(vis=self.msfile,imagename=imnat,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='natural', parallel=self.parallel)
          ret0 = tclean(vis=self.msfile,imagename=imbriggs0,imsize=100,cell='8.0arcsec',specmode='cube', perchanweightdensity=True,deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='briggs', robust=0, parallel=self.parallel)
          ret_2=tclean(vis=self.msfile,imagename=imbriggs_2,imsize=100,cell='8.0arcsec',specmode='cube', perchanweightdensity=True,deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='briggs', robust=-2.0, parallel=self.parallel)
          ret_3=tclean(vis=self.msfile,imagename=imbriggs_3,imsize=100,cell='8.0arcsec',specmode='cube', perchanweightdensity=True,deconvolver='hogbom',niter=1,threshold='0Jy',interactive=0, weighting='briggs', robust=0.0, uvtaper=['50arcsec'], parallel=self.parallel)
          self.assertTrue(os.path.exists(imnat+'.image') and os.path.exists(imbriggs0+'.image') and os.path.exists(imbriggs_2+'.image') and  os.path.exists(imbriggs_3+'.image') )
          self.assertTrue(self.th.check_beam_compare(imbriggs0+'.image', imnat+'.image', operator.lt))
          self.assertTrue(self.th.check_beam_compare(imbriggs_2+'.image', imbriggs0+'.image', operator.lt))
          self.assertTrue(self.th.check_beam_compare(imbriggs0+'.image', imbriggs_3+'.image', operator.lt))
#     def test_cube_D3(self):
#          """ EMPTY : [cube] Test_Cube_D3 : specmode cubesrc - Doppler correct to a SOURCE ephemeris"""
#          ret = tclean(vis=self.msfile,field='1',spw='0:105~135',specmode='cubesrc',nchan=30,start=105,width=1,veltype='radio',imagename=self.img,imsize=256,cell='0.01arcmin',phasecenter=1,deconvolver='hogbom',niter=10)
#          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )

     @unittest.skipIf(is_CASA6, "Skip because plotms is not available in casatasks")
     def test_cube_continuum_subtract_uvsub(self):
          """ [cube] Test_Cube_continuum_subtract :  Using uvsub """
          self.prepData('refim_point_withline.ms')
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='reset0')
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='data',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step0data.png',title="original data")
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step0model.png',title="empty model")

          # Let's include a subdir in the output image name. This could cause failures, at
          # least in parallel mode (CAS-10937).
          imagename = os.path.join(self.img_subdir, self.img)
          print("IMAGENAME=",imagename)
          ret = tclean(vis=self.msfile,imagename=imagename,imsize=100,cell='8.0arcsec', spw='0:12~19',niter=50,gain=0.2,savemodel='modelcolumn',
                       deconvolver='mtmfs',parallel=self.parallel)
#          self.assertTrue(self.th.exists(self.img+'.model') )
#          self.assertTrue( self.th.check_modelchan(self.msfile,10) == 0.0 and self.th.check_modelchan(self.msfile,3) > 0.0 )
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step1.png',title="model after partial mtmfs on some channels")

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='reset0')

#          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=[imagename+'.model.tt0',imagename+'.model.tt1'], 
#                       spw='0',niter=0,savemodel='modelcolumn',deconvolver='mtmfs',parallel=self.parallel)
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=[imagename+'.model.tt0',imagename+'.model.tt1'], 
                       spw='0',niter=0,savemodel='modelcolumn',deconvolver='mtmfs',parallel=self.parallel)

#          self.assertTrue( self.th.check_modelchan(self.msfile,10) > 0.0 and self.th.check_modelchan(self.msfile,3) > 0.0 
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step2.png',title="model after mtmfs predict on full spw" )

          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='data-model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step3data.png',title="data-model")
          

#     def test_cube_continuum_subtract_otf(self):
#          """ EMPTY : [cube] Test_Cube_continuum_subtract :  On-The-Fly using multifield """
#          self.prepData('refim_point_withline.ms')

     def test_cube_badchannel_restoringbeam(self):
          """ [cube] Test auto restoring beam with a bad edge channel """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
          report1=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',0.889,[54,50,0,0]) , (self.img+'.image',0.0602,[54,50,0,19]) , (self.img+'.residual',0.033942,[54,50,0,19]) ])
          # first channel's psf is 'bad' and wider along one axis. This offcenter location is higher in value

          #  For parallel mode, to get common beam, need to run anoter tclean run with serial
          if self.parallel:
              ret = tclean(vis=self.msfile,imagename=self.img+'1',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',
                       restoration=False, parallel=self.parallel)
              ret2 = tclean(vis=self.msfile,imagename=self.img+'1',specmode='cube',imsize=100,cell='10.0arcsec',niter=0,deconvolver='hogbom',
                       restoration=True, restoringbeam='common', calcres=False, calcpsf=False, parallel=False)
          else:    
              ret = tclean(vis=self.msfile,imagename=self.img+'1',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',
                       restoringbeam='common',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'1.psf') and os.path.exists(self.img+'1.image') )
          report2=self.th.checkall(imgexist=[self.img+'1.image'],imgval=[(self.img+'1.image',0.8906,[54,50,0,0]), (self.img+'1.image',0.51977,[54,50,0,19]) , (self.img+'1.residual',0.033942,[54,50,0,19]) ])
          # OLD - first channel has been restored by a 'common' beam picked from channel 2
          self.assertTrue(self.check_final(report1+report2))

#  def test_cube_explicit_restoringbeam(self):
#          """ [cube] Test explicit restoring beams : Test peak flux and off source value for smoothed residuals"""



     def test_cube_common_restoringbeam(self):
          """ [cube] Test_cube_restoringbeam (cas10849/10946) : Test parallel and serial run on same refconcat images  """
          
          self.prepData('refim_point.ms')
          
          # Imaging run - no restoration (serial or parallel)
          ret = tclean(vis=self.msfile,imagename=self.img,
                       imsize=100,cell='10.0arcsec',
                       interactive=0,niter=10,specmode='cube',
                       restoration=False, parallel=self.parallel)
          
          # Serial restart for restoration only (serial only)
          retpar = tclean(vis=self.msfile,imagename=self.img,
                          imsize=100,cell='10.0arcsec',
                          interactive=0,niter=0,specmode='cube',
                          restoration=True, restoringbeam='common',parallel=False, #### always False. 
                          calcres=False, calcpsf=False)
          
          header = imhead(self.img+'.image',verbose=False)
               
          estr = "["+inspect.stack()[1][3]+"] Has single restoring beam ? : " + self.th.verdict( 'restoringbeam' in header) + "\n"

          report2 = self.th.checkall(imgexist=[self.img+'.image'], 
                                     imgval=[(self.img+'.image',0.770450,[54,50,0,1]),
                                            (self.img+'.image',0.567246,[54,50,0,15])  ])
          
          ## Pass or Fail (and why) ?
          self.assertTrue(self.check_final(estr+report2))

     # Now works for parallel  with CAS-9396 
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily for 5.5")
#     def test_cube_chanchunks(self):
#          """ [cube] Test channel chunking for large cubes """
#          self.prepData('refim_point.ms')
##          ret = tclean(vis=self.msfile,imagename=self.img,specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom', savemodel='modelcolumn')
##          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
##          report=self.th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.5002,[50,50,0,0]) , (self.img+'.image',0.769,[50,50,0,19]) ])
#
#          ret = tclean(vis=self.msfile,imagename=self.img+'cc',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',parallel=self.parallel)
#          self.assertTrue(os.path.exists(self.img+'cc.psf') and os.path.exists(self.img+'cc.image') )
#          report=self.th.checkall(imgexist=[self.img+'cc.image'],imgval=[(self.img+'cc.image',1.5002,[50,50,0,0]) , (self.img+'cc.image',0.769,[50,50,0,19]) ])
#          self.assertTrue(self.check_final(report))

     def test_cube_chanchunks_auto(self):
          """ [cube] Test channel chunking for large cubes : automatic calc of nchanchunks """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img+'cc',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'cc.psf') and os.path.exists(self.img+'cc.image') )
          report=self.th.checkall(imgexist=[self.img+'cc.image'],imgval=[(self.img+'cc.image',1.5002,[50,50,0,0]) , (self.img+'cc.image',0.769,[50,50,0,19]) ])
          self.assertTrue(self.check_final(report))


     # Now passes for parallel
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_cube_chanchunks_savemodel(self):
          """ [cube] Test channel chunking for large cubes and save model """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img+'cc',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom', savemodel='modelcolumn',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'cc.psf') and os.path.exists(self.img+'cc.image') )
          report=self.th.checkall(imgexist=[self.img+'cc.image'],imgval=[(self.img+'cc.image',1.5002,[50,50,0,0]) , (self.img+'cc.image',0.769,[50,50,0,19]) ])
          self.assertTrue( self.th.check_modelchan(self.msfile,5) > 0.0 and self.th.check_modelchan(self.msfile,18) > 0.0 )
          self.assertTrue(self.check_final(report))
          
#     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     @unittest.skip('Skip. This test deprecated. no longer need mtmfs for cube use msclean')
     def test_cube_mtmfs_nterms1(self):		
          """ [cube] Test mtmfs with cube and nterms = 1 """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img+'cc', specmode='cube', imsize=200,cell='8.0arcsec',niter=10,deconvolver='mtmfs',nterms=1,interactive=0,parallel=self.parallel,scales=[0,20,40,100])
          report=self.th.checkall(ret=ret, imgexist=[self.img+'cc.psf.tt0', self.img+'cc.residual.tt0', self.img+'cc.image.tt0', self.img+'cc.model.tt0'],imgval=[(self.img+'cc.image.tt0',1.0,[100,100,0,0]),(self.img+'cc.image.tt0',0.492,[100,100,0,1]),(self.img+'cc.image.tt0',0.281,[100,100,0,2])])		
          self.assertTrue(self.check_final(report))
          
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     @unittest.skip('Skip. This test deprecated. no longer need mtmfs for cube use msclean')
     def test_cubedata_mtmfs_nterms1(self):
          """ [cube] Test mtmfs with cube data and nterms = 1 """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img+'cc', specmode='cubedata', imsize=200,cell='8.0arcsec',niter=10,deconvolver='mtmfs',nterms=1,interactive=0,parallel=self.parallel,scales=[0,20,40,100])
          report=self.th.checkall(ret=ret, imgexist=[self.img+'cc.psf.tt0', self.img+'cc.residual.tt0', self.img+'cc.image.tt0', self.img+'cc.model.tt0'],imgval=[(self.img+'cc.image.tt0',1.0,[100,100,0,0]),(self.img+'cc.image.tt0',0.492,[100,100,0,1]),(self.img+'cc.image.tt0',0.281,[100,100,0,2])])		
          self.assertTrue(self.check_final(report)) 

     def test_cube_twoMS_startfreq(self):
          """ [cube] Test cube with list of two MSs with start in frequency specified (test CAS-12877 fix) """
          # The two tclean runs should produce identical cubes with the same image spectral coordinates
          ms1 = 'refim_point_first11chans.ms'
          ms2 = 'refim_point_last10chans.ms'
          self.prepData(ms1)
          self.prepData(ms2)
          # start f is at chan4 (1.2GHz in TOPO) 
          ret = tclean(vis=[ms1, ms2],field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       deconvolver='hogbom',\
                       spw=['0','0'], start='1.199989GHz', imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.199989e9)

          ret2 = tclean(vis=[ms2, ms1],field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       deconvolver='hogbom',\
                       spw=['0','0'], start='1.199989GHz', imagename=self.img+'_reverse',veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'_reverse.psf') and os.path.exists(self.img+'_reverse.image') )
          report2 = self.th.check_spec_frame(self.img+'_reverse.image', 'LSRK', 1.199989e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_twoMS_startvel(self):
          """ [cube] Test cube with list of two MSs with start in velocity specified (test CAS-12877 fix) """
          # The two tclean runs should produce identical cubes with the same image spectral coordinates
          ms1 = 'refim_point_first11chans.ms'
          ms2 = 'refim_point_last10chans.ms'
          self.prepData(ms1)
          self.prepData(ms2)
          ret = tclean(vis=[ms1, ms2],field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       deconvolver='hogbom',\
                       spw=['0','0'], start='11994.3km/s', width='-11991.7km/s',imagename=self.img,veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image'))
          report = self.th.check_spec_frame(self.img+'.image', 'LSRK', 1.199989e9)

          ret2 = tclean(vis=[ms2, ms1],field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       deconvolver='hogbom',\
                       spw=['0','0'], start='11994.3km/s',width='-11991.7km/s',  imagename=self.img+'_reverse',veltype='radio',outframe='LSRK',parallel=self.parallel)
          self.assertTrue(os.path.exists(self.img+'_reverse.psf') and os.path.exists(self.img+'_reverse.image'))
          report2 = self.th.check_spec_frame(self.img+'_reverse.image', 'LSRK', 1.199989e9)
          self.assertTrue(self.check_final(report+report2))

     def test_cube_flagged_mosaic_hogbom(self):
          """CAS-12957: 0-value channels aren't skipped with gridder=mosaic and initial channels are flagged"""
          # These tests are mainly here as regression test. The bug related to CAS-12957 was only known to affect multiscale clean, and here we test for similar bugs in hogbom.
          self.prepData('refim_twochan.ms')
          if is_CASA6 or not ParallelTaskHelper.isMPIEnabled():
               logstart = self.th.get_log_length()
          flagdata(self.msfile, spw='*:0')
          ret = tclean(self.msfile, imagename=self.img, specmode='cube', imsize=20, cell='8.0arcsec', scales=[0,5,10], niter=10, cycleniter=10, threshold=0, nchan=2, spw='0', interactive=0, \
                       deconvolver='hogbom', gridder='mosaic')
          report=self.th.checkall(imgexist=[self.img+'.model'], imgval=[(self.img+'.model', 0.01324, [10,10,0,1])], \
                                  imgvalexact=[(self.img+'.model', 0, [1,1,0,0]), (self.img+'.model', 0, [10,10,0,0])]   , tfmask=[(self.img+'.image',['mask0'])])#, epsilon=0.2)

          prefix = r"::deconvolve::MPIServer-[0-9]+" if self.parallel else "::deconvolve"
          # the channel output is a questionmark, because the iterators could give each MPI process a single channel (which is also why we can't count on there being a ":C1")
          if is_CASA6 or not ParallelTaskHelper.isMPIEnabled(): # the Casa 6 Casa5 Subtree Pull Request bamboo plan creates multiple log files -> this test doesn't work there
               report += self.th.check_logs(logstart, expected=[r"Processing channels in range \[[0-9]+, [0-9]+\]", prefix+r"[\t ]+\[tst(:C0)?\] iters="])#, prefix+r"[\t ]+\[tst(:C1)?\] iters="])

          self.assertTrue(self.check_final(pstr=report))

     def test_cube_flagged_mosaic_clark(self):
          """CAS-12957: 0-value channels aren't skipped with gridder=mosaic and initial channels are flagged"""
          # These tests are mainly here as regression test. The bug related to CAS-12957 was only known to affect multiscale clean, and here we test for similar bugs in clark.
          self.prepData('refim_twochan.ms')
          flagdata(self.msfile, spw='*:0')
          ret = tclean(self.msfile, imagename=self.img, specmode='cube', imsize=20, cell='8.0arcsec', scales=[0,5,10], niter=10, cycleniter=10, threshold=0, nchan=2, spw='0', interactive=0, \
                       deconvolver='clark', gridder='mosaic')
          report=self.th.checkall(imgexist=[self.img+'.model'], imgval=[(self.img+'.model', 0.01252, [10,10,0,1])], \
                                  imgvalexact=[(self.img+'.model', 0, [1,1,0,0]), (self.img+'.model', 0, [10,10,0,0])])#, epsilon=0.2)
          self.assertTrue(self.check_final(pstr=report))

     def test_cube_flagged_mosaic_multiscale(self):
          """CAS-12957: 0-value channels aren't skipped with gridder=mosaic and initial channels are flagged"""
          self.prepData('refim_twochan.ms')
          flagdata(self.msfile, spw='*:0')
          ret = tclean(self.msfile, imagename=self.img, specmode='cube', imsize=20, cell='8.0arcsec', scales=[0,5,10], niter=10, cycleniter=10, threshold=0, nchan=2, spw='0', interactive=0, \
                       deconvolver='multiscale', gridder='mosaic')
          report=self.th.checkall(imgexist=[self.img+'.model'], imgval=[(self.img+'.model', 0.01086, [10,10,0,1])], \
                                  imgvalexact=[(self.img+'.model', 0, [1,1,0,0]), (self.img+'.model', 0, [10,10,0,0])])#, epsilon=0.2)
          self.assertTrue(self.check_final(pstr=report))

     @unittest.skip('Skip test.')
     def test_cube_flagged_mosaic_mtmfs(self):
          """CAS-12957: 0-value channels aren't skipped with gridder=mosaic and initial channels are flagged"""
          # These tests are mainly here as regression test. The bug related to CAS-12957 was only known to affect multiscale clean, and here we test for similar bugs in mtmfs.
          self.prepData('refim_twochan.ms')
          flagdata(self.msfile, spw='*:0')
          ret = tclean(self.msfile, imagename=self.img, imsize=20, cell='8.0arcsec', scales=[0,5,10], niter=10, cycleniter=10, threshold=0, nchan=2, spw='0', interactive=0, \
                       deconvolver='mtmfs', nterms=1, gridder='mosaic')
          report=self.th.checkall(imgexist=[self.img+'.model.tt0'], imgval=[(self.img+'.model.tt0', 0.00530, [10,10,0,1])], \
                                  imgvalexact=[(self.img+'.model.tt0', 0, [1,1,0,0]), (self.img+'.model.tt0', 0, [10,10,0,0])])#, epsilon=0.2)
          self.assertTrue(self.check_final(pstr=report))

##############################################
##############################################

##Task level tests : masks and clean boxes.
class test_mask(testref_base):

     def test_mask_1(self):
          """ [mask] test_mask_1 : Input mask as file and string : mfs """
          self.prepData('refim_twochan.ms')
          mstr = 'circle[[50pix,80pix],10pix]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',
                        mask=self.img+'.mask.txt',parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=mstr,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.mask', self.img+'2.mask'], imgval=[(self.img+'1.mask',0.0,[50,50,0,0]),(self.img+'1.mask',1.0,[50,80,0,0]),(self.img+'2.mask',0.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,80,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_2(self):
          """ [mask] test_mask_2 :  Input mask as file and string : cube (few channels) """
          self.prepData('refim_point.ms')
          mstr =  'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',
                        interactive=0,usemask='user',mask=self.img+'.mask.txt',parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',
                        interactive=0,usemask='user',mask=mstr,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.mask', self.img+'2.mask'], imgval=[(self.img+'1.mask',0.0,[50,50,0,1]),(self.img+'1.mask',1.0,[50,50,0,2]),(self.img+'1.mask',1.0,[50,50,0,10]),(self.img+'1.mask',0.0,[50,50,0,11]),(self.img+'2.mask',0.0,[50,50,0,1]),(self.img+'2.mask',1.0,[50,50,0,2]),(self.img+'2.mask',1.0,[50,50,0,10]),(self.img+'2.mask',0.0,[50,50,0,11])])
          self.assertTrue(self.check_final(report))

     def test_mask_3(self):
          """ [mask] test_mask_3 : Input mask as image-to-be-regridded (ra/dec) : mfs """
          self.prepData('refim_twochan.ms')
          mstr = 'circle[[50pix,50pix],10pix]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,
                        usemask='user',mask=self.img+'.mask.txt',parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,
                        usemask='user',mask=self.img+'1.mask',phasecenter='J2000 19h59m57.5s +40d49m00.077s',parallel=self.parallel) # shift phasecenter
          report=self.th.checkall(imgexist=[self.img+'1.mask', self.img+'2.mask'], imgval=[(self.img+'1.mask',1.0,[50,50,0,0]),(self.img+'2.mask',1.0,[91,13,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_4(self):
          """ [mask] test_mask_4 :  Input mask as image-to-be-regridded(ra/dec/specframe) : cube """
          self.prepData('refim_point.ms')
          mstr =  'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',
                        interactive=0,usemask='user',mask=self.img+'.mask.txt',parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',
                        start='1.3GHz',interactive=0,usemask='user',mask=self.img+'1.mask',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.mask', self.img+'2.mask'], imgval=[(self.img+'1.mask',0.0,[50,50,0,1]),(self.img+'1.mask',1.0,[50,50,0,2]),(self.img+'1.mask',1.0,[50,50,0,10]),(self.img+'1.mask',0.0,[50,50,0,11]),(self.img+'2.mask',1.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,50,0,4]),(self.img+'2.mask',0.0,[50,50,0,10])])
          self.assertTrue(self.check_final(report))

     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily for 5.5")
     # parallel mode issue was fixed in imageanalysis 2019.05.23
     def test_mask_5(self):
          """ [mask] test_mask_5 : Input cube mask that has different chan
          ranges (use mask from the 1st tclean with a different channel range in the 2nd tclean run)"""
          self.prepData('refim_point.ms')
          mstr = 'circle[[50pix,50pix],10pix]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=1,deconvolver='hogbom',specmode='cube',
                        start=0,nchan=10,interactive=0,usemask='user',mask=self.img+'.mask.txt',parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=1,deconvolver='hogbom',specmode='cube',
                        start=5,nchan=10,interactive=0,usemask='user',mask=self.img+'1.mask',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.mask', self.img+'2.mask'], imgval=[(self.img+'1.mask',1.0,[50,50,0,1]),(self.img+'1.mask',1.0,[50,50,0,2]),(self.img+'1.mask',1.0,[50,50,0,9]),(self.img+'2.mask',1.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,50,0,4]),(self.img+'2.mask',0.0,[50,50,0,5])])
          self.assertTrue(self.check_final(report))

# the option, auto-thresh removed
#     def test_mask_autobox(self):
#         # changed to use threshold based automasking 
#          """ [mask] test_mask_autobox :  Autobox """
#          self.prepData('refim_twochan.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh')
#          # temporarily change value test to make it pass until extra masking in final minor cycle is resolved....
#          #report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,80,0,0])])
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
#          self.assertTrue(self.check_final(report))

# the option, auto-thresh removed
#     def test_mask_autobox_redraw(self):
#         # changed to use threshold based automasking 
#          """ [mask] test_mask_autobox_redraw :  Autoboxing with a redraw after each major cycle """
#          self.prepData('refim_eptwochan.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh',maskthreshold=0.5)
#          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=20,cycleniter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh',maskthreshold=0.5)
#          # tweak in automask threshold in the code changed masking extent 2016-03-21
#          #report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[60,30,0,0]),(self.img+'2.mask',1.0,[60,30,0,0])])
          # temporarily change the value test for unmasked region to make it pass (replace with the above when the extra masking issue is resolved...)
#          #report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[60,85,0,0]),(self.img+'2.mask',1.0,[60,30,0,0])])
#          #change in behavior due to automask code modification on July 1st,2016
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[60,85,0,0]),(self.img+'2.mask',0.0,[60,30,0,0])])
#          self.assertTrue(self.check_final(report))

# the option, auto-thresh removed
#     def test_mask_autobox_nmask(self):
#          """ [mask] test_mask_autobox_nmask : Autoboxing with nmask """
#          # this won't be triggering actual pruning but just to check going into write places
#          self.prepData('refim_point.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',
#                       interactive=0,usemask='auto-thresh',nmask=3)
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
#          self.assertTrue(self.check_final(report)
           
# the option, auto-thresh2 removed
#     def test_mask_autobox2_nmask(self):
#          """ [mask] test_mask_autobox2_nmask : Autoboxing (no binning) with nmask"""
#          # this won't be triggering actual pruning but just to check going into write places
#          self.prepData('refim_point.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',
#                       interactive=0,usemask='auto-thresh2',nmask=3)
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
#          self.assertTrue(self.check_final(report))

# the option, auto-thresh removed
#     def test_mask_autobox_pbmask(self):
#         """ [mask] test_mask_autobox_nmask : Autoboxing  with pbmask"""
#          # this won't be triggering actual pruning but just to check going into write places
#          self.prepData('refim_point.ms')
#          # change imsize to see the pbmask boundary
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=500,cell='8.0arcsec',niter=10,deconvolver='hogbom',
#                       interactive=0,usemask='auto-thresh', pbmask=0.2)
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[250,250,0,0]),(self.img+'.mask',0.0,[250,285,0,0]),(self.img+'.mask',0.0,[360,360])])
#          self.assertTrue(self.check_final(report))

#     @unittest.skip('Skip. This test deprecated. removed autoadjust param.')
#     def test_mask_autobox_autoadjust(self):
#          """ [mask] test_mask_autobox_autoadjust : Autoboxing with autoadjust=T """
#          self.prepData('refim_point.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',
#                       interactive=0,usemask='auto-thresh',autoadjust=True)
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
#          self.assertTrue(self.check_final(report)
     @unittest.skip('Skip test.')
     def test_mask_pbmask(self):
          """ [mask] test_mask_pbmask :  pb mask """
          pass

     @unittest.skip('Skip test.')
     def test_mask_combined_1(self):
          """ [mask] test_mask_combined_1 :  string + pbmask """
          pass

     @unittest.skip('Skip test.')
     def test_mask_combined_2(self):
          """ [mask] test_mask_combined_2 :  Autobox + pbmask """
          pass

     @unittest.skip('Skip test.')
     def test_mask_outlier(self):
          """ [mask] test_mask_outlier : With outlier fields """
          pass

# the option, auto-thresh removed
#     def test_mask_restart(self):
#          """ [mask] test_mask_restart : Test that mask reloads upon restart """
#          self.prepData('refim_twochan.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh')
#          ret2 = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          #report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,80,0,0])])
          # temporarily change the value test for unmasked region to make it pass (replace with the above when the extra masking issue is resolved...)
#          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
#          self.assertTrue(self.check_final(report))

     # AUTOMASK TESTS
     def test_mask_autobox_multithresh(self):
          """ [mask] test_mask__autobox_multithresh :  multi-threshold Autobox (default)"""
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_newnoise(self):
          """ [mask] test_mask__autobox_multithresh_newnoise :  multi-threshold Autobox invoking the new noise calc."""
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh', fastnoise=False)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_with_nsigma(self):
          """ [mask] test_mask__autobox_multithresh :  multi-threshold Autobox invoking the new noise calc."""
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh', nsigma=3.0)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_with_nsigma_newnoise(self):
          """ [mask] test_mask__autobox_multithresh :  multi-threshold Autobox invoking the new noise calc."""
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh', nsigma=3.0, fastnoise=False)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_with_prune(self):
          """ [mask] test_mask__autobox_multithresh_with_prune :  multi-threshold Autobox (minbeamfrac=0.3)"""
          # also test for a bug fix to the new pruneRegions (only caused the failure when image size large
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=1000,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh',
          minbeamfrac=0.3,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[500,500,0,0]),(self.img+'.mask',0.0,[500,510,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_with_stopmask(self):
          """ [mask] test_mask__autobox_multithresh_with_stopmask :  multi-threshold Autobox (minbeamfrac=0.3) with stop mask on """
          # will trigger stop mask condition for the last cycle (Cycle 4) - does not change output mask but can be checked on the log 
          #  
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=100,deconvolver='hogbom',interactive=0,
           usemask='auto-multithresh', minbeamfrac=0.3, minpercentchange=0.2,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[63,50,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_with_absorption(self):
          """ [mask] test_mask__autobox_multithresh_on_absorption :  multi-threshold Autobox (minbeamfrac=0.3) on the data with both emission and absorption  """
          # data with a emission pt and absorption pt.
          self.prepData('refim_point_pos_neg.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=100,deconvolver='hogbom',interactive=0,
                       usemask='auto-multithresh', negativethreshold=5.0,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[60,40,0,0]),(self.img+'.mask',0.0,[65,50,0,0])])
          self.assertTrue(self.check_final(report))

     def test_mask_autobox_multithresh_mfs_IQUV(self):
          """ [mask] test_mask__autobox_multithresh_mtmfs_IQUV :  multi-threshold Autobox (minbeamfrac=0.3) with cube full polarizaiton (IQUV) imaging """
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='mfs',interpolation='nearest',usemask="auto-multithresh", verbose=True, parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[40,60,0,0]),(self.img+'.mask',0.0,[65,50,0,0]), (self.img+'.mask', 1.0,[40,60,3,0])])
          self.assertTrue(self.check_final(report))
     
     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "MPI and serial results for the .mask images vary widely. To be fixed in CAS-13582.")
     def test_mask_autobox_multithresh_cube_IQUV(self):
          """ [mask] test_mask__autobox_multithresh_cube_IQUV :  multi-threshold Autobox (minbeamfrac=0.05) with cube full polarizaiton (IQUV) imaging """
          self.prepData('refim_point_linXY.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube',interpolation='nearest',usemask="auto-multithresh", minbeamfrac=0.05,  verbose=True, parallel=self.parallel)
          # test values that pass for serial runs:
          report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[35,75,0,0]),(self.img+'.mask',0.0,[32,80,1,1]), (self.img+'.mask',1.0,[35,60,1,1]), (self.img+'.mask',1.0,[60,30,3,0])])
          # test values that pass for mpi runs:
          #report=self.th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[35,75,0,0]),(self.img+'.mask',0.0,[36,74,1,1]), (self.img+'.mask',1.0,[41,67,1,1]), (self.img+'.mask',1.0,[60,30,3,0])])

          self.assertTrue(self.check_final(report))

#     def test_mask_outregion(self):
#          """ [mask] test_mask_outregion : Input mask has region that goes outside the image """
#          self.prepData('refim_twochan.ms')
#          mstr = 'circle[[50pix,110pix],20pix]'
#          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=mstr)
#          report=self.th.checkall(imgexist=[self.img+'2.mask'], imgval=[(self.img+'2.mask',0.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,95,0,0])])

     def test_mask_zerostart(self):
          """ [mask] test_mask_zerostart : Test that a zero starting mask is caught  """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='10.0arcsec',niter=0,interactive=0,parallel=self.parallel)
          os.system('cp -r ' + self.img + '.residual '+ self.img+'2.inpmask')
          _ia.open(self.img+'2.inpmask')
          pix =_ia.getchunk()
          pix.fill(0.0)
          _ia.putchunk(pix)
          _ia.close()

          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='10.0arcsec',niter=10,interactive=0,mask=self.img+'2.inpmask')

          report=self.th.checkall(ret=ret, imgexist=[self.img+'2.mask'], imgval=[(self.img+'2.model',0.0,[50,50,0,0]),(self.img+'2.mask',0.0,[50,50,0,0])], stopcode=7)

          self.assertTrue(self.check_final(report))

# the option, auto-thresh removed
#     def test_mask_zeroauto(self):
#          """ [mask] test_mask_zeroauto : Test that an automask-generated zero mask is caught  """
#          self.prepData('refim_point.ms')
#          ret = tclean(vis=self.msfile, imagename=self.img,niter=0,interactive=0,usemask='auto-thresh',maskthreshold='40.0Jy')
#          ret = tclean(vis=self.msfile, imagename=self.img,niter=10,interactive=0,usemask='auto-thresh',maskthreshold='40.0Jy')
#
#          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'], imgval=[(self.img+'.model',0.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,50,0,0])], stopcode=7)
#
#          self.assertTrue(self.check_final(report))

     def test_mask_expand_contstokesImask_to_cube(self):
          """ [mask] test_mask_expand_contstokesImask_to_cube : Test for
          expanding input continuum Stokes I mask to cube imaging  """
          self.prepData('refim_point_linRL.ms')
          self.prepInputmask('refim_cont_stokesI_input.mask')
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec', niter=10,interactive=0,interpolation='nearest', usemask='user', mask=self.maskname)

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,50,0,1]),(self.img+'.mask',1.0,[50,50,0,2]), (self.img+'.mask',0.0,[65,65,0,1])])

          self.assertTrue(self.check_final(report))

     def test_mask_expand_contstokesImask_nodegen_to_cube(self):
          """ [mask] test_mask_expand_contstokesImask_nodegen_to_cube : Test for
          expanding input continuum Stokes I mask with its degenerate axes removed to cube imaging  """
          self.prepData('refim_point_linRL.ms')
          self.prepInputmask('refim_cont_stokesI_input.mask')
          imsubimage(imagename=self.maskname, outfile=self.maskname+"_dropdeg",dropdeg=True, overwrite=True)
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', usemask='user',
          mask=self.maskname+"_dropdeg")
          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
                                  imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,50,0,1]),(self.img+'.mask',1.0,[50,50,0,2]), (self.img+'.mask',0.0,[65,65,0,1])])

          self.assertTrue(self.check_final(report))

        
     def test_mask_expand_contstokesImask_to_IQUV(self):
          """ [mask] test_mask_expand_contstokesImask_to_IQUV : Test for expanding
          input continuum Stokes I mask to continuum multi-stokes imaging  """
          self.prepData('refim_point_linRL.ms')
          self.prepInputmask('refim_cont_stokesI_input.mask')
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="mfs", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0, stokes='IQUV', usemask='user', mask=self.maskname)

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,50,1,0]),(self.img+'.mask',1.0,[50,50,2,0]),(self.img+'.mask',1.0,[50,50,3,0]), (self.img+'.mask',0.0,[65,65,2,0])])

          self.assertTrue(self.check_final(report))

     def test_mask_expand_contstokesImask_nodegen_to_IQUV(self):
          """ [mask] test_mask_expand_contstokesImask_nodegen_to_IQUV : Test for expanding
          input continuum Stokes I mask with its degenerate axes removed to continuum multi-stokes imaging  """
          self.prepData('refim_point_linRL.ms')
          self.prepInputmask('refim_cont_stokesI_input.mask')
          imsubimage(imagename=self.maskname, outfile=self.maskname+"_dropdeg", dropdeg=True, overwrite=True)
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="mfs", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0, stokes='IQUV', usemask='user',
          mask=self.maskname+"_dropdeg")

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,50,1,0]),(self.img+'.mask',1.0,[50,50,2,0]),(self.img+'.mask',1.0,[50,50,3,0]), (self.img+'.mask',0.0,[65,65,2,0])])

          self.assertTrue(self.check_final(report))


     def test_mask_expand_contstokesImask_to_cube_IQUV(self):
          """ [mask] test_mask_extend_contstokesImask_to_cube_IQUV : Test for extending
          input continuum Stokes I mask to cube multi-stokes imaging  """
          self.prepData('refim_point_linRL.ms')
          self.prepInputmask('refim_cont_stokesI_input.mask')
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', stokes='IQUV', usemask='user', mask=self.maskname)

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),
                 (self.img+'.mask',1.0,[50,50,1,0]),
                 (self.img+'.mask',1.0,[50,50,2,0]),
                 (self.img+'.mask',1.0,[50,50,3,0]), 
                 (self.img+'.mask',1.0,[50,50,0,1]),
                 (self.img+'.mask',1.0,[50,50,1,1]),
                 (self.img+'.mask',1.0,[50,50,2,1]),
                 (self.img+'.mask',1.0,[50,50,3,1]),
                 (self.img+'.mask',1.0,[50,50,1,2]),
                 (self.img+'.mask',0.0,[65,65,0,0]),
                 (self.img+'.mask',0.0,[65,65,2,1]),
                 ])

          self.assertTrue(self.check_final(report))


     def test_mask_expand_contstokesImask_nodegen_to_cube_IQUV(self):
          """ [mask] test_mask_extend_contstokesImask_nodegen_to_cube_IQUV : Test for extending
          input continuum Stokes I mask with its denenerate axes removed to cube multi-stokes imaging  """
          self.prepData('refim_point_linRL.ms')
          self.prepInputmask('refim_cont_stokesI_input.mask')
          imsubimage(imagename=self.maskname, outfile=self.maskname+"_dropdeg",dropdeg=True, overwrite=True)
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', stokes='IQUV',
          usemask='user', mask=self.maskname+"_dropdeg")

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),
                 (self.img+'.mask',1.0,[50,50,1,0]),
                 (self.img+'.mask',1.0,[50,50,2,0]),
                 (self.img+'.mask',1.0,[50,50,3,0]), 
                 (self.img+'.mask',1.0,[50,50,0,1]),
                 (self.img+'.mask',1.0,[50,50,1,1]),
                 (self.img+'.mask',1.0,[50,50,2,1]),
                 (self.img+'.mask',1.0,[50,50,3,1]),
                 (self.img+'.mask',1.0,[50,50,1,2]),
                 (self.img+'.mask',0.0,[65,65,0,0]),
                 (self.img+'.mask',0.0,[65,65,2,1]),
                 ])

          self.assertTrue(self.check_final(report))

     def test_mask_expand_contstokesIQUVmask_to_cube_IQUV(self):
          """ [mask] test_mask_expand_contstokesIQUVmask_to_cube_IQUV : Test for expanding
          input continuum Stokes IQUV mask to cube IQUV imaging  """
          # extending to all channels and preserving mask of each stokes 
          self.prepData('refim_point_linRL.ms') 
          # input mask will different for different stokes plane
          self.prepInputmask('refim_cont_stokesIQUV_input.mask')
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', stokes='IQUV', usemask='user', mask=self.maskname)

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),
                 (self.img+'.mask',1.0,[50,50,0,1]),
                 (self.img+'.mask',1.0,[50,50,1,0]),
                 (self.img+'.mask',1.0,[50,50,1,2]), 
                 (self.img+'.mask',1.0,[50,50,2,0]),
                 (self.img+'.mask',1.0,[50,50,2,1]),
                 (self.img+'.mask',1.0,[43,31,3,0]),
                 (self.img+'.mask',1.0,[43,31,3,2]),
                 (self.img+'.mask',0.0,[61,51,0,0]),
                 (self.img+'.mask',0.0,[50,63,1,1]),
                 (self.img+'.mask',0.0,[37,65,2,2]),
                 (self.img+'.mask',0.0,[34,70,0,1]),
                 ])

          self.assertTrue(self.check_final(report))


     def test_mask_expand_contstokesIQUVmask_nodegen_to_cube_IQUV(self):
          """ [mask] test_mask_expand_contstokesIQUVmask_nodegen_to_cube_IQUV : Test for expanding
          input continuum Stokes IQUV mask with its degenerate axes removed to cube IQUV imaging  """
          # extending to all channels and preserving mask of each stokes 
          self.prepData('refim_point_linRL.ms') 
          # input mask will different for different stokes plane
          self.prepInputmask('refim_cont_stokesIQUV_input.mask')
          imsubimage(self.maskname, outfile=self.maskname+"_dropdeg",dropdeg=True, overwrite=True);
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', stokes='IQUV',
          usemask='user', mask=self.maskname+"_dropdeg")

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),
                 (self.img+'.mask',1.0,[50,50,0,1]),
                 (self.img+'.mask',1.0,[50,50,1,0]),
                 (self.img+'.mask',1.0,[50,50,1,2]), 
                 (self.img+'.mask',1.0,[50,50,2,0]),
                 (self.img+'.mask',1.0,[50,50,2,1]),
                 (self.img+'.mask',1.0,[43,31,3,0]),
                 (self.img+'.mask',1.0,[43,31,3,2]),
                 (self.img+'.mask',0.0,[61,51,0,0]),
                 (self.img+'.mask',0.0,[50,63,1,1]),
                 (self.img+'.mask',0.0,[37,65,2,2]),
                 (self.img+'.mask',0.0,[34,70,0,1]),
                 ])

          self.assertTrue(self.check_final(report))

     def test_mask_expand_cubestokesImask_to_cube_IQUV(self):
          """ [mask] test_mask_expand_contstokesIQUVmask_to_cube_IQUV : Test for expanding
          input cube Stokes I mask to cube (of the same spectral coordinates)  IQUV imaging  """
          # extending to all channels and preserving mask of each stokes 
          self.prepData('refim_point_linRL.ms') 
          # input mask will different for different stokes plane
          self.prepInputmask('refim_cube_StokesI_input.mask')
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', stokes='IQUV',
          usemask='user', mask=self.maskname)

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),
                 (self.img+'.mask',1.0,[50,50,0,1]),
                 (self.img+'.mask',1.0,[50,50,1,0]),
                 (self.img+'.mask',1.0,[50,50,1,2]), 
                 (self.img+'.mask',1.0,[50,50,2,0]),
                 (self.img+'.mask',1.0,[50,50,2,1]),
                 (self.img+'.mask',1.0,[37,63,3,0]),
                 (self.img+'.mask',0.0,[46,49,3,2]),
                 (self.img+'.mask',0.0,[63,51,0,0]),
                 (self.img+'.mask',0.0,[50,63,1,1]),
                 (self.img+'.mask',0.0,[43,59,2,2]),
                 (self.img+'.mask',1.0,[43,57,2,2]),
                 (self.img+'.mask',0.0,[43,70,0,1]),
                 ])

          self.assertTrue(self.check_final(report))

     def test_mask_expand_cubestokesImask_nodegen_to_cube_IQUV(self):
          """ [mask] test_mask_expand_contstokesIQUVmask_nodegen_to_cube_IQUV : Test for expanding
          input cube Stokes I mask with its degenerate axes removed to cube (of the same spectral coordinates)  IQUV imaging  """
          # extending to all channels and preserving mask of each stokes 
          self.prepData('refim_point_linRL.ms') 
          # input mask will different for different stokes plane
          self.prepInputmask('refim_cube_StokesI_input.mask')
          imsubimage(self.maskname, outfile=self.maskname+"_dropdeg",dropdeg=True, overwrite=True);
          ret = tclean(vis=self.msfile,
          imagename=self.img, specmode="cube", imsize=100, cell='8.0arcsec',
          niter=10,interactive=0,interpolation='nearest', stokes='IQUV',
          usemask='user', mask=self.maskname+'_dropdeg')

          report=self.th.checkall(ret=ret, imgexist=[self.img+'.mask'],
          imgval=[(self.img+'.mask',1.0,[50,50,0,0]),
                 (self.img+'.mask',1.0,[50,50,0,1]),
                 (self.img+'.mask',1.0,[50,50,1,0]),
                 (self.img+'.mask',1.0,[50,50,1,2]), 
                 (self.img+'.mask',1.0,[50,50,2,0]),
                 (self.img+'.mask',1.0,[50,50,2,1]),
                 (self.img+'.mask',1.0,[37,63,3,0]),
                 (self.img+'.mask',0.0,[46,49,3,2]),
                 (self.img+'.mask',0.0,[63,51,0,0]),
                 (self.img+'.mask',0.0,[50,63,1,1]),
                 (self.img+'.mask',0.0,[43,59,2,2]),
                 (self.img+'.mask',1.0,[43,57,2,2]),
                 (self.img+'.mask',0.0,[43,70,0,1]),
                 ])

          self.assertTrue(self.check_final(report))

##############################################
##############################################

class test_wproject(testref_base):

     def test_wterm_wproject(self):
          """ [wproject] Test_Widefield_wproj : W-Projection """ 
          self.prepData("refim_point_wterm_vlad.ms")
          msname = self.msfile
          #msname = '/home/vega/rurvashi/TestCASA/VerificationTests/WProjection/refim_point_wterm_vlad.ms'

          ## Without w-term corrections, the source peak will be 0.768
          #tclean(vis=msname, imagename=self.img+'wno', imsize=2048, cell='10.0arcsec',niter=0, weighting='uniform', gridder='standard', pblimit=-0.1)

          tclean(vis=msname, imagename=self.img+'.wyes',  imsize=2048, cell='10.0arcsec',niter=0, weighting='uniform', gridder='wproject', wprojplanes=16, pblimit=-0.1,parallel=self.parallel)

          report=self.th.checkall(imgexist=[self.img+'.wyes.image'],imgval=[(self.img+'.wyes.psf',1.0,[1024,1024,0,0]),(self.img+'.wyes.image',1.0,[1158,1384,0,0]) ], tfmask=[(self.img+'.wyes.image',['T'])] )
          self.assertTrue(self.check_final(report))

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Facetted imaging tests parallel are skipped temporarily until a fix is found. ")
     def test_wterm_facets(self):
          """ [wproject] Test_Widefield_wproj : Facets """ 
          self.prepData("refim_point_wterm_vlad.ms")
          msname = self.msfile
          #msname = '/home/vega/rurvashi/TestCASA/VerificationTests/WProjection/refim_point_wterm_vlad.ms'

          tclean(vis=msname, imagename=self.img+'.facet',  imsize=2048, cell='10.0arcsec',niter=0, weighting='uniform', gridder='widefield', wprojplanes=1,facets=4, pblimit=-0.1,parallel=self.parallel)
          
          ## Current value with facets=4 is 0.988. 
          report=self.th.checkall(imgexist=[self.img+'.facet.image'],imgval=[(self.img+'.facet.psf',1.0,[1024,1024,0,0]),(self.img+'.facet.image',1.0,[1158,1384,0,0]) ] )
          self.assertTrue(self.check_final(report))


     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Facetted imaging tests in parallel are skipped temporarily until a fix is found. ")
     def test_wterm_wproject_facets(self):
          """ [wproject] Test_Widefield_wproj : Facets with wprojection per facet""" 
          self.prepData("refim_point_wterm_vlad.ms")
          msname = self.msfile
          #msname = '/home/vega/rurvashi/TestCASA/VerificationTests/WProjection/refim_point_wterm_vlad.ms'

          tclean(vis=msname, imagename=self.img+'.wp.facet',  imsize=2048, cell='10.0arcsec',niter=0, weighting='uniform', gridder='widefield', wprojplanes=4,facets=4, pblimit=-0.1,parallel=self.parallel)
          
          ## Current value with facets=4 is 0.988. 
          report=self.th.checkall(imgexist=[self.img+'.wp.facet.image'],imgval=[(self.img+'.wp.facet.psf',1.0,[1024,1024,0,0]),(self.img+'.wp.facet.image',1.0,[1158,1384,0,0]) ] )
          self.assertTrue(self.check_final(report))

  
     @unittest.skip('Skip test for wterm imaging with awproject until the numerical error has been addressed in CAS-13191')
     def test_wterm_awproject(self):
          """ [wproject] Test_Widefield_wproj : W-Projection using the AWProject gridder """ 
          self.prepData("refim_point_wterm_vlad.ms")
          msname = self.msfile
          #msname = '/home/vega/rurvashi/TestCASA/VerificationTests/WProjection/refim_point_wterm_vlad.ms'

          ### Peak value with gridder='awproject' comes out as 0.85479 instead of the 1.0 (0.998) that is made by gridder='wproject'. 
          tclean(vis=msname, imagename=self.img+'.awp',  imsize=2048, cell='10.0arcsec',niter=0, weighting='uniform', gridder='awproject', wprojplanes=16, pblimit=-0.1, psterm=True, aterm=False, wbawp=False,cfcache=self.img+'_use_awp.cf',parallel=self.parallel)

          report=self.th.checkall(imgexist=[self.img+'.awp.image'],imgval=[(self.img+'.awp.psf',1.0,[1024,1024,0,0]),(self.img+'.awp.image',1.0,[1158,1384,0,0]) ] )
          self.assertTrue(self.check_final(report))

          
  



##############################################
##############################################

##Task level tests : awproject and mosaics
class test_widefield(testref_base):
     def test_widefield_aproj_mfs(self):
          """ [widefield] Test_Widefield_aproj : MFS with narrowband AWProjection (wbawp=F, 1spw)  stokes I """
          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='1',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=30,gridder='awproject',cfcache='',wbawp=False,conjbeams=True,psterm=False,computepastep=360.0,
                       rotatepastep=360.0,deconvolver='hogbom',savemodel='modelcolumn',parallel=self.parallel)
         ## ret = tclean(vis=self.msfile,spw='2',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',wbawp=False,conjbeams=True,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom')
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.96,[256,256,0,0]),(self.img+'.pb',0.96,[256,256,0,0]),(self.img+'.weight',0.496,[256,256,0,0]) ] )
          #
          # Changed to the following for 5.5.0 release of AWP.  Will revisit and replace the test MS later.
          #
          #report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.705,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
          ## weight is pbsq which is 0.7^2 = 0.49 (approx).
          self.assertTrue(self.check_final(report))

          #do stokes V too.....
##     @unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_widefield_aproj_cube(self):
          """ [widefield] Test_Widefield_aproj_cube_aproj : Cube with AW-Projection  and rotation off """

          casalog.post("EMPTY TEST")
          return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       specmode='cube',niter=1,gain=1.0,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,
                       conjbeams=False,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.11,[256,256,0,0]),(self.img+'.weight',0.34,[256,256,0,0]) ] )
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
          self.assertTrue(self.check_final(report))

     ## Test normtype too somewhere..

     def test_widefield_wbaproj_mfs(self):
          """ [widefield] Test_Widefield_wbaproj_mfs : MFS with wideband AWProjection (wbawp=T, allspw) and nt=1 stokes I  """

          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=30,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=True,psterm=False,computepastep=360.0,
                       rotatepastep=360.0,deconvolver='hogbom',pblimit=0.3,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
          #
          # Changed to the following for 5.5.0 release of AWP.  Will revisit and replace the test MS later.
          #
          #report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.698,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
          self.assertTrue(self.check_final(report))


          #do stokes V too..

     def test_widefield_aproj_mtmfs(self):
          """ [widefield] Test_Widefield_aproj_mtmfs : MFS with AWProjection (wbawp=T,conjbeams=F, allspw) and nt=2 stokes I  """

          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=0,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=False,psterm=False,computepastep=360.0,
                       rotatepastep=360.0,deconvolver='mtmfs',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.96,[256,256,0,0]),(self.img+'.weight.tt0',0.48,[256,256,0,0]),(self.img+'.alpha',0.04,[256,256,0,0]) ] )
          #
          # Changed to the following for 5.5.0 release of AWP.  Will revisit and replace the test MS later.
          #
          #report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.705,[256,256,0,0]),(self.img+'.weight.tt0',0.48,[256,256,0,0]),(self.img+'.alpha',0.04,[256,256,0,0]) ] )
          self.assertTrue(self.check_final(report))
          ## alpha should represent that of the mosaic PB (twice).. -0.1 doesn't look right. Sigh.... well.. it should converge to zero.
          ## alpha keeps increasing in magnitude with niter.... not right.
          ## restricting this check to niter=0 only. Alpha should represent that of the PB.

     def test_widefield_wbaproj_mtmfs(self):
          """ [widefield] Test_Widefield_wbaproj_mtmfs : MFS with wideband AWProjection (wbawp=T,conjbeams=T, allspw) and nt=2 stokes I  """

          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=30,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=True,psterm=False,
                       computepastep=360.0,rotatepastep=360.0,deconvolver='mtmfs',pblimit=0.1,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.96,[256,256,0,0]),(self.img+'.weight.tt0',0.486,[256,256,0,0]),(self.img+'.alpha',0.0,[256,256,0,0]) ] )
          #
          # Changed to the following for 5.5.0 release of AWP.  Will revisit and replace the test MS later.
          #
          #report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.696,[256,256,0,0]),(self.img+'.weight.tt0',0.486,[256,256,0,0]),(self.img+'.alpha',0.0,[256,256,0,0]) ] )
          ## alpha should be ZERO as the pb spectrum has been taken out.
          self.assertTrue(self.check_final(report))


#     def test_widefield_wbaproj_subsets(self):
#          """ [widefield] Test_Widefield_wbaproj_subsets : MFS with the AWProjection gridder and A,W turned off  """
#          self.prepData("refim_mawproject.ms")
#          ## PS only
#          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',psterm=True,aterm=False,wprojplanes=1,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',pblimit=0.3)
#          #report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
#
#          ## W and PS only
#          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',psterm=True,aterm=False,wprojplanes=16,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',pblimit=0.3)
#          #report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )


#     def test_widefield_multispws(self):
#          """ [widefield] Test_Widefield_multispws : Test cube imaging with mosaicft and awproj  """

     ## CHECK NORMALIZATION OF WEIGHTIMAGE = normed to peak=1
     ## TODO : make vpman recognize EVLA in addition to VLA.
     def test_widefield_mosaicft_mfs(self):
          """ [widefield] Test_Widefield_mosaic : MFS with mosaicft  stokes I """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='1',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=30,gridder='mosaicft',deconvolver='hogbom',pblimit=0.3,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.961231,[256,256,0,0]),(self.img+'.weight',0.50576,[256,256,0,0]) ] , tfmask=[(self.img+'.image',['mask0'])] )
          #ret = clean(vis=self.msfile,spw='1',field='*',imagename=self.img+'.old',imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,imagermode='mosaic',psfmode='hogbom')
          self.assertTrue(self.check_final(report))

          #do stokes V too..

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. mosaic, Briggs weighting with mosweight=True. Enable this after fixing CAS-11978")
     def test_widefield_mosaicft_mfs_mosweightTrue(self):
          """ [widefield] Test_Widefield_mosaic : MFS with mosaicft  stokes I briggs mosweight=True(default)"""
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='1',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=30,gridder='mosaicft',deconvolver='hogbom',pblimit=0.3,weighting='briggs', parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.962813, [256,256,0,0]),(self.img+'.weight',0.50520, [256,256,0,0]) ] )
          #ret = clean(vis=self.msfile,spw='1',field='*',imagename=self.img+'.old',imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,imagermode='mosaic',psfmode='hogbom')
          self.assertTrue(self.check_final(report))

     def test_widefield_mosaicft_mtmfs(self):
          """ [widefield] Test_Widefield_mosaicft_mtmfs : MT-MFS with mosaicft  stokes I, alpha """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=60,gridder='mosaicft',deconvolver='mtmfs', conjbeams=False,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.9413,[256,256,0,0]),(self.img+'.weight.tt0',0.50546,[256,256,0,0]),(self.img+'.alpha', 0.07367,[256,256,0,0]) ] , tfmask=[(self.img+'.image.tt0',['mask0'])] )
          ## alpha should represent that of the mosaic PB (twice)... and should then converge to zero
          self.assertTrue(self.check_final(report))
          
     def test_widefield_mosaicft_mtmfs_conj(self):
          """ [widefield] Test_Widefield_mosaicft_mtmfs : MT-MFS with mosaicft  stokes I, alpha """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=60,gridder='mosaicft',deconvolver='mtmfs', conjbeams=True,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.9638,[256,256,0,0]),(self.img+'.weight.tt0',0.49804,[256,256,0,0]),(self.img+'.alpha',-0.03692,[256,256,0,0]) ] )
          ## alpha should represent that of the mosaic PB (twice)... and should then converge to zero
          self.assertTrue(self.check_final(report))
     
     def test_widefield_mosaicft_mtmfs_pbsquare(self):
          """ [widefield] Test_Widefield_mosaicft_mtmfs : MT-MFS with mosaicft  stokes I, alpha """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       niter=60,gridder='mosaic',deconvolver='mtmfs', conjbeams=False, normtype='pbsquare',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imgval=[(self.img+'.image.tt0',0.9194,[256,256,0,0]),(self.img+'.weight.tt0',0.5059,[256,256,0,0]),(self.img+'.alpha',0.021195,[256,256,0,0]) ] )
          ## alpha should represent that of the mosaic PB (twice)... and should then converge to zero
          self.assertTrue(self.check_final(report))

     def test_widefield_mosaicft_cube(self):
          """ [widefield] Test_Widefield_mosaicft_cube : MFS with mosaicft  stokes I """
          self.prepData("refim_mawproject.ms")
#          _vp.setpbpoly(telescope='EVLA', coeff=[1.0, -1.529e-3, 8.69e-7, -1.88e-10]) 
#          _vp.saveastable('evlavp.tab')
          ret = tclean(vis=self.msfile,spw='*',field='0',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",
                       specmode='cube',niter=10,gridder='mosaicft',deconvolver='hogbom',gain=0.1,stokes='I',parallel=self.parallel) #,vptable='evlavp.tab')
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.7987,[256,256,0,0]),(self.img+'.weight',0.6528,[256,256,0,0]) ] )
          self.assertTrue(self.check_final(report))

          #do stokes V too..
     def test_mosaicft_newpsfphasecenter(self):
          """
          test_mosaicft_newpsfphasecenter : different phasecenter for psf
          """
          self.prepData("refim_mawproject.ms")
          ret=tclean(vis="refim_mawproject.ms",field="*",spw="1",datacolumn="corrected",imagename=self.img,imsize=512,cell="10.0arcsec",phasecenter="J2000 19:59:28.500 +40.44.01.50",stokes="I",specmode="mfs",gridder="mosaic",psfphasecenter="J2000 19:59:28.520 +40.44.01.51",vptable="",pblimit=0.3,normtype="flatnoise",deconvolver="hogbom",restoration=True,weighting="natural", niter=30,gain=0.1, usemask="user",mask="",restart=True,savemodel="none",calcres=True,calcpsf=True, parallel=self.parallel)

          # need to add check of the actual coordinates of the peak of psf. That would match with psfphasecenter value...
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imgval=[(self.img+'.image',0.96,[256,256,0,0]), (self.img+'.psf',1.0,[256,256,0,0])])

          self.assertTrue(self.check_final(report))

     def test_mosaicft_newpsfphasecenter_cube(self):
          """
          test_mosaicft_newpsfphasecenter_cube : different phasecenter for psf
          """
          self.prepData("refim_mawproject.ms")
          ret=tclean(vis="refim_mawproject.ms",field="*",spw="*",datacolumn="corrected",imagename=self.img,imsize=512,cell="10.0arcsec",phasecenter="J2000 19:59:28.500 +40.44.01.50",stokes="I",specmode="cube",gridder="mosaic",psfphasecenter="J2000 19:59:28.520 +40.44.01.51",vptable="",pblimit=0.3,normtype="flatnoise",deconvolver="hogbom",restoration=True,weighting="natural", niter=30,gain=0.1, usemask="user",mask="",restart=True,savemodel="none",calcres=True,calcpsf=True, parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'], imgval=[(self.img+'.image',0.99,[256,256,0,0]), (self.img+'.psf',1.0,[256,256,0,0])])
          self.assertTrue(self.check_final(report))


     
##############################################
##############################################

##Task level tests : model prediction.
class test_modelvis(testref_base):
     
     def test_modelvis_1(self):
          """ [modelpredict] Test_modelvis_1 : mfs with no save model """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='none',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==False )

     #Passes for parallel after CAS-9386 
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_2(self):
          """ [modelpredict] Test_modelvis_2 : mfs with save model column """
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='modelcolumn',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

          ##Predict from input model image (startmodel)
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', niter=0,
                       savemodel='modelcolumn',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     #Passes with parallel after CAS-9386
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_3(self):
          """ [modelpredict] Test_modelvis_3 : mfs with save virtual model """
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

          ##Predict from input model image (startmodel)
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', niter=0,
                       savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_4(self):
          """ [modelpredict] Test_modelvis_4 : mt-mfs with no save model """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',savemodel='none',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==False )

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled() or is_CASA6, "Skip the test temporarily, plotms unavailable in casatasks")
     def test_modelvis_5(self):
          """ [modelpredict] Test_modelvis_5 : mt-mfs with save model column """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',
                       savemodel='modelcolumn',parallel=self.parallel)
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='data',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.data.png',title="original data")
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.model.png',title="empty model")
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'],
                       niter=0,deconvolver='mtmfs',savemodel='modelcolumn',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     #tablelock issue with n=2 (mpicasa -n 3 or 4 worked) on CAS-9386
     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_6(self):
          """ [modelpredict] Test_modelvis_6 : mt-mfs with save virtual model """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'],
                       niter=0,deconvolver='mtmfs',savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     # Passes after CAS-9386
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_7(self):
          """ [modelpredict] Test_modelvis_7 : cube with chan selection and save model column """
          ## check explicit channels ...
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0:5~12',imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,savemodel='modelcolumn',
                       start=5,nchan=8,interpolation='nearest',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )
          reportcv=self.th.check_chanvals(self.msfile, [(10,">",0.0),(3,"==",1.0)])
          self.assertTrue(self.check_final(reportcv))
          
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0',imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',specmode='cube',
                       niter=0,savemodel='modelcolumn',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )
          reportcv=self.th.check_chanvals(self.msfile,[(10,">",0.0),(3,"==",self.th.check_modelchan(self.msfile,1))])
          self.assertTrue(self.check_final(reportcv))


     #Passes after CAS-9386
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_8(self):
          """ [modelpredict] Test_modelvis_8 : cube with chan selection and save virtual model """
          ## check explicit channels ...
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0:5~12',imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,
                       savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0',imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',specmode='cube',
                       niter=0,savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     # tested on CAS-9386, now seems to work...
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_9(self):
          """ [modelpredict] Test_modelvis_9 : Don't de-grid channels with zero model. Also test limited-freq mask """
          self.prepData("refim_point.ms")
          masklist=''  # choose only a few channels here.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,mask=masklist,parallel=self.parallel)
          self.assertTrue(self.th.exists(self.img+'.model') )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',niter=0,
                       savemodel='modelcolumn',parallel=self.parallel)

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',niter=0,
                       savemodel='virtual',parallel=self.parallel)

     # tested on CAS-9386, now seems to work...
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_10(self):
          """ [modelpredict] Test_modelvis_10 : Use input model of different (narrower) freq range than data """
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', spw='0:5~12',niter=10,savemodel='modelcolumn',parallel=self.parallel)
          self.assertTrue(self.th.exists(self.img+'.model') )
#          self.assertTrue( self.th.check_modelchan(self.msfile,10) > 0.0 and self.th.check_modelchan(self.msfile,3) == 1.0 )
          reportcv=self.th.check_chanvals(self.msfile,[(10,">",0.0),(3,"==",1.0)])
          self.assertTrue(self.check_final(reportcv))


          ## add model expansion parameter
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,
                       savemodel='modelcolumn',parallel=self.parallel)
#          self.assertTrue( self.th.check_modelchan(self.msfile,10) > 0.0 and self.th.check_modelchan(self.msfile,3) > 0.0 )
          reportcv=self.th.check_chanvals(self.msfile,[(10,">",0.0),(3,">",0.0)])
          self.assertTrue(self.check_final(reportcv))
                    

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,
                       savemodel='virtual',parallel=self.parallel)
          ## cannot check anything here....  just that it runs without error

     # mfs case but fixed? with CAS-9386
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_modelvis_11(self):
          """ [modelpredict] Test_modelvis_11 : Predict model image over channel gaps not included in imaging """
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', spw='0:0~8;12~19',niter=10,savemodel='modelcolumn',parallel=self.parallel)
          self.assertTrue(self.th.exists(self.img+'.model') )
          ###vi2 leave unselected channel as is so it will be 1.0
          self.assertTrue( (self.th.check_modelchan(self.msfile,10) == 0.0) or (np.abs(self.th.check_modelchan(self.msfile,10)-1) < 1.0e-12)   and self.th.check_modelchan(self.msfile,3) > 0.0 )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,
                       savemodel='modelcolumn',parallel=self.parallel)
          self.assertTrue( self.th.check_modelchan(self.msfile,10) > 0.0 and self.th.check_modelchan(self.msfile,3) > 0.0 )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,
                       savemodel='virtual',parallel=self.parallel)
          ## cannot check anything here....  just that it runs without error

          
     def test_modelvis_12(self):
          """ [modelpredict] Test_modelvis_12 : (CAS-12618) mfs with automask and save model column (saving model via a separate predict model step)"""
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='modelcolumn',usemask='auto-multithresh', parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     def test_modelvis_13(self):
          """ [modelpredict] Test_modelvis_13 : (CAS-12618) mfs with automask and save virtual model (saving model via a separate predict model step)"""
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='virtual',usemask='auto-multithresh', parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_14(self):
          """ [modelpredict] Test_modelvis_14 : (cas-12618) mfs with nsigma and save model column (saving model via a separate predict model step)"""
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='modelcolumn',nsigma=1.0, parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     def test_modelvis_15(self):
          """ [modelpredict] Test_modelvis_15 : (CAS-12618) mfs with nsigma and save model column (saving model via a separate predict model step)"""
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='virtual',nsigma=1.0, parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_16(self):
          """ [modelpredict] Test_modelvis_16 : (CAS-12618) cube with  and save  model column for auto-multithresh"""
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,
                       usemask='auto-multithresh', savemodel='modelcolumn',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     def test_modelvis_17(self):
          """ [modelpredict] Test_modelvis_17 : (CAS-12618) cube with  and save virtual model for auto-multithreseh"""
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,
                       usemask='auto-multithresh', savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_18(self):
          """ [modelpredict] Test_modelvis_18 : (CAS-12618) cube with  and save model column for nsgima >0.0"""
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,
                       nsigma=1.0, savemodel='modelcolumn',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     def test_modelvis_19(self):
          """ [modelpredict] Test_modelvis_19 : (CAS-12618) cube with  and save virtual model for nsima >0.0"""
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,
                       nsigma=1.0, savemodel='virtual',parallel=self.parallel)
          hasmodcol, modsum, hasvirmod = self.th.check_model(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

class test_startmodel(testref_base):
     def test_startmodel_regrid_mfs(self):
          """ [modelpredict] Test_startmodel_regrid_mfs : Regrid input model onto new image grid : mfs (ra/dec) """
          self.prepData('refim_twopoints_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=50,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,
                        phasecenter='J2000 19h58m40.801s +40d55m59.863s',parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=200,cell='8.0arcsec',niter=0,deconvolver='hogbom',interactive=0,
                        startmodel=self.img+'1.model',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual', self.img+'2.residual'], imgval=[(self.img+'1.residual',1.7963,[25,25,0,0]),(self.img+'2.residual',1.910,[168,190,0,0])])
          self.assertTrue(self.check_final(report))

     def test_startmodel_regrid_cube(self):
          """ [modelpredict] Test_startmodel_regrid_cube : Regrid input model onto new image grid : cube (ra/dec/specframe)"""
          self.prepData('refim_point.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=50,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',
                        start='1.05GHz',width='50MHz',nchan=20,parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,deconvolver='hogbom',interactive=0,
                        startmodel=self.img+'1.model',specmode='cube',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual', self.img+'2.residual'], imgval=[(self.img+'1.residual',0.362,[25,25,0,5]),(self.img+'2.residual',0.362,[50,50,0,6])])
          self.assertTrue(self.check_final(report))



#     def test_startmodel_14(self):
#          """ [modelpredict] Test_startmodel_14 : Regrid input model onto new image grid : mtmfs (ra/dec/terms)"""

#     def test_startmodel_15(self):
#          """ [modelpredict] Test_startmodel_15 : Regrid input model onto new image grid : mfs (imsize/cell)"""

     def test_startmodel_mfs_continue(self):
          """ [startmodel] test_startmodel_mfs_continue : Restart a run with no parameter changes"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual'], imgval=[(self.img+'1.residual',0.35304,[50,50,0,0])])
          ret2 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual'], imgval=[(self.img+'1.residual',0.1259,[50,50,0,0])])
          self.assertTrue(self.check_final(report))

     def test_startmodel_mfs_restart(self):
          """ [startmodel] test_startmodel_mfs_restart : Restart a run using 'startmodel' and changed imagename"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,
                        startmodel=self.img+'1.model',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual', self.img+'2.residual'], imgval=[(self.img+'1.residual',0.35304,[50,50,0,0]),(self.img+'2.residual',0.1259,[50,50,0,0])])
          self.assertTrue(self.check_final(report))

     def test_startmodel_mfs_changeshape_1(self):
          """ [startmodel] test_startmodel_mfs_changeshape_1 : Restart a run but change shape only (cas-6937)"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual'], imgval=[(self.img+'1.residual',0.35304,[50,50,0,0])])

          try:
               ## This run should fail with an exception (if __rethrow_exceptions = True )
               ret2 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=120,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
               correct=False
          except Exception as e:
               correct=True
          #self.assertTrue(correct)
          
          ## Check that there is no change in output value.... 
          ## i.e. the second run should have failed.     
          report=self.th.checkall(imgval=[(self.img+'1.residual',0.35304,[50,50,0,0])])
          self.assertTrue(self.check_final(report))

     def test_startmodel_mfs_changeshape_2(self):
          """ [startmodel] test_startmodel_mfs_changeshape_2 : Restart a run using 'startmodel' and change shape and imagename"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,parallel=self.parallel)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=120,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,
                        startmodel=self.img+'1.model',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'1.residual', self.img+'2.residual'], imgval=[(self.img+'1.residual',0.35304,[50,50,0,0]),(self.img+'2.residual',0.1259,[60,60,0,0])])
          self.assertTrue(self.check_final(report))

     def test_startmodel_mtmfs_restart(self):
          """ [startmodel] test_startmodel_mtmfs_restart : Restart a multi-term run using 'startmodel' and changed imagename"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,parallel=self.parallel)
          # start with full model
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,
                        startmodel=[self.img+'1.model.tt0',self.img+'1.model.tt1'],parallel=self.parallel)

          # start with model only for tt0
          ret3 = tclean(vis=self.msfile,imagename=self.img+'3',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,
                        startmodel=self.img+'1.model.tt0',parallel=self.parallel)
          # start with model only for tt1
          ret3 = tclean(vis=self.msfile,imagename=self.img+'4',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,
                        startmodel=['',self.img+'1.model.tt1'],parallel=self.parallel)

          report=self.th.checkall(imgexist=[self.img+'1.residual.tt0', self.img+'2.residual.tt0', self.img+'3.residual.tt0', self.img+'4.residual.tt0', self.img+'1.residual.tt1', self.img+'2.residual.tt1', self.img+'3.residual.tt1', self.img+'4.residual.tt1'], imgval=[  (self.img+'1.residual.tt0',0.39226,[50,50,0,0]),
                             (self.img+'2.residual.tt0',0.13677,[50,50,0,0]),
                             (self.img+'3.residual.tt0',0.13677,[50,50,0,0]),
                             (self.img+'4.residual.tt0',0.39226,[50,50,0,0]),  
                             (self.img+'1.residual.tt1',-0.04358,[50,50,0,0]),
                             (self.img+'2.residual.tt1',-0.01519,[50,50,0,0]),
                             (self.img+'3.residual.tt1',-0.04358,[50,50,0,0]),
                             (self.img+'4.residual.tt1',-0.01519,[50,50,0,0])     ] )

          self.assertTrue(self.check_final(report))

     ## Run this test only in parallel mode as it tests combinations of serial and parallel runs. 
     @unittest.skipIf(not ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_csys_startmodel_restart_mfs(self):
          """ [startmodel] test_csys_startmodel_restart_cube

          Run a sequence of tclean runs to trigger a complicated situation of restarts, mixing serial/parallel and model writes. 
          This sequence, coupled with the algorithm options listed below in tests #1 through #6 trigger three different errors that
          have been fixed in this branch, and one that will be addressed via CAS-9386 (cube refactor). 

          tclean call sequence : 
          --- (a) Parallel run for niter=0
          --- (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
          --- (c) Serial model-predict run (without/with startmodel) : in one case it reuses prev image-set. in other case it reuses only 'model'. 
          --- (d) Impbcor on the output of (b)

          Note that this is not a full fix of the various instances of the 'latpole' inconsistency, but only a workaround. 
          Hence it needs a test to ensure this keeps working. 
          """
          self.prepData('refim_point.ms')

          specmode='mfs'
          runpar=self.parallel

          # Test 1 :  mfs + hogbom/mtmfs + usestartmod=False : Triggers both the latpole errors 
          #                                                                      (a) when the serial modelpredict runs on the output of the prev parallel run.
          #                                                                      (b) when the ouput of a parallel run is used by impbcor.
          # Test 2 : mfs + hogbom/mtmfs + usestartmod=True : Triggers the problem of 'model' being insufficient to create a valid ImStore. 
          # Run for deconvolver='hogbom' and 'mtmfs' because it exercises both SIImageStore and SIImageStoreMultiTerm 
          #      where some infrastructure code is re-implemented for multi-term vs single-term. 

          report = ''

          for runtype in ['serial','parallel']:
               for deconvolver in ['hogbom','mtmfs']:
                    for usestartmod in [False,True]:
                         infostr = 'Run with '+specmode +' - ' + deconvolver + ' - usestartmodel = ' + str(usestartmod) + ' - imaging in ' + runtype 
                         print(infostr)
                         report = report + infostr+'\n'

                         # Clear the model column.
                         self.th.delmodels(msname=self.msfile,modcol='reset0')

                         os.system('rm -rf savemod*')

                         ## (a) Always parallel run for niter=0
                         tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19',   niter=0,  imsize=100, cell='10.0arcsec', deconvolver=deconvolver,  specmode=specmode,  parallel=True)
                         
                         ## (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
                         if runtype=='serial':
                              runpar=False
                         else:
                              runpar=True
                         tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19', niter=10, imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode,  calcpsf=False, calcres=False,  parallel=runpar)
                         
                         # (c) Serial model-predict run (without/with startmodel) : 
                         if usestartmod==False:    # Do the restart by re-using the same image name as before.
                              imname2 = 'savemod.par'
                              startmodel=[]
                         else:
                              imname2 = 'savemod.ser'  # Do the restart with a new imagename, and using 'startmodel'
                              if deconvolver=='mtmfs':
                                   startmodel=['savemod.par.model.tt0', 'savemod.par.model.tt1']
                              else:
                                   startmodel= ['savemod.par.model']

                         tclean( vis=self.msfile, imagename=imname2, spw='0:0~5,0:15~19', niter=0,  imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode, savemodel='modelcolumn', startmodel=startmodel, calcres=False,  calcpsf= False, restoration=False, parallel=False)
                         
                         # Check the values of the model column
                         report = report + self.th.check_chanvals(self.msfile,[(0,">",0.5),(19,"<",0.9)])
                         
                         ## (d) Impbcor on the output of (b) + check the output of pbcor
                         if deconvolver=='hogbom':
                              impbcor(imagename='savemod.par'+'.image', pbimage='savemod.par'+'.pb', outfile='savemod.par'+'.impbcor',overwrite=True)
                              report=report +self.th.checkall(imgexist=['savemod.par.impbcor'], imgval=[  ('savemod.par.impbcor',1.1,[50,50,0,0]) ])
                         else:
                              impbcor(imagename='savemod.par'+'.image.tt0', pbimage='savemod.par'+'.pb.tt0', outfile='savemod.par'+'.impbcor.tt0',overwrite=True)
                              report=report +self.th.checkall(imgexist=['savemod.par.impbcor.tt0'], imgval=[  ('savemod.par.impbcor.tt0',1.1,[50,50,0,0]) ])

          self.assertTrue(self.check_final(report))


     def test_csys_startmodel_restart_cube(self):
          """ [startmodel] test_csys_startmodel_restart_cube : Check that csys differences w.r.to latpoles for parallel vs serial runs are appropriately squashed. 

          Run a sequence of tclean runs to trigger a complicated situation of restarts, mixing serial/parallel and model writes. 
          This sequence, coupled with the algorithm options listed below in tests #1 through #6 trigger three different errors that
          have been fixed in this branch, and one that will be addressed via CAS-9386 (cube refactor). 

          tclean call sequence : 
          --- (a) Parallel run for niter=0
          --- (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
          --- (c) Serial model-predict run (without/with startmodel) : in one case it reuses prev image-set. in other case it reuses only 'model'. 
          --- (d) Impbcor on the output of (b)

          Note that this is not a full fix of the various instances of the 'latpole' inconsistency, but only a workaround. 
          Hence it needs a test to ensure this keeps working. 
          """
          self.prepData('refim_point.ms')

          specmode='cube'
          deconvolver='hogbom'
          runpar=self.parallel

          report=''

          # Test  : cube + hogbom + usestartmod=True/False :  Triggers the problem of refconcat image outputs being incompatible with restarts in serial later.
          ####for runtype in ['serial','parallel']:     ## This will fail in 'serial' because we cannot mix and match refconcat images with regular ones. 
          ####                                                                 This must be revisited after CAS-9386. 
          for runtype in ['parallel']:
               for usestartmod in [False,True]:
                    print('Run with '+specmode +' - ' + deconvolver + ' - ' + str(usestartmod))
                    
                    infostr = 'Run with '+specmode +' - ' + deconvolver + ' - usestartmodel = ' + str(usestartmod) + ' - imaging in ' + runtype 
                    print(infostr)
                    report = report + infostr+'\n'
                    
                    # Clear the model column.
                    self.th.delmodels(msname=self.msfile,modcol='reset0')
                    
                    os.system('rm -rf savemod*')
                    
                    ## (a) Always parallel run for niter=0
                    tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19',   niter=0,  imsize=100, cell='10.0arcsec', deconvolver=deconvolver,  specmode=specmode,  parallel=True)
                    
                    ## (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
                    if runtype=='serial':
                         runpar=False
                    else:
                         runpar=True
                    tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19', niter=10, imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode,  calcpsf=False, calcres=False,  parallel=runpar)
                         
                    # (c) Serial model-predict run (without/with startmodel) : 
                    if usestartmod==False:    # Do the restart by re-using the same image name as before.
                         imname2 = 'savemod.par'
                         startmodel=[]
                    else:
                         imname2 = 'savemod.ser'  # Do the restart with a new imagename, and using 'startmodel'
                         if deconvolver=='mtmfs':
                              startmodel=['savemod.par.model.tt0', 'savemod.par.model.tt1']
                         else:
                              startmodel= ['savemod.par.model']
                              
                    tclean( vis=self.msfile, imagename=imname2, spw='0:0~5,0:15~19', niter=0,  imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode, savemodel='modelcolumn', startmodel=startmodel, calcres=False,  calcpsf= False, restoration=False, parallel=False)
                              
                    # Check the values of the model column
                    report = report + self.th.check_chanvals(self.msfile,[(0,">",0.5),(19,"<",0.9)])
                    
                    ## (d) Impbcor on the output of (b) + check the output of pbcor (first channel)
                    impbcor(imagename='savemod.par'+'.image', pbimage='savemod.par'+'.pb', outfile='savemod.par'+'.impbcor',overwrite=True)
                    report=report +self.th.checkall(imgexist=['savemod.par.impbcor'], imgval=[  ('savemod.par.impbcor',1.5,[50,50,0,0]) ])


          self.assertTrue(self.check_final(report))


     # Now passes for n=3, 4, with CAS-9386 but table lock issue in model write for n=2

     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_startmodel_with_mask_mfs(self):
          """ [startmodel] test_startmodel_with_mask_mfs : Mask out some regions in the startmodel, before prediction """
          self.prepData("refim_twopoints_twochan.ms")
          ## image both sources
          tclean(vis=self.msfile,spw='0:0',imagename=self.img,imsize=200,cell='10.0arcsec',niter=100,parallel=self.parallel)
          ## mask out all but the brightest outlier source
          _ia.open(self.img+'.model')
          _ia.calcmask(mask='"'+self.img+'.model">2.0') # keep only high values. i.e. the far out bright source.
          _ia.close()
          ## predict only that outlier source
          tclean(vis=self.msfile,spw='0:0',imagename=self.img+'.2',niter=0,savemodel='modelcolumn',imsize=200,cell='10.0arcsec',
                 startmodel=self.img+'.model',parallel=self.parallel)
          ## subtract it out
          uvsub(vis=self.msfile)
          ## image the rest of the field.
          tclean(vis=self.msfile,spw='0:0',imagename=self.img+'.3',imsize=200,cell='10.0arcsec',niter=100,parallel=self.parallel)
          
          report=self.th.checkall(imgexist=[self.img+'.model',self.img+'.2.model',self.img+'.3.model'],
                        imgval=[ (self.img+'.model',1.497,[100,100,0,0]),
                                (self.img+'.model',7.474,[154,172,0,0]),
                                (self.img+'.3.model',1.497,[100,100,0,0]), 
                                (self.img+'.3.model',0.024,[154,172,0,0])   ] )
          self.assertTrue(self.check_final(report))
          
     #@unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_startmodel_with_mask_mtmfs(self):
          """ [startmodel] test_startmodel_with_mask_mtmfs : Mask out some regions in the startmodel, before prediction """
          self.prepData("refim_twopoints_twochan.ms")
          ## image both sources
          tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='10.0arcsec',niter=100,deconvolver='mtmfs',parallel=self.parallel)
          ## mask out all but the brightest outlier source
          _ia.open(self.img+'.model.tt0')
          _ia.calcmask(mask='"'+self.img+'.model.tt0">2.0') # keep only high values. i.e. the far out bright source.
          _ia.close()
          _ia.open(self.img+'.model.tt1')
          _ia.calcmask(mask='"'+self.img+'.model.tt0">2.0') # keep only high values. i.e. the far out bright source.
          _ia.close()
          ## predict only that outlier source
          tclean(vis=self.msfile,imagename=self.img+'.2',niter=0,savemodel='modelcolumn',imsize=200,cell='10.0arcsec',
                 startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'],deconvolver='mtmfs',parallel=self.parallel)
          ## subtract it out
          uvsub(vis=self.msfile)
          ## image the rest of the field.
          tclean(vis=self.msfile,imagename=self.img+'.3',imsize=200,cell='10.0arcsec',niter=100,deconvolver='mtmfs',parallel=self.parallel)
          
          report=self.th.checkall(imgexist=[self.img+'.model.tt0',self.img+'.2.model.tt0',self.img+'.3.model.tt0'],
                        imgval=[ 
                                (self.img+'.model.tt0',1.11,[100,100,0,0]),
                                (self.img+'.model.tt0',5.56,[154,172,0,0]),
                                (self.img+'.3.model.tt0',1.12,[100,100,0,0]), 
                                (self.img+'.3.model.tt0',0.013,[154,172,0,0]),
                                (self.img+'.model.tt1',-1.1027,[100,100,0,0]),
                                (self.img+'.model.tt1',-5.602,[154,172,0,0]),
                                (self.img+'.3.model.tt1',-1.124,[100,100,0,0]), 
                                (self.img+'.3.model.tt1',-0.021,[154,172,0,0])   ] )
          self.assertTrue(self.check_final(report))
          
 ##############################################
class test_pbcor(testref_base):
     def setUp(self):
          super(test_pbcor, self).setUp()

          _vp.setpbpoly(telescope='EVLA', coeff=[1.0, -1.529e-3, 8.69e-7, -1.88e-10]) 
          _vp.saveastable('evlavp.tab')

     def test_pbcor_mfs(self):
          """ [pbcor] Test pbcor with mfs"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.7,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0])])
          self.assertTrue(self.check_final(report))

     def test_pbcor_mtmfs(self):
          """ [pbcor] Test pbcor with mtmfs"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True, deconvolver='mtmfs',parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.pb.tt0'], imgexistnot=[self.img+'.image.tt0.pbcor', self.img+'.alpha.pbcor'], imgval=[(self.img+'.pb.tt0',0.7,[256,256,0,0])])  
          #report=self.th.checkall(imgexist=[self.img+'.image.tt0', self.img+'.pb.tt0', self.img+'.image.tt0.pbcor', self.img+'.alpha.pbcor'], imgval=[(self.img+'.pb.tt0',0.7,[256,256,0,0]),(self.img+'.image.tt0.pbcor',1.0,[256,256,0,0]),(self.img+'.alpha',-0.7,[256,256,0,0]), (self.img+'.alpha.pbcor',-0.7,[256,256,0,0]) ])  
          # uncorrected alpha, for now. 
          self.assertTrue(self.check_final(report))

     def test_pbcor_cube_basic(self):
          """ [pbcor] Test pbcor with cube"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='cube', vptable='evlavp.tab', pbcor=True,parallel=self.parallel)
          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.79,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2]),(self.img+'.image.pbcor',1.0,[256,256,0,2])])
          self.assertTrue(self.check_final(report))

     def test_pbcor_cube_twosteps(self):
          """ [pbcor] Test pbcor with cube with imaging and pbcor separately"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='cube',vptable='evlavp.tab', pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgexist=[self.img+'.image',self.img+'.pb'],imgexistnot=[self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.79,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2])])
          ret2 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=0,calcres=False,calcpsf=False, specmode='cube', vptable='evlavp.tab', pbcor=True,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.79,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2]),(self.img+'.image.pbcor',1.0,[256,256,0,2])])
          self.assertTrue(self.check_final(report1+report2))

     def test_pbcor_cube_mosaicft(self):
          """ [pbcor] Test pbcor with cube with mosaicft"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='cube', vptable='evlavp.tab',pbcor=True, gridder='mosaic',pblimit=-0.2,parallel=self.parallel)

          report=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.79,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2]),(self.img+'.image.pbcor',1.0,[256,256,0,2])], tfmask=[(self.img+'.image',['T']), (self.img+'.pb',['mask0']),(self.img+'.image.pbcor',['mask0']) ] )
          self.assertTrue(self.check_final(report))

     def test_pbcor_mfs_restart(self):
          """ [pbcor] Test pbcor with mfs and a restart"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=0, specmode='mfs', vptable='evlavp.tab', pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb'], imgexistnot=[self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.7,[256,256,0,0])])

          ret2 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imgval=[(self.img+'.pb',0.7,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0])])
          self.assertTrue(self.check_final(report1+report2))

     def test_pbcor_turn_off_pbmask(self):
          """ [pbcor] Test pbcor with mfs where the internal T/F mask is turned off"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=0, specmode='mfs', vptable='evlavp.tab', pbcor=True,parallel=self.parallel)
          report1=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb'], imgval=[(self.img+'.pb',0.7,[256,256,0,0])], imgmask=[(self.img+'.pb',False,[10,10,0,0]), (self.img+'.image',False,[10,10,0,0])] )

          ret2 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", 
                        niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True, calcpsf=False, calcres=False, pblimit=-0.2,parallel=self.parallel)
          report2=self.th.checkall(imgexist=[self.img+'.image', self.img+'.pb'], imgval=[(self.img+'.pb',0.7,[256,256,0,0])] , imgmask=[(self.img+'.pb',False,[10,10,0,0]), (self.img+'.image',True,[10,10,0,0])] , tfmask=[(self.img+'.image',['T']), (self.img+'.pb',['mask0']),(self.img+'.image.pbcor',['mask0'])] )

          self.assertTrue(self.check_final(report1+report2))


#####################################################
#####################################################
#####################################################
#####################################################
### Heterogeneous array imaging
class test_hetarray_imaging(testref_base):
     '''
     Tests for all kinds of heterogeneous imaging.

     Type 1 :  Antennas of different shapes and/or sizes :  gridder='mosaic'.  
     Type 2 :  Antennas have different pointing offsets (groups of antennas and time-dependence ) :  gridder='awproject'. [ Later, via CAS-11191, for 'mosaic' too ]

     Current Test list : 
     test_het_pointing_offsets_awproject_cube :  With CAS-12617  :  Test antenna-dependent and time-dependent pointing offset corrections
     test_het_pointing_offsets_awproject_mtmfs :  With CAS-12617 :  Test antenna-dependent and time-dependent pointing offset correct

     Tests to add later : 
     test_het_pointing_offsets_mosaic_cube :   With CAS-11191  :  Test antenna-dependent and time-dependent pointing offset correct
     test_het_pointing_offsets_mosaic_mtmfs :    With CAS-11191  :  Test antenna-dependent and time-dependent pointing offset correct
     test_het_antenna_mosaic :   Test ALMA 7m+12m dataset with and without cross-baselines.  
     test_het_antenna_mosaic_simulate :  With CAS-11464 : Test model prediction for a generic het array with dish diameter modified in the ANTENNA subtable. 
     '''     
#     @unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test till awproject works with CAS-9386")
     def test_het_pointing_offsets_awproject_cube(self):
          '''
          This dataset has two groups of antennas and two timesteps, with pointing centers forming the corners of a square around the source (and MS phasecenter). 
          Cube imaging with awproject :  For all three channels, check that the source and PB are the same such that pbcorrected intensity is 1.0 Jy. 
          '''
          self.prepData('refim_heterogeneous_pointings.ms')
          msname = self.msfile
          #msname = '/home/vega/rurvashi/TestCASA/VerificationTests/PointingCorrection/simulation_pointing/sim_data.ms'
          self.baselines  ={'grp1':'0,2,4,6,8,10,12,14,16,18,20,22,24,26', 
                            'grp2':'1,3,5,7,9,11,13,15,17,19,21,23,25' }

          os.environ['ATerm_OVERSAMPLING'] = '5'
          os.environ['ATerm_CONVSIZE'] = '512'
          os.environ['PO_DEBUG'] = '0'

          ## No corrections :  usepointing=False : PB goes to location of MS field phasecenter (middle of the image)
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=self.img+'_pcorr0_uspF', niter=0, specmode='cube', nchan=3,start='1.9GHz', width='0.4GHz', interpolation='nearest', pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=False)
          report1=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          (self.img+'_pcorr0_uspF.image' ,0.40,[1024,1024,0,0]), 
                                          (self.img+'_pcorr0_uspF.image' ,0.22,[1024,1024,0,1]), 
                                          (self.img+'_pcorr0_uspF.image' ,0.09,[1024,1024,0,2]), 
                                          ## Check PB at source location
                                          ## Check that PB peak is at the expected location
                                          (self.img+'_pcorr0_uspF.pb' ,1.0,[1024,1024,0,0]), 
                                          (self.img+'_pcorr0_uspF.pb' ,1.0,[1024,1024,0,1]), 
                                          (self.img+'_pcorr0_uspF.pb' ,1.0,[1024,1024,0,2]) ],
                                          tfmask=[(self.img+'_pcorr0_uspF.image',['T'])] )

          ## Note : Test for PB location in the following tests. Pick the expected location (for a single PB) and test that the value is 1.0

          ## Average correction : usepointing=True, pointingoffsetsigdev=[2000,2000] :  PB goes to the average location of all antennas, for first timestep. PB[-100,0] 
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=self.img+'_pcorr0_uspT', niter=0, specmode='cube', nchan=3,start='1.9GHz', width='0.4GHz', interpolation='nearest', pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[2000.0,2000.0])
          report2=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          (self.img+'_pcorr0_uspT.image' ,0.40,[1024,1024,0,0]), 
                                          (self.img+'_pcorr0_uspT.image' ,0.22,[1024,1024,0,1]), 
                                          (self.img+'_pcorr0_uspT.image' ,0.09,[1024,1024,0,2]), 
                                          ## Check PB at source location
                                          (self.img+'_pcorr0_uspT.pb' ,0.65,[1024,1024,0,0]), 
                                          (self.img+'_pcorr0_uspT.pb' ,0.50,[1024,1024,0,1]), 
                                          (self.img+'_pcorr0_uspT.pb' ,0.35,[1024,1024,0,2]), 
                                          ## Check that PB peak is at the expected location
                                          (self.img+'_pcorr0_uspT.pb' ,1.0,[924,1024,0,0]) ] )   



          ## Antenna/time correction : usepointing=True, pointingoffsetsigdev=[20,2000], timerange='time1', antenna='grp1' : PB = PB[-100,+100] in the top left corner. 
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=self.img+'_pcorr1_time1_grp1', niter=0, specmode='cube', nchan=3,start='1.9GHz', width='0.4GHz', interpolation='nearest', pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,2000.0], timerange='<21:40:00.0', antenna=self.baselines['grp1']+ ' & ' )
          report3=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          (self.img+'_pcorr1_time1_grp1.image' ,0.40,[1024,1024,0,0]), 
                                          (self.img+'_pcorr1_time1_grp1.image' ,0.22,[1024,1024,0,1]), 
                                          (self.img+'_pcorr1_time1_grp1.image' ,0.09,[1024,1024,0,2]), 
                                          ## Check PB at source location
                                          (self.img+'_pcorr1_time1_grp1.pb' ,0.40,[1024,1024,0,0]), 
                                          (self.img+'_pcorr1_time1_grp1.pb' ,0.22,[1024,1024,0,1]), 
                                          (self.img+'_pcorr1_time1_grp1.pb' ,0.09,[1024,1024,0,2]), 
                                          ## Check that PB peak is at the expected location
                                          (self.img+'_pcorr1_time1_grp1.pb' ,1.0,[924,1124,0,0]) ] )   
          
          
          ## Cross baselines only : usepointing=True, pointingoffsetsigdev=[20,2000], timerange='time1', antenna='grp1 & grp2' : PB = PB[-100,0].  Note : THIS IS WRONG. 
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=self.img+'_pcorr1_time1_cross', niter=0, specmode='cube', nchan=3,start='1.9GHz', width='0.4GHz', interpolation='nearest', pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,2000.0], timerange='<21:40:00.0', antenna=self.baselines['grp1'] + ' & ' + self.baselines['grp2'])
          report4=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          (self.img+'_pcorr1_time1_cross.image' ,0.40,[1024,1024,0,0]), 
                                          (self.img+'_pcorr1_time1_cross.image' ,0.22,[1024,1024,0,1]), 
                                          (self.img+'_pcorr1_time1_cross.image' ,0.09,[1024,1024,0,2]), 
                                          ## Check PB at source location
                                          (self.img+'_pcorr1_time1_cross.pb' ,0.65,[1024,1024,0,0]), 
                                          (self.img+'_pcorr1_time1_cross.pb' ,0.50,[1024,1024,0,1]), 
                                          (self.img+'_pcorr1_time1_cross.pb' ,0.35,[1024,1024,0,2]), 
                                          ## Check that PB peak is at the expected location
                                          (self.img+'_pcorr1_time1_cross.pb' ,1.0,[924,1024,0,0]) ] )   
          report4 = report4 + "This test checks for cross-baseline PB values that are known to be incorrect. Edit these values/test once the algorithm for cross-baseline PBs is fixed.\n"
          
          ## Four corners : usepointing=True, pointingoffsetsigdev=[20,20], timerange='*', antenna='grp1,grp2' : PB = Sum of PB in all 4 corners (with no cross-terms). Flux/alpha are correct. 
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=self.img+'_pcorr2_4corners', niter=0, specmode='cube', nchan=3,start='1.9GHz', width='0.4GHz', interpolation='nearest', pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0], antenna=self.baselines['grp1']+' & ; '+self.baselines['grp2']+ ' &')
          report5=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          (self.img+'_pcorr2_4corners.image' ,0.77,[1024,1024,0,0]), 
                                          (self.img+'_pcorr2_4corners.image' ,0.42,[1024,1024,0,1]), 
                                          (self.img+'_pcorr2_4corners.image' ,0.17,[1024,1024,0,2]), 
                                          ## Check PB at source location
                                          (self.img+'_pcorr2_4corners.pb' ,0.77,[1024,1024,0,0]), 
                                          (self.img+'_pcorr2_4corners.pb' ,0.42,[1024,1024,0,1]), 
                                          (self.img+'_pcorr2_4corners.pb' ,0.17,[1024,1024,0,2]), 
                                          ## Check that PB peak is at the expected location
                                          (self.img+'_pcorr2_4corners.pb' ,0.925,[924,1124,0,0]),  
                                          (self.img+'_pcorr2_4corners.pb' ,1.0,[924,924,0,0]),   
                                          (self.img+'_pcorr2_4corners.pb' ,0.925,[1124,1124,0,0]),   
                                          (self.img+'_pcorr2_4corners.pb' ,1.0,[1124,924,0,0]) ] )   
          report5 = report5 + "This test leaves out cross-baselines. Edit later to include them, once the algorithm for cross-baseline PBs is fixed.\n"
        
          #### Note : Add a run with all antennas ONLY after the cross-baselines imaging is correct.  
          #### grp1 has 14 ants. grp2 has 13.  But, the PBs for grp1 have the peak of 1.0 whereas grp2 has 0.93.  Needs to be understood. But, image and pb values match, so flux is ok. 

          ### Set these back to the default values encoded in  code/synthesis/TransformMachines2/ATerm.h 
          os.environ['ATerm_OVERSAMPLING'] = '20'
          os.environ['ATerm_CONVSIZE'] = '2048'


          self.assertTrue(self.check_final(report1+report2+report3+report4+report5))

     ###########################

     def test_het_pointing_offsets_awproject_mtmfs(self):
          '''
          This dataset has two groups of antennas and two timesteps, with pointing centers forming the corners of a square around the source (and MS phasecenter). 
          MTMFS imaging with awproject : Check that source and PB are the same. Check that alpha is 0.0 (with conjbeams=True). 
          '''
          self.prepData('refim_heterogeneous_pointings.ms')
          msname = self.msfile
          #msname = '/home/vega/rurvashi/TestCASA/VerificationTests/PointingCorrection/simulation_pointing/sim_data.ms'

          self.baselines  ={'grp1':'0,2,4,6,8,10,12,14,16,18,20,22,24,26', 
                            'grp2':'1,3,5,7,9,11,13,15,17,19,21,23,25' }

          os.environ['ATerm_OVERSAMPLING'] = '5'
          os.environ['ATerm_CONVSIZE'] = '512'
          os.environ['PO_DEBUG'] = '0'

          ######################
          ## Top Left : Grp 1 and Time 1.  

          im_pcorr2 = 'try_mt_pcorr1_grp1_time1'  # Apply antenna grouping corrections and time-dependence corrections
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=im_pcorr2+'.std', niter=0, specmode='mfs', deconvolver='mtmfs', conjbeams=True,  pblimit=-0.01,gridder='standard',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0],   antenna=self.baselines['grp1']+' &' , timerange='<21:40:00.0')
          
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=im_pcorr2, niter=0, specmode='mfs', deconvolver='mtmfs', conjbeams=True,  pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0],   antenna=self.baselines['grp1'] +' &' , timerange='<21:40:00.0')
          
          os.system('rm -rf '+im_pcorr2+'.psf*')
          os.system('cp -r ' + im_pcorr2+'.std.psf.tt0 ' + im_pcorr2+'.psf.tt0')
          os.system('cp -r ' + im_pcorr2+'.std.psf.tt1 ' + im_pcorr2+'.psf.tt1')
          os.system('cp -r ' + im_pcorr2+'.std.psf.tt2 ' + im_pcorr2+'.psf.tt2')
          
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=im_pcorr2, niter=10, specmode='mfs', deconvolver='mtmfs', conjbeams=True,  pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0],   antenna=self.baselines['grp1']+' &'  , timerange='<21:40:00.0', calcpsf=False, calcres=False)
          
          report1=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          ('try_mt_pcorr1_grp1_time1.image.tt0' ,0.2,[1024,1024,0,0]), ## Average of 0.4,0.22,0.09 for the 3 chans.
                                          ## Check PB at source location
                                          ## Check that PB peak is at the expected location
                                          ('try_mt_pcorr1_grp1_time1.pb.tt0' ,0.2,[1024,1024,0,0]), 
                                          ## Check alpha
                                          ('try_mt_pcorr1_grp1_time1.alpha' ,0.0,[1024,1024,0,0]) ] )


          ######################
          ## 4 corners : All baselines and both timesteps ( leave out cross baselines for now ).
          
          im_pcorr2 = 'try_mt_pcorr2_4corners'  # Apply antenna grouping corrections and time-dependence corrections
         
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=im_pcorr2+'.std', niter=0, specmode='mfs', deconvolver='mtmfs', conjbeams=True,  pblimit=-0.01,gridder='standard',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0],   antenna=self.baselines['grp1'] + ' & ; ' + self.baselines['grp2'] + ' &')
         
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=im_pcorr2, niter=0, specmode='mfs', deconvolver='mtmfs', conjbeams=True,  pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0],   antenna=self.baselines['grp1'] + ' & ; ' + self.baselines['grp2'] + ' &')
         
          os.system('rm -rf '+im_pcorr2+'.psf*')
          os.system('cp -r ' + im_pcorr2+'.std.psf.tt0 ' + im_pcorr2+'.psf.tt0')
          os.system('cp -r ' + im_pcorr2+'.std.psf.tt1 ' + im_pcorr2+'.psf.tt1')
          os.system('cp -r ' + im_pcorr2+'.std.psf.tt2 ' + im_pcorr2+'.psf.tt2')
          
          tclean(vis=msname, datacolumn='observed', imsize=2048,cell=5.0, imagename=im_pcorr2, niter=10, specmode='mfs', deconvolver='mtmfs', conjbeams=True,  pblimit=-0.01,gridder='awproject',wbawp=True,psterm=False, usepointing=True, pointingoffsetsigdev=[20.0,20.0],   antenna=self.baselines['grp1'] + ' & ; ' + self.baselines['grp2'] + ' &', calcres=False, calcpsf=False)
          
          report2=self.th.checkall(imgval=[
                                          ## Check source intensity
                                          ('try_mt_pcorr2_4corners.image.tt0' ,0.38,[1024,1024,0,0]), ## Average of 0.4,0.22,0.09 for the 3 chans.
                                          ## Check PB at source location
                                          ## Check that PB peak is at the expected location
                                          ('try_mt_pcorr2_4corners.pb.tt0' ,0.38,[1024,1024,0,0]), 
                                          ## Check alpha
                                          ('try_mt_pcorr2_4corners.alpha' ,0.0,[1024,1024,0,0]) ] )


          ### Set these back to the default values encoded in  code/synthesis/TransformMachines2/ATerm.h 
          os.environ['ATerm_OVERSAMPLING'] = '20'
          os.environ['ATerm_CONVSIZE'] = '2048'

          self.assertTrue(self.check_final(report1+report2))

     ###########################

#####################################################
#####################################################
#####################################################
#####################################################
class test_mosaic_mtmfs(testref_base):
     def test_mtmfs_standard_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = ''
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          widebandpbcor(vis=self.msfile, imagename=self.img,nterms=2,reffreq='1.5GHz',pbmin=0.1, field=field, spwlist=[0,1,2], chanlist=[0,0,0], weightlist=[1.0,1.0,1.0])
          report1=self.th.checkall(imgval=[(self.img+'.pbcor.image.tt0',1.021635,[512,596,0,0]),(self.img+'.pbcor.image.alpha', -0.4811849,[512,596,0,0]),(self.img+'.pbcor.workdirectory/' + self.img + '.pb.alpha',-1.526915,[512,596,0,0])])
         
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          widebandpbcor(vis=self.msfile, imagename=self.img,nterms=2,reffreq='1.5GHz',pbmin=0.1, field=field, spwlist=[0,1,2], chanlist=[0,0,0], weightlist=[1.0,1.0,1.0])
          report2=self.th.checkall(imgval=[(self.img+'.pbcor.image.tt0',1.02163529396,[512,596,0,0]),(self.img+'.pbcor.image.alpha', -0.481181353331 ,[512,596,0,0]),(self.img+'.pbcor.workdirectory/' + self.img + '.pb.alpha',-1.52691471577,[512,596,0,0])])
          
          self.assertTrue(self.check_final(report1+report2))
     
     ###########################
     #warning
     def test_mtmfs_standard_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.696767389774,[512,596,0,0]),(self.img+'.alpha', -1.22704184055,[512,596,0,0])])
          widebandpbcor(vis=self.msfile, imagename=self.img,nterms=2,reffreq='1.5GHz',pbmin=0.1, field=field, spwlist=[0,1,2], chanlist=[0,0,0], weightlist=[1.0,1.0,1.0])
          report2=self.th.checkall(imgval=[(self.img+'.pbcor.image.tt0',0.7115,[512,596,0,0]),(self.img+'.pbcor.image.alpha',-1.2540572,[512,596,0,0]),(self.img+'.pbcor.workdirectory/' + self.img + '.pb.alpha',0.027015,[512,596,0,0])])
          
          os.system('rm -rf ' + self.img+'*')
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.696766972542 ,[512,596,0,0]),(self.img+'.alpha',-1.2270386219,[512,596,0,0])])
          widebandpbcor(vis=self.msfile, imagename=self.img,nterms=2,reffreq='1.5GHz',pbmin=0.1, field=field, spwlist=[0,1,2], chanlist=[0,0,0], weightlist=[1.0,1.0,1.0])
          report4=self.th.checkall(imgval=[(self.img+'.pbcor.image.tt0',0.7115,[512,596,0,0]),(self.img+'.pbcor.image.alpha',-1.2540572,[512,596,0,0]),(self.img+'.pbcor.workdirectory/' + self.img + '.pb.alpha',0.027015,[512,596,0,0])])
          
          self.assertTrue(self.check_final(report1 + report2 + report3 + report4 + '\n Warning: values must be theoretically validated'))
          
#############################
     #
     def test_mtmfs_mosaic_cbFalse_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.5302894,[512,596,0,0]),(self.img+'.pb.tt0',0.5144197,[512,596,0,0]),(self.img+'.alpha',-3.29580,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.450807213783 ,[512,596,0,0]),(self.img+'.pb.tt0',0.514419734478,[512,596,0,0]),(self.img+'.alpha',-1.6693893671,[512,596,0,0])])
          self.assertTrue(self.check_final(report1+report2))

     def test_mtmfs_mosaic_cbFalse_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.97402,[512,596,0,0]),(self.img+'.pb.tt0', 0.982510,[512,596,0,0]),(self.img+'.alpha', -0.9880,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.982747793198,[512,596,0,0]),(self.img+'.pb.tt0', 0.982510685921,[512,596,0,0]),(self.img+'.alpha', -0.69006639719,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2 + '\n Warning: values must be theoretically validated'))
      
     def test_mtmfs_mosaic_cbTrue_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.492861807346,[512,596,0,0]),(self.img+'.pb.tt0', 0.492951124907,[512,596,0,0]),(self.img+'.alpha',-0.586531162262,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.494806170464,[512,596,0,0]),(self.img+'.pb.tt0', 0.492951124907,[512,596,0,0]),(self.img+'.alpha',-0.548104763031,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2))
     
     def test_mtmfs_mosaic_cbTrue_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',1.004476,[512,596,0,0]),(self.img+'.pb.tt0',0.98724,[512,596,0,0]),(self.img+'.alpha',-1.216693,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.989150226116,[512,596,0,0]),(self.img+'.pb.tt0',0.987243115902,[512,596,0,0]),(self.img+'.alpha',-0.771613717079,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2 + '\n Warning: values must be theoretically validated'))
          
          ###################    
     def test_mtmfs_mosaic_cbTrue_onefield_use_standard_psf(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbTrue")
          phasecenter = '' 
          field='0'
          
          ## Run the standard gridder imaging.
          tclean(vis=self.msfile, imagename='std',niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='standard',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          
          ## Make the residual images using mosaic/awproject
          os.system('rm -rf '+self.img+'.*'); 
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
           
          ## Replace the PSFs with those from the standard gridder
          os.system('rm -rf '+self.img+'.psf*')
          os.system('cp -r std.psf.tt0 '+self.img+'.psf.tt0')
          os.system('cp -r std.psf.tt1 '+self.img+'.psf.tt1')
          os.system('cp -r std.psf.tt2 '+self.img+'.psf.tt2')
          
          ## Restart tclean with calcres=False and calcpsf=False
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,calcres=False, calcpsf=False,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)

          report=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.492902100086 ,[512,596,0,0]),(self.img+'.pb.tt0',   0.492951124907 ,[512,596,0,0]),(self.img+'.alpha',    -0.591838240623,[512,596,0,0])])
          self.assertTrue(self.check_final(report))
             
     def test_mtmfs_mosaic_cbTrue_twofield_use_standard_psf(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbTrue")
          phasecenter = '' 
          field='0,1'
          
          ## Run the standard gridder imaging.
          tclean(vis=self.msfile, imagename='std',niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='standard',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          
          ## Make the residual images using mosaic/awproject
          os.system('rm -rf '+self.img+'.*'); 
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
           
          ## Replace the PSFs with those from the standard gridder
          os.system('rm -rf '+self.img+'.psf*')
          os.system('cp -r std.psf.tt0 '+self.img+'.psf.tt0')
          os.system('cp -r std.psf.tt1 '+self.img+'.psf.tt1')
          os.system('cp -r std.psf.tt2 '+self.img+'.psf.tt2')
          
          ## Restart tclean with calcres=False and calcpsf=False
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,calcres=False, calcpsf=False,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          report=self.th.checkall(imgval=[(self.img+'.image.tt0',  0.982702195644,[512,596,0,0]),(self.img+'.pb.tt0', 0.987243115902,[512,596,0,0]),(self.img+'.alpha',  -0.55044066906,[512,596,0,0])])
          self.assertTrue(self.check_final(report))

#############################        

     def test_mtmfs_awproject_cbFalse_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter =''
          field='0'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field, cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.5269711,[512,596,0,0]),(self.img+'.pb.tt0',0.50752753,[512,596,0,0]),(self.img+'.alpha',-3.24132061,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field, cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.451119929552,[512,596,0,0]),(self.img+'.pb.tt0',0.507527530193,[512,596,0,0]),(self.img+'.alpha',-1.65221953392,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2))


     def test_mtmfs_awproject_cbFalse_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,  cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.97862583398,[512,596,0,0]),(self.img+'.pb.tt0',0.979142010212,[512,596,0,0]),(self.img+'.alpha',-1.24368548393,[512,596,0,0])])
          
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject', field=field,  cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.975049376488,[512,596,0,0]),(self.img+'.pb.tt0',0.979141950607,[512,596,0,0]),(self.img+'.alpha',-0.80353230238,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2 + '\n Warning: values must be theoretically validated'))
          
     def test_mtmfs_awproject_cbTrue_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbTrue")
          phasecenter = '' 
          field='0'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,  cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.477538466454,[512,596,0,0]),(self.img+'.pb.tt0', 0.479197412729,[512,596,0,0]),(self.img+'.alpha',  -0.562356948853,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,  cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.482177525759,[512,596,0,0]),(self.img+'.pb.tt0', 0.479197442532,[512,596,0,0]),(self.img+'.alpha',  -0.568624258041,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2))
             
     def test_mtmfs_awproject_cbTrue_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbTrue")
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,  cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.974482476711 ,[512,596,0,0]),(self.img+'.pb.tt0',0.979797422886,[512,596,0,0]),(self.img+'.alpha', -0.538577735424 ,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,  cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.97661459446 ,[512,596,0,0]),(self.img+'.pb.tt0',0.979797422886,[512,596,0,0]),(self.img+'.alpha', -0.538577854633 ,[512,596,0,0])])
          self.assertTrue(self.check_final(report1+report2))
###################    
     def test_mtmfs_awproject_cbTrue_onefield_use_standard_psf(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbTrue")
          phasecenter = '' 
          field='0'
          
          ## Run the standard gridder imaging.
          tclean(vis=self.msfile, imagename='std',niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='standard',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          
          ## Make the residual images using mosaic/awproject
          os.system('rm -rf '+self.img+'.*'); 
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
           
          ## Replace the PSFs with those from the standard gridder
          os.system('rm -rf '+self.img+'.psf*')
          os.system('cp -r std.psf.tt0 '+self.img+'.psf.tt0')
          os.system('cp -r std.psf.tt1 '+self.img+'.psf.tt1')
          os.system('cp -r std.psf.tt2 '+self.img+'.psf.tt2')
          
          ## Restart tclean with calcres=False and calcpsf=False
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,calcres=False, calcpsf=False,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)

          report=self.th.checkall(imgval=[(self.img+'.image.tt0',0.477585822344,[512,596,0,0]),(self.img+'.pb.tt0',  0.479197412729 ,[512,596,0,0]),(self.img+'.alpha',   -0.569523513317,[512,596,0,0])])
          self.assertTrue(self.check_final(report))
             
     def test_mtmfs_awproject_cbTrue_twofield_use_standard_psf(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbTrue")
          phasecenter = '' 
          field='0,1'
          
          ## Run the standard gridder imaging.
          tclean(vis=self.msfile, imagename='std',niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='standard',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          
          ## Make the residual images using mosaic/awproject
          os.system('rm -rf '+self.img+'.*'); 
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache,conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
           
          ## Replace the PSFs with those from the standard gridder
          os.system('rm -rf '+self.img+'.psf*')
          os.system('cp -r std.psf.tt0 '+self.img+'.psf.tt0')
          os.system('cp -r std.psf.tt1 '+self.img+'.psf.tt1')
          os.system('cp -r std.psf.tt2 '+self.img+'.psf.tt2')
          
          ## Restart tclean with calcres=False and calcpsf=False
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=True, wbawp=True, psterm=False,pblimit=0.1,calcres=False, calcpsf=False,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',parallel=self.parallel)
          report=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.974559485912,[512,596,0,0]),(self.img+'.pb.tt0', 0.979797422886,[512,596,0,0]),(self.img+'.alpha',  -0.542647540569,[512,596,0,0])])
          self.assertTrue(self.check_final(report))
          
          
     def test_mtmfs_mosaic_cbFalse_onefield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, usepointing=True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.5302894,[512,596,0,0]),(self.img+'.pb.tt0',0.5144197,[512,596,0,0]),(self.img+'.alpha',-3.29580,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, usepointing=True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.450807213783 ,[512,596,0,0]),(self.img+'.pb.tt0',0.514419734478,[512,596,0,0]),(self.img+'.alpha',-1.6693893671,[512,596,0,0])])
          self.assertTrue(self.check_final(report1+report2))

     def test_mtmfs_mosaic_cbFalse_twofield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, usepointing=True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.97402,[512,596,0,0]),(self.img+'.pb.tt0', 0.982510,[512,596,0,0]),(self.img+'.alpha', -0.9880,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, usepointing=True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.982747793198,[512,596,0,0]),(self.img+'.pb.tt0', 0.982510685921,[512,596,0,0]),(self.img+'.alpha', -0.69006639719,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2 + '\n Warning: values must be theoretically validated'))
          
          
          
     def test_mtmfs_awproject_cbFalse_onefield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter =''
          field='0'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field, usepointing=True, cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.5269711,[512,596,0,0]),(self.img+'.pb.tt0',0.50752753,[512,596,0,0]),(self.img+'.alpha',-3.24132061,[512,596,0,0])])
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field, usepointing=True, cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.451119929552,[512,596,0,0]),(self.img+'.pb.tt0',0.507527530193,[512,596,0,0]),(self.img+'.alpha',-1.65221953392,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2))


     def test_mtmfs_awproject_cbFalse_twofield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,  usepointing=True, cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.tt0',0.97862583398,[512,596,0,0]),(self.img+'.pb.tt0',0.979142010212,[512,596,0,0]),(self.img+'.alpha',-1.24368548393,[512,596,0,0])])
          
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject', field=field,  usepointing=True, cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel)
          report2=self.th.checkall(imgval=[(self.img+'.image.tt0',0.975049376488,[512,596,0,0]),(self.img+'.pb.tt0',0.979141950607,[512,596,0,0]),(self.img+'.alpha',-0.80353230238,[512,596,0,0])])
          self.assertTrue(self.check_final(report1 + report2 + '\n Warning: values must be theoretically validated'))
          

     def test_mtmfs_mosaic_mwFalse_briggs_twofield(self):  ## Added in CAS-13438 to catch failing case (pb os zero). 
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=5,specmode='mfs',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,deconvolver='mtmfs',nterms=2,reffreq='1.5GHz',pbcor=False,parallel=self.parallel,mosweight=False,weighting='briggs')
          report=self.th.checkall(imgval=[(self.img+'.image.tt0', 0.9794,[512,596,0,0]),(self.img+'.pb.tt0', 0.9817,[512,596,0,0]),(self.img+'.alpha', -0.797,[512,596,0,0])])

          self.assertTrue(self.check_final(report))
          
###########################################################
###########################################################
###########################################################
class test_mosaic_cube(testref_base):
     
     @unittest.skip('Skip test.') #Skip until CAS-13420 is resolved
     def test_mosaic_briggsbwtaper(self):
          self.prepData('refim_alma_mosaic.ms')
          
          tclean(vis=self.msfile,imagename=self.img+'1',imsize=[350,280],cell=[0.06,0.06 ],specmode='cube',niter=0,gridder='mosaic',phasecenter='J2000 12:01:52.430856 -18.51.49.94369',weighting='briggsbwtaper',robust=0.5 ,perchanweightdensity=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'1.image', 1,[175,140,0,0])])
          
          _ia.open(self.img+'1.image')
          nchan = _ia.shape()[3]
          briggsbwtaper_beamarea = np.zeros(nchan)
          for k in range(nchan):
            briggsbwtaper_beamarea[k]= _ia.beamarea(k,0)['arcsec2']
          _ia.close()
          
          _, report2 = self.th.check_val(briggsbwtaper_beamarea[0], 0.31155618 , valname='beam_area_chan_0', exact=False)
          _, report3 = self.th.check_val(briggsbwtaper_beamarea[1], 0.30887386 , valname='beam_area_chan_1', exact=False)
          _, report4 = self.th.check_val(briggsbwtaper_beamarea[2], 0.30639076 , valname='beam_area_chan_2', exact=False)
 
          tclean(vis=self.msfile,imagename=self.img+'2',imsize=[350,280],cell=[0.06,0.06 ],specmode='cube',niter=0,gridder='mosaic',phasecenter='J2000 12:01:52.430856 -18.51.49.94369',weighting='briggs',robust=0.5 ,perchanweightdensity=True,parallel=self.parallel)
          report5=self.th.checkall(imgval=[(self.img+'2.image', 1,[175,140,0,0])])
          
          self.assertTrue(self.th.check_beam_compare(self.img+'1.image', self.img+'2.image'))
          
          self.assertTrue(self.check_final(pstr = report1+report2+report3+report4+report5))


     def test_cube_standard_onefield(self):
          """ [mosaic_cube] Test_cube_standard_onefield : One field, standard gridder  """
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = ''
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor', 1.10922563076,[512,596,0,0]),(self.img+'.image.pbcor',0.989520609379,[512,596,0,1]),(self.img+'.image.pbcor',0.90361648798,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.505622766021, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor', 1.10922467709,[512,596,0,0]),(self.img+'.image.pbcor', 0.989521086216,[512,596,0,1]),(self.img+'.image.pbcor',0.903617084026,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.505622766021, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
     
     def test_cube_standard_twofield(self):
          """ [mosaic_cube] Test_cube_standard_twofield : two field, cube imaging with standard gridder  """
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter ='J2000 19h59m28.5 +40d40m01.5' # pointing center of field0
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',0.895364701748,[512,596,0,0]),(self.img+'.image.pbcor',0.701696813107,[512,596,0,1]),(self.img+'.image.pbcor',0.539412796497,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -1.249799253, valname='Spectral flux', exact=False)
          
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec', gridder='standard',field=field,pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor', 0.895363807678,[512,596,0,0]),(self.img+'.image.pbcor',0.701696515083,[512,596,0,1]),(self.img+'.image.pbcor',0.539413094521,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -1.24979542771, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
     
###########################################################
     def test_cube_mosaic_cbFalse_mwFalse_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.1021887064,[512,596,0,0]),(self.img+'.image.pbcor',0.983152985573,[512,596,0,1]),(self.img+'.image.pbcor',0.905956506729,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index,  -0.483567460558, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11281049252,[512,596,0,0]),(self.img+'.image.pbcor',0.994237959385,[512,596,0,1]),(self.img+'.image.pbcor',0.912590324879,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index,  -0.48920856272, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
          
     def test_cube_mosaic_cbFalse_mwFalse_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.10354316235,[512,596,0,0]),(self.img+'.image.pbcor',0.981979846954,[512,596,0,1]),(self.img+'.image.pbcor', 0.895015060902,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.516544552, valname='Spectral flux', exact=False)
          
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11373543739,[512,596,0,0]),(self.img+'.image.pbcor',0.994349777699,[512,596,0,1]),(self.img+'.image.pbcor', 0.909038066864,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.500876470767, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
           
     def test_cube_mosaic_cbFalse_mwTrue_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.1021887064,[512,596,0,0]),(self.img+'.image.pbcor',0.983152985573,[512,596,0,1]),(self.img+'.image.pbcor',0.905956506729,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index,  -0.483548182038, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11281049252,[512,596,0,0]),(self.img+'.image.pbcor',0.994237959385,[512,596,0,1]),(self.img+'.image.pbcor',0.912590324879,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index,  -0.48920856272, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))

     def test_cube_mosaic_cbFalse_mwTrue_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.10354316235,[512,596,0,0]),(self.img+'.image.pbcor', 0.981979727745,[512,596,0,1]),(self.img+'.image.pbcor', 0.895014822483,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.51654520997, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11373543739,[512,596,0,0]),(self.img+'.image.pbcor', 0.994349777699,[512,596,0,1]),(self.img+'.image.pbcor', 0.909038066864,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.50087647076, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
#####################################################  
     #@unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_cube_awproject_cbFalse_mwFalse_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter = '' 
          field='0'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor', 1.1262229681,[512,596,0,0]),(self.img+'.image.pbcor', 0.996681272984,[512,596,0,1]),(self.img+'.image.pbcor', 0.879481077194,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index,  -0.641964, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)      
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor', 1.13023257256,[512,596,0,0]),(self.img+'.image.pbcor', 1.00326132774,[512,596,0,1]),(self.img+'.image.pbcor',0.8681204319,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index,  -0.618663982179, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
          
     #@unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_cube_awproject_cbFalse_mwFalse_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11765027046,[512,596,0,0]),(self.img+'.image.pbcor',0.990656971931,[512,596,0,1]),(self.img+'.image.pbcor', 0.879846811295,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.590028509558, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.12036144733,[512,596,0,0]),(self.img+'.image.pbcor',0.994982719421,[512,596,0,1]),(self.img+'.image.pbcor', 0.889532327652,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.569002802902, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
          
     #@unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_cube_awproject_cbFalse_mwTrue_onefield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter = '' 
          field='0'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.1262229681,[512,596,0,0]),(self.img+'.image.pbcor', 0.996681272984,[512,596,0,1]),(self.img+'.image.pbcor', 0.8681204319,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.641964870168, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.13023257256,[512,596,0,0]),(self.img+'.image.pbcor',  1.00326132774,[512,596,0,1]),(self.img+'.image.pbcor', 0.879481077194,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index,  -0.61866398217, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
          
     #@unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_cube_awproject_cbFalse_mwTrue_twofield(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11765027046,[512,596,0,0]),(self.img+'.image.pbcor',0.990656971931,[512,596,0,1]),(self.img+'.image.pbcor', 0.879846811295,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.590028509558, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.12036144733,[512,596,0,0]),(self.img+'.image.pbcor',0.994982719421,[512,596,0,1]),(self.img+'.image.pbcor', 0.889532327652,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.569002802902, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4)) 
          
#####################################################  
##################################################### 
##################################################### 
     def test_cube_mosaic_cbFalse_mwFalse_onefield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,  usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.1021887064,[512,596,0,0]),(self.img+'.image.pbcor',0.983152985573,[512,596,0,1]),(self.img+'.image.pbcor',0.905956506729,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index,  -0.483567460558, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,  usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11281049252,[512,596,0,0]),(self.img+'.image.pbcor',0.994237959385,[512,596,0,1]),(self.img+'.image.pbcor',0.912590324879,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index,  -0.48920856272, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
          
     def test_cube_mosaic_cbFalse_mwFalse_twofield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,  usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.10354316235,[512,596,0,0]),(self.img+'.image.pbcor',0.981979846954,[512,596,0,1]),(self.img+'.image.pbcor', 0.895015060902,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.516544552, valname='Spectral flux', exact=False)
          
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='mosaic',field=field,  usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11373543739,[512,596,0,0]),(self.img+'.image.pbcor',0.994349777699,[512,596,0,1]),(self.img+'.image.pbcor', 0.909038066864,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.500876470767, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))

     #@unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_cube_awproject_cbFalse_mwFalse_onefield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter = '' 
          field='0'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor', 1.1262229681,[512,596,0,0]),(self.img+'.image.pbcor', 0.996681272984,[512,596,0,1]),(self.img+'.image.pbcor', 0.879481077194,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index,  -0.641964, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)      
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor', 1.13023257256,[512,596,0,0]),(self.img+'.image.pbcor', 1.00326132774,[512,596,0,1]),(self.img+'.image.pbcor',0.8681204319,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index,  -0.618663982179, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))
          
     #@unittest.skipIf(True, "The awproject gridder does not currently work with specmode='cube'.")
     def test_cube_awproject_cbFalse_mwFalse_twofield_upTrue(self):
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          self.prepCfcache("cfcache_oneshiftpoint_mosaic_cbFalse")
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.11765027046,[512,596,0,0]),(self.img+'.image.pbcor',0.990656971931,[512,596,0,1]),(self.img+'.image.pbcor', 0.879846811295,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report2 = self.th.check_val(spectral_index, -0.590028509558, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='awproject',field=field,cfcache=self.cfcache, usepointing = True, conjbeams=False, wbawp=True, psterm=False,pblimit=0.1,reffreq='1.5GHz',pbcor=True,mosweight=False,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.12036144733,[512,596,0,0]),(self.img+'.image.pbcor',0.994982719421,[512,596,0,1]),(self.img+'.image.pbcor', 0.889532327652,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          _, report4 = self.th.check_val(spectral_index, -0.569002802902, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(report1+report2+report3+report4))

######
# wprojecton tests 
     def test_cube_wproj_onefield(self):
          """ [test_mosaic_cube] test_cube_wproj_onefield : One field, widefield cube imaging with wprojection gridder """
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = ''
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='wproject',field=field, wprojplanes=4, pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.10409939289,[512,596,0,0]),(self.img+'.image.pbcor',0.983174622059,[512,596,0,1]),(self.img+'.image.pbcor',0.896172583103,[512,596,0,2])])

          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          report2 = self.th.check_val(spectral_index,  -0.51459974954, valname='Spectral flux', exact=False)

          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='wproject',field=field, wprojplanes=4, pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.10441792011,[512,596,0,0]),(self.img+'.image.pbcor',0.98375672102,[512,596,0,1]),(self.img+'.image.pbcor',0.897617280483,[512,596,0,2])])

          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          report4 = self.th.check_val(spectral_index,  -0.511338498497, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(str(report1)+str(report2)+str(report3)+str(report4)))

     def test_cube_wproj_onefield_autowprojplanes(self):
          """ [test_mosaic_cube] test_cube_wproj_onefield_autowprojplanes : One field, widefield cube imaging, gridder='wproject', automaticalluy calculate wprojplanes  """
          # reduce the imsize to 512 so for faster excution
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = '' 
          field='0'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=512, phasecenter=phasecenter,cell='10.0arcsec',gridder='wproject',field=field, wprojplanes=-1, pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.09958565,[256,340,0,0]),(self.img+'.image.pbcor',0.97893494,[256,340,0,1]),(self.img+'.image.pbcor',0.89169014,[256,340,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[256,340,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[256,340,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          report2 = self.th.check_val(spectral_index,  -0.51459974954, valname='Spectral flux', exact=False)
          
          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=512, phasecenter=phasecenter,cell='10.0arcsec',gridder='wproject',field=field, wprojplanes=-1, pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',1.10441792011,[256,340,0,0]),(self.img+'.image.pbcor',0.98375672102,[256,340,0,1]),(self.img+'.image.pbcor',0.897617280483,[256,340,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[256,340,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[256,340,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          report4 = self.th.check_val(spectral_index,  -0.511338498497, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(str(report1)+str(report2)+str(report3)+str(report4)))
          
     def test_cube_wproj_twofield(self):
          """ [test_mosaic_cube] test_cube_wproj_twofield : Two fields, widefield cube imaging with wprojection """
          self.prepData('refim_oneshiftpoint.mosaic.ms')
          phasecenter = 'J2000 19h59m28.5 +40d40m01.5' # pointing center of field0 
          field='0,1'
          tclean(vis=self.msfile, imagename=self.img,niter=0,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='wproject',field=field, wprojplanes=4,  pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
          report1=self.th.checkall(imgval=[(self.img+'.image.pbcor',0.894496798515,[512,596,0,0]),(self.img+'.image.pbcor',0.701376080513,[512,596,0,1]),(self.img+'.image.pbcor',0.539906442165,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          report2 = self.th.check_val(spectral_index,-1.24515141863, valname='Spectral flux', exact=False)

          tclean(vis=self.msfile, imagename=self.img,niter=10,specmode='cube',spw='*',imsize=1024, phasecenter=phasecenter,cell='10.0arcsec',gridder='wproject',field=field, wprojplanes=4,pblimit=0.1,reffreq='1.5GHz',pbcor=True,parallel=self.parallel)
 
          report3=self.th.checkall(imgval=[(self.img+'.image.pbcor',0.894754827023,[512,596,0,0]),(self.img+'.image.pbcor',0.701884686947,[512,596,0,1]),(self.img+'.image.pbcor',0.540905594826,[512,596,0,2])])
          
          source_flux_v0 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,0])
          source_flux_v2 = self.th.get_pix(self.img+'.image.pbcor',[512,596,0,2])
          v0 = 1.2 #In GHz
          v2 = 1.8 #In GHz
          spectral_index = np.log(source_flux_v0/source_flux_v2)/np.log(v0/v2)
          report4 = self.th.check_val(spectral_index, -1.24130281995, valname='Spectral flux', exact=False)
          self.assertTrue(self.check_final(str(report1)+str(report2)+str(report3)+str(report4)))

###########################################################
###########################################################
###########################################################
class test_ephemeris(testref_base):

     def test_onefield_mfs_eph(self):
          " [ephemeris] test_onefield_mfs_eph : single field (standard gridder), mfs mode "

          self.prepData('venus_ephem_test.ms')
          ret = tclean(vis=self.msfile, field='0', imagename=self.img, imsize=[288, 288], cell=['0.14arcsec'], phasecenter='TRACKFIELD', specmode='mfs', gridder='standard', niter=0, interactive=0, parallel=self.parallel)

          # Retrieve original image and test image statistics
          _ia.open(refdatapath+'venus_sf_ephem_test.residual')
          orig_stats = _ia.statistics()
          orig_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          _ia.open(self.img+'.residual')
          test_stats = _ia.statistics()
          test_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          # Determine metrics for testing
          # Check 1: tests flux stays within 1% of original image
          if (test_stats['sum'] - orig_stats['sum'])/orig_stats['sum'] < 0.01:
               result = True
          else:
               result = False
          _, report1 = self.th.check_val(result, True, valname='Flux within 1% of original', exact=True)

          # Check 2: tests positions shifts stays within 1% of original image
          if np.sum(np.absolute(test_freqavg - orig_freqavg)) / np.sum(np.absolute(orig_freqavg)) < 0.01:
               result = True
          else:
               result = False
          _, report2 = self.th.check_val(result, True, valname='Position shift within 1% of original', exact=True)

          # Check 3: tests position shifts are less than 10% of angular resolution; distance in pixels multiplied by cell size in arcsecs; PSF beam width calculated using lambda/max_baseline
          psf_beam_width = 1.176
          distance = np.sqrt((orig_stats['maxpos'][0] - test_stats['maxpos'][0])**2 + (orig_stats['maxpos'][1] - test_stats['maxpos'][1])**2)*0.14
          if distance/psf_beam_width < 0.1:
               result = True
          else:
               result = False
          _, report3 = self.th.check_val(result, True, valname='Position shift lass than 10% of angular resolution', exact=True)

          report = report1 + report2 + report3
          self.assertTrue(self.check_final(pstr=report))


     def test_onefield_cube_eph(self):
          " [ephemeris] test_onefield_cube_eph : single field (standard gridder), cubesource mode "

          self.prepData('venus_ephem_test.ms')
          ret = tclean(vis=self.msfile, field='0', imagename=self.img, imsize=[288, 288], cell=['0.14arcsec'], phasecenter='TRACKFIELD', specmode='cubesource', gridder='standard', niter=0, interactive=0, parallel=False)

          # Retrieve original image and test image statistics
          _ia.open(refdatapath+'venus_sf_ephem_test.residual')
          orig_stats = _ia.statistics()
          orig_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          _ia.open(self.img+'.residual')
          test_stats = _ia.statistics()
          test_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          # Determine metrics for testing
          # Check 1: tests flux stays within 1% of original image
          if (test_stats['sum'] - orig_stats['sum'])/orig_stats['sum'] < 0.01:
               result = True
          else:
               result = False
          _, report1 = self.th.check_val(result, True, valname='Flux within 1% of original', exact=True)

          # Check 2: tests positions shifts stays within 1% of original image
          if np.sum(np.absolute(test_freqavg - orig_freqavg)) / np.sum(np.absolute(orig_freqavg)) < 0.01:
               result = True
          else:
               result = False
          _, report2 = self.th.check_val(result, True, valname='Position shift within 1% of original', exact=True)

          # Check 3: tests position shifts are less than 10% of angular resolution; distance in pixels multiplied by cell size in arcsecs; PSF beam width calculated using lambda/max_baseline
          psf_beam_width = 1.176
          distance = np.sqrt((orig_stats['maxpos'][0] - test_stats['maxpos'][0])**2 + (orig_stats['maxpos'][1] - test_stats['maxpos'][1])**2)*0.14
          if distance/psf_beam_width < 0.1:
               result = True
          else:
               result = False
          _, report3 = self.th.check_val(result, True, valname='Position shift lass than 10% of angular resolution', exact=True)

          report = report1 + report2 + report3
          self.assertTrue(self.check_final(pstr=report))


     def test_multifield_mfs_eph(self):
          " [ephemeris] test_multifield_mfs_eph : multifield (mosaic gridder), mfs mode "

          self.prepData('venus_ephem_test.ms')
          ret = tclean(vis=self.msfile, imagename=self.img, imsize=[480, 420], cell=['0.14arcsec'], phasecenter='TRACKFIELD', specmode='mfs', gridder='mosaic', niter=0, interactive=0, parallel=self.parallel)

          # Retrieve original image and test image statistics
          _ia.open(refdatapath+'venus_mos_ephem_test.residual')
          orig_stats = _ia.statistics()
          orig_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          _ia.open(self.img+'.residual')
          test_stats = _ia.statistics()
          test_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          # Determine metrics for testing
          # Check 1: tests flux stays within 1% of original image
          if (test_stats['sum'] - orig_stats['sum'])/orig_stats['sum'] < 0.01:
               result = True
          else:
               result = False
          _, report1 = self.th.check_val(result, True, valname='Flux within 1% of original', exact=True)

          # Check 2: tests positions shifts stays within 1% of original image
          if np.sum(np.absolute(test_freqavg - orig_freqavg)) / np.sum(np.absolute(orig_freqavg)) < 0.01:
               result = True
          else:
               result = False
          _, report2 = self.th.check_val(result, True, valname='Position shift within 1% of original', exact=True)

          # Check 3: tests position shifts are less than 10% of angular resolution; distance in pixels multiplied by cell size in arcsecs; PSF beam width calculated using lambda/max_baseline
          psf_beam_width = 1.176
          distance = np.sqrt((orig_stats['maxpos'][0] - test_stats['maxpos'][0])**2 + (orig_stats['maxpos'][1] - test_stats['maxpos'][1])**2)*0.14
          if distance/psf_beam_width < 0.1:
               result = True
          else:
               result = False
          _, report3 = self.th.check_val(result, True, valname='Position shift lass than 10% of angular resolution', exact=True)

          report = report1 + report2 + report3
          self.assertTrue(self.check_final(pstr=report))


     def test_multifield_cube_eph(self):
          " [ephemeris] test_multifield_cube_eph : multifield (mosaic gridder), cubesource mode "

          self.prepData('venus_ephem_test.ms')
          ret = tclean(vis=self.msfile, imagename=self.img, imsize=[480, 420], cell=['0.14arcsec'], phasecenter='TRACKFIELD', specmode='cubesource', gridder='mosaic', niter=0, interactive=0, parallel=False)

          # Retrieve original image and test image statistics
          _ia.open(refdatapath+'venus_mos_ephem_test.residual')
          orig_stats = _ia.statistics()
          orig_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          _ia.open(self.img+'.residual')
          test_stats = _ia.statistics()
          test_freqavg = _ia.statistics(axes=[2])['sum']
          _ia.close()

          # Determine metrics for testing
          # Check 1: tests flux stays within 1% of original image
          if (test_stats['sum'] - orig_stats['sum'])/orig_stats['sum'] < 0.01:
               result = True
          else:
               result = False
          _, report1 = self.th.check_val(result, True, valname='Flux within 1% of original', exact=True)

          # Check 2: tests positions shifts stays within 1% of original image
          if np.sum(np.absolute(test_freqavg - orig_freqavg)) / np.sum(np.absolute(orig_freqavg)) < 0.01:
               result = True
          else:
               result = False
          _, report2 = self.th.check_val(result, True, valname='Position shift within 1% of original', exact=True)

          # Check 3: tests position shifts are less than 10% of angular resolution; distance in pixels multiplied by cell size in arcsecs; PSF beam width calculated using lambda/max_baseline
          psf_beam_width = 1.176
          distance = np.sqrt((orig_stats['maxpos'][0] - test_stats['maxpos'][0])**2 + (orig_stats['maxpos'][1] - test_stats['maxpos'][1])**2)*0.14
          if distance/psf_beam_width < 0.1:
               result = True
          else:
               result = False
          _, report3 = self.th.check_val(result, True, valname='Position shift lass than 10% of angular resolution', exact=True)

          report = report1 + report2 + report3
          self.assertTrue(self.check_final(pstr=report))


class test_errors_failures(testref_base):

     def test_wrong_spw_select_data(self):
          """
          test_wrong_spw: should produce exception in selectData.

          Should throw something like
           RuntimeError: Error in selectData() : Spw Expression: No match found for 33,
          """
          self.prepData('refim_twochan.ms')

          with self.assertRaises(RuntimeError):
               ret = tclean(vis=self.msfile, imagename=self.img+'wrong_spw', imsize=98,
                            cell='8.0arcsec', niter=2, specmode='cube', deconvolver='hogbom',
                            nchan=10, restfreq=['1.25GHz'],
                            field='0', phasecenter="J2000 19:59:28.500 +40.44.01.50",
                            spw='33:5~19', start=0, width=1,
                            veltype='radio', outframe='TOPO',
                            parallel=self.parallel)

     def test_wrong_field_select_data(self):
          """
          test_wrong_field: should produce exception in selectData.

          Exception like:
            RuntimeError: Error in selectData() : Field Expression: Partial or no match for Field ID list [33]
          """
          self.prepData('refim_twochan.ms')

          with self.assertRaises(RuntimeError):
               ret = tclean(vis=self.msfile, imagename=self.img+'wrong', imsize=98,
                            cell='8.0arcsec', niter=2, specmode='cube', deconvolver='hogbom',
                            nchan=10, restfreq=['1.25GHz'],
                            field='33', phasecenter="J2000 19:59:28.500 +40.44.01.50",
                            spw='0', start=0, width=1,
                            veltype='radio', outframe='TOPO',
                            parallel=self.parallel)

     def test_bad_freqframe_define_image(self):
          """
          test_bad_freqframe: produce exception in initializeImagers/defineImage

          File ".../__casac__/synthesisimager.py", line 180, in defineimage
            return _synthesisimager.synthesisimager_defineimage(self, *args, **kwargs)
          RuntimeError: Invalid Image Parameter set : Invalid Frequency Frame fail_!
          """
          self.prepData('refim_twochan.ms')

          with self.assertRaises(RuntimeError):
               ret = tclean(vis=self.msfile, imagename=self.img+'wrong', imsize=98,
                            cell='8.0arcsec', niter=2, specmode='cube', deconvolver='hogbom',
                            nchan=10, restfreq=['1.25GHz'],
                            field='0', phasecenter="J2000 02:59:28.500 -40.44.01.50",
                            spw='0', start=22, width=1,
                            veltype='radio', outframe='fail_!',
                            parallel=self.parallel)

     def test_error_gridding(self):
          """ test_error_gridding: produces exception in gridding

          One way to make gridding fail is to give imsize=3,4, etc. This fails like:
           RuntimeError: Error in making PSF : One or more  of the cube section failed in de/gridding. Return values for the sections: [0]
          With a preceding SEVERE error messages:

          SEVERE  task_tclean::FTMachine::initMaps  number of pixels 6 on x axis is smaller that the gridding support 6 Please use a larger value
          SEVERE  task_tclean::CubeMajorCycleAlgorithm::task (file .../code/synthesis/ImagerObjects/CubeMajorCycleAlgorithm.cc, line 136) Exception: Error in making PSF :
          SEVERE  task_tclean::FTMachine::initMapsnumber of pixels 6 on x axis is smaller that the gridding support 6 Please use a larger value
          """
          self.prepData('refim_twochan.ms')

          with self.assertRaises(RuntimeError):
               ret = tclean(vis=self.msfile, imagename=self.img+'wrong', imsize=4,
                            cell='8.0arcsec', niter=2, specmode='cube', deconvolver='hogbom',
                            nchan=10, restfreq=['1.25GHz'],
                            field='0', phasecenter="J2000 19:59:28.500 +40.44.01.50",
                            spw='0', start=22, width=1,
                            veltype='radio', outframe='LSRK',
                            parallel=self.parallel)


if is_CASA6:
     if __name__ == '__main__':
          unittest.main()

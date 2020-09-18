##########################################################################
##########################################################################
#
# Test programs for the minor cycle deconvolve task :  test_req_task_deconvolve
#
# Each of the following categories (classes) has a set of tests within it.
#
#  test_onefield                 # basic tests, deconvolution algorithms
#  test_iterbot                  # iteration control options for mfs and cube
#  test_multifield               # multiple fields of same type and with different shapes/deconvolvers/gridders
#  test_stokes                   # multiple stokes planes, imaging with flagged correlations..
#  test_cube                     # all things cube. Spectral frame setup, handling empty channels, etc
#  test_widefield                # facets, wprojection, imagemosaic, mosaicft, awproject
#  test_mask                     # input mask options : regridding, mask file, automasking, etc
#  test_modelvis                 # saving models (column/otf), using starting models, predict-only (setjy)
#  test_ephemeris                # ephemeris tests for gridder standard and mosaic, mode mfs and cubesource
#
# To run from within casapy :  
#
#  runUnitTest.main(['test_req_task_deconvolve'])                                              # Run all tests
#  runUnitTest.main(['test_req_task_deconvolve[test_onefield]'])                               # Run tests from test_onefield
#  runUnitTest.main(['test_req_task_deconvolve[test_onefield_mtmfs]'])                         # Run one specific test
#  runUnitTest.main(['test_req_task_deconvolve[test_onefield_mtmfs,test_onefield_hogbom]'])    # Multiple specific tests
#
# To see the full list of tests :   grep "\"\"\" \[" test_req_task_deconvolve.py
#
#  These tests need data stored in data/regression/unittest/clean/refimager
#
#  For a developer build, to get the datasets locally 
#
#  --- Get the basic data repo :  svn co https://svn.cv.nrao.edu/svn/casa-data/distro data
#  --- Make directories : mkdir -p data/regression/unittest/clean; cd data/regression/unittest/clean
#  --- Get test datasets :  svn co https://svn.cv.nrao.edu/svn/casa-data/trunk/regression/unittest/clean/refimager
#
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
#  task_tclean is used to generate the input tables to task_deconvolve from these data, then
#  the test is run on task_deconvolve. It is probable, therefore, that if task_tclean tests
#  start failing, then task_deconvolve tests will also start failing.
#
##########################################################################
#
#  Tests
#
#Onefield tests (mimicing tclean)
#1. Onefield Multiscale: Should produce the same results as the test by the same name for tclean.
#imsize:200, cell:'8.0arcsec', deconvolver:'multiscale', scales:[0,20,40,100], niter:10
#testname: test_onefield_multiscale
#
#2. Onefield mfmfs: Should produce the same results as the test by the same name for tclean.
#imsize:100, cell:'8.0arcsec', deconvolver:'mtmfs', niter:10
#testname: test_onefield_mtmfs
#
#3. onefield_rectangular_pixels(self): Should produce the same results as the test by the same name for tclean.
#imsize:100, cell:['10.0arcsec','30.0arcsec'], niter:10
#testname: test_onefield_rectangular_pixels
#
#
#
#Iterbot tests (mimicing tclean)
#4. Iterbot Clark mfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'clark',gain:0.15, niter=20
#testname: test_iterbot_mfs_4
#
#5. Iterbot Clark mfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'clark',threshold:'0.5Jy',gain:0.15, niter=10, threshold='0.5Jy'
#testname: test_iterbot_mfs_5
#
#
#
#Stokes tests (mimicing tclean)
#6. Stokes I mfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',stokes:'I', niter=10
#testname: test_stokes_mfs_I
#
#7. Stokes IQUV mtmfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',stokes:'IQUV',deconvolver:'mtmfs',nterms:2, niter=10
#testname: test_stokes_mtmfs_IQUV
#
#
#
#Cube tests (mimicing tclean)
#8. Cube: Should produce the same results as the test by the same name for tclean.
#field:'0', imsize:100, cell:'8.0arcsec', specmode:'cube', nchan:10, restfreq:['1.25GHz'], phasecenter:"J2000 19:59:28.500 +40.44.01.50", deconvolver:'hogbom', spw:0, start:0, width:1, veltype:'radio', outframe:'LSRK', interpolation:'linear', niter=10
#testname: test_cube_0
#
#9. Cube, 'chanchunks' Auto: Should produce the same results as the test by the same name for tclean.
#specmode:'cube',imsize:100,cell:'10.0arcsec',deconvolver:'hogbom',chanchunks:-1, niter=10, deconvolver='hogbom'
#testname: test_cube_chanchunks_auto
#
#
#
#Masking tests (mimicing tclean)
#10. User Mask: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',usemask:'user',mask:self.img+'.mask.txt', niter=10
#deconvolve runs: mask='tst.mask.txt', mask=mstr
#testname: test_mask_1
#
#11. User Mask: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',specmode:'cube',interactive:0,usemask:'user', niter=10
#two mask runs: mask:'tst.mask.txt', mask:mstr
#testname: test_mask_2
#
#12. Missing Mask File: tst.mask is sometimes required
#mask='tst.model.txt'
#rm 'tst.model.txt', mask='tst.model.txt'
#testname: test_mask_missingfile
#
#13. Auto Mask: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',interactive:0,usemask:'auto-multithresh',niter=10
#testname: test_mask_autobox_multithresh
#
#
#
#Show multiple executions of deconvolve get a probably-correct answer
#14. Multiple Deconvolve Runs: execute deconvolve once and compare the results to those of running deconvolve twice in a row.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',threshold:'1mJy'
#first run to get expected value: niter=20
#second/third runs to get actual value: niter=10, niter=10
#testname: test_multirun_hoghog
#
#15. Deconvolve Multiscale + Hogbom: Tests the example use case of using hogbom to speed up cleaning after the benefits of multiscale have worn off.
#Note: we only test for task completion, don't know what the value should be at the end. (TODO needs validation)
#imsize:200, cell:'8.0arcsec', deconvolver:'multiscale', scales:[0,20,40,100],niter=10
#imsize:200, cell:'8.0arcsec', deconvolver:'hogbom', niter=10
#testname: test_multirun_multiscalehog
#
#16. Run then Restore: Tests the example use case of using task_deconvolve for just the restoration step.
#imsize:100, cell:['10.0arcsec','30.0arcsec'], restoration=False, niter=10
#imsize:100, cell:['10.0arcsec','30.0arcsec'], restoration=True, niter=0
#testname: test_multirun_norestore_restore
#
#
#
#'img' value input checking
#17. Missing Table: tst.residual is always required
#rm 'tst.residual', niter=10
#testname: test_imgval_missingimgs_residual
#
#18. Missing Table: tst.psf is always required
#rm 'tst.psf', niter=10
#testname: test_imgval_missingimgs_psf
#
#19. Missing Table: tst.model is used to continue deconvolution, but is not required.
#rm 'tst.model', niter=10
#testname: test_imgval_missingimgs_model
#
#20. Missing Table: tst.sumwt is never required
#rm 'tst.sumwt', niter=10
#testname: test_imgval_missingimgs_sumwt
#
#21. Reorder Image Axes: tst.residual must have the axes as is given in tclean
#imtrans("tst.residual", order="3012"), niter=10
#testname: test_imgval_axesmismatch_residual
#
#22. Reorder Image Axes: tst.psf must have the axes as is given in tclean
#imtrans("tst.psf", order="3012"), niter=10
#testname: test_imgval_axesmismatch_psf
#
#23. Reorder Image Axes: tst.model must have the axes as is given in tclean
#imtrans("tst.model", order="3012"), niter=10
#testname: test_imgval_axesmismatch_model
#
#24. Reorder Image Axes: tst.pb must have the axes as is given in tclean
#imtrans("tst.pb", order="3012"), niter=10
#testname: test_imgval_axesmismatch_pb
#
#25. Rebin To Smaller Image: everything else must have the same shape as tst.residual
#imrebin("tst.residual", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_residual
#
#26. Rebin To Smaller Image: tst.psf must have the same shape as tst.residual
#imrebin("tst.psf", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_psf
#
#27. Rebin To Smaller Image: tst.model must have the same shape as tst.residual
#imrebin("tst.model", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_model
#
#28. Rebin To Smaller Image: tst.pb must have the same shape as tst.residual
#Note: cpp code does not mind that the .pb image has a weird shape
#imrebin("tst.pb", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_pb
#
#29. Empty 'startmodel' Parameter String: Deconvolve should ignore all empty strings entered for the startmodel
#startmodel='', startmodel=['', '', ''], startmodel=['', '', self.mname2, '', '']
#testname: test_imgval_startmodel_empty
#
#30. Parameter 'startmodel' Does Not Exist: Throws an error if startmodel is set but does not exist
#startmodel='doesnotexists.model'
#testname: test_imgval_startmodel_dne
#
#31. Both 'startmodel' And tst.model Exist: Throws an error if startmodel is set and tst.model exists (must be one or the other, not both)
#startmodel='tst_2.model'
#testname: test_imgval_startmodel_model_exists
#
#32. Parameter 'startmodel' Set: Tests ability of deconvolve to copy startmodel to tst.model before starting deconvolution
#startmodel='tst_2.model'
#testname: test_imgval_startmodel_basic_copy
#
#33. Reordered Startmodel Axes: Try to deconvolve with mismatched axes between startmodel and psf (should fail).
#startmodel='tst_2.model', imtrans(order="3012")
#testname:test_imgval_startmodel_axesmismatch
#
#34. Coordinate System Change: Task deconvolve should regrid the csys of the startmodel to that of tst.residual
#set_crval0(51), niter=10
#testname:test_imgval_startmodel_csysmismatch
#
#35. Image Shape Change: Task deconvolve should regrid the shape of the startmodel to that of tst.residual
#imrebin(factor=[2,2]), niter=10
#testname:test_imgval_startmodel_shapemismatch
#
#
#
#Multiple deconvolves update the .residual
#36. Hogbom Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':hogbom, 'niter':10
#testname: test_residual_update_hogbom
#
#37. Clark Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':clark, 'niter':10
#testname: test_residual_update_clark
#
#38. Clarkstokes Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':clarkstokes, 'niter':10
#testname: test_residual_update_clarkstokes
#
#39. Multiscale Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':multiscale, 'niter':10
#testname: test_residual_update_multiscale
#
#40. MTMFS Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':mtmfs, 'niter':10
#testname: test_residual_update_mtmfs
#
#41. Mem Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':mem, 'niter':10
#testname: test_residual_update_mem
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
from casatestutils.imagerhelpers import TestHelpers 

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import quanta, measures, image, vpmanager, calibrater
    from casatasks import casalog, deconvolve, tclean, imtrans, imrebin, imregrid
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    # sys.path.append(os.path.abspath(os.path.basename(__file__)))
    # from testhelper_imager import TestHelpers

    _ia = image( )
    _vp = vpmanager( )
    _cb = calibrater( )
    _qa = quanta( )
    _me = measures( )

    # refdatapath = ctsys.resolve('regression/unittest/clean/refimager/')
    refdatapath = "/export/home/figs/bbean/casadata/master/regression/unittest/clean/refimager/"
else:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from parallel.parallel_task_helper import ParallelTaskHelper
    from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
    # from imagerhelpers.testhelper_imager import TestHelpers

    _ia = iatool( )
    _vp = vptool( )
    _cb = cbtool( )
    # not local tools
    _qa = qa
    _me = me

    refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/clean/refimager/'
    #refdatapath = "/export/home/riya/rurvashi/Work/ImagerRefactor/Runs/UnitData"
    #refdatapath = "/home/vega/rurvashi/TestCASA/ImagerRefactor/Runs/WFtests"

th = TestHelpers()

## List to be run
def suite():
    return [test_onefield, test_iterbot, test_multifield, test_stokes, test_cube, test_mask, test_multirun, test_imgval]
 
## Base Test class with Utility functions
class testref_base(unittest.TestCase):
    cachedir = "test_req_task_deconvolve_cache"

    @classmethod
    def setUpClass(cls):
        casalog.setlogfile('testlog.log')

    def setUp(self):
        self.epsilon = 0.05
        self.msfile = ""
        self.imgsrc = ""
        self.cachedir = testref_base.cachedir
        self.inptbls = [".residual", ".psf", ".pb", ".model"]
        self.img = "tst"
        # To use subdir in the output image names in some tests (CAS-10937)
        self.img_subdir = 'refimager_tst_subdir'
        self.parallel = False
        self.nnode = 0
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True
        self.loglen = len(open('testlog.log').read())

        # self.th = TestHelpers()

    def tearDown(self):
        # Default: delete all (input and output data)
        # self.delData()
        # leave for input and output (e.g. for debugging)
        self.delData(delinput=False, deloutput=False)

    @classmethod
    def tearDownClass(cls):
        # if os.path.exists(cls.cachedir):
        #     os.system('rm -rf '+cls.cachedir)
        pass

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self,imgsrc="",cache=False,stopearly=False,tclean_args={},onlycopyms=False,delold=True):
        """
        Copies the .ms file to the local directory and
        prepares the tclean(niter=0) tst.* tables for use with deconvolve.
        - cache is mostly useful when running the same test repeatedly during development.
        - onlycopyms is mostly useful when comparing the deconvolve test against the original tclean test method
        - delold is mostly useful when creating mask files, so that the newly created mask isn't removed

        1 stop if onlycopyms is True, otherwise:
        2 if test_req_task_deconvolve_residuals/imgsrc.residual exists, skip to (6)
        3 execute tclean with the given tclean_args
        4 stop if cache=False, otherwise:
        5 copy all tst.* tables to test_req_task_deconvolve_residuals/ and remove all tst.* tables
        6 copy all test_req_task_deconvolve_residuals/imgsrc.* to the working directory
        """
        if imgsrc.endswith(".ms"):
            self.msfile = imgsrc
            imgsrc = imgsrc[:-3]
        if imgsrc != "":
            self.imgsrc = imgsrc
        if delold:
            self.delData()
        if stopearly:
            return
        if not os.path.exists(self.cachedir):
            os.system('mkdir ' + self.cachedir)
        print(refdatapath)

        # set some default values for tclean
        defs = { 'vis'          : self.imgsrc+'.ms',
                 'imagename'    : self.img,
                 'niter'        : 0,
                 'restoration'  : False,
                 'calcres'      : True,
                 'pbcor'        : True,
                 'parallel'     : self.parallel }
        for k in defs.keys():
            tclean_args[k] = defs[k] if k not in tclean_args else tclean_args[k]
        
        if (onlycopyms):
            shutil.copytree(os.path.join(refdatapath,self.imgsrc+'.ms'), self.imgsrc+'.ms')
        elif (cache):
            resfiles = [os.path.join(self.cachedir, self.imgsrc + '.residual')]
            resfiles.append(resfiles[0] + '.tt0')
            print(resfiles)

            # create the residuals and other feeder tables, as necessary
            if not os.path.exists(resfiles[0]) and not os.path.exists(resfiles[1]):
                shutil.copytree(os.path.join(refdatapath,self.imgsrc+'.ms'), self.imgsrc+'.ms')
                # return
                tclean(**tclean_args)
                for tbl in self.inptbls:
                    srctbl = self.img + tbl
                    dsttbl = os.path.join(self.cachedir, self.imgsrc + tbl)
                    if os.path.exists(srctbl):
                        shutil.copytree(srctbl, dsttbl)
                os.system('rm -rf ' + self.imgsrc + '.*')
                os.system('rm -rf ' + self.img + '.*')

             # copy necessary files to the working directory
            for tbl in self.inptbls:
                if os.path.exists(os.path.join(self.cachedir, self.imgsrc+tbl)):
                    shutil.copytree(os.path.join(self.cachedir, self.imgsrc+tbl), self.img+tbl)
        else:
            print(self.imgsrc)

            # create the residuals and other feeder tables
            shutil.copytree(os.path.join(refdatapath,self.imgsrc+'.ms'), self.imgsrc+'.ms')
            tclean(**tclean_args)
            os.system('rm -rf ' + self.imgsrc + '*')

    def getLogStr(self):
        # return only the part of the log that applies to the current test method
        logstr = open('testlog.log').read()
        return logstr[self.loglen:]

    def delData(self,msname="",delinput=True,deloutput=True):
        if msname != "":
            self.msfile=msname
        if (os.path.exists(self.msfile) and delinput):
            os.system('rm -rf ' + self.msfile)
        if (os.path.exists('usermask.mask.txt') and delinput):
            os.system('rm -rf usermask.mask.txt')
        if deloutput:
            os.system('rm -rf ' + self.img_subdir)
            os.system('rm -rf ' + self.img+'*')

    def checkfinal(self,pstr=""):
        #pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  casa -c `echo $CASAPATH | awk '{print $1}'`/gcwrap/python/scripts/regressions/admin/runUnitTest.py test_refimager["+ inspect.stack()[1][3] +"]"
        pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  runUnitTest.main(['test_req_task_deconvolve["+ inspect.stack()[1][3] +"]'])"
        casalog.post(pstr,'INFO')
        if( pstr.count("(Fail") > 0 ):
            self.fail("\n"+pstr)

##############################################
##############################################

##Task level tests : one field, 2chan.
class test_onefield(testref_base):

    # Test 1
    def test_onefield_multiscale(self):
        """ [onefield] test_onefield_multiscale """
        ######################################################################################
        # Test mfs with multiscale minor cycle. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_eptwochan.ms', tclean_args={'imsize':200, 'cell':'8.0arcsec', 'deconvolver':'multiscale', 'scales':[0,20,40,100]})
        results = deconvolve(imagename=self.img, niter=10, deconvolver='multiscale', scales=[0,20,40,100], interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.823, modflux=3.816, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'],
                           imgval=[(self.img+'.model',0.234,[94,107,0,0])])
        self.checkfinal(pstr=report)

    # Test 2
    def test_onefield_mtmfs(self):
        """ [onefield] test_onefield_mtmfs """
        ######################################################################################
        # Test mt-mfs with minor cycle iterations . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'mtmfs'})
        results = deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs', interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.392, modflux=0.732, iterdone=10, imgexist=[self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0', self.img+'.model.tt0',self.img+'.model.tt1',self.img+'.alpha'],
                           imgval=[(self.img+'.model.tt0',0.733,[50,50,0,0]),(self.img+'.image.tt1',0.019,[2,94,0,0])])
        ## iterdone=11 only because of the return (iterdone_p+1) in MultiTermMatrixCleaner::mtclean() !
        self.checkfinal(pstr=report)

    # Test 3
    def test_onefield_rectangular_pixels(self):
        """ [onefield] test_onefield_rectangular_pixels """
        ######################################################################################
        # Test restoration with rectangular pixels (cas-7171). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point.ms', tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec']})
        deconvolve(imagename=self.img, niter=10)
        report=th.checkall(imgexist=[self.img+'.image'], imgval=[(self.img+'.image',0.482,[50,49,0,0])] )
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : iteration controls
class test_iterbot(testref_base):

    # Test 4
    def test_iterbot_mfs_4(self):
        """ [iterbot] test_iterbot_mfs_4 """
        ######################################################################################
        # Test Iterations with high gain . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'clark','gain':0.15})
        results = deconvolve(imagename=self.img, deconvolver='clark', niter=20, gain=0.15,interactive=0)
        report=th.checkall(ret=results['retrec'], iterdone=14, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                           imgval=[(self.img+'.model',0.937,[50,50,0,0])])

        self.checkfinal(report)

    # Test 5
    def test_iterbot_mfs_5(self):
        """ [iterbot] test_iterbot_mfs_5 """
        ######################################################################################
        # Threshold test . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'clark','threshold':'0.5Jy','gain':0.15})
        results = deconvolve(imagename=self.img, deconvolver='clark', niter=10, threshold='0.5Jy', gain=0.15, interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.499, modflux=0.626, iterdone=5, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                           imgval=[(self.img+'.model',0.626,[50,50,0,0])])

        self.checkfinal(report)

##############################################
##############################################

##Task level tests : multi-field, 2chan.
class test_multifield(testref_base):
    # n/a, outlier tests are for running a combined major cycle on many images
    # just run deconvolve on each of the individual images
    pass

##############################################
##############################################

##Task level tests : Stokes imaging options
class test_stokes(testref_base):

    # Test 6
    def test_stokes_mfs_I(self):
        """ [stokes] test_stokes_mfs_I """
        ######################################################################################
        # Test_Stokes_I_mfs mfs with stokes I. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point_linRL.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','stokes':'I'})
        deconvolve(imagename=self.img, niter=10)
        report=th.checkall(imgexist=[self.img+'.image'],
                           imgval=[(self.img+'.image',1.0,[50,50,0,0])])
        self.checkfinal(report)

    # Test 7
    def test_stokes_mtmfs_IQUV(self):
        """ [stokes] test_stokes_mtmfs_IQUV """
        ######################################################################################
        # Test mtmfs with stokes IQUV. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point_linRL.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','stokes':'IQUV','deconvolver':'mtmfs','nterms':2})
        deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs', nterms=2)
        report=th.checkall(imgexist=[self.img+'.image.tt0'],imgexistnot=[self.img+'.image.alpha'],
                           imgval=[(self.img+'.image.tt0',1.0,[50,50,0,0]),(self.img+'.image.tt0',2.0,[50,50,1,0]), (self.img+'.image.tt0',3.0,[50,50,2,0]),(self.img+'.image.tt0',4.0,[50,50,3,0]) ])
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : cube.
class test_cube(testref_base):

    # def __init__(self,methodName='runTest'):
    #      testref_base.__init__(self,methodName)
    #      self.test_cube_0.__func__.__doc__ %="aaaa"

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
        
        # self.test_cube_0.__func__.__doc__ %=self.testList[0]['desc']
    

    def get_cubetclean_args(self, testid):
        """ core function to execute a cube tclean """
        if 'interpolation' in self.testList[testid]:
            interpolation = self.testList[testid]['interpolation']
        else:
            interpolation = 'linear'

        tclean_args = { 'field':'0', 'imsize':100, 'cell':'8.0arcsec',
                        'specmode':'cube', 'nchan':10, 'restfreq':['1.25GHz'],
                        'phasecenter':"J2000 19:59:28.500 +40.44.01.50", 'deconvolver':'hogbom',
                        'spw':self.testList[testid]['spw'],
                        #'imagename':self.img+self.testList[testid]['imagename'],
                        'start':self.testList[testid]['start'],
                        'width':self.testList[testid]['width'],
                        'veltype':self.testList[testid]['veltype'],
                        'outframe':self.testList[testid]['outframe'],
                        'interpolation':interpolation }
        return tclean_args

    # Test 8
    @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_0(self):
        """ [cube] test_cube_0 """
        ######################################################################################
        # Test_Cube_0 new . Should produce the same results as tclean.
        ######################################################################################
        testid=0
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50002,[50,50,0,0])])
        report2 = th.checkspecframe(self.img+'.image','LSRK',999988750)
        self.checkfinal(report)#+report2)

    # Test 9
    def test_cube_chanchunks_auto(self):
        """ [cube] test_cube_chanchunks_auto """
        ######################################################################################
        # Test channel chunking for large cubes : automatic calc of nchanchunks . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point.ms', tclean_args={'specmode':'cube','imsize':100,'cell':'10.0arcsec','deconvolver':'hogbom','chanchunks':-1})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom')
        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
        report=th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.5002,[50,50,0,0]) , (self.img+'.image',0.769,[50,50,0,19]) ])
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : masks and clean boxes.
class test_mask(testref_base):

    # Test 10
    def test_mask_1(self):
        """ [mask] test_mask_1 """
        ######################################################################################
        # Test Input mask as file and string : mfs . Should produce the same results as tclean.
        ######################################################################################
        mstr = 'circle[[50pix,80pix],10pix]'
        self.delData('refim_twochan.ms') # delete data here, since we're not doing that in prepData
        
        th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
        # delold=False -> don't delete the mask file
        self.prepData('refim_twochan.ms', delold=False, tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','usemask':'user','mask':self.img+'.mask.txt'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=self.img+'.mask.txt')
        report1=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,80,0,0])])

        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','usemask':'user','mask':mstr})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=mstr)
        report2=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,80,0,0])])

        self.checkfinal(report1+report2)

    # Test 11
    def test_mask_2(self):
        """ [mask] test_mask_2 """
        ######################################################################################
        # Test  Input mask as file and string : cube (few channels) . Should produce the same results as tclean.
        ######################################################################################
        mstr =  'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
        self.delData('refim_point.ms') # delete data here, since we're not doing that in prepData

        th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
        # delold=False -> don't delete the mask file
        self.prepData('refim_point.ms', delold=False, tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','specmode':'cube','interactive':0,'usemask':'user','mask':self.img+'.mask.txt'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=self.img+'.mask.txt')
        report1=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,1]),(self.img+'.mask',1.0,[50,50,0,2]),(self.img+'.mask',1.0,[50,50,0,10]),(self.img+'.mask',0.0,[50,50,0,11])])

        self.prepData('refim_point.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','specmode':'cube','interactive':0,'usemask':'user','mask':mstr})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=mstr)
        report2=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,1]),(self.img+'.mask',1.0,[50,50,0,2]),(self.img+'.mask',1.0,[50,50,0,10]),(self.img+'.mask',0.0,[50,50,0,11])])

        self.checkfinal(report1+report2)

    # Test 12
    def test_mask_missingfile(self):
        """ [imgval] test_mask_missingfile """
        ######################################################################################
        # tst.mask is sometimes required
        ######################################################################################
        mstr = 'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
        mstr = '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n'
        mname = self.img+'.mask.txt'
        th.write_file(mname, mstr)
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','usemask':'user','mask':mname})
        deconvolve_args={'imagename':self.img, 'niter':10, 'deconvolver':'hogbom', 'usemask':'user', 'mask':mname}

        # because tclean apparently deletes the mask file... ugh
        self.assertFalse(os.path.exists(mname))
        th.write_file(mname, mstr)

        # this should be fine ('tst.mask.txt' is present and 'tst.mask' is not)
        os.system("rm -rf "+self.img+".mask")
        deconvolve(**deconvolve_args)

        # this should error out ('tst.mask.txt' has been deleted)
        os.system("rm -rf "+self.img+".mask")
        os.system("rm -rf "+mname)
        strcheck = "'mask' parameter specified as a filename '"+mname+"', but no such file exists"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(**deconvolve_args)
        else:
            try:
                deconvolve(**deconvolve_args)
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # AUTOMASK TESTS
    # Test 13
    def test_mask_autobox_multithresh(self):
        """ [mask] test_mask_autobox_multithresh """
        ######################################################################################
        # Test multi-threshold Autobox (default). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','interactive':0,'usemask':'auto-multithresh'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='auto-multithresh')
        report=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : run deconvolve multiple times in a row
class test_multirun(testref_base):

    # Test 14
    def test_multirun_hoghog(self):
        """ [onefield] test_multirun_hoghog """
        ######################################################################################
        # Test running hogbom twice in a row and show that it gets the same value as one run with twice the iterations . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','threshold':'1mJy'})
        results = deconvolve(imagename=self.img, deconvolver='hogbom', niter=20, threshold='1mJy', interactive=0)
        report1=th.checkall(ret=results['retrec'], peakres=0.148, modflux=1.008, iterdone=20, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                            imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','threshold':'1mJy'})
        results2 = deconvolve(imagename=self.img, deconvolver='hogbom', niter=10, threshold='1mJy', interactive=0, restoration=False)
        report2=th.checkall(ret=results2['retrec'], peakres=0.353, modflux=0.772, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'],
                            imgval=[(self.img+'.model',0.772,[50,50,0,0])])
        results3 = deconvolve(imagename=self.img, deconvolver='hogbom', niter=10, threshold='1mJy', interactive=0)
        report3=th.checkall(ret=results3['retrec'], peakres=0.148, modflux=1.008, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                            imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.checkfinal(report1 + report2 + report3)

    # Test 15
    def test_multirun_multiscalehog(self):
        """ [onefield] test_multirun_multiscalehog """
        ######################################################################################
        # Test running multiscale clean followed by hogbom . Should produce the same results as tclean.
        # Note: we only test for task completion, don't know what the value should be at the end. (TODO needs validation)
        ######################################################################################
        self.prepData('refim_eptwochan.ms', tclean_args={'imsize':200, 'cell':'8.0arcsec', 'deconvolver':'multiscale', 'scales':[0,20,40,100]})
        results1 = deconvolve(imagename=self.img, niter=10, deconvolver='multiscale', scales=[0,20,40,100], interactive=0, restoration=False)
        report1 = th.checkall(ret=results1['retrec'], peakres=0.823, modflux=3.816, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.model'], imgexistnot=[self.img+'.image'])
        results2 = deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0)
        report2 = th.checkall(ret=results2['retrec'], iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image', self.img+'.model'])

        self.checkfinal(report1 + report2)

    # Test 16
    def test_multirun_norestore_restore(self):
        """ [onefield] test_multirun_norestore_restore """
        ######################################################################################
        # Test to test the retore-only feature . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point.ms', tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec']})
        results1 = deconvolve(imagename=self.img, niter=10, restoration=False)
        report1=th.checkall(ret=results1['retrec'], iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.model'], imgexistnot=[self.img+'.image'])
        results2 = deconvolve(imagename=self.img, niter=0, restoration=True)
        report2=th.checkall(ret=results2['retrec'], iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.model', self.img+'.image'],
                            imgval=[(self.img+'.image',0.482,[50,49,0,0])] )
        self.checkfinal(report1 + report2)

##############################################
##############################################

##Task level tests : test validation of input images and copying/converting startmodel images
class test_imgval(testref_base):
    def ivsetup(self, stopearly=False):
        self.prepData('refim_point.ms', cache=True, stopearly=stopearly, tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec']})
        self.mname = self.img + ".model"
        self.mname2 = self.img + "_2.model"

    def tearDown(self):
        super(test_imgval, self).tearDown()

        # remove all tst_2.* and tst_bak.* images
        for uname in ["_bak", "_2"]:
            fn = self.img + uname + ".*"
            if os.path.exists(fn):
                os.system('rm -rf '+fn)

    def get_csys_crval0(self, img):
        _ia.open(img)
        csys = _ia.coordsys()
        # pnt = csys.torecord()['direction0']['crval'][0]
        pnt = csys.torecord()['direction0']['crpix'][0]
        _ia.close()
        return csys, pnt

    def set_crval0(self, img, crval0):
        _ia.open(img)
        csys = _ia.coordsys()
        rec = csys.torecord()
        # rec['direction0']['crval'][0] = crval0
        rec['direction0']['crpix'][0] = crval0
        csys.fromrecord(rec)
        _ia.setcoordsys(csys.torecord())
        _ia.close()

    def get_shape(self, img):
        _ia.open(img)
        ret = _ia.shape()
        _ia.close()
        return ret

    def tst_imgval_missingimgs(self, ext):
        # helper method for test_imgval_missingimgs_*
        self.ivsetup()
        os.system("mv {0}{1} {0}_bak{1}".format(self.img, ext))

        strcheck = "missing one or more of the required images"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(imagename=self.img, niter=10)
        else:
            try:
                deconvolve(imagename=self.img, niter=10)
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # Test 17
    def test_imgval_missingimgs_residual(self):
        """ [imgval] test_imgval_missingimgs_residual """
        ######################################################################################
        # tst.residual and tst.psf are always required
        ######################################################################################
        self.tst_imgval_missingimgs(".residual")

    # Test 18
    def test_imgval_missingimgs_psf(self):
        """ [imgval] test_imgval_missingimgs_psf """
        ######################################################################################
        # tst.residual and tst.psf are always required
        ######################################################################################
        # Note: cpp code doesn't throw an exception when psf is missing, just prints a warning
        self.tst_imgval_missingimgs(".psf")

    # Test 19
    def test_imgval_missingimgs_model(self):
        """ [imgval] test_imgval_missingimgs_models """
        ######################################################################################
        # tst.model is used to continue deconvolution, but is not required.
        ######################################################################################
        self.ivsetup()
        if (os.path.exists(self.img+'.model')):
            os.system("rm -rf "+self.img+".model")
        deconvolve(imagename=self.img, niter=10)

    # Test 20
    def test_imgval_missingimgs_sumwt(self):
        """ [imgval] test_imgval_missingimgs_sumwt """
        ######################################################################################
        # tst.sumwt is never required
        ######################################################################################
        self.ivsetup()
        if (os.path.exists(self.img+'.sumwt')):
            os.system("rm -rf "+self.img+".sumwt")

        # Should be fine. Sumwt should not be required for task deconvolve.
        deconvolve(imagename=self.img, niter=10)

    def tst_imgval_axesmismatch(self, ext):
        # helper method for test_imgval_axesmismatch_*
        self.ivsetup()

        # cause the tst.model image to exist prior to looking for it
        deconvolve(imagename=self.img, niter=1, restoration=False)

        # verify the images have the wrong format and an exception is thrown
        # Note: cpp code does not mind that the .pb image has a weird axes order
        fn1 = self.img + ext
        fn2 = self.img + "_bak" + ext
        os.system("mv {0} {1}".format(fn1, fn2))

        imtrans(imagename=self.img+"_bak"+ext, outfile=self.img+ext, order="3012")
        strcheck = "There is a shape mismatch between existing images"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(imagename=self.img, niter=10)
        else:
            try:
                deconvolve(imagename=self.img, niter=10)
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # Test 21
    def test_imgval_axesmismatch_residual(self):
        """ [imgval] test_imgval_axesmismatch_residual """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.tst_imgval_axesmismatch(".residual")

    # Test 22
    def test_imgval_axesmismatch_psf(self):
        """ [imgval] test_imgval_axesmismatch_psf """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.tst_imgval_axesmismatch(".psf")

    # Test 23
    def test_imgval_axesmismatch_model(self):
        """ [imgval] test_imgval_axesmismatch_model """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.tst_imgval_axesmismatch(".model")

    # Test 24
    def test_imgval_axesmismatch_pb(self):
        """ [imgval] test_imgval_axesmismatch_pb """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        # Note: cpp code does not mind that the .pb image has a weird axes order
        # self.tst_imgval_axesmismatch(".pb")
        pass

    def tst_imgval_shapemismatch(self, ext):
        """All input images should have the same shape as the tst.residual image"""
        self.ivsetup()

        # cause the tst.model image to exist prior to looking for it
        deconvolve(imagename=self.img, niter=1, restoration=False)

        # verify the images have the wrong format and an exception is thrown
        # Note: cpp code does not mind that the .pb image has a weird shape
        os.system("mv {0}{1} {0}_bak{1}".format(self.img, ext))

        imrebin(imagename=self.img+"_bak"+ext, outfile=self.img+ext, factor=[50,50])
        strcheck = "There is a shape mismatch between existing images"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(imagename=self.img, niter=10)
        else:
            try:
                deconvolve(imagename=self.img, niter=10)
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # Test 25
    def test_imgval_shapemismatch_residual(self):
        """ [imgval] test_imgval_shapemismatch_residual """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.tst_imgval_shapemismatch(".residual")

    # Test 26
    def test_imgval_shapemismatch_psf(self):
        """ [imgval] test_imgval_shapemismatch_psf """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.tst_imgval_shapemismatch(".psf")

    # Test 27
    def test_imgval_shapemismatch_model(self):
        """ [imgval] test_imgval_shapemismatch_model """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.tst_imgval_shapemismatch(".model")

    # Test 28
    def test_imgval_shapemismatch_pb(self):
        """ [imgval] test_imgval_shapemismatch_pb """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        # Note: cpp code does not mind that the .pb image has a weird shape
        # self.tst_imgval_shapemismatch(".pb")
        pass
    
    # TODO figure out why running the startmodel_axesmismatch test immediately before this test causes an exception to be thrown
    # Test 29
    @unittest.skip("if test_imgval_startmodel_axesmismatch executes immediately before this test then it fails")
    def test_imgval_startmodel_empty(self):
        """ [imgval] test_imgval_startmodel_empty """
        ######################################################################################
        # Deconvolve should ignore all empty strings entered for the startmodel
        ######################################################################################
        # self.test_imgval_startmodel_axesmismatch()
        self.ivsetup()

        # basic test with empty string as startmodel
        deconvolve(imagename=self.img, niter=10, startmodel='')

        # basic test with list of empty strings as startmodel
        deconvolve(imagename=self.img, niter=10, startmodel=['', '', ''])

        # basic copy test where empty string is discarded from list
        self.mname = self.img + ".model"
        self.mname2 = self.img + "_2.model"
        os.system("mv {0} {1}".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=['', '', self.mname2, '', ''])

    # Test 30
    def test_imgval_startmodel_dne(self):
        """ [imgval] test_imgval_startmodel_dne """
        ######################################################################################
        # Throws an error if startmodel is set but does not exist
        ######################################################################################
        self.ivsetup()
        strcheck = "does not exist"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(imagename=self.img, niter=10, startmodel='doesnotexists.model')
        else:
            try:
                deconvolve(imagename=self.img, niter=10, startmodel='doesnotexists.model')
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # Test 31
    def test_imgval_startmodel_model_exists(self):
        """ [imgval] test_imgval_startmodel_model_exists """
        ######################################################################################
        # Throws an error if startmodel is set and tst.model exists (must be one or the other, not both)
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10)
        shutil.copytree(self.mname, self.mname2)
        
        strcheck = "exists"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        else:
            try:
                deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # TODO figure out why running the startmodel_axesmismatch test immediately before this test causes an exception to be thrown
    # Test 32
    @unittest.skip("if test_imgval_startmodel_axesmismatch executes immediately before this test then it fails")
    def test_imgval_startmodel_basic_copy(self):
        """ [imgval] test_imgval_startmodel_basic_copy """
        ######################################################################################
        # Tests ability of deconvolve to copy startmodel to tst.model before starting deconvolution
        ######################################################################################
        # self.test_imgval_startmodel_axesmismatch()
        self.ivsetup(stopearly=False)
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of
        os.system("mv {0} {1}".format(self.mname, self.mname2))
            
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get moved to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get copied!".format(self.mname))

    # Test 33
    def test_imgval_startmodel_axesmismatch(self):
        """ [imgval] test_imgval_startmodel_axesmismatch """
        ######################################################################################
        # Tests the existing functionality. If in the future the logic is added to auto-translate images, this test can be removed.
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of
        imtrans(imagename=self.mname, outfile=self.mname2, order="3012")
            
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get translated to {1}!".format(self.mname, self.mname2))
        strcheck = "Error in setting"
        if is_CASA6:
            with self.assertRaisesRegex(RuntimeError, strcheck):
                deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        else:
            try:
                deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
                self.assertFalse("Error case failed to throw exception")
            except RuntimeError as e:
                self.assertTrue(strcheck in str(e))

    # Test 34
    def test_imgval_startmodel_csysmismatch(self):
        """ [imgval] test_imgval_startmodel_csysmismatch """
        ######################################################################################
        # Tests the ability of the deconvolve regrid the csys of the startmodel to that of tst.residual
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of
        os.system("mv {0} {1}".format(self.mname, self.mname2))

        # change csys
        csys, oldpnt = self.get_csys_crval0(self.mname2)
        newpnt = 51#oldpnt * 0.9999
        self.assertNotEqual(oldpnt, newpnt, "Change amount not big enough")
        self.set_crval0(self.mname2, newpnt)
        csys2, newpnt2 = self.get_csys_crval0(self.mname2)
        self.assertEqual(newpnt2, newpnt, "Image {0} did not get its csys.direction0.crval[0] value updated properly from {1} to the expected {2}! (actual value is {3})".format(self.mname2, oldpnt, newpnt, newpnt2))

        # test that deconvolve regrids the image
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get rebinned to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get regridded back!".format(self.mname))
        csys3, regridpnt = self.get_csys_crval0(self.mname)
        self.assertAlmostEqual(regridpnt, oldpnt, "Image {0} did not get its csys.direction0.crval[0] value regridded properly from {1} to {2}! (actual value is {3})".format(self.mname2, newpnt, oldpnt, regridpnt))

    # Test 35
    def test_imgval_startmodel_shapemismatch(self):
        """ [imgval] test_imgval_startmodel_shapemismatch """
        ######################################################################################
        # Tests the ability of the deconvolve regrid the shape of the startmodel to that of tst.residual
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of

        # change shape
        s = self.get_shape(self.mname)
        self.assertEqual([s[0], s[1]], [100, 100], "Image shape ({0}) didn't start as expected ({1})".format(s, [100, 100]))
        imrebin(imagename=self.mname, outfile=self.mname2, factor=[2,2])

        # sanity: make sure the shape changed
        s = self.get_shape(self.mname2)
        self.assertEqual([s[0], s[1]], [50, 50], "Image shape ({0}) didn't get rebinned as expected ({1})".format(s, [50, 50]))

        # run deconvolve and make sure the shape gets regridded back in
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get rebinned to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get regridded back!".format(self.mname))
        s = self.get_shape(self.mname)
        self.assertEqual([s[0], s[1]], [100, 100], "Image shape ({0}) didn't get regridded as expected ({1})".format(s, [100, 100]))


##############################################
##############################################

##Task level tests : verify that the .residual image gets updated between consecutive runs of deconvolve
class test_residual_update(testref_base):
    def get_stats(self):
        _ia.open(self.img+".residual")
        ret = _ia.statistics()
        _ia.close()
        return ret

    def cmp_stats(self, stats1, stats2):
        errmsg = "Difference between residual statistics in first and second run didn't change.\n" +\
                 "First run: {0}\n\nSecond run: {1}".format(stats1, stats2)
        self.assertNotAlmostEqual(stats1['max'][0], stats2['max'][0], delta=0.00001, msg=errmsg)
        self.assertNotAlmostEqual(stats1['min'][0], stats2['min'][0], delta=0.00001, msg=errmsg)
        self.assertNotAlmostEqual(stats1['mean'][0], stats2['mean'][0], delta=0.00001, msg=errmsg)

    def tst_residual_update(self, deconvolver):
        """Helper method to execute the residual update tests for non-multiterm images"""
        self.prepData('refim_point.ms', stopearly=stopearly, tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec'], 'deconvolver':deconvolver})

        # run and get the statistics
        deconvolve(imagename=self.img, deconvolver=deconvolver, niter=10)
        stats1 = self.get_stats()
        deconvolve(imagename=self.img, deconvolver=deconvolver, niter=10)
        stats2 = self.get_stats()

        # verify there's a difference in the stats
        self.cmp_stats(stats1, stats2)

    # Test 36
    def test_residual_update_hogbom(self):
        """ [residual_update] test_residual_update_hogbom """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that hogbom does this correctly.
        ######################################################################################
        self.tst_residual_update('hogbom')

    # Test 37
    def test_residual_update_clark(self):
        """ [residual_update] test_residual_update_clark """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that clark does this correctly.
        ######################################################################################
        self.tst_residual_update('clark')

    # Test 38
    def test_residual_update_clarkstokes(self):
        """ [residual_update] test_residual_update_clarkstokes """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that clarkstokes does this correctly.
        ######################################################################################
        self.tst_residual_update('clarkstokes')

    # Test 39
    def test_residual_update_multiscale(self):
        """ [residual_update] test_residual_update_multiscale """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that multiscale does this correctly.
        ######################################################################################
        self.tst_residual_update('multiscale')

    # Test 40
    def test_residual_update_mtmfs(self):
        """ [residual_update] test_residual_update_mtmfs """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that mtmfs does this correctly.
        ######################################################################################
        self.tst_residual_update('mtmfs')

    # Test 41
    def test_residual_update_mem(self):
        """ [residual_update] test_residual_update_mem """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that mem does this correctly.
        ######################################################################################
        self.tst_residual_update('mem')

if __name__ == '__main__':
    unittest.main()
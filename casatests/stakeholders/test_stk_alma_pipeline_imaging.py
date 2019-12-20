##########################################################################
##########################################################################
# test_tclean_ALMA_pipeline
#
# Author: Shawn Thomas Booth


'''
Datasets (MOUS)
E2E6.1.00034.S (uid://A002/Xcff05c/X1ec)
2018.1.00879.S (uid://A001/X133d/X169f)
E2E6.1.00020.S (uid://A002/Xcff05c/Xe5)
2017.1.00750.T (uid://A001/X131b/X57)

Test list
- Single field cube - E2E6.1.00034.S
- SF MFS - E2E6.1.00020.S
- SF mtmfs - E2E6.1.00020.S
- SF ephemeris cube (multi-EB) - 2017.1.00750.T
- SF ephemeris MFS - 2018.1.00879.S
- SF ephemeris mtmfs - 2018.1.00879.S
- SF Calibrator - E2E6.1.00034.S
- Mosaic cube - E2E6.1.00034.S
- Mosaic MFS - E2E6.1.00020.S
- Mosaic mtmfs - E2E6.1.00020.S
- Mosaic ephemeris cube - 2018.1.00879.S
- Mosaic ephemeris MFS - 2018.1.00879.S
- Mosaic ephemeris mtmfs - 2018.1.00879.S
'''
##########################################################################
##########################################################################

# Imports #
import os
import subprocess
import glob
import time
import sys
import unittest
from __main__ import default  # reset given task to its default values
from tasks import *  # Imports all casa tasks
from taskinit import *  # Imports all casa tools
import numpy
import shutil
import inspect
from parallel.parallel_task_helper import ParallelTaskHelper

import casaTestHelper
#data_path = '/lustre/naasc/sciops/comm/sbooth/CASA_ALMA_pipeline/data_dir/'
data_path = os.environ.get('CASAPATH').split()[0] + '/data/stakeholder/alma'


## Base Test class with Utility functions
class test_tclean_base(unittest.TestCase):

    def setUp(self):
        self._myia = iatool()
        self.epsilon = 0.01 # sets epsilon as a percentage (1%)
        self.msfile = ""
        self.img = "tst"
        # To use subdir in the output image names in some tests (CAS-10937)
        self.img_subdir = 'refimager_tst_subdir'
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

#    def rename(self, file_name, file_list):
#        for item in file_list:
#            sp_item = item.split('iter0')
#            process = subprocess.Popen(['mv', item, \
#                                        sp_item[0]+'iter1'+sp_item[1]])
#            process.wait()
#        if self.parallel:
#            subfile_list = glob.glob(file_name+'1.workdirectory/*')
#            for item in subfile_list:
#                sp_item = item.split('iter0')
#                process = subprocess.Popen(['mv', item, \
#                                            sp_item[0]+'iter1'+sp_item[1]])
#                process.wait()

    def rename(self, file_name, file_list):
        for item in file_list:
            sp_item = item.split('iter0')
            shutil.copytree(item, sp_item[0]+'iter1'+sp_item[1])
        if self.parallel:
            subfile_list = glob.glob(file_name+'1.workdirectory/*')
            for item in subfile_list:
                sp_item = item.split('iter0')
                shutil.copytree(item, sp_item[0]+'iter1'+sp_item[1])

    def tearDown(self):
        casaTestHelper.generate_weblog("tclean_ALMA_pipeline",test_dict)
        self._myia.done()
        """ don't delete it all """
#        self.delData()

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self, msname=""):
        if msname != "":
             self.msfile=msname

    def delData(self, msname=""):
        if msname != "":
             self.msfile=msname
#        if (os.path.exists(self.msfile)):
#             os.popen('rm -rf ' + self.msfile)
#        os.popen('rm -rf ' + self.img_subdir)
#        os.popen('rm -rf ' + self.img+'*')

    def prepInputmask(self, maskname=""):
        if maskname!="":
            self.maskname=maskname
        if (os.path.exists(self.maskname)):
            os.popen('rm -rf ' + self.maskname)
        shutil.copytree(refdatapath+self.maskname, self.maskname)

    def checkfinal(self, pstr=""):
#        pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  casa "
#                "-c `echo $CASAPATH | awk '{print $1}'`/gcwrap/python/"
#                "scripts/regressions/admin/runUnitTest.py test_refimager["
#                + inspect.stack()[1][3] +"]"
        pstr += "["+inspect.stack()[1][3]+"] : To re-run this test : " \
                "runUnitTest.main(['test_tclean["+ inspect.stack()[1][3] \
                +"]'])"
        casalog.post(pstr,'INFO')
        if( pstr.count("(Fail") > 0 ):
             self.fail("\n"+pstr)

    def image_stats(self, img, suf):
        self._myia.open(img+suf)
        rms = self._myia.statistics()['rms'][0]  # reads RMS value for image
        int_flux = self._myia.statistics()['sum'][0]  # reads integrated 
                                                      # flux value for image
        major = self._myia.commonbeam()['major']['value']
        minor = self._myia.commonbeam()['minor']['value']
        pa = self._myia.commonbeam()['pa']['value']
        self._myia.close()

        return rms, int_flux, major, minor, pa

    def image_list(self, img, mode):
        standard = [img+'.psf', img+'.residual', img+'.image', \
                    img+'.image.pbcor', img+'.mask', img+'.pb', \
                    img+'.model', img+'.sumwt']
        mosaic = [img+'.weight']
        mtmfs = [img+'.alpha', img+'.alpha.error', img+'.alpha.pbcor', \
                 img+'.psf.tt0', img+'.psf.tt1', img+'.psf.tt2', \
                 img+'.residual.tt0', img+'.residual.tt1', img+'.image.tt0',\
                 img+'.image.tt1', img+'.image.tt0.pbcor', \
                 img+'.image.tt1.pbcor', img+'.mask', img+'.pb.tt0', \
                 img+'.model.tt0', img+'.model.tt1', img+'.sumwt.tt0', \
                 img+'.sumwt.tt1', img+'.sumwt.tt2']
        mos_mtmfs = [img+'.weight.tt0', img+'.weight.tt1', img+'.weight.tt2']

        if mode == 'standard':
            img_list = standard
        if mode == 'mosaic':
            img_list = standard+mosaic
        if mode == 'mtmfs':
            img_list = mtmfs
        if mode == 'mos_mtmfs':
            img_list = mtmfs+mos_mtmfs

        return img_list


##############################################
##############################################
test_dict = {}
class Test_standard(test_tclean_base):


    @casaTestHelper.stats_dict(test_dict)
    def test_standard_cube(self):
        '''Standard (single field) cube imaging'''

        file_name = 'standard_cube.iter'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], \
               imagename=file_name+'0', field='1', spw=['0'], \
               imsize=[108, 108], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS'
               ' 00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
               nchan=508, start='220.2526744GHz', width='0.244174087287MHz', \
               outframe='LSRK', pblimit=0.2, perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               deconvolver='hogbom', usepointing=False, restoration=False, \
               pbcor=False, weighting='briggs', restoringbeam='common', \
               robust=0.5, npixels=0, niter=0, threshold='0.0mJy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=self.parallel)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], \
               imagename=file_name+'1', field='1', spw=['0'], \
               imsize=[108, 108], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS'
               ' 00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
               nchan=508, start='220.2526744GHz', width='0.244174087287MHz', \
               outframe='LSRK', perchanweightdensity=False, \
               usepointing=False, pblimit=0.2, nsigma=0.0, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               deconvolver='hogbom', restoration=True, \
               pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
               niter=20000, threshold='0.354Jy', interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, calcres=False, calcpsf=False, \
               savemodel='none', parallel=self.parallel)

        if self.parallel:
            tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], \
                   imagename=file_name+'1', spw=['0'], field='1', \
                   imsize=[108, 108], cell=['1.1arcsec'], phasecenter='ICRS'
                   ' 00:45:54.3836 -073.15.29.413', stokes='I', \
                   antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                   intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                   specmode='cube', nchan=508, start='220.2526744GHz', \
                   width='0.244174087287MHz', outframe='LSRK', \
                   perchanweightdensity=False, gridder='mosaic', \
                   chanchunks=-1, mosweight=False, usepointing=False, \
                   pblimit=0.2, deconvolver='hogbom', restoration=True, \
                   restoringbeam='common', pbcor=True, weighting='briggs', \
                   robust=0.5, npixels=0, niter=0, threshold='0.354Jy', \
                   nsigma=0.0, interactive=0, usemask='auto-multithresh', \
                   sidelobethreshold=1.25, noisethreshold=5.0, \
                   lownoisethreshold=2.0, negativethreshold=0.0, \
                   minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                   minpercentchange=1.0, fastnoise=False, restart=True, \
                   calcres=False, calcpsf=False, savemodel='none', \
                   parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.143988420605,
                           'Integrated Flux':-48.1222390321,
                           'Beam - Major Axis':8.43566099983,
                           'Beam - Minor Axis':5.96328384076,
                           'Beam - Position Angle':69.882947913}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'standard'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks that a pixel value at edge of pb is not flagged &
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [52, 84, 0, 0]), \
                       (img+'.image', False, [52, 85, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', \
            epsilon = abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "E2E6.1.00034.S_tclean.ms")

        test_dict["test_standard_cube"]['self.parallel'] = self.parallel
        test_dict["test_standard_cube"]['report'] = report
        immoments(imagename = img+'.image', moments = 8, outfile = \
                  img+'.image.moment8')
        imview(img+'.image.moment8', out = img+'.image.moment8.png')
        test_dict["test_standard_cube"]['images'] = \
            [img+'.image.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mfs(self):
        '''Standard (single field) MFS imaging'''

        file_name = 'standard_mfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='2', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'0', stokes='I', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS'
               ' 13:56:07.2100 +005.15.17.200', specmode='mfs', nchan=-1, \
               outframe='LSRK', perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               nterms=2, restoration=False, restoringbeam='common', \
               pbcor=False, weighting='briggs', robust=0.5, npixels=0, \
               niter=0, threshold='0.0mJy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='2', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'1', stokes='I', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS'
               ' 13:56:07.2100 +005.15.17.200', specmode='mfs', nchan=-1, \
               outframe='LSRK', perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               nterms=2, restoration=True, restoringbeam='common', \
               pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
               niter=30000, threshold='0.00723Jy', nsigma=0.0, \
               interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, restart=True, \
               calcres=False, calcpsf=False, savemodel='none', \
               parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.00419614745362, 
                           'Integrated Flux':1.78775458621,
                           'Beam - Major Axis':17.1119709015,
                           'Beam - Minor Axis':9.91198253632,
                           'Beam - Position Angle':87.5630874634}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'standard'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 0.0407012216747, [64, 63, 0, 0])], 
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [63, 95, 0, 0]), \
                       (img+'.image', False, [64, 96, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', \
            epsilon = abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "E2E6.1.00020.S_tclean.ms")

        test_dict["test_standard_mfs"]['self.parallel'] = self.parallel
        test_dict["test_standard_mfs"]['report'] = report
        imview(img+'.image', out = img+'.image.png')
        test_dict["test_standard_mfs"]['images'] = [img+'.image.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mtmfs(self):
        '''Single field mtmfs imaging'''

        file_name = 'standard_mtmfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='2', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'0', stokes='I', \
               imsize=[126, 126], cell=['2arcsec'],phasecenter='ICRS'
               ' 13:56:07.2100 +005.15.17.200', specmode='mfs', nchan=-1, \
               outframe='LSRK', perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               nterms=2, restoration=False, restoringbeam='common', \
               pbcor=False, weighting='briggs', robust=0.5, npixels=0, \
               niter=0, threshold='0.0mJy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='2', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'1', stokes='I', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS'
               ' 13:56:07.2100 +005.15.17.200', specmode='mfs', nchan=-1, \
               outframe='LSRK', perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               nterms=2, restoration=True, restoringbeam='common', \
               pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
               niter=30000, threshold='0.00723Jy', nsigma=0.0, \
               interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, restart=True, \
               calcres=False, calcpsf=False, savemodel='none', \
               parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image.tt0')
        expected_values = {'RMS':0.00417779723296, 
                           'Integrated Flux':1.78041236983,
                           'Beam - Major Axis':17.7945289612,
                           'Beam - Minor Axis':10.0446176529,
                           'Beam - Position Angle':86.6707839966}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mtmfs'),
            # checks the peak flux value and location
            imgval = [(img+'.image.tt0', 0.0372908785939, [63, 63, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image.tt0', True, [62, 95, 0, 0]), \
                       (img+'.image.tt0', False, [62, 96, 0, 0])]) 

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = 
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "E2E6.1.00020.S_tclean.ms")

        test_dict["test_standard_mtmfs"]['self.parallel'] = self.parallel
        test_dict["test_standard_mtmfs"]['report'] = report
        imview(img+'.image.tt0', out = img+'.image.tt0.png')
        test_dict["test_standard_mtmfs"]['images'] = [img+'.image.tt0.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_cube_eph(self):
        '''Single field multi-EB ephemeris cube imaging'''

        file_name = 'standard_cube_eph.iter'
        self.prepData(data_path+'2017.1.00750.T_tclean_exe1.ms')
        self.prepData(data_path+'2017.1.00750.T_tclean_exe2.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'2017.1.00750.T_tclean_exe1.ms',
               data_path+'2017.1.00750.T_tclean_exe2.ms'], \
               field='21PGiacobini-Zinner', spw=['0', '0'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11', \
               '0,1,2,3,4,5,6,7,8,9'], scan=['7,11,15,19,23',
               '8,12,16,20,24'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'0', imsize=[80, 80], \
               cell=['0.66arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='cubesource', nchan=1000, start=1550, width=1, \
               perchanweightdensity=False, gridder='standard', \
               chanchunks=-1, mosweight=False, usepointing=False, \
               pblimit=0.2, deconvolver='hogbom', restoration=False, \
               restoringbeam='common', pbcor=False, weighting='briggs', \
               robust=0.5, npixels=0, niter=0, threshold='0.0mJy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'2017.1.00750.T_tclean_exe1.ms', \
               data_path+'2017.1.00750.T_tclean_exe2.ms'], \
               field='21PGiacobini-Zinner', spw=['0', '0'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11', \
               '0,1,2,3,4,5,6,7,8,9'], scan=['7,11,15,19,23', \
               '8,12,16,20,24'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'1', imsize=[80, 80], \
               cell=['0.66arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='cubesource', nchan=1000, start=1550, width=1, \
               perchanweightdensity=False, gridder='standard', \
               chanchunks=-1, mosweight=False, usepointing=False, \
               pblimit=0.2, deconvolver='hogbom', restoration=True, \
               restoringbeam='common', pbcor=True, weighting='briggs', \
               robust=0.5, npixels=0, niter=30000, threshold='0.274Jy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, restart=True, \
               calcres=False, calcpsf=False, savemodel='none', \
               parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.0560069684214,
                           'Integrated Flux':2285.025188,
                           'Beam - Major Axis':4.49769604687,
                           'Beam - Minor Axis':3.3237527868,
                           'Beam - Position Angle':87.0964067383}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'standard'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 3.21681857109, [46, 41, 0, 489])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                       (img+'.image', False, [40, 73, 0, 0])]) 

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                    ["2017.1.00750.T_tclean_exe1.ms", \
                                     "2017.1.00750.T_tclean_exe2.ms"])

        test_dict["test_standard_cube_eph"]['self.parallel'] = self.parallel
        test_dict["test_standard_cube_eph"]['report'] = report
        immoments(imagename = img+'.image', moments = 8, \
                  outfile = img+'.image.moment8')
        imview(img+'.image.moment8', out = img+'.image.moment8.png')
        test_dict["test_standard_cube_eph"]['images'] = \
            [img+'.image.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mfs_eph(self):
        '''Standard (single field) ephemeris mfs imaging'''

        file_name = 'standard_mfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='0', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'0', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=2.0, \
               noisethreshold=4.25, lownoisethreshold=1.5, \
               negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               savemodel='none', parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='0', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'1', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=True, restoringbeam='common', pbcor=True, \
               weighting='briggs', robust=0.5, npixels=0, niter=7000000, \
               threshold='0.0316Jy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=2.0, \
               noisethreshold=4.25, lownoisethreshold=1.5, \
               negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, calcres=False, calcpsf=False, \
               savemodel='none', parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.470603939015,
                           'Integrated Flux':-1482.12338717,
                           'Beam - Major Axis':0.925863981247,
                           'Beam - Minor Axis':0.716557502747,
                           'Beam - Position Angle':-87.4016647339}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'standard'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 1.37355649471, [318, 222, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [240, 328, 0, 0]), \
                       (img+'.image', False, [240, 329, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "2018.1.00879.S_tclean.ms")

        test_dict["test_standard_mfs_eph"]['self.parallel'] = self.parallel
        test_dict["test_standard_mfs_eph"]['report'] = report
        imview(img+'.image', out = img+'.image.png')
        test_dict["test_standard_mfs_eph"]['images'] = [img+'.image.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mtmfs_eph(self):
        '''Standard (single field) ephemeris mtmfs imaging'''

        file_name = 'standard_mtmfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='0', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'0', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=2.0, \
               noisethreshold=4.25, lownoisethreshold=1.5, \
               negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               savemodel='none', parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='0', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'1', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='standard', chanchunks=-1, mosweight=False, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               restoration=True, restoringbeam='common', pbcor=True, \
               weighting='briggs', robust=0.5, npixels=0, niter=7000000, \
               threshold='0.0316Jy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=2.0, \
               noisethreshold=4.25, lownoisethreshold=1.5, \
               negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, calcres=False, calcpsf=False, \
               savemodel='none', parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image.tt0')
        expected_values = {'RMS':0.470603939015,
                           'Integrated Flux':-1482.12338717,
                           'Beam - Major Axis':0.925863981247,
                           'Beam - Minor Axis':0.716557502747,
                           'Beam - Position Angle':-87.4016647339}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mtmfs'),
            # checks the peak flux value and location
            imgval = [(img+'.image.tt0', 1.37355649471, [318, 222, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image.tt0', True, [240, 328, 0, 0]), \
                       (img+'.image.tt0', False, [240, 329, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "2018.1.00879.S_tclean.ms")

        test_dict["test_standard_mtmfs_eph"]['self.parallel'] = self.parallel
        test_dict["test_standard_mtmfs_eph"]['report'] = report
        imview(img+'.image.tt0', out = img+'.image.tt0.png')
        test_dict["test_standard_mtmfs_eph"]['images'] = \
           [img+'.image.tt0.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_cal(self):
        '''Calibrator image'''

        file_name = 'standard_cal.iter'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], field='0', \
               spw=['0'], antenna=['0,1,2,3,4,5,6,7,8'], scan=['3'], \
               intent='CALIBRATE_BANDPASS#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'0', imsize=[90, 90], stokes='I', \
               cell=['0.85arcsec'], phasecenter='ICRS 22:58:05.9629 '
               '-027.58.21.257', specmode='mfs', nchan=-1, outframe='LSRK', \
               perchanweightdensity=False, gridder='standard', \
               chanchunks=-1, mosweight=False, usepointing=False, \
               pblimit=0.2, deconvolver='hogbom', restoration=False, \
               restoringbeam='common', pbcor=False, weighting='briggs', \
               robust=0.5, npixels=0, niter=0, threshold='0.0mJy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.5, noisethreshold=6.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], field='0', \
               spw=['0'], antenna=['0,1,2,3,4,5,6,7,8'], scan=['3'], \
               intent='CALIBRATE_BANDPASS#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'1', imsize=[90, 90], stokes='I', \
               cell=['0.85arcsec'], phasecenter='ICRS 22:58:05.9629 '
               '-027.58.21.257', specmode='mfs', nchan=-1, outframe='LSRK', \
               perchanweightdensity=False, gridder='standard', \
               chanchunks=-1, mosweight=False, usepointing=False, \
               pblimit=0.2, deconvolver='hogbom', restoration=True, \
               restoringbeam='common', pbcor=True, weighting='briggs', \
               robust=0.5, npixels=0, niter=300000, threshold='0.0241Jy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.5, noisethreshold=6.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, restart=True, \
               calcres=False, calcpsf=False, savemodel='none', \
               parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.203889568665,
                           'Integrated Flux':172.689694168,
                           'Beam - Major Axis':9.98569583893,
                           'Beam - Minor Axis':4.62464284897,
                           'Beam - Position Angle':-86.3871307373}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'standard'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 2.40606951714, [45, 45, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [45, 85, 0, 0]), \
                       (img+'.image', False, [45, 86, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
                                   "E2E6.1.00034.S_tclean.ms")

        test_dict["test_standard_cal"]['self.parallel'] = self.parallel
        test_dict["test_standard_cal"]['report'] = report
        imview(img+'.image', out = img+'.image.png')
        test_dict["test_standard_cal"]['images'] = [img+'.image.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)

##############################################
##############################################

class Test_mosaic(test_tclean_base):


    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_cube(self):
        '''Mosaic cube imaging'''

        file_name = 'mosaic_cube.iter'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], \
               field='SMIDGE_NWCloud', spw=['0'], antenna=['0,1,2,3,4,5,'
               '6,7,8'], scan=['8,12,16'], intent='OBSERVE_TARGET'
               '#ON_SOURCE', datacolumn='data', imagename=file_name+'0', \
               imsize=[108, 108], cell=['1.1arcsec'], phasecenter='ICRS '
               '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
               nchan=508, start='220.2526744GHz', width='0.244174087287MHz', \
               outframe='LSRK', perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', interactive=0, usemask='auto-'
               'multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=self.parallel)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], \
               field='SMIDGE_NWCloud', spw=['0'], antenna=['0,1,2,3,4,5,'
               '6,7,8'], scan=['8,12,16'], intent='OBSERVE_TARGET'
               '#ON_SOURCE', datacolumn='data',  imagename=file_name+'1', \
               imsize=[108, 108], cell=['1.1arcsec'], phasecenter='ICRS '
               '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
               nchan=508, start='220.2526744GHz', width='0.244174087287MHz', \
               outframe='LSRK', perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=True, pbcor=True, weighting='briggs', \
               robust=0.5, npixels=0, niter=20000, threshold='0.354Jy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=1.25, noisethreshold=5.0, \
               lownoisethreshold=2.0, negativethreshold=0.0, \
               minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, restart=True, \
               savemodel='none', calcres=False, calcpsf=False, \
               parallel=self.parallel)

        if self.parallel:
            tclean(vis=[data_path+'E2E6.1.00034.S_tclean.ms'], \
                   field='SMIDGE_NWCloud', spw=['0'], antenna=['0,1,2,'
                   '3,4,5,6,7,8'], scan=['8,12,16'], intent='OBSERVE'
                   '_TARGET#ON_SOURCE', datacolumn='data', \
                   imagename=file_name+'1', imsize=[108, 108], \
                   cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                   ' -073.15.29.413', stokes='I', specmode='cube', \
                   nchan=508, start='220.2526744GHz', width= \
                   '0.244174087287MHz', outframe='LSRK', \
                   perchanweightdensity=False, gridder='mosaic', \
                   chanchunks=-1, mosweight=True, usepointing=False, \
                   pblimit=0.2, deconvolver='hogbom', restoration=True, \
                   restoringbeam='common', pbcor=True, weighting='briggs', \
                   robust=0.5, npixels=0, niter=0, threshold='0.354Jy', \
                   nsigma=0.0, interactive=0, usemask='auto-multithresh', \
                   sidelobethreshold=1.25, noisethreshold=5.0, \
                   lownoisethreshold=2.0, negativethreshold=0.0, \
                   minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                   minpercentchange=1.0, fastnoise=False, restart=True, \
                   savemodel='none', calcres=False, calcpsf=False, \
                   parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.0878815938952,
                           'Integrated Flux':306.734229447,
                           'Beam - Major Axis':8.3559513092,
                           'Beam - Minor Axis':5.94227743149,
                           'Beam - Position Angle':69.9616699219}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mosaic'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 1.19889616966, [45, 38, 0, 4])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [52, 99, 0, 0]), \
                       (img+'.image', False, [52, 100, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
                                   "E2E6.1.00034.S_tclean.ms")

        test_dict["test_mosaic_cube"]['self.parallel'] = self.parallel
        test_dict["test_mosaic_cube"]['report'] = report
        immoments(imagename = img+'.image', moments = 8, outfile = \
                  img+'.image.moment8')
        imview(img+'.image.moment8', out = img+'.image.moment8.png')
        test_dict["test_mosaic_cube"]['images'] = [img+'.image.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mfs(self):
        '''Mosaic MFS imaging'''

        file_name = 'mosaic_mfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='NGC5363', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'0', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS '
               '13:56:07.2100 +005.15.17.200', stokes='I', specmode='mfs', \
               nchan=-1, outframe='LSRK', perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='NGC5363', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'1', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS '
               '13:56:07.2100 +005.15.17.200', stokes='I', specmode='mfs', \
               nchan=-1, outframe='LSRK', perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=True, restoringbeam='common', pbcor=True, \
               weighting='briggs', robust=0.5, npixels=0, niter=30000, \
               threshold='0.00723Jy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, savemodel='none', calcres=False, \
               calcpsf=False, parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.00253072066577,
                           'Integrated Flux':1.46365389783,
                           'Beam - Major Axis':16.9057273865,
                           'Beam - Minor Axis':9.88278484344,
                           'Beam - Position Angle':87.210609436}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mosaic'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 0.0363491363823, [63, 63, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [64, 112, 0, 0]), \
                       (img+'.image', False, [64, 113, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
                                   "E2E6.1.00020.S_tclean.ms")

        test_dict["test_mosaic_mfs"]['self.parallel'] = self.parallel
        test_dict["test_mosaic_mfs"]['report'] = report
        imview(img+'.image',out=img+'.image.png')
        test_dict["test_mosaic_mfs"]['images'] = [img+'.image.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mtmfs(self):
        '''Mosaic mtmfs imaging'''

        file_name = 'mosaic_mtmfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='NGC5363', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'0', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS '
               '13:56:07.2100 +005.15.17.200', stokes='I', specmode='mfs', \
               nchan=-1, outframe='LSRK', perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'E2E6.1.00020.S_tclean.ms'], field='NGC5363', \
               spw=['0:113.893653412~114.629981537GHz;114.8809581'
               '~115.758887787GHz'], antenna=['0,1,2,3,4,5,6,7,8'], \
               scan=['6,9'], intent='OBSERVE_TARGET#ON_SOURCE', \
               datacolumn='data', imagename=file_name+'1', \
               imsize=[126, 126], cell=['2arcsec'], phasecenter='ICRS '
               '13:56:07.2100 +005.15.17.200', stokes='I', specmode='mfs', \
               nchan=-1, outframe='LSRK', perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               restoration=True, restoringbeam='common', pbcor=True, \
               weighting='briggs', robust=0.5, npixels=0, niter=30000, \
               threshold='0.00723Jy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=1.25, \
               noisethreshold=5.0, lownoisethreshold=2.0, \
               negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, savemodel='none', calcres=False, \
               calcpsf=False, parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image.tt0')
        expected_values = {'RMS':0.0025180800487,
                           'Integrated Flux':1.46911711973,
                           'Beam - Major Axis':17.6737804413,
                           'Beam - Minor Axis':10.0601730347,
                           'Beam - Position Angle':86.6785964966}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mos_mtmfs'),
            # checks the peak flux value and location
            imgval = [(img+'.image.tt0', 0.0350823290646, [63, 63, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image.tt0', True, [64, 113, 0, 0]), \
                       (img+'.image.tt0', False, [64, 114, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "E2E6.1.00020.S_tclean.ms")

        test_dict["test_mosaic_mtmfs"]['self.parallel'] = self.parallel
        test_dict["test_mosaic_mtmfs"]['report'] = report
        imview(img+'.image.tt0',out = img+'.image.tt0.png')
        test_dict["test_mosaic_mtmfs"]['images'] = [img+'.image.tt0.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)

#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_cube_eph(self):
        '''Mosaic ephemeris cube imaging'''

        file_name = 'mosaic_cube_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='Venus', \
               spw=['0'], antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,'
               '16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,'
               '35,36,37,38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'0', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='cubesource', nchan=948, start=5, width=1, \
               perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
               mosweight=True, usepointing=False, pblimit=0.2, \
               deconvolver='hogbom', restoration=False, \
               restoringbeam='common', pbcor=False, weighting='briggs', \
               robust=0.5, npixels=0, niter=0, threshold='0.0mJy', \
               nsigma=0.0, interactive=0, usemask='auto-multithresh', \
               sidelobethreshold=2.0, noisethreshold=4.25, \
               lownoisethreshold=1.5, negativethreshold=15.0, \
               minbeamfrac=0.3, growiterations=50, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='Venus', \
               spw=['0'], antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,'
              '16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,'
              '35,36,37,38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
              intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
              imagename=file_name+'1', imsize=[480, 420], \
              cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
              specmode='cubesource', nchan=948, start=5, width=1, \
              perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
              mosweight=True, usepointing=False, pblimit=0.2, \
              deconvolver='hogbom', restoration=True, \
              restoringbeam='common', pbcor=True, weighting='briggs', \
              robust=0.5, npixels=0, niter=700000, threshold='0.0106Jy', \
              nsigma=0.0, interactive=0, usemask='auto-multithresh', \
              sidelobethreshold=2.0, noisethreshold=4.25, \
              lownoisethreshold=1.5, negativethreshold=15.0, \
              minbeamfrac=0.3, growiterations=50, dogrowprune=True, \
              minpercentchange=1.0, fastnoise=False, restart=True, \
              savemodel='none', calcres=False, calcpsf=False, parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.85061036,
                           'Integrated Flux':6958649.90461285,
                           'Beam - Major Axis':0.9308618304028009,
                           'Beam - Minor Axis':0.719328746862788,
                           'Beam - Position Angle':-88.12278704464944}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mosaic'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 2.42139792, [289, 216, 0, 357])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [209, 387, 0, 0]), \
                       (img+'.image', False, [209, 388, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                  ["2017.1.00750.T_tclean_exe1.ms", \
                                   "2017.1.00750.T_tclean_exe2.ms"])

        test_dict["test_mosaic_cube_eph"]['self.parallel'] = self.parallel
        test_dict["test_mosaic_cube_eph"]['report'] = report
        immoments(imagename = img+'.image', moments = 8, outfile = \
                  img+'.image.moment8')
        imview(img+'.image.moment8',out = img+'.image.moment8.png')
        test_dict["test_mosaic_cube_eph"]['images'] = \
            [img+'.image.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mfs_eph(self):
        '''Mosaic ephemeris mfs imaging'''

        file_name = 'mosaic_mfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='Venus', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'0', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto' \
               '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
               lownoisethreshold=1.5, negativethreshold=0.0, \
               minbeamfrac=0.3, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='Venus', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'1', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='hogbom', \
               restoration=True, restoringbeam='common', pbcor=True, \
               weighting='briggs', robust=0.5, npixels=0, niter=7000000, \
               threshold='0.0316Jy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=2.0, \
               noisethreshold=4.25, lownoisethreshold=1.5, \
               negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, calcres=False, calcpsf=False, parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image')
        expected_values = {'RMS':0.849691489705,
                           'Integrated Flux':7342.75576754,
                           'Beam - Major Axis':0.929236471653,
                           'Beam - Minor Axis':0.717731595039,
                           'Beam - Position Angle':-88.23828125}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mosaic'),
            # checks the peak flux value and location
            imgval = [(img+'.image', 2.36644268036, [153, 151, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image', True, [209, 387, 0, 0]), \
                       (img+'.image', False, [209, 388, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
                                   "2018.1.00879.S_tclean.ms")

        test_dict["test_mosaic_mfs_eph"]['self.parallel'] = self.parallel
        test_dict["test_mosaic_mfs_eph"]['report'] = report
        imview(img+'.image',out = img+'.image.png')
        test_dict["test_mosaic_mfs_eph"]['images'] = [img+'.image.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mtmfs_eph(self):
        '''Mosaic ephemeris mtmfs imaging'''

        file_name = 'mosaic_mtmfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='Venus', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'0', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               restoration=False, restoringbeam='common', pbcor=False, \
               weighting='briggs', robust=0.5, npixels=0, niter=0, \
               threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto' \
               '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
               lownoisethreshold=1.5, negativethreshold=0.0, \
               minbeamfrac=0.3, growiterations=75, dogrowprune=True, \
               minpercentchange=1.0, fastnoise=False, savemodel='none', \
               parallel=False)

        # move files to iter1
        print('Renaming iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.rename(file_name, file_list)
        time.sleep(5)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=[data_path+'2018.1.00879.S_tclean.ms'], field='Venus', \
               spw=['0:261.752937691~261.774177925GHz;261.783699409'
               '~261.837898628GHz;261.958504097~261.984871284GHz'], \
               antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,'
               '19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,'
               '38,39,40,41,42,43,44,45,46'], scan=['7,11'], \
               intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
               imagename=file_name+'1', imsize=[480, 420], \
               cell=['0.14arcsec'], phasecenter='TRACKFIELD', \
               specmode='mfs', nchan=-1, perchanweightdensity=False, \
               gridder='mosaic', chanchunks=-1, mosweight=True, \
               usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
               restoration=True, restoringbeam='common', pbcor=True, \
               weighting='briggs', robust=0.5, npixels=0, niter=7000000, \
               threshold='0.0316Jy', nsigma=0.0, interactive=0, \
               usemask='auto-multithresh', sidelobethreshold=2.0, \
               noisethreshold=4.25, lownoisethreshold=1.5, \
               negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
               dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
               restart=True, calcres=False, calcpsf=False, parallel=False)

        img = os.getcwd()+'/'+file_name+'1'

        # restoring image values retrieve
        rms, int_flux, major, minor, pa = self.image_stats(img,'.image.tt0')
        expected_values = {'RMS':0.849691489705,
                           'Integrated Flux':7342.75576754,
                           'Beam - Major Axis':0.929236471653,
                           'Beam - Minor Axis':0.717731595039,
                           'Beam - Position Angle':-88.23828125}

        report1 = casaTestHelper.checkall( \
            # checks if the expected files are created
            imgexist = self.image_list(img,'mos_mtmfs'),
            # checks the peak flux value and location
            imgval = [(img+'.image.tt0', 2.36644268036, [153, 151, 0, 0])],
            # checks that a pixel value at edge of pb is not flagged
            # checks that a pixel value at edge of pb is flagged
            imgmask = [(img+'.image.tt0', True, [209, 387, 0, 0]), \
                       (img+'.image.tt0', False, [209, 388, 0, 0])])

        out, report2 = casaTestHelper.check_val( \
            rms, expected_values['RMS'], valname = 'RMS', epsilon = \
            abs(expected_values['RMS']*self.epsilon))
        out, report3 = casaTestHelper.check_val( \
            int_flux, expected_values['Integrated Flux'], valname = \
            'Integrated Flux', epsilon = \
            abs(expected_values['Integrated Flux']*self.epsilon))
        out, report4 = casaTestHelper.check_val( \
            major, expected_values['Beam - Major Axis'], valname = \
            'Beam - Major Axis', epsilon = \
            expected_values['Beam - Major Axis']*self.epsilon)
        out, report5 = casaTestHelper.check_val( \
            minor, expected_values['Beam - Minor Axis'], valname = \
            'Beam - Minor Axis', epsilon = \
            expected_values['Beam - Minor Axis']*self.epsilon)
        out, report6 = casaTestHelper.check_val( \
            pa, expected_values['Beam - Position Angle'], valname = \
            'Beam - Position Angle', epsilon = \
            abs(expected_values['Beam - Position Angle'])*self.epsilon)

        report = report1 + report2 + report3 + report4 + report5 + report6

        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
                                   "2018.1.00879.S_tclean.ms")

        test_dict["test_mosaic_mtmfs_eph"]['self.parallel'] = self.parallel
        test_dict["test_mosaic_mtmfs_eph"]['report'] = report
        imview(img+'.image.tt0',out = img+'.image.tt0.png')
        test_dict["test_mosaic_mtmfs_eph"]['images'] = [img+'.image.tt0.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
                        msg = report)


def suite():
     return [Test_standard, Test_mosaic]

# Main #
if __name__ == '__main__':
    unittest.main()



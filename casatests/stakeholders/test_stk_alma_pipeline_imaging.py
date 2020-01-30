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
import scipy
from parallel.parallel_task_helper import ParallelTaskHelper

# export PYTHONPATH=/lustre/naasc/sciops/comm/sbooth/CASA_ALMA_pipeline/tclean_test_dir/casaTestHelper
import casaTestHelper
data_path = '/lustre/naasc/sciops/comm/sbooth/CASA_ALMA_pipeline/data_dir/'
#data_path = os.environ.get('CASAPATH').split()[0] + '/casa-data-vt/vlass/'


## Base Test class with Utility functions
class test_tclean_base(unittest.TestCase):

    def setUp(self):
        self._myia = iatool()
#        self.epsilon = 0.01 # sets epsilon as a percentage (1%)
        self.msfile = ""
#        self.img = "tst"
        # To use subdir in the output image names in some tests (CAS-10937)
#        self.img_subdir = 'refimager_tst_subdir'
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

    def tearDown(self):
        casaTestHelper.generate_weblog("tclean_ALMA_pipeline",test_dict)
        self._myia.done()
        """ don't delete it all """
#        self.delData()

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self, msname=[""]):
        if msname != [""]:
             self.msfile=msname

    def delData(self, msname=[""]):
        if msname != [""]:
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

    def getNameDoc(self):
#        print inspect.stack()[1]
        testname=inspect.stack()[1][3]
        print("Test name  : " + testname)
#        tname = inspect.getframeinfo(inspect.currentframe()).function
        doc = eval('self.'+testname + '.__doc__')
        print("Doc : " +  doc)
        return testname, doc

    # function to copy iter0 images to iter1 images
    def copy_images(self, file_name):
        file_list = glob.glob(file_name+'0*')
        os.popen('rm -rf '+file_name+'1*')
        for item in file_list:
            sp_item = item.split('iter0')
            shutil.copytree(item, sp_item[0]+'iter1'+sp_item[1])
        if os.path.isdir(file_name+'1.workdirectory/'):
            work0_list = glob.glob(file_name+'0.workdirectory/*')
            for item in work0_list:
                sp_item = item.split('iter0')
                new_item = sp_item[0]+'iter1'+sp_item[1]+'iter1'+sp_item[2]
                shutil.copytree(item, new_item)

    # function to return per-channel beam statistics
    # will be deprecated and combined into image_stats 
    # once CASA beam issue is fixed
    def cube_beam_stats(self, img):
        self._myia.open(img)

        res_bmin_dict = {}; res_bmaj_dict = {}; res_pa_dict = {}
        beam_dict = self._myia.restoringbeam()['beams']
        for item in beam_dict.keys():
            res_bmin_dict[item] = beam_dict[item]['*0']['minor']['value']
            res_bmaj_dict[item] = beam_dict[item]['*0']['major']['value']
            res_pa_dict[item] = \
                beam_dict[item]['*0']['positionangle']['value']

        self._myia.close()

        return res_bmin_dict, res_bmaj_dict, res_pa_dict

    # function that takes an image file and returns a statistics dictionary
    def image_stats(self, img, suf, region_file=None):
        self._myia.open(img+suf)

        im_size = self._myia.boundingbox()['imageShape'].tolist()
        npts = im_size[0]*im_size[1]*im_size[3]
        start = self._myia.summary()['refval'][3]
        end = start + (self._myia.summary()['incr'][3] * \
                       self._myia.summary()['shape'][3])

        if '.image' in suf:
            stats_dict = {'com_bmin': \
                self._myia.commonbeam()['minor']['value'], 'com_bmaj': \
                self._myia.commonbeam()['major']['value'], 'com_pa': \
                self._myia.commonbeam()['pa']['value'], 'npts': npts, \
                'npts_unmasked': self._myia.statistics()['npts'][0], \
                'npts_real': numpy.count_nonzero(~numpy.isnan( \
                self._myia.getchunk())), 'freq_bin': \
                self._myia.summary()['incr'][3], 'start': start, 'end': end, \
                'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist(), 'im_rms': \
                self._myia.statistics()['rms'][0],'rms_per_chan': 1, \
                'im_sum': self._myia.statistics()['sum'][0], 'regn_sum': \
                self._myia.statistics(region=region_file)['sum'][0], \
                'im_fit':0}
            if 'mosaic' in  img:
                stats_dict['rms_per_field'] = 1

#                stats_dict['rms_per_chan']
#                stats_dict['rms_per_field']

#                fit_dict = self._myia.fitcomponents( \
#                    region=region_file)['results']['component0']
#                stats_dict['im_fit'] = [fit_dict['pixelcoords'].tolist(), \
#                    fit_dict['peak']['value'], fit_dict, fit_dict]

        if '.mask' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'mask_pix': \
                numpy.count_nonzero(self._myia.getchunk()), 'mask_regns': \
                scipy.ndimage.label(self._myia.getchunk())[1]}

        if '.pb' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist(),'npts_0.2': \
                numpy.count_nonzero(self._myia.getchunk()>0.2), 'npts_0.5': \
                numpy.count_nonzero(self._myia.getchunk()>0.5), 'pb_im_fit':0}

#                stats_dict['pb_im_fit']: \
#                    self._myia.fitcomponents(region=region_file)}

        if '.psf' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist(), 'im_rms': \
                self._myia.statistics()['rms'][0], 'im_sum': \
                self._myia.statistics()['sum'][0], 'cen_im_fit':0}

#                stats_dict['cen_im_fit']: \
#                    self._myia.fitcomponents(region=region_file)}

        if '.residual' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist(), 'im_rms': \
                self._myia.statistics()['rms'][0], 'im_sum': \
                self._myia.statistics()['sum'][0], 'regn_sum': \
                self._myia.statistics(region=region_file)['sum'][0]}

        if '.model' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist(), 'im_rms': \
                self._myia.statistics()['rms'][0], 'im_sum': \
                self._myia.statistics()['sum'][0], 'regn_sum': \
                self._myia.statistics(region=region_file)['sum'][0], \
                'mask_non0': numpy.count_nonzero(self._myia.getchunk())}

        if '.sumwt' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist()}

        if '.weight' in suf:
            stats_dict = {'npts': npts, 'npts_unmasked': \
                self._myia.statistics()['npts'][0], 'npts_real': \
                numpy.count_nonzero(~numpy.isnan(self._myia.getchunk())), \
                'freq_bin': self._myia.summary()['incr'][3], 'start': start, \
                'end': end, 'nchan': im_size[3], 'max_val': \
                self._myia.statistics()['max'][0], 'max_val_pos': \
                self._myia.statistics()['maxpos'].tolist(), 'min_val': \
                self._myia.statistics()['min'][0], 'min_val_pos': \
                self._myia.statistics()['minpos'].tolist(), 'im_rms': \
                self._myia.statistics()['rms'][0], 'im_sum': \
                self._myia.statistics()['sum'][0], 'npts_0.2': \
                numpy.count_nonzero(self._myia.getchunk()>0.2), 'npts_0.3': \
                numpy.count_nonzero(self._myia.getchunk()>0.3)}

        self._myia.close()

        return stats_dict

    # function used to return expected imaging output files
    def image_list(self, img, mode):
        standard = [img+'.psf', img+'.residual', img+'.image', \
            img+'.image.pbcor', img+'.mask', img+'.pb', img+'.model', \
            img+'.sumwt']
        mosaic = [img+'.weight']
        mtmfs = [img+'.alpha', img+'.alpha.error', img+'.alpha.pbcor', \
           img+'.psf.tt0', img+'.psf.tt1', img+'.psf.tt2', \
           img+'.residual.tt0', img+'.residual.tt1', img+'.image.tt0',\
           img+'.image.tt1', img+'.image.tt0.pbcor', img+'.image.tt1.pbcor', \
           img+'.mask', img+'.pb.tt0', img+'.model.tt0', img+'.model.tt1', \
           img+'.sumwt.tt0', img+'.sumwt.tt1', img+'.sumwt.tt2']
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

    # function that takes and image and turns it into a .png for weblog
    def png_creator(self, img, range_list):
        immoments(imagename = img, moments = 8, outfile = img+'.moment8')
        imview(raster={'file': img+'.moment8', 'range': range_list}, \
            out = {'file': img+'.moment8.png'})
        os.popen('mogrify -trim '+img+'.moment8.png')


##############################################
##############################################
test_dict = {}
class Test_standard(test_tclean_base):


    @casaTestHelper.stats_dict(test_dict)
    def test_standard_cube(self):
        '''
        Standard (single field) cube imaging
        central field of SMIDGE_NWCloud (field 3), spw 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_cube.iter'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, imagename=file_name+'0', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS'
            ' 00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526744GHz', width='0.244174087287MHz',  \
            outframe='LSRK', pblimit=0.2, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            deconvolver='hogbom', usepointing=False, restoration=False, \
            pbcor=False, weighting='briggs', restoringbeam='common', \
            robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        if self.parallel:
            tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
                spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
                scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
                datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
                '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
                nchan=508, start='220.2526744GHz', width='0.244174087287MHz',\
                outframe='LSRK', perchanweightdensity=False, \
                usepointing=False, pblimit=0.2, nsigma=0.0, \
                gridder='standard', chanchunks=-1, mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggs', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none', \
                parallel=True)

            # retrieve per-channel beam statistics (only in parallel)
            img = os.getcwd()+'/'+file_name+'1'
            res_bmin_dict, res_bmaj_dict, res_pa_dict = \
                self.cube_beam_stats(img+'.image')

            tclean(vis=self.msfile, imagename=file_name+'1', spw=['0'], \
                field='1', imsize=[80, 80], cell=['1.1arcsec'], \
                phasecenter='ICRS 00:45:54.3836 -073.15.29.413', stokes='I', \
                antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                specmode='cube', nchan=508, start='220.2526744GHz', \
                width='0.244174087287MHz', outframe='LSRK', \
                perchanweightdensity=False, gridder='standard', \
                chanchunks=-1, mosweight=False, usepointing=False, \
                pblimit=0.2, deconvolver='hogbom', restoration=True, \
                restoringbeam='common', pbcor=True, weighting='briggs', \
                robust=0.5, npixels=0, niter=0, threshold='0.354Jy', \
                nsigma=0.0, interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none',  \
                parallel=False)

        else:
            tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
                spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
                scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
                datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
                '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
                nchan=508, start='220.2526744GHz', width='0.244174087287MHz',\
                outframe='LSRK', perchanweightdensity=False, \
                usepointing=False, pblimit=0.2, nsigma=0.0, \
                gridder='standard', chanchunks=-1, mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggs', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none', \
                restoringbeam='common', parallel=False)

            img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/standard_cube.image.crtf')

        exp_im_stats = {'com_bmaj': 8.509892605313942,
            'com_bmin': 5.950050676606115,
            'com_pa': 72.54607919421503,
            'npts': 3251200,
            'npts_unmasked': 1522476.0,
            'freq_bin': 244174.08728027344,
            'start': 220252674399.99997,
            'end': 0,
            'nchan': 508,
            'max_val': 0.94608676433563232,
            'max_val_pos':[38, 36, 0, 254],
            'min_val': -0.70467984676361084,
            'min_val_pos':[18, 57, 0, 374],
            'im_rms': 0.143986095161,
            'rms_per_chan': 1,
            'im_sum': 168.242297979,
            'regn_sum': 72.3549563158,
            'npts_real': 3251200}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 70, 0, 0]), \
                (img+'.image', False, [40, 71, 0, 0]), \
                (img+'.image', True, [10, 40, 0, 0]), \
                (img+'.image', False, [9, 40, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 3251200,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'mask_pix': 437,
            'mask_regns': 1,
            'npts_real': 3251200}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/standard_cube.pb.crtf')

        exp_pb_stats = {'npts': 3251200,
            'npts_unmasked': 1522476.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.200896695256,
            'min_val_pos':[25, 13, 0, 396],
            'npts_0.2': 1522476,
            'npts_0.5': 736092,
            'npts_real': 3251200}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/standard_cube.psf.crtf')

        exp_psf_stats = {'npts': 3251200,
            'npts_unmasked': 3251200.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.218764916062,
            'min_val_pos':[1, 16, 0, 503],
            'im_rms':  0.136036099793,
            'im_sum': 7472.57665916,
            'npts': 3251200}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/standard_cube.residual.crtf')

        exp_resid_stats = {'npts': 3251200,
            'npts_unmasked': 1522476.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 0.785612404346,
            'max_val_pos':[42, 43, 0, 256],
            'min_val': -0.704679846764,
            'min_val_pos':[18, 57, 0, 374],
            'im_rms': 0.143918523224,
            'im_sum': 124.317946204,
            'npts': 3251200}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/standard_cube.model.crtf')

        exp_model_stats = {'npts': 3251200,
            'npts_unmasked': 3251200.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 0.286023736,
            'max_val_pos':[38, 36, 0, 254],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.000249846621096,
            'im_sum': 0.92636379227,
            'regn_sum': 0.92636379227,
            'mask_non0': 6,
            'npts': 3251200}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 508,
            'npts_unmasked': 508.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 94.4766769409,
            'max_val_pos':[0, 0, 0, 17],
            'min_val': 94.4766464233,
            'min_val_pos':[0, 0, 0, 449],
            'npts': 508}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[0.3, 1.0])
        self.png_creator(img+'.residual', range_list=[0.3, 1.0])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mfs(self):
        '''
        Standard (single field) MFS imaging
        central field of NGC5363 (field 2), spw 16 & 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_mfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='2', spw=['0:113.893653412~'
            '114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', stokes='I', imsize=[80, 80], \
            cell=['2arcsec'], phasecenter='ICRS 13:56:07.2100 '
            '+005.15.17.200', specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', nterms=2, restoration=False, \
            restoringbeam='common', pbcor=False, weighting='briggs', \
            robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='2',  spw=['0:113.893653412~'
            '114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', stokes='I', imsize=[80, 80], \
            cell=['2arcsec'], phasecenter='ICRS 13:56:07.2100 '
            '+005.15.17.200', specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', nterms=2, restoration=True, \
            restoringbeam='common', pbcor=True, weighting='briggs', \
            robust=0.5, npixels=0, niter=30000, threshold='0.00723Jy', \
            nsigma=0.0, interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, restart=True, calcres=False, calcpsf=False, \
            savemodel='none', parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/standard_mfs.image.crtf')

        exp_im_stats = {'com_bmaj': 18.0537223816,
            'com_bmin': 10.3130550385,
            'com_pa': 86.4389877319,
            'npts': 6400,
            'npts_unmasked': 3793.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0362096913159,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00253091799095,
            'min_val_pos':[51, 36, 0, 0],
            'im_rms': 0.00317099603729,
            'rms_per_chan': 1,
            'im_sum': 1.72629857491,
            'regn_sum': 1.70481428877,
            'npts_real': 6400}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 74, 0, 0]), \
                (img+'.image', False, [40, 75, 0, 0]), \
                (img+'.image', True, [6, 40, 0, 0]), \
                (img+'.image', False, [5, 40, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 6400,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 334,
            'mask_regns': 1,
            'npts_real': 6400}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/standard_mfs.pb.crtf')

        exp_pb_stats = {'npts': 6400,
            'npts_unmasked': 3793.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.200896695256,
            'min_val_pos':[33, 6, 0, 0],
            'npts_0.2': 3793,
            'npts_0.5': 1813,
            'npts_real': 6400}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/standard_mfs.psf.crtf')

        exp_psf_stats = {'npts': 6400,
            'npts_unmasked': 6400.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.170877560973,
            'min_val_pos':[45, 14, 0, 0],
            'im_rms':  0.11175375188,
            'im_sum': 13.431971317,
            'npts_real': 6400}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/standard_mfs.residual.crtf')

        exp_resid_stats = {'npts': 6400,
            'npts_unmasked': 3793.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.00680132163689,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00267858314328,
            'min_val_pos':[51, 36, 0, 0],
            'im_rms': 0.00119958583186,
            'im_sum': 0.167714220393,
            'npts_real': 6400}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/standard_mfs.model.crtf')

        exp_model_stats = {'npts': 6400,
            'npts_unmasked': 6400.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0253213103861,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.000320901758782,
            'im_sum': 0.029550973326,
            'regn_sum': 0.029550973326,
            'mask_non0': 2,
            'npts_real': 6400}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 3208318.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 3208318.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[-0.003, 0.04])
        self.png_creator(img+'.residual', range_list=[-0.003, 0.04])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mtmfs(self):
        '''
        Single field mtmfs imaging
        central field of NGC5363 (field 2), spw 16 & 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_mtmfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='2', spw=['0:113.893653412'
            '~114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', stokes='I', imsize=[80, 80], \
            cell=['2arcsec'], phasecenter='ICRS 13:56:07.2100 '
            '+005.15.17.200', specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='mtmfs', nterms=2, restoration=False, \
            restoringbeam='common', pbcor=False, weighting='briggs', \
            robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='2', spw=['0:113.893653412'
            '~114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', stokes='I', imsize=[80, 80], \
            cell=['2arcsec'], phasecenter='ICRS 13:56:07.2100 '
            '+005.15.17.200', specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='mtmfs', nterms=2, restoration=True, \
            restoringbeam='common', pbcor=True, weighting='briggs', \
            robust=0.5, npixels=0, niter=30000, threshold='0.00723Jy', \
            nsigma=0.0, interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, restart=True, calcres=False, calcpsf=False, \
            savemodel='none', parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.image.tt0.crtf')

        exp_im_stats = {'com_bmaj': 18.0537223816,
            'com_bmin': 10.3130550385,
            'com_pa': 86.4389877319,
            'npts': 6400,
            'npts_unmasked': 3793.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0375871881843,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00482273567468,
            'min_val_pos':[40, 62, 0, 0],
            'im_rms': 0.00341863379853,
            'rms_per_chan': 1,
            'im_sum': 1.82207681158,
            'regn_sum': 1.77460362448,
            'npts_real': 6400}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [40, 74, 0, 0]), \
                (img+'.image.tt0', False, [40, 75, 0, 0]), \
                (img+'.image.tt0', True, [6, 40, 0, 0]), \
                (img+'.image.tt0', False, [5, 40, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 6400,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 332,
            'mask_regns': 1,
            'npts_real': 6400}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.pb.tt0.crtf')

        exp_pb_stats = {'npts': 6400,
            'npts_unmasked': 3793.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.200896695256,
            'min_val_pos':[33, 6, 0, 0],
            'npts_0.2': 3793,
            'npts_0.5': 1813,
            'npts_real': 6400}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.psf.tt0.crtf')

        exp_psf_stats = {'npts': 6400,
            'npts_unmasked': 6400.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.170877560973,
            'min_val_pos':[45, 14, 0, 0],
            'im_rms':  0.11175375188,
            'im_sum': 13.431971317,
            'npts_real': 6400}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual.tt0', \
            region_file = data_path+
            'region_files/standard_mtmfs.residual.tt0.crtf')

        exp_resid_stats = {'npts': 6400,
            'npts_unmasked': 3793.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.00681467168033,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00269496999681,
            'min_val_pos':[51, 36, 0, 0],
            'im_rms': 0.00121173798212,
            'im_sum': 0.146451243258,
            'npts_real': 6400}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.model.tt0.crtf')

        exp_model_stats = {'npts': 6400,
            'npts_unmasked': 6400.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0257838983089,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.00032823232983,
            'im_sum': 0.0307542048395,
            'regn_sum': 0.0307542048395,
            'mask_non0': 2,
            'npts_real': 6400}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt.tt0')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 15849921197.895538,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 3208318.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 3208318.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image.tt0', range_list=[-0.005, 0.04])
        self.png_creator(img+'.residual.tt0', range_list=[-0.005, 0.04])

        test_dict[testname]['images'] = \
            [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_cube_eph(self):
        '''
        Single field multi-EB ephemeris cube imaging
        field 21PGiacobini-Zinner, spw 20
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_cube_eph.iter'
        self.prepData([data_path+'2017.1.00750.T_tclean_exe1.ms', 
            data_path+'2017.1.00750.T_tclean_exe2.ms'])

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='21PGiacobini-Zinner', spw=['0', '0'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11', '0,1,2,3,4,5,6,7,8,9'], \
            scan=['7,11,15,19,23','8,12,16,20,24'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[80, 80], cell=['0.66arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='cubesource', \
            nchan=1000, start=1550, width=1, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            usepointing=False, pblimit=0.2, deconvolver='hogbom', \
            restoration=False, restoringbeam='common', pbcor=False, \
            weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, savemodel='none', \
            parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='21PGiacobini-Zinner', spw=['0', '0'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11', '0,1,2,3,4,5,6,7,8,9'], \
            scan=['7,11,15,19,23', '8,12,16,20,24'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[80, 80], cell=['0.66arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='cubesource', \
            nchan=1000, start=1550, width=1, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            usepointing=False, pblimit=0.2, deconvolver='hogbom', \
            restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggs', robust=0.5, npixels=0, niter=30000, \
            threshold='0.274Jy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            calcres=False, calcpsf=False, savemodel='none', \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/standard_cube_eph.image.crtf')

        exp_im_stats = {'com_bmaj': 4.50223609831,
            'com_bmin': 3.32703591993,
            'com_pa': 87.0952374381,
            'npts': 6400000,
            'npts_unmasked': 3233000.0,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'max_val': 3.06201076508,
            'max_val_pos':[46, 41, 0, 491],
            'min_val': -0.362504482269,
            'min_val_pos':[59, 39, 0, 494],
            'im_rms': 0.055973852238,
            'rms_per_chan': 1,
            'im_sum': 2249.5560029,
            'regn_sum': 253.368211955,
            'npts_real': 6400000}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 6400000,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'mask_pix': 9932,
            'mask_regns': 1,
            'npts_real': 6400000}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/standard_cube_eph.pb.crtf')

        exp_pb_stats = {'npts': 6400000,
            'npts_unmasked': 3233000.0,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.20036059618,
            'min_val_pos':[39, 8, 0, 810],
            'npts_0.2': 3233000,
            'npts_0.5': 1549000,
            'npts_real': 6400000}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/standard_cube_eph.psf.crtf')

        exp_psf_stats = {'npts': 6400000,
            'npts_unmasked': 6400000.0,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.164189130068,
            'min_val_pos':[36, 35, 0, 993],
            'im_rms':  0.0871161935921,
            'im_sum': 2742.74484326,
            'npts_real': 6400000}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/standard_cube_eph.residual.crtf')

        exp_resid_stats = {'npts': 6400000,
            'npts_unmasked': 3233000.0,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'max_val': 0.366728395224,
            'max_val_pos':[39, 68, 0, 502],
            'min_val': -0.338401287794,
            'min_val_pos':[32, 47, 0, 493],
            'im_rms': 0.0469751110358,
            'im_sum': 242.908819752,
            'npts_real': 6400000}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/standard_cube_eph.model.crtf')

        exp_model_stats = {'npts': 6400000,
            'npts_unmasked': 6400000.0,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'max_val': 1.27968358994,
            'max_val_pos':[46, 41, 0, 490],
            'min_val': -0.0999799370766,
            'min_val_pos':[58, 39, 0, 494],
            'im_rms': 0.00170565742213,
            'im_sum': 51.6511641908,
            'regn_sum': 5.50620770082,
            'mask_non0': 417,
            'npts_real': 6400000}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 1000,
            'npts_unmasked': 1000.0,
            'freq_bin': 122071.64398193359,
            'start': 0,
            'end': 0,
            'nchan': 1000,
            'max_val': 1009.50134277,
            'max_val_pos':[0, 0, 0, 980],
            'min_val': 1007.94287109,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1000}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[0.0, 3.25])
        self.png_creator(img+'.residual', range_list=[0.0, 3.25])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png', img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mfs_eph(self):
        '''
        Standard (single field) ephemeris mfs imaging
        central field of Venus (field 2), spw 25 & 45
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_mfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='2', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'],  \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[288, 288], cell=['0.14arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='mfs', nchan=-1, \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='2', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[288, 288], cell=['0.14arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='mfs', nchan=-1, \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=7000000, threshold='0.0316Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=2.0, \
            noisethreshold=4.25, lownoisethreshold=1.5, \
            negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, calcres=False, calcpsf=False, \
            savemodel='none', parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/standard_mfs_eph.image.crtf')

        exp_im_stats = {'com_bmaj': 0.875921189785,
            'com_bmin': 0.673674583435,
            'com_pa': 88.5387649536,
            'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.03113353252,
            'max_val_pos':[224, 153, 0, 0],
            'min_val': -1.01794064045,
            'min_val_pos':[222, 93, 0, 0],
            'im_rms': 0.359352011299,
            'rms_per_chan': 1,
            'im_sum': -1491.198136,
            'regn_sum': 3362.95355159,
            'npts_real': 82944}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [144, 266, 0, 0]), \
                (img+'.image', False, [144, 267, 0, 0]), \
                (img+'.image', True, [22, 145, 0, 0]), \
                (img+'.image', False, [21, 145, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 82944,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 82944}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/standard_mfs_eph.pb.crtf')

        exp_pb_stats = {'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': 0.200061768293,
            'min_val_pos':[114, 25, 0, 0],
            'npts_0.2': 47325,
            'npts_0.5': 22362,
            'npts_real': 82944}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/standard_mfs_eph.psf.crtf')

        exp_psf_stats = {'npts': 82944,
            'npts_unmasked': 82944.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': -0.0609973333776,
            'min_val_pos':[140, 137, 0, 0],
            'im_rms':  0.0198315591613,
            'im_sum': 16.3422700718,
            'npts_real': 82944}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/standard_mfs_eph.residual.crtf')

        exp_resid_stats = {'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.03113353252,
            'max_val_pos':[224, 93, 0, 0],
            'min_val': -0.676309525967,
            'min_val_pos':[54, 126, 0, 0],
            'im_rms': 0.359352011299,
            'im_sum': -1491.198136,
            'npts_real': 82944}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/standard_mfs_eph.model.crtf')

        exp_model_stats = {'npts': 82944,
            'npts_unmasked': 82944.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.0,
            'im_sum': 0.0,
            'regn_sum': 0.0,
            'mask_non0': 0,
            'npts_real': 82944}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 23234590.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 23234590.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7

        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[-1.05, 1.05])
        self.png_creator(img+'.residual', range_list=[-1.05, 1.05])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_mtmfs_eph(self):
        '''
        Standard (single field) ephemeris mtmfs imaging
        central field of Venus (field 2), spw 25 & 45
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_mtmfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='2', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[288, 288], nterms=2, \
            cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
            specmode='mfs', nchan=-1, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
            restoration=False, restoringbeam='common', pbcor=False, \
            weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='2', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[288, 288], nterms=2, \
            cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
            specmode='mfs', nchan=-1, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
            restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggs', robust=0.5, npixels=0, niter=7000000, \
            threshold='0.0316Jy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, restart=True, calcres=False, calcpsf=False, \
            savemodel='none', parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.image.tt0.crtf')

        exp_im_stats = {'com_bmaj': 0.875921189785,
            'com_bmin': 0.673674583435,
            'com_pa': 88.5387649536,
            'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.01515209675,
            'max_val_pos':[224, 93, 0, 0],
            'min_val': -1.01363265514,
            'min_val_pos':[222, 93, 0, 0],
            'im_rms': 0.361662749039,
            'rms_per_chan': 1,
            'im_sum': -1558.29933771,
            'regn_sum': 3388.57268598,
            'npts_real': 82944}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [144, 266, 0, 0]), \
                (img+'.image.tt0', False, [144, 267, 0, 0]), \
                (img+'.image.tt0', True, [22, 145, 0, 0]), \
                (img+'.image.tt0', False, [21, 145, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 82944,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 82944}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.pb.tt0.crtf')

        exp_pb_stats = {'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': 0.200061768293,
            'min_val_pos':[114, 25, 0, 0],
            'npts_0.2': 47325,
            'npts_0.5': 22362,
            'npts_real': 82944}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.psf.tt0.crtf')

        exp_psf_stats = {'npts': 82944,
            'npts_unmasked': 82944.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': -0.0609973333776,
            'min_val_pos':[140, 137, 0, 0],
            'im_rms':  0.0198315591613,
            'im_sum': 16.3422700718,
            'npts_real': 82944}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual.tt0', \
            region_file = data_path+
            'region_files/standard_mtmfs_eph.residual.tt0.crtf')

        exp_resid_stats = {'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.03113353252,
            'max_val_pos':[224, 153, 0, 0],
            'min_val': -1.01794064045,
            'min_val_pos':[222, 93, 0, 0],
            'im_rms': 0.359352011299,
            'im_sum': -1491.198136,
            'npts_real': 82944}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.model.tt0.crtf')

        exp_model_stats = {'npts': 82944,
            'npts_unmasked': 82944.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.0,
            'im_sum': 0.0,
            'regn_sum': 0.0,
            'mask_non0': 0,
            'npts_real': 82944}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt.tt0')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 16762501225.396851,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 23234590.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 23234590.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image.tt0', range_list=[-1.05, 1.05])
        self.png_creator(img+'.residual.tt0', range_list=[-1.05, 1.05])

        test_dict[testname]['images'] = \
            [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_standard_cal(self):
        '''
        Calibrator image
        field J2258-2758, spw 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'standard_cal.iter'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='0', spw=['0'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['3'], \
            intent='CALIBRATE_BANDPASS#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[90, 90], stokes='I',  \
            cell=['0.85arcsec'], phasecenter='ICRS 22:58:05.9629 '
            '-027.58.21.257', specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.5, noisethreshold=6.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, savemodel='none', \
            parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='0', spw=['0'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['3'], \
            intent='CALIBRATE_BANDPASS#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[90, 90], stokes='I', \
            cell=['0.85arcsec'], phasecenter='ICRS 22:58:05.9629 '
            '-027.58.21.257', specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='standard', chanchunks=-1, \
            mosweight=False, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=300000, threshold='0.0241Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=1.5, \
            noisethreshold=6.0, lownoisethreshold=2.0, \
            negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, calcres=False, calcpsf=False, savemodel='none', \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/standard_cal.image.crtf')

        exp_im_stats = {'com_bmaj': 9.98560905457,
            'com_bmin': 4.6246509552,
            'com_pa': -86.3877105713,
            'npts': 8100,
            'npts_unmasked': 5041.0,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 2.40606927872,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': -0.0115160318092,
            'min_val_pos':[50, 16, 0, 0],
            'im_rms': 0.203888800435,
            'rms_per_chan': 1,
            'im_sum': 172.688482511,
            'regn_sum': 171.572869264,
            'npts_real': 8100}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [45, 85, 0, 0]), \
                (img+'.image', False, [45, 86, 0, 0]), \
                (img+'.image', True, [5, 45, 0, 0]), \
                (img+'.image', False, [4, 45, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 8100,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 404,
            'mask_regns': 1,
            'npts_real': 8100}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/standard_cal.pb.crtf')

        exp_pb_stats = {'npts': 8100,
            'npts_unmasked': 5041.0,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': 0.200092822313,
            'min_val_pos':[36, 6, 0, 0],
            'npts_0.2': 5041,
            'npts_0.5': 2409,
            'npts_real': 8100}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/standard_cal.psf.crtf')

        exp_psf_stats = {'npts': 8100,
            'npts_unmasked': 8100.0,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': -0.186672970653,
            'min_val_pos':[31, 41, 0, 0],
            'im_rms':  0.125651854223,
            'im_sum': 40.2140428556,
            'npts_real': 8100}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/standard_cal.residual.crtf')

        exp_resid_stats = {'npts': 8100,
            'npts_unmasked': 5041.0,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0233333036304,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': -0.0115160560235,
            'min_val_pos':[50, 16, 0, 0],
            'im_rms': 0.00410021113948,
            'im_sum': 0.122163831366,
            'npts_real': 8100}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/standard_cal.model.crtf')

        exp_model_stats = {'npts': 8100,
            'npts_unmasked': 8100.0,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 2.38273620605,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.0264748467339,
            'im_sum': 2.38273620605,
            'regn_sum': 2.38273620605,
            'mask_non0': 1,
            'npts_real': 8100}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 125009872.91876221,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 201537.8125,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 201537.8125,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[-0.015, 2.5])
        self.png_creator(img+'.residual', range_list=[-0.015, 2.5])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


###############################################
###############################################

class Test_mosaic(test_tclean_base):


    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_cube(self):
        '''
        Mosaic cube imaging
        field SMIDGE_NWCloud, spw 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'mosaic_cube.iter'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[108, 108], cell=['1.1arcsec'], \
            phasecenter='ICRS 00:45:54.3836 -073.15.29.413', stokes='I', \
            specmode='cube', nchan=508, start='220.2526744GHz', \
            width='0.244174087287MHz', outframe='LSRK', \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        if self.parallel:
            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'],scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526744GHz', width='0.244174087287MHz', \
                outframe='LSRK', perchanweightdensity=False, \
                gridder='mosaic', chanchunks=-1, mosweight=True, \
                usepointing=False, pblimit=0.2, deconvolver='hogbom', \
                restoration=True, pbcor=True, weighting='briggs', robust=0.5,\
                npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
                interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                savemodel='none', calcres=False, calcpsf=False, \
                parallel=True)

            # retrieve per-channel beam statistics
            img = os.getcwd()+'/'+file_name+'1'
            res_bmin_dict, res_bmaj_dict, res_pa_dict = \
                self.cube_beam_stats(img+'.image')

            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526744GHz', width='0.244174087287MHz', \
                outframe='LSRK', perchanweightdensity=False, \
                gridder='mosaic', chanchunks=-1, mosweight=True, \
                usepointing=False, pblimit=0.2, deconvolver='hogbom', \
                restoration=True, restoringbeam='common', pbcor=True, \
                weighting='briggs', robust=0.5, npixels=0, niter=0, \
                threshold='0.354Jy', nsigma=0.0, interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                savemodel='none', calcres=False, calcpsf=False, 
                parallel=False)

        else:
            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'],scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526744GHz', width='0.244174087287MHz', \
                outframe='LSRK', perchanweightdensity=False, \
                gridder='mosaic', chanchunks=-1, mosweight=True, \
                usepointing=False, pblimit=0.2, deconvolver='hogbom', \
                restoration=True, pbcor=True, weighting='briggs', robust=0.5,\
                npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
                interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                savemodel='none', calcres=False, calcpsf=False, \
                restoringbeam='common', parallel=False)

            img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/mosaic_cube.image.crtf')

        exp_im_stats = {'com_bmaj': 8.79758391563,
            'com_bmin': 6.10500958355,
            'com_pa': 64.9303736341,
            'npts': 5925312,
            'npts_unmasked': 3338068.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 1.18474924564,
            'max_val_pos':[44, 38, 0, 252],
            'min_val': -0.417252391577,
            'min_val_pos':[24, 31, 0, 498],
            'im_rms': 0.0878815938952,
            'rms_per_chan': 1,
            'im_sum': 365.155811516,
            'regn_sum': 91.3425360965,
            'npts_real': 5925312,
            'rms_per_field': 1}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = casaTestHelper.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 5925312,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'mask_pix': 3929,
            'mask_regns': 1,
            'npts_real': 5925312}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/mosaic_cube.pb.crtf')

        exp_pb_stats = {'npts': 5925312,
            'npts_unmasked': 3338068.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[54, 54, 0, 0],
            'min_val': 0.200047940016,
            'min_val_pos':[98, 64, 0, 339],
            'npts_0.2': 3338068,
            'npts_0.5': 1825336,
            'npts_real': 5925312}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/mosaic_cube.psf.crtf')

        exp_psf_stats = {'npts': 5925312,
            'npts_unmasked': 5925312.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[54, 54, 0, 0],
            'min_val': -0.168375626206,
            'min_val_pos':[63, 54, 0, 12],
            'im_rms':  0.0604957837605,
            'im_sum': 67.0230730532,
            'npts_real': 5925312}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/mosaic_cube.residual.crtf')

        exp_resid_stats = {'npts': 5925312,
            'npts_unmasked': 3338068.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 0.490073978901,
            'max_val_pos':[49, 60, 0, 249],
            'min_val': -0.417252391577,
            'min_val_pos':[24, 31, 0, 498],
            'im_rms': 0.087537475151,
            'im_sum': 68.1217086753,
            'npts_real': 5925312}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/mosaic_cube.model.crtf')

        exp_model_stats = {'npts': 5925312,
            'npts_unmasked': 5925312.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 0.510256052017,
            'max_val_pos':[55, 57, 0, 255],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.000534824373949,
            'im_sum': 5.9057832174,
            'regn_sum': 1.13614080101,
            'mask_non0': 32,
            'npts_real': 5925312}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 508,
            'npts_unmasked': 508.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 120.900192261,
            'max_val_pos':[0, 0, 0, 447],
            'min_val': 120.665283203,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 508}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # .weight report
        wt_stats_dict = self.image_stats(img, '.weight')

        exp_wt_stats = {'npts': 5925312,
            'npts_unmasked': 5925312.0,
            'freq_bin': 244174.08728027344,
            'start': 0,
            'end': 0,
            'nchan': 508,
            'max_val': 0.393758654594,
            'max_val_pos':[54, 54, 0, 0],
            'min_val': 7.45326979086e-05,
            'min_val_pos':[97, 107, 0, 170],
            'im_rms': 0.140904168376,
            'im_sum': 506774.321109,
            'npts_0.2': 1058307,
            'npts_0.3': 504152,
            'npts_real': 5925312}

        out, report8_a = casaTestHelper.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = casaTestHelper.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = casaTestHelper.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = casaTestHelper.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = casaTestHelper.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = casaTestHelper.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = casaTestHelper.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = casaTestHelper.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = casaTestHelper.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = casaTestHelper.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = casaTestHelper.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = casaTestHelper.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[0.15, 1.2])
        self.png_creator(img+'.residual', range_list=[0.15, 1.2])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mfs(self):
        '''
        Mosaic MFS imaging
        field NGC5363, spw 16 & 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'mosaic_mfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='NGC5363', spw=['0:113.893653412'
            '~114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[126, 126], cell=['2arcsec'], \
            phasecenter='ICRS 13:56:07.2100 +005.15.17.200', stokes='I', \
            specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom',  restoration=False, restoringbeam='common',\
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        file_list = glob.glob(file_name+'0*')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='NGC5363', spw=['0:113.893653412'
            '~114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[126, 126], cell=['2arcsec'], \
            phasecenter='ICRS 13:56:07.2100 +005.15.17.200', stokes='I', \
            specmode='mfs', nchan=-1, outframe='LSRK', \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=30000, threshold='0.00723Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=1.25, \
            noisethreshold=5.0, lownoisethreshold=2.0, \
            negativethreshold=0.0, minbeamfrac=0.1, growiterations=75, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/mosaic_mfs.image.crtf')

        exp_im_stats = {'com_bmaj': 17.6737785339,
            'com_bmin': 10.060172081,
            'com_pa': 86.6785964966,
            'npts': 15876,
            'npts_unmasked': 8454.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0345157124102,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00193646573462,
            'min_val_pos':[91, 52, 0, 0],
            'im_rms':  0.00202865568675,
            'rms_per_chan': 1,
            'im_sum': 1.51296002558,
            'regn_sum': 1.58850855171,
            'npts_real': 15876,
            'rms_per_field': 1}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [64, 114, 0, 0]), \
                      (img+'.image', False, [64, 115, 0, 0]), \
                      (img+'.image', True, [11, 60, 0, 0]), \
                      (img+'.image', False, [10, 60, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = casaTestHelper.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 15876,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 360,
            'mask_regns': 1,
            'npts_real': 15876}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/mosaic_mfs.pb.crtf')

        exp_pb_stats = {'npts': 15876,
            'npts_unmasked': 8454.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 0.200080409646,
            'min_val_pos':[102, 28, 0, 0],
            'npts_0.2': 8454,
            'npts_0.5': 4497,
            'npts_real': 15876}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/mosaic_mfs.psf.crtf')

        exp_psf_stats = {'npts': 15876,
            'npts_unmasked': 15876.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.169576942921,
            'min_val_pos':[61, 57, 0, 0],
            'im_rms':  0.0501544137568,
            'im_sum': 0.00266315032371,
            'npts_real': 15876}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/mosaic_mfs.residual.crtf')

        exp_resid_stats = {'npts': 15876,
            'npts_unmasked': 8454.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.00713972514495,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00197902694345,
            'min_val_pos':[72, 58, 0, 0],
            'im_rms': 0.000868739852256,
            'im_sum': 0.134135425201,
            'npts_real': 15876}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/mosaic_mfs.model.crtf')

        exp_model_stats = {'npts': 15876,
            'npts_unmasked': 15876.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0273759867996,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.000217269736505,
            'im_sum': 0.0273759867996,
            'regn_sum': 0.0273759867996,
            'mask_non0': 1,
            'npts_real': 15876}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 4396210.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 4396210.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # .weight report
        wt_stats_dict = self.image_stats(img, '.weight')

        exp_wt_stats = {'npts': 15876,
            'npts_unmasked': 15876.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.426215529442,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 4.20721153205e-05,
            'min_val_pos':[1, 24, 0, 0],
            'im_rms': 0.144357241177,
            'im_sum': 1347.62643852,
            'npts_0.2': 2778,
            'npts_0.3': 1468,
            'npts_real': 15876}

        out, report8_a = casaTestHelper.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = casaTestHelper.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = casaTestHelper.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = casaTestHelper.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = casaTestHelper.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = casaTestHelper.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = casaTestHelper.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = casaTestHelper.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = casaTestHelper.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = casaTestHelper.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = casaTestHelper.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = casaTestHelper.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[-0.002, 0.035])
        self.png_creator(img+'.residual', range_list=[-0.002, 0.035])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mtmfs(self):
        '''
        Mosaic mtmfs imaging
        field NGC5363, spw 16 & 22
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'mosaic_mtmfs.iter'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='NGC5363', spw=['0:113.893653412'
            '~114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', nterms=2, imsize=[126, 126], \
            cell=['2arcsec'], phasecenter='ICRS 13:56:07.2100 '
            '+005.15.17.200', stokes='I', specmode='mfs', nchan=-1, \
            outframe='LSRK', perchanweightdensity=False, gridder='mosaic', \
            chanchunks=-1, mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='mtmfs', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='NGC5363', spw=['0:113.893653412'
            '~114.629981537GHz;114.8809581~115.758887787GHz,1:99.90983902'
            '~100.494799957GHz;100.68327652~101.77312027GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['6,9'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', nterms=2, imsize=[126, 126], \
            cell=['2arcsec'], phasecenter='ICRS 13:56:07.2100 '
            '+005.15.17.200', stokes='I', specmode='mfs', nchan=-1, \
            outframe='LSRK', perchanweightdensity=False, gridder='mosaic', \
            chanchunks=-1, mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='mtmfs', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=30000, threshold='0.00723Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=1.25, \
            noisethreshold=5.0, lownoisethreshold=2.0, negativethreshold=0.0,\
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mos_mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.image.tt0.crtf')

        exp_im_stats = {'com_bmaj': 17.6737785339,
            'com_bmin': 10.060172081,
            'com_pa': 86.6785964966,
            'npts': 15876,
            'npts_unmasked': 8454.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0349145904183,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00258040986955,
            'min_val_pos':[70, 97, 0, 0],
            'im_rms':  0.00207815366878,
            'rms_per_chan': 1,
            'im_sum': 1.47797238623,
            'regn_sum': 1.58228412218,
            'npts_real': 15876,
            'rms_per_field': 1}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [64, 114, 0, 0]), \
                      (img+'.image.tt0', False, [64, 115, 0, 0]), \
                      (img+'.image.tt0', True, [11, 60, 0, 0]), \
                      (img+'.image.tt0', False, [10, 60, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = casaTestHelper.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 15876,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 360,
            'mask_regns': 1,
            'npts_real': 15876}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.pb.tt0.crtf')

        exp_pb_stats = {'npts': 15876,
            'npts_unmasked': 8454.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 0.200080409646,
            'min_val_pos':[102, 28, 0, 0],
            'npts_0.2': 8454,
            'npts_0.5': 4497,
            'npts_real': 15876}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.psf.tt0.crtf')

        exp_psf_stats = {'npts': 15876,
            'npts_unmasked': 15876.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.169576942921,
            'min_val_pos':[61, 57, 0, 0],
            'im_rms':  0.0501544137568,
            'im_sum': 0.00266315032371,
            'npts_real': 15876}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual.tt0', \
            region_file = data_path+
            'region_files/mosaic_mtmfs.residual.tt0.crtf')

        exp_resid_stats = {'npts': 15876,
            'npts_unmasked': 8454.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.00713019352406,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00197783415206,
            'min_val_pos':[72, 58, 0, 0],
            'im_rms': 0.000865575412319,
            'im_sum': 0.130538275658,
            'npts_real': 15876}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.model.tt0.crtf')

        exp_model_stats = {'npts': 15876,
            'npts_unmasked': 15876.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0272394195199,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.000216185869206,
            'im_sum': 0.0272394195199,
            'regn_sum': 0.0272394195199,
            'mask_non0': 1,
            'npts_real': 15876}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt.tt0')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 4396210.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 4396210.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # .weight report
        wt_stats_dict = self.image_stats(img, '.weight.tt0')

        exp_wt_stats = {'npts': 15876,
            'npts_unmasked': 15876.0,
            'freq_bin': 15849925874.83342,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.426215529442,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 4.20721153205e-05,
            'min_val_pos':[1, 24, 0, 0],
            'im_rms': 0.144357241177,
            'im_sum': 1347.62643852,
            'npts_0.2': 2778,
            'npts_0.3': 1468,
            'npts_real': 15876}

        out, report8_a = casaTestHelper.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = casaTestHelper.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = casaTestHelper.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = casaTestHelper.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = casaTestHelper.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = casaTestHelper.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = casaTestHelper.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = casaTestHelper.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = casaTestHelper.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = casaTestHelper.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = casaTestHelper.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = casaTestHelper.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image.tt0', range_list=[-0.003, 0.035])
        self.png_creator(img+'.residual.tt0', range_list=[-0.003, 0.035])

        test_dict[testname]['images'] = \
            [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)

#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_cube_eph(self):
        '''
        Mosaic ephemeris cube imaging
        field Venus, spw 45
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'mosaic_cube_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='Venus', spw=['1'], antenna=['0,1,2,3,'
            '4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,'
            '27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46'], \
            scan=['7,11'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='corrected', imagename=file_name+'0', \
            imsize=[480, 420], cell=['0.14arcsec'], phasecenter='TRACKFIELD',\
            stokes='I', specmode='cubesource', nchan=948, start=5, width=1, \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=15.0, minbeamfrac=0.3, \
            growiterations=50, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='Venus', spw=['1'], antenna=['0,1,2,3,'
            '4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,'
            '27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46'], \
            scan=['7,11'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='corrected', imagename=file_name+'1', \
            imsize=[480, 420], cell=['0.14arcsec'], phasecenter='TRACKFIELD',\
            stokes='I', specmode='cubesource', nchan=948, start=5, width=1, \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=700000, threshold='0.0106Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=2.0, \
            noisethreshold=4.25, lownoisethreshold=1.5, \
            negativethreshold=15.0, minbeamfrac=0.3, growiterations=50, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/mosaic_cube_eph.image.crtf')

        exp_im_stats = {'com_bmaj': 0.930861414189,
            'com_bmin': 0.71932871452,
            'com_pa': -88.1227929847,
            'npts': 191116800,
            'npts_unmasked': 105006865.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 0.116102233529,
            'max_val_pos':[286, 233, 0, 568],
            'min_val': -0.0639033019543,
            'min_val_pos':[211, 186, 0, 592],
            'im_rms':  0.0106607721423,
            'rms_per_chan': 1,
            'im_sum': 7294.44604131,
            'regn_sum': 28.4447036739,
            'npts_real': 191116800,
            'rms_per_field': 1}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [209, 387, 0, 0]), \
                      (img+'.image', False, [209, 388, 0, 0]), \
                      (img+'.image', True, [20, 204, 0, 0]), \
                      (img+'.image', False, [19, 204, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = casaTestHelper.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 191116800,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'mask_pix': 156228,
            'mask_regns': 39,
            'npts_real': 191116800}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/mosaic_cube_eph.pb.crtf')

        exp_pb_stats = {'npts': 191116800,
            'npts_unmasked': 105006865.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 1.0,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 0.200000017881,
            'min_val_pos':[56, 302, 0, 567],
            'npts_0.2': 105006529,
            'npts_0.5': 60643408,
            'npts_real': 191116800}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/mosaic_cube_eph.psf.crtf')

        exp_psf_stats = {'npts': 191116800,
            'npts_unmasked': 191116800.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 1.0,
            'max_val_pos':[240, 210, 0, 0],
            'min_val': -0.0456548333168,
            'min_val_pos':[230, 216, 0, 14],
            'im_rms':  0.0122779442959,
            'im_sum': 130.930138815,
            'npts_real': 191116800}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/mosaic_cube_eph.residual.crtf')

        exp_resid_stats = {'npts': 191116800,
            'npts_unmasked': 105006865.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 0.0589463338256,
            'max_val_pos':[269, 264, 0, 765],
            'min_val': -0.0639033019543,
            'min_val_pos':[211, 186, 0, 592],
            'im_rms': 0.0105758073519,
            'im_sum': 3583.31516361,
            'npts_real': 191116800}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/mosaic_cube_eph.model.crtf')

        exp_model_stats = {'npts': 191116800,
            'npts_unmasked': 191116800.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 0.0502305738628,
            'max_val_pos':[258, 261, 0, 504],
            'min_val': -0.0168753024191,
            'min_val_pos':[304, 224, 0, 652],
            'im_rms': 7.7618298009e-05,
            'im_sum': 95.9311681981,
            'regn_sum': 0.698479007697,
            'mask_non0': 15139,
            'npts_real': 191116800}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 948,
            'npts_unmasked': 948.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 45510.7695312,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 45195.4453125,
            'min_val_pos':[0, 0, 0, 594],
            'npts_real': 948}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # .weight report
        wt_stats_dict = self.image_stats(img, '.weight')

        exp_wt_stats = {'npts': 191116800,
            'npts_unmasked': 191116800.0,
            'freq_bin': 244151.1796875,
            'start': 0,
            'end': 0,
            'nchan': 948,
            'max_val': 0.317849487066,
            'max_val_pos':[211, 203, 0, 10],
            'min_val': 6.55390249449e-05,
            'min_val_pos':[451, 65, 0, 947],
            'im_rms': 0.119764405387,
            'im_sum': 13810530.8091,
            'npts_0.2': 28707625,
            'npts_0.3': 6190533,
            'npts_real': 191116800}

        out, report8_a = casaTestHelper.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = casaTestHelper.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = casaTestHelper.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = casaTestHelper.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = casaTestHelper.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = casaTestHelper.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = casaTestHelper.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = casaTestHelper.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = casaTestHelper.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = casaTestHelper.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = casaTestHelper.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = casaTestHelper.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[-0.01, 0.1])
        self.png_creator(img+'.residual', range_list=[-0.01, 0.1])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mfs_eph(self):
        '''
        Mosaic ephemeris mfs imaging
        field Venus, spw 25 & 45
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'mosaic_mfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='Venus', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,'
            '42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[480, 420], cell=['0.14arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='mfs', nchan=-1, \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto' \
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='Venus', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,'
            '42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[480, 420], cell=['0.14arcsec'], \
            phasecenter='TRACKFIELD', specmode='mfs', nchan=-1, \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=7000000, threshold='0.0316Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=2.0, \
            noisethreshold=4.25, lownoisethreshold=1.5, \
            negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, calcres=False, calcpsf=False, \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image', region_file = \
            data_path+'region_files/mosaic_mfs_eph.image.crtf')

        exp_im_stats = {'com_bmaj': 0.914226949215,
            'com_bmin': 0.708592534065,
            'com_pa': -89.3612976074,
            'npts': 201600,
            'npts_unmasked': 113589.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 2.06048989296,
            'max_val_pos':[291, 212, 0, 0],
            'min_val': -2.1858522892,
            'min_val_pos':[290, 152, 0, 0],
            'im_rms':  0.676557465791,
            'rms_per_chan': 1,
            'im_sum': 5498.32523989,
            'regn_sum': 8725.50744967,
            'npts_real': 201600,
            'rms_per_field': 1}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [211, 390, 0, 0]), \
                      (img+'.image', False, [211, 391, 0, 0]), \
                      (img+'.image', True, [18, 205, 0, 0]), \
                      (img+'.image', False, [17, 205, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = casaTestHelper.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 201600,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 201600}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb', region_file = \
            data_path+'region_files/mosaic_mfs_eph.pb.crtf')

        exp_pb_stats = {'npts': 201600,
            'npts_unmasked': 113589.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 0.200001135468,
            'min_val_pos':[81, 343, 0, 0],
            'npts_0.2': 113589,
            'npts_0.5': 64549,
            'npts_real': 201600}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf', region_file = \
            data_path+'region_files/mosaic_mfs_eph.psf.crtf')

        exp_psf_stats = {'npts': 201600,
            'npts_unmasked': 201600.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[240, 210, 0, 0],
            'min_val': -0.0525086708367,
            'min_val_pos':[230, 216, 0, 0],
            'im_rms':  0.0111846421981,
            'im_sum': 0.100949260701,
            'npts_real': 201600}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual', region_file = \
            data_path+'region_files/mosaic_mfs_eph.residual.crtf')

        exp_resid_stats = {'npts': 201600,
            'npts_unmasked': 113589.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 2.06048989296,
            'max_val_pos':[291, 212, 0, 0],
            'min_val': -2.1858522892,
            'min_val_pos':[290, 152, 0, 0],
            'im_rms': 0.676557465791,
            'im_sum': 5498.32523989,
            'npts_real': 201600}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model', region_file = \
            data_path+'region_files/mosaic_mfs_eph.model.crtf')

        exp_model_stats = {'npts': 201600,
            'npts_unmasked': 201600.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.0,
            'im_sum': 0.0,
            'regn_sum': 0.0,
            'mask_non0': 0,
            'npts_real': 201600}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 30068706.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 30068706.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # .weight report
        wt_stats_dict = self.image_stats(img, '.weight')

        exp_wt_stats = {'npts': 201600,
            'npts_unmasked': 201600.0,
            'freq_bin': 16762504556.453735,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.333539396524,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 9.39047822612e-05,
            'min_val_pos':[451, 64, 0, 0],
            'im_rms': 0.12506836881,
            'im_sum': 15366.9703442,
            'npts_0.2': 32025,
            'npts_0.3': 9855,
            'npts_real': 201600}

        out, report8_a = casaTestHelper.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = casaTestHelper.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = casaTestHelper.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = casaTestHelper.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = casaTestHelper.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = casaTestHelper.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = casaTestHelper.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = casaTestHelper.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = casaTestHelper.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = casaTestHelper.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = casaTestHelper.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = casaTestHelper.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        casaTestHelper.add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image', range_list=[-2.2, 2.1])
        self.png_creator(img+'.residual', range_list=[-2.2, 2.1])

        test_dict[testname]['images'] = \
            [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @casaTestHelper.stats_dict(test_dict)
    def test_mosaic_mtmfs_eph(self):
        '''
        Mosaic ephemeris mtmfs imaging
        field Venus, spw 25 & 45
        '''

        testname, testdoc = self.getNameDoc()
        file_name = 'mosaic_mtmfs_eph.iter'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='Venus', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[480, 420], nterms=2, \
            cell=['0.14arcsec'], phasecenter='TRACKFIELD', stokes='I', \
            specmode='mfs', nchan=-1, perchanweightdensity=False, \
            gridder='mosaic', chanchunks=-1, mosweight=True, \
            usepointing=False, pblimit=0.2, deconvolver='mtmfs', \
            restoration=False, restoringbeam='common', pbcor=False, \
            weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto' \
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_images(file_name)

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='Venus', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'1', imsize=[480, 420], nterms=2, \
            cell=['0.14arcsec'], phasecenter='TRACKFIELD', specmode='mfs', \
            nchan=-1, perchanweightdensity=False, gridder='mosaic', \
            chanchunks=-1, mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='mtmfs', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggs', robust=0.5, npixels=0, \
            niter=7000000, threshold='0.0316Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=2.0, \
            noisethreshold=4.25, lownoisethreshold=1.5, \
            negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, calcres=False, calcpsf=False,  \
            parallel=self.parallel)

        img = os.getcwd()+'/'+file_name+'1'

        report0 = casaTestHelper.checkall( \
            imgexist = self.image_list(img, 'mos_mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img, '.image.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.image.tt0.crtf')

        exp_im_stats = {'com_bmaj': 0.914226949215,
            'com_bmin': 0.708592534065,
            'com_pa': -89.3612976074,
            'npts': 201600,
            'npts_unmasked': 113589.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 2.05563259125,
            'max_val_pos':[291, 212, 0, 0],
            'min_val': -2.19232487679,
            'min_val_pos':[290, 152, 0, 0],
            'im_rms':  0.672551533853,
            'rms_per_chan': 1,
            'im_sum': 5410.97024339,
            'regn_sum': 8700.32470433,
            'npts_real': 201600,
            'rms_per_field': 1}

        report1_a = casaTestHelper.checkall( \
#            # checks the peak flux value and location
#            imgval = [(img+'.image', 1.01393294334, [52, 50, 0, 6])],
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [211, 390, 0, 0]), \
                      (img+'.image.tt0', False, [211, 391, 0, 0]), \
                      (img+'.image.tt0', True, [18, 205, 0, 0]), \
                      (img+'.image.tt0', False, [17, 205, 0, 0])])

        out, report1_b = casaTestHelper.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = casaTestHelper.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = casaTestHelper.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = casaTestHelper.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = casaTestHelper.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = casaTestHelper.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = casaTestHelper.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = casaTestHelper.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = casaTestHelper.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = casaTestHelper.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = casaTestHelper.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = casaTestHelper.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = casaTestHelper.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = casaTestHelper.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = casaTestHelper.check_val( \
            im_stats_dict['rms_per_chan'], exp_im_stats['rms_per_chan'], \
            valname='RMS per channel of .image', exact=False, \
            epsilon=0.01)
        out, report1_u = casaTestHelper.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = casaTestHelper.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = casaTestHelper.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 201600,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 201600}

        out, report2_a = casaTestHelper.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = casaTestHelper.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = casaTestHelper.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = casaTestHelper.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = casaTestHelper.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = casaTestHelper.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = casaTestHelper.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report2 = report2_a + report2_b + report2_c + report2_d + \
            report2_e + report2_f + report2_g + report2_h

        # .pb report
        pb_stats_dict = self.image_stats(img, '.pb.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.pb.tt0.crtf')

        exp_pb_stats = {'npts': 201600,
            'npts_unmasked': 113589.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 0.200001135468,
            'min_val_pos':[81, 343, 0, 0],
            'npts_0.2': 113589,
            'npts_0.5': 64549,
            'npts_real': 201600}

        out, report3_a = casaTestHelper.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = casaTestHelper.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = casaTestHelper.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = casaTestHelper.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = casaTestHelper.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = casaTestHelper.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = casaTestHelper.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = casaTestHelper.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = casaTestHelper.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = casaTestHelper.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = casaTestHelper.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = casaTestHelper.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.psf.tt0.crtf')

        exp_psf_stats = {'npts': 201600,
            'npts_unmasked': 201600.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[240, 210, 0, 0],
            'min_val': -0.0525086708367,
            'min_val_pos':[230, 216, 0, 0],
            'im_rms':  0.0111846421981,
            'im_sum': 0.100949260701,
            'npts_real': 201600}

        out, report4_a = casaTestHelper.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = casaTestHelper.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = casaTestHelper.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = casaTestHelper.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = casaTestHelper.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = casaTestHelper.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = casaTestHelper.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = casaTestHelper.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = casaTestHelper.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = casaTestHelper.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = casaTestHelper.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = casaTestHelper.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual.tt0', \
            region_file = data_path+
            'region_files/mosaic_mtmfs_eph.residual.tt0.crtf')

        exp_resid_stats = {'npts': 201600,
            'npts_unmasked': 113589.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 2.06048989296,
            'max_val_pos':[291, 212, 0, 0],
            'min_val': -2.1858522892,
            'min_val_pos':[290, 152, 0, 0],
            'im_rms': 0.676557465791,
            'im_sum': 5498.32523989,
            'npts_real': 201600}

        out, report5_a = casaTestHelper.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = casaTestHelper.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = casaTestHelper.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = casaTestHelper.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = casaTestHelper.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = casaTestHelper.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = casaTestHelper.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = casaTestHelper.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = casaTestHelper.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = casaTestHelper.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = casaTestHelper.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = casaTestHelper.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report5 = report5_a + report5_b + report5_c + report5_d + \
            report5_e + report5_f + report5_g + report5_h + report5_i + \
            report5_j + report5_k + report5_l + report5_m + report5_n + \
            report5_o + report5_p + report5_q

        # .model report
        model_stats_dict = self.image_stats(img, '.model.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.model.tt0.crtf')

        exp_model_stats = {'npts': 201600,
            'npts_unmasked': 201600.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 0.0,
            'min_val_pos':[0, 0, 0, 0],
            'im_rms': 0.0,
            'im_sum': 0.0,
            'regn_sum': 0.0,
            'mask_non0': 0,
            'npts_real': 201600}

        out, report6_a = casaTestHelper.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = casaTestHelper.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = casaTestHelper.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = casaTestHelper.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = casaTestHelper.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = casaTestHelper.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = casaTestHelper.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = casaTestHelper.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = casaTestHelper.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = casaTestHelper.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = casaTestHelper.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = casaTestHelper.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = casaTestHelper.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = casaTestHelper.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report6 = report6_a + report6_b + report6_c + report6_d + \
            report6_e + report6_f + report6_g + report6_h + report6_i + \
            report6_j + report6_k + report6_l + report6_m + report6_n + \
            report6_o + report6_p + report6_q + report6_r + report6_s

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img, '.sumwt.tt0')

        exp_sumwt_stats = {'npts': 1,
            'npts_unmasked': 1.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 30068706.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 30068706.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = casaTestHelper.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = casaTestHelper.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = casaTestHelper.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = casaTestHelper.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = casaTestHelper.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = casaTestHelper.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = casaTestHelper.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = casaTestHelper.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # .weight report
        wt_stats_dict = self.image_stats(img, '.weight.tt0')

        exp_wt_stats = {'npts': 201600,
            'npts_unmasked': 201600.0,
            'freq_bin': 16762504556.453705,
            'start': 0,
            'end': 0,
            'nchan': 1,
            'max_val': 0.333539396524,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 9.39047822612e-05,
            'min_val_pos':[451, 64, 0, 0],
            'im_rms': 0.12506836881,
            'im_sum': 15366.9703442,
            'npts_0.2': 32025,
            'npts_0.3': 9855,
            'npts_real': 201600}

        out, report8_a = casaTestHelper.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = casaTestHelper.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = casaTestHelper.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = casaTestHelper.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = casaTestHelper.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = casaTestHelper.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = casaTestHelper.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = casaTestHelper.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = casaTestHelper.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = casaTestHelper.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = casaTestHelper.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = casaTestHelper.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = casaTestHelper.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = casaTestHelper.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = casaTestHelper.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        casaTestHelper.add_to_dict(self, output=test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        self.png_creator(img+'.image.tt0', range_list=[-2.2, 2.1])
        self.png_creator(img+'.residual.tt0', range_list=[-2.2, 2.1])

        test_dict[testname]['images'] = \
            [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(casaTestHelper.check_final(pstr = report), \
            msg = report)


def suite():
     return [Test_standard, Test_mosaic]

# Main #
if __name__ == '__main__':
    unittest.main()



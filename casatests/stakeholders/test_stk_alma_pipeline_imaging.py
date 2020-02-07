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
import numpy
import shutil
import inspect
import scipy

from casatestutils.imagerhelpers import TestHelpers
th = TestHelpers()
from casatestutils import generate_weblog
from casatestutils import add_to_dict
from casatestutils import stats_dict
#from distutils.dir_util import copy_tree


CASA6 = False
try:
    from casatools import ctsys, quanta, measures, image, vpmanager, calibrater
    from casatasks import casalog, delmod, imsubimage, tclean, uvsub, imhead, imsmooth, immath, widebandpbcor, immoments#, imview
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    CASA6 = True
    _ia = image( )
#    _vp = vpmanager( )
#    _cb = calibrater( )
#    _qa = quanta( )
#    _me = measures( )

except ImportError:
    from __main__ import default  # reset given task to its default values
    from tasks import *  # Imports all casa tasks
    from taskinit import *  # Imports all casa tools
    from parallel.parallel_task_helper import ParallelTaskHelper

    _ia = iatool( )
#    _vp = vptool( )
#    _cb = cbtool( )
#    # not local tools
#    _qa = qa
#    _me = me

# location of data
data_path = '/lustre/naasc/sciops/comm/sbooth/CASA_ALMA_pipeline/data_dir/'
#data_path = os.environ.get('CASAPATH').split()[0] + '/casa-data-vt/vlass/'


## Base Test class with Utility functions
class test_tclean_base(unittest.TestCase):

    def setUp(self):
        self._myia = _ia
#        self.epsilon = 0.01 # sets epsilon as a percentage (1%)
        self.msfile = ""
#        self.img = "tst"
        self.img_subdir = 'testdir'
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

    def tearDown(self):
        generate_weblog("tclean_ALMA_pipeline",test_dict)
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
        if (os.path.exists(self.msfile)):
             os.popen('rm -rf ' + self.msfile)
        os.popen('rm -rf ' + self.img_subdir)
        os.popen('rm -rf ' + self.img+'*')

    def prepInputmask(self, maskname=""):
        if maskname!="":
            self.maskname=maskname
        if (os.path.exists(self.maskname)):
            os.popen('rm -rf ' + self.maskname)
        shutil.copytree(refdatapath+self.maskname, self.maskname, symlinks=True)

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

    # compares 2 lists and returns if they are equivalent (within error)
    def check_list_vals(self, list1, list2, epsilon=None):
        if len(list1) == len(list2) and epsilon is None:
            if list1 == list2:
                result = True
            else:
                result = False
        elif len(list1) == len(list2) and epsilon is not None:
            i = 0
            while i < len(list1):
                result, pstr = th.check_val(list1[i], list2[i], \
                    epsilon=epsilon)
                if result == False:
                    break
                i += 1
        else:
            result = False
        return result

    def getNameDoc(self):
#        print inspect.stack()[1]
        testname=inspect.stack()[1][3]
        print("Test name  : " + testname)
#        tname = inspect.getframeinfo(inspect.currentframe()).function
        doc = eval('self.'+testname + '.__doc__')
        print("Doc : " +  doc)
        return testname, doc

    # function to copy iter0 images to iter1 images
#    def copy_images(self, file_name):
#        file_list = glob.glob(file_name+'0*')
#        os.popen('rm -rf '+file_name+'1*')
#        for item in file_list:
#            sp_item = item.split('iter0')
#            copy_tree(item, sp_item[0]+'iter1'+sp_item[1])
#        if os.path.isdir(file_name+'1.workdirectory/'):
#            work0_list = glob.glob(file_name+'0.workdirectory/*')
#            for item in work0_list:
#                sp_item = item.split('iter0')
#                new_item = sp_item[0]+'iter1'+sp_item[1]+'iter1'+sp_item[2]
#                copy_tree(item, new_item)

    # function to copy iter0 images to iter1 images (taken from pipeline)
    def copy_products(self, old_pname, new_pname, ignore=None):
        imlist = glob.glob('%s.*' % (old_pname))
        imlist = [xx for xx in imlist if ignore is None or ignore not in xx]
        for image_name in imlist:
            newname = image_name.replace(old_pname, new_pname)
            if image_name == old_pname + '.workdirectory':
                mkcmd = 'mkdir '+ newname
                os.system(mkcmd)
                self.copy_products(os.path.join(image_name, old_pname), os.path.join(newname, new_pname))
            else:
                shutil.copytree(image_name, newname, symlinks=True)


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
    def image_stats(self, img, suf, region_file=None):#, field_regions=None):
        self._myia.open(img+suf)

        im_size = self._myia.boundingbox()['imageShape'].tolist()
        npts = im_size[0]*im_size[1]*im_size[3]
        start = float( \
            self._myia.statistics()['blcf'].split(', ')[3].split('Hz')[0])
        end = float( \
            self._myia.statistics()['trcf'].split(', ')[3].split('Hz')[0])

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
                self._myia.statistics()['rms'][0], 'im_sum': \
                self._myia.statistics()['sum'][0], 'regn_sum': \
                self._myia.statistics(region=region_file)['sum'][0]}
            try:
                fit_dict = self._myia.fitcomponents( \
                    region=region_file)['results']['component0']
                stats_dict['im_fit'] = [fit_dict['pixelcoords'].tolist(),\
                    fit_dict['spectrum']['channel'], \
                    fit_dict['spectrum']['frequency']['m0']['value'], \
                    fit_dict['peak']['value'], \
                    fit_dict['shape']['majoraxis']['value'], \
                    fit_dict['shape']['minoraxis']['value']]
            except KeyError:
                 stats_dict['im_fit'] = [[1.0, 1.0], 1, 1.0, 1.0, 1.0, 1.0]
            if 'cube' in  img:
                stats_dict['rms_per_chan'] = \
                    self._myia.statistics(axes=[0,1])['rms'].tolist()
            if 'mosaic' in  img:
                #field_regions = []
                stats_dict['rms_per_field'] = 1

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
                numpy.count_nonzero(self._myia.getchunk()>0.5)}
            try:
                fit_dict = self._myia.fitcomponents( \
                    region=region_file)['results']['component0']
                stats_dict['pb_im_fit'] = \
                    [fit_dict['pixelcoords'].tolist(), \
                     fit_dict['spectrum']['channel'], \
                     fit_dict['spectrum']['frequency']['m0']['value'], \
                     fit_dict['peak']['value'], \
                     fit_dict['shape']['majoraxis']['value'], \
                     fit_dict['shape']['minoraxis']['value']]
            except KeyError:
                stats_dict['pb_im_fit'] = [[1.0, 1.0], 1, 1.0, 1.0, 1.0, 1.0]

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
                self._myia.statistics()['sum'][0]}
            try:
                fit_dict = self._myia.fitcomponents( \
                    region=region_file)['results']['component0']
                stats_dict['cen_im_fit'] = \
                    [fit_dict['pixelcoords'].tolist(), \
                     fit_dict['spectrum']['channel'], \
                     fit_dict['spectrum']['frequency']['m0']['value'], \
                     fit_dict['peak']['value'], \
                     fit_dict['shape']['majoraxis']['value'], \
                     fit_dict['shape']['minoraxis']['value']]
            except KeyError:
                stats_dict['cen_im_fit'] = [[1.0, 1.0], 1, 1.0, 1.0, 1.0, 1.0]

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


    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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
            img = file_name+'1'
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

            img = file_name+'1'

        report0 = th.checkall( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 0.94608676433563232,
            'max_val_pos': [38, 36, 0, 254],
            'min_val': -0.70467984676361084,
            'min_val_pos': [18, 57, 0, 374],
            'im_rms': 0.143986095161,
            'rms_per_chan': [0.12634143789549535, 0.14184103443455054, 0.146980848491423, 0.11560714721775701, 0.16661204089962214, 0.14476185721954618, 0.1316304465304373, 0.11433992978005528, 0.1442336908847237, 0.16543686095641913, 0.11873031602382604, 0.15284339493903618, 0.17195923304927352, 0.15046376352374186, 0.14070242684409134, 0.13744696370230045, 0.11158401597402363, 0.12687194824894843, 0.14295610599459097, 0.16862585818094997, 0.13008023644212954, 0.14186757490541813, 0.1169917541756216, 0.13693402712045005, 0.16194773183534902, 0.13634870573065122, 0.13860445389090703, 0.1309492701136035, 0.14819650092974662, 0.15030436484700252, 0.14931368127692368, 0.11984058396768074, 0.13119726238629406, 0.1278997483823513, 0.15680618961802364, 0.14782343113879803, 0.1452811146145065, 0.14962350388870774, 0.12727392822661138, 0.12403611951801675, 0.13971310504808565, 0.14443747442976043, 0.13947457857066817, 0.14603448593891352, 0.1487357653330162, 0.13728792717834695, 0.12754218448487697, 0.13472363429296094, 0.17318897000268654, 0.15007875971445414, 0.1210452212469422, 0.15977553256440455, 0.13077138200444186, 0.12679151267047647, 0.12091204082027505, 0.1338966333695089, 0.13556991771575277, 0.15456345918134376, 0.12044465611280492, 0.14836982561861023, 0.1349896116866282, 0.15311214064438922, 0.11251497655504887, 0.134867796496014, 0.13574313457554038, 0.14582224580240324, 0.12753531221719416, 0.15335445312643003, 0.13482732612307324, 0.1622050903445585, 0.13260306174268546, 0.1345326608100535, 0.16404765102131583, 0.13449430188802702, 0.14543809289295098, 0.1606584196112734, 0.12484651484486906, 0.16251383851634701, 0.13756025624117688, 0.13165353467440083, 0.1308248320448295, 0.14752778635690292, 0.1274645256107852, 0.16421712463271607, 0.15255317243782812, 0.1497707840063393, 0.11911825364867326, 0.14541033702618353, 0.1659723426787793, 0.1554971226410762, 0.14703675741501698, 0.12325846980328654, 0.15070706791866434, 0.1243073669840061, 0.13646642844468243, 0.1301143392639293, 0.12734602178400867, 0.1553600823593344, 0.15035594210430997, 0.11530605847413075, 0.1611567346343003, 0.12221832384850957, 0.14207389319672978, 0.14522516033398006, 0.1345322788758837, 0.1486176245373929, 0.15765848896613346, 0.1308440759384876, 0.1466820831226493, 0.13598865468593319, 0.15187538855740168, 0.1478468013010444, 0.1383732515889412, 0.1276861625889527, 0.11697230161534232, 0.13739607388184524, 0.11303259344169146, 0.1361001584583741, 0.12857356426667815, 0.1437570752313611, 0.13169397143643052, 0.15326431411050365, 0.12383180315967929, 0.1526310794015497, 0.14746177769245866, 0.15194893390457265, 0.1421630320154613, 0.15662308690272084, 0.12239198421329735, 0.12071542153915982, 0.14268554321174182, 0.13489697242976567, 0.15127855443293006, 0.1542443819690316, 0.15752918577920158, 0.11478434733366248, 0.17298964180575135, 0.13177526480150695, 0.12236732291938952, 0.15625856947990782, 0.13687165189461548, 0.1536631153928859, 0.14669563803395924, 0.1277170908624889, 0.14966567842171496, 0.12823515897560267, 0.13577828413547297, 0.16140169123660877, 0.13133284404676335, 0.14223570583416104, 0.1603292311222728, 0.10759630495294702, 0.15787039978749143, 0.1327200609847152, 0.14655899389809018, 0.14008820956915727, 0.1442107348583108, 0.1317943450568934, 0.12972989243424995, 0.1625036947147829, 0.12241712383574781, 0.14998173521745944, 0.13352731228428555, 0.1741676258276787, 0.15545996482656257, 0.13121844421079562, 0.1389256768353536, 0.1475992903718036, 0.14205849908080379, 0.14975427804440275, 0.1532491403618113, 0.12531915969323904, 0.14153689035122896, 0.16741877503811964, 0.1355536447212321, 0.12548585056941425, 0.16334800417248366, 0.14220841606737944, 0.1376802362928535, 0.1394159389365598, 0.1533008119644231, 0.12568227593323275, 0.14138024496799945, 0.14688836279261966, 0.12037367892758656, 0.12335138886587714, 0.16740640885840646, 0.11756235238942149, 0.13221931449560975, 0.14605469946826174, 0.12287649136200192, 0.13900407591276098, 0.1477935699475207, 0.14723640198504923, 0.12637771862286276, 0.14264989851200444, 0.14188497863070984, 0.1517498748029243, 0.1745550071541481, 0.14693061119966988, 0.12180541963696558, 0.17178472812899895, 0.134842796032342, 0.1587769050427257, 0.16022475326023228, 0.12598385136025822, 0.12173065475536829, 0.1358700032273519, 0.12249230371601251, 0.1320416693266833, 0.1380195667444624, 0.17036819494074398, 0.14449179298441997, 0.1363579047545357, 0.15814587607932587, 0.1387404461979901, 0.13421674959986293, 0.1221729254232071, 0.15007074873391474, 0.1519841002224019, 0.17405910974305452, 0.10810253208919626, 0.14404509620995673, 0.12925102011532486, 0.13284702789448985, 0.16316517507291742, 0.18004246985230368, 0.12352109323053732, 0.13438971701846103, 0.14110722423724795, 0.15240505247738928, 0.16299890474660164, 0.13862726296963718, 0.13653417057201506, 0.15748574227626927, 0.13330507817933285, 0.11630210517195279, 0.14310200319865532, 0.16947649357417122, 0.19276632648628003, 0.1442624172150719, 0.12588748723036136, 0.13766261231222038, 0.15501574319393477, 0.1467664214746981, 0.1437631603776764, 0.13281269178755073, 0.1499498907469051, 0.1547831067839161, 0.1650203851926907, 0.19221241068202002, 0.17354684080024252, 0.1914610395870755, 0.1870175530216679, 0.17581044678251068, 0.1969886211075383, 0.1793297247200161, 0.1650936078172776, 0.14492844003488278, 0.14169417048994912, 0.13913741690930143, 0.15199069574737809, 0.14420138013804454, 0.13427068237336662, 0.14854837369704055, 0.16108337901230768, 0.14392226756737841, 0.17357790177064894, 0.12582205913761169, 0.14917348406521977, 0.14570283814332685, 0.1494709027791571, 0.15333214078415874, 0.13788344445589443, 0.15113305293127033, 0.15910328490835218, 0.1258386682524315, 0.14889037473182778, 0.12206410303883597, 0.14151660747434333, 0.12015317625283857, 0.13619775353885638, 0.1360718460404539, 0.12638979179451687, 0.13793807234055996, 0.11510797437711365, 0.151683286728862, 0.12186947178372404, 0.1334587198282186, 0.1416652191079726, 0.17466019895415022, 0.1378517218379097, 0.11249994464540782, 0.1279851283688596, 0.14018299952556243, 0.14434793746130153, 0.16650070810701043, 0.16748683720689606, 0.11940475115442413, 0.1276077403883665, 0.12359063458201924, 0.1412078454825856, 0.14551852706758697, 0.14986530033150053, 0.1325683624037802, 0.13820026200226213, 0.12812701294141524, 0.12753530671589502, 0.13497760358639382, 0.170481828949963, 0.13917782659135283, 0.14415530157163367, 0.17676595928739716, 0.16642659928559336, 0.1488633024931477, 0.13762287309163865, 0.13173759984975023, 0.13142953909930827, 0.13550297890119667, 0.1378187724239619, 0.13522211241434626, 0.16692597968173054, 0.12681333395848324, 0.14653010389914287, 0.12049158054527685, 0.1394775351353454, 0.1308254468074764, 0.15490863570934446, 0.13846726479055535, 0.12889021958687075, 0.1475440113411787, 0.1317061942800524, 0.13705790461652367, 0.13974036959646274, 0.14395769504564201, 0.16390950740528484, 0.14702614449178092, 0.11790240659365156, 0.15825207098676003, 0.1572387243129431, 0.16049823864649457, 0.15795518128891234, 0.12433544601523717, 0.12148637251440078, 0.17047767727369215, 0.1542909493696503, 0.14083894029456667, 0.1273419049433489, 0.12328954463396018, 0.14383278211432127, 0.16043058197358323, 0.12709149940547382, 0.14015718181807055, 0.13290234816222832, 0.14241483672006272, 0.13977682956445464, 0.19082240791944033, 0.1362551319469732, 0.13878314241029005, 0.13024580725200274, 0.17584969348391408, 0.1585665364841856, 0.1753088252644526, 0.14469549936190076, 0.12561131070797885, 0.1275068394660626, 0.1246509966110453, 0.1353771632037516, 0.14968376152813326, 0.1306367587605417, 0.17509696143235529, 0.12303417070683348, 0.12984385988761357, 0.15655156396975772, 0.15158283415352666, 0.14190571663235674, 0.16013330804491582, 0.15701423529031192, 0.1437377861978616, 0.13178714397605581, 0.14896642830351547, 0.1543110828211962, 0.12729972586137084, 0.13258135641591057, 0.1297547997545033, 0.13623621902113375, 0.16832161541183457, 0.1515993553463645, 0.14230780146154384, 0.15752128344423338, 0.16153801377050137, 0.17727249054646335, 0.17081959421469814, 0.132415345138367, 0.1911207789325767, 0.12167933736058796, 0.16199351646871693, 0.12220003129636688, 0.12276819690997134, 0.1615915263092078, 0.12453980969317847, 0.13583553980995322, 0.15247696319780732, 0.1588409406770753, 0.12523094487161238, 0.14442835547547206, 0.14678401935124402, 0.16646235670431703, 0.13318109022977143, 0.1296318756165084, 0.17474714467039837, 0.14996240468941135, 0.152459961330088, 0.14174866849523865, 0.1531711479592131, 0.13140647178793072, 0.12876054234189516, 0.16073620618951742, 0.1551915897989922, 0.13258798290709237, 0.12079795445447011, 0.17410998739129824, 0.15215648462348455, 0.11414501724048388, 0.1583626102822897, 0.14397582688024274, 0.14206199694234994, 0.1301244502774394, 0.14832273898073517, 0.14965582403580915, 0.11974680479280232, 0.130301576062527, 0.1501180642039695, 0.144040093741809, 0.14644461407160658, 0.16466625212890748, 0.12145879983950111, 0.15515579044778377, 0.11925841619372717, 0.14869156838984582, 0.14241952971976954, 0.14769053467507798, 0.1574041746549612, 0.13729987849339584, 0.13259542009008363, 0.12816033699683535, 0.14196238422903054, 0.14271636018596165, 0.1303221375620509, 0.16261332868157302, 0.16340560251777758, 0.1496448141003338, 0.13870981454297002, 0.12581695569707782, 0.10706545850873779, 0.16477345617357922, 0.13501062540902478, 0.15224041747140973, 0.14304328175572253, 0.1730320818474802, 0.13262178336064356, 0.13626925589130173, 0.19930734468055353, 0.12840465325131772, 0.14565284364647132, 0.14342338225449644, 0.14414011525387105, 0.1314767224882364, 0.14667809040912927, 0.15534831496680646, 0.14770322269613237, 0.18150555408581706, 0.14175110723675505, 0.13126842096913788, 0.1821782599052581, 0.1230212008372949, 0.13879304863372663, 0.1539616899583395, 0.14735461585329224, 0.1103039067487409, 0.1400589610059953, 0.12936178694337058, 0.1185505506087342, 0.15904169231653986, 0.12089109358577753, 0.1479414523691717, 0.13846627238497242, 0.15232750215647736, 0.1318383069114147, 0.13798313113493177, 0.12800661214057107, 0.15270008178936095, 0.12254175320261147, 0.13897505123792747, 0.13607974358604047, 0.15255689754961757, 0.11221218682645963, 0.14990854525350425, 0.15160985722773732, 0.13575646404498873, 0.14576400908648826, 0.15836555616934264, 0.13583817797589506, 0.12589022941650596, 0.15590488932502036, 0.16965183711831752, 0.1483222833838575, 0.14620965767678082],
            'im_sum': 168.242297979,
            'regn_sum': 72.3549563158,
            'npts_real': 3251200,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 70, 0, 0]), \
                (img+'.image', False, [40, 71, 0, 0]), \
                (img+'.image', True, [10, 40, 0, 0]), \
                (img+'.image', False, [9, 40, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( \
            self.check_list_vals(im_stats_dict['rms_per_chan'], \
                exp_im_stats['rms_per_chan'], epsilon=0.01), True, \
            valname='RMS per channel of .image', exact=True)
        out, report1_u = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 3251200,
            'freq_bin': 244174.08728027344,
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'mask_pix': 437,
            'mask_regns': 1,
            'npts_real': 3251200}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.200896695256,
            'min_val_pos':[25, 13, 0, 396],
            'npts_0.2': 1522476,
            'npts_0.5': 736092,
            'npts_real': 3251200,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.218764916062,
            'min_val_pos':[1, 16, 0, 503],
            'im_rms':  0.136036099793,
            'im_sum': 7472.57665916,
            'npts': 3251200,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 0.785612404346,
            'max_val_pos':[42, 43, 0, 256],
            'min_val': -0.704679846764,
            'min_val_pos':[18, 57, 0, 374],
            'im_rms': 0.143918523224,
            'im_sum': 124.317946204,
            'npts': 3251200}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 94.4766769409,
            'max_val_pos':[0, 0, 0, 17],
            'min_val': 94.4766464233,
            'min_val_pos':[0, 0, 0, 449],
            'npts': 508}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[0.3, 1.0])
            self.png_creator(img+'.residual', range_list=[0.3, 1.0])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.0362096913159,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00253091799095,
            'min_val_pos':[51, 36, 0, 0],
            'im_rms': 0.00317099603729,
            'im_sum': 1.72629857491,
            'regn_sum': 1.70481428877,
            'npts_real': 6400,
            'im_fit': [[40.2022706573, 40.0784833662],
                        0, 107.840244976, 0.0368173095435,
                        17.888484296, 9.90872728645]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 74, 0, 0]), \
                (img+'.image', False, [40, 75, 0, 0]), \
                (img+'.image', True, [6, 40, 0, 0]), \
                (img+'.image', False, [5, 40, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_a1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 6400,
            'freq_bin': 15849921197.895538,
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'mask_pix': 334,
            'mask_regns': 1,
            'npts_real': 6400}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.200896695256,
            'min_val_pos':[33, 6, 0, 0],
            'npts_0.2': 3793,
            'npts_0.5': 1813,
            'npts_real': 6400,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.170877560973,
            'min_val_pos':[45, 14, 0, 0],
            'im_rms':  0.11175375188,
            'im_sum': 13.431971317,
            'npts_real': 6400,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.00680132163689,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00267858314328,
            'min_val_pos':[51, 36, 0, 0],
            'im_rms': 0.00119958583186,
            'im_sum': 0.167714220393,
            'npts_real': 6400}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 3208318.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 3208318.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[-0.003, 0.04])
            self.png_creator(img+'.residual', range_list=[-0.003, 0.04])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.0375871881843,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00482273567468,
            'min_val_pos':[40, 62, 0, 0],
            'im_rms': 0.00341863379853,
            'im_sum': 1.82207681158,
            'regn_sum': 1.77460362448,
            'npts_real': 6400,
            'im_fit': [[40.4007657574, 39.9910814072],
                        0, 107.840244976, 0.0382050339329,
                        18.0250373345, 9.89487712402]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [40, 74, 0, 0]), \
                (img+'.image.tt0', False, [40, 75, 0, 0]), \
                (img+'.image.tt0', True, [6, 40, 0, 0]), \
                (img+'.image.tt0', False, [5, 40, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_a1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 6400,
            'freq_bin': 15849921197.895538,
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'mask_pix': 332,
            'mask_regns': 1,
            'npts_real': 6400}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.200896695256,
            'min_val_pos':[33, 6, 0, 0],
            'npts_0.2': 3793,
            'npts_0.5': 1813,
            'npts_real': 6400,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.170877560973,
            'min_val_pos':[45, 14, 0, 0],
            'im_rms':  0.11175375188,
            'im_sum': 13.431971317,
            'npts_real': 6400,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.00681467168033,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.00269496999681,
            'min_val_pos':[51, 36, 0, 0],
            'im_rms': 0.00121173798212,
            'im_sum': 0.146451243258,
            'npts_real': 6400}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 3208318.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 3208318.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image.tt0', range_list=[-0.005, 0.04])
            self.png_creator(img+'.residual.tt0', range_list=[-0.005, 0.04])

            test_dict[testname]['images'] = \
                [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 3.544453e+11,
            'end': 3.545672e+11,
            'nchan': 1000,
            'max_val': 3.06201076508,
            'max_val_pos':[46, 41, 0, 491],
            'min_val': -0.362504482269,
            'min_val_pos':[59, 39, 0, 494],
            'im_rms': 0.055973852238,
            'rms_per_chan': [0.04147284042825574, 0.04254928243180794, 0.046807795799524314, 0.05027397223185615, 0.04726460222419306, 0.041559297830375534, 0.04236561722381218, 0.043036468623728126, 0.04185685692088316, 0.045531459154953734, 0.0502720842227908, 0.046580946962388896, 0.04349004222357965, 0.04088516731949348, 0.04032681058965125, 0.04004370630764985, 0.04282532717338722, 0.04243048499252981, 0.04449161973003974, 0.04324216366829654, 0.04616183427621934, 0.045546850381536906, 0.04568917172164622, 0.04878374160870331, 0.04532779077099853, 0.0437936219389287, 0.04390807637740492, 0.048466189860678285, 0.046529104851415704, 0.04044114715260885, 0.04105299548884187, 0.04718021385024875, 0.048635734842170555, 0.04778352173735742, 0.048864397284943815, 0.04947867522897103, 0.05141910763698068, 0.0512721608744214, 0.0473725195650241, 0.043381703225800525, 0.04474353375743355, 0.04225543218841734, 0.04182061034187776, 0.04143552743108357, 0.03986439693964313, 0.045028854960487516, 0.04438820259973242, 0.042863745521477134, 0.04849035569380246, 0.04884082602653869, 0.044756041531294664, 0.04599423395626361, 0.049482043151084176, 0.0496152044643088, 0.045281948947307755, 0.043645237105021745, 0.04347131672008097, 0.044574875397560484, 0.04510197212386486, 0.0473453186056419, 0.046313552258314065, 0.0472346150744253, 0.043175615839884204, 0.041733352765603186, 0.048210122453321784, 0.04932141826597058, 0.048629347551728284, 0.048537272402429395, 0.047245553989599516, 0.04574888610875129, 0.04270406001206756, 0.046658923904596826, 0.05038738811398297, 0.050753471866074394, 0.05332664694276341, 0.0500212682167146, 0.04005595830599356, 0.04256708608554199, 0.04674317187094298, 0.04220196590717436, 0.042417250547235955, 0.04638325160660877, 0.04977411965178434, 0.048010039596647354, 0.04729628659810954, 0.04552976930765748, 0.045897452383889176, 0.05011848864926129, 0.049265748843579454, 0.044939308431170936, 0.04173942757602199, 0.043309993723027386, 0.04481787437052452, 0.04593964591705647, 0.042168150769080266, 0.043064963635774126, 0.044709873996171136, 0.04705284280702873, 0.0462879042533439, 0.04793783471906006, 0.04795045974722719, 0.043234046665818066, 0.04657259303546827, 0.05060673222029783, 0.04982265694059221, 0.04945810522631542, 0.04718691129436864, 0.04599097717885053, 0.04652832978367185, 0.052507053711751205, 0.051181773441787363, 0.04511396852603947, 0.04586263498786632, 0.040426992884568934, 0.04110525614812427, 0.04356376644251651, 0.04452582713745292, 0.045080577140126664, 0.04593602103694487, 0.044815607997532654, 0.04513071706324825, 0.04809051889261438, 0.05037930616846981, 0.05192179480783372, 0.04521341852721519, 0.03962756432350773, 0.04590142198822248, 0.046515903084460436, 0.04488842262996963, 0.04761712356190545, 0.04794787134003834, 0.045634077003998914, 0.04401564779461644, 0.04075926253790344, 0.046011153342612625, 0.052943227415454486, 0.04586891075704417, 0.041053053549788235, 0.04741848461012887, 0.04760539104537539, 0.04251590135895046, 0.041517475889717256, 0.041530348823189914, 0.04181852835766343, 0.04435229347677239, 0.041494787553196834, 0.04029171963642991, 0.04171074628726203, 0.03944305739155072, 0.03930109753701833, 0.04150470464803849, 0.04361220599575606, 0.04108243710715513, 0.04468045730535825, 0.046349270421574155, 0.04739880380909585, 0.045192504763623735, 0.040344656318779815, 0.04519294126647845, 0.049069650042015946, 0.04867861195647122, 0.05150697532031088, 0.05310825022071752, 0.044183347666423387, 0.04385582347461477, 0.04832622326670505, 0.05039773328491312, 0.04995486556079216, 0.044406026526960116, 0.04224896017060205, 0.04339237711388555, 0.043054962255623365, 0.03815805334382629, 0.039538862002337145, 0.03972032675929341, 0.03983988379595894, 0.04067903646705838, 0.042816131958255746, 0.044500702340071294, 0.046921846157623126, 0.05137318219477987, 0.049149740418630204, 0.043796797416599566, 0.045302958655315764, 0.04174656753685317, 0.04673517759601561, 0.04394365100006975, 0.04337034201775957, 0.03993686958213624, 0.03600549752919795, 0.043866442563853844, 0.04967750426554719, 0.04328006231133492, 0.04429873465092799, 0.04539662341787065, 0.04208124262819621, 0.0419930913030589, 0.043390621818758704, 0.04031071173345182, 0.040404111183998906, 0.045092347850751044, 0.04461721051697477, 0.047339081036115664, 0.04902937369564867, 0.04788705109593518, 0.044187432479277686, 0.04743470609004175, 0.045923903982013736, 0.044487812095947744, 0.04569301579760048, 0.04588298557818862, 0.048879884521756484, 0.042517514688907354, 0.044736166501443096, 0.049966263636027984, 0.05055485639716179, 0.0486507675230859, 0.045900183028439706, 0.04597463115002017, 0.04874118410854506, 0.0492627400756792, 0.04553659078435011, 0.04308247909546547, 0.04661898391374173, 0.046133985226344924, 0.04326784460986761, 0.04608934609263339, 0.04625450203821501, 0.044862609177145935, 0.04011425405357339, 0.041241389971566234, 0.03875820612078441, 0.04402089540797002, 0.04758335721211567, 0.0494907576224672, 0.04651213583236676, 0.04263585154588318, 0.04172307447971044, 0.04349353240515943, 0.04139629311641164, 0.044157035731005885, 0.04641142701759145, 0.04776945914094206, 0.04658322019166245, 0.04808416923554729, 0.04936108762937629, 0.04552844004064125, 0.04525718576147768, 0.050084592913668284, 0.050497649572359075, 0.04819990615260935, 0.043508765026292744, 0.04331538798826631, 0.04182089077294588, 0.04060874225805633, 0.04479726493754151, 0.047444774373354884, 0.04987173103462332, 0.050304956312335815, 0.042820208031155285, 0.046200529213358, 0.049547291267549205, 0.0439794169450193, 0.041240659042780806, 0.043480246063883995, 0.0436023938949218, 0.04327850965717666, 0.04376974285821181, 0.04364660222066811, 0.04336487246668942, 0.04202446629715709, 0.0446809233205639, 0.04429657923257662, 0.04586834962102086, 0.04777632595404274, 0.04711156800879183, 0.04588720887454971, 0.04171637723240844, 0.04258842902985198, 0.043791244321268845, 0.04444364218656609, 0.04703530623267712, 0.046967155407046726, 0.04423379398880278, 0.046597813073495226, 0.04515926828622885, 0.04865114788979305, 0.04959718488366777, 0.04775289058645823, 0.048222751666638024, 0.045807848861385864, 0.04258228402458693, 0.043295652945759516, 0.040739713165586305, 0.04112738264439764, 0.045910299621454445, 0.04599175115680861, 0.044287782495768244, 0.04662480181390475, 0.04529554134253019, 0.046560456501458534, 0.04506625139233557, 0.04502184582085042, 0.042166746741141986, 0.0411152754577657, 0.03926189149501848, 0.039107527777057026, 0.04681587804048637, 0.0467105611576896, 0.046323540451818135, 0.047918802925478454, 0.04657382422156662, 0.04444849776991001, 0.04346268753310575, 0.044898397078793936, 0.046510495726545836, 0.04619833742205952, 0.0443555071018013, 0.04635832632481224, 0.041976743688728786, 0.04197602077685787, 0.04478463724471771, 0.047207294982882224, 0.05136573873750946, 0.04973287522222192, 0.045711977533150296, 0.040456014664830654, 0.03924539842104608, 0.046609134010290224, 0.04733368546103009, 0.04334939035917689, 0.04290796909827518, 0.045276995722941284, 0.046520486617919486, 0.04815054390965692, 0.04350658880383937, 0.04483676528017031, 0.04077733169883272, 0.041984761989770365, 0.040793755345831954, 0.04472020759595963, 0.04204742196256374, 0.04600734863246967, 0.05089399077936206, 0.0476489534013462, 0.0456432877701943, 0.044223470648491386, 0.04530694593376634, 0.04337233068228447, 0.042694682962867124, 0.044726289013241535, 0.043300741131810816, 0.04210560761472504, 0.041739307544568935, 0.04189755601252967, 0.040777149672434856, 0.04220437809752992, 0.04419949695012377, 0.04359320445482362, 0.04527755744850593, 0.046697574147745516, 0.04642331922810677, 0.043483738542462264, 0.04368515330767118, 0.04908300908360263, 0.04250379298118959, 0.04768927918228166, 0.0486006267244303, 0.049653949354931635, 0.04641633776412558, 0.04737868673319071, 0.046764416207862085, 0.04461670467967515, 0.04868254895425987, 0.0476770300145769, 0.04107924608973292, 0.04188106955545214, 0.041590035586901046, 0.04214482902672925, 0.045136435370418934, 0.04573554214962054, 0.043710136784232044, 0.04634140816704903, 0.04404580814065792, 0.038977016608844756, 0.04117627295447821, 0.04243400270559211, 0.04363716717613951, 0.04419958014197212, 0.04782353317334333, 0.04906799379537854, 0.0512348649962232, 0.046385535687935825, 0.0403796970773063, 0.04050700815126366, 0.040063995010656335, 0.04031612066398897, 0.04164336103806869, 0.046497487361936596, 0.04562422571208214, 0.04263755686426088, 0.04483032698151463, 0.04587078156764805, 0.04362771229290711, 0.04381177070732172, 0.04620332641463225, 0.04248238362177694, 0.04620226614206356, 0.04879854913727185, 0.04526347063166013, 0.04465446848896567, 0.049905242523157714, 0.044088661998105946, 0.04208975459552379, 0.04509748114630728, 0.047764989805243976, 0.047269720611939715, 0.04815713354608546, 0.04861013081395479, 0.047923545468142856, 0.043222858729962114, 0.04395516586340558, 0.04633226277208989, 0.04633822508846244, 0.04466181106665416, 0.044429807483860266, 0.043064439123369755, 0.0406185164405471, 0.04336445242364316, 0.04935831428699568, 0.051896146917270104, 0.04800247160267496, 0.0451928435647402, 0.04857713796305899, 0.04545063178306572, 0.04804396178280432, 0.05062919060884726, 0.04786223876583842, 0.04566827171942874, 0.043626694816840776, 0.041250182796215465, 0.04709180830104477, 0.04772551420079806, 0.044386072616289254, 0.03820143940422228, 0.04413829408333633, 0.04230936833406113, 0.04910812334129778, 0.049843192718712584, 0.048737041668602155, 0.04588920688875621, 0.03628094221396932, 0.0371987923390054, 0.04471470429430708, 0.050663574304905326, 0.043588998319089856, 0.04221773793026161, 0.04144727099823885, 0.044023383956817666, 0.04470813033701908, 0.04269963975357676, 0.042875526590775696, 0.04317741059571664, 0.042984637475138296, 0.04319761839249863, 0.042861390095655066, 0.0452468522274471, 0.049152419395045284, 0.0469133600403256, 0.04228455356894668, 0.0390540855759181, 0.04197843679175804, 0.049840238850523244, 0.04706824869796128, 0.042537650152546964, 0.047732161147428794, 0.04341518757549293, 0.0476977950606753, 0.04549419041632738, 0.03875515499440533, 0.042797424424739626, 0.04596345837561548, 0.0438363052201628, 0.04685885716066778, 0.04283950332858225, 0.044488182626806164, 0.05995010602561655, 0.10557374827025408, 0.18872692387137646, 0.28302762974690815, 0.3498003625351968, 0.3672106885356358, 0.33194227396054365, 0.290440560615709, 0.26012873399840286, 0.24803777375051594, 0.24978939300892694, 0.25141736578190277, 0.25421665255928555, 0.26251178802202624, 0.25890228013598743, 0.23200587051623178, 0.18629342623098122, 0.1303196902675577, 0.07996660345798334, 0.04899539873406976, 0.04090759649654674, 0.043628462817536906, 0.05009979123148561, 0.04956125914220415, 0.044080405061950774, 0.04478967661369681, 0.04570876170701758, 0.04409156219707164, 0.04304988904995974, 0.04781018261476629, 0.04835789663918392, 0.04507151460403614, 0.046165905299855176, 0.04519394844018466, 0.04332556952868729, 0.04057131126983464, 0.04298739652919473, 0.04432500568724349, 0.04972848221531476, 0.049019669301815824, 0.04668382236769563, 0.04472111573456095, 0.045103852453247786, 0.046491820979471625, 0.05029817544174287, 0.050276552117727094, 0.04607120553921124, 0.046587396318693006, 0.04654908671975277, 0.0477396384141432, 0.05256856375387284, 0.05309433276083172, 0.04936056085087188, 0.04465264498510182, 0.04066864955986269, 0.038536803593856646, 0.044266520348999756, 0.044096677361870014, 0.04390677053924258, 0.04488540119711562, 0.04541089658069488, 0.04699651590931532, 0.04734119093292678, 0.04514869425580925, 0.04047704621692836, 0.044688492688345544, 0.04221500509145915, 0.04450735469144093, 0.042870304238128605, 0.045817794483138954, 0.04708641472548149, 0.040032566176785585, 0.04085404222004716, 0.04308076248218408, 0.04352225643436645, 0.04455733495127243, 0.04627376965000834, 0.04544360293708542, 0.04476633455817276, 0.04422718765041479, 0.043247403149073264, 0.04420384687399638, 0.04400101820916508, 0.044306358418837066, 0.04445515506908776, 0.04650940955011436, 0.042975867175620554, 0.04506894253575688, 0.04815973017230859, 0.046770673106161904, 0.04723785961435044, 0.04722458534100488, 0.04635888252639177, 0.04754855183877157, 0.04541021786961315, 0.045993677503893035, 0.04879607967477977, 0.043012398343036544, 0.04339187856746779, 0.04716829159128152, 0.044042929942474544, 0.04060802127590302, 0.03933483270377363, 0.0421084592357695, 0.045799139539569794, 0.04728564817775729, 0.048709802287630354, 0.04625098222367802, 0.04142378728333036, 0.04654172737077701, 0.04695454123756478, 0.04254219696670264, 0.04167825362567383, 0.043590015524733454, 0.04517446609996377, 0.04538394511111757, 0.0411803245006179, 0.04162769198992844, 0.04566795976989248, 0.0469352857577231, 0.045475823563583345, 0.04670207165848767, 0.04944871386138451, 0.04261661221332014, 0.04350784548879662, 0.04563064656476053, 0.04645374564617096, 0.047011364136469676, 0.04626821704856244, 0.04601759529369392, 0.04692011789955231, 0.046935547183383414, 0.04811322379382826, 0.05092280831305335, 0.05169333173453074, 0.050615414866423104, 0.04965990304886755, 0.04985520414774309, 0.049556207960304625, 0.045668964894689947, 0.04091929343081231, 0.04134479306355263, 0.04851193747165114, 0.04945249500604705, 0.04719240158553954, 0.04446873081259328, 0.04328173884244462, 0.04263930758124629, 0.04538798376013595, 0.04733862004263131, 0.042818994443588124, 0.03991049287167168, 0.04003339147973695, 0.045131965116429994, 0.044555886972228564, 0.03904228167420214, 0.03964296000593127, 0.0398357586726478, 0.04462944027508994, 0.047207297513162834, 0.04720091701217782, 0.04573250672526995, 0.050619550625776444, 0.05167748038945954, 0.05037746570705263, 0.04136834989720388, 0.0395474086935746, 0.04275837492827358, 0.0439212218316328, 0.0471906556974699, 0.050993019281882726, 0.049124147702145744, 0.048244397599512094, 0.04506197168417649, 0.04568375276302797, 0.04812036783271446, 0.04445936319441292, 0.041293535244715976, 0.043179887675183666, 0.046055289083475205, 0.04706507397502584, 0.04447967606146786, 0.04222015492968105, 0.0439674462328621, 0.043896529077583894, 0.04266544335030779, 0.04266644623124096, 0.04357712294408853, 0.05231398531202364, 0.053370700541464616, 0.04857282762175127, 0.045258380741403106, 0.045303357380300946, 0.04101809933111672, 0.041189074434246936, 0.04633904453435081, 0.049717255912949114, 0.05010614817105909, 0.04393734175641571, 0.03867286537447114, 0.0429230714573404, 0.04488217279877273, 0.04522370048445455, 0.04336725660321416, 0.04356264989657724, 0.04146827380306933, 0.03665312851533432, 0.039980260984851516, 0.04228302535435613, 0.04389068048602057, 0.03997240613607379, 0.04241841231637553, 0.04160034918348827, 0.04520530739905513, 0.04133848236544517, 0.04417032186418197, 0.04732307233198403, 0.046066585583925186, 0.04794749014534085, 0.0430060722593525, 0.04177280157123228, 0.04380811047918848, 0.040779536741100096, 0.044002495409217926, 0.042339960518415154, 0.04284779608368661, 0.04144915330440946, 0.04175329012568266, 0.04468389809649514, 0.04262292048772433, 0.04101753248574084, 0.03833202726106456, 0.041432154426680995, 0.04545721669426245, 0.043500403222287314, 0.04247370233471467, 0.042711987628776316, 0.045197003435745016, 0.046061750636820546, 0.04279507295911485, 0.04239066668592003, 0.04714143730922104, 0.04661255729559933, 0.04238422932818252, 0.04179816251339896, 0.04481205604751281, 0.04862619399232404, 0.04290828113838176, 0.04221396938127627, 0.04242662994290929, 0.04450718119590875, 0.04471420159102595, 0.041561459460195016, 0.04061336003877967, 0.041745476588896754, 0.04286446946892797, 0.04403000667460554, 0.042516346869012976, 0.04391550291853197, 0.04147492526437, 0.044390159913787586, 0.04858618926975879, 0.047057979274585916, 0.044837316309386074, 0.045456397939683736, 0.046873185359702536, 0.04545495694981475, 0.04359060042114437, 0.045989007837274644, 0.049686754595291746, 0.04536022945930599, 0.04770465331936742, 0.0510614831423978, 0.04590996154799734, 0.04601963242407712, 0.04619596080750725, 0.045306980051320944, 0.04432076077286921, 0.038777354135178126, 0.039785511964539784, 0.042794096518758616, 0.044378135486347496, 0.040378706840431426, 0.04270144914894412, 0.04485723781386858, 0.04092189619423321, 0.044567443137470936, 0.046683929341345125, 0.0446591792962899, 0.04663066745026505, 0.04778417013088139, 0.0479986312062676, 0.04691279172857023, 0.04550992586868622, 0.04481415989096017, 0.046264411558386266, 0.04684048149792687, 0.04508317803023017, 0.04418840392406166, 0.04514330368780914, 0.04331853163851515, 0.041418043475806104, 0.044335444134957894, 0.043813612310493855, 0.042237129336115634, 0.04484883967809909, 0.04719279649342131, 0.048484876497386795, 0.04999861337401947, 0.047139021664050894, 0.04624424651115784, 0.04946026305108584, 0.048753425401869756, 0.04795785902297371, 0.04548142114132507, 0.04597234810107064, 0.04578669626370894, 0.044150493323153066, 0.04347549735022207, 0.045012961497875516, 0.042469833670680926, 0.043123635922132565, 0.0444188259911608, 0.044045252309259125, 0.041347516913099874, 0.04423993898811506, 0.04833449411104301, 0.04399536175530263, 0.04193561651797477, 0.047193999173540044, 0.04640949391019046, 0.04377876735567927, 0.04784563639430407, 0.04990155675893534, 0.04705787638590471, 0.04776353982103787, 0.05032860054279513, 0.04713583458148935, 0.041586710934546633, 0.04081836062004791, 0.046730461342141934, 0.04770022113801226, 0.041521491950066124, 0.04322288824535129, 0.03871971212467864, 0.03877418394769848, 0.04189861595321477, 0.04379541870491918, 0.045818379140407366, 0.04553142481385117, 0.04557339175358304, 0.040221668383312335, 0.03847824536287222, 0.04226084888225819, 0.04278014077365339, 0.04682600952614384, 0.048543746308069145, 0.04672005940927673, 0.04372755168922572, 0.043360368499897335, 0.0406401455261877, 0.04106608036483906, 0.03901295012577516, 0.04432469216665788, 0.045229796180029196, 0.043080105650654996, 0.04000644096010989, 0.041872795599633035, 0.04367030740721623, 0.04062843867079192, 0.04046255587943874, 0.04189940339714176, 0.042461327527919236, 0.04502541852078475, 0.04357770467605546, 0.041638598996523384, 0.04244333432938209, 0.04597283261976164, 0.04870528030372, 0.046562672932188146, 0.04048613501112604, 0.043194916055069865, 0.04705216547916007, 0.04707229565946042, 0.04308733577727289, 0.042290353014862815, 0.04191029758476797, 0.043622005420342246, 0.04629853652312584, 0.04812334763810347, 0.0454554723964637, 0.044883953341833926, 0.04616194471799996, 0.04146985976458868, 0.04493084061523062, 0.04648662940955314, 0.03901313717168857, 0.04093255386975347, 0.04531066917587298, 0.04790801736759378, 0.043300807696865916, 0.03795010478256023, 0.03937132546536799, 0.03998005844569648, 0.044712329482512296, 0.04344663738412382, 0.04384624678585034, 0.04189198462079245, 0.0447850884961132, 0.04360860246046017, 0.046212492689080946, 0.04662574746696418, 0.040175151217488075, 0.042548892150321385, 0.04273708134110358, 0.0436652394233731, 0.0446567273402838, 0.041504006631368794, 0.04180817535232774, 0.04180189242542836, 0.04464597692689543, 0.04051048022670634, 0.0402043786014598, 0.044332536061165495, 0.04451349302294642, 0.04721775825825243, 0.05034142046308038, 0.04914557159307879, 0.04305314404299413, 0.04501919188309524, 0.04385981707446499, 0.04026426160488087, 0.043707857906778516, 0.0434888399935484, 0.03994276460009276, 0.040291017821495335, 0.04551819051577132, 0.04785682015877868, 0.048275626437799854, 0.04494963910113805, 0.042300285689887605, 0.045716179659770935, 0.046907677550574234, 0.042150397450181565, 0.0417516956499194, 0.04703357874043769, 0.04323223727854463, 0.048300785507211216, 0.051514165729642306, 0.045485370078021274, 0.03691851661101719, 0.035216417938379464, 0.04073572180872651, 0.03847618248085107, 0.03924093325651294, 0.04169824170116127, 0.041857460173997246, 0.042763224613731665, 0.04363214657035114, 0.04600106990967965, 0.051028150668357815, 0.045746850061159285, 0.04206677455142946, 0.048022086691811045, 0.044914251190311984, 0.042020771247839775, 0.04626404511121505, 0.04497765856342021, 0.04481991912280334, 0.04506179626666786, 0.04985039232547424, 0.05054033486622929, 0.043881117212786314, 0.045088403501057986, 0.041061363556784206, 0.03711119731324331, 0.039768631455522105, 0.0477101845505575, 0.04633349558816296, 0.045643402391453314, 0.04548923525790922, 0.04656460145077493, 0.04989878555412256, 0.05033863381352162, 0.046023336890052224, 0.045066109334492205, 0.04464246708556863, 0.04342978956885089, 0.04410532977158723, 0.04664016435292943, 0.04798874070999172, 0.044607643307464015, 0.043805868432151375, 0.045124615570696955, 0.04565951172721549, 0.046829926257196015, 0.044668867909129516, 0.05062380133045267, 0.051219538639712014, 0.04842948233786371, 0.04432788396188177, 0.04266337961965428, 0.041294265588724774, 0.043132576757287024, 0.044563631847001746, 0.04507019653789561, 0.04458375243067548, 0.04125790045420121, 0.04146933310219249, 0.04236381054760312, 0.04492870140983278, 0.04750546903562345, 0.04511795490590616, 0.04460467237817201],
            'im_sum': 2249.5560029,
            'regn_sum': 253.368211955,
            'npts_real': 6400000,
            'im_fit': [[45.802306052, 41.1650012396],
                        491, 354.505209348, 2.89344750697,
                        6.25839481605, 5.74536585146]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( \
            self.check_list_vals(im_stats_dict['rms_per_chan'], \
                exp_im_stats['rms_per_chan'], epsilon=0.01), True, \
            valname='RMS per channel of .image', exact=True)
        out, report1_u = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 6400000,
            'freq_bin': 122071.64398193359,
            'start': 3.544453e+11,
            'end': 3.545672e+11,
            'nchan': 1000,
            'mask_pix': 9932,
            'mask_regns': 1,
            'npts_real': 6400000}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 3.544453e+11,
            'end': 3.545672e+11,
            'nchan': 1000,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': 0.20036059618,
            'min_val_pos':[39, 8, 0, 810],
            'npts_0.2': 3233000,
            'npts_0.5': 1549000,
            'npts_real': 6400000,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 3.544453e+11,
            'end': 3.545672e+11,
            'nchan': 1000,
            'max_val': 1.0,
            'max_val_pos':[40, 40, 0, 0],
            'min_val': -0.164189130068,
            'min_val_pos':[36, 35, 0, 993],
            'im_rms':  0.0871161935921,
            'im_sum': 2742.74484326,
            'npts_real': 6400000,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 3.544453e+11,
            'end': 3.545672e+11,
            'nchan': 1000,
            'max_val': 0.366728395224,
            'max_val_pos':[39, 68, 0, 502],
            'min_val': -0.338401287794,
            'min_val_pos':[32, 47, 0, 493],
            'im_rms': 0.0469751110358,
            'im_sum': 242.908819752,
            'npts_real': 6400000}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 3.544453e+11,
            'end': 3.545672e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 3.544453e+11,
            'end': 3.545672e+11,
            'nchan': 1000,
            'max_val': 1009.50134277,
            'max_val_pos':[0, 0, 0, 980],
            'min_val': 1007.94287109,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1000}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[0.0, 3.25])
            self.png_creator(img+'.residual', range_list=[0.0, 3.25])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png', img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        tclean(vis=self.msfile, field='0', spw=['0:245.220516619'
            '~245.273983416GHz,1:261.752937691~261.774177925GHz;261.783699409'
            '~261.837898628GHz;261.958504097~261.984871284GHz'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'
            '21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,'
            '41,42,43,44,45,46'], scan=['7,11'], \
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
        self.copy_products(file_name+'0', file_name+'1')

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='0', spw=['0:245.220516619'
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

        report0 = th.checkall( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.03113353252,
            'max_val_pos':[224, 153, 0, 0],
            'min_val': -1.01794064045,
            'min_val_pos':[222, 93, 0, 0],
            'im_rms': 0.359352011299,
            'im_sum': -1491.198136,
            'regn_sum': 3362.95355159,
            'npts_real': 82944,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        0, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [144, 266, 0, 0]), \
                (img+'.image', False, [144, 267, 0, 0]), \
                (img+'.image', True, [22, 145, 0, 0]), \
                (img+'.image', False, [21, 145, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_a1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 82944,
            'freq_bin': 16762501225.396851,
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 82944}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': 0.200061768293,
            'min_val_pos':[114, 25, 0, 0],
            'npts_0.2': 47325,
            'npts_0.5': 22362,
            'npts_real': 82944,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        0, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': -0.0609973333776,
            'min_val_pos':[140, 137, 0, 0],
            'im_rms':  0.0198315591613,
            'im_sum': 16.3422700718,
            'npts_real': 82944,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        0, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.03113353252,
            'max_val_pos':[224, 93, 0, 0],
            'min_val': -0.676309525967,
            'min_val_pos':[54, 126, 0, 0],
            'im_rms': 0.359352011299,
            'im_sum': -1491.198136,
            'npts_real': 82944}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 23234590.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 23234590.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[-1.05, 1.05])
            self.png_creator(img+'.residual', range_list=[-1.05, 1.05])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        tclean(vis=self.msfile, field='0', spw=['0:245.220516619'
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
        self.copy_products(file_name+'0', file_name+'1')

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='0', spw=['0:245.220516619'
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

        report0 = th.checkall( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.01515209675,
            'max_val_pos':[224, 153, 0, 0],
            'min_val': -1.01363265514,
            'min_val_pos':[222, 93, 0, 0],
            'im_rms': 0.361662749039,
            'im_sum': -1558.29933771,
            'regn_sum': 3388.57268598,
            'npts_real': 82944,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        0, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [144, 266, 0, 0]), \
                (img+'.image.tt0', False, [144, 267, 0, 0]), \
                (img+'.image.tt0', True, [22, 145, 0, 0]), \
                (img+'.image.tt0', False, [21, 145, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_a1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 82944,
            'freq_bin': 16762501225.396851,
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 82944}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': 0.200061768293,
            'min_val_pos':[25, 114, 0, 0],
            'npts_0.2': 47325,
            'npts_0.5': 22362,
            'npts_real': 82944,
            'pb_im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report3_r = th.check_val(im_stats_dict['pb_im_fit'][0][0], \
            exp_im_stats['pb_im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report3_s = th.check_val(im_stats_dict['pb_im_fit'][0][1], \
            exp_im_stats['pb_im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report3_t = th.check_val( \
            im_stats_dict['pb_im_fit'][1], exp_im_stats['pb_im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report3_u = th.check_val( \
            im_stats_dict['pb_im_fit'][2], exp_im_stats['pb_im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report3_v = th.check_val(im_stats_dict['pb_im_fit'][3], \
            exp_im_stats['pb_im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report3_w = th.check_val(im_stats_dict['pb_im_fit'][4], \
            exp_im_stats['pb_im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report3_x = th.check_val(im_stats_dict['pb_im_fit'][5], \
            exp_im_stats['pb_im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report3 = report3_a + report3_b + report3_c + report3_d + \
            report3_e + report3_f + report3_g + report3_h + report3_i + \
            report3_j + report3_k + report3_l + report3_m + report3_n + \
            report3_o + report3_p + report3_q + report3_r + report3_s + \
            report3_t + report3_u + report3_v + report3_w + report3_x

        # .psf report
        psf_stats_dict = self.image_stats(img, '.psf.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.psf.tt0.crtf')

        exp_psf_stats = {'npts': 82944,
            'npts_unmasked': 82944.0,
            'freq_bin': 16762501225.396851,
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[144, 144, 0, 0],
            'min_val': -0.0609973333776,
            'min_val_pos':[140, 137, 0, 0],
            'im_rms':  0.0198315591613,
            'im_sum': 16.3422700718,
            'npts_real': 82944,
            'cen_im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report4_r = th.check_val(im_stats_dict['cen_im_fit'][0][0], \
            exp_im_stats[cen_'im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report4_s = th.check_val(im_stats_dict['cen_im_fit'][0][1], \
            exp_im_stats['cen_im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report4_t = th.check_val( \
            im_stats_dict['cen_im_fit'][1], exp_im_stats['cen_im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report4_u = th.check_val( \
            im_stats_dict['cen_im_fit'][2], exp_im_stats['cen_im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report4_v = th.check_val(im_stats_dict['cen_im_fit'][3], \
            exp_im_stats['cen_im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report4_w = th.check_val(im_stats_dict['cen_im_fit'][4], \
            exp_im_stats['cen_im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report4_x = th.check_val(im_stats_dict['cen_im_fit'][5], \
            exp_im_stats['cen_im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report4 = report4_a + report4_b + report4_c + report4_d + \
            report4_e + report4_f + report4_g + report4_h + report4_i + \
            report4_j + report4_k + report4_l + report4_m + report4_n + \
            report4_o + report4_p + report4_q + report4_r + report4_s + \
            report4_t + report4_u + report4_v + report4_w + report4_x

        # .residual report
        resid_stats_dict = self.image_stats(img, '.residual.tt0', \
            region_file = data_path+
            'region_files/standard_mtmfs_eph.residual.tt0.crtf')

        exp_resid_stats = {'npts': 82944,
            'npts_unmasked': 47329.0,
            'freq_bin': 16762501225.396851,
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.03113353252,
            'max_val_pos':[224, 153, 0, 0],
            'min_val': -1.01794064045,
            'min_val_pos':[222, 93, 0, 0],
            'im_rms': 0.359352011299,
            'im_sum': -1491.198136,
            'npts_real': 82944}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 23234590.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 23234590.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image.tt0', range_list=[-1.05, 1.05])
            self.png_creator(img+'.residual.tt0', range_list=[-1.05, 1.05])

            test_dict[testname]['images'] = \
                [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 2.20301e+11,
            'end': 2.20301e+11,
            'nchan': 1,
            'max_val': 2.40606927872,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': -0.0115160318092,
            'min_val_pos':[50, 16, 0, 0],
            'im_rms': 0.203888800435,
            'im_sum': 172.688482511,
            'regn_sum': 171.572869264,
            'npts_real': 8100,
            'im_fit': [[45.000405766, 45.0014155577],
                        0, 220.300765422, 2.40974849537,
                        9.96002749264, 4.61946099469]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [45, 85, 0, 0]), \
                (img+'.image', False, [45, 86, 0, 0]), \
                (img+'.image', True, [5, 45, 0, 0]), \
                (img+'.image', False, [4, 45, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_a1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 8100,
            'freq_bin': 125009872.91876221,
            'start': 2.20301e+11,
            'end': 2.20301e+11,
            'nchan': 1,
            'mask_pix': 404,
            'mask_regns': 1,
            'npts_real': 8100}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 2.20301e+11,
            'end': 2.20301e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': 0.200092822313,
            'min_val_pos':[36, 6, 0, 0],
            'npts_0.2': 5041,
            'npts_0.5': 2409,
            'npts_real': 8100,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 2.20301e+11,
            'end': 2.20301e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': -0.186672970653,
            'min_val_pos':[31, 41, 0, 0],
            'im_rms':  0.125651854223,
            'im_sum': 40.2140428556,
            'npts_real': 8100,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 2.20301e+11,
            'end': 2.20301e+11,
            'nchan': 1,
            'max_val': 0.0233333036304,
            'max_val_pos':[45, 45, 0, 0],
            'min_val': -0.0115160560235,
            'min_val_pos':[50, 16, 0, 0],
            'im_rms': 0.00410021113948,
            'im_sum': 0.122163831366,
            'npts_real': 8100}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 2.20301e+11,
            'end': 2.20301e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 2.20301e+11,
            'end': 2.20301e+11,
            'nchan': 1,
            'max_val': 201537.8125,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 201537.8125,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report7 = report7_a + report7_b + report7_c + report7_d + \
            report7_e + report7_f + report7_g + report7_h + report7_i + \
            report7_j + report7_k + report7_l + report7_m + report7_n + \
            report7_o

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7


        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[-0.015, 2.5])
            self.png_creator(img+'.residual', range_list=[-0.015, 2.5])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


###############################################
###############################################

class Test_mosaic(test_tclean_base):


    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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
            img = file_name+'1'
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

            img = file_name+'1'

        report0 = th.checkall( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 1.18474924564,
            'max_val_pos':[44, 38, 0, 252],
            'min_val': -0.417252391577,
            'min_val_pos':[24, 31, 0, 498],
            'im_rms': 0.0878815938952,
            'rms_per_chan': [0.08352891852269921, 0.08135269305379152, 0.08395260930289478, 0.08279304060774786, 0.08596507857262789, 0.08278415992796316, 0.08377357600543471, 0.08394694956313518, 0.08294962386707326, 0.08607043136961237, 0.07680081900457504, 0.09496464679558769, 0.08214392421353917, 0.0857648739112381, 0.08281633996127867, 0.0894806463001662, 0.081589405081739, 0.07924880491236708, 0.07453379448699567, 0.09379627665607704, 0.08157315058137742, 0.09166424964725516, 0.08450533906882468, 0.08106955380637881, 0.08737684423508374, 0.09073897650010625, 0.09084105654348026, 0.07784857687669093, 0.08502491422043724, 0.08973744200475807, 0.08225747757890416, 0.08302380899013398, 0.0770200580373865, 0.08881158355253989, 0.08865534449297038, 0.09159308234217008, 0.08115203022056021, 0.08341851059515137, 0.08987845777172865, 0.08508065238113408, 0.07778307167083158, 0.086558262984569, 0.08338769817916242, 0.08071653499000843, 0.08865177850010854, 0.08534950297770955, 0.07616257958327323, 0.09190545934190404, 0.08895841881770002, 0.0896694261175686, 0.07893346400403167, 0.07455014243908499, 0.09047688550066521, 0.0869112050227715, 0.09350825779882399, 0.07646177066867972, 0.08249014753514353, 0.09101549816292659, 0.08639556025249603, 0.08038246669695737, 0.08824074457422779, 0.08618380562988853, 0.08205108347689954, 0.07771920892371241, 0.08768241197771068, 0.08334623972314524, 0.08815925813978431, 0.08212951592340229, 0.0823396048959934, 0.08234068229313843, 0.08884264540759013, 0.07862615163124824, 0.08655155235861448, 0.08979383171683271, 0.08809736503013246, 0.07192760633210062, 0.08230751478411082, 0.07857712439720876, 0.07948313531326559, 0.09225583910089932, 0.0827786011949029, 0.08972941404057228, 0.0887804253347499, 0.088568302489554, 0.08542522263383671, 0.08705651177174693, 0.09503784914334118, 0.09211196203585155, 0.09124273060482811, 0.0894655349173044, 0.08578196967967916, 0.08600001037669011, 0.08725679761935697, 0.09250336443228556, 0.09141071794136318, 0.08257340474387778, 0.08489860170832744, 0.09095982550662123, 0.08386307951475633, 0.08245963060222646, 0.08082452670747917, 0.08421465818141162, 0.07716197318352062, 0.08491216770880909, 0.09219997011780534, 0.08999904369748837, 0.08917413213129059, 0.08669554914308757, 0.08478598632777916, 0.08272607334416665, 0.09204243226438258, 0.08766378811832916, 0.0865825253580717, 0.07758524805681824, 0.09085091409160705, 0.08521446839560726, 0.0884995697340635, 0.0884733008608976, 0.08577677843356687, 0.09306189230759793, 0.08612785329237876, 0.08258623029884975, 0.07941031533209075, 0.08826196881912282, 0.09840895220275017, 0.08933425191997195, 0.08745304230772397, 0.08566154826880104, 0.09365349891994729, 0.09053316121385116, 0.08810691082194215, 0.08630703760915545, 0.08753186529189085, 0.0937019174346115, 0.0773479036533113, 0.08029497918496241, 0.08724287089313901, 0.08790671029158231, 0.08773912124792194, 0.08266091420346706, 0.08132362905839345, 0.09216634845393681, 0.08762278965147091, 0.09332173318661345, 0.07778672535312933, 0.08408442636695551, 0.08553415947786376, 0.08232058120295321, 0.07780988703074035, 0.08409905782923954, 0.07296439116382543, 0.08925116729628191, 0.09008277668363449, 0.09520991051126101, 0.08553135515808247, 0.084718775948537, 0.08883634984779534, 0.08643599841065604, 0.08339630810303196, 0.08862968769223899, 0.0784883898433635, 0.08099618435890771, 0.08565713860724568, 0.08556682581312276, 0.09128248240593151, 0.08827156303897946, 0.0880790536420939, 0.08578294204187285, 0.08362128063154757, 0.09123831674563974, 0.08722115387496024, 0.09621944328499646, 0.09149116203081945, 0.09118103925268266, 0.08347804034781846, 0.0855692865777019, 0.0897847427742151, 0.07186203727120243, 0.08802972450167627, 0.08442626094668859, 0.08723051889453723, 0.08159640903835322, 0.07548809021347838, 0.0915760655325942, 0.1000199147822218, 0.09787395465748581, 0.08693633490200812, 0.08442264879607819, 0.08721652263357231, 0.0880380635040301, 0.08729925477779406, 0.08519856552892391, 0.08975256932981747, 0.07958208410639706, 0.0914096789556792, 0.09811570777169866, 0.07949196007023566, 0.09006036869795747, 0.0912952536965593, 0.09805249071833734, 0.0897004902132416, 0.08952989036422254, 0.0790849762337038, 0.09347465569139526, 0.09000611959279463, 0.08239101076565238, 0.08373667485591837, 0.0871320704983511, 0.08767396777909775, 0.0965335008897428, 0.08633964662006433, 0.08735933635150735, 0.09035634185538645, 0.08468159828070221, 0.09116890211273068, 0.08804541661464536, 0.08484392727135703, 0.08203979237254262, 0.0889071541447407, 0.0889633730061124, 0.08339973306782386, 0.08330574319001677, 0.08871583444712416, 0.09420962372690272, 0.08521243909929632, 0.08355595121547242, 0.0933236552392962, 0.08719802034432139, 0.09062561297486302, 0.08535346593907196, 0.08813919625097401, 0.092913551212774, 0.09227103971954637, 0.08161432828410638, 0.09072439427447547, 0.09244590715316298, 0.08709339445075073, 0.09924553460801067, 0.1002690868580879, 0.08840557566191516, 0.08133715093336598, 0.09149908227822648, 0.07670073484069043, 0.09249672315752175, 0.0782222639219109, 0.0990254821246598, 0.09280158172995194, 0.08730248696252146, 0.09456558538062404, 0.10081497167443383, 0.10403855973688014, 0.12877330926044459, 0.1386207992219812, 0.13938585037552162, 0.1235720046054832, 0.12067870584152206, 0.10645742986583119, 0.08552529466321122, 0.08268301454232765, 0.08542822040617518, 0.09335660100821405, 0.08885775208502826, 0.08644558884485783, 0.08440970509099167, 0.091904188949323, 0.09288900294073746, 0.09184545798056735, 0.0978475829246572, 0.09536587471324205, 0.09325738676721759, 0.08191614651111201, 0.09130572131132277, 0.09356730058206153, 0.09150948292317468, 0.08930212344785793, 0.08973008889149876, 0.08468876678845778, 0.09439900459047333, 0.08340248011080888, 0.09426625948231673, 0.0846651688286838, 0.08220153827772422, 0.09338524135684238, 0.08541949594391877, 0.07720477952979755, 0.08539185143996587, 0.09204753294875682, 0.08278104837254098, 0.07886155262558238, 0.08341737313272374, 0.08757850027055837, 0.08360325848864149, 0.09166550957557691, 0.09185533846656913, 0.08026617289995454, 0.09106087061249621, 0.08802259852925379, 0.09053259952518417, 0.08328401754435075, 0.08888144950799276, 0.0810576402470788, 0.08479306981872027, 0.08580610257468782, 0.09179491887497933, 0.08841796481950012, 0.0845148944552744, 0.08053481062011857, 0.08862797798354732, 0.09462089152970354, 0.08582787706368114, 0.08814348149204695, 0.08504514772616785, 0.08371224832082137, 0.08678333110045997, 0.09353126562315753, 0.09129854527725832, 0.08116659030415371, 0.09015138809948552, 0.08710081730619218, 0.09437955706280836, 0.08963858464777974, 0.09313459101197238, 0.08856120416820244, 0.08552675606885124, 0.08351176926318361, 0.08411701581891168, 0.08427020020603929, 0.09163946448881417, 0.0805306218916976, 0.08160963806132211, 0.08804552292687956, 0.09626408912773088, 0.08913709670428199, 0.09096834064650154, 0.0851269228240773, 0.09017277009104614, 0.08476290074193281, 0.07632278322336213, 0.08385737538890878, 0.08700039219503956, 0.08866885268530736, 0.08466059774774966, 0.0845759814329557, 0.0790621930867725, 0.08771807117918605, 0.08893473780006282, 0.09310980223541039, 0.09306479580057536, 0.09147420625138089, 0.08949657274281757, 0.08192815006108203, 0.07905463600626796, 0.09666550899459639, 0.0808647476478242, 0.08495044094490133, 0.0916137566838688, 0.09649894883996268, 0.0920414733457368, 0.08356993363476055, 0.09422414928552236, 0.08778457312089456, 0.08987693020950831, 0.09777964977624237, 0.09060848800058789, 0.09087547326073596, 0.09065288590043372, 0.09815595961513002, 0.08630801892602018, 0.08960594520751539, 0.09100452485239247, 0.09259096682095072, 0.09364434916529175, 0.0853051896352123, 0.08391849112594159, 0.08978397560355508, 0.08356817274105656, 0.08639129387321305, 0.07641054760007339, 0.08566749942250246, 0.09110851912136013, 0.08938769699436465, 0.08943111829547572, 0.08605404789230106, 0.08796239383347541, 0.08454535587717901, 0.0929903335078121, 0.08246802178760436, 0.08817700985492735, 0.0820807551149751, 0.08511665149857307, 0.0914822776922275, 0.08779917531640212, 0.0779155145245507, 0.08062449848958794, 0.09151321230973979, 0.08251138633527932, 0.08314391480095229, 0.09660065726688243, 0.09161786920726972, 0.09195890233721427, 0.09484458463309785, 0.08672723967704118, 0.09056091164304954, 0.08950078278038455, 0.08453380213065807, 0.08621663260330842, 0.0903504892183676, 0.08888219947767556, 0.09691781310403562, 0.0829510997618925, 0.08538905047156051, 0.08536033187251706, 0.09253646586834798, 0.08719827400559628, 0.08741965478896539, 0.0908875936865952, 0.08650855583109175, 0.0911287851112432, 0.0870327992529023, 0.09334187790207615, 0.08691128023529447, 0.0829607319773311, 0.08561452123819384, 0.09416699659894374, 0.09865975100963004, 0.08059372543252491, 0.08162290581093083, 0.07969201254174872, 0.09014664917727015, 0.07748434709443736, 0.09115822285795372, 0.0874199097979386, 0.08331094704142918, 0.08450373137759747, 0.0873987436304484, 0.0792090383862253, 0.08682449919890456, 0.08898017363528224, 0.0891014307981212, 0.08578455417679541, 0.09612343808420232, 0.07718957370637691, 0.08963596680253917, 0.09053358289784386, 0.08104077369182848, 0.08805192318672665, 0.09036158074767288, 0.09733534898340712, 0.08642234534217835, 0.0873180423896324, 0.08509746809331085, 0.0927045997121519, 0.08985111493399206, 0.09486348772674255, 0.09000267315635847, 0.08474935185150062, 0.08247809724159613, 0.08802043461258492, 0.0865126082695781, 0.08138258449164787, 0.08795575893476618, 0.09456070784224847, 0.09680636657810467, 0.08476700362040308, 0.08670105726809373, 0.08636724547506389, 0.08412716463066074, 0.08428810620773179, 0.08964151944607554, 0.08689624513108678, 0.08902965594822292, 0.09221339424501332, 0.08726043474359556, 0.08607641544478577, 0.09100554185059015, 0.08492870794208009, 0.08529837352577493, 0.09521158569562824, 0.08914943856267118, 0.087555731639101, 0.0862048336688618, 0.08984721078423315, 0.08217617292259297, 0.08824966695006062, 0.07486261467473272, 0.08753387468056906, 0.08379545796577004, 0.09274777146757962, 0.09220642715156253, 0.0792962124445207, 0.09090807566463247, 0.08751737767113807, 0.07961706129268199, 0.0941224640615722, 0.0795895706794336, 0.09562104758697967, 0.08020847225726233, 0.09417989583892716, 0.09061167014269772, 0.08898710965217106, 0.0897447948736654, 0.08398102899483291, 0.08684215184345169, 0.096630024031122, 0.08473098919932259, 0.09179580145778438, 0.07887094371255345, 0.08638286938225163],
            'im_sum': 365.155811516,
            'regn_sum': 91.3425360965,
            'npts_real': 5925312,
            'rms_per_field': 1,
            'im_fit': [[44.5674514546, 38.386085988],
                        252, 220.31420627, 1.15058925098,
                        9.43393089916, 8.09228180871]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( \
            self.check_list_vals(im_stats_dict['rms_per_chan'], \
                exp_im_stats['rms_per_chan'], epsilon=0.01), True, \
            valname='RMS per channel of .image', exact=True)
        out, report1_u = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = th.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_b1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_e1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1 + report1_e1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 5925312,
            'freq_bin': 244174.08728027344,
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'mask_pix': 3929,
            'mask_regns': 1,
            'npts_real': 5925312}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[54, 54, 0, 0],
            'min_val': 0.200047940016,
            'min_val_pos':[98, 64, 0, 339],
            'npts_0.2': 3338068,
            'npts_0.5': 1825336,
            'npts_real': 5925312,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 1.0,
            'max_val_pos':[54, 54, 0, 0],
            'min_val': -0.168375626206,
            'min_val_pos':[63, 54, 0, 12],
            'im_rms':  0.0604957837605,
            'im_sum': 67.0230730532,
            'npts_real': 5925312,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 0.490073978901,
            'max_val_pos':[49, 60, 0, 249],
            'min_val': -0.417252391577,
            'min_val_pos':[24, 31, 0, 498],
            'im_rms': 0.087537475151,
            'im_sum': 68.1217086753,
            'npts_real': 5925312}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
            'nchan': 508,
            'max_val': 120.900192261,
            'max_val_pos':[0, 0, 0, 447],
            'min_val': 120.665283203,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 508}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
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
            'start': 2.202527e+11,
            'end': 2.203765e+11,
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

        out, report8_a = th.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = th.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = th.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = th.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = th.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = th.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = th.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = th.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = th.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = th.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = th.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = th.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = th.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = th.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = th.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = th.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = th.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = th.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[0.15, 1.2])
            self.png_creator(img+'.residual', range_list=[0.15, 1.2])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.0345157124102,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00193646573462,
            'min_val_pos':[91, 52, 0, 0],
            'im_rms':  0.00202865568675,
            'im_sum': 1.51296002558,
            'regn_sum': 1.58850855171,
            'npts_real': 15876,
            'rms_per_field': 1,
            'im_fit': [[62.9942562846, 62.995885097],
                        0, 107.840245142, 0.0352250058226,
                        17.4609357952, 9.70983031045]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [64, 114, 0, 0]), \
                      (img+'.image', False, [64, 115, 0, 0]), \
                      (img+'.image', True, [11, 60, 0, 0]), \
                      (img+'.image', False, [10, 60, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 15876,
            'freq_bin': 15849925874.83342,
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'mask_pix': 360,
            'mask_regns': 1,
            'npts_real': 15876}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 0.200080409646,
            'min_val_pos':[102, 28, 0, 0],
            'npts_0.2': 8454,
            'npts_0.5': 4497,
            'npts_real': 15876,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.169576942921,
            'min_val_pos':[61, 57, 0, 0],
            'im_rms':  0.0501544137568,
            'im_sum': 0.00266315032371,
            'npts_real': 15876,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.00713972514495,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00197902694345,
            'min_val_pos':[72, 58, 0, 0],
            'im_rms': 0.000868739852256,
            'im_sum': 0.134135425201,
            'npts_real': 15876}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 4396210.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 4396210.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
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

        out, report8_a = th.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = th.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = th.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = th.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = th.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = th.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = th.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = th.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = th.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = th.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = th.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = th.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = th.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = th.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = th.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = th.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = th.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = th.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[-0.002, 0.035])
            self.png_creator(img+'.residual', range_list=[-0.002, 0.035])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.0349145904183,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00258040986955,
            'min_val_pos':[70, 97, 0, 0],
            'im_rms':  0.00207815366878,
            'im_sum': 1.47797238623,
            'regn_sum': 1.58228412218,
            'npts_real': 15876,
            'rms_per_field': 1,
            'im_fit': [[63.0904967336, 62.9480581202],
                        0, 107.840245142, 0.0358052067244,
                        17.191876841, 9.68274896612]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [64, 114, 0, 0]), \
                      (img+'.image.tt0', False, [64, 115, 0, 0]), \
                      (img+'.image.tt0', True, [11, 60, 0, 0]), \
                      (img+'.image.tt0', False, [10, 60, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 15876,
            'freq_bin': 15849925874.83342,
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'mask_pix': 360,
            'mask_regns': 1,
            'npts_real': 15876}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': 0.200080409646,
            'min_val_pos':[102, 28, 0, 0],
            'npts_0.2': 8454,
            'npts_0.5': 4497,
            'npts_real': 15876,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.169576942921,
            'min_val_pos':[61, 57, 0, 0],
            'im_rms':  0.0501544137568,
            'im_sum': 0.00266315032371,
            'npts_real': 15876,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 0.00713019352406,
            'max_val_pos':[63, 63, 0, 0],
            'min_val': -0.00197783415206,
            'min_val_pos':[72, 58, 0, 0],
            'im_rms': 0.000865575412319,
            'im_sum': 0.130538275658,
            'npts_real': 15876}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
            'nchan': 1,
            'max_val': 4396210.5,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 4396210.5,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
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
            'start': 1.0784e+11,
            'end': 1.0784e+11,
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

        out, report8_a = th.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = th.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = th.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = th.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = th.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = th.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = th.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = th.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = th.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = th.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = th.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = th.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = th.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = th.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = th.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = th.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = th.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = th.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image.tt0', range_list=[-0.003, 0.035])
            self.png_creator(img+'.residual.tt0', range_list=[-0.003, 0.035])

            test_dict[testname]['images'] = \
                [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)

#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
            'nchan': 948,
            'max_val': 0.116102233529,
            'max_val_pos':[286, 233, 0, 568],
            'min_val': -0.0639033019543,
            'min_val_pos':[211, 186, 0, 592],
            'im_rms':  0.0106607721423,
            'rms_per_chan': [0.010544357743534185, 0.012164182832671495, 0.012074626408193847, 0.009608773051277117, 0.010081045610929583, 0.010510081502495537, 0.010527319732380157, 0.012514230024815057, 0.012753066817019756, 0.011681796537235542, 0.010530484090355572, 0.009501404486316086, 0.009490010520613795, 0.011376546145212042, 0.010273170623339859, 0.010443578473712106, 0.01014520940946224, 0.009211688240177103, 0.010328996596638062, 0.008857390893316648, 0.00969819905194494, 0.010097891889316984, 0.010715532081831087, 0.010566538223521662, 0.011550578029757013, 0.011835921152101766, 0.009523790740444496, 0.009610787542876992, 0.00889374579987403, 0.011338637905980761, 0.01160108504073822, 0.010099260357842827, 0.010038389403518005, 0.009466211294483254, 0.009647570939463142, 0.009353805624141558, 0.009747882683864631, 0.010480292978630698, 0.009020740823907015, 0.010563402579460455, 0.00935619375374926, 0.010208074753991778, 0.012220975043057088, 0.01114448147947678, 0.012098410856094012, 0.011752245125698849, 0.011295086034410555, 0.008843893127631184, 0.009017742131154766, 0.009294017517839195, 0.009439753643632231, 0.010298112798281545, 0.00983020859606145, 0.009457164872774293, 0.00912083668936835, 0.00939036093928175, 0.010078765634384814, 0.009624914878293646, 0.010992622185683143, 0.01216552068938102, 0.009417394615374051, 0.009573237173795658, 0.009836088837253487, 0.010043608815207505, 0.010220718217145255, 0.009479141768911016, 0.009419449299275788, 0.009420613091044862, 0.010629038754111685, 0.009649729908379494, 0.00859903466986099, 0.01035027268097479, 0.009671533858829157, 0.009735338494338334, 0.008349141951145837, 0.009526646655757486, 0.010130986474178332, 0.01273518138792606, 0.013225776820163676, 0.01208756113589589, 0.009598698796406378, 0.008617137137373982, 0.009399195723952727, 0.008975335184425445, 0.010881867424952852, 0.010133162330386693, 0.009899938735784194, 0.008983615752536919, 0.00949018251262457, 0.008751000010565207, 0.00911715549919336, 0.01001691482386163, 0.011819454830932457, 0.01109435639469309, 0.009557013726682832, 0.009709700220863355, 0.00985575414918664, 0.00998599406831754, 0.010091692782471956, 0.009572553930785084, 0.009344805379320887, 0.010135569827288124, 0.010056596352430012, 0.010832596842626574, 0.011605177042901175, 0.010443478320509977, 0.011356906939813928, 0.013042351451522865, 0.011139121213493433, 0.011820354361659095, 0.010041036239125105, 0.009512866701197682, 0.008612372264544444, 0.011002754164077793, 0.01053947338995949, 0.009506280721006897, 0.010188693032720478, 0.010354320149696072, 0.008838743955308245, 0.009906391825595074, 0.010290362507151294, 0.009280775710372087, 0.010735599452561117, 0.00835892430872177, 0.009964775347837304, 0.009121537140106377, 0.010162452525963437, 0.009692431760002418, 0.010050320410351348, 0.00907025900862872, 0.008130321688587665, 0.008942863823954227, 0.008146927578559456, 0.008819099589515166, 0.009665321994575399, 0.011361786727881892, 0.010731646979190946, 0.009549049512605124, 0.009248255030069584, 0.008275531984229586, 0.008997829132744339, 0.008966784777533397, 0.009388067977054564, 0.009972295171472296, 0.009735405654831794, 0.008485573682731, 0.009884558661097958, 0.01136664155742694, 0.008559108886390115, 0.009269856861011708, 0.008936470207152974, 0.009895498390789828, 0.009940648833550225, 0.009422436343932283, 0.01035819562040688, 0.01141470440757288, 0.010892989271589057, 0.010986451080061113, 0.01097053729839708, 0.009701094872393484, 0.01187355095514777, 0.012151980932532249, 0.009675848445614829, 0.00958223979015258, 0.008792361818022029, 0.008410285919927433, 0.009153095144803381, 0.008790964318579916, 0.008385975655907417, 0.00851976304347304, 0.009015694116101187, 0.00835006411696345, 0.008642582934538689, 0.00928710258838381, 0.01011067254224128, 0.011346564051039329, 0.010691953731919317, 0.009399566699808988, 0.009405531094879179, 0.010660394515107623, 0.010943244210507325, 0.00981095454240694, 0.009177624842325642, 0.010028030036168279, 0.011925784455447959, 0.009799696329063506, 0.008578421999458309, 0.008813715404176527, 0.009018413653366789, 0.009635099991548523, 0.009328979148105208, 0.008936577895979679, 0.00865068779495704, 0.011714015409683263, 0.009143298376525971, 0.0093722810186909, 0.010313452344660377, 0.008834954958669152, 0.010253318668779522, 0.009394638743027534, 0.010273425578342389, 0.011646966102790951, 0.012083569546386952, 0.012821615929802169, 0.01052987662527084, 0.010668361229993557, 0.008969315179641104, 0.009878504584890173, 0.009083355155260162, 0.009585038959164524, 0.01028527457111657, 0.011225333943569707, 0.009544898534207502, 0.010404987829551536, 0.009053075213135307, 0.009650786360307993, 0.010729754441548991, 0.011605132587502206, 0.009889033707096453, 0.009781390989644614, 0.011094611168044211, 0.012109054097268448, 0.010700499779381503, 0.01003792070939281, 0.011936605596990372, 0.011782773204662209, 0.00902961949098768, 0.010523704170008573, 0.009595105693410697, 0.010549837480142371, 0.010263110511920222, 0.009057636657348751, 0.01003624492957981, 0.011148296021869435, 0.009292331119195234, 0.008568711151010893, 0.009286797026264382, 0.01104598145381877, 0.010708854544863955, 0.010133324959161678, 0.012099744538605444, 0.012550124524717517, 0.008724364514235848, 0.008432250211062139, 0.009543766794389941, 0.009179955842670515, 0.009372585011702034, 0.009984965976952535, 0.008916048975671567, 0.009531643766691949, 0.00973942530515858, 0.01155893736638345, 0.009960265126137914, 0.008745760727547035, 0.009761776201659514, 0.01165928993115887, 0.010717212417944754, 0.011452037184819215, 0.010123049584187018, 0.00993441625515861, 0.009020593800774006, 0.009744067962277569, 0.009045716844558526, 0.010409172682813921, 0.010404044896126233, 0.009443969673974904, 0.01053770196088235, 0.011304803231166756, 0.00891148598331849, 0.0098978355019775, 0.009997024546928806, 0.009529657038869916, 0.009132380264567346, 0.010029448833730532, 0.010382881585698646, 0.009897208722102207, 0.009173659257114417, 0.009951231438726353, 0.01119932254374386, 0.011238295376514357, 0.01019615977334354, 0.009454245320281105, 0.009159666432461626, 0.009516838818912762, 0.008526334557564018, 0.009532637310721296, 0.0100845317703788, 0.010350639749445592, 0.009988597763990931, 0.010297622180831748, 0.01054102173760143, 0.010165435633922918, 0.008824701125215799, 0.00929833359003822, 0.009970561992630996, 0.010441506742680622, 0.008956932734484425, 0.0096514541102629, 0.009914471240254193, 0.010434250188652615, 0.01008036328298848, 0.010190819434328414, 0.009972681972539957, 0.010909137402546235, 0.013307048044779377, 0.01531369837138347, 0.01500782754653681, 0.010933753203318824, 0.00989177816589834, 0.009755611724437025, 0.011861803340913515, 0.009742804396682332, 0.009648183781161846, 0.011273772535601327, 0.011178953294584638, 0.010092910322705055, 0.009742215852868137, 0.009412600922294893, 0.011739646011811268, 0.011447443333342305, 0.010210912752285169, 0.010571829916096837, 0.012608339250828814, 0.013429207981686387, 0.010219208191188569, 0.0102535199365424, 0.011182581218080596, 0.012605602929318725, 0.012178550835134635, 0.010575806669458634, 0.01002954521194429, 0.009673264605422331, 0.010998378306806791, 0.012029893543896105, 0.010413268982625557, 0.010345477969553553, 0.010573977486438197, 0.012546954893696205, 0.011930234015695566, 0.010846313902970462, 0.011028267399350523, 0.010746512656062558, 0.01139699819869516, 0.010072069764373028, 0.009821142143793416, 0.011103098419875952, 0.012001135706067574, 0.011271820179113977, 0.01101392198086214, 0.010070536855658731, 0.010286201247195452, 0.010744364052413354, 0.012280220469572589, 0.011766010900169847, 0.011267844321872443, 0.010700172065778115, 0.00994925423698614, 0.010592208789229089, 0.0107102409232148, 0.010908709803928817, 0.010167156820648308, 0.009894447324525211, 0.010248123990740118, 0.01111849264679597, 0.012685849012986548, 0.011007547190594696, 0.00985806937224629, 0.010179497177305983, 0.01028098465015516, 0.011341065512633674, 0.010680276073247395, 0.009915179056407288, 0.0136551603473631, 0.013949462320439256, 0.010070230603335612, 0.010188485879220622, 0.010130233667421397, 0.009725768026021004, 0.010377466126136131, 0.009900726063704273, 0.009694481054639385, 0.009889371475791948, 0.010348629326478943, 0.011200587311736272, 0.010394943138605752, 0.011211652210645877, 0.010398454260753361, 0.010484057077021847, 0.011841168665806675, 0.009559723149813782, 0.010662785514774637, 0.011022018920913198, 0.010817023627469266, 0.009773502272303091, 0.010080396884424254, 0.01011651082652112, 0.009684927864918041, 0.009740614318067492, 0.010498062520519658, 0.010232104753925952, 0.012726303747665138, 0.009784364506650206, 0.010363682617365616, 0.010402214484748035, 0.011535827951877641, 0.010260452361780402, 0.009539657080016932, 0.009580809579446535, 0.013653094599638642, 0.013736790600000876, 0.011166365637421464, 0.01025133925493026, 0.011469928172074343, 0.01217997501951802, 0.010227602614939019, 0.011038370480602322, 0.011309919147343288, 0.01223123344997119, 0.010501539757993731, 0.011402965590321959, 0.010633220088281636, 0.01022400444841367, 0.009042635684091464, 0.010136158531394041, 0.010445046108694655, 0.010734748421322677, 0.010655949801969078, 0.010856986048873123, 0.01145987627223437, 0.01118386152130654, 0.010432606490126728, 0.010927140525502452, 0.011106983241414145, 0.009580214309081271, 0.011058833941938202, 0.010345010851262552, 0.011329825022648473, 0.011177725691948037, 0.012835355262928456, 0.011562555056677574, 0.010994165533931842, 0.011325552812055372, 0.012615387798011564, 0.013235599589657601, 0.012547479558920332, 0.01125155382236267, 0.011226333781524142, 0.009867099926252274, 0.010150976916901701, 0.00981456892640829, 0.010714013420052564, 0.012100414836000127, 0.013482305562984867, 0.012101608805487574, 0.011406733386446304, 0.011865020706158967, 0.010192404876503065, 0.010379782500127232, 0.010288816861602905, 0.011073875574562582, 0.011720944887620786, 0.011660526270432651, 0.011958079143296477, 0.010372963634161838, 0.010986196217407736, 0.010002087405069806, 0.010538518418827155, 0.009704941123596739, 0.01072720678974811, 0.010712087520487025, 0.011460224410312943, 0.011066250051779282, 0.011052989415309224, 0.011792111340277907, 0.010461854625871915, 0.010940461844210448, 0.010633699307646589, 0.010624022094149118, 0.011272242032483613, 0.010749869822808143, 0.011597062804568674, 0.01110807968738095, 0.010639725479597148, 0.010908395584711409, 0.013349933044667018, 0.0123243550483051, 0.01190271898274429, 0.012490664839553555, 0.012177466069438828, 0.011310780630467789, 0.01140940358525683, 0.010873455550563021, 0.010788849773109501, 0.011220843005138393, 0.011463926910727748, 0.010524467595075394, 0.010552208382105969, 0.011241882980229972, 0.011374143162653274, 0.012908776430391664, 0.01321474049935015, 0.011057229033692072, 0.011967901597247985, 0.010634871974387526, 0.012847393919979678, 0.011640242470566048, 0.011764354794187262, 0.011304207829268558, 0.01102205702407454, 0.011534653134776302, 0.011215772462419877, 0.01135268923219837, 0.013785467940882474, 0.01276952440553062, 0.010324791311942041, 0.011742973970237897, 0.01239586472863265, 0.011820481236013234, 0.011435139644675858, 0.01251611692344489, 0.01181323543372642, 0.010259393448015694, 0.011032029278022465, 0.01138934125996457, 0.011995068459853424, 0.01253808441029633, 0.010578386704919134, 0.01094522258091802, 0.012023277021817103, 0.011983116519653754, 0.010726480423429016, 0.010276658120803209, 0.010088385037851429, 0.010749839256216968, 0.010689062089672774, 0.010843887916758405, 0.011330782706137311, 0.011705587331553243, 0.01078000585040528, 0.010178986357696613, 0.011618165219610762, 0.011328541439778008, 0.010569799899248606, 0.01026219816096713, 0.010753014189736407, 0.010277625362334897, 0.011523218875299335, 0.012303865209747848, 0.013293662800221375, 0.012412719348663672, 0.011435203330629742, 0.01042380482199826, 0.011658856837661085, 0.011710373960274486, 0.012691089005450436, 0.013248739273267862, 0.012153365143317933, 0.010525964427932811, 0.01296496282370707, 0.011463993021550007, 0.011076028570385438, 0.012261069450156316, 0.012477345372605819, 0.01226426795633862, 0.011606690695066016, 0.012143643295226977, 0.011312259757459552, 0.010992315957294285, 0.01532943533106626, 0.011092024194712127, 0.011358704668808766, 0.010694287202125683, 0.010878015428811366, 0.011018691580814973, 0.01155699106864085, 0.011618045609910352, 0.010675703668991806, 0.011333729498651912, 0.011866470967968508, 0.010675876192637095, 0.010828976214518499, 0.011561959603500675, 0.012470936294951428, 0.01245237169975303, 0.011408547219405325, 0.013920173570889711, 0.011122172612362303, 0.011236706234314497, 0.011720912402096465, 0.01191176160235763, 0.013553146389120613, 0.015318038165449256, 0.012672532647506268, 0.011305885443388376, 0.012750520103345639, 0.012389453325464502, 0.013207008246510289, 0.011086890499330761, 0.01098436381363253, 0.011103446616269198, 0.011398117165826218, 0.011732621911866645, 0.013128667659727181, 0.011091147599992579, 0.010766761434001697, 0.01083985050764076, 0.012238486678992926, 0.01094365956136194, 0.010698246509702924, 0.01057389471760903, 0.011229440083161086, 0.011883817906017858, 0.01190745663644984, 0.013542478866497072, 0.011197611080423948, 0.01146965810623239, 0.011649422173178298, 0.010376195634184132, 0.013337478973473264, 0.01197910038576815, 0.012784823664355818, 0.012015337660711407, 0.012304555746398167, 0.012812880568715579, 0.011655837565537766, 0.010859847432719553, 0.010762540627138698, 0.011476999313503817, 0.012024669186394674, 0.010324907533403371, 0.011546371472360498, 0.012115046621070774, 0.01099964975371205, 0.010849575032497491, 0.01114553213835079, 0.009807210731195264, 0.011080608708596127, 0.01201945659989251, 0.011022773830302809, 0.011436605238807905, 0.010857167354394917, 0.010591362626392799, 0.0115255128187247, 0.011692115744482277, 0.011529151340701394, 0.014546531051519426, 0.012343550629432947, 0.012459213299429967, 0.011571824556046613, 0.011685125860296047, 0.01103394931287493, 0.010259184236371994, 0.013641962227648085, 0.010321747947849359, 0.011727616781904605, 0.011250622217955376, 0.01041199616493244, 0.010669919945036114, 0.01258526956270449, 0.011311246370266953, 0.010734114755387245, 0.011639490573190551, 0.011108436929631612, 0.013200806248911713, 0.013496042494565829, 0.01263957078900326, 0.011033512373956906, 0.011465822657641444, 0.010904556041386378, 0.013599431652120868, 0.010494701065814625, 0.009791714244609558, 0.008964002347221375, 0.010294884814437069, 0.01029810507291961, 0.013871652574021487, 0.011472848640171904, 0.010402558630102332, 0.010557678551486387, 0.009852673673925049, 0.012398389293026047, 0.011621507282032663, 0.009967140563868136, 0.01280606309090768, 0.009330114333904703, 0.011560788782769582, 0.01178744196876202, 0.010988048679879849, 0.010201127564532547, 0.01134249736787167, 0.0112425701365602, 0.009867586204808782, 0.011750882025482367, 0.012177731615191933, 0.010404844723897366, 0.010418253188478511, 0.010502939935957061, 0.010729987468609769, 0.011754797421181766, 0.011650378921476334, 0.011996508944877037, 0.010463821794374323, 0.011196107180678032, 0.011977020564712098, 0.011647565795514793, 0.011667368202315304, 0.010968175921477986, 0.011953650960032411, 0.010752291589534856, 0.013116454301674716, 0.013736783965398751, 0.012341451475925004, 0.012270776323448042, 0.010340757877337107, 0.010888267913550308, 0.011224134834619379, 0.009625939319606883, 0.009925185847972355, 0.009482847873369228, 0.010986612600040985, 0.011111298972473446, 0.010908854783352807, 0.010461443328959671, 0.010672619491429456, 0.010430146689624591, 0.01048379433666622, 0.011649538699094615, 0.011147832115376993, 0.012018384032169337, 0.01106759447869807, 0.009612710254814293, 0.011704904587770964, 0.010687630427532768, 0.00977516337390196, 0.009024562841962082, 0.009423225036468887, 0.010119454525294596, 0.010635528565256807, 0.010127006601197426, 0.009939977835643725, 0.009345031850815238, 0.012190187121965084, 0.010505037696602245, 0.010111623680303187, 0.009392067339766443, 0.009598288193694876, 0.010185604002435787, 0.012016813288971662, 0.011440126912197184, 0.01113106520286566, 0.013267667301903222, 0.010274530795273512, 0.00919477747332055, 0.008477021036017681, 0.009880232093632027, 0.010951146029138608, 0.010497290941526653, 0.0101716980211132, 0.009817587320741077, 0.011145982789163423, 0.010386744956456397, 0.010405913581197197, 0.011799570023060412, 0.009959634440018265, 0.010776640330465702, 0.014795886500233542, 0.010412884156670894, 0.00946176633444524, 0.010404398154873385, 0.010569348306726706, 0.01105674783544768, 0.011642019259874332, 0.010869528182031853, 0.00964747204271998, 0.010349378025902125, 0.01166077592649701, 0.009868153806603148, 0.011171345663126206, 0.009778759650023791, 0.00933177560738003, 0.008950225790739364, 0.012115801604401157, 0.01126339538430799, 0.011600988570399887, 0.010382586795258353, 0.009854484795689965, 0.00995869241937479, 0.010101066712418692, 0.009517649770286157, 0.00956091319387346, 0.00995987388676793, 0.01048802924636113, 0.009375408553236436, 0.010436849637018556, 0.009571260898819529, 0.010330123615380856, 0.010231236592148914, 0.009068717732484721, 0.012595809050203199, 0.011709465834549918, 0.010539977695652392, 0.009511595225627547, 0.009840731370430966, 0.009848102777870583, 0.009608981914888091, 0.010055363855187233, 0.010272034031323568, 0.009757275399578366, 0.00883997274911326, 0.01098936713047078, 0.010257393159455496, 0.008866838143674737, 0.009219709731706395, 0.009311583955426563, 0.008674466811696816, 0.0096540493127476, 0.010538075791420543, 0.010042582180017924, 0.009909704620295333, 0.010894836793617874, 0.010087624096615516, 0.009062855774670009, 0.00812474147026903, 0.009976151867812316, 0.010284316215532786, 0.009803365346022868, 0.009042753754190724, 0.008677203077829744, 0.010820179361645397, 0.011956743158821254, 0.011828637976101555, 0.010064940500071068, 0.012140840189519836, 0.011529829256948418, 0.009769740360379742, 0.009798044818356054, 0.011163036229953963, 0.010162791378669255, 0.009562504630899395, 0.01018206270241318, 0.01065164838099367, 0.008050497733057918, 0.010564249500390884, 0.00913150261764931, 0.01189821258736614, 0.010281986565400374, 0.008722411331043001, 0.011243438676865784, 0.01021610726013466, 0.008446091356968478, 0.008802457384278576, 0.008755990819858067, 0.00942145168545389, 0.0103427513442704, 0.009146632808879647, 0.010156712119463762, 0.009226573078603733, 0.008643679346235627, 0.009621798172188535, 0.00945240939004691, 0.007987468625400975, 0.00956493398918754, 0.010967775187150761, 0.010778074787163187, 0.009783134215189442, 0.00849674192273554, 0.008523265892321852, 0.010682347563223654, 0.01103103376005863, 0.012920569032402155, 0.012576980073156378, 0.011663880550671054, 0.010389746148365865, 0.010974010265535003, 0.01027064270445675, 0.00910778953628219, 0.009717706992309581, 0.009204976889737556, 0.008937822399494208, 0.009159958099142034, 0.010185995284263286, 0.01037940029857688, 0.010350445653442427, 0.011018941047741851, 0.009398759692134868, 0.011269752581613057, 0.012120975503486077, 0.010278812005534693, 0.009059875210870283, 0.008017090351726633, 0.008980427560123044, 0.00929248698563944, 0.010446507414092424, 0.008751499674850765, 0.009878823323884735, 0.008833954421548829, 0.009960604179491224, 0.011695771161034712, 0.010756479586952724, 0.008128773838808048, 0.0098897419700091, 0.007881070830252409, 0.00862996356881199, 0.009719473281425455, 0.009495391836243962, 0.009408180456787592, 0.008584712788125356, 0.00850299458181752, 0.008971239806717447, 0.009936199265547629, 0.009180681201659328, 0.009223184334341562, 0.011107303967276301, 0.009883751151564106, 0.007907321565336023, 0.0091262066720501, 0.009216752861686254, 0.009178905103799978, 0.009403308191776394, 0.01218332281191191, 0.011385888752999493, 0.009074419664166824, 0.009246354843707325, 0.0076941547167105094, 0.008952745121255394, 0.00925779245226104, 0.009884907644872749, 0.010255639476036957, 0.008037014749273852, 0.00941973765990461, 0.01002110887152415, 0.008457981115620032, 0.010194642285385696, 0.01281761748930397, 0.011057454468532681, 0.010187400989895463, 0.00863257461266235, 0.009555793520219453, 0.008919144122653173, 0.008304604711662937, 0.008522730427382655, 0.01013926741089452, 0.009828267951840651, 0.010160592946727233, 0.010312307637638811, 0.01043479509090218, 0.009167526409330535, 0.008890760078922063],
            'im_sum': 7294.44604131,
            'regn_sum': 28.4447036739,
            'npts_real': 191116800,
            'rms_per_field': 1,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [209, 387, 0, 0]), \
                      (img+'.image', False, [209, 388, 0, 0]), \
                      (img+'.image', True, [20, 204, 0, 0]), \
                      (img+'.image', False, [19, 204, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( \
            self.check_list_vals(im_stats_dict['rms_per_chan'], \
                exp_im_stats['rms_per_chan'], epsilon=0.01), True, \
            valname='RMS per channel of .image', exact=True)
        out, report1_u = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_w = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_x = th.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 191116800,
            'freq_bin': 244151.1796875,
            'start': 261764375854.0,
            'end': 261995831173.0,
            'nchan': 948,
            'mask_pix': 156228,
            'mask_regns': 39,
            'npts_real': 191116800}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
            'nchan': 948,
            'max_val': 1.0,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 0.200000017881,
            'min_val_pos':[56, 302, 0, 567],
            'npts_0.2': 105006529,
            'npts_0.5': 60643408,
            'npts_real': 191116800,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
            'nchan': 948,
            'max_val': 1.0,
            'max_val_pos':[240, 210, 0, 0],
            'min_val': -0.0456548333168,
            'min_val_pos':[230, 216, 0, 14],
            'im_rms':  0.0122779442959,
            'im_sum': 130.930138815,
            'npts_real': 191116800,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
            'nchan': 948,
            'max_val': 0.0589463338256,
            'max_val_pos':[269, 264, 0, 765],
            'min_val': -0.0639033019543,
            'min_val_pos':[211, 186, 0, 592],
            'im_rms': 0.0105758073519,
            'im_sum': 3583.31516361,
            'npts_real': 191116800}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
            'nchan': 948,
            'max_val': 45510.7695312,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 45195.4453125,
            'min_val_pos':[0, 0, 0, 594],
            'npts_real': 948}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
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
            'start': 261764375854.0,
            'end': 261995831173.0,
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

        out, report8_a = th.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = th.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = th.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = th.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = th.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = th.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = th.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = th.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = th.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = th.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = th.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = th.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = th.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = th.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = th.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = th.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = th.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = th.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[-0.01, 0.1])
            self.png_creator(img+'.residual', range_list=[-0.01, 0.1])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 2.06048989296,
            'max_val_pos':[291, 212, 0, 0],
            'min_val': -2.1858522892,
            'min_val_pos':[290, 152, 0, 0],
            'im_rms':  0.676557465791,
            'im_sum': 5498.32523989,
            'regn_sum': 8725.50744967,
            'npts_real': 201600,
            'rms_per_field': 1,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [211, 390, 0, 0]), \
                      (img+'.image', False, [211, 391, 0, 0]), \
                      (img+'.image', True, [18, 205, 0, 0]), \
                      (img+'.image', False, [17, 205, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_v = th.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)
        out, report1_w = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_a1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1

        # .mask report
        mask_stats_dict = self.image_stats(img, '.mask')

        exp_mask_stats = {'npts': 201600,
            'freq_bin': 16762504556.453735,
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'mask_pix': 0,
            'mask_regns': 0,
            'npts_real': 201600}

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[211, 203, 0, 0],
            'min_val': 0.200001135468,
            'min_val_pos':[81, 343, 0, 0],
            'npts_0.2': 113589,
            'npts_0.5': 64549,
            'npts_real': 201600,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 1.0,
            'max_val_pos':[240, 210, 0, 0],
            'min_val': -0.0525086708367,
            'min_val_pos':[230, 216, 0, 0],
            'im_rms':  0.0111846421981,
            'im_sum': 0.100949260701,
            'npts_real': 201600,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 2.06048989296,
            'max_val_pos':[291, 212, 0, 0],
            'min_val': -2.1858522892,
            'min_val_pos':[290, 152, 0, 0],
            'im_rms': 0.676557465791,
            'im_sum': 5498.32523989,
            'npts_real': 201600}

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
            'nchan': 1,
            'max_val': 30068706.0,
            'max_val_pos':[0, 0, 0, 0],
            'min_val': 30068706.0,
            'min_val_pos':[0, 0, 0, 0],
            'npts_real': 1}

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
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
            'start': 2.53574e+11,
            'end': 2.53574e+11,
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

        out, report8_a = th.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = th.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = th.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = th.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = th.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = th.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = th.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = th.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = th.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = th.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = th.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = th.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = th.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = th.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = th.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = th.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = th.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = th.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image', range_list=[-2.2, 2.1])
            self.png_creator(img+'.residual', range_list=[-2.2, 2.1])

            test_dict[testname]['images'] = \
                [img+'.image.moment8.png',img+'.residual.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


#-------------------------------------------------#
    @stats_dict(test_dict)
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
        self.copy_products(file_name+'0', file_name+'1')

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

        report0 = th.checkall( \
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
            'im_sum': 5410.97024339,
            'regn_sum': 8700.32470433,
            'npts_real': 201600,
            'rms_per_field': 1,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        report1_a = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [211, 390, 0, 0]), \
                      (img+'.image.tt0', False, [211, 391, 0, 0]), \
                      (img+'.image.tt0', True, [18, 205, 0, 0]), \
                      (img+'.image.tt0', False, [17, 205, 0, 0])])

        out, report1_b = th.check_val(im_stats_dict['com_bmaj'], \
            exp_im_stats['com_bmaj'], valname='Common beam major axis', \
            exact=False, epsilon=0.01)
        out, report1_c = th.check_val(im_stats_dict['com_bmin'], \
            exp_im_stats['com_bmin'], valname='Common beam minor axis', \
            exact=False, epsilon=0.01)
        out, report1_d = th.check_val(im_stats_dict['com_pa'], \
            exp_im_stats['com_pa'], valname='Common beam position angle', \
            exact=False, epsilon=0.01)
        out, report1_e = th.check_val(im_stats_dict['npts'], \
            exp_im_stats['npts'], valname='Number of pixels in .image', \
            exact=True)
        out, report1_f = th.check_val( \
            im_stats_dict['npts_unmasked'], exp_im_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .image', exact=True)
        out, report1_g = th.check_val(im_stats_dict['freq_bin'], \
            exp_im_stats['freq_bin'], valname='Frequency bin of .image', \
            exact=True)
        out, report1_h = th.check_val(im_stats_dict['start'], \
            exp_im_stats['start'], valname='Start channel of .image', \
            exact=True)
        out, report1_i = th.check_val(im_stats_dict['end'], \
            exp_im_stats['end'], valname='End channel of .image', exact=True)
        out, report1_j = th.check_val(im_stats_dict['nchan'], \
            exp_im_stats['nchan'], valname='Number of channels of .image', \
            exact=True)
        out, report1_k = th.check_val(im_stats_dict['max_val'], \
            exp_im_stats['max_val'], valname='Peak .image value', \
            exact=False, epsilon=0.01)
        out, report1_l = th.check_val( \
            im_stats_dict['max_val_pos'][0], exp_im_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .image', exact=True)
        out, report1_m = th.check_val( \
            im_stats_dict['max_val_pos'][1], exp_im_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .image', exact=True)
        out, report1_n = th.check_val( \
            im_stats_dict['max_val_pos'][3], exp_im_stats['max_val_pos'][3], \
            valname='Channel of peak value of .image', exact=True)
        out, report1_o = th.check_val(im_stats_dict['min_val'], \
            exp_im_stats['min_val'], valname='Minimum .image value', \
            exact=False, epsilon=0.01)
        out, report1_p = th.check_val( \
            im_stats_dict['min_val_pos'][0], exp_im_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .image', \
            exact=True)
        out, report1_q = th.check_val( \
            im_stats_dict['min_val_pos'][1], exp_im_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .image', \
            exact=True)
        out, report1_r = th.check_val( \
            im_stats_dict['min_val_pos'][3], exp_im_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .image', exact=True)
        out, report1_s = th.check_val(im_stats_dict['im_rms'], \
            exp_im_stats['im_rms'], valname='RMS of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_t = th.check_val( im_stats_dict['im_sum'], \
            exp_im_stats['im_sum'], valname='Sum of the whole .image', \
            exact=False, epsilon=0.01)
        out, report1_u = th.check_val(im_stats_dict['regn_sum'], \
            exp_im_stats['regn_sum'], valname='Sum of a .image region', \
            exact=False, epsilon=0.01)
        out, report1_v = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)
        out, report1_w = th.check_val( \
            im_stats_dict['rms_per_field'], exp_im_stats['rms_per_field'], \
            valname='RMS per field of .image', exact=False, \
            epsilon=0.01)
        out, report1_x = th.check_val(im_stats_dict['im_fit'][0][0], \
            exp_im_stats['im_fit'][0][0], valname='Fit center x coord', \
            exact=False, epsilon=0.01)
        out, report1_y = th.check_val(im_stats_dict['im_fit'][0][1], \
            exp_im_stats['im_fit'][0][1], valname='Fit center y coord', \
            exact=False, epsilon=0.01)
        out, report1_z = th.check_val( \
            im_stats_dict['im_fit'][1], exp_im_stats['im_fit'][1], \
            valname='Fit channel location', exact=True)
        out, report1_a1 = th.check_val( \
            im_stats_dict['im_fit'][2], exp_im_stats['im_fit'][2], \
            valname='Frequency of fit', exact=True)
        out, report1_b1 = th.check_val(im_stats_dict['im_fit'][3], \
            exp_im_stats['im_fit'][3], valname='Peak of the fit', \
            exact=False, epsilon=0.01)
        out, report1_c1 = th.check_val(im_stats_dict['im_fit'][4], \
            exp_im_stats['im_fit'][4], valname='Major axis of fit', \
            exact=False, epsilon=0.01)
        out, report1_d1 = th.check_val(im_stats_dict['im_fit'][5], \
            exp_im_stats['im_fit'][5], valname='Minor axis of fit', \
            exact=False, epsilon=0.01)

        report1 = report1_a + report1_b + report1_c + report1_d + \
            report1_e + report1_f + report1_g + report1_h + report1_i + \
            report1_j + report1_k + report1_l + report1_m + report1_n + \
            report1_o + report1_p + report1_q + report1_r + report1_s + \
            report1_t + report1_u + report1_v + report1_w + report1_x + \
            report1_y + report1_z + report1_a1 + report1_b1 + report1_c1 + \
            report1_d1

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

        out, report2_a = th.check_val(mask_stats_dict['npts'], \
            exp_mask_stats['npts'], valname='Number of pixels in .mask', \
            exact=True)
        out, report2_b = th.check_val( \
            mask_stats_dict['freq_bin'], exp_mask_stats['freq_bin'], \
            valname='Frequency bin of .mask', exact=True)
        out, report2_c = th.check_val(mask_stats_dict['start'], \
            exp_mask_stats['start'], valname='Start channel of .mask', \
            exact=True)
        out, report2_d = th.check_val(mask_stats_dict['end'], \
            exp_mask_stats['end'], valname='End channel of .mask', exact=True)
        out, report2_e = th.check_val(mask_stats_dict['nchan'], \
            exp_mask_stats['nchan'], valname='Number of channels in .mask', \
            exact=True)
        out, report2_f = th.check_val( \
            mask_stats_dict['mask_pix'], exp_mask_stats['mask_pix'], \
            valname='Number of pixels masked', exact=True)
        out, report2_g = th.check_val( \
            mask_stats_dict['mask_regns'], exp_mask_stats['mask_regns'], \
            valname='Number of regions in .mask', exact=True)
        out, report2_h = th.check_val( \
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
            'npts_real': 201600,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report3_a = th.check_val(pb_stats_dict['npts'], \
            exp_pb_stats['npts'], valname='Number of pixels in .pb', \
            exact=True)
        out, report3_b = th.check_val( \
            pb_stats_dict['npts_unmasked'], exp_pb_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .pb', exact=True)
        out, report3_c = th.check_val(pb_stats_dict['freq_bin'], \
            exp_pb_stats['freq_bin'], valname='Frequency bin of .pb', \
            exact=True)
        out, report3_d = th.check_val(pb_stats_dict['start'], \
            exp_pb_stats['start'], valname='Start channel of .pb', exact=True)
        out, report3_e = th.check_val(pb_stats_dict['end'], \
            exp_pb_stats['end'], valname='End channel of .pb', exact=True)
        out, report3_f = th.check_val(pb_stats_dict['nchan'], \
            exp_pb_stats['nchan'], valname='Number of channels of .pb', \
            exact=True)
        out, report3_g = th.check_val(pb_stats_dict['max_val'], \
            exp_pb_stats['max_val'], valname='Maximum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_h = th.check_val( \
            pb_stats_dict['max_val_pos'][0], exp_pb_stats['max_val_pos'][0], \
            valname='RA pixel location of peak value of .pb', exact=True)
        out, report3_i = th.check_val( \
            pb_stats_dict['max_val_pos'][1], exp_pb_stats['max_val_pos'][1], \
            valname='Dec pixel location of peak value of .pb', exact=True)
        out, report3_j = th.check_val( \
            pb_stats_dict['max_val_pos'][3], exp_pb_stats['max_val_pos'][3], \
            valname='Channel of peak value of .pb', exact=True)
        out, report3_k = th.check_val(pb_stats_dict['min_val'], \
            exp_pb_stats['min_val'], valname='Minimum .pb value', \
            exact=False, epsilon=0.01)
        out, report3_l = th.check_val( \
            pb_stats_dict['min_val_pos'][0], exp_pb_stats['min_val_pos'][0], \
            valname='RA pixel location of minimum value of .pb', exact=True)
        out, report3_m = th.check_val( \
            pb_stats_dict['min_val_pos'][1], exp_pb_stats['min_val_pos'][1], \
            valname='Dec pixel location of minimum value of .pb', exact=True)
        out, report3_n = th.check_val( \
            pb_stats_dict['min_val_pos'][3], exp_pb_stats['min_val_pos'][3], \
            valname='Channel of minimum value of .pb', exact=True)
        out, report3_o = th.check_val(pb_stats_dict['npts_0.2'], \
            exp_pb_stats['npts_0.2'], valname='Number of points above .pb '
            '0.2', exact=False, epsilon=0.01)
        out, report3_p = th.check_val(pb_stats_dict['npts_0.5'], \
            exp_pb_stats['npts_0.5'], valname='Number of points above .pb '
            '0.5', exact=False, epsilon=0.01)
        out, report3_q = th.check_val( \
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
            'npts_real': 201600,
            'im_fit': [[38.303527208322556, 37.46382702762046],
                        254, 220.31469461816917, 0.909677723621573,
                        13.925317200593012, 7.225723641827008]}

        out, report4_a = th.check_val(psf_stats_dict['npts'], \
            exp_psf_stats['npts'], valname='Number of pixels in .psf', \
            exact=True)
        out, report4_b = th.check_val( \
            psf_stats_dict['npts_unmasked'], exp_psf_stats['npts_unmasked'], \
            valname='Number of unmasked pixels in .psf', exact=True)
        out, report4_c = th.check_val( \
            psf_stats_dict['freq_bin'], exp_psf_stats['freq_bin'], \
            valname='Frequency bin of .psf', exact=True)
        out, report4_d = th.check_val(psf_stats_dict['start'], \
            exp_psf_stats['start'], valname='Start channel of .psf', \
            exact=True)
        out, report4_e = th.check_val(psf_stats_dict['end'], \
            exp_psf_stats['end'], valname='End channel of .psf', exact=True)
        out, report4_f = th.check_val(psf_stats_dict['nchan'], \
            exp_psf_stats['nchan'], valname='Number of channels of .psf', \
            exact=True)
        out, report4_g = th.check_val(psf_stats_dict['max_val'], \
            exp_psf_stats['max_val'], valname='Maximum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_h = th.check_val( \
            psf_stats_dict['max_val_pos'][0], \
            exp_psf_stats['max_val_pos'][0], valname='RA pixel location of '
            'peak value of .psf', exact=True)
        out, report4_i = th.check_val( \
            psf_stats_dict['max_val_pos'][1], \
            exp_psf_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .psf', exact=True)
        out, report4_j = th.check_val( \
            psf_stats_dict['max_val_pos'][3], \
            exp_psf_stats['max_val_pos'][3], valname='Channel of peak value'
            ' of .psf', exact=True)
        out, report4_k = th.check_val(psf_stats_dict['min_val'], \
            exp_psf_stats['min_val'], valname='Minimum .psf value', \
            exact=False, epsilon=0.01)
        out, report4_l = th.check_val( \
            psf_stats_dict['min_val_pos'][0], \
            exp_psf_stats['min_val_pos'][0], valname='RA pixel location of '
            'minimum value of .psf', exact=True)
        out, report4_m = th.check_val( \
            psf_stats_dict['min_val_pos'][1], \
            exp_psf_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .psf', exact=True)
        out, report4_n = th.check_val( \
            psf_stats_dict['min_val_pos'][3], \
            exp_psf_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .psf', exact=True)
        out, report4_o = th.check_val(psf_stats_dict['im_rms'], \
            exp_psf_stats['im_rms'], valname='RMS of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_p = th.check_val(psf_stats_dict['im_sum'], \
            exp_psf_stats['im_sum'], valname='Sum of the whole .psf', \
            exact=False, epsilon=0.01)
        out, report4_q = th.check_val( \
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

        out, report5_a = th.check_val(resid_stats_dict['npts'], \
            exp_resid_stats['npts'], valname='Number of pixels in '
            '.residual', exact=True)
        out, report5_b = th.check_val( \
            resid_stats_dict['npts_unmasked'], \
            exp_resid_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .residual', exact=True)
        out, report5_c = th.check_val( \
            resid_stats_dict['freq_bin'], exp_resid_stats['freq_bin'], \
            valname='Frequency bin of .residual', exact=True)
        out, report5_d = th.check_val(resid_stats_dict['start'], \
            exp_resid_stats['start'], valname='Start channel of .residual', \
            exact=True)
        out, report5_e = th.check_val(resid_stats_dict['end'], \
            exp_resid_stats['end'], valname='End channel of .residual', \
            exact=True)
        out, report5_f = th.check_val(resid_stats_dict['nchan'], \
            exp_resid_stats['nchan'], valname='Number of channels of '
            '.residual', exact=True)
        out, report5_g = th.check_val( \
            resid_stats_dict['max_val'], exp_resid_stats['max_val'], \
            valname='Maximum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_h = th.check_val( \
            resid_stats_dict['max_val_pos'][0], \
            exp_resid_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_i = th.check_val( \
            resid_stats_dict['max_val_pos'][1], \
            exp_resid_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .residual', exact=True)
        out, report5_j = th.check_val( \
            resid_stats_dict['max_val_pos'][3], \
            exp_resid_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .residual', exact=True)
        out, report5_k = th.check_val( \
            resid_stats_dict['min_val'], exp_resid_stats['min_val'], \
            valname='Minimum .residual value', exact=False, \
            epsilon=0.01)
        out, report5_l = th.check_val( \
            resid_stats_dict['min_val_pos'][0], \
            exp_resid_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_m = th.check_val( \
            resid_stats_dict['min_val_pos'][1], \
            exp_resid_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .residual', exact=True)
        out, report5_n = th.check_val( \
            resid_stats_dict['min_val_pos'][3], \
            exp_resid_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .residual', exact=True)
        out, report5_o = th.check_val( \
            resid_stats_dict['im_rms'], exp_resid_stats['im_rms'], \
            valname='RMS of the whole .residual', exact=False, epsilon=0.01)
        out, report5_p = th.check_val( \
            resid_stats_dict['im_sum'], exp_resid_stats['im_sum'], \
            valname='Sum of the whole .residual', exact=False, epsilon=0.01)
        out, report5_q = th.check_val( \
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

        out, report6_a = th.check_val(model_stats_dict['npts'], \
            exp_model_stats['npts'], valname='Number of pixels in .model', \
            exact=True)
        out, report6_b = th.check_val( \
            model_stats_dict['npts_unmasked'], \
            exp_model_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .model', exact=True)
        out, report6_c = th.check_val( \
            model_stats_dict['freq_bin'], exp_model_stats['freq_bin'], \
            valname='Frequency bin of .model', exact=True)
        out, report6_d = th.check_val(model_stats_dict['start'], \
            exp_model_stats['start'], valname='Start channel of .model', \
            exact=True)
        out, report6_e = th.check_val(model_stats_dict['end'], \
            exp_model_stats['end'], valname='End channel of .model', \
            exact=True)
        out, report6_f = th.check_val(model_stats_dict['nchan'], \
            exp_model_stats['nchan'], valname='Number of channels of '
            '.model', exact=True)
        out, report6_g = th.check_val( \
            model_stats_dict['max_val'], exp_model_stats['max_val'], \
            valname='Maximum .model value', exact=False, epsilon=0.01)
        out, report6_h = th.check_val( \
            model_stats_dict['max_val_pos'][0], \
            exp_model_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .model', exact=True)
        out, report6_i = th.check_val( \
            model_stats_dict['max_val_pos'][1], \
            exp_model_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .model', exact=True)
        out, report6_j = th.check_val( \
            model_stats_dict['max_val_pos'][3], \
            exp_model_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .model', exact=True)
        out, report6_k = th.check_val( \
            model_stats_dict['min_val'], exp_model_stats['min_val'], \
            valname='Minimum .model value', exact=False, epsilon=0.01)
        out, report6_l = th.check_val( \
            model_stats_dict['min_val_pos'][0], \
            exp_model_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_m = th.check_val( \
            model_stats_dict['min_val_pos'][1], \
            exp_model_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .model', exact=True)
        out, report6_n = th.check_val( \
            model_stats_dict['min_val_pos'][3], \
            exp_model_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .model', exact=True)
        out, report6_o = th.check_val( \
            model_stats_dict['im_rms'], exp_model_stats['im_rms'], \
            valname='RMS of the whole .model', exact=False, epsilon=0.01)
        out, report6_p = th.check_val( \
            model_stats_dict['im_sum'], exp_model_stats['im_sum'], \
            valname='Sum of the whole .model', exact=False, epsilon=0.01)
        out, report6_q = th.check_val( \
            model_stats_dict['regn_sum'], exp_model_stats['regn_sum'], \
            valname='Sum of a region of .model', exact=False, epsilon=0.01)
        out, report6_r = th.check_val( \
            model_stats_dict['mask_non0'], \
            exp_model_stats['mask_non0'], valname='Non zero values in masked'
            ' regions of .model', exact=True)
        out, report6_s = th.check_val( \
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

        out, report7_a = th.check_val(sumwt_stats_dict['npts'], \
            exp_sumwt_stats['npts'], valname='Number of pixels in .sumwt', \
            exact=True)
        out, report7_b = th.check_val( \
            sumwt_stats_dict['npts_unmasked'], \
            exp_sumwt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .sumwt', exact=True)
        out, report7_c = th.check_val( \
            sumwt_stats_dict['freq_bin'], exp_sumwt_stats['freq_bin'], \
            valname='Frequency bin of .sumwt', exact=True)
        out, report7_d = th.check_val(sumwt_stats_dict['start'], \
            exp_sumwt_stats['start'], valname='Start channel of .sumwt', \
            exact=True)
        out, report7_e = th.check_val(sumwt_stats_dict['end'], \
            exp_sumwt_stats['end'], valname='End channel of .sumwt', \
            exact=True)
        out, report7_f = th.check_val(sumwt_stats_dict['nchan'], \
            exp_sumwt_stats['nchan'], valname='Number of channels of '
            '.sumwt', exact=True)
        out, report7_g = th.check_val( \
            sumwt_stats_dict['max_val'], exp_sumwt_stats['max_val'], \
            valname='Maximum .sumwt value', exact=False, epsilon=0.01)
        out, report7_h = th.check_val( \
            sumwt_stats_dict['max_val_pos'][0], \
            exp_sumwt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_i = th.check_val( \
            sumwt_stats_dict['max_val_pos'][1], \
            exp_sumwt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .sumwt', exact=True)
        out, report7_j = th.check_val( \
            sumwt_stats_dict['max_val_pos'][3], \
            exp_sumwt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .sumwt', exact=True)
        out, report7_k = th.check_val( \
            sumwt_stats_dict['min_val'], exp_sumwt_stats['min_val'], \
            valname='Minimum .sumwt value', exact=False, epsilon=0.01)
        out, report7_l = th.check_val( \
            sumwt_stats_dict['min_val_pos'][0], \
            exp_sumwt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_m = th.check_val( \
            sumwt_stats_dict['min_val_pos'][1], \
            exp_sumwt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .sumwt', exact=True)
        out, report7_n = th.check_val( \
            sumwt_stats_dict['min_val_pos'][3], \
            exp_sumwt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .sumwt', exact=True)
        out, report7_o = th.check_val( \
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

        out, report8_a = th.check_val(wt_stats_dict['npts'], \
            exp_wt_stats['npts'], valname='Number of pixels in .weight', \
            exact=True)
        out, report8_b = th.check_val( \
            wt_stats_dict['npts_unmasked'], \
            exp_wt_stats['npts_unmasked'], valname='Number of unmasked '
            'pixels in .weight', exact=True)
        out, report8_c = th.check_val( \
            wt_stats_dict['freq_bin'], exp_wt_stats['freq_bin'], \
            valname='Frequency bin of .weight', exact=True)
        out, report8_d = th.check_val(wt_stats_dict['start'], \
            exp_wt_stats['start'], valname='Start channel of .weight', \
            exact=True)
        out, report8_e = th.check_val(wt_stats_dict['end'], \
            exp_wt_stats['end'], valname='End channel of .weight', \
            exact=True)
        out, report8_f = th.check_val(wt_stats_dict['nchan'], \
            exp_wt_stats['nchan'], valname='Number of channels of '
            '.weight', exact=True)
        out, report8_g = th.check_val( \
            wt_stats_dict['max_val'], exp_wt_stats['max_val'], \
            valname='Maximum .weight value', exact=False, epsilon=0.01)
        out, report8_h = th.check_val( \
            wt_stats_dict['max_val_pos'][0], \
            exp_wt_stats['max_val_pos'][0], valname='RA pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_i = th.check_val( \
            wt_stats_dict['max_val_pos'][1], \
            exp_wt_stats['max_val_pos'][1], valname='Dec pixel location of'
            ' peak value of .weight', exact=True)
        out, report8_j = th.check_val( \
            wt_stats_dict['max_val_pos'][3], \
            exp_wt_stats['max_val_pos'][3], valname='Channel of peak '
            'value of .weight', exact=True)
        out, report8_k = th.check_val( \
            wt_stats_dict['min_val'], exp_wt_stats['min_val'], \
            valname='Minimum .weight value', exact=False, epsilon=0.01)
        out, report8_l = th.check_val( \
            wt_stats_dict['min_val_pos'][0], \
            exp_wt_stats['min_val_pos'][0], valname='RA pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_m = th.check_val( \
            wt_stats_dict['min_val_pos'][1], \
            exp_wt_stats['min_val_pos'][1], valname='Dec pixel location of'
            ' minimum value of .weight', exact=True)
        out, report8_n = th.check_val( \
            wt_stats_dict['min_val_pos'][3], \
            exp_wt_stats['min_val_pos'][3], valname='Channel of minimum '
            'value of .weight', exact=True)
        out, report8_o = th.check_val( \
            wt_stats_dict['im_rms'], exp_wt_stats['im_rms'], \
            valname='RMS of the whole .weight', exact=False, epsilon=0.01)
        out, report8_p = th.check_val( \
            wt_stats_dict['im_sum'], exp_wt_stats['im_sum'], \
            valname='Sum of the whole .weight', exact=False, epsilon=0.01)
        out, report8_q = th.check_val(wt_stats_dict['npts_0.2'], \
            exp_wt_stats['npts_0.2'], valname='Number of points above .wt '
            '0.2', exact=False, epsilon=0.01)
        out, report8_r = th.check_val(wt_stats_dict['npts_0.3'], \
            exp_wt_stats['npts_0.3'], valname='Number of points above .wt '
            '0.3', exact=False, epsilon=0.01)
        out, report8_s = th.check_val( \
            im_stats_dict['npts_real'], exp_im_stats['npts_real'], \
            valname='Number of real pixels in .image', exact=True)

        report8 = report8_a + report8_b + report8_c + report8_d + \
            report8_e + report8_f + report8_g + report8_h + report8_i + \
            report8_j + report8_k + report8_l + report8_m + report8_n + \
            report8_o + report8_p + report8_q + report8_r + report8_s

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output=test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report

        if not CASA6:
            self.png_creator(img+'.image.tt0', range_list=[-2.2, 2.1])
            self.png_creator(img+'.residual.tt0', range_list=[-2.2, 2.1])

            test_dict[testname]['images'] = \
                [img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png']

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


def suite():
     return [Test_standard, Test_mosaic]

# Main #
if __name__ == '__main__':
    unittest.main()



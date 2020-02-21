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
import matplotlib.pyplot as pyplot

from casatestutils.imagerhelpers import TestHelpers
th = TestHelpers()
from casatestutils import generate_weblog
from casatestutils import add_to_dict
from casatestutils import stats_dict


CASA6 = False
try:
    from casatools import ctsys, quanta, measures, image, vpmanager, calibrater
    from casatasks import casalog, delmod, imsubimage, tclean, uvsub, imhead, imsmooth, immath, widebandpbcor, immoments#, imview
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    CASA6 = True
    _ia = image()

except ImportError:
    from __main__ import default  # reset given task to its default values
    from tasks import *  # Imports all casa tasks
    from taskinit import *  # Imports all casa tools
    from parallel.parallel_task_helper import ParallelTaskHelper

    _ia = iatool()

# location of data
data_path = '/lustre/naasc/sciops/comm/sbooth/CASA_ALMA_pipeline/data_dir/'
#data_path = os.environ.get('CASAPATH').split()[0] + '/casa-data-vt/vlass/'


## Base Test class with Utility functions
class test_tclean_base(unittest.TestCase):

    def setUp(self):
        self._myia = _ia
        self.epsilon = 0.01 # sets epsilon as a percentage (1%)
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
        if(pstr.count("(Fail") > 0 ):
             self.fail("\n"+pstr)

    def check_list_vals(self, list1, list2, test, epsilon=None):
        """ compares 2 lists and returns if they are equivalent (within error) 
        """
        report = ''
        if len(list1) == len(list2) and epsilon is None:
            i = 0
            while i < len(list1):
                result, pstr = th.check_val(list1[i], list2[i], \
                    valname=test+' index '+str(i), exact=True)
                if result == False:
                    report = pstr
                    break
                i += 1
        elif len(list1) == len(list2) and epsilon is not None:
            i = 0
            while i < len(list1):
                result, pstr = th.check_val(list1[i], list2[i], \
                    valname=test+' index '+str(i), exact=False, epsilon=epsilon)
                if result == False:
                    report = pstr
                    break
                i += 1
        else:
            result = False

        return result, report

    def getNameDoc(self):
        testname = inspect.stack()[1][3]
        print("Test name  : " + testname)
        doc = eval('self.'+testname + '.__doc__')
        print("Doc : " +  doc)

        return testname, doc

    def copy_products(self, old_pname, new_pname, ignore=None):
        """ function to copy iter0 images to iter1 images
            (taken from pipeline)
        """
        imlist = glob.glob('%s.*' % (old_pname))
        imlist = [xx for xx in imlist if ignore is None or ignore not in xx]
        for image_name in imlist:
            newname = image_name.replace(old_pname, new_pname)
            if image_name == old_pname + '.workdirectory':
                mkcmd = 'mkdir '+ newname
                os.system(mkcmd)
                self.copy_products(os.path.join(image_name, old_pname), \
                    os.path.join(newname, new_pname))
            else:
                shutil.copytree(image_name, newname, symlinks=True)

    def cube_beam_stats(self, img):
        """ function to return per-channel beam statistics
            will be deprecated and combined into image_stats 
            once CASA beam issue is fixed
        """
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

    def image_stats(self, img, region_file=None, field_regions=[]):
        """ function that takes an image file and returns a statistics
            dictionary
        """
        self._myia.open(img)
        stats_dict = {}

        statistics = self._myia.statistics()
        chunk = self._myia.getchunk()

        # stats returned for all images
        im_size = self._myia.boundingbox()['imageShape'].tolist()
        stats_dict['npts'] = im_size[0]*im_size[1]*im_size[3]
        stats_dict['npts_unmasked'] = statistics['npts'][0]
        stats_dict['npts_real'] = numpy.count_nonzero(~numpy.isnan(chunk))
        stats_dict['freq_bin'] = self._myia.summary()['incr'][3]
        stats_dict['start'] = float( \
            statistics['blcf'].split(', ')[3].split('Hz')[0])
        stats_dict['end'] = float( \
            statistics['trcf'].split(', ')[3].split('Hz')[0])
        stats_dict['start_delta'] = stats_dict['start']
        stats_dict['end_delta'] = stats_dict['end']
        stats_dict['nchan'] = im_size[3]

        # stats returned for all images except .mask
        if not img.endswith('.mask'):
            stats_dict['max_val'] = statistics['max'][0]
            stats_dict['max_val_pos'] = statistics['maxpos'].tolist()
            max_loc = [stats_dict['max_val_pos'][0], \
                stats_dict['max_val_pos'][1]]
            stats_dict['min_val'] = statistics['min'][0]
            stats_dict['min_val_pos'] = statistics['minpos'].tolist()
            stats_dict['im_rms'] = statistics['rms'][0]

        # stats returned if a region file is given
        if region_file != None:
            region_stats = self._myia.statistics(region=region_file)
            stats_dict['regn_sum'] = region_stats['sum'][0]

        # stats returned for .image(.tt0)
        if img.endswith('.image') or img.endswith('.image.tt0'):
            commonbeam = self._myia.commonbeam()
            stats_dict['com_bmin'] = commonbeam['minor']['value']
            stats_dict['com_bmaj'] = commonbeam['major']['value']
            stats_dict['com_pa'] = commonbeam['pa']['value']
            if 'mfs_eph' not in img:
                try:
                    fit_dict = self._myia.fitcomponents( \
                        region=region_file)['results']['component0']
                    stats_dict['im_fit'] = [fit_dict['peak']['value'], \
                        fit_dict['shape']['majoraxis']['value'], \
                        fit_dict['shape']['minoraxis']['value']]
                    stats_dict['im_fit_loc'] = [ \
                        fit_dict['spectrum']['channel'], \
                        fit_dict['spectrum']['frequency']['m0']['value']]
                    stats_dict['im_fit_pix'] = fit_dict['pixelcoords'].tolist()
                except KeyError:
                    stats_dict['im_fit'] = [1.0, 1.0, 1.0]
                    stats_dict['im_fit_loc'] = [1.0, 1.0]
                    stats_dict['im_fit_pix'] = [1.0, 1.0]
            if 'cube' in img:
                stats_dict['rms_per_chan'] = \
                    self._myia.statistics(axes=[0,1])['rms'].tolist()
                stats_dict['profile'] = self.cube_profile_fit( \
                    img, max_loc, stats_dict['nchan'])
            if 'mosaic' in img:
                stats_dict['rms_per_field'] = []
                for region in field_regions:
                    stats_dict['rms_per_field'].append( \
                        self._myia.statistics(region=region)['rms'][0])

        # stats returned if not .pb(.tt0), .sumwt(.tt0), or .mask
        if ('.pb' not in img and '.sumwt' not in img and not
                img.endswith('.mask')):
            stats_dict['im_sum'] = statistics['sum'][0]

        if img.endswith('.mask'):
            stats_dict['mask_pix'] = numpy.count_nonzero(chunk)
            stats_dict['mask_regns'] = scipy.ndimage.label(chunk)[1]

        if img.endswith('.pb') or img.endswith('.pb.tt0'):
            stats_dict['npts_0.2'] = numpy.count_nonzero(chunk>0.2)
            stats_dict['npts_0.5'] = numpy.count_nonzero(chunk>0.5)
            try:
                fit_dict = self._myia.fitcomponents( \
                    region=region_file)['results']['component0']
                stats_dict['pb_fit'] = [fit_dict['peak']['value'], \
                    fit_dict['shape']['majoraxis']['value'], \
                    fit_dict['shape']['minoraxis']['value']]
                stats_dict['pb_fit_loc'] = [ \
                    fit_dict['spectrum']['channel'], \
                    fit_dict['spectrum']['frequency']['m0']['value']]
                stats_dict['pb_fit_pix'] = fit_dict['pixelcoords'].tolist()
            except KeyError:
                stats_dict['pb_fit'] = [1.0, 1.0, 1.0]
                stats_dict['pb_fit_loc'] = [1.0, 1.0]
                stats_dict['pb_fit_pix'] = [1.0, 1.0]

        if img.endswith('.psf') or img.endswith('.psf.tt0'):
            try:
                fit_dict = self._myia.fitcomponents( \
                    region=region_file)['results']['component0']
                stats_dict['psf_fit'] = [fit_dict['peak']['value'], \
                    fit_dict['shape']['majoraxis']['value'], \
                    fit_dict['shape']['minoraxis']['value']]
                stats_dict['psf_fit_loc'] = [ \
                    fit_dict['spectrum']['channel'], \
                    fit_dict['spectrum']['frequency']['m0']['value']]
                stats_dict['psf_fit_pix'] = \
                    fit_dict['pixelcoords'].tolist()
            except KeyError:
                stats_dict['psf_fit'] = [1.0, 1.0, 1.0]
                stats_dict['psf_fit_loc'] = [1.0, 1.0]
                stats_dict['psf_fit_pix'] = [1.0, 1.0]

        if (img.endswith('.model') or img.endswith('.model.tt0') or
                img.endswith('.alpha')):
            stats_dict['mask_non0'] = numpy.count_nonzero(chunk)

        if img.endswith('.weight') or img.endswith('.weight.tt0'):
            stats_dict['npts_0.2'] = numpy.count_nonzero(chunk>0.2)
            stats_dict['npts_0.3'] = numpy.count_nonzero(chunk>0.3)

        self._myia.close()

        return stats_dict

    def image_list(self, img, mode):
        """ function used to return expected imaging output files """
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

    def stats_compare(self, exp_dict, stats_dict, suffix):
        """ function to compare expected dictionary with returned
            dictionary
        """
        report = ''
        for key in exp_dict:
            if type(exp_dict[key][1]) == list:
                if exp_dict[key][0] == True:
                    result, pstr = self.check_list_vals(stats_dict[key], 
                        exp_dict[key][1], test=suffix+' '+key)
                    report += th.check_val(result, True, \
                        valname=suffix+' '+key, exact=True)[1]
                    report += pstr
                else:
                    result, pstr = self.check_list_vals(stats_dict[key], \
                        exp_dict[key][1], test=suffix+' '+key, \
                        epsilon=self.epsilon)
                    report += th.check_val(result, True, \
                        valname=suffix+' '+key, exact=True)[1]
                    report += pstr
            else:
                if exp_dict[key][0] == True:
                    report += th.check_val(stats_dict[key], \
                        exp_dict[key][1], valname=suffix+' '+key, exact=True)[1]
                elif exp_dict[key][0] == False and exp_dict[key][1] == 0.0:
                    report += th.check_val(stats_dict[key], \
                        exp_dict[key][1], valname=suffix+' '+key, exact=True)[1]
                else:
                    report += th.check_val(stats_dict[key], \
                        exp_dict[key][1], valname=suffix+' '+key, exact=False, \
                        epsilon=self.epsilon)[1]

        return report

    def mom8_creator(self, img, range_list):
        """ function that takes and image and turns it into a .png for
            weblog
        """
        immoments(imagename = img, moments = 8, outfile = img+'.moment8')
        imview(raster={'file': img+'.moment8', 'range': range_list}, \
            out = {'file': img+'.moment8.png'})
        os.popen('mogrify -trim '+img+'.moment8.png')

    def cube_profile_fit(self, img, max_loc, nchan):
        """ function that will retrieve a profile for cubes at the max position
            and create a png showing the profile plot; must be called with
            image already opened
        """
        pyplot.clf()
        box = str(max_loc[0])+','+str(max_loc[1])+','+str(max_loc[0])+','+str(max_loc[1])
        profile = self._myia.fitprofile(box=box)['gs']['amp'][0][0][0][0][0]
        X = self._myia.getchunk(blc=max_loc, trc=max_loc, axes=[0,1])[0][0][0]
        pyplot.title('Frequency Profile at Max Value Position')
        pyplot.xlabel('Channel Number')
        pyplot.xlim(0,(nchan+1))
        pyplot.ylabel('Amplitude (Jy/Beam)')
        pyplot.plot(X)
        pyplot.savefig(img+'.profile.png')
        pyplot.clf()

        return profile


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
        img = file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, imagename=file_name+'0', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS'
            ' 00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526744GHz', width='0.244174087287MHz', \
            outframe='LSRK', pblimit=0.2, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            deconvolver='hogbom', usepointing=False, restoration=False, \
            pbcor=False, weighting='briggs', restoringbeam='common', \
            robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=True)

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
                calcres=False, calcpsf=False, savemodel='none', \
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


        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/standard_cube.image.crtf')

        exp_im_stats = {'com_bmaj': [False, 8.509892605313942],
            'com_bmin': [False, 5.950050676606115],
            'com_pa': [False, 72.54607919421503],
            'npts': [True, 3251200],
            'npts_unmasked': [True, 1522476.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.94608676433563232],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, -0.70467984676361084],
            'min_val_pos': [True, [18, 57, 0, 374]],
            'im_rms': [False, 0.143986095161],
            'rms_per_chan': [False, [0.12634143789549535, 0.14184103443455054, 0.146980848491423, 0.11560714721775701, 0.16661204089962214, 0.14476185721954618, 0.1316304465304373, 0.11433992978005528, 0.1442336908847237, 0.16543686095641913, 0.11873031602382604, 0.15284339493903618, 0.17195923304927352, 0.15046376352374186, 0.14070242684409134, 0.13744696370230045, 0.11158401597402363, 0.12687194824894843, 0.14295610599459097, 0.16862585818094997, 0.13008023644212954, 0.14186757490541813, 0.1169917541756216, 0.13693402712045005, 0.16194773183534902, 0.13634870573065122, 0.13860445389090703, 0.1309492701136035, 0.14819650092974662, 0.15030436484700252, 0.14931368127692368, 0.11984058396768074, 0.13119726238629406, 0.1278997483823513, 0.15680618961802364, 0.14782343113879803, 0.1452811146145065, 0.14962350388870774, 0.12727392822661138, 0.12403611951801675, 0.13971310504808565, 0.14443747442976043, 0.13947457857066817, 0.14603448593891352, 0.1487357653330162, 0.13728792717834695, 0.12754218448487697, 0.13472363429296094, 0.17318897000268654, 0.15007875971445414, 0.1210452212469422, 0.15977553256440455, 0.13077138200444186, 0.12679151267047647, 0.12091204082027505, 0.1338966333695089, 0.13556991771575277, 0.15456345918134376, 0.12044465611280492, 0.14836982561861023, 0.1349896116866282, 0.15311214064438922, 0.11251497655504887, 0.134867796496014, 0.13574313457554038, 0.14582224580240324, 0.12753531221719416, 0.15335445312643003, 0.13482732612307324, 0.1622050903445585, 0.13260306174268546, 0.1345326608100535, 0.16404765102131583, 0.13449430188802702, 0.14543809289295098, 0.1606584196112734, 0.12484651484486906, 0.16251383851634701, 0.13756025624117688, 0.13165353467440083, 0.1308248320448295, 0.14752778635690292, 0.1274645256107852, 0.16421712463271607, 0.15255317243782812, 0.1497707840063393, 0.11911825364867326, 0.14541033702618353, 0.1659723426787793, 0.1554971226410762, 0.14703675741501698, 0.12325846980328654, 0.15070706791866434, 0.1243073669840061, 0.13646642844468243, 0.1301143392639293, 0.12734602178400867, 0.1553600823593344, 0.15035594210430997, 0.11530605847413075, 0.1611567346343003, 0.12221832384850957, 0.14207389319672978, 0.14522516033398006, 0.1345322788758837, 0.1486176245373929, 0.15765848896613346, 0.1308440759384876, 0.1466820831226493, 0.13598865468593319, 0.15187538855740168, 0.1478468013010444, 0.1383732515889412, 0.1276861625889527, 0.11697230161534232, 0.13739607388184524, 0.11303259344169146, 0.1361001584583741, 0.12857356426667815, 0.1437570752313611, 0.13169397143643052, 0.15326431411050365, 0.12383180315967929, 0.1526310794015497, 0.14746177769245866, 0.15194893390457265, 0.1421630320154613, 0.15662308690272084, 0.12239198421329735, 0.12071542153915982, 0.14268554321174182, 0.13489697242976567, 0.15127855443293006, 0.1542443819690316, 0.15752918577920158, 0.11478434733366248, 0.17298964180575135, 0.13177526480150695, 0.12236732291938952, 0.15625856947990782, 0.13687165189461548, 0.1536631153928859, 0.14669563803395924, 0.1277170908624889, 0.14966567842171496, 0.12823515897560267, 0.13577828413547297, 0.16140169123660877, 0.13133284404676335, 0.14223570583416104, 0.1603292311222728, 0.10759630495294702, 0.15787039978749143, 0.1327200609847152, 0.14655899389809018, 0.14008820956915727, 0.1442107348583108, 0.1317943450568934, 0.12972989243424995, 0.1625036947147829, 0.12241712383574781, 0.14998173521745944, 0.13352731228428555, 0.1741676258276787, 0.15545996482656257, 0.13121844421079562, 0.1389256768353536, 0.1475992903718036, 0.14205849908080379, 0.14975427804440275, 0.1532491403618113, 0.12531915969323904, 0.14153689035122896, 0.16741877503811964, 0.1355536447212321, 0.12548585056941425, 0.16334800417248366, 0.14220841606737944, 0.1376802362928535, 0.1394159389365598, 0.1533008119644231, 0.12568227593323275, 0.14138024496799945, 0.14688836279261966, 0.12037367892758656, 0.12335138886587714, 0.16740640885840646, 0.11756235238942149, 0.13221931449560975, 0.14605469946826174, 0.12287649136200192, 0.13900407591276098, 0.1477935699475207, 0.14723640198504923, 0.12637771862286276, 0.14264989851200444, 0.14188497863070984, 0.1517498748029243, 0.1745550071541481, 0.14693061119966988, 0.12180541963696558, 0.17178472812899895, 0.134842796032342, 0.1587769050427257, 0.16022475326023228, 0.12598385136025822, 0.12173065475536829, 0.1358700032273519, 0.12249230371601251, 0.1320416693266833, 0.1380195667444624, 0.17036819494074398, 0.14449179298441997, 0.1363579047545357, 0.15814587607932587, 0.1387404461979901, 0.13421674959986293, 0.1221729254232071, 0.15007074873391474, 0.1519841002224019, 0.17405910974305452, 0.10810253208919626, 0.14404509620995673, 0.12925102011532486, 0.13284702789448985, 0.16316517507291742, 0.18004246985230368, 0.12352109323053732, 0.13438971701846103, 0.14110722423724795, 0.15240505247738928, 0.16299890474660164, 0.13862726296963718, 0.13653417057201506, 0.15748574227626927, 0.13330507817933285, 0.11630210517195279, 0.14310200319865532, 0.16947649357417122, 0.19276632648628003, 0.1442624172150719, 0.12588748723036136, 0.13766261231222038, 0.15501574319393477, 0.1467664214746981, 0.1437631603776764, 0.13281269178755073, 0.1499498907469051, 0.1547831067839161, 0.1650203851926907, 0.19221241068202002, 0.17354684080024252, 0.1914610395870755, 0.1870175530216679, 0.17581044678251068, 0.1969886211075383, 0.1793297247200161, 0.1650936078172776, 0.14492844003488278, 0.14169417048994912, 0.13913741690930143, 0.15199069574737809, 0.14420138013804454, 0.13427068237336662, 0.14854837369704055, 0.16108337901230768, 0.14392226756737841, 0.17357790177064894, 0.12582205913761169, 0.14917348406521977, 0.14570283814332685, 0.1494709027791571, 0.15333214078415874, 0.13788344445589443, 0.15113305293127033, 0.15910328490835218, 0.1258386682524315, 0.14889037473182778, 0.12206410303883597, 0.14151660747434333, 0.12015317625283857, 0.13619775353885638, 0.1360718460404539, 0.12638979179451687, 0.13793807234055996, 0.11510797437711365, 0.151683286728862, 0.12186947178372404, 0.1334587198282186, 0.1416652191079726, 0.17466019895415022, 0.1378517218379097, 0.11249994464540782, 0.1279851283688596, 0.14018299952556243, 0.14434793746130153, 0.16650070810701043, 0.16748683720689606, 0.11940475115442413, 0.1276077403883665, 0.12359063458201924, 0.1412078454825856, 0.14551852706758697, 0.14986530033150053, 0.1325683624037802, 0.13820026200226213, 0.12812701294141524, 0.12753530671589502, 0.13497760358639382, 0.170481828949963, 0.13917782659135283, 0.14415530157163367, 0.17676595928739716, 0.16642659928559336, 0.1488633024931477, 0.13762287309163865, 0.13173759984975023, 0.13142953909930827, 0.13550297890119667, 0.1378187724239619, 0.13522211241434626, 0.16692597968173054, 0.12681333395848324, 0.14653010389914287, 0.12049158054527685, 0.1394775351353454, 0.1308254468074764, 0.15490863570934446, 0.13846726479055535, 0.12889021958687075, 0.1475440113411787, 0.1317061942800524, 0.13705790461652367, 0.13974036959646274, 0.14395769504564201, 0.16390950740528484, 0.14702614449178092, 0.11790240659365156, 0.15825207098676003, 0.1572387243129431, 0.16049823864649457, 0.15795518128891234, 0.12433544601523717, 0.12148637251440078, 0.17047767727369215, 0.1542909493696503, 0.14083894029456667, 0.1273419049433489, 0.12328954463396018, 0.14383278211432127, 0.16043058197358323, 0.12709149940547382, 0.14015718181807055, 0.13290234816222832, 0.14241483672006272, 0.13977682956445464, 0.19082240791944033, 0.1362551319469732, 0.13878314241029005, 0.13024580725200274, 0.17584969348391408, 0.1585665364841856, 0.1753088252644526, 0.14469549936190076, 0.12561131070797885, 0.1275068394660626, 0.1246509966110453, 0.1353771632037516, 0.14968376152813326, 0.1306367587605417, 0.17509696143235529, 0.12303417070683348, 0.12984385988761357, 0.15655156396975772, 0.15158283415352666, 0.14190571663235674, 0.16013330804491582, 0.15701423529031192, 0.1437377861978616, 0.13178714397605581, 0.14896642830351547, 0.1543110828211962, 0.12729972586137084, 0.13258135641591057, 0.1297547997545033, 0.13623621902113375, 0.16832161541183457, 0.1515993553463645, 0.14230780146154384, 0.15752128344423338, 0.16153801377050137, 0.17727249054646335, 0.17081959421469814, 0.132415345138367, 0.1911207789325767, 0.12167933736058796, 0.16199351646871693, 0.12220003129636688, 0.12276819690997134, 0.1615915263092078, 0.12453980969317847, 0.13583553980995322, 0.15247696319780732, 0.1588409406770753, 0.12523094487161238, 0.14442835547547206, 0.14678401935124402, 0.16646235670431703, 0.13318109022977143, 0.1296318756165084, 0.17474714467039837, 0.14996240468941135, 0.152459961330088, 0.14174866849523865, 0.1531711479592131, 0.13140647178793072, 0.12876054234189516, 0.16073620618951742, 0.1551915897989922, 0.13258798290709237, 0.12079795445447011, 0.17410998739129824, 0.15215648462348455, 0.11414501724048388, 0.1583626102822897, 0.14397582688024274, 0.14206199694234994, 0.1301244502774394, 0.14832273898073517, 0.14965582403580915, 0.11974680479280232, 0.130301576062527, 0.1501180642039695, 0.144040093741809, 0.14644461407160658, 0.16466625212890748, 0.12145879983950111, 0.15515579044778377, 0.11925841619372717, 0.14869156838984582, 0.14241952971976954, 0.14769053467507798, 0.1574041746549612, 0.13729987849339584, 0.13259542009008363, 0.12816033699683535, 0.14196238422903054, 0.14271636018596165, 0.1303221375620509, 0.16261332868157302, 0.16340560251777758, 0.1496448141003338, 0.13870981454297002, 0.12581695569707782, 0.10706545850873779, 0.16477345617357922, 0.13501062540902478, 0.15224041747140973, 0.14304328175572253, 0.1730320818474802, 0.13262178336064356, 0.13626925589130173, 0.19930734468055353, 0.12840465325131772, 0.14565284364647132, 0.14342338225449644, 0.14414011525387105, 0.1314767224882364, 0.14667809040912927, 0.15534831496680646, 0.14770322269613237, 0.18150555408581706, 0.14175110723675505, 0.13126842096913788, 0.1821782599052581, 0.1230212008372949, 0.13879304863372663, 0.1539616899583395, 0.14735461585329224, 0.1103039067487409, 0.1400589610059953, 0.12936178694337058, 0.1185505506087342, 0.15904169231653986, 0.12089109358577753, 0.1479414523691717, 0.13846627238497242, 0.15232750215647736, 0.1318383069114147, 0.13798313113493177, 0.12800661214057107, 0.15270008178936095, 0.12254175320261147, 0.13897505123792747, 0.13607974358604047, 0.15255689754961757, 0.11221218682645963, 0.14990854525350425, 0.15160985722773732, 0.13575646404498873, 0.14576400908648826, 0.15836555616934264, 0.13583817797589506, 0.12589022941650596, 0.15590488932502036, 0.16965183711831752, 0.1483222833838575, 0.14620965767678082]],
            'im_sum': [False, 168.242297979],
            'regn_sum': [False, 72.3549563158],
            'npts_real': [True, 3251200],
            'profile': [False, 0.90354546054915941],
            'im_fit': [False, [0.909677723621573, 13.925317200593012, \
                       7.225723641827008]],
            'im_fit_loc': [True, [254, 220.31469461816917]],
            'im_fit_pix': [False, [38.303527208322556, 37.46382702762046]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 70, 0, 0]), \
                (img+'.image', False, [40, 71, 0, 0]), \
                (img+'.image', True, [10, 40, 0, 0]), \
                (img+'.image', False, [9, 40, 0, 0])])

        # .image report
        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')


        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 3251200],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'mask_pix': [True, 437],
            'mask_regns': [True, 1],
            'npts_real': [True, 3251200]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/standard_cube.pb.crtf')

        exp_pb_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 1522476.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'min_val_pos': [True, [25, 13, 0, 396]],
            'im_rms': [False, 0.578238326026],
            'npts_0.2': [True, 1522476],
            'npts_0.5': [True, 736092],
            'npts_real': [True, 3251200],
            'pb_fit': [False, [1.0308127949041446, 46.61751391582679, \
                       46.61253844001269]],
            'pb_fit_loc': [True, [254, 220.31469461816917]],
            'pb_fit_pix': [False, [40.00032808200995, 40.00099739969875]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/standard_cube.psf.crtf')

        exp_psf_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.218764916062],
            'min_val_pos': [True, [1, 16, 0, 503]],
            'im_rms': [False, 0.136036099793],
            'im_sum': [False, 7472.57665916],
            'npts_real': [True, 3251200],
            'psf_fit': [False, [1.0959863390592945, 7.672871552789668, \
                        5.141790170376213]],
            'psf_fit_loc': [True, [254, 220.31469461816917]],
            'psf_fit_pix': [False, [39.99880225653659, 39.99524870969922]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/standard_cube.residual.crtf')

        exp_resid_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 1522476.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.785612404346],
            'max_val_pos': [True, [42, 43, 0, 256]],
            'min_val': [False, -0.704679846764],
            'min_val_pos': [True, [18, 57, 0, 374]],
            'im_rms': [False, 0.143918523224],
            'im_sum': [False, 124.317946204],
            'regn_sum': [False, 33.355380173],
            'npts_real': [True, 3251200]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/standard_cube.model.crtf')

        exp_model_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.286023736],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000249846621096],
            'im_sum': [False, 0.92636379227],
            'regn_sum': [False, 0.92636379227],
            'mask_non0': [True, 6],
            'npts_real': [True, 3251200]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 94.4766769409],
            'max_val_pos': [True, [0, 0, 0, 17]],
            'min_val': [False, 94.4766464233],
            'min_val_pos': [True, [0, 0, 0, 449]],
            'im_rms': [False, 94.476661707],
            'npts_real': [True, 508]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[0.3, 1.0])
            self.mom8_creator(img+'.residual', range_list=[0.3, 1.0])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[testname]['images'].append(img+'.image.profile.png')

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=True)

        # move files to iter1
        print('Copying iter0 files to iter1')
        self.copy_products(file_name+'0', file_name+'1')

        print("STARTING: iter1 routine")

        # iter1 (restart)
        tclean(vis=self.msfile, field='2', spw=['0:113.893653412~'
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
            savemodel='none', parallel=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/standard_mfs.image.crtf')

        exp_im_stats = {'com_bmaj': [False, 18.0536975861],
            'com_bmin': [False, 10.3130340576],
            'com_pa': [False, 86.4390563965],
            'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.03620968014],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.0025309077464],
            'min_val_pos': [True, [51, 36, 0, 0]],
            'im_rms': [False, 0.00317099023963],
            'im_sum': [False, 1.72629247889],
            'regn_sum': [False, 1.70480877149],
            'npts_real': [True, 6400],
            'im_fit': [False, [0.0368173095435, 17.888484296, 9.90872728645]],
            'im_fit_loc': [True, [0, 107.84024497577988]],
            'im_fit_pix': [False, [40.2022706573, 40.0784833662]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 74, 0, 0]), \
                (img+'.image', False, [40, 75, 0, 0]), \
                (img+'.image', True, [6, 40, 0, 0]), \
                (img+'.image', False, [5, 40, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 6400],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 334],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/standard_mfs.pb.crtf')

        exp_pb_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'min_val_pos': [True, [33, 6, 0, 0]],
            'im_rms': [False, 0.577629540243],
            'npts_0.2': [True, 3793],
            'npts_0.5': [True, 1813],
            'npts_real': [True, 6400],
            'pb_fit': [False, [1.0467417343495562, 92.30725376920157, \
                       92.30671415384658]],
            'pb_fit_loc': [True, [0, 107.84024497577988]],
            'pb_fit_pix': [False, [39.99973335198128, 40.00036927599604]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/standard_mfs.psf.crtf')

        exp_psf_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.170876815915],
            'min_val_pos': [True, [45, 14, 0, 0]],
            'im_rms': [False, 0.111753598174],
            'im_sum': [False, 13.4319936319],
            'npts_real': [True, 6400],
            'psf_fit': [False, [1.097246906267534, 15.626704258596684, \
                        9.180460042245928]],
            'psf_fit_loc': [True, [0, 107.84024497577988]],
            'psf_fit_pix': [False, [40.01095621317507, 39.995429898147734]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/standard_mfs.residual.crtf')

        exp_resid_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.00680131698027],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.00267856917344],
            'min_val_pos': [True, [51, 36, 0, 0]],
            'im_rms': [False, 0.00119958416595],
            'im_sum': [False, 0.167713507041],
            'regn_sum': [False, 0.263926472682],
            'npts_real': [True, 6400]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.resid')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/standard_mfs.model.crtf')

        exp_model_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0253213085234],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000320901735817],
            'im_sum': [False, 0.0295509714633],
            'regn_sum': [False, 0.0295509714633],
            'mask_non0': [True, 2],
            'npts_real': [True, 6400]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 3208318.25],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 3208318.25],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 3208318.32244],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[-0.003, 0.04])
            self.mom8_creator(img+'.residual', range_list=[-0.003, 0.04])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=True)

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
            savemodel='none', parallel=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.image.tt0.crtf')

        exp_im_stats = {'com_bmaj': [False, 18.0536975861],
            'com_bmin': [False, 10.3130340576],
            'com_pa': [False, 86.4390563965],
            'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0375871770084],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.00482273893431],
            'min_val_pos': [True, [40, 62, 0, 0]],
            'im_rms': [False, 0.00341862729159],
            'im_sum': [False, 1.82207041922],
            'regn_sum': [False, 1.7745978582],
            'npts_real': [True, 6400],
            'im_fit': [False, [0.03820503393292664, 18.02503733453906, \
                       9.894877124019276]],
            'im_fit_loc': [True, [0, 107.84024497577988]],
            'im_fit_pix': [False, [40.2022706573, 40.0784833662]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [40, 74, 0, 0]), \
                (img+'.image.tt0', False, [40, 75, 0, 0]), \
                (img+'.image.tt0', True, [6, 40, 0, 0]), \
                (img+'.image.tt0', False, [5, 40, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image.tt0')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 6400],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 332],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.pb.tt0.crtf')

        exp_pb_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'min_val_pos': [True, [33, 6, 0, 0]],
            'im_rms': [False, 0.577629540243],
            'npts_0.2': [True, 3793],
            'npts_0.5': [True, 1813],
            'npts_real': [True, 6400],
            'pb_fit': [False, [1.0467417343495562, 92.30725376920157, \
                       92.30671415384658]],
            'pb_fit_loc': [True, [0, 107.84024497577988]],
            'pb_fit_pix': [False, [39.99973335198128, 40.00036927599604]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb.tt0')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.psf.tt0.crtf')

        exp_psf_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.170876815915],
            'min_val_pos': [True, [45, 14, 0, 0]],
            'im_rms': [False, 0.111753598174],
            'im_sum': [False, 13.4319936319],
            'npts_real': [True, 6400],
            'psf_fit': [False, [1.097246906267534, 15.626704258596684, \
                        9.180460042245928]],
            'psf_fit_loc': [True, [0, 107.84024497577988]],
            'psf_fit_pix': [False, [40.01095621317507, 39.995429898147734]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf.tt0')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            region_file = data_path+
            'region_files/standard_mtmfs.residual.tt0.crtf')

        exp_resid_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.00681467074901],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.00269495742396],
            'min_val_pos': [True, [51, 36, 0, 0]],
            'im_rms': [False, 0.00121173651371],
            'im_sum': [False, 0.146450882651],
            'regn_sum': [False, 0.279778897976],
            'npts_real': [True, 6400]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0')

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', region_file = \
            data_path+'region_files/standard_mtmfs.model.tt0.crtf')

        exp_model_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0257838927209],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000328232262345],
            'im_sum': [False, 0.0307541997172],
            'regn_sum': [False, 0.0307541997172],
            'mask_non0': [True, 2],
            'npts_real': [True, 6400]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model.tt0')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 3208318.25],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 3208318.25],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 3208318.32244],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0')

        # .alpha report
        alpha_stats_dict = self.image_stats(img+'.alpha', region_file = \
            data_path+'region_files/standard_mtmfs.alpha.crtf')

        exp_alpha_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 2009.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 55.2813644409],
            'max_val_pos': [True, [28, 15, 0, 0]],
            'min_val': [False, -55.0754890442],
            'min_val_pos': [True, [58, 63, 0, 0]],
            'im_rms': [False, 16.165969128],
            'im_sum': [False, 19177.0033636],
            'regn_sum': [False, 2.50601291656],
            'mask_non0': [True, 2009],
            'npts_real': [True, 6400]}

        report9 = self.stats_compare(exp_alpha_stats, alpha_stats_dict, \
            '.alpha')

        # .alpha.error report
        error_stats_dict = self.image_stats(img+'.alpha.error', region_file = \
            data_path+'region_files/standard_mtmfs.alpha.error.crtf')

        exp_error_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 2009.0],
            'freq_bin': [True, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 78.1796569824],
            'max_val_pos': [True, [28, 15, 0, 0]],
            'min_val': [False, 0.0368266552687],
            'min_val_pos': [True, [20, 3, 0, 0]],
            'im_rms': [False, 22.7563186515],
            'im_sum': [False, 35943.3947913],
            'regn_sum': [False, 0.745403200388],
            'npts_real': [True, 6400]}

        report10 = self.stats_compare(exp_error_stats, error_stats_dict, \
            '.alpha.error')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10

        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image.tt0', range_list=[-0.005, 0.04])
            self.mom8_creator(img+'.residual.tt0', range_list=[-0.005, 0.04])
            test_dict[testname]['images'].extend( \
                (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            parallel=False)

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
            parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/standard_cube_eph.image.crtf')

        exp_im_stats = {'com_bmaj': [False, 4.49769604687],
            'com_bmin': [False, 3.3237527868],
            'com_pa': [False, 87.0964067383],
            'npts': [True, 6400000],
            'npts_unmasked': [True, 3233000.0],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 3.2168185710906982],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.35654923319816589],
            'min_val_pos': [True, [58, 39, 0, 490]],
            'im_rms': [False, 0.0560069684214],
            'rms_per_chan': [False, [0.04643669957034992, 0.04934490831784589, 0.04896325882011928, 0.045246930831289396, 0.04076479281458919, 0.04066843326216516, 0.043761578500556515, 0.04370275108520549, 0.0445865859431914, 0.04802272981355104, 0.046277171156772115, 0.04421698263536416, 0.04156787733109207, 0.04051142659268897, 0.039050378679608116, 0.04274874247664772, 0.04252468925473175, 0.043598688189381424, 0.04212993795874748, 0.04291563807365895, 0.04663133269893015, 0.044995182687279166, 0.046652174816757166, 0.04680608296946909, 0.04689052908184202, 0.04461073330219135, 0.043414456008460174, 0.04697870918788608, 0.04121882663735321, 0.04252376733369552, 0.04487817651244527, 0.046179300427229446, 0.05082631575222008, 0.05096767735717168, 0.05157594918063651, 0.053087774183373154, 0.04951326648446268, 0.045951663538404505, 0.04496180664310007, 0.04043775122677042, 0.03863681189715942, 0.04198142992416419, 0.046355293794327235, 0.045532964002632854, 0.04173551104543082, 0.04007043615342306, 0.043528066688108526, 0.04885149532839892, 0.050238176414982176, 0.04596874411684181, 0.04195383438441166, 0.04626081402776205, 0.0492666486971535, 0.05150304313267975, 0.04462335509905232, 0.04118038686235944, 0.04780955467331849, 0.04752930157265889, 0.044205684120537145, 0.042692262207786574, 0.04261267884394058, 0.041558523019057145, 0.04598449726235426, 0.05011069650686886, 0.04700991930989118, 0.04754284812449084, 0.04777701297333353, 0.044937587455719685, 0.04545473401615083, 0.045628517495069715, 0.046974171924376536, 0.04791987552049449, 0.051103967401651365, 0.05078302116118002, 0.04916354347812644, 0.045185515685354324, 0.04149383314091068, 0.04252371812661069, 0.04659269834497681, 0.0513482193106561, 0.04888969235184415, 0.04687727401908186, 0.04365857821009889, 0.041131483603749484, 0.043500582322170855, 0.0462475887996055, 0.050947374919630015, 0.04740160325995145, 0.0470193902756731, 0.04820730198011588, 0.044227183958923326, 0.042122827379657904, 0.042941860821600106, 0.041389301431217086, 0.04059296483437863, 0.04341043774321849, 0.0465241236562757, 0.050343051630505134, 0.05028274999938501, 0.04155552969131352, 0.039886981070265566, 0.046645263962396155, 0.051459870958774946, 0.04958652396902705, 0.045265090075850355, 0.04434772609236033, 0.049116532851638736, 0.05221305702489155, 0.04706088732472824, 0.04602886320332381, 0.04613588584535028, 0.04371228016495952, 0.03982330273063204, 0.04193092175965313, 0.049198482375183476, 0.04559050484741412, 0.04460361880216926, 0.04353058920316075, 0.045899887482949125, 0.046768575899857695, 0.04587087214138045, 0.04843732516362743, 0.04694015149126692, 0.04388274962922123, 0.043307500885936795, 0.046105744604161274, 0.050346868624568816, 0.047200642319448524, 0.043589961639616794, 0.04335470016064548, 0.04324714935672587, 0.04454941674749486, 0.04737927686773757, 0.05082163506683548, 0.04848585994107168, 0.0411642812817369, 0.04576066235499968, 0.04473033253026546, 0.04278733014801023, 0.04134775293581778, 0.043486359980697865, 0.046006054192836494, 0.043096775414862926, 0.04208176804672694, 0.04049787676913971, 0.03776413600971469, 0.04096302522240697, 0.04113701125252589, 0.04084318844319109, 0.03993326431649316, 0.04205014226137875, 0.04165528446233006, 0.04743968297025062, 0.0475521877300109, 0.04475605845619125, 0.04401328458943738, 0.041007266936161024, 0.04599262159243273, 0.046865999859658745, 0.04827691420906647, 0.05228095627985328, 0.05380087529459753, 0.045382540014794075, 0.04408051931220572, 0.050880525362058296, 0.04633788067452101, 0.04452905578109357, 0.04791010864887247, 0.04529597872689918, 0.04525854238251303, 0.046747443761763714, 0.03933074371811037, 0.03598068749641659, 0.041265993728334616, 0.0430330960969908, 0.04067461852978585, 0.041035175357556804, 0.045119814280047195, 0.047730924757603155, 0.04759626186493997, 0.04727124180084044, 0.04424468462136652, 0.04251731984028063, 0.04176012565602086, 0.04875414512653576, 0.049819399667323835, 0.04249500530271126, 0.04048631054645446, 0.04019127484662825, 0.04447260323303356, 0.04733720183513169, 0.04759559578370847, 0.046415899414605034, 0.04081800248138881, 0.040835322610187075, 0.04078854436028061, 0.04084601596595794, 0.04280266371654101, 0.04384684014451019, 0.04357951587924562, 0.04493228494689293, 0.044139830914717695, 0.0435869145723245, 0.04764017721689463, 0.048284337537651045, 0.04666980732346283, 0.04733310508129286, 0.044116863618867234, 0.04559536238825363, 0.04448530889797417, 0.04612795674019243, 0.04340643325575028, 0.045288094093386654, 0.04878446017187089, 0.04885753033555586, 0.04941586336470389, 0.046422340562798904, 0.044323117513462226, 0.047708631750805794, 0.04795955693138798, 0.04403887278767426, 0.043242816075875455, 0.04505113360417497, 0.045971733801389685, 0.04692591382520756, 0.04526334465868224, 0.0402667998958263, 0.041599834702646075, 0.04137294130862181, 0.04270784407047202, 0.041988761534425746, 0.04231406077286222, 0.051940796708851544, 0.05046188340262736, 0.04423891883374399, 0.04311071331069119, 0.041819134036779615, 0.04229742503352504, 0.04170398273505288, 0.04562529179144891, 0.04874608112356234, 0.045505038077315624, 0.04507038545756651, 0.04516625412224887, 0.04788075845520247, 0.04825376763517949, 0.04906145926913566, 0.049139822772542675, 0.04557852751238228, 0.04720599191117688, 0.043884980387716814, 0.04140054881324252, 0.04211615786203093, 0.043127698114725284, 0.04911141042084732, 0.05107081928037268, 0.044796984464178013, 0.04495833810462226, 0.04730374501586064, 0.04667756324575635, 0.045159866333060114, 0.04042904164023826, 0.04090576051606706, 0.046863276603866556, 0.04659173239577013, 0.04376744319773422, 0.04522289686498856, 0.04473657323218306, 0.042819092227464356, 0.04113471729537815, 0.04591263892335558, 0.04798013561676214, 0.04731129275067027, 0.04827915688307468, 0.04467434152154626, 0.043046478969549565, 0.042103651285131743, 0.04240517430258763, 0.04373979060557242, 0.045668840553286565, 0.046285137568665945, 0.04211919000123673, 0.043169481337927564, 0.05065536603759789, 0.04964137715560091, 0.04872082434854346, 0.048823130419572945, 0.044647556725311885, 0.04686071279804561, 0.04898506218686401, 0.04248479479575401, 0.041238959229725446, 0.041711815950868036, 0.04448585314572631, 0.045569453998905156, 0.043782848170480226, 0.03714516951598358, 0.041388339824126735, 0.04794368187941998, 0.04774424802037788, 0.04185035411590102, 0.04187212720629191, 0.04061630982271867, 0.04014055898470098, 0.03755902905378686, 0.041731525656284275, 0.04552077472723898, 0.045384557347911, 0.04561425607077976, 0.04396240289328813, 0.04560866514676667, 0.046018010543353076, 0.041758975644607906, 0.044584517396196746, 0.044393563966261325, 0.044515088147448055, 0.048423700561444925, 0.045279573442018585, 0.04361878407496164, 0.0468957317021233, 0.0498702556814145, 0.04923825866076725, 0.04933686700621225, 0.04709859751278079, 0.04167128843031174, 0.04061805103901589, 0.03885135890432584, 0.04321653772782609, 0.044339345524985654, 0.04449007623693989, 0.0480083250328015, 0.04728131076450639, 0.042495683485785384, 0.04335727837754768, 0.042470766390855536, 0.04273025309548049, 0.044332202260225044, 0.042474455866522876, 0.04016055605615106, 0.042300011432415185, 0.04416254354402529, 0.04563914350851595, 0.04686015983173452, 0.04631455813837692, 0.04498015818897669, 0.043175559477057215, 0.04383425133770565, 0.04159113410245179, 0.04489350716614797, 0.04650094370947378, 0.043072941586427806, 0.042205330071801495, 0.039883474726589395, 0.044825203490698275, 0.04528892805777943, 0.04282610154191627, 0.04349512611998853, 0.04247412091291035, 0.04415929041944228, 0.048223187057763434, 0.04726879001889971, 0.04394868050355618, 0.04415635728956364, 0.04501345441129184, 0.04418666629872543, 0.04668576976453467, 0.04556394565391764, 0.04741009982528531, 0.04740140452150073, 0.04583037551978373, 0.046567557496669576, 0.04677469032755722, 0.046169211332555686, 0.043724071095036704, 0.04287974465665924, 0.04304122053501154, 0.04102254270029318, 0.04361672048216947, 0.046468143734660614, 0.04638227729706821, 0.04188419989140981, 0.04141905130776169, 0.04601534424565727, 0.04315761018707787, 0.04178065317591785, 0.04181189318236191, 0.04433161553697957, 0.04835683069281614, 0.05122027651258199, 0.04888875491912803, 0.04709415480984611, 0.03982308750616121, 0.03931970337124684, 0.041135584388753595, 0.04103474060272379, 0.042688688753807644, 0.045098100375485284, 0.044081581919921765, 0.044867634085172664, 0.04376194005125083, 0.04112353809689816, 0.04165256642058881, 0.04113864076272384, 0.04194089691466448, 0.0426357193530805, 0.044169091773530005, 0.044706911723422046, 0.045676954240245315, 0.0475735363640219, 0.045252129087844414, 0.04327843818779502, 0.040787748524881914, 0.04349982070782621, 0.04623701927905451, 0.05003818830094732, 0.0492637202304341, 0.05030172011167053, 0.05041856505148443, 0.04865025938675271, 0.04242114836177711, 0.04332913791842073, 0.04417555276786257, 0.046511419107420855, 0.047855191840534814, 0.0447112493727955, 0.04291434192759193, 0.042000052238535396, 0.047208554560098065, 0.04991329683989256, 0.043522094691402766, 0.04728530679991288, 0.043409909773413066, 0.04413453137291575, 0.04378692865928274, 0.04318986388883555, 0.04768566840485644, 0.0493043820754111, 0.045915063051690456, 0.04372051075206559, 0.04177309296445956, 0.04594702153236263, 0.043136187115026786, 0.03989584010763349, 0.042052831627502514, 0.04395385380777047, 0.044490797911686523, 0.05050721479049005, 0.05280016432973104, 0.045817643009786586, 0.043342516330605346, 0.03668099527899864, 0.037354424964758944, 0.04417339767971324, 0.048705771265993306, 0.04417246001889241, 0.03920079992318, 0.04236443562241198, 0.04569530642633145, 0.04585207146091095, 0.04393512941738649, 0.04442675712801038, 0.04313956775505375, 0.03960015974712282, 0.04595448951857105, 0.046253067292557175, 0.04764265636037589, 0.048664342796493544, 0.044490016511759545, 0.04446073291073993, 0.04325106910170598, 0.04371067608468382, 0.045815516334865726, 0.04364942129237882, 0.040614872275104554, 0.04925116326612862, 0.04622908354230324, 0.04187263674394814, 0.04246865124807005, 0.042407043107492326, 0.04395209375510593, 0.046255182065641, 0.04466090683862347, 0.043865504621226147, 0.042107514869371046, 0.045106163590565594, 0.07187955608053229, 0.1473963504142755, 0.2590806857156774, 0.35083385144133233, 0.372380544172661, 0.3385626044480845, 0.300050027009694, 0.27503491394149054, 0.2542837886157723, 0.2419519658614859, 0.2481632024383586, 0.2588597980725786, 0.2698371638935117, 0.2849747471499766, 0.26991800949889033, 0.20471858845828023, 0.12442205394588181, 0.06371996453205579, 0.04474790491119776, 0.04331778024846207, 0.043425178103955415, 0.04141722416093238, 0.04346704298140081, 0.04842483936808344, 0.04576883917836126, 0.041948629698781825, 0.041139433859260756, 0.041359269306182375, 0.04476797275140064, 0.04730401097466933, 0.04854807973917275, 0.04776793906389466, 0.04695441386976931, 0.04561123621247733, 0.04328810557431784, 0.0421941040364879, 0.04567863392551492, 0.04607553504600313, 0.04745289236281092, 0.050799363620603805, 0.047009742333265186, 0.04197312616246828, 0.045568658799004055, 0.04950290466374596, 0.05047313776113717, 0.04866930525486824, 0.041709352811006434, 0.04346991804565062, 0.04847793967620536, 0.0483407353418824, 0.04893693944483718, 0.04964036128679813, 0.04585043053261588, 0.042137232515707786, 0.04318024712159876, 0.03882668894551708, 0.03992657190878028, 0.04341166932492751, 0.04522723382472872, 0.04602373201322939, 0.043574401205969185, 0.04347286511713428, 0.04389442975326638, 0.046474225185014734, 0.04468946526873493, 0.04521254048699986, 0.045612321205862956, 0.045254081348552126, 0.04373868382604877, 0.044731104709681925, 0.04419661974812566, 0.03971135024435756, 0.040276425161050844, 0.04538473582431024, 0.04695881226018986, 0.04384551885057054, 0.04374323253548989, 0.04378309423401829, 0.042320665449719776, 0.04504178403394185, 0.04346549189106467, 0.04465447999006352, 0.045680549745101176, 0.046689972040828386, 0.04171381599865176, 0.04088252520815791, 0.04502160504474559, 0.047372363873766914, 0.05073673375295184, 0.0470445713112432, 0.04961697023359589, 0.04824373439850615, 0.0400694047756427, 0.04493917216262888, 0.05025685507158057, 0.04784519980053747, 0.045304927551867505, 0.04581202904349647, 0.044853469155612266, 0.04537774246156747, 0.0416293398608068, 0.04012314740565847, 0.039889998035148665, 0.043182397253171426, 0.048325503442644566, 0.0500434886759839, 0.04499255258371695, 0.04305837266895443, 0.041098643216759476, 0.04484789503370013, 0.04943351137544058, 0.042150912815641435, 0.04326682724906145, 0.046785962527478776, 0.042390627529364476, 0.041870969796478384, 0.04575952382122322, 0.04464385769860984, 0.04468575848348887, 0.04600236348234317, 0.04910839704277886, 0.04499056162828212, 0.0440582963307661, 0.04417858953950231, 0.045987527884238946, 0.045880564612362756, 0.04507177207972796, 0.046831333948118364, 0.04785959151143914, 0.04414535812583369, 0.04466612718088915, 0.04749075630750218, 0.04939464543074291, 0.04922078962789247, 0.04975050201601396, 0.04898735973778251, 0.04883496451383182, 0.046957028246326936, 0.04603417847246742, 0.042974412285528056, 0.03734283029861541, 0.04290257841202429, 0.0491248511597453, 0.047062472661873855, 0.04764013621951441, 0.04410154416697243, 0.04301389109429321, 0.044186292324747155, 0.04091910249679947, 0.04155120567845409, 0.04476902014742051, 0.04108120421202448, 0.037906657399788365, 0.043727289913506266, 0.042871601202105925, 0.03900637961843521, 0.036193845297806844, 0.04184542629394857, 0.04701315036410403, 0.04471123658120445, 0.04535258675302806, 0.050623824743081816, 0.05364527249445001, 0.046842276346725456, 0.038842274435487, 0.039585603998699086, 0.04186078559942754, 0.04158055169968772, 0.04788888006918037, 0.04878977348407639, 0.047802142754339305, 0.043335450547719104, 0.04465722940976265, 0.04500927942866699, 0.046016621948327174, 0.045449259413764626, 0.045445092988202714, 0.04273032326669586, 0.043450268414430906, 0.0452201565785205, 0.043983216026702174, 0.042267470283976184, 0.03913285551073241, 0.04535011019030258, 0.047601998718010685, 0.04236531816603007, 0.041294137653438645, 0.04703683659909296, 0.05060791140142113, 0.0490695284000884, 0.04659290053276082, 0.045697348140200045, 0.04077968590537448, 0.042542229405807035, 0.04833345266185202, 0.04694899703990895, 0.04753129321637921, 0.045508370850503334, 0.04245549456493364, 0.03906240157284455, 0.04077469386005198, 0.04569938738853675, 0.04440478762482793, 0.04458027978559779, 0.041291274184566136, 0.040327004579782603, 0.03772530113394665, 0.04083635914799946, 0.04480978869822854, 0.04629559979344282, 0.04138521838261721, 0.04302891111392147, 0.038370248128620034, 0.04209043287611531, 0.046217035613680665, 0.045800208835282394, 0.046851059570329256, 0.04183024348314623, 0.04036128158434568, 0.04154121832239926, 0.0413199319399299, 0.045272944823609315, 0.04445619088339592, 0.043584303677589226, 0.04009361595030041, 0.04333911757881482, 0.04634525839423941, 0.04634810861597266, 0.04732622386379959, 0.04224104408669772, 0.03908553395350114, 0.04176098735093898, 0.04365687872838128, 0.04525577395865094, 0.04636682286835371, 0.042718943926862415, 0.04378330826555344, 0.041424912747139844, 0.04158788349112557, 0.04253349156312666, 0.045348277990678455, 0.0470667024556114, 0.04534775298446029, 0.042046294836338385, 0.0411202805610577, 0.0447037980612116, 0.046461880529271776, 0.042166059117019254, 0.04066365805660356, 0.046661351791754235, 0.04804066024995413, 0.04473367304366499, 0.04118334887302664, 0.03687109983877782, 0.040931664331577064, 0.04099064458130503, 0.04277329900971529, 0.04022544876582107, 0.03842619559805958, 0.03955790365944766, 0.048853977038304616, 0.052350159708389035, 0.043366490118314324, 0.042284043375617086, 0.046783215575199044, 0.045916450282866, 0.04430257999109638, 0.04406340741696814, 0.045304210724048936, 0.04943023454037068, 0.04720492054704553, 0.046689832402328775, 0.04818384407209958, 0.0476481653943768, 0.0453458612262087, 0.040030114060635696, 0.04200276911815852, 0.0425060220901035, 0.040335283679017245, 0.04246410750799881, 0.04210916420131698, 0.040953013589040944, 0.040394321415695204, 0.04003555612982741, 0.044162707889056144, 0.0478552122030664, 0.0446140069212816, 0.04707846756265828, 0.047873583302387744, 0.047578558370327005, 0.04637937855240745, 0.0471908192573982, 0.04821753327034351, 0.04908841396554945, 0.046608834469953696, 0.044578939926042345, 0.04570005131701801, 0.042183998353329234, 0.03986806719065372, 0.041386458533410865, 0.044472471009591064, 0.04641284219986197, 0.04305330530314275, 0.04350602307903746, 0.04826561883855598, 0.04711496736718839, 0.050660924340931085, 0.05049390547388511, 0.04647009178482704, 0.04722936558353204, 0.04483411288299754, 0.04761332083261617, 0.046953461575535954, 0.04567243876512055, 0.04581749461306325, 0.04435481490424548, 0.042242332368740136, 0.04275183069544311, 0.04508095276399907, 0.04254990506220978, 0.040790627358346654, 0.04028581784318738, 0.04422786974798521, 0.04567626616109576, 0.04527964741315845, 0.04446371266660218, 0.04349445071655706, 0.04122914864839567, 0.03986658201593119, 0.042328231176402935, 0.046343793953189175, 0.04838952376368143, 0.04752407407902626, 0.0474961806565066, 0.04921827790106632, 0.0461377516515905, 0.04619995255644665, 0.04544427868597971, 0.041700931095521476, 0.04294287611850517, 0.04766903736560473, 0.04654348486860325, 0.039713977300771586, 0.039356751756926496, 0.03913220192515639, 0.039131362684817775, 0.039753132787580955, 0.042272765954301884, 0.04626841895646131, 0.044354353817838095, 0.041321187216488565, 0.042235087230061676, 0.04553720900537041, 0.04347648632966383, 0.0413496924476828, 0.04656233345688962, 0.04515312599195977, 0.04195368109573263, 0.041978131864340025, 0.04192552913436148, 0.03945625491737432, 0.04396846601436579, 0.043606450770260395, 0.04187681162413345, 0.04173631142978733, 0.04035408491782999, 0.04056327326050339, 0.0395087085528608, 0.040300122224325774, 0.04132111896757895, 0.04192198691810278, 0.04083790031822528, 0.039537336324875376, 0.04390551354786044, 0.04300962686921604, 0.04228082909060399, 0.0417565132505904, 0.0454002603419375, 0.04956280704504699, 0.04799419864665441, 0.045962094381381415, 0.0467473484273944, 0.044792070505431394, 0.04109111560252778, 0.04218760049990035, 0.04299570423703939, 0.040379807727104054, 0.04329150734012735, 0.04625143635561658, 0.04698562949105908, 0.050376418896302615, 0.04631249530889829, 0.04166623484418304, 0.040619462020136686, 0.0415672974517897, 0.04513067594110064, 0.045756965299573, 0.043529754631348574, 0.04340680883579754, 0.04370621525120933, 0.04266492708218574, 0.039113762928348804, 0.0380401026990651, 0.039976651574658253, 0.043445633528437735, 0.04549565098121662, 0.04831654339511471, 0.04572863462146055, 0.03953925284100334, 0.04053883676878803, 0.043876981147326, 0.04141695355051718, 0.041887461127441634, 0.04056677256034229, 0.04186037444981721, 0.04201413302871151, 0.04573851274855546, 0.041440183918970946, 0.03850311313601368, 0.037626510852843485, 0.04094179197839728, 0.04580007732330971, 0.04460738621891311, 0.04675278645179763, 0.0478605348944842, 0.04306252971921845, 0.046171581575171455, 0.05050074884315606, 0.042798874627819575, 0.04156016348168418, 0.044862605539110534, 0.04285813852906759, 0.04032873819014946, 0.0412731534837442, 0.04350566852020168, 0.04225372734410245, 0.04634286632670015, 0.04550365056295158, 0.045465956532854325, 0.05065921696811653, 0.04830588294157991, 0.04330708334672608, 0.04064806228847181, 0.03914909641618058, 0.041654448170194834, 0.04615595205922822, 0.04700632494550399, 0.04765357505656815, 0.04343840127677389, 0.0355150293972138, 0.03858834128735246, 0.038564138772430145, 0.043522168150675826, 0.043839343130339124, 0.03817820557819878, 0.04137383994518887, 0.04412145086969157, 0.04622839892712683, 0.04926430171308205, 0.047914247587726576, 0.04572585602454942, 0.04401998823089076, 0.04404035191781479, 0.04769100526666846, 0.048987866360454185, 0.04317200864218212, 0.042478490338744235, 0.04206542906810126, 0.04202805035461332, 0.04502553794599071, 0.04895049138390894, 0.049372283427571675, 0.04436728410362847, 0.042676705096514625, 0.0388932106584164, 0.03893482966119943, 0.04151332630361281, 0.04524423881588369, 0.04506116503171938, 0.04532681103261335, 0.0466493567116118, 0.048090064899859696, 0.04802299070388316, 0.04980025022656742, 0.04389577911463016, 0.04241263878375371, 0.04724930900263318, 0.04747258403601618, 0.042306567956874014, 0.04638399027666991, 0.04749861770039475, 0.047082672401003424, 0.04362764515742035, 0.04541169690222177, 0.04596341709748536, 0.048390183305535114, 0.04560789490195511, 0.04696830554394434, 0.05086854951618504, 0.04647093018422276, 0.03634809104465572, 0.04124641086518968, 0.04241076035782614, 0.04709708747155942, 0.04778492899706399, 0.04420316303483561, 0.04496193347190997, 0.04466609027599868, 0.04118875739174372, 0.042672531495944965, 0.044356981785933176, 0.04867307382129206, 0.046866106238546514, 0.039972418432447705, 0.04638912575140644]],
            'im_sum': [False, 2285.0251880029814],
            'regn_sum': [False, 202.79424607008696],
            'npts_real': [True, 6400000],
            'profile': [False, 2.73209990935],
            'im_fit': [False, [2.3181817983302015, 6.557786339807762, \
                       5.498814625695614]],
            'im_fit_loc': [True, [491, 354.50520934823834]],
            'im_fit_pix': [False, [46.17680391271608, 41.12935673096737]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 6400000],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'mask_pix': [True, 9362],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/standard_cube_eph.pb.crtf')

        exp_pb_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 3233000.0],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20036059618],
            'min_val_pos': [True, [20, 15, 0, 891]],
            'im_rms': [False, 0.57684102],
            'npts_0.2': [True, 3233000],
            'npts_0.5': [True, 1549000],
            'npts_real': [True, 6400000],
            'pb_fit': [False, [1.0468474913487789, 28.075184571586952, \
                       28.075184571520158]],
            'pb_fit_loc': [True, [491, 354.50520934823834]],
            'pb_fit_pix': [False, [40.0, 40.0]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/standard_cube_eph.psf.crtf')

        exp_psf_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.164189130068],
            'min_val_pos': [True, [36, 35, 0, 993]],
            'im_rms': [False, 0.0871161935921],
            'im_sum': [False, 2742.74484326],
            'npts_real': [True, 6400000],
            'psf_fit': [False, [1.1041899877875863, 3.9440139645863947, \
                        2.8581115435698514]],
            'psf_fit_loc': [True, [491, 354.50520934823834]],
            'psf_fit_pix': [False, [39.99991479119429, 39.99762629804452]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/standard_cube_eph.residual.crtf')

        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 3233000.0],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 0.313477724791],
            'max_val_pos': [True, [59, 65, 0, 492]],
            'min_val': [False, -0.335443556309],
            'min_val_pos': [True, [32, 47, 0, 491]],
            'im_rms': [False, 0.0465509555972],
            'im_sum': [False, 236.741826246],
            'regn_sum': [False, 51.1043745019],
            'npts_real': [True, 6400000]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/standard_cube_eph.model.crtf')

        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 1.3696539402],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.0673163980246],
            'min_val_pos': [True, [58, 39, 0, 490]],
            'im_rms': [False, 0.00172803556487],
            'im_sum': [False, 52.6738070883],
            'regn_sum': [False, 4.03587769717],
            'mask_non0': [True, 436],
            'npts_real': [True, 6400000]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [True, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 1009.50134277],
            'max_val_pos': [True, [0, 0, 0, 979]],
            'min_val': [False, 1007.97528076],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 1008.71815136],
            'npts_real': [True, 1000]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[0.0, 3.25])
            self.mom8_creator(img+'.residual', range_list=[0.0, 3.25])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[testname]['images'].append(img+'.image.profile.png')

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=False)

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
            savemodel='none', parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/standard_mfs_eph.image.crtf')

        exp_im_stats = {'com_bmaj': [False, 0.875946879387],
            'com_bmin': [False, 0.673672378063],
            'com_pa': [False, 88.5368652344],
            'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.03113353252],
            'max_val_pos': [True, [224, 153, 0, 0]],
            'min_val': [False, -1.01794064045],
            'min_val_pos': [True, [222, 93, 0, 0]],
            'im_rms': [False, 0.359352011299],
            'im_sum': [False, -1491.198136],
            'regn_sum': [False, 3362.95355159],
            'npts_real': [True, 82944]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [144, 266, 0, 0]), \
                (img+'.image', False, [144, 267, 0, 0]), \
                (img+'.image', True, [22, 145, 0, 0]), \
                (img+'.image', False, [21, 145, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 82944],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 82944]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/standard_mfs_eph.pb.crtf')

        exp_pb_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, 0.200061768293],
            'min_val_pos': [True, [114, 25, 0, 0]],
            'im_rms': [False, 0.569403001285],
            'npts_0.2': [True, 47329],
            'npts_0.5': [True, 22365],
            'npts_real': [True, 82944],
            'pb_fit': [False, [1.0286609217550002, 22.907692916947163, \
                       22.90769291676479]],
            'pb_fit_loc': [True, [0, 253.57442221593894]],
            'pb_fit_pix': [False, [144.0, 144.0]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/standard_mfs_eph.psf.crtf')

        exp_psf_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 82944.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, -0.0609973631799],
            'min_val_pos': [True, [140, 137, 0, 0]],
            'im_rms': [False, 0.019837926364],
            'im_sum': [False, 16.3427572285],
            'npts_real': [True, 82944],
            'psf_fit': [False, [0.9200466881631709, 0.9746655722260728, \
                        0.7626550313652652]],
            'psf_fit_loc': [True, [0, 253.57442221593894]],
            'psf_fit_pix': [False, [144.00051463175717, 144.00004766689185]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/standard_mfs_eph.residual.crtf')

        exp_resid_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.03113353252],
            'max_val_pos': [True, [224, 153, 0, 0]],
            'min_val': [False, -1.01794064045],
            'min_val_pos': [True, [222, 93, 0, 0]],
            'im_rms': [False, 0.359352011299],
            'im_sum': [False, -1491.198136],
            'regn_sum': [False, 3362.95355159],
            'npts_real': [True, 82944]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/standard_mfs_eph.model.crtf')

        exp_model_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 82944.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.0],
            'im_sum': [False, 0.0],
            'regn_sum': [False, 0.0],
            'mask_non0': [True, 0],
            'npts_real': [True, 82944]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 23234454.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 23234454.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 23234453.7637],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[-1.05, 1.05])
            self.mom8_creator(img+'.residual', range_list=[-1.05, 1.05])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=False)

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
            savemodel='none', parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.image.tt0.crtf')

        exp_im_stats = {'com_bmaj': [False, 0.875946879387],
            'com_bmin': [False, 0.673672378063],
            'com_pa': [False, 88.5368652344],
            'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.01515197754],
            'max_val_pos': [True, [224, 153, 0, 0]],
            'min_val': [False, -1.01363253593],
            'min_val_pos': [True, [222, 93, 0, 0]],
            'im_rms': [False, 0.361662708984],
            'im_sum': [False, -1558.29916138],
            'regn_sum': [False, 3388.5723138],
            'npts_real': [True, 82944]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [144, 266, 0, 0]), \
                (img+'.image.tt0', False, [144, 267, 0, 0]), \
                (img+'.image.tt0', True, [22, 145, 0, 0]), \
                (img+'.image.tt0', False, [21, 145, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 82944],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 82944]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.pb.tt0.crtf')

        exp_pb_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, 0.200061768293],
            'min_val_pos': [True, [114, 25, 0, 0]],
            'im_rms': [False, 0.569403001285],
            'npts_0.2': [True, 47329],
            'npts_0.5': [True, 22365],
            'npts_real': [True, 82944],
            'pb_fit': [False, [1.0286609217550002, 22.907692916947163, \
                       22.90769291676479]],
            'pb_fit_loc': [True, [0, 253.57442221593894]],
            'pb_fit_pix': [False, [144.0, 144.0]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.psf.tt0.crtf')

        exp_psf_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 82944.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, -0.0609973631799],
            'min_val_pos': [True, [140, 137, 0, 0]],
            'im_rms': [False, 0.019837926364],
            'im_sum': [False, 16.3427572285],
            'npts_real': [True, 82944],
            'psf_fit': [False, [0.9200466881631709, 0.9746655722260728, \
                        0.7626550313652652]],
            'psf_fit_loc': [True, [0, 253.57442221593894]],
            'psf_fit_pix': [False, [144.00051463175717, 144.00004766689185]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            region_file = data_path+
            'region_files/standard_mtmfs_eph.residual.tt0.crtf')

        exp_resid_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.03113353252],
            'max_val_pos': [True, [224, 153, 0, 0]],
            'min_val': [False, -1.01794064045],
            'min_val_pos': [True, [222, 93, 0, 0]],
            'im_rms': [False, 0.359352011299],
            'im_sum': [False, -1491.198136],
            'regn_sum': [False, 3362.95355159],
            'npts_real': [True, 82944]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', region_file = \
            data_path+'region_files/standard_mtmfs_eph.model.tt0.crtf')

        exp_model_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 82944.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.0],
            'im_sum': [False, 0.0],
            'regn_sum': [False, 0.0],
            'mask_non0': [True, 0],
            'npts_real': [True, 82944]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 23234454.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 23234454.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 23234453.7637],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # .alpha report
        alpha_stats_dict = self.image_stats(img+'.alpha', region_file = \
            data_path+'region_files/standard_mtmfs_eph.alpha.crtf')

        exp_alpha_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 29656.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 93.7298202515],
            'max_val_pos': [True, [245, 3, 0, 0]],
            'min_val': [False, -90.3992004395],
            'min_val_pos': [True, [42, 48, 0, 0]],
            'im_rms': [False, 15.9100777789],
            'im_sum': [False, -60401.9121475],
            'regn_sum': [False, 2535.505823],
            'mask_non0': [True, 29656],
            'npts_real': [True, 82944]}

        report9 = self.stats_compare(exp_alpha_stats, alpha_stats_dict, \
            '.alpha')

        # .alpha.error report
        error_stats_dict = self.image_stats(img+'.alpha.error', region_file = \
            data_path+'region_files/standard_mtmfs_eph.alpha.error.crtf')

        exp_error_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 29656.0],
            'freq_bin': [True, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 132.553985596],
            'max_val_pos': [True, [245, 3, 0, 0]],
            'min_val': [False, 0.0016036973102],
            'min_val_pos': [True, [2, 146, 0, 0]],
            'im_rms': [False, 22.5002473865],
            'im_sum': [False, 472110.530338],
            'regn_sum': [False, 179381.029198],
            'npts_real': [True, 82944]}

        report10 = self.stats_compare(exp_error_stats, error_stats_dict, \
            '.alpha.error')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image.tt0', range_list=[-1.05, 1.05])
            self.mom8_creator(img+'.residual.tt0', range_list=[-1.05, 1.05])
            test_dict[testname]['images'].extend( \
                (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='0', spw=['0'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['3'], \
            intent='CALIBRATE_BANDPASS#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[90, 90], stokes='I', \
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
            parallel=False)

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
            parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/standard_cal.image.crtf')

        exp_im_stats = {'com_bmaj': [False, 9.98569583893],
            'com_bmin': [False, 4.62464284897],
            'com_pa': [False, -86.3871307373],
            'npts': [True, 8100],
            'npts_unmasked': [True, 5041.0],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.40606927872],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, -0.0115236453712],
            'min_val_pos': [True, [50, 16, 0, 0]],
            'im_rms': [False, 0.203889562304],
            'im_sum': [False, 172.68969442],
            'regn_sum': [False, 171.573849684],
            'npts_real': [True, 8100],
            'im_fit': [False, [2.40974849537, 9.96002749264, 4.61946099469]],
            'im_fit_loc': [True, [0, 220.30076542192973]],
            'im_fit_pix': [False, [45.000405766, 45.0014155577]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [45, 85, 0, 0]), \
                (img+'.image', False, [45, 86, 0, 0]), \
                (img+'.image', True, [5, 45, 0, 0]), \
                (img+'.image', False, [4, 45, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 8100],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 407],
            'mask_regns': [True, 1],
            'npts_real': [True, 8100]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/standard_cal.pb.crtf')

        exp_pb_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 5041.0],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, 0.200092822313],
            'min_val_pos': [True, [36, 6, 0, 0]],
            'im_rms': [False, 0.577170049921],
            'npts_0.2': [True, 5041],
            'npts_0.5': [True, 2409],
            'npts_real': [True, 8100],
            'pb_fit': [False, [1.0468035426303963, 45.181424068122176, \
                       45.18134398951289]],
            'pb_fit_loc': [True, [0, 220.30076542192973]],
            'pb_fit_pix': [False, [45.000270927482546, 45.00030384048325]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/standard_cal.psf.crtf')

        exp_psf_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 8100.0],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, -0.186668172479],
            'min_val_pos': [True, [31, 41, 0, 0]],
            'im_rms': [False, 0.125651062092],
            'im_sum': [False, 40.2072189807],
            'npts_real': [True, 8100],
            'psf_fit': [False, [1.0640200932648511, 8.801094080240267, \
                        4.303338569406158]],
            'psf_fit_loc': [True, [0, 220.30076542192973]],
            'psf_fit_pix': [False, [44.99810399006913, 44.996587647973605]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/standard_cal.residual.crtf')

        exp_resid_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 5041.0],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0233334042132],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, -0.01152363047],
            'min_val_pos': [True, [50, 16, 0, 0]],
            'im_rms': [False, 0.00410389073592],
            'im_sum': [False, 0.122164499838],
            'regn_sum': [False, 1.0518577156],
            'npts_real': [True, 8100]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/standard_cal.model.crtf')

        exp_model_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 8100.0],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.38273620605],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.0264748467339],
            'im_sum': [False, 2.38273620605],
            'regn_sum': [False, 2.38273620605],
            'mask_non0': [True, 1],
            'npts_real': [True, 8100]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 201522.109375],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 201522.109375],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 201522.108127],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[-0.015, 2.5])
            self.mom8_creator(img+'.residual', range_list=[-0.015, 2.5])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=True)

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

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/mosaic_cube.image.crtf', field_regions = \
            ['circle[[00:45:54.383559, -73.15.29.41306], 22.45arcsec]',
             'circle[[00:45:49.435664, -73.15.35.13742], 22.45arcsec]',
             'circle[[00:45:53.057440, -73.15.50.79016], 22.45arcsec]',
             'circle[[00:45:50.762696, -73.15.13.76177], 22.45arcsec]',
             'circle[[00:45:58.006248, -73.15.45.06040], 22.45arcsec]',
             'circle[[00:45:55.708764, -73.15.08.03543], 22.45arcsec]',
             'circle[[00:45:59.330540, -73.15.23.68133], 22.45arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 8.79758467135],
            'com_bmin': [False, 6.10500627626],
            'com_pa': [False, 64.9303836873],
            'npts': [True, 5925312],
            'npts_unmasked': [True, 3338068.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.18467354774],
            'max_val_pos': [True, [44, 38, 0, 252]],
            'min_val': [False, -0.417252391577],
            'min_val_pos': [True, [24, 31, 0, 498]],
            'im_rms': [False, 0.0878779753708],
            'rms_per_chan': [False, [0.08352891852269921, 0.08135269305379152, 0.08395260930289478, 0.08279304060774786, 0.08596507857262789, 0.08278415992796316, 0.08377357600543471, 0.08394694956313518, 0.08294962386707326, 0.08607043136961237, 0.07680081900457504, 0.09496464679558769, 0.08214392421353917, 0.0857648739112381, 0.08281633996127867, 0.0894806463001662, 0.081589405081739, 0.07924880491236708, 0.07453379448699567, 0.09379627665607704, 0.08157315058137742, 0.09166424964725516, 0.08450533906882468, 0.08106955380637881, 0.08737684423508374, 0.09073897650010625, 0.09084105654348026, 0.07784857687669093, 0.08502491422043724, 0.08973744200475807, 0.08225747757890416, 0.08302380899013398, 0.0770200580373865, 0.08881158355253989, 0.08865534449297038, 0.09159308234217008, 0.08115203022056021, 0.08341851059515137, 0.08987845777172865, 0.08508065238113408, 0.07778307167083158, 0.086558262984569, 0.08338769817916242, 0.08071653499000843, 0.08865177850010854, 0.08534950297770955, 0.07616257958327323, 0.09190545934190404, 0.08895841881770002, 0.0896694261175686, 0.07893346400403167, 0.07455014243908499, 0.09047688550066521, 0.0869112050227715, 0.09350825779882399, 0.07646177066867972, 0.08249014753514353, 0.09101549816292659, 0.08639556025249603, 0.08038246669695737, 0.08824074457422779, 0.08618380562988853, 0.08205108347689954, 0.07771920892371241, 0.08768241197771068, 0.08334623972314524, 0.08815925813978431, 0.08212951592340229, 0.0823396048959934, 0.08234068229313843, 0.08884264540759013, 0.07862615163124824, 0.08655155235861448, 0.08979383171683271, 0.08809736503013246, 0.07192760633210062, 0.08230751478411082, 0.07857712439720876, 0.07948313531326559, 0.09225583910089932, 0.0827786011949029, 0.08972941404057228, 0.0887804253347499, 0.088568302489554, 0.08542522263383671, 0.08705651177174693, 0.09503784914334118, 0.09211196203585155, 0.09124273060482811, 0.0894655349173044, 0.08578196967967916, 0.08600001037669011, 0.08725679761935697, 0.09250336443228556, 0.09141071794136318, 0.08257340474387778, 0.08489860170832744, 0.09095982550662123, 0.08386307951475633, 0.08245963060222646, 0.08082452670747917, 0.08421465818141162, 0.07716197318352062, 0.08491216770880909, 0.09219997011780534, 0.08999904369748837, 0.08917413213129059, 0.08669554914308757, 0.08478598632777916, 0.08272607334416665, 0.09204243226438258, 0.08766378811832916, 0.0865825253580717, 0.07758524805681824, 0.09085091409160705, 0.08521446839560726, 0.0884995697340635, 0.0884733008608976, 0.08577677843356687, 0.09306189230759793, 0.08612785329237876, 0.08258623029884975, 0.07941031533209075, 0.08826196881912282, 0.09840895220275017, 0.08933425191997195, 0.08745304230772397, 0.08566154826880104, 0.09365349891994729, 0.09053316121385116, 0.08810691082194215, 0.08630703760915545, 0.08753186529189085, 0.0937019174346115, 0.0773479036533113, 0.08029497918496241, 0.08724287089313901, 0.08790671029158231, 0.08773912124792194, 0.08266091420346706, 0.08132362905839345, 0.09216634845393681, 0.08762278965147091, 0.09332173318661345, 0.07778672535312933, 0.08408442636695551, 0.08553415947786376, 0.08232058120295321, 0.07780988703074035, 0.08409905782923954, 0.07296439116382543, 0.08925116729628191, 0.09008277668363449, 0.09520991051126101, 0.08553135515808247, 0.084718775948537, 0.08883634984779534, 0.08643599841065604, 0.08339630810303196, 0.08862968769223899, 0.0784883898433635, 0.08099618435890771, 0.08565713860724568, 0.08556682581312276, 0.09128248240593151, 0.08827156303897946, 0.0880790536420939, 0.08578294204187285, 0.08362128063154757, 0.09123831674563974, 0.08722115387496024, 0.09621944328499646, 0.09149116203081945, 0.09118103925268266, 0.08347804034781846, 0.0855692865777019, 0.0897847427742151, 0.07186203727120243, 0.08802972450167627, 0.08442626094668859, 0.08723051889453723, 0.08159640903835322, 0.07548809021347838, 0.0915760655325942, 0.1000199147822218, 0.09787395465748581, 0.08693633490200812, 0.08442264879607819, 0.08721652263357231, 0.0880380635040301, 0.08729925477779406, 0.08519856552892391, 0.08975256932981747, 0.07958208410639706, 0.0914096789556792, 0.09811570777169866, 0.07949196007023566, 0.09006036869795747, 0.0912952536965593, 0.09805249071833734, 0.0897004902132416, 0.08952989036422254, 0.0790849762337038, 0.09347465569139526, 0.09000611959279463, 0.08239101076565238, 0.08373667485591837, 0.0871320704983511, 0.08767396777909775, 0.0965335008897428, 0.08633964662006433, 0.08735933635150735, 0.09035634185538645, 0.08468159828070221, 0.09116890211273068, 0.08804541661464536, 0.08484392727135703, 0.08203979237254262, 0.0889071541447407, 0.0889633730061124, 0.08339973306782386, 0.08330574319001677, 0.08871583444712416, 0.09420962372690272, 0.08521243909929632, 0.08355595121547242, 0.0933236552392962, 0.08719802034432139, 0.09062561297486302, 0.08535346593907196, 0.08813919625097401, 0.092913551212774, 0.09227103971954637, 0.08161432828410638, 0.09072439427447547, 0.09244590715316298, 0.08709339445075073, 0.09924553460801067, 0.1002690868580879, 0.08840557566191516, 0.08133715093336598, 0.09149908227822648, 0.07670073484069043, 0.09249672315752175, 0.0782222639219109, 0.0990254821246598, 0.09280158172995194, 0.08730248696252146, 0.09456558538062404, 0.10081497167443383, 0.10403855973688014, 0.12877330926044459, 0.1386207992219812, 0.13938585037552162, 0.1235720046054832, 0.12067870584152206, 0.10645742986583119, 0.08552529466321122, 0.08268301454232765, 0.08542822040617518, 0.09335660100821405, 0.08885775208502826, 0.08644558884485783, 0.08440970509099167, 0.091904188949323, 0.09288900294073746, 0.09184545798056735, 0.0978475829246572, 0.09536587471324205, 0.09325738676721759, 0.08191614651111201, 0.09130572131132277, 0.09356730058206153, 0.09150948292317468, 0.08930212344785793, 0.08973008889149876, 0.08468876678845778, 0.09439900459047333, 0.08340248011080888, 0.09426625948231673, 0.0846651688286838, 0.08220153827772422, 0.09338524135684238, 0.08541949594391877, 0.07720477952979755, 0.08539185143996587, 0.09204753294875682, 0.08278104837254098, 0.07886155262558238, 0.08341737313272374, 0.08757850027055837, 0.08360325848864149, 0.09166550957557691, 0.09185533846656913, 0.08026617289995454, 0.09106087061249621, 0.08802259852925379, 0.09053259952518417, 0.08328401754435075, 0.08888144950799276, 0.0810576402470788, 0.08479306981872027, 0.08580610257468782, 0.09179491887497933, 0.08841796481950012, 0.0845148944552744, 0.08053481062011857, 0.08862797798354732, 0.09462089152970354, 0.08582787706368114, 0.08814348149204695, 0.08504514772616785, 0.08371224832082137, 0.08678333110045997, 0.09353126562315753, 0.09129854527725832, 0.08116659030415371, 0.09015138809948552, 0.08710081730619218, 0.09437955706280836, 0.08963858464777974, 0.09313459101197238, 0.08856120416820244, 0.08552675606885124, 0.08351176926318361, 0.08411701581891168, 0.08427020020603929, 0.09163946448881417, 0.0805306218916976, 0.08160963806132211, 0.08804552292687956, 0.09626408912773088, 0.08913709670428199, 0.09096834064650154, 0.0851269228240773, 0.09017277009104614, 0.08476290074193281, 0.07632278322336213, 0.08385737538890878, 0.08700039219503956, 0.08866885268530736, 0.08466059774774966, 0.0845759814329557, 0.0790621930867725, 0.08771807117918605, 0.08893473780006282, 0.09310980223541039, 0.09306479580057536, 0.09147420625138089, 0.08949657274281757, 0.08192815006108203, 0.07905463600626796, 0.09666550899459639, 0.0808647476478242, 0.08495044094490133, 0.0916137566838688, 0.09649894883996268, 0.0920414733457368, 0.08356993363476055, 0.09422414928552236, 0.08778457312089456, 0.08987693020950831, 0.09777964977624237, 0.09060848800058789, 0.09087547326073596, 0.09065288590043372, 0.09815595961513002, 0.08630801892602018, 0.08960594520751539, 0.09100452485239247, 0.09259096682095072, 0.09364434916529175, 0.0853051896352123, 0.08391849112594159, 0.08978397560355508, 0.08356817274105656, 0.08639129387321305, 0.07641054760007339, 0.08566749942250246, 0.09110851912136013, 0.08938769699436465, 0.08943111829547572, 0.08605404789230106, 0.08796239383347541, 0.08454535587717901, 0.0929903335078121, 0.08246802178760436, 0.08817700985492735, 0.0820807551149751, 0.08511665149857307, 0.0914822776922275, 0.08779917531640212, 0.0779155145245507, 0.08062449848958794, 0.09151321230973979, 0.08251138633527932, 0.08314391480095229, 0.09660065726688243, 0.09161786920726972, 0.09195890233721427, 0.09484458463309785, 0.08672723967704118, 0.09056091164304954, 0.08950078278038455, 0.08453380213065807, 0.08621663260330842, 0.0903504892183676, 0.08888219947767556, 0.09691781310403562, 0.0829510997618925, 0.08538905047156051, 0.08536033187251706, 0.09253646586834798, 0.08719827400559628, 0.08741965478896539, 0.0908875936865952, 0.08650855583109175, 0.0911287851112432, 0.0870327992529023, 0.09334187790207615, 0.08691128023529447, 0.0829607319773311, 0.08561452123819384, 0.09416699659894374, 0.09865975100963004, 0.08059372543252491, 0.08162290581093083, 0.07969201254174872, 0.09014664917727015, 0.07748434709443736, 0.09115822285795372, 0.0874199097979386, 0.08331094704142918, 0.08450373137759747, 0.0873987436304484, 0.0792090383862253, 0.08682449919890456, 0.08898017363528224, 0.0891014307981212, 0.08578455417679541, 0.09612343808420232, 0.07718957370637691, 0.08963596680253917, 0.09053358289784386, 0.08104077369182848, 0.08805192318672665, 0.09036158074767288, 0.09733534898340712, 0.08642234534217835, 0.0873180423896324, 0.08509746809331085, 0.0927045997121519, 0.08985111493399206, 0.09486348772674255, 0.09000267315635847, 0.08474935185150062, 0.08247809724159613, 0.08802043461258492, 0.0865126082695781, 0.08138258449164787, 0.08795575893476618, 0.09456070784224847, 0.09680636657810467, 0.08476700362040308, 0.08670105726809373, 0.08636724547506389, 0.08412716463066074, 0.08428810620773179, 0.08964151944607554, 0.08689624513108678, 0.08902965594822292, 0.09221339424501332, 0.08726043474359556, 0.08607641544478577, 0.09100554185059015, 0.08492870794208009, 0.08529837352577493, 0.09521158569562824, 0.08914943856267118, 0.087555731639101, 0.0862048336688618, 0.08984721078423315, 0.08217617292259297, 0.08824966695006062, 0.07486261467473272, 0.08753387468056906, 0.08379545796577004, 0.09274777146757962, 0.09220642715156253, 0.0792962124445207, 0.09090807566463247, 0.08751737767113807, 0.07961706129268199, 0.0941224640615722, 0.0795895706794336, 0.09562104758697967, 0.08020847225726233, 0.09417989583892716, 0.09061167014269772, 0.08898710965217106, 0.0897447948736654, 0.08398102899483291, 0.08684215184345169, 0.096630024031122, 0.08473098919932259, 0.09179580145778438, 0.07887094371255345, 0.08638286938225163]],
            'im_sum': [False, 365.155811516],
            'regn_sum': [False, 91.3382614758],
            'npts_real': [True, 5925312],
            'rms_per_field': [False, [0.089055932035104174, 0.087593317854421537, 0.088754854977895078, 0.087673584891401882, 0.088859674046269987, 0.087793556001629414, 0.087858029812576927]],
            'profile': [False, 1.2198751],
            'im_fit': [False, [1.150589250979386, 9.433930899157085, \
                       8.092281808705982]],
            'im_fit_loc': [True, [252, 220.3142062699946]],
            'im_fit_pix': [False, [44.567451454649344, 38.38608598796125]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 5925312],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'mask_pix': [True, 3928],
            'mask_regns': [True, 1],
            'npts_real': [True, 5925312]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/mosaic_cube.pb.crtf')

        exp_pb_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 3338068.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, 0.200100466609],
            'min_val_pos': [True, [43, 99, 0, 447]],
            'im_rms': [False, 0.615051728939],
            'npts_0.2': [True, 3338068],
            'npts_0.5': [True, 1825336],
            'npts_real': [True, 5925312],
            'pb_fit': [False, [1.081101854419883, 69.22423084266158, \
                       69.16483065383713]],
            'pb_fit_loc': [True, [252, 220.3142062699946]],
            'pb_fit_pix': [False, [54.07005299848095, 54.01350317818875]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/mosaic_cube.psf.crtf')

        exp_psf_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, -0.168374374509],
            'min_val_pos': [True, [63, 54, 0, 12]],
            'im_rms': [False, 0.0604965318073],
            'im_sum': [False, 66.9960417559],
            'npts_real': [True, 5925312],
            'psf_fit': [False, [1.1012741250775981, 7.879923449290915, \
                        5.24248657397286]],
            'psf_fit_loc': [True, [252, 220.3142062699946]],
            'psf_fit_pix': [False, [53.987524754275434, 54.00459088886782]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/mosaic_cube.residual.crtf')

        exp_resid_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 3338068.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.490056365728],
            'max_val_pos': [True, [49, 60, 0, 249]],
            'min_val': [False, -0.417252391577],
            'min_val_pos': [True, [24, 31, 0, 498]],
            'im_rms': [False, 0.0875338790102],
            'im_sum': [False, 68.0894097835],
            'regn_sum': [False, 37.8514072159],
            'npts_real': [True, 5925312]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/mosaic_cube.model.crtf')

        exp_model_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.510233163834],
            'max_val_pos': [True, [55, 57, 0, 255]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000534801027393],
            'im_sum': [False, 5.90552844107],
            'regn_sum': [False, 1.13609256595],
            'mask_non0': [True, 32],
            'npts_real': [True, 5925312]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 120.900192261],
            'max_val_pos': [True, [0, 0, 0, 447]],
            'min_val': [False, 120.665283203],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 120.761758344],
            'npts_real': [True, 508]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight')

        exp_wt_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [True, 244174.08728027344],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.39385086298],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, 7.44352073525e-05],
            'min_val_pos': [True, [97, 107, 0, 224]],
            'im_rms': [False, 0.140916352052],
            'im_sum': [False, 506828.976527],
            'npts_0.2': [True, 1058348],
            'npts_0.3': [True, 504298],
            'npts_real': [True, 5925312]}

        report9 = self.stats_compare(exp_wt_stats, wt_stats_dict, '.weight')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9


        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[0.15, 1.2])
            self.mom8_creator(img+'.residual', range_list=[0.15, 1.2])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[testname]['images'].append(img+'.image.profile.png')

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
        img = os.getcwd()+'/'+file_name+'1'
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
            deconvolver='hogbom', restoration=False, restoringbeam='common',\
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=True)

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
            parallel=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/mosaic_mfs.image.crtf', field_regions = \
            ['circle[[13:56:07.210000, +05.15.17.20000], 45.9arcsec]',
             'circle[[13:56:06.355525, +05.15.59.74129], 45.9arcsec]',
             'circle[[13:56:04.316267, +05.15.27.41716], 45.9arcsec]',
             'circle[[13:56:09.249291, +05.15.49.52355], 45.9arcsec]',
             'circle[[13:56:05.170768, +05.14.44.87604], 45.9arcsec]',
             'circle[[13:56:10.103706, +05.15.06.98201], 45.9arcsec]',
             'circle[[13:56:08.064442, +05.14.34.65864], 45.9arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 17.6737785339],
            'com_bmin': [False, 10.060172081],
            'com_pa': [False, 86.6785964966],
            'npts': [True, 15876],
            'npts_unmasked': [True, 8454.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0345157124102],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, -0.0019364656182],
            'min_val_pos': [True, [91, 52, 0, 0]],
            'im_rms': [False, 0.00202865566602],
            'im_sum': [False, 1.51295999996],
            'regn_sum': [False, 1.58850853646],
            'npts_real': [True, 15876],
            'rms_per_field': [False, [0.0043183333329864697, 0.0035324567542234214, 0.0033643881411162453, 0.0034484378708067886, 0.0035013779149007367, 0.0033359473616749375, 0.0035589233954835845]],
            'im_fit': [False, [0.03522500582263719, 17.46093579518058, 
                       9.709830310449933]],
            'im_fit_loc': [True, [0, 107.8402451422565]],
            'im_fit_pix': [False, [62.9942562846151, 62.995885097033394]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [64, 114, 0, 0]), \
                      (img+'.image', False, [64, 115, 0, 0]), \
                      (img+'.image', True, [11, 60, 0, 0]), \
                      (img+'.image', False, [10, 60, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 15876],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 360],
            'mask_regns': [True, 1],
            'npts_real': [True, 15876]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/mosaic_mfs.pb.crtf')

        exp_pb_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 8454.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 0.200080364943],
            'min_val_pos': [True, [102, 28, 0, 0]],
            'im_rms': [False, 0.604913988836],
            'npts_0.2': [True, 8454],
            'npts_0.5': [True, 4497],
            'npts_real': [True, 15876],
            'pb_fit': [False, [1.0693559655652305, 141.80580479462876, \
                       141.74549135472637]],
            'pb_fit_loc': [True, [0, 107.8402451422565]],
            'pb_fit_pix': [False, [62.975154097364715, 62.94725116661756]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/mosaic_mfs.psf.crtf')

        exp_psf_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, -0.169576942921],
            'min_val_pos': [True, [61, 57, 0, 0]],
            'im_rms': [False, 0.0501544145677],
            'im_sum': [False, 0.00266324029255],
            'npts_real': [True, 15876],
            'psf_fit': [False, [1.088207720785799, 15.893701850875548, \
                        8.795192549423799]],
            'psf_fit_loc': [True, [0, 107.8402451422565]],
            'psf_fit_pix': [False, [62.98416058527938, 63.00086190688355]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, \
            '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/mosaic_mfs.residual.crtf')

        exp_resid_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 8454.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.00713972467929],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, -0.00197902671061],
            'min_val_pos': [True, [72, 58, 0, 0]],
            'im_rms': [False, 0.000868739833233],
            'im_sum': [False, 0.134135401952],
            'regn_sum': [False, 0.291504465893],
            'npts_real': [True, 15876]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/mosaic_mfs.model.crtf')

        exp_model_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0273759867996],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000217269736505],
            'im_sum': [False, 0.0273759867996],
            'regn_sum': [False, 0.0273759867996],
            'mask_non0': [True, 1],
            'npts_real': [True, 15876]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 4396210.5],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 4396210.5],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 4396210.53446],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight')

        exp_wt_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.426215469837],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 4.20771575591e-05],
            'min_val_pos': [True, [1, 24, 0, 0]],
            'im_rms': [False, 0.144357241275],
            'im_sum': [False, 1347.6264633],
            'npts_0.2': [True, 2778],
            'npts_0.3': [True, 1468],
            'npts_real': [True, 15876]}

        report9 = self.stats_compare(exp_wt_stats, wt_stats_dict, '.weight')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9


        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[-0.002, 0.035])
            self.mom8_creator(img+'.residual', range_list=[-0.002, 0.035])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=True)

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
            parallel=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mos_mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.image.tt0.crtf', \
            field_regions = \
            ['circle[[13:56:07.210000, +05.15.17.20000], 45.9arcsec]',
             'circle[[13:56:06.355525, +05.15.59.74129], 45.9arcsec]',
             'circle[[13:56:04.316267, +05.15.27.41716], 45.9arcsec]',
             'circle[[13:56:09.249291, +05.15.49.52355], 45.9arcsec]',
             'circle[[13:56:05.170768, +05.14.44.87604], 45.9arcsec]',
             'circle[[13:56:10.103706, +05.15.06.98201], 45.9arcsec]',
             'circle[[13:56:08.064442, +05.14.34.65864], 45.9arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 17.6737785339],
            'com_bmin': [False, 10.060172081],
            'com_pa': [False, 86.6785964966],
            'npts': [True, 15876],
            'npts_unmasked': [True, 8454.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.034914586693],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, -0.00258040963672],
            'min_val_pos': [True, [70, 97, 0, 0]],
            'im_rms': [False, 0.00207815358165],
            'im_sum': [False, 1.47797251044],
            'regn_sum': [False, 1.58228409861],
            'npts_real': [True, 15876],
            'rms_per_field': [False, [0.0043714750541124468, 0.003595955649075896, 0.003430812546940405, 0.0034721618952933994, 0.0035884838221589368, 0.0033594003986614097, 0.0036156061407187608]],
            'im_fit': [False, [0.03580520672441999, 17.19187684101627, \
                       9.68274896612347]],
            'im_fit_loc': [True, [0, 107.8402451422565]],
            'im_fit_pix': [False, [63.09049673358014, 62.94805812018937]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [64, 114, 0, 0]), \
                      (img+'.image.tt0', False, [64, 115, 0, 0]), \
                      (img+'.image.tt0', True, [11, 60, 0, 0]), \
                      (img+'.image.tt0', False, [10, 60, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image.tt0')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 15876],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 360],
            'mask_regns': [True, 1],
            'npts_real': [True, 15876]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.pb.tt0.crtf')

        exp_pb_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 8454.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 0.200080364943],
            'min_val_pos': [True, [102, 28, 0, 0]],
            'im_rms': [False, 0.604913988836],
            'npts_0.2': [True, 8454],
            'npts_0.5': [True, 4497],
            'npts_real': [True, 15876],
            'pb_fit': [False, [1.0693559655651996, 141.80580479464936, \
                       141.74549135470988]],
            'pb_fit_loc': [True, [0, 107.8402451422565]],
            'pb_fit_pix': [False, [62.975154097364715, 62.94725116661756]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb.tt0')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.psf.tt0.crtf')

        exp_psf_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, -0.169576942921],
            'min_val_pos': [True, [61, 57, 0, 0]],
            'im_rms': [False, 0.0501544145677],
            'im_sum': [False, 0.00266324029255],
            'npts_real': [True, 15876],
            'psf_fit': [False, [1.0781857293103545, 15.898196388608632, \
                        8.995969894587292]],
            'psf_fit_loc': [True, [0, 107.8402451422565]],
            'psf_fit_pix': [False, [62.991298508308404, 63.00339664380328]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf.tt0')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            region_file = data_path+
            'region_files/mosaic_mtmfs.residual.tt0.crtf')

        exp_resid_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 8454.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.00713019119576],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, -0.00197783415206],
            'min_val_pos': [True, [72, 58, 0, 0]],
            'im_rms': [False, 0.000865575345209],
            'im_sum': [False, 0.130538249029],
            'regn_sum': [False, 0.283332834988],
            'npts_real': [True, 15876]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0')

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs.model.tt0.crtf')

        exp_model_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0272394213825],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000216185883988],
            'im_sum': [False, 0.0272394213825],
            'regn_sum': [False, 0.0272394213825],
            'mask_non0': [True, 1],
            'npts_real': [True, 15876]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model.tt0')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 4396210.5],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 4396210.5],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 4396210.53446],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0')

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight.tt0')

        exp_wt_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.426215469837],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 4.20771575591e-05],
            'min_val_pos': [True, [1, 24, 0, 0]],
            'im_rms': [False, 0.144357241275],
            'im_sum': [False, 1347.6264633],
            'npts_0.2': [True, 2778],
            'npts_0.3': [True, 1468],
            'npts_real': [True, 15876]}

        report9 = self.stats_compare(exp_wt_stats, wt_stats_dict, \
            '.weight.tt0')

        # .alpha report
        alpha_stats_dict = self.image_stats(img+'.alpha', region_file = \
            data_path+'region_files/mosaic_mtmfs.alpha.crtf')

        exp_alpha_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 1628.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 39.2035827637],
            'max_val_pos': [True, [52, 71, 0, 0]],
            'min_val': [False, -34.6691818237],
            'min_val_pos': [True, [82, 85, 0, 0]],
            'im_rms': [False, 11.2641311632],
            'im_sum': [False, 8782.46312377],
            'regn_sum': [False, -46.7235065232],
            'mask_non0': [True, 1628],
            'npts_real': [True, 15876]}

        report10 = self.stats_compare(exp_alpha_stats, alpha_stats_dict, \
            '.alpha')

        # .alpha.error report
        error_stats_dict = self.image_stats(img+'.alpha.error', region_file = \
            data_path+'region_files/mosaic_mtmfs.alpha.error.crtf')

        exp_error_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 1628.0],
            'freq_bin': [True, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 55.4204788208],
            'max_val_pos': [True, [52, 71, 0, 0]],
            'min_val': [False, 0.0116222631186],
            'min_val_pos': [True, [15, 49, 0, 0]],
            'im_rms': [False, 15.8834614404],
            'im_sum': [False, 19545.8077074],
            'regn_sum': [False, 289.468003012],
            'npts_real': [True, 15876]}

        report11 = self.stats_compare(exp_error_stats, error_stats_dict, \
            '.alpha.error')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11


        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image.tt0', range_list=[-0.003, 0.035])
            self.mom8_creator(img+'.residual.tt0', range_list=[-0.003, 0.035])
            test_dict[testname]['images'].extend( \
                (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=False)

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
            parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/mosaic_cube_eph.image.crtf', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 0.930861830403],
            'com_bmin': [False, 0.719328746863],
            'com_pa': [False, -88.1227870446],
            'npts': [True, 191116800],
            'npts_unmasked': [True, 105033551.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.0847059637308],
            'max_val_pos': [True, [290, 236, 0, 511]],
            'min_val': [False, -0.0710606276989],
            'min_val_pos': [True, [192, 206, 0, 674]],
            'im_rms': [False, 0.0122806235373],
            'rms_per_chan': [False, [0.015667012104297096, 0.01249556647636659, 0.011119753312996314, 0.011559686287211745, 0.012295549107925172, 0.013497595553035495, 0.013890985018415562, 0.01537521583867029, 0.012560026473467654, 0.01211174756347237, 0.01068202599475987, 0.011446394296400014, 0.013360170834984883, 0.011779590786641968, 0.011813658887463094, 0.01105808607636477, 0.011211227704164385, 0.012672658897254005, 0.009644489317411001, 0.012306186479311533, 0.011820657820002498, 0.012528421200108672, 0.012097912099003863, 0.014149728955396576, 0.013040286513860739, 0.010926749613896915, 0.01135630004195933, 0.011450114870829052, 0.013874579045960853, 0.012713343335685372, 0.012365580200560802, 0.01112892658338281, 0.011735154726597849, 0.010271529404584005, 0.012034388520273378, 0.011825840449510846, 0.01201477612491181, 0.011271561871280695, 0.01263051888325016, 0.010917661123920426, 0.012304591002675239, 0.01382691832578588, 0.012025131483665455, 0.01471038786901529, 0.01268901239406491, 0.012269947721949661, 0.010339439002704098, 0.010588455121686411, 0.010663853975518148, 0.011927246586817075, 0.011971607513601665, 0.011740843559651189, 0.010860788064116878, 0.0113500230752135, 0.011498619788115836, 0.011798959247821511, 0.010929419432007908, 0.014579761324742322, 0.012976422990339583, 0.010218616946376293, 0.011098072571178239, 0.01168771322025399, 0.011600913193876674, 0.011684300942923293, 0.01071321921174995, 0.011492474041242522, 0.012011160150216459, 0.01290685298415816, 0.010431320458890848, 0.010472939825276044, 0.013447112178611359, 0.011687328995254547, 0.010898530945791702, 0.009752776379182097, 0.01164870192026142, 0.012676083477754277, 0.01431102043820875, 0.014295093842479985, 0.012626217708334785, 0.011231466919543581, 0.010535901939971674, 0.01062302128541249, 0.011298369163697381, 0.012502624763702992, 0.011305836954063448, 0.01215581953177045, 0.010607739795551263, 0.01112311151599714, 0.010136484643153825, 0.01128990114770225, 0.012014173589023637, 0.013193805876756995, 0.012330109734653474, 0.010766011891193288, 0.01153437959726003, 0.011978851287972363, 0.011109003021191307, 0.011833746981365703, 0.011880050178226903, 0.009862018971003887, 0.012800148141077113, 0.012007057726628611, 0.012985090676833535, 0.0129505958602198, 0.011696207783179233, 0.013661835268823226, 0.015246672300909788, 0.012047095989451426, 0.014346897791128894, 0.011754599415304922, 0.011009518999420173, 0.011052019581176067, 0.012820563697322922, 0.01167880516927239, 0.010768042038144418, 0.012407448804067968, 0.011258173984772548, 0.010006383872394846, 0.012053714508583965, 0.011282069296445192, 0.01185568667180287, 0.01224398735628035, 0.009954101011210293, 0.01224524095389172, 0.01083068385084185, 0.012042693023557548, 0.010466482413183785, 0.01270824749297588, 0.009858469333879292, 0.009616370018819147, 0.011251488506345967, 0.010944010326885257, 0.009937336043765989, 0.011811563673806756, 0.01374671433132631, 0.011885072616527947, 0.01128936295560062, 0.010745268810504142, 0.010783213461157924, 0.010401671249270496, 0.010546985468598878, 0.011704982925431526, 0.01151082357877418, 0.011389394982234517, 0.009891395767145527, 0.012676962217854196, 0.01239854403608683, 0.01012692292366004, 0.01208213727619496, 0.010426410721875919, 0.011489717695650697, 0.011351311722094469, 0.010788345620860689, 0.012413824258742542, 0.012869284594202364, 0.01254719377424318, 0.012493275011557376, 0.012026588643842862, 0.011444080682226833, 0.014504798148840915, 0.012783237416754493, 0.010978226098473022, 0.011226830327894566, 0.010386044707730815, 0.01179048796956224, 0.010466753647787678, 0.010520561042133916, 0.009806686930384046, 0.009825709113667067, 0.011163605933804482, 0.009898306765781937, 0.010639624156931048, 0.010765250926759881, 0.011871167097798471, 0.012699735217972713, 0.012271717028350141, 0.01011927965625519, 0.01185201591720369, 0.01187934596239281, 0.012860096441123013, 0.01045831490108293, 0.01168302685217849, 0.01231538733094323, 0.013711623609843774, 0.010509256728676152, 0.009725104065778598, 0.010742348970319717, 0.010497654900964383, 0.01226776948490034, 0.009968096362512743, 0.01065210210757194, 0.0113790694870219, 0.013691499504725839, 0.009709697345128918, 0.011668556206361638, 0.01203851403100267, 0.011090927358203469, 0.012322941164683295, 0.010763675962361283, 0.01243680774921738, 0.01273630413642944, 0.014925176287167791, 0.013343151602031603, 0.012184525933592252, 0.011868284483931634, 0.010612397890136424, 0.010925259212727632, 0.010736722739232834, 0.01128555748033192, 0.012106425521907454, 0.012536567241883659, 0.011811387511674406, 0.011410855057029651, 0.010610454803475318, 0.010985570534714531, 0.013420524115250722, 0.01237714866693324, 0.011400339204092964, 0.011530006934858634, 0.013404558644802023, 0.01356996738906636, 0.011809592810750743, 0.012195836458918021, 0.013983595304231499, 0.012723543201613906, 0.011350835725055812, 0.012149584287404331, 0.011272979457253907, 0.012262167105146599, 0.011175511659374244, 0.0102697228260645, 0.01334158429888008, 0.012408056314862907, 0.011491512834255597, 0.00928563088957559, 0.011723585674628282, 0.012562330134652718, 0.012321060378936478, 0.011175930032929103, 0.014444827741007454, 0.013240412166233477, 0.00954358182816961, 0.009886456935597223, 0.011315336508350314, 0.011160043368747817, 0.011084649991455213, 0.011426321774805799, 0.010549594186710465, 0.011211051291167739, 0.01176768646746374, 0.014686912922251879, 0.011246937183158924, 0.009773020704868391, 0.012859597369260864, 0.012513039230848268, 0.013313363039583453, 0.012077863600292971, 0.011426654372094657, 0.011508622055300996, 0.011247566432148776, 0.011623431851054248, 0.011014612057401632, 0.013022506823681756, 0.011485416194269382, 0.01154829448377298, 0.01202691068493326, 0.012642809433009043, 0.010836720661823337, 0.012182321211669447, 0.012349974771377626, 0.01129390722457532, 0.011135976780555478, 0.011531711369332466, 0.012456469996333023, 0.01066628966176333, 0.011106915433366485, 0.012535628404137928, 0.01329621583757732, 0.01309891488346235, 0.011823826809276102, 0.010766426258224856, 0.010769239154148987, 0.0110549819438351, 0.00941967982610333, 0.012016650037610828, 0.011327327654288295, 0.01209848828634942, 0.011506245997467066, 0.011684994479376746, 0.013606364448973455, 0.010446608584558138, 0.01066351605331571, 0.010707363943853028, 0.011928301857376082, 0.01294521124420252, 0.010131421709150719, 0.011849062523828911, 0.011443882324740128, 0.012461616479499914, 0.011076627457455922, 0.012813631803261182, 0.011813942249706206, 0.012869872838703512, 0.016215007268429744, 0.016455835790109574, 0.016195858609782746, 0.010586943752781703, 0.011714040295832838, 0.011002274392087834, 0.013100914878598957, 0.011252873658270539, 0.011303009258586938, 0.012859801115861635, 0.012461620893746536, 0.011357224642469422, 0.011637071400518138, 0.011468136930703879, 0.01411435236184562, 0.01213930177391437, 0.012012513247874964, 0.012399215107416643, 0.017256708026613384, 0.011747981333506227, 0.012212040710131029, 0.011702161891262641, 0.01357489047444866, 0.0142883693783086, 0.013345754616212652, 0.0111940675739781, 0.011675373377606931, 0.011441040364084332, 0.013384698628456768, 0.013142898529307593, 0.01225725868583698, 0.012226995314265549, 0.011976385838259086, 0.014630107009854095, 0.01393710853346633, 0.011607705747133181, 0.013218275149292351, 0.011873190322311846, 0.013747319388144195, 0.010499983139191808, 0.011624206615490757, 0.013225475425745141, 0.013680190321206091, 0.012046704878481942, 0.012471069574344537, 0.010953227174093479, 0.011926261656122828, 0.012459493726885675, 0.014416014101405097, 0.013092512809823969, 0.012391343965965609, 0.012032749220004132, 0.011174480504353115, 0.0128020365473887, 0.012373815528334766, 0.013220253396707947, 0.011397076938533486, 0.011737376857836166, 0.011339198156054645, 0.013683087656476136, 0.014055337259697667, 0.011872209955622893, 0.011489327319776421, 0.011356744920908554, 0.011924222252823101, 0.012924011814257761, 0.01156472654485118, 0.012452508623193576, 0.01615067926057454, 0.014393758763013523, 0.012021023556503968, 0.011686927548314128, 0.012558578572353821, 0.011698219106403014, 0.011850872070781452, 0.011244104935548218, 0.010803851710396785, 0.011326571988274907, 0.012046055166178125, 0.012789254887232475, 0.012054328352918125, 0.012531458975574155, 0.011667295333375608, 0.012715182062018394, 0.01416008745695674, 0.011590350023313635, 0.013280049474709302, 0.012317202369770713, 0.011350237164077244, 0.011468547414280923, 0.011839276627787746, 0.011984505075773167, 0.010498556329113671, 0.011607479477631402, 0.01146628697832332, 0.012454745876293334, 0.015324136554616518, 0.012648464631195741, 0.011476812708577248, 0.012388293945433912, 0.013097915532805426, 0.010844862222727058, 0.010849704992030628, 0.011857198994433613, 0.014369812497106239, 0.01614444176168378, 0.010859408246343039, 0.012355198379138267, 0.013369345496468774, 0.01352324109813424, 0.01249093576911017, 0.011839983100254041, 0.013790424278145797, 0.013615888137811373, 0.012137597782120882, 0.011920659584678292, 0.01198453326792445, 0.011455625248390569, 0.00998759798749584, 0.011907441379145489, 0.011445931803455656, 0.012802399352945356, 0.011807120236149467, 0.01226251626371192, 0.013803196660688412, 0.011396570477306092, 0.011871861677445866, 0.012409145934442605, 0.012190261729970274, 0.011355705245344998, 0.013059478167556213, 0.011340084762450894, 0.014094393707136091, 0.012227696710476979, 0.015109377241227427, 0.012328532759057612, 0.012948495238519125, 0.012448990941930516, 0.014909025593587475, 0.013966538678131317, 0.013951304003408212, 0.012643902341794237, 0.012381886725880696, 0.011608607970324055, 0.011524217824005454, 0.011753643583578194, 0.011909335216658386, 0.01500036992819057, 0.0146564784502679, 0.013039737284216963, 0.013751934845130839, 0.013710006202952523, 0.01146836988848032, 0.011545597482743033, 0.011534520901890178, 0.013313753578876611, 0.012923196186255153, 0.013864645252987274, 0.01291016687946222, 0.012187645424073405, 0.0123954851702774, 0.011090370304889066, 0.012569939183393495, 0.011171977695060169, 0.01193123241184618, 0.012223023550711422, 0.012575288535573641, 0.013397233030543118, 0.012677733376390488, 0.012824707036730977, 0.012077443785414215, 0.011917723071670982, 0.012440011533195431, 0.012737137257114249, 0.013226824784925779, 0.011705373001984695, 0.013236156564942205, 0.012428547205061067, 0.011897634394487568, 0.012340276631541067, 0.016514178673390664, 0.01279768810975663, 0.013295330804538989, 0.013923981704086203, 0.013293372296425762, 0.012807871537395961, 0.01233101894981078, 0.012328987872006852, 0.01163784902034462, 0.013313092897728503, 0.012985108364251114, 0.011142998656392029, 0.012369903540412508, 0.012947299676825685, 0.012653423845945748, 0.01562996207353602, 0.014022980315620364, 0.013079373939625572, 0.013669501711869474, 0.011451239023725265, 0.01387333519961978, 0.013272974615061776, 0.013010075621525907, 0.012057515236229667, 0.013348290271015285, 0.012474986304014272, 0.012822235425304088, 0.013469504747600882, 0.015049934910563677, 0.013108662650721938, 0.01169414937155411, 0.014209631933301799, 0.014029595420562672, 0.01313572956489975, 0.013768070645688276, 0.013669663522466149, 0.012789385397857027, 0.011759854640416133, 0.013381506917531618, 0.01317316470287087, 0.013731014632491167, 0.014096287340039539, 0.011977139109191822, 0.01269262588325782, 0.014216384246730137, 0.012780653296103527, 0.012899528621024064, 0.011032967652830057, 0.012102612001078233, 0.011976078228056453, 0.011883895621375998, 0.012802365574674025, 0.01324801385166456, 0.013229524374842769, 0.012151160553514988, 0.011697133232144788, 0.013819175149287335, 0.012140935480252437, 0.011981813853002501, 0.012084722992447795, 0.012585444087037377, 0.012414152508510147, 0.013316552564554114, 0.013253620957310662, 0.014878797637568952, 0.013615894050495254, 0.012709069629370266, 0.012427103871389584, 0.013255141575299076, 0.012738711652745057, 0.014904432456876815, 0.014210867966534892, 0.013509920011482207, 0.012425165695341685, 0.012533528802502332, 0.014336868162538341, 0.011833811236181477, 0.015193744677447416, 0.01431279129775192, 0.013259559570353551, 0.012994100393289547, 0.013864774950178729, 0.012159929964791984, 0.01419863835690063, 0.014190789497478196, 0.013672387203935446, 0.012055278077236548, 0.012085479493499343, 0.012708973627136714, 0.012057165753153366, 0.012948753136076728, 0.01300072865900603, 0.011737981749916887, 0.013527304652929072, 0.012929913822303431, 0.011886809618695744, 0.012833849789732805, 0.01302949373163326, 0.013895403620873681, 0.013171820320213534, 0.013953961887679132, 0.013734000181833348, 0.01262587625087098, 0.012199103159950037, 0.013643673456292997, 0.013397836614489634, 0.01594172388898437, 0.015573035287868315, 0.013917856223549286, 0.012536353641427949, 0.014860466731552354, 0.013552708900045911, 0.013629466127257123, 0.012378335224334382, 0.013116351705612704, 0.012177988274798253, 0.012805373309334908, 0.013766331009848752, 0.014681878707956165, 0.011926702054932887, 0.012452273266241426, 0.012579633618734222, 0.012271260658847647, 0.013434178141045432, 0.011256170201273603, 0.012241164545083469, 0.012715311381559832, 0.014158217761420192, 0.013710431625047755, 0.015544655166662385, 0.012479494691650066, 0.012794382968093331, 0.013723467058383875, 0.01129790706237436, 0.01386400204442386, 0.013948937970667924, 0.01474167932185738, 0.012394940495298034, 0.014413979480373103, 0.013415306676650376, 0.013015296980135697, 0.012460061661220363, 0.012148683531904852, 0.014565058530881577, 0.013354317640191604, 0.010932991047176424, 0.013472405580926574, 0.014090167896930112, 0.011289132690347727, 0.013030955001855788, 0.011930998072436357, 0.01119650545759846, 0.01336252769201265, 0.013985875731024684, 0.011744402770443252, 0.012987363600946232, 0.012570373512332105, 0.011967267431977555, 0.013963824316748901, 0.013113917131760142, 0.015175164104330932, 0.01327121004907913, 0.01480994279997818, 0.013017833326488197, 0.012590284721083084, 0.013636613698566152, 0.013023427264598512, 0.0118689366333651, 0.014431815130736951, 0.011510423632609004, 0.013751161061425392, 0.011946926261948976, 0.011574213816363992, 0.013103519491424586, 0.01256833652316557, 0.013717201694781961, 0.011799548188880763, 0.014221219755103955, 0.013122869137854539, 0.015681024673163967, 0.013977036504021368, 0.01385158873128627, 0.01192090212998216, 0.013405064782284188, 0.012370543487842847, 0.013543641566948098, 0.012190038354089203, 0.010485951880655316, 0.01065474654892235, 0.011894241281652475, 0.012235290499779875, 0.016341337277583624, 0.011951837129874415, 0.012786343847113055, 0.011700169381836972, 0.0116987237015556, 0.01584129766486376, 0.012934100451129409, 0.012421043404934572, 0.013319782161544239, 0.010822586447971657, 0.013768225063031962, 0.01256900129001914, 0.013155623546818948, 0.011216305081099302, 0.01342997193575067, 0.01202191464374732, 0.010697624321778235, 0.015231405911477284, 0.012707984996807302, 0.011509494099144266, 0.012370090248740612, 0.012333106398064926, 0.012355927978977287, 0.013324826178046834, 0.012937849744614949, 0.013251091849170686, 0.012241393042655245, 0.012809410411935642, 0.013415437736483378, 0.013549261827461793, 0.013063192366571214, 0.012121670719955952, 0.01336030963204456, 0.012030706877818232, 0.01620522326989027, 0.013898156174912801, 0.013845480988651687, 0.013887455990710647, 0.012068390405162546, 0.011755939925068324, 0.013275604498608711, 0.011216611629757777, 0.011677478386439837, 0.010605866598758461, 0.0123599433807257, 0.013049424843940767, 0.011751712640564415, 0.011750132281128709, 0.012538068443295915, 0.01178746192810137, 0.012416680391644138, 0.013967663494309064, 0.011784356183646988, 0.014604214590101753, 0.011833948152033797, 0.011315804095294313, 0.01332743466825628, 0.013053434708846068, 0.011059804195895667, 0.010564587906126596, 0.01102267731863153, 0.011342558989477891, 0.01202689018525376, 0.012453421183126498, 0.011102185306584089, 0.010463866575165995, 0.015269163260545289, 0.011918415549432885, 0.011369128293761672, 0.010371364768147193, 0.011736371996752575, 0.012797487595901094, 0.013259851973223793, 0.012380797309777503, 0.013308425109731127, 0.015320758463204735, 0.010225926749762294, 0.011258787297291056, 0.009966991432626141, 0.011479673309463388, 0.012747496824004115, 0.012375931475739925, 0.011058517478639826, 0.012762924931169286, 0.012477944671249878, 0.01260839664096796, 0.011204724607765226, 0.013902208290096557, 0.010776228597758947, 0.013309499022509929, 0.017424540106428448, 0.011666803285392461, 0.010612250155377485, 0.01225365093143365, 0.012565477652135646, 0.012560022064203515, 0.013311882508317791, 0.012505545505753059, 0.010306455499657464, 0.012676663075725451, 0.012498771477844908, 0.01058998021853162, 0.01405727829447646, 0.011586866156371556, 0.010518281562552243, 0.010731082611461988, 0.013668577521770654, 0.012951237133989603, 0.013297250697740166, 0.011256633699887369, 0.011410614744018555, 0.011920105512413305, 0.010903972882515475, 0.01151646075083865, 0.011652106943818302, 0.011279191661803575, 0.013096542735610034, 0.010571055541644716, 0.0125850059607434, 0.010687465761303071, 0.012287345913328777, 0.011395341842353034, 0.01066365582500683, 0.014598748750369274, 0.012549443442399444, 0.011441826793593216, 0.010884568047403273, 0.01180816627630304, 0.01120499414253317, 0.012397216994592853, 0.011586660529077834, 0.012109595113175883, 0.01012537646617585, 0.010614883001281113, 0.013062279023311165, 0.011156508307384487, 0.011360974428249597, 0.010662569680858604, 0.011096498962291788, 0.010554891464213529, 0.011990213710704774, 0.01199714967099034, 0.011854981439770644, 0.011862683712584867, 0.012802582494057203, 0.011577453698377574, 0.01019052971268619, 0.01019901528332309, 0.012281126786588159, 0.011415812682763855, 0.01105581294331331, 0.01066694467332525, 0.010788622197805227, 0.012827960070295203, 0.014827747643843274, 0.01228615008936392, 0.012695947954523156, 0.01372504743362642, 0.012856030967871494, 0.010597007817762974, 0.012304902393004517, 0.012065110220915125, 0.01195668425986831, 0.01086087593260838, 0.01231769839828131, 0.012128147369877246, 0.010817010237205124, 0.012643034637853885, 0.011662533330489484, 0.01332388398778439, 0.011379239016704942, 0.011338092027712373, 0.013532638696355816, 0.01088929059364428, 0.009833061982280674, 0.01039296581185181, 0.010551371059577607, 0.012800002623802556, 0.011294950650386104, 0.010968921330141941, 0.011815483268140185, 0.010991027109824806, 0.010101876381734976, 0.011413636304184781, 0.01168187248586187, 0.010061068733763744, 0.01197014119071967, 0.012698376660930412, 0.011322714727807322, 0.01090459130971119, 0.009559547835004872, 0.010498329644215184, 0.013744539638296114, 0.01241662294231113, 0.013920642752687732, 0.01423709054803628, 0.012246450050510711, 0.011798522946392928, 0.012688307869426999, 0.011701411381744752, 0.010963335503348588, 0.011647664158072897, 0.010302975974959723, 0.011032126020690486, 0.011010171786400392, 0.012070857382656947, 0.012159491609797332, 0.012081814922729835, 0.01247264839908568, 0.011200668427178526, 0.013522937854091454, 0.013105319145488747, 0.011579740082727632, 0.010161886039664997, 0.010813245667125665, 0.010532288644280609, 0.011266459095452083, 0.012340346869037254, 0.010811434812069995, 0.012335213065123588, 0.00964337831495196, 0.012187168579263661, 0.013251917772334918, 0.011498371091243765, 0.010161965220729475, 0.01181737461439899, 0.009784767253164278, 0.009961151956167761, 0.011528718678809174, 0.010742033001898285, 0.010879391169694889, 0.010316128340493388, 0.010247573903152401, 0.011062856172903247, 0.011417489025016405, 0.01119193099686049, 0.01098276222674334, 0.014026370736078514, 0.01043225916128798, 0.010194459622871846, 0.010815687677297334, 0.011134576557489298, 0.01142961024001326, 0.010585867874587658, 0.015266592121768712, 0.011920413664227218, 0.009755732807066618, 0.012196390497070086, 0.009319626263621858, 0.011464227492400532, 0.01051867055357182, 0.012952488572362067, 0.010904039384865052, 0.009939432112089783, 0.011965844054709109, 0.010384358639967691, 0.010084010975718798, 0.013177624524319486, 0.014828138853365192, 0.011173145152672995, 0.011590247202021274, 0.010687161551733566, 0.011618327067404854, 0.010073717357638476, 0.01029100516208372, 0.010542481107447248, 0.011770135056296347, 0.011999996525704881, 0.011631265879266437, 0.012435408590750038, 0.012434157423771368, 0.010597499708746861, 0.010475004983816016, 0.013007106929861926]],
            'im_sum': [False, 4153.80461124],
            'regn_sum': [False, 22.0496815369],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.013743433433999516, 0.012734920080011661, 0.012977344646935137, 0.012849700732137936, 0.012869319430463244, 0.012786113424398437, 0.012763233601525502]],
            'profile': [False, 0.0247435603272],
            'im_fit': [False, [0.07905756034384644, 12.886359417267892, \
                       5.310724943044715]],
            'im_fit_loc': [True, [511, 261.889137107219]],
            'im_fit_pix': [False, [299.5756533902605, 255.6464223194929]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [209, 387, 0, 0]), \
                      (img+'.image', False, [209, 388, 0, 0]), \
                      (img+'.image', True, [20, 204, 0, 0]), \
                      (img+'.image', False, [19, 204, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 191116800],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'mask_pix': [True, 9091],
            'mask_regns': [True, 33],
            'npts_real': [True, 191116800]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/mosaic_cube_eph.pb.crtf')

        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 105033551.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [211, 203, 0, 0]],
            'min_val': [False, 0.200000017881],
            'min_val_pos': [True, [380, 278, 0, 194]],
            'im_rms': [False, 0.638066212618],
            'npts_0.2': [True, 105033551],
            'npts_0.5': [True, 60649291],
            'npts_real': [True, 191116800],
            'pb_fit': [False, [1.1107636648844048, 36.56192593982793, \
                       36.54997859671651]],
            'pb_fit_loc': [True, [511, 261.889137107219]],
            'pb_fit_pix': [False, [211.18281521819424, 203.29195188890523]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/mosaic_cube_eph.psf.crtf')

        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.045659199357],
            'min_val_pos': [True, [230, 216, 0, 13]],
            'im_rms': [False, 0.0122782213796],
            'im_sum': [False, 130.74164209],
            'npts_real': [True, 191116800],
            'psf_fit': [False, [0.8855960185478, 1.0832399168917628, \
                        0.85106729999774]],
            'psf_fit_loc': [True, [511, 261.889137107219]],
            'psf_fit_pix': [False, [239.96614077346965, 209.99207884772346]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/mosaic_cube_eph.residual.crtf')

        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 105033551.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.0698709115386],
            'max_val_pos': [True, [276, 250, 0, 321]],
            'min_val': [False, -0.0710606276989],
            'min_val_pos': [True, [192, 206, 0, 674]],
            'im_rms': [False, 0.012272523599],
            'im_sum': [False, 3748.45107633],
            'regn_sum': [False, 5.95724748177],
            'npts_real': [True, 191116800]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/mosaic_cube_eph.model.crtf')

        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.0327642522752],
            'max_val_pos': [True, [288, 214, 0, 511]],
            'min_val': [False, -0.00106700242031],
            'min_val_pos': [True, [295, 238, 0, 589]],
            'im_rms': [False, 2.51895960856e-05],
            'im_sum': [False, 10.4715985146],
            'regn_sum': [False, 0.46631132951],
            'mask_non0': [True, 1459],
            'npts_real': [True, 191116800]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, 
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 45509.8398438],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 45195.1835938],
            'min_val_pos': [True, [0, 0, 0, 593]],
            'im_rms': [False, 45284.18444383],
            'npts_real': [True, 948]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight')

        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [True, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.3178614676],
            'max_val_pos': [True, [211, 202, 0, 944]],
            'min_val': [False, 6.66096239001e-05],
            'min_val_pos': [True, [451, 65, 0, 947]],
            'im_rms': [False, 0.119822765111],
            'im_sum': [False, 13818688.5072],
            'npts_0.2': [True, 28725785],
            'npts_0.3': [True, 6228582],
            'npts_real': [True, 191116800]}

        report9 = self.stats_compare(exp_wt_stats, wt_stats_dict, '.weight')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9


        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[-0.01, 0.1])
            self.mom8_creator(img+'.residual', range_list=[-0.01, 0.1])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[testname]['images'].append(img+'.image.profile.png')

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=False)

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
            parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', region_file = \
            data_path+'region_files/mosaic_mfs_eph.image.crtf', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 0.914226949215],
            'com_bmin': [False, 0.708592534065],
            'com_pa': [False, -89.3612976074],
            'npts': [True, 201600],
            'npts_unmasked': [True, 113589.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.06048989296],
            'max_val_pos': [True, [291, 212, 0, 0]],
            'min_val': [False, -2.1858522892],
            'min_val_pos': [True, [290, 152, 0, 0]],
            'im_rms': [False, 0.676557465791],
            'im_sum': [False, 5498.32523989],
            'regn_sum': [False, 8725.50744961],
            'npts_real': [True, 201600],
            'rms_per_field': [False, [0.93211966613880237, 0.75587343134288443, 0.89153649532739099, 0.78815835277889079, 0.9296910564725368, 0.78753090736250964, 0.81230003610334034]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [211, 390, 0, 0]), \
                      (img+'.image', False, [211, 391, 0, 0]), \
                      (img+'.image', True, [18, 205, 0, 0]), \
                      (img+'.image', False, [17, 205, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 201600],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 201600]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', region_file = \
            data_path+'region_files/mosaic_mfs_eph.pb.crtf')

        exp_pb_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 113589.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [211, 203, 0, 0]],
            'min_val': [False, 0.200001135468],
            'min_val_pos': [True, [81, 343, 0, 0]],
            'im_rms': [False, 0.63128956829],
            'npts_0.2': [True, 113589],
            'npts_0.5': [True, 64549],
            'npts_real': [True, 201600],
            'pb_fit': [False, [1.1028920489851448, 36.87369130938714, \
                       36.84614469310595]],
            'pb_fit_loc': [True, [0, 253.5744226264628]],
            'pb_fit_pix': [False, [211.18238054179585, 203.38210680708025]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', region_file = \
            data_path+'region_files/mosaic_mfs_eph.psf.crtf')

        exp_psf_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 201600.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.0525086708367],
            'min_val_pos': [True, [230, 216, 0, 0]],
            'im_rms': [False, 0.0111846421981],
            'im_sum': [False, 0.100949261998],
            'npts_real': [True, 201600],
            'psf_fit': [False, [0.8980212570855989, 1.0458854777504984, \
                        0.8222593788495552]],
            'psf_fit_loc': [True, [0, 253.5744226264628]],
            'psf_fit_pix': [False, [239.96621779301014, 209.99390876796625]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', region_file = \
            data_path+'region_files/mosaic_mfs_eph.residual.crtf')

        exp_resid_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 113589.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.06048989296],
            'max_val_pos': [True, [291, 212, 0, 0]],
            'min_val': [False, -2.1858522892],
            'min_val_pos': [True, [290, 152, 0, 0]],
            'im_rms': [False, 0.676557465791],
            'im_sum': [False, 5498.32523989],
            'regn_sum': [False, 8725.50744961],
            'npts_real': [True, 201600]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual')

        # .model report
        model_stats_dict = self.image_stats(img+'.model', region_file = \
            data_path+'region_files/mosaic_mfs_eph.model.crtf')

        exp_model_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 201600.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.0],
            'im_sum': [False, 0.0],
            'regn_sum': [False, 0.0],
            'mask_non0': [True, 0],
            'npts_real': [True, 201600]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 30068706.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 30068706.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 30068705.591],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt')

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight')

        exp_wt_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 201600.0],
            'freq_bin': [True, 16762504556.453735],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.333539396524],
            'max_val_pos': [True, [211, 203, 0, 0]],
            'min_val': [False, 9.39047822612e-05],
            'min_val_pos': [True, [451, 64, 0, 0]],
            'im_rms': [False, 0.12506836881],
            'im_sum': [False, 15366.9703442],
            'npts_0.2': [True, 32025],
            'npts_0.3': [True, 9855],
            'npts_real': [True, 201600]}

        report9 = self.stats_compare(exp_wt_stats, wt_stats_dict, '.weight')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9


        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image', range_list=[-2.2, 2.1])
            self.mom8_creator(img+'.residual', range_list=[-2.2, 2.1])
            test_dict[testname]['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

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
        img = os.getcwd()+'/'+file_name+'1'
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
            fastnoise=False, savemodel='none', parallel=False)

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
            restart=True, calcres=False, calcpsf=False, \
            parallel=False)

        report0 = th.checkall(imgexist = self.image_list(img, 'mos_mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.image.tt0.crtf', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 0.914226949215],
            'com_bmin': [False, 0.708592534065],
            'com_pa': [False, -89.3612976074],
            'npts': [True, 201600],
            'npts_unmasked': [True, 113589.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.05563259125],
            'max_val_pos': [True, [291, 212, 0, 0]],
            'min_val': [False, -2.19232487679],
            'min_val_pos': [True, [290, 152, 0, 0]],
            'im_rms': [False, 0.672551536032],
            'im_sum': [False, 5410.97026619],
            'regn_sum': [False, 8700.32473019],
            'npts_real': [True, 201600],
            'rms_per_field': [False, [0.92450562082276955, 0.74815259370356468, 0.88345818716087821, 0.78563887415834988, 0.92110354253354942, 0.78525487979592601, 0.80567239544155445]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [211, 390, 0, 0]), \
                      (img+'.image.tt0', False, [211, 391, 0, 0]), \
                      (img+'.image.tt0', True, [18, 205, 0, 0]), \
                      (img+'.image.tt0', False, [17, 205, 0, 0])])

        report2 = self.stats_compare(exp_im_stats, im_stats_dict, '.image.tt0')

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 201600],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 201600]}

        report3 = self.stats_compare(exp_mask_stats, mask_stats_dict, '.mask')

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.pb.tt0.crtf')

        exp_pb_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 113589.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [211, 203, 0, 0]],
            'min_val': [False, 0.200001135468],
            'min_val_pos': [True, [81, 343, 0, 0]],
            'im_rms': [False, 0.63128956829],
            'npts_0.2': [True, 113589],
            'npts_0.5': [True, 64549],
            'npts_real': [True, 201600],
            'pb_fit': [False, [1.1028920489851448, 36.87369130938714, \
                       36.84614469310595]],
            'pb_fit_loc': [True, [0, 253.5744226264628]],
            'pb_fit_pix': [False, [211.18238054179585, 203.38210680708025]]}

        report4 = self.stats_compare(exp_pb_stats, pb_stats_dict, '.pb.tt0')

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.psf.tt0.crtf')

        exp_psf_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 201600.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.0525086708367],
            'min_val_pos': [True, [230, 216, 0, 0]],
            'im_rms': [False, 0.0111846421981],
            'im_sum': [False, 0.100949258623],
            'npts_real': [True, 201600],
            'psf_fit': [False, [0.8980212570855989, 1.0458854777504984, \
                        0.8222593788495552]],
            'psf_fit_loc': [True, [0, 253.5744226264628]],
            'psf_fit_pix': [False, [239.96621779301014, 209.99390876796625]]}

        report5 = self.stats_compare(exp_psf_stats, psf_stats_dict, '.psf.tt0')

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            region_file = data_path+
            'region_files/mosaic_mtmfs_eph.residual.tt0.crtf')

        exp_resid_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 113589.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.06048989296],
            'max_val_pos': [True, [291, 212, 0, 0]],
            'min_val': [False, -2.1858522892],
            'min_val_pos': [True, [290, 152, 0, 0]],
            'im_rms': [False, 0.676557465789],
            'im_sum': [False, 5498.32524009],
            'regn_sum': [False, 8725.50744971],
            'npts_real': [True, 201600]}

        report6 = self.stats_compare(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0')

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.model.tt0.crtf')

        exp_model_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 201600.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.0],
            'im_sum': [False, 0.0],
            'regn_sum': [False, 0.0],
            'mask_non0': [True, 0],
            'npts_real': [True, 201600]}

        report7 = self.stats_compare(exp_model_stats, model_stats_dict, \
            '.model.tt0')

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 30068706.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 30068706.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 30068705.591],
            'npts_real': [True, 1]}

        report8 = self.stats_compare(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0')

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight.tt0')

        exp_wt_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 201600.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.333539396524],
            'max_val_pos': [True, [211, 203, 0, 0]],
            'min_val': [False, 9.39047822612e-05],
            'min_val_pos': [True, [451, 64, 0, 0]],
            'im_rms': [False, 0.12506836881],
            'im_sum': [False, 15366.9703442],
            'npts_0.2': [True, 32025],
            'npts_0.3': [True, 9855],
            'npts_real': [True, 201600]}

        report9 = self.stats_compare(exp_wt_stats, wt_stats_dict, \
            '.weight.tt0')

        # .alpha report
        alpha_stats_dict = self.image_stats(img+'.alpha', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.alpha.crtf')

        exp_alpha_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 49203.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 62.1377220154],
            'max_val_pos': [True, [173, 232, 0, 0]],
            'min_val': [False, -58.0292701721],
            'min_val_pos': [True, [204, 163, 0, 0]],
            'im_rms': [False, 9.26911023158],
            'im_sum': [False, 40833.0406524],
            'regn_sum': [False, 29081.164461],
            'mask_non0': [True, 49203],
            'npts_real': [True, 201600]}

        report10 = self.stats_compare(exp_alpha_stats, alpha_stats_dict, \
            '.alpha')

        # .alpha.error report
        error_stats_dict = self.image_stats(img+'.alpha.error', region_file = \
            data_path+'region_files/mosaic_mtmfs_eph.alpha.error.crtf')

        exp_error_stats = {'npts': [True, 201600],
            'npts_unmasked': [True, 49203.0],
            'freq_bin': [True, 16762504556.453705],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 87.8760070801],
            'max_val_pos': [True, [173, 232, 0, 0]],
            'min_val': [False, 0.000175367502379],
            'min_val_pos': [True, [224, 36, 0, 0]],
            'im_rms': [False, 13.1085011733],
            'im_sum': [False, 463643.715791],
            'regn_sum': [False, 105848.469828],
            'npts_real': [True, 201600]}

        report11 = self.stats_compare(exp_error_stats, error_stats_dict, \
            '.alpha.error')

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11

        add_to_dict(self, output=test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['report'] = report
        test_dict[testname]['images'] = []

        if not CASA6:
            self.mom8_creator(img+'.image.tt0', range_list=[-2.2, 2.1])
            self.mom8_creator(img+'.residual.tt0', range_list=[-2.2, 2.1])
            test_dict[testname]['images'].extend( \
                (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = report)


def suite():
     return [Test_standard, Test_mosaic]

# Main #
if __name__ == '__main__':
    unittest.main()



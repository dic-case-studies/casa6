##########################################################################
##########################################################################
# test_stk_alma_pipeline_imaging.py
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
# [https://open-jira.nrao.edu/browse/CAS-12428]
#
#
##########################################################################

'''
Datasets (MOUS)
E2E6.1.00034.S (uid://A002/Xcff05c/X1ec)
2018.1.00879.S (uid://A001/X133d/X169f)
E2E6.1.00020.S (uid://A002/Xcff05c/Xe5)
2017.1.00750.T (uid://A001/X131b/X57)

Test list - 18 total
- Single field(SF) cube - E2E6.1.00034.S
- SF cube with perchanweightdensity(PCWD)+briggsbwtaper - E2E6.1.00034.S
- SF MFS - E2E6.1.00020.S
- SF mtmfs - E2E6.1.00020.S
- SF ephemeris cube (multi-EB) - 2017.1.00750.T
- SF ephemeris cube (multi-EB) with PCWD+briggsbwtaper - 2017.1.00750.T
- SF ephemeris MFS - 2018.1.00879.S
- SF ephemeris mtmfs - 2018.1.00879.S
- SF Calibrator - E2E6.1.00034.S
- SF ephemeris Calibrator - 2018.1.00879.S  
- Mosaic cube - E2E6.1.00034.S
- Mosaic cube with PCWD+briggsbwtaper- E2E6.1.00034.S
- Mosaic MFS - E2E6.1.00020.S
- Mosaic mtmfs - E2E6.1.00020.S
- Mosaic ephemeris cube - 2018.1.00879.S
- Mosaic ephemeris cube with PCWD+briggsbwtaper - 2018.1.00879.S
- Mosaic ephemeris MFS - 2018.1.00879.S
- Mosaic ephemeris mtmfs - 2018.1.00879.S

Each test stores reference values in dictionaries for the metrics 
to be tested:
The following variable names are used:
the standard default sets are,
exp_im_stats, exp_mask_stats, exp_pb_stats, exp_psf_stats, 
exp_model_stats, exp_resid_stats, exp_sumwt_stats
exp_wt_stats (for mosaic)
Addtionally, for cube imaging (if self.parallel=True), 
exp_bmin_dict, exp_bmaj_dict, exp_pa_dict
And for mtmfs
exp_im1_stats, exp_model1_stats, exp_resid1_stats, exp_sumwt1_stats

'''

##########################################################################
##########################################################################

# Imports #
import os
import glob
import sys
import subprocess
import unittest
import numpy
import shutil
import inspect
import scipy
import matplotlib.pyplot as pyplot
import json
import pickle
import copy

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
    ctsys_resolve = ctsys.resolve

except ImportError:
    from __main__ import default  # reset given task to its default values
    from tasks import *  # Imports all casa tasks
    from taskinit import *  # Imports all casa tools
    from parallel.parallel_task_helper import ParallelTaskHelper

    _ia = iatool()
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0], 'casatestdata/')
        return os.path.join(dataPath,apath)

# location of data
data_path = ctsys_resolve('stakeholder/alma/')

# save the dictionaries of the metrics to files (per tests) 
# mostly useful for the maintenace (updating the expected metric i
# parameters based
# on the current metrics)
savemetricdict=True

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
        #print("TEST_DICT=",test_dict)
        generate_weblog("tclean_ALMA_pipeline",test_dict)
        print("Closing ia tool")
        self._myia.done()
        """ don't delete it all """
#        self.delData()

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self, msname=[""]):
        if msname != [""]:
             self.msfile=msname

    def delData(self, msname=[""]):
        del_files = [self.img_subdir]
        if msname != [""]:
             self.msfile=msname
        if (os.path.exists(self.msfile)):
            del_files.append(self.msfile)
        img_files = glob.glob(self.img+'*')
        del_files += img_files
        for f in del_files:
            shutil.rmtree(f)

    def prepInputmask(self, maskname=""):
        if maskname!="":
            self.maskname=maskname
        if (os.path.exists(self.maskname)):
            shutil.rmtree(self.maskname)
        shutil.copytree(refdatapath+self.maskname, self.maskname, symlinks=True)

    def checkfinal(self, pstr=""):
        pstr += "["+inspect.stack()[1][3]+"] : To re-run this test : " \
                "runUnitTest.main(['test_tclean["+ inspect.stack()[1][3] \
                +"]'])"
        casalog.post(pstr,'INFO')
        if(pstr.count("(Fail") > 0 ):
             self.fail("\n"+pstr)

    def check_dict_vals_beam(self, exp_dict, act_dict, suffix, epsilon=0.01):
        """ Compares expected dictionary with actual dictionary. Useful for comparing the restoring beam.

            Parameters
            ----------
            exp_dict: dictionary
                Expected values, as key:value pairs.
                Keys must match between exp_dict and act_dict.
                Values are compared between exp_dict and act_dict. A summary
                line is returned with the first mismatched value, or the last
                successfully matched value.
            act_dict: dictionary
                Actual values to compare to exp_dict (and just the values).
            suffix: string
                For use with summary print statements.
        """
        report = ''
        eps = epsilon
        passed = True
        chans = 0
        for key in exp_dict:
            result = th.check_val(act_dict[key], exp_dict[key],
                valname=suffix+' chan'+str(chans), epsilon=eps)[1]
            chans += 1
            if 'Fail' in result:
                passed = False
                break
        report += th.check_val(passed, True, valname=suffix+' chan'+str(chans), exact=True)[1]

        return report

    def copy_products(self, old_pname, new_pname, ignore=None):
        """ function to copy iter0 images to iter1 images
            (taken from pipeline)
        """
        imlist = glob.glob('%s.*' % old_pname)
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

    def cube_beam_stats(self, image):
        """ function to return per-channel beam statistics
            will be deprecated and combined into image_stats 
            once CASA beam issue is fixed
        """
        self._myia.open(image)

        bmin_dict = {}; bmaj_dict = {}; pa_dict = {}
        beam_dict = self._myia.restoringbeam()['beams']
        for item in beam_dict.keys():
            bmin_dict[item] = beam_dict[item]['*0']['minor']['value']
            bmaj_dict[item] = beam_dict[item]['*0']['major']['value']
            pa_dict[item] = beam_dict[item]['*0']['positionangle']['value']

        self._myia.close()

        return bmin_dict, bmaj_dict, pa_dict

    def image_stats(self, image, fit_region=None, field_regions=None, masks=None):
        """ function that takes an image file and returns a statistics
            dictionary
        """
        self._myia.open(image)
        imagename=os.path.basename(image)
        stats_dict = {}

        statistics = self._myia.statistics()
        # Return data chunk; transpose to make channel selection easier
        chunk = numpy.transpose(self._myia.getchunk(dropdeg=True))

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
        if not image.endswith('.mask'):
            stats_dict['max_val'] = statistics['max'][0]
            stats_dict['max_val_pos'] = statistics['maxpos'].tolist()
            max_loc = [stats_dict['max_val_pos'][0], \
                stats_dict['max_val_pos'][1]]
            stats_dict['min_val'] = statistics['min'][0]
            stats_dict['min_val_pos'] = statistics['minpos'].tolist()
            stats_dict['im_rms'] = statistics['rms'][0]

        # stats returned if a region file is given
        if fit_region != None:
            if '_cube' in imagename:
                if '.pb' in imagename:
                    fit_region = fit_region + ', range=[%schan,%schan]'\
                        % (int(im_size[3]/2), int(im_size[3]/2))
                if '.psf' in imagename:
                    # using chan 1 as first because ia.fitcomponents fits
                    # every channel if chan=0
                    fit_regions = [(fit_region + ', range=[%schan,%schan]' \
                                   % (1, 1)), \
                                  (fit_region + ', range=[%schan,%schan]' \
                                   % (int(im_size[3]/2), int(im_size[3]/2))), \
                                  (fit_region + ', range=[%schan,%schan]' \
                                   % ((im_size[3]-1), (im_size[3]-1)))]
                    i = 0
                    for region in fit_regions:
                        try:
                            fit_dict = self._myia.fitcomponents( \
                                region=region)['results']['component0']
                            stats_dict['fit_'+str(i)] = [ \
                                fit_dict['peak']['value'], \
                                fit_dict['shape']['majoraxis']['value'], \
                                fit_dict['shape']['minoraxis']['value']]
                            stats_dict['fit_loc_chan_'+str(i)] = fit_dict['spectrum']['channel']
                            stats_dict['fit_loc_freq_'+str(i)] = \
                                fit_dict['spectrum']['frequency']['m0']['value']
                            stats_dict['fit_pix_'+str(i)] = \
                                fit_dict['pixelcoords'].tolist()
                        except KeyError:
                            stats_dict['fit_'+str(i)] = [1.0, 1.0, 1.0]
                            stats_dict['fit_loc_chan_'+str(i)] = 1.0
                            stats_dict['fit_loc_freq_'+str(i)] = 1.0
                            stats_dict['fit_pix_'+str(i)] = [1.0, 1.0]
                        i += 1
                if '.model' in imagename:
                    fit_region = fit_region
                if '.model' not in imagename and '.pb' not in imagename and '.psf' not in imagename:
                    # WARN: If max value channel is 0, tool fits all channels
                    fit_region = fit_region + ', range=[%schan,%schan]' \
                        % (stats_dict['max_val_pos'][3], \
                        stats_dict['max_val_pos'][3])
            if '.psf' in imagename and '_cube' in imagename:
                stats_dict['regn_sum'] = self._myia.statistics( \
                    region=fit_regions[1])['sum'][0]
            else:
                stats_dict['regn_sum'] = self._myia.statistics( \
                    region=fit_region)['sum'][0]
            if ('image' in imagename and 'mosaic_cube_eph' not in imagename) or 'pb' in imagename or ('psf' in imagename and 'cube' not in imagename):
                try:
                    fit_dict = self._myia.fitcomponents( \
                        region=fit_region)['results']['component0']
                    stats_dict['fit'] = [fit_dict['peak']['value'], \
                        fit_dict['shape']['majoraxis']['value'], \
                        fit_dict['shape']['minoraxis']['value']]
                    stats_dict['fit_loc_chan'] = fit_dict['spectrum']['channel']
                    stats_dict['fit_loc_freq'] = fit_dict['spectrum']['frequency']['m0']['value']
                    stats_dict['fit_pix'] = fit_dict['pixelcoords'].tolist()
                except KeyError:
                    stats_dict['fit'] = [1.0, 1.0, 1.0]
                    stats_dict['fit_loc_chan'] = 1.0
                    stats_dict['fit_loc_freq'] = 1.0
                    stats_dict['fit_pix'] = [1.0, 1.0]

        # stats returned for .image(.tt0)
        if 'image' in imagename:
            commonbeam = self._myia.commonbeam()
            stats_dict['com_bmin'] = commonbeam['minor']['value']
            stats_dict['com_bmaj'] = commonbeam['major']['value']
            stats_dict['com_pa'] = commonbeam['pa']['value']
            if 'cube' in imagename:
                stats_dict['rms_per_chan'] = \
                    self._myia.statistics(axes=[0,1])['rms'].tolist()
                stats_dict['profile'] = self.cube_profile_fit( \
                    image, max_loc, stats_dict['nchan'])
            if 'mosaic' in imagename:
                stats_dict['rms_per_field'] = []
                for region in field_regions:
                    stats_dict['rms_per_field'].append( \
                        self._myia.statistics(region=region)['rms'][0])

        # stats returned if not .pb(.tt0), .sumwt(.tt0), or .mask
        # if 'pb' not in image and 'sumwt' not in image and not image.endswith('.mask'):
        stats_dict['im_sum'] = statistics['sum'][0]

        if image.endswith('.mask'):
            stats_dict['mask_pix'] = numpy.count_nonzero(chunk)
            stats_dict['mask_regns'] = scipy.ndimage.label(chunk)[1]
            stats_dict['mask'] = ~numpy.array(chunk, dtype=bool)

        if 'pb' in imagename:
            pb_mask_02 = chunk>0.2
            pb_mask_05 = chunk>0.5
            if 'cube' in image:
                pb_02_list = []
                pb_05_list = []
                i = 0
                for chan in chunk:
                    pb_02_list.append(numpy.count_nonzero(chan*pb_mask_02[i]))
                    pb_05_list.append(numpy.count_nonzero(chan*pb_mask_05[i]))
                    i += 1
                stats_dict['npts_0.2'] = pb_02_list
                stats_dict['npts_0.5'] = pb_05_list
            else:
                stats_dict['npts_0.2'] = numpy.count_nonzero(pb_mask_02)
                stats_dict['npts_0.5'] = numpy.count_nonzero(pb_mask_05)
            if 'mosaic' in imagename:
                stats_dict['pb_mask_0.2'] = pb_mask_02
                stats_dict['pb_mask_0.5'] = pb_mask_05

        if 'model' in imagename or image.endswith('.alpha'):
            stats_dict['mask_non0'] = numpy.count_nonzero(chunk*masks)

        if 'weight' in imagename:
            if 'cube' in imagename:
                wt_02_list = []
                wt_05_list = []
                i = 0
                for chan in chunk:
                    wt_02_list.append(numpy.count_nonzero(chan*masks[0][i]))
                    wt_05_list.append(numpy.count_nonzero(chan*masks[1][i]))
                    i += 1
                stats_dict['npts_0.2'] = wt_02_list
                stats_dict['npts_0.5'] = wt_05_list
            else:
                stats_dict['npts_0.2'] = numpy.count_nonzero(chunk*masks[0])
                stats_dict['npts_0.5'] = numpy.count_nonzero(chunk*masks[1])

        self._myia.close()

        return stats_dict

    def image_list(self, image, mode):
        """ function used to return expected imaging output files """
        standard = [image+'.psf', image+'.residual', image+'.image', \
            image+'.image.pbcor', image+'.mask', image+'.pb', image+'.model', \
            image+'.sumwt']
        mosaic = [image+'.weight']
        mtmfs = [image+'.alpha', image+'.alpha.error', image+'.alpha.pbcor', \
           image+'.psf.tt0', image+'.psf.tt1', image+'.psf.tt2', \
           image+'.residual.tt0', image+'.residual.tt1', image+'.image.tt0',\
           image+'.image.tt1', image+'.image.tt0.pbcor', image+'.image.tt1.pbcor', \
           image+'.mask', image+'.pb.tt0', image+'.model.tt0', image+'.model.tt1', \
           image+'.sumwt.tt0', image+'.sumwt.tt1', image+'.sumwt.tt2']
        mos_mtmfs = [image+'.weight.tt0', image+'.weight.tt1', image+'.weight.tt2']

        if mode == 'standard':
            img_list = standard
        if mode == 'mosaic':
            img_list = standard+mosaic
        if mode == 'mtmfs':
            img_list = mtmfs
        if mode == 'mos_mtmfs':
            img_list = mtmfs+mos_mtmfs

        return img_list

    def mom8_creator(self, image, range_list):
        """ function that takes and image and turns it into a .png for
            weblog
        """
        immoments(imagename = image, moments = 8, outfile = image+'.moment8')
        imview(raster={'file': image+'.moment8', 'range': range_list}, \
            out = {'file': image+'.moment8.png'})
        subprocess.call('mogrify -trim '+image+'.moment8.png', shell=True)

    def cube_profile_fit(self, image, max_loc, nchan):
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
        pyplot.savefig(image+'.profile.png')
        pyplot.clf()

        return profile

    def filter_report(self, report, showonlyfail=True):
        """ function to filter the test report, the input report is expected to be a string with the newline code """
        ret = ''
        if showonlyfail:
            filter='Fail'
        else:
            filter='Pass' 
        
        if report!='':
            testItems = report.split('\n')
            retitems=[]
            for testitem in testItems:
                if '[ check_ims ]' in testitem or '[ check_pixmask ]' in testitem or '[ check_val ]' in testitem:
                    if '( '+filter in testitem:
                        retitems.append(testitem)
            nfail = len(retitems)
            msg = str(nfail)+' individual test failure(s) '
            ret = '\n' + '\n'.join(retitems)
            ret += '\n' + msg
        return ret

     
    def save_dict_to_file(self, topkey, indict, outfilename, outformat='JSON'):
        """ function that will save input Python dictionaries to a JSON file (default)
            or pickle file. 
            Create a separate file with outfilename_topkey for now
        """
        if type(indict)!=dict:
            print("indict is not a dict. Saved file may not be in correct format") 
        nestedDict={}
        nestedDict[topkey]=indict 
        print("Saving %s dictionaries", len(indict))
        if outformat == 'pickle':
            # writing to pickle: note if writing this way (without protocol=2) 
            # in casa6 and read in casa5 it will fail 
            with open(outfilename+'_'+topkey+'.pickle', 'wb') as outf:
                pickle.dump(nestedDict, outf)
        elif outformat== 'JSON':
            with open(outfilename+'.json', 'w') as outf:
                json.dump(nestedDict, outf)
        else:
            print("no saving with format:", outformat)

    def modify_dict(self, output=None, testname=None, parallel=None):
        ''' Modified test_dict costructed by casatestutils add_to_dict to include only 
            the task commands executed and also add self.parallel value to the dictionary.
            The cube imaging cases usually have if-else conditional based on parallel mode is on or not
            to trigger different set of tclean commands.
            Assumption: self.parallel is used to trigger different tclean commands at iter1 step.
            For self.parallel=True, iter1 has two tclean commands (2nd and 3rd tclean commands within
            each relevante test(cube) and so in test_dict['taskcall'], 1st(iter0) and 2nd and 3rd commands
            are the ones acutually executed and should remove 4th (self.parallel=False) case.
        '''
        if testname in output:
            if 'taskcall' in output[testname] and len(output[testname]['taskcall'])==4: 
                if parallel:
                    # 0,1,2th in the list are used pop last one
                    output[testname]['taskcall'].pop()
                else:
                    output[testname]['taskcall'].pop(1)
                    output[testname]['taskcall'].pop(1)
            output[testname]['self.parallel']=parallel



##############################################
# TESTS
##############################################
test_dict = {}
class Test_standard(test_tclean_base):


    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cube(self):
        ''' Standard (single field) cube imaging - central field of SMIDGE_NWCloud (field 3), spw 22 '''

        file_name = 'standard_cube.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, imagename=file_name+'0', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS'
            ' 00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526743594GHz', width='0.2441741MHz', \
            outframe='LSRK', pblimit=0.2, perchanweightdensity=False, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            deconvolver='hogbom', usepointing=False, restoration=False, \
            pbcor=False, weighting='briggs', restoringbeam='common', \
            robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
                nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
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
                parallel=True, verbose=True)

            # retrieve per-channel beam statistics (only in parallel)
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(image=img+'.image')

            tclean(vis=self.msfile, imagename=file_name+'1', spw=['0'], \
                field='1', imsize=[80, 80], cell=['1.1arcsec'], \
                phasecenter='ICRS 00:45:54.3836 -073.15.29.413', stokes='I', \
                antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                specmode='cube', nchan=508, start='220.2526743594GHz', \
                width='0.2441741MHz', outframe='LSRK', \
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
                parallel=False, verbose=True)

        else:
            tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
                spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
                scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
                datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
                '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
                nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
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
                restoringbeam='common', parallel=False, verbose=True)


        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report(test_standard_cube)
        im_stats_dict = self.image_stats(image=img+'.image', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 8.509892605313942],
            'com_bmin': [False, 5.950050676606115],
            'com_pa': [False, 72.54607919421503],
            'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 66.089371562004089],
            'npts_real': [True, 3251200],
            'profile': [False, 0.90354546054915941],
            'fit': [False, [0.9222914423989481, 14.289411729097703, \
                    6.881365508161659]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [38.263402177385942, 37.306443753086633]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 70, 0, 0]), \
                (img+'.image', False, [40, 71, 0, 0]), \
                (img+'.image', True, [10, 40, 0, 0]), \
                (img+'.image', False, [9, 40, 0, 0])])

        # .image report
        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(image=img+'.mask')

        # test_standard_cube.exp_mask_stats
        exp_mask_stats = {'npts': [True, 3251200],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'mask_pix': [False, 437],
            'mask_regns': [True, 1],
            'npts_real': [True, 3251200]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube.exp_mask_stats
        exp_pb_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'im_rms': [False, 0.578238326026],
            'npts_0.2': [False, [2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997]],
            'npts_0.5': [False, [1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449]],
            'npts_real': [True, 3251200],
            'fit': [False, [1.0308127949041446, 46.61751391582679, \
                       46.61253844001269]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [40.00032808200995, 40.00099739969875]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(image=img+'.psf', fit_region = \
            'ellipse[[11.47648725deg, -73.25812003deg], [8.0291arcsec, 6.8080arcsec], 90.00000000deg]')

        # test_standard_cube.exp_psf_stats
        exp_psf_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.1],
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
            'fit_0': [False, [1.0959520385253885, 7.675969776744627, \
                      5.143545685744538]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [39.998799449215596, 39.99523672953224]],
            'fit_1': [False, [1.0959863390592945, 7.672871552789668, \
                      5.141790170376213]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [39.99880225653659, 39.99524870969922]],
            'fit_2': [False, [1.0960422882714267, 7.669928861314425, \
                      5.140004109591353]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [39.9988051116427, 39.995258207738075]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube.exp_resid_stats
        exp_resid_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 19.0307790947],
            'npts_real': [True, 3251200]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube.exp_model_stats
        exp_model_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.1],
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
            'mask_non0': [True, 0],
            'npts_real': [True, 3251200]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        # test_standard_cube.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.1],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube.exp_bmin_dict
            exp_bmin_dict = {'*0': 5.944103717803955, '*1': 5.94410514831543, '*10': 5.944015026092529, '*100': 5.94343900680542, '*101': 5.943426609039307, '*102': 5.943421840667725, '*103': 5.943414211273193, '*104': 5.943387985229492, '*105': 5.943387985229492, '*106': 5.943348407745361, '*107': 5.943348407745361, '*108': 5.94331693649292, '*109': 5.94331693649292, '*11': 5.943998336791992, '*110': 5.94331693649292, '*111': 5.94331693649292, '*112': 5.94331693649292, '*113': 5.94331693649292, '*114': 5.94331693649292, '*115': 5.94331693649292, '*116': 5.94331693649292, '*117': 5.943300247192383, '*118': 5.9432783126831055, '*119': 5.9432783126831055, '*12': 5.94395637512207, '*120': 5.943252086639404, '*121': 5.943252086639404, '*122': 5.943222999572754, '*123': 5.943180561065674, '*124': 5.943145275115967, '*125': 5.943161487579346, '*126': 5.943154811859131, '*127': 5.943154811859131, '*128': 5.943154811859131, '*129': 5.9431471824646, '*13': 5.94395637512207, '*130': 5.9431471824646, '*131': 5.9431471824646, '*132': 5.9431471824646, '*133': 5.943125247955322, '*134': 5.943125247955322, '*135': 5.94309663772583, '*136': 5.94309663772583, '*137': 5.943094253540039, '*138': 5.94307279586792, '*139': 5.94307279586792, '*14': 5.94395637512207, '*140': 5.943065166473389, '*141': 5.943065166473389, '*142': 5.943065166473389, '*143': 5.943065166473389, '*144': 5.943065166473389, '*145': 5.943065166473389, '*146': 5.943028926849365, '*147': 5.943028926849365, '*148': 5.943004131317139, '*149': 5.942957878112793, '*15': 5.94395637512207, '*150': 5.942957878112793, '*151': 5.942957878112793, '*152': 5.94295072555542, '*153': 5.94295072555542, '*154': 5.942929267883301, '*155': 5.942929267883301, '*156': 5.94288444519043, '*157': 5.94288444519043, '*158': 5.94288444519043, '*159': 5.94288444519043, '*16': 5.943938255310059, '*160': 5.94288444519043, '*161': 5.942834854125977, '*162': 5.942834854125977, '*163': 5.942787170410156, '*164': 5.942765712738037, '*165': 5.942748546600342, '*166': 5.942732810974121, '*167': 5.942732810974121, '*168': 5.94273042678833, '*169': 5.94273042678833, '*17': 5.94391393661499, '*170': 5.94273042678833, '*171': 5.94273042678833, '*172': 5.942715644836426, '*173': 5.942707538604736, '*174': 5.942705154418945, '*175': 5.942681789398193, '*176': 5.942651271820068, '*177': 5.9426374435424805, '*178': 5.942581653594971, '*179': 5.942581653594971, '*18': 5.943907260894775, '*180': 5.942537307739258, '*181': 5.942537307739258, '*182': 5.942537307739258, '*183': 5.942534923553467, '*184': 5.942537784576416, '*185': 5.942537784576416, '*186': 5.942508697509766, '*187': 5.942508697509766, '*188': 5.942508697509766, '*189': 5.942508697509766, '*19': 5.943899154663086, '*190': 5.942508697509766, '*191': 5.942508697509766, '*192': 5.942503929138184, '*193': 5.942503929138184, '*194': 5.942503929138184, '*195': 5.942503929138184, '*196': 5.942503929138184, '*197': 5.942503929138184, '*198': 5.942503452301025, '*199': 5.9424729347229, '*2': 5.9440813064575195, '*20': 5.943899154663086, '*200': 5.9424729347229, '*201': 5.9424729347229, '*202': 5.94244384765625, '*203': 5.942445278167725, '*204': 5.942437171936035, '*205': 5.942437171936035, '*206': 5.942437171936035, '*207': 5.942437171936035, '*208': 5.942437171936035, '*209': 5.942437171936035, '*21': 5.943899154663086, '*210': 5.9424333572387695, '*211': 5.9424333572387695, '*212': 5.9424333572387695, '*213': 5.9424333572387695, '*214': 5.9424333572387695, '*215': 5.9424333572387695, '*216': 5.942417144775391, '*217': 5.942417144775391, '*218': 5.942412376403809, '*219': 5.942412376403809, '*22': 5.943899154663086, '*220': 5.942407131195068, '*221': 5.942407131195068, '*222': 5.942407131195068, '*223': 5.942407131195068, '*224': 5.942375659942627, '*225': 5.94236946105957, '*226': 5.94236946105957, '*227': 5.942365646362305, '*228': 5.942312717437744, '*229': 5.942312717437744, '*23': 5.943899154663086, '*230': 5.942312717437744, '*231': 5.942311763763428, '*232': 5.942311763763428, '*233': 5.942311763763428, '*234': 5.942311763763428, '*235': 5.942282676696777, '*236': 5.942282676696777, '*237': 5.942237854003906, '*238': 5.942237854003906, '*239': 5.942210674285889, '*24': 5.9438910484313965, '*240': 5.942195415496826, '*241': 5.942195415496826, '*242': 5.942195415496826, '*243': 5.942178249359131, '*244': 5.942124366760254, '*245': 5.942124366760254, '*246': 5.942124366760254, '*247': 5.942124366760254, '*248': 5.942124366760254, '*249': 5.942124366760254, '*25': 5.943882465362549, '*250': 5.942102909088135, '*251': 5.942102909088135, '*252': 5.942080497741699, '*253': 5.942091464996338, '*254': 5.942091464996338, '*255': 5.94207763671875, '*256': 5.942057132720947, '*257': 5.942053318023682, '*258': 5.942053318023682, '*259': 5.942053318023682, '*26': 5.943882465362549, '*260': 5.942053318023682, '*261': 5.94204568862915, '*262': 5.94204568862915, '*263': 5.94204568862915, '*264': 5.94204568862915, '*265': 5.94204568862915, '*266': 5.941995143890381, '*267': 5.941966533660889, '*268': 5.941969394683838, '*269': 5.941969394683838, '*27': 5.943882465362549, '*270': 5.94195032119751, '*271': 5.941928386688232, '*272': 5.941928386688232, '*273': 5.941908359527588, '*274': 5.941908359527588, '*275': 5.941908359527588, '*276': 5.941908359527588, '*277': 5.941891670227051, '*278': 5.941891670227051, '*279': 5.941891193389893, '*28': 5.943875789642334, '*280': 5.941891193389893, '*281': 5.941891193389893, '*282': 5.941891193389893, '*283': 5.941833972930908, '*284': 5.941833972930908, '*285': 5.941833972930908, '*286': 5.941833972930908, '*287': 5.941833972930908, '*288': 5.941832542419434, '*289': 5.941824913024902, '*29': 5.943875789642334, '*290': 5.941824913024902, '*291': 5.941812038421631, '*292': 5.941808223724365, '*293': 5.941808223724365, '*294': 5.941808223724365, '*295': 5.941808223724365, '*296': 5.941801071166992, '*297': 5.941802501678467, '*298': 5.941776752471924, '*299': 5.941746711730957, '*3': 5.9440813064575195, '*30': 5.943875789642334, '*300': 5.9417266845703125, '*301': 5.9417266845703125, '*302': 5.941720008850098, '*303': 5.941673278808594, '*304': 5.941673278808594, '*305': 5.9416584968566895, '*306': 5.9416584968566895, '*307': 5.9416584968566895, '*308': 5.9416584968566895, '*309': 5.9416584968566895, '*31': 5.943875789642334, '*310': 5.941618919372559, '*311': 5.941612243652344, '*312': 5.941611289978027, '*313': 5.941611289978027, '*314': 5.941610336303711, '*315': 5.941610336303711, '*316': 5.941610336303711, '*317': 5.941610336303711, '*318': 5.941595077514648, '*319': 5.941595077514648, '*32': 5.943851947784424, '*320': 5.941573143005371, '*321': 5.941573143005371, '*322': 5.941573143005371, '*323': 5.94158935546875, '*324': 5.94158935546875, '*325': 5.941571235656738, '*326': 5.941502571105957, '*327': 5.941502571105957, '*328': 5.941502571105957, '*329': 5.941502571105957, '*33': 5.943851947784424, '*330': 5.941502571105957, '*331': 5.941502571105957, '*332': 5.941502571105957, '*333': 5.9414777755737305, '*334': 5.9414777755737305, '*335': 5.941437721252441, '*336': 5.941437721252441, '*337': 5.941437721252441, '*338': 5.94136381149292, '*339': 5.94136381149292, '*34': 5.943851947784424, '*340': 5.94134521484375, '*341': 5.941339015960693, '*342': 5.941339015960693, '*343': 5.941314220428467, '*344': 5.941311836242676, '*345': 5.941303253173828, '*346': 5.941303253173828, '*347': 5.941303253173828, '*348': 5.941303253173828, '*349': 5.941303253173828, '*35': 5.943843364715576, '*350': 5.941303253173828, '*351': 5.941265106201172, '*352': 5.941249847412109, '*353': 5.941249847412109, '*354': 5.941246509552002, '*355': 5.941246509552002, '*356': 5.941233158111572, '*357': 5.941233158111572, '*358': 5.941233158111572, '*359': 5.941233158111572, '*36': 5.943843364715576, '*360': 5.941233158111572, '*361': 5.9412522315979, '*362': 5.9412522315979, '*363': 5.9412522315979, '*364': 5.941198825836182, '*365': 5.941179275512695, '*366': 5.941154479980469, '*367': 5.941154479980469, '*368': 5.941154479980469, '*369': 5.9411540031433105, '*37': 5.943843364715576, '*370': 5.941145896911621, '*371': 5.941145896911621, '*372': 5.941145896911621, '*373': 5.941145896911621, '*374': 5.941145896911621, '*375': 5.941142559051514, '*376': 5.941142559051514, '*377': 5.941142559051514, '*378': 5.941142559051514, '*379': 5.941142559051514, '*38': 5.943840980529785, '*380': 5.941142559051514, '*381': 5.941142559051514, '*382': 5.941142559051514, '*383': 5.941102504730225, '*384': 5.941102504730225, '*385': 5.941101551055908, '*386': 5.941097259521484, '*387': 5.94110107421875, '*388': 5.94110107421875, '*389': 5.941096305847168, '*39': 5.943840980529785, '*390': 5.941096305847168, '*391': 5.941074848175049, '*392': 5.941074848175049, '*393': 5.941074848175049, '*394': 5.94106912612915, '*395': 5.941043376922607, '*396': 5.941043376922607, '*397': 5.941043376922607, '*398': 5.941043376922607, '*399': 5.941043376922607, '*4': 5.9440813064575195, '*40': 5.943840980529785, '*400': 5.941017150878906, '*401': 5.940994739532471, '*402': 5.940994739532471, '*403': 5.940974712371826, '*404': 5.940969467163086, '*405': 5.940962791442871, '*406': 5.940962791442871, '*407': 5.940912246704102, '*408': 5.940912246704102, '*409': 5.940912246704102, '*41': 5.943840980529785, '*410': 5.940889835357666, '*411': 5.940889835357666, '*412': 5.940889835357666, '*413': 5.940889835357666, '*414': 5.940862655639648, '*415': 5.940829277038574, '*416': 5.940829277038574, '*417': 5.940829277038574, '*418': 5.940817356109619, '*419': 5.9407854080200195, '*42': 5.943844795227051, '*420': 5.940733909606934, '*421': 5.940733909606934, '*422': 5.940733909606934, '*423': 5.940733909606934, '*424': 5.940709114074707, '*425': 5.940709114074707, '*426': 5.940701961517334, '*427': 5.940701961517334, '*428': 5.940701961517334, '*429': 5.940701961517334, '*43': 5.943844795227051, '*430': 5.940701961517334, '*431': 5.940697193145752, '*432': 5.940697193145752, '*433': 5.940697193145752, '*434': 5.940697193145752, '*435': 5.940667629241943, '*436': 5.940667629241943, '*437': 5.940667629241943, '*438': 5.940667629241943, '*439': 5.940667629241943, '*44': 5.943831920623779, '*440': 5.940667629241943, '*441': 5.9406633377075195, '*442': 5.9406633377075195, '*443': 5.9406633377075195, '*444': 5.9406633377075195, '*445': 5.9406633377075195, '*446': 5.9406633377075195, '*447': 5.9406633377075195, '*448': 5.9406633377075195, '*449': 5.940667629241943, '*45': 5.94381856918335, '*450': 5.94066858291626, '*451': 5.94066858291626, '*452': 5.940656661987305, '*453': 5.940609931945801, '*454': 5.940609931945801, '*455': 5.940609931945801, '*456': 5.940608024597168, '*457': 5.940608024597168, '*458': 5.940586566925049, '*459': 5.940587520599365, '*46': 5.943813800811768, '*460': 5.940587520599365, '*461': 5.940557479858398, '*462': 5.940558910369873, '*463': 5.940558910369873, '*464': 5.940558910369873, '*465': 5.940556049346924, '*466': 5.940556049346924, '*467': 5.940556049346924, '*468': 5.940556049346924, '*469': 5.940556049346924, '*47': 5.943813800811768, '*470': 5.940556049346924, '*471': 5.940539836883545, '*472': 5.940539836883545, '*473': 5.940534591674805, '*474': 5.940478324890137, '*475': 5.940452575683594, '*476': 5.940426826477051, '*477': 5.940426826477051, '*478': 5.9403886795043945, '*479': 5.940366744995117, '*48': 5.943813800811768, '*480': 5.940366744995117, '*481': 5.940366744995117, '*482': 5.940362453460693, '*483': 5.940362453460693, '*484': 5.940362453460693, '*485': 5.940317630767822, '*486': 5.940317630767822, '*487': 5.940317630767822, '*488': 5.940317630767822, '*489': 5.940317630767822, '*49': 5.943813800811768, '*490': 5.940311908721924, '*491': 5.940311908721924, '*492': 5.940311908721924, '*493': 5.940311908721924, '*494': 5.940311908721924, '*495': 5.940311908721924, '*496': 5.940311908721924, '*497': 5.940313816070557, '*498': 5.940313816070557, '*499': 5.940313816070557, '*5': 5.944080829620361, '*50': 5.943783283233643, '*500': 5.940237045288086, '*501': 5.940237045288086, '*502': 5.940237045288086, '*503': 5.9401655197143555, '*504': 5.9401655197143555, '*505': 5.940146446228027, '*506': 5.940146446228027, '*507': 5.940146446228027, '*51': 5.943783283233643, '*52': 5.943783283233643, '*53': 5.943783283233643, '*54': 5.943783283233643, '*55': 5.943783283233643, '*56': 5.943750381469727, '*57': 5.943746089935303, '*58': 5.943742275238037, '*59': 5.943742275238037, '*6': 5.944068908691406, '*60': 5.943739414215088, '*61': 5.943734169006348, '*62': 5.943727970123291, '*63': 5.943719863891602, '*64': 5.943719863891602, '*65': 5.943719863891602, '*66': 5.943719863891602, '*67': 5.943719863891602, '*68': 5.943719863891602, '*69': 5.943702697753906, '*7': 5.9440226554870605, '*70': 5.943700313568115, '*71': 5.943700313568115, '*72': 5.943695545196533, '*73': 5.943671703338623, '*74': 5.943672180175781, '*75': 5.943672180175781, '*76': 5.943665504455566, '*77': 5.943638801574707, '*78': 5.943638801574707, '*79': 5.943638801574707, '*8': 5.9440226554870605, '*80': 5.943630695343018, '*81': 5.943630695343018, '*82': 5.943585395812988, '*83': 5.943577289581299, '*84': 5.943577289581299, '*85': 5.943577289581299, '*86': 5.94357442855835, '*87': 5.94357442855835, '*88': 5.94357442855835, '*89': 5.943568229675293, '*9': 5.944015026092529, '*90': 5.943568229675293, '*91': 5.943568229675293, '*92': 5.943528175354004, '*93': 5.94349479675293, '*94': 5.94349479675293, '*95': 5.943480014801025, '*96': 5.943466663360596, '*97': 5.943450450897217, '*98': 5.943450450897217, '*99': 5.943450450897217}
            # test_standard_cube.exp_bmaj_dict
            exp_bmaj_dict = {'*0': 8.501383781433105, '*1': 8.5013427734375, '*10': 8.501267433166504, '*100': 8.49976634979248, '*101': 8.499771118164062, '*102': 8.499725341796875, '*103': 8.499642372131348, '*104': 8.499631881713867, '*105': 8.499631881713867, '*106': 8.49958610534668, '*107': 8.49958610534668, '*108': 8.499568939208984, '*109': 8.499568939208984, '*11': 8.50125789642334, '*110': 8.499568939208984, '*111': 8.499568939208984, '*112': 8.499568939208984, '*113': 8.499568939208984, '*114': 8.499568939208984, '*115': 8.499568939208984, '*116': 8.499568939208984, '*117': 8.499425888061523, '*118': 8.499419212341309, '*119': 8.499419212341309, '*12': 8.501240730285645, '*120': 8.49941635131836, '*121': 8.49941635131836, '*122': 8.49940013885498, '*123': 8.499324798583984, '*124': 8.499231338500977, '*125': 8.499188423156738, '*126': 8.499133110046387, '*127': 8.499133110046387, '*128': 8.499133110046387, '*129': 8.499103546142578, '*13': 8.501240730285645, '*130': 8.499103546142578, '*131': 8.499103546142578, '*132': 8.499103546142578, '*133': 8.499098777770996, '*134': 8.499098777770996, '*135': 8.49912166595459, '*136': 8.49912166595459, '*137': 8.499106407165527, '*138': 8.498995780944824, '*139': 8.498995780944824, '*14': 8.501240730285645, '*140': 8.499002456665039, '*141': 8.499002456665039, '*142': 8.499002456665039, '*143': 8.499002456665039, '*144': 8.499002456665039, '*145': 8.499002456665039, '*146': 8.498955726623535, '*147': 8.498955726623535, '*148': 8.49891185760498, '*149': 8.498912811279297, '*15': 8.501240730285645, '*150': 8.498912811279297, '*151': 8.498912811279297, '*152': 8.498917579650879, '*153': 8.498917579650879, '*154': 8.498931884765625, '*155': 8.498931884765625, '*156': 8.49893569946289, '*157': 8.49893569946289, '*158': 8.49893569946289, '*159': 8.49893569946289, '*16': 8.501228332519531, '*160': 8.49893569946289, '*161': 8.498952865600586, '*162': 8.498952865600586, '*163': 8.498932838439941, '*164': 8.498906135559082, '*165': 8.498762130737305, '*166': 8.498750686645508, '*167': 8.498750686645508, '*168': 8.498736381530762, '*169': 8.498736381530762, '*17': 8.501107215881348, '*170': 8.498736381530762, '*171': 8.498736381530762, '*172': 8.498745918273926, '*173': 8.498677253723145, '*174': 8.498661994934082, '*175': 8.498669624328613, '*176': 8.498653411865234, '*177': 8.49864387512207, '*178': 8.498629570007324, '*179': 8.498629570007324, '*18': 8.50111198425293, '*180': 8.498604774475098, '*181': 8.498604774475098, '*182': 8.498604774475098, '*183': 8.498589515686035, '*184': 8.498554229736328, '*185': 8.498554229736328, '*186': 8.498551368713379, '*187': 8.498551368713379, '*188': 8.498551368713379, '*189': 8.498551368713379, '*19': 8.501042366027832, '*190': 8.498551368713379, '*191': 8.498551368713379, '*192': 8.498496055603027, '*193': 8.498496055603027, '*194': 8.498496055603027, '*195': 8.498496055603027, '*196': 8.498496055603027, '*197': 8.498496055603027, '*198': 8.498467445373535, '*199': 8.49845027923584, '*2': 8.501371383666992, '*20': 8.501042366027832, '*200': 8.49845027923584, '*201': 8.49845027923584, '*202': 8.498435020446777, '*203': 8.498387336730957, '*204': 8.498392105102539, '*205': 8.498392105102539, '*206': 8.498392105102539, '*207': 8.498392105102539, '*208': 8.498392105102539, '*209': 8.498392105102539, '*21': 8.501042366027832, '*210': 8.498343467712402, '*211': 8.498343467712402, '*212': 8.498343467712402, '*213': 8.498343467712402, '*214': 8.498343467712402, '*215': 8.498343467712402, '*216': 8.49832534790039, '*217': 8.49832534790039, '*218': 8.498236656188965, '*219': 8.498236656188965, '*22': 8.501042366027832, '*220': 8.498186111450195, '*221': 8.498186111450195, '*222': 8.498186111450195, '*223': 8.498186111450195, '*224': 8.49821949005127, '*225': 8.498163223266602, '*226': 8.498163223266602, '*227': 8.498127937316895, '*228': 8.498111724853516, '*229': 8.498111724853516, '*23': 8.501042366027832, '*230': 8.498111724853516, '*231': 8.498093605041504, '*232': 8.498093605041504, '*233': 8.49803638458252, '*234': 8.49803638458252, '*235': 8.498059272766113, '*236': 8.498059272766113, '*237': 8.498035430908203, '*238': 8.498035430908203, '*239': 8.498000144958496, '*24': 8.500977516174316, '*240': 8.497992515563965, '*241': 8.497992515563965, '*242': 8.497992515563965, '*243': 8.498001098632812, '*244': 8.497981071472168, '*245': 8.497981071472168, '*246': 8.497981071472168, '*247': 8.497981071472168, '*248': 8.497981071472168, '*249': 8.497981071472168, '*25': 8.500937461853027, '*250': 8.497970581054688, '*251': 8.497970581054688, '*252': 8.497958183288574, '*253': 8.497900009155273, '*254': 8.497900009155273, '*255': 8.497868537902832, '*256': 8.497882843017578, '*257': 8.497847557067871, '*258': 8.497847557067871, '*259': 8.497847557067871, '*26': 8.500937461853027, '*260': 8.497847557067871, '*261': 8.497831344604492, '*262': 8.497831344604492, '*263': 8.497831344604492, '*264': 8.497831344604492, '*265': 8.497831344604492, '*266': 8.49786376953125, '*267': 8.497773170471191, '*268': 8.497733116149902, '*269': 8.497733116149902, '*27': 8.500937461853027, '*270': 8.4977445602417, '*271': 8.497739791870117, '*272': 8.497739791870117, '*273': 8.497754096984863, '*274': 8.497754096984863, '*275': 8.497754096984863, '*276': 8.497712135314941, '*277': 8.497693061828613, '*278': 8.497693061828613, '*279': 8.497665405273438, '*28': 8.500882148742676, '*280': 8.497665405273438, '*281': 8.497665405273438, '*282': 8.497665405273438, '*283': 8.49767017364502, '*284': 8.49767017364502, '*285': 8.49767017364502, '*286': 8.49767017364502, '*287': 8.49767017364502, '*288': 8.497586250305176, '*289': 8.497503280639648, '*29': 8.500882148742676, '*290': 8.497503280639648, '*291': 8.4975004196167, '*292': 8.497406005859375, '*293': 8.497406005859375, '*294': 8.497406005859375, '*295': 8.497406005859375, '*296': 8.497410774230957, '*297': 8.497359275817871, '*298': 8.497355461120605, '*299': 8.497365951538086, '*3': 8.501371383666992, '*30': 8.500882148742676, '*300': 8.49736213684082, '*301': 8.49736213684082, '*302': 8.49730396270752, '*303': 8.497294425964355, '*304': 8.497294425964355, '*305': 8.497184753417969, '*306': 8.497184753417969, '*307': 8.497184753417969, '*308': 8.497184753417969, '*309': 8.497184753417969, '*31': 8.500882148742676, '*310': 8.497190475463867, '*311': 8.497152328491211, '*312': 8.49712085723877, '*313': 8.49712085723877, '*314': 8.497085571289062, '*315': 8.497085571289062, '*316': 8.497085571289062, '*317': 8.497085571289062, '*318': 8.497077941894531, '*319': 8.497077941894531, '*32': 8.500880241394043, '*320': 8.49706745147705, '*321': 8.49706745147705, '*322': 8.49706745147705, '*323': 8.496978759765625, '*324': 8.496978759765625, '*325': 8.496988296508789, '*326': 8.49695873260498, '*327': 8.49695873260498, '*328': 8.49695873260498, '*329': 8.49695873260498, '*33': 8.500880241394043, '*330': 8.49695873260498, '*331': 8.49695873260498, '*332': 8.49695873260498, '*333': 8.496955871582031, '*334': 8.496955871582031, '*335': 8.496960639953613, '*336': 8.496960639953613, '*337': 8.496960639953613, '*338': 8.496953964233398, '*339': 8.496953964233398, '*34': 8.500880241394043, '*340': 8.496964454650879, '*341': 8.496933937072754, '*342': 8.496933937072754, '*343': 8.496898651123047, '*344': 8.49687671661377, '*345': 8.496785163879395, '*346': 8.496785163879395, '*347': 8.496785163879395, '*348': 8.496785163879395, '*349': 8.496785163879395, '*35': 8.500884056091309, '*350': 8.496785163879395, '*351': 8.49678897857666, '*352': 8.496676445007324, '*353': 8.496676445007324, '*354': 8.496588706970215, '*355': 8.496588706970215, '*356': 8.496580123901367, '*357': 8.496580123901367, '*358': 8.496580123901367, '*359': 8.496580123901367, '*36': 8.500884056091309, '*360': 8.496580123901367, '*361': 8.49642276763916, '*362': 8.49642276763916, '*363': 8.49642276763916, '*364': 8.496406555175781, '*365': 8.496430397033691, '*366': 8.496448516845703, '*367': 8.496448516845703, '*368': 8.496448516845703, '*369': 8.496413230895996, '*37': 8.500884056091309, '*370': 8.49639892578125, '*371': 8.49639892578125, '*372': 8.49639892578125, '*373': 8.49639892578125, '*374': 8.49639892578125, '*375': 8.496352195739746, '*376': 8.496352195739746, '*377': 8.496352195739746, '*378': 8.496352195739746, '*379': 8.496352195739746, '*38': 8.500860214233398, '*380': 8.496352195739746, '*381': 8.496352195739746, '*382': 8.496352195739746, '*383': 8.496353149414062, '*384': 8.496353149414062, '*385': 8.496321678161621, '*386': 8.496297836303711, '*387': 8.496237754821777, '*388': 8.496237754821777, '*389': 8.49616813659668, '*39': 8.500860214233398, '*390': 8.49616813659668, '*391': 8.4961576461792, '*392': 8.4961576461792, '*393': 8.4961576461792, '*394': 8.49609375, '*395': 8.496075630187988, '*396': 8.496075630187988, '*397': 8.496075630187988, '*398': 8.496075630187988, '*399': 8.496075630187988, '*4': 8.501371383666992, '*40': 8.500860214233398, '*400': 8.496064186096191, '*401': 8.496079444885254, '*402': 8.496079444885254, '*403': 8.496061325073242, '*404': 8.49602222442627, '*405': 8.49592113494873, '*406': 8.49592113494873, '*407': 8.495946884155273, '*408': 8.495946884155273, '*409': 8.495946884155273, '*41': 8.500860214233398, '*410': 8.495942115783691, '*411': 8.495942115783691, '*412': 8.495942115783691, '*413': 8.495942115783691, '*414': 8.495928764343262, '*415': 8.495914459228516, '*416': 8.495914459228516, '*417': 8.495914459228516, '*418': 8.495862007141113, '*419': 8.495800971984863, '*42': 8.500804901123047, '*420': 8.495829582214355, '*421': 8.495829582214355, '*422': 8.495829582214355, '*423': 8.495829582214355, '*424': 8.495783805847168, '*425': 8.495783805847168, '*426': 8.495759010314941, '*427': 8.495759010314941, '*428': 8.495759010314941, '*429': 8.495759010314941, '*43': 8.500804901123047, '*430': 8.495759010314941, '*431': 8.495683670043945, '*432': 8.495683670043945, '*433': 8.495683670043945, '*434': 8.495683670043945, '*435': 8.49563980102539, '*436': 8.49563980102539, '*437': 8.49563980102539, '*438': 8.49563980102539, '*439': 8.49563980102539, '*44': 8.500791549682617, '*440': 8.49563980102539, '*441': 8.495579719543457, '*442': 8.495579719543457, '*443': 8.495579719543457, '*444': 8.495579719543457, '*445': 8.495579719543457, '*446': 8.495579719543457, '*447': 8.495579719543457, '*448': 8.495579719543457, '*449': 8.495530128479004, '*45': 8.500797271728516, '*450': 8.495490074157715, '*451': 8.495490074157715, '*452': 8.495403289794922, '*453': 8.495404243469238, '*454': 8.495404243469238, '*455': 8.495404243469238, '*456': 8.495386123657227, '*457': 8.495386123657227, '*458': 8.495345115661621, '*459': 8.495296478271484, '*46': 8.50074291229248, '*460': 8.495296478271484, '*461': 8.495307922363281, '*462': 8.495274543762207, '*463': 8.495274543762207, '*464': 8.495274543762207, '*465': 8.495259284973145, '*466': 8.495259284973145, '*467': 8.495259284973145, '*468': 8.495259284973145, '*469': 8.495259284973145, '*47': 8.50074291229248, '*470': 8.495259284973145, '*471': 8.495250701904297, '*472': 8.495250701904297, '*473': 8.495205879211426, '*474': 8.495189666748047, '*475': 8.495183944702148, '*476': 8.495180130004883, '*477': 8.495180130004883, '*478': 8.495171546936035, '*479': 8.495190620422363, '*48': 8.50074291229248, '*480': 8.495190620422363, '*481': 8.495190620422363, '*482': 8.495136260986328, '*483': 8.495136260986328, '*484': 8.495136260986328, '*485': 8.495139122009277, '*486': 8.495139122009277, '*487': 8.495139122009277, '*488': 8.495139122009277, '*489': 8.495139122009277, '*49': 8.50074291229248, '*490': 8.495108604431152, '*491': 8.495108604431152, '*492': 8.495108604431152, '*493': 8.495108604431152, '*494': 8.495108604431152, '*495': 8.495108604431152, '*496': 8.495108604431152, '*497': 8.495058059692383, '*498': 8.495058059692383, '*499': 8.495058059692383, '*5': 8.501351356506348, '*50': 8.500726699829102, '*500': 8.49502944946289, '*501': 8.49502944946289, '*502': 8.49502944946289, '*503': 8.495006561279297, '*504': 8.495006561279297, '*505': 8.494993209838867, '*506': 8.494993209838867, '*507': 8.494993209838867, '*51': 8.500726699829102, '*52': 8.50068473815918, '*53': 8.50068473815918, '*54': 8.50068473815918, '*55': 8.50068473815918, '*56': 8.500676155090332, '*57': 8.500635147094727, '*58': 8.500604629516602, '*59': 8.500604629516602, '*6': 8.501274108886719, '*60': 8.500577926635742, '*61': 8.500565528869629, '*62': 8.50053882598877, '*63': 8.500545501708984, '*64': 8.500545501708984, '*65': 8.500545501708984, '*66': 8.500545501708984, '*67': 8.500545501708984, '*68': 8.500545501708984, '*69': 8.500402450561523, '*7': 8.501262664794922, '*70': 8.500380516052246, '*71': 8.500380516052246, '*72': 8.500340461730957, '*73': 8.500334739685059, '*74': 8.500287055969238, '*75': 8.500287055969238, '*76': 8.500186920166016, '*77': 8.500128746032715, '*78': 8.500128746032715, '*79': 8.500128746032715, '*8': 8.501262664794922, '*80': 8.500129699707031, '*81': 8.500129699707031, '*82': 8.500100135803223, '*83': 8.500035285949707, '*84': 8.500035285949707, '*85': 8.500035285949707, '*86': 8.499926567077637, '*87': 8.499926567077637, '*88': 8.499926567077637, '*89': 8.499892234802246, '*9': 8.501267433166504, '*90': 8.499892234802246, '*91': 8.499892234802246, '*92': 8.499894142150879, '*93': 8.499887466430664, '*94': 8.499887466430664, '*95': 8.499871253967285, '*96': 8.499770164489746, '*97': 8.499777793884277, '*98': 8.499777793884277, '*99': 8.499777793884277}
            # test_standard_cube.exp_pa_dict
            exp_pa_dict = {'*0': 72.54618072509766, '*1': 72.5458984375, '*10': 72.54621124267578, '*100': 72.54096984863281, '*101': 72.54096221923828, '*102': 72.54052734375, '*103': 72.54006958007812, '*104': 72.54045867919922, '*105': 72.54045867919922, '*106': 72.54037475585938, '*107': 72.54037475585938, '*108': 72.54067993164062, '*109': 72.54067993164062, '*11': 72.54644775390625, '*110': 72.54067993164062, '*111': 72.54067993164062, '*112': 72.54067993164062, '*113': 72.54067993164062, '*114': 72.54067993164062, '*115': 72.54067993164062, '*116': 72.54067993164062, '*117': 72.54016876220703, '*118': 72.54037475585938, '*119': 72.54037475585938, '*12': 72.547119140625, '*120': 72.54056549072266, '*121': 72.54056549072266, '*122': 72.54084777832031, '*123': 72.54094696044922, '*124': 72.54084014892578, '*125': 72.54124450683594, '*126': 72.54103088378906, '*127': 72.54103088378906, '*128': 72.54103088378906, '*129': 72.5400161743164, '*13': 72.547119140625, '*130': 72.5400161743164, '*131': 72.5400161743164, '*132': 72.5400161743164, '*133': 72.54019927978516, '*134': 72.54019927978516, '*135': 72.54015350341797, '*136': 72.54015350341797, '*137': 72.53997039794922, '*138': 72.54032897949219, '*139': 72.54032897949219, '*14': 72.547119140625, '*140': 72.54021453857422, '*141': 72.54021453857422, '*142': 72.54021453857422, '*143': 72.54021453857422, '*144': 72.54021453857422, '*145': 72.54021453857422, '*146': 72.53984069824219, '*147': 72.53984069824219, '*148': 72.53958129882812, '*149': 72.53987884521484, '*15': 72.547119140625, '*150': 72.53987884521484, '*151': 72.53987884521484, '*152': 72.5394287109375, '*153': 72.5394287109375, '*154': 72.53932189941406, '*155': 72.53932189941406, '*156': 72.53958892822266, '*157': 72.53958892822266, '*158': 72.53958892822266, '*159': 72.53958892822266, '*16': 72.54747009277344, '*160': 72.53958892822266, '*161': 72.53949737548828, '*162': 72.53949737548828, '*163': 72.54027557373047, '*164': 72.54119110107422, '*165': 72.54067993164062, '*166': 72.54100036621094, '*167': 72.54100036621094, '*168': 72.54080963134766, '*169': 72.54080963134766, '*17': 72.54749298095703, '*170': 72.54080963134766, '*171': 72.54080963134766, '*172': 72.53984832763672, '*173': 72.53959655761719, '*174': 72.53941345214844, '*175': 72.53883361816406, '*176': 72.53913116455078, '*177': 72.53937530517578, '*178': 72.53985595703125, '*179': 72.53985595703125, '*18': 72.54743957519531, '*180': 72.54026794433594, '*181': 72.54026794433594, '*182': 72.54026794433594, '*183': 72.54010009765625, '*184': 72.54019165039062, '*185': 72.54019165039062, '*186': 72.54039001464844, '*187': 72.54039001464844, '*188': 72.54039001464844, '*189': 72.54039001464844, '*19': 72.54718780517578, '*190': 72.54039001464844, '*191': 72.54039001464844, '*192': 72.54012298583984, '*193': 72.54012298583984, '*194': 72.54012298583984, '*195': 72.54012298583984, '*196': 72.54012298583984, '*197': 72.54012298583984, '*198': 72.54006958007812, '*199': 72.54035949707031, '*2': 72.54558563232422, '*20': 72.54718780517578, '*200': 72.54035949707031, '*201': 72.54035949707031, '*202': 72.54064178466797, '*203': 72.54032135009766, '*204': 72.54026794433594, '*205': 72.54026794433594, '*206': 72.54026794433594, '*207': 72.54026794433594, '*208': 72.54026794433594, '*209': 72.54026794433594, '*21': 72.54718780517578, '*210': 72.54014587402344, '*211': 72.54014587402344, '*212': 72.54014587402344, '*213': 72.54014587402344, '*214': 72.54014587402344, '*215': 72.54014587402344, '*216': 72.54064178466797, '*217': 72.54064178466797, '*218': 72.54141235351562, '*219': 72.54141235351562, '*22': 72.54718780517578, '*220': 72.54093933105469, '*221': 72.54093933105469, '*222': 72.54093933105469, '*223': 72.54093933105469, '*224': 72.54061889648438, '*225': 72.54041290283203, '*226': 72.54041290283203, '*227': 72.54000854492188, '*228': 72.54046630859375, '*229': 72.54046630859375, '*23': 72.54718780517578, '*230': 72.54046630859375, '*231': 72.54094696044922, '*232': 72.54094696044922, '*233': 72.54092407226562, '*234': 72.54092407226562, '*235': 72.5408935546875, '*236': 72.5408935546875, '*237': 72.54131317138672, '*238': 72.54131317138672, '*239': 72.54127502441406, '*24': 72.54695129394531, '*240': 72.54150390625, '*241': 72.54150390625, '*242': 72.54150390625, '*243': 72.54146575927734, '*244': 72.54193115234375, '*245': 72.54193115234375, '*246': 72.54193115234375, '*247': 72.54193115234375, '*248': 72.54193115234375, '*249': 72.54193115234375, '*25': 72.54627990722656, '*250': 72.54222869873047, '*251': 72.54222869873047, '*252': 72.54244995117188, '*253': 72.54251098632812, '*254': 72.54251098632812, '*255': 72.5427474975586, '*256': 72.5426254272461, '*257': 72.54222106933594, '*258': 72.54222106933594, '*259': 72.54222106933594, '*26': 72.54627990722656, '*260': 72.54222106933594, '*261': 72.54263305664062, '*262': 72.54263305664062, '*263': 72.54263305664062, '*264': 72.54263305664062, '*265': 72.54263305664062, '*266': 72.54269409179688, '*267': 72.54336547851562, '*268': 72.54347229003906, '*269': 72.54347229003906, '*27': 72.54627990722656, '*270': 72.54338836669922, '*271': 72.5435791015625, '*272': 72.5435791015625, '*273': 72.54358673095703, '*274': 72.54358673095703, '*275': 72.54358673095703, '*276': 72.54357147216797, '*277': 72.54421997070312, '*278': 72.54421997070312, '*279': 72.54415893554688, '*28': 72.54607391357422, '*280': 72.54415893554688, '*281': 72.54415893554688, '*282': 72.54415893554688, '*283': 72.54443359375, '*284': 72.54443359375, '*285': 72.54443359375, '*286': 72.54443359375, '*287': 72.54443359375, '*288': 72.54430389404297, '*289': 72.54385375976562, '*29': 72.54607391357422, '*290': 72.54385375976562, '*291': 72.54399108886719, '*292': 72.54380798339844, '*293': 72.54380798339844, '*294': 72.54380798339844, '*295': 72.54380798339844, '*296': 72.54374694824219, '*297': 72.54345703125, '*298': 72.54364776611328, '*299': 72.54367065429688, '*3': 72.54558563232422, '*30': 72.54607391357422, '*300': 72.54383087158203, '*301': 72.54383087158203, '*302': 72.54360961914062, '*303': 72.54415130615234, '*304': 72.54415130615234, '*305': 72.54395294189453, '*306': 72.54395294189453, '*307': 72.54395294189453, '*308': 72.54395294189453, '*309': 72.54395294189453, '*31': 72.54607391357422, '*310': 72.54415893554688, '*311': 72.54348754882812, '*312': 72.54342651367188, '*313': 72.54342651367188, '*314': 72.54335021972656, '*315': 72.54335021972656, '*316': 72.54335021972656, '*317': 72.54335021972656, '*318': 72.54357147216797, '*319': 72.54357147216797, '*32': 72.54622650146484, '*320': 72.54387664794922, '*321': 72.54387664794922, '*322': 72.54387664794922, '*323': 72.54393005371094, '*324': 72.54393005371094, '*325': 72.54389190673828, '*326': 72.54483795166016, '*327': 72.54483795166016, '*328': 72.54483795166016, '*329': 72.54483795166016, '*33': 72.54622650146484, '*330': 72.54483795166016, '*331': 72.54483795166016, '*332': 72.54483795166016, '*333': 72.54502868652344, '*334': 72.54502868652344, '*335': 72.54524993896484, '*336': 72.54524993896484, '*337': 72.54524993896484, '*338': 72.54592895507812, '*339': 72.54592895507812, '*34': 72.54622650146484, '*340': 72.5458755493164, '*341': 72.54499053955078, '*342': 72.54499053955078, '*343': 72.5450210571289, '*344': 72.54485321044922, '*345': 72.54437255859375, '*346': 72.54437255859375, '*347': 72.54437255859375, '*348': 72.54437255859375, '*349': 72.54437255859375, '*35': 72.54571533203125, '*350': 72.54437255859375, '*351': 72.54459381103516, '*352': 72.54454803466797, '*353': 72.54454803466797, '*354': 72.54415130615234, '*355': 72.54415130615234, '*356': 72.54439544677734, '*357': 72.54439544677734, '*358': 72.54439544677734, '*359': 72.54439544677734, '*36': 72.54571533203125, '*360': 72.54439544677734, '*361': 72.54432678222656, '*362': 72.54432678222656, '*363': 72.54432678222656, '*364': 72.54479217529297, '*365': 72.54451751708984, '*366': 72.54434967041016, '*367': 72.54434967041016, '*368': 72.54434967041016, '*369': 72.54428100585938, '*37': 72.54571533203125, '*370': 72.5446548461914, '*371': 72.5446548461914, '*372': 72.5446548461914, '*373': 72.5446548461914, '*374': 72.5446548461914, '*375': 72.5445327758789, '*376': 72.5445327758789, '*377': 72.5445327758789, '*378': 72.5445327758789, '*379': 72.5445327758789, '*38': 72.5455322265625, '*380': 72.5445327758789, '*381': 72.5445327758789, '*382': 72.5445327758789, '*383': 72.5447769165039, '*384': 72.5447769165039, '*385': 72.54471588134766, '*386': 72.54437255859375, '*387': 72.54450225830078, '*388': 72.54450225830078, '*389': 72.54415893554688, '*39': 72.5455322265625, '*390': 72.54415893554688, '*391': 72.54446411132812, '*392': 72.54446411132812, '*393': 72.54446411132812, '*394': 72.54413604736328, '*395': 72.54474639892578, '*396': 72.54474639892578, '*397': 72.54474639892578, '*398': 72.54474639892578, '*399': 72.54474639892578, '*4': 72.54558563232422, '*40': 72.5455322265625, '*400': 72.54512786865234, '*401': 72.54499053955078, '*402': 72.54499053955078, '*403': 72.54547119140625, '*404': 72.5450210571289, '*405': 72.54446411132812, '*406': 72.54446411132812, '*407': 72.54449462890625, '*408': 72.54449462890625, '*409': 72.54449462890625, '*41': 72.5455322265625, '*410': 72.54468536376953, '*411': 72.54468536376953, '*412': 72.54468536376953, '*413': 72.54468536376953, '*414': 72.54508972167969, '*415': 72.54559326171875, '*416': 72.54559326171875, '*417': 72.54559326171875, '*418': 72.54500579833984, '*419': 72.54496765136719, '*42': 72.545654296875, '*420': 72.54496002197266, '*421': 72.54496002197266, '*422': 72.54496002197266, '*423': 72.54496002197266, '*424': 72.54521942138672, '*425': 72.54521942138672, '*426': 72.54429626464844, '*427': 72.54429626464844, '*428': 72.54429626464844, '*429': 72.54429626464844, '*43': 72.545654296875, '*430': 72.54429626464844, '*431': 72.54389953613281, '*432': 72.54389953613281, '*433': 72.54389953613281, '*434': 72.54389953613281, '*435': 72.54237365722656, '*436': 72.54237365722656, '*437': 72.54237365722656, '*438': 72.54237365722656, '*439': 72.54237365722656, '*44': 72.54613494873047, '*440': 72.54237365722656, '*441': 72.54206848144531, '*442': 72.54206848144531, '*443': 72.54206848144531, '*444': 72.54206848144531, '*445': 72.54206848144531, '*446': 72.54206848144531, '*447': 72.54206848144531, '*448': 72.54206848144531, '*449': 72.54217529296875, '*45': 72.54609680175781, '*450': 72.5418930053711, '*451': 72.5418930053711, '*452': 72.54092407226562, '*453': 72.54122161865234, '*454': 72.54122161865234, '*455': 72.54122161865234, '*456': 72.5416259765625, '*457': 72.5416259765625, '*458': 72.54183959960938, '*459': 72.5418472290039, '*46': 72.54585266113281, '*460': 72.5418472290039, '*461': 72.54179382324219, '*462': 72.54154205322266, '*463': 72.54154205322266, '*464': 72.54154205322266, '*465': 72.54129028320312, '*466': 72.54129028320312, '*467': 72.54129028320312, '*468': 72.54129028320312, '*469': 72.54129028320312, '*47': 72.54585266113281, '*470': 72.54129028320312, '*471': 72.5415267944336, '*472': 72.5415267944336, '*473': 72.54109191894531, '*474': 72.54157257080078, '*475': 72.54180145263672, '*476': 72.5419921875, '*477': 72.5419921875, '*478': 72.54230499267578, '*479': 72.54096984863281, '*48': 72.54585266113281, '*480': 72.54096984863281, '*481': 72.54096984863281, '*482': 72.54082489013672, '*483': 72.54082489013672, '*484': 72.54082489013672, '*485': 72.54109191894531, '*486': 72.54109191894531, '*487': 72.54109191894531, '*488': 72.54109191894531, '*489': 72.54109191894531, '*49': 72.54585266113281, '*490': 72.54019927978516, '*491': 72.54019927978516, '*492': 72.54019927978516, '*493': 72.54019927978516, '*494': 72.54019927978516, '*495': 72.54019927978516, '*496': 72.54019927978516, '*497': 72.53990936279297, '*498': 72.53990936279297, '*499': 72.53990936279297, '*5': 72.54608154296875, '*50': 72.546142578125, '*500': 72.54069519042969, '*501': 72.54069519042969, '*502': 72.54069519042969, '*503': 72.5416488647461, '*504': 72.5416488647461, '*505': 72.54203033447266, '*506': 72.54203033447266, '*507': 72.54203033447266, '*51': 72.546142578125, '*52': 72.5461196899414, '*53': 72.5461196899414, '*54': 72.5461196899414, '*55': 72.5461196899414, '*56': 72.54659271240234, '*57': 72.54618835449219, '*58': 72.54581451416016, '*59': 72.54581451416016, '*6': 72.54618072509766, '*60': 72.54548645019531, '*61': 72.54580688476562, '*62': 72.5449447631836, '*63': 72.54444122314453, '*64': 72.54444122314453, '*65': 72.54444122314453, '*66': 72.54444122314453, '*67': 72.54444122314453, '*68': 72.54444122314453, '*69': 72.54393768310547, '*7': 72.54672241210938, '*70': 72.54375457763672, '*71': 72.54375457763672, '*72': 72.54331970214844, '*73': 72.54352569580078, '*74': 72.54351806640625, '*75': 72.54351806640625, '*76': 72.54296875, '*77': 72.5428695678711, '*78': 72.5428695678711, '*79': 72.5428695678711, '*8': 72.54672241210938, '*80': 72.54236602783203, '*81': 72.54236602783203, '*82': 72.54192352294922, '*83': 72.54167938232422, '*84': 72.54167938232422, '*85': 72.54167938232422, '*86': 72.54149627685547, '*87': 72.54149627685547, '*88': 72.54149627685547, '*89': 72.5408706665039, '*9': 72.54621124267578, '*90': 72.5408706665039, '*91': 72.5408706665039, '*92': 72.5411148071289, '*93': 72.5413818359375, '*94': 72.5413818359375, '*95': 72.54183197021484, '*96': 72.54070281982422, '*97': 72.54064178466797, '*98': 72.54064178466797, '*99': 72.54064178466797}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed=self.filter_report(report)
        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, 'test_standard_cube', self.parallel)

        test_dict['test_standard_cube']['report'] = report
        test_dict['test_standard_cube']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(image=img+'.image', range_list=[0.3, 1.0])
        self.mom8_creator(image=img+'.residual', range_list=[0.3, 1.0])
        test_dict['test_standard_cube']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_standard_cube']['images'].append(img+'.image.profile.png')

        if savemetricdict:
            ### serialize ndarray in mask_stats_dcit
            mask_stats_mod_dict = copy.deepcopy(mask_stats_dict)
            mask_stats_mod_dict['mask'] = mask_stats_dict['mask'].tolist()
            #create a nested dictionary containing exp dictionaries to save
            savedict = {}
            #list of stats to save
            # im_stats, mask_stats, pb_stats, psf_stats,\
            # model_stats, resid_stats, sumwt_stats]
            savedict['im_stats_dict']=im_stats_dict
            savedict['mask_stats_dict']=mask_stats_mod_dict
            savedict['pb_stats_dict']=pb_stats_dict
            savedict['psf_stats_dict']=psf_stats_dict
            savedict['model_stats_dict']=model_stats_dict
            savedict['resid_stats_dict']=resid_stats_dict
            savedict['sumwt_stats_dict']=sumwt_stats_dict
            
            if self.parallel:
                savedict['bmin_dict']=bmin_dict
                savedict['bmaj_dict']=bmaj_dict
                savedict['pa_dict']=pa_dict

            self.save_dict_to_file('standard_cube',savedict, 'standard_cube_current_metrics')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cube
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cube_briggsbwtaper(self):
        ''' Standard (single field) cube imaging with briggsbwtaper - central field of SMIDGE_NWCloud (field 3), spw 22 '''
        
        test_name = self._testMethodName
        file_name = 'standard_cube_briggsbwtaper.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, imagename=file_name+'0', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS'
            ' 00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526743594GHz', width='0.2441741MHz', \
            outframe='LSRK', pblimit=0.2, perchanweightdensity=True,\
            gridder='standard', chanchunks=-1, mosweight=False, \
            deconvolver='hogbom', usepointing=False, restoration=False, \
            pbcor=False, weighting='briggsbwtaper', restoringbeam='common', \
            robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
                nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
                outframe='LSRK', perchanweightdensity=True, \
                usepointing=False, pblimit=0.2, nsigma=0.0, \
                gridder='standard', chanchunks=-1, mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggsbwtaper', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none', \
                parallel=True, verbose=True)

            # retrieve per-channel beam statistics 
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(image=img+'.image')

            tclean(vis=self.msfile, imagename=file_name+'1', spw=['0'], \
                field='1', imsize=[80, 80], cell=['1.1arcsec'], \
                phasecenter='ICRS 00:45:54.3836 -073.15.29.413', stokes='I', \
                antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                specmode='cube', nchan=508, start='220.2526743594GHz', \
                width='0.2441741MHz', outframe='LSRK', \
                perchanweightdensity=True, gridder='standard', \
                chanchunks=-1, mosweight=False, usepointing=False, \
                pblimit=0.2, deconvolver='hogbom', restoration=True, \
                restoringbeam='common', pbcor=True, weighting='briggsbwtaper', \
                robust=0.5, npixels=0, niter=0, threshold='0.354Jy', \
                nsigma=0.0, interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none', \
                parallel=False, verbose=True)

        else:
            tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
                spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
                scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
                datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
                '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
                nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
                outframe='LSRK', perchanweightdensity=True, \
                usepointing=False, pblimit=0.2, nsigma=0.0, \
                gridder='standard', chanchunks=-1, mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggsbwtaper', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none', \
                restoringbeam='common', parallel=False, verbose=True)

            # retrieve per-channel beam statistics 
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(image=img+'.psf')


        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report(test_standard_cube_briggsbwtaper)
        im_stats_dict = self.image_stats(image=img+'.image', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 8.900085299026527],
            'com_bmin': [False, 6.193896355268416],
            'com_pa': [False, 67.88903485967856],
            'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0724563598632812],
            'max_val_pos': [True, [37, 36, 0, 254]],
            'min_val': [False, -0.6416710615158081],
            'min_val_pos': [True, [21, 27, 0, 163]],
            'im_rms': [False, 0.1410881918873013],
            'rms_per_chan': [False, [0.12290898085061072, 0.14257689228073617, 0.12671417592671352, 0.11172644319283263, 0.14934046219727876, 0.1233022963004575, 0.1316617075885439, 0.10687715500507994, 0.15440277261843877, 0.15923512337364315, 0.11261330252885768, 0.14019044996404106, 0.16725342914945754, 0.13082002809405885, 0.15717325857778172, 0.13676716196926497, 0.10348059645324566, 0.12580125461063632, 0.14055468413599767, 0.1760971440523076, 0.13231902566151393, 0.15615432219739606, 0.11358763822766985, 0.11163603242896791, 0.1693573388568213, 0.13672584704600213, 0.12994425281122007, 0.11473852955818571, 0.142020732342772, 0.13874814584221792, 0.13613511688061466, 0.11921928506172898, 0.12695356702893465, 0.12230085119899346, 0.15655691669204613, 0.15605048925812473, 0.15797374641067014, 0.14376616156759608, 0.14173606556109017, 0.13120065249352203, 0.1318030853954161, 0.15100264053794282, 0.1436096861820802, 0.13537857172023654, 0.13946825736059518, 0.14121042075276188, 0.12627069016144918, 0.13207764149477344, 0.16521686719234963, 0.13933029682048714, 0.12131217458425789, 0.15558850064116342, 0.13480354316441812, 0.12706503426856744, 0.11836229576232353, 0.132689723664371, 0.14928922875072054, 0.17390335924051184, 0.1092401543197615, 0.14026683218962221, 0.14748731051595246, 0.1364589517289056, 0.11889854371099778, 0.12379782754125412, 0.12539928407638745, 0.13752761848185477, 0.1273665204254799, 0.14758429600783765, 0.13281798078835313, 0.15519223821410325, 0.12790003959049917, 0.13521019720568053, 0.14849112259696942, 0.1387760086699735, 0.17177836770016536, 0.1638507135605381, 0.1176411129524291, 0.15015600238308188, 0.14643014972275534, 0.12108409399419816, 0.12374091398285429, 0.16404649343289684, 0.11272305338423198, 0.17919523358197306, 0.13279110674294275, 0.14060781659541516, 0.11297028401625303, 0.13175733566982806, 0.15864676649364556, 0.13302407490398102, 0.1481985455324882, 0.13089260233543767, 0.1372059292952956, 0.14102983233928806, 0.1220417217334515, 0.1426333709833995, 0.12368184736821582, 0.16439302084521634, 0.1459633864090159, 0.12872364003940287, 0.15385481190277425, 0.11530116401993136, 0.1399014814090996, 0.15108650078046265, 0.13163828622886412, 0.13012072641263844, 0.1552045845207469, 0.1414108649296818, 0.15066556681410845, 0.12723762438104666, 0.14272569631181775, 0.12939676717125606, 0.12646083232388844, 0.13677927294840872, 0.12962505468458058, 0.1311267359331575, 0.11422832927804552, 0.15544387820117134, 0.1263821785038708, 0.15246375483272553, 0.13229857256860048, 0.14271875167726572, 0.1365588911184279, 0.14706019722813085, 0.12233854359244917, 0.15551604126741714, 0.13949406092889993, 0.14687632039060197, 0.1284839673034992, 0.12048579865672468, 0.14386156211949575, 0.12301608111671818, 0.15633816770695802, 0.152685823688177, 0.15616188431107525, 0.12034648579669188, 0.1600003030577954, 0.13632843570617187, 0.1358342661382129, 0.14932409256263407, 0.11427223515942367, 0.14402441494552998, 0.14543940938082892, 0.11584322136769924, 0.16021733082403353, 0.12312756250651857, 0.1429147478677048, 0.16275926891647502, 0.12722449762548585, 0.1353235846302639, 0.14031805327698588, 0.10244717267795887, 0.14750901728799104, 0.12391837658241082, 0.16384804706901723, 0.14698556489729708, 0.14634371131379295, 0.12823369020590178, 0.11607963704471314, 0.14646514573566127, 0.11859529977750853, 0.1533009033735843, 0.15336123085541178, 0.17443566969044008, 0.1555579778464484, 0.13826371457422384, 0.11937072266086565, 0.12316574686386468, 0.12534163692198463, 0.1398177372337716, 0.1365838635984754, 0.12786115306289006, 0.13358163640872325, 0.16595092356031152, 0.12566456317623703, 0.13374773305225734, 0.18307306969371642, 0.14357620956527756, 0.15029078821943345, 0.1441831126912868, 0.1286293264882945, 0.13796961831311216, 0.13655707738120842, 0.14898989652043723, 0.10487188237865079, 0.12110149636336076, 0.1732258165345567, 0.1258294635110184, 0.1293292682830415, 0.14835522460567896, 0.11714215623512433, 0.12336058127049529, 0.14850957290685526, 0.16245145333354116, 0.12567893517161385, 0.13945285055643794, 0.1326859978378357, 0.1501689237016373, 0.1691912274339341, 0.15546082480381557, 0.1477211560178025, 0.16453661789524993, 0.1296846736678987, 0.15113428267252244, 0.16287548117997902, 0.12530420684085453, 0.13715954980291806, 0.126287655972795, 0.12477887384768895, 0.1328796002994287, 0.15465673064748145, 0.152248855597909, 0.14567271683798488, 0.1447253085311366, 0.14718154934940425, 0.12002317412259286, 0.11994038184848867, 0.1225248444063849, 0.1606576330330884, 0.1575838334025464, 0.14794883712334206, 0.1198054930530698, 0.13093245582330149, 0.1193792195905996, 0.12872417714115242, 0.15919141537195755, 0.18863723322738327, 0.121344016789415, 0.12774182706619602, 0.1324124143853652, 0.1587027875381588, 0.1707034109763738, 0.13406855193238387, 0.16393863526552424, 0.14378843863622814, 0.1276297098368369, 0.12008561845486326, 0.14780126749090602, 0.15817158027819092, 0.1786780394222075, 0.13726383233839543, 0.12569766757859185, 0.1596675592417091, 0.16200854176505602, 0.12889566050040083, 0.15395982551893328, 0.12398253328925761, 0.16334089334278282, 0.1573077882947232, 0.170278289772687, 0.1666365467136698, 0.1671250251718543, 0.20094385619045316, 0.20868728889081525, 0.18798507042676246, 0.20115970102455497, 0.1969195747003099, 0.16750740331618408, 0.13464917184693717, 0.12186140004433678, 0.13016520313062654, 0.15062514933380772, 0.13319730810106403, 0.12619178774612666, 0.15374280490389886, 0.1648069702400786, 0.15804410763272553, 0.14341178171965052, 0.1189070708644023, 0.13209812168686544, 0.13882846984218505, 0.13208821168106688, 0.13791261074756286, 0.1453722035587563, 0.13825645591288435, 0.1516176969593033, 0.13011226211946536, 0.14358664748656805, 0.1303609224393637, 0.15077768552902335, 0.11141425054556503, 0.12638892584818023, 0.14539209390111735, 0.12672085174296246, 0.1376905126761013, 0.11847701938516829, 0.1493730341986286, 0.11881681796106171, 0.11524070672956355, 0.13626068268069916, 0.1599080181716221, 0.14914461232679782, 0.12046711164162706, 0.13358857007614966, 0.1427676196526599, 0.13022404399407686, 0.17035263383318347, 0.1490074914838121, 0.1108352158203021, 0.10605700116090709, 0.11644847544054283, 0.14747264605319355, 0.14798831629351297, 0.13857134855887465, 0.1254343257457443, 0.13533826261335155, 0.12349769953995753, 0.12369173353800561, 0.14159323775793625, 0.16141269592176555, 0.15027150207838086, 0.13528185931436473, 0.17177910100841462, 0.16826814945487623, 0.15318178780409825, 0.12352182208401302, 0.14072943729028406, 0.13555999358770096, 0.1318537126585016, 0.14300431721908874, 0.1338096674083135, 0.15435867612490678, 0.13382476750670413, 0.14071234769239446, 0.1298025807991066, 0.13367157152957698, 0.1360502830532836, 0.1286917199710627, 0.13588835476750216, 0.12912777669894607, 0.15391289894940063, 0.12971942700508174, 0.1424283386465245, 0.15262061324904286, 0.13986033561204655, 0.1585869329264211, 0.1306428353030295, 0.11774463633838317, 0.15687562694714619, 0.13970948459254875, 0.14744333107759136, 0.1629349106803633, 0.13554728893996187, 0.11500732988026163, 0.15541370308713684, 0.13442487800956282, 0.14279933203008277, 0.11841773858040139, 0.10870069867018264, 0.14002291416599788, 0.1360447208654598, 0.124799034714854, 0.14762145655003475, 0.12934922185777842, 0.14655087036345735, 0.1546010202478706, 0.1445847521700007, 0.12788201616824202, 0.13186944260510575, 0.13296541840724005, 0.1602869274805873, 0.13716135160222698, 0.15556461028481286, 0.1456293758406546, 0.13242711567023355, 0.12797086053543882, 0.11551070966785679, 0.12243135390015093, 0.14637372710341348, 0.1343558614100962, 0.1598336049756603, 0.12471996121575767, 0.1288330843188406, 0.15862914261693611, 0.14674188655158835, 0.12650161852733688, 0.14112542813017817, 0.1371463889530532, 0.13878559149751493, 0.1302282646206918, 0.15417551030646207, 0.16222337296089528, 0.1303365009634827, 0.1331633473842199, 0.1388438852458591, 0.13287020966561522, 0.16439303231555155, 0.15056164511677264, 0.14024125440378682, 0.15720272106754057, 0.15355559536592728, 0.12804250001584644, 0.1406720940126104, 0.1359784432597291, 0.17919362636908276, 0.11627777766964895, 0.15136889746768192, 0.1201429059678373, 0.13310125549985022, 0.16090526540064085, 0.11309208550065801, 0.1372077803261786, 0.12212380922432872, 0.16534890768313149, 0.1270983494080957, 0.146260444784059, 0.13556246334582084, 0.16565031211054604, 0.13163678342508658, 0.1252684586299344, 0.16403250616956142, 0.15260359126111875, 0.14671952083405962, 0.13897185124528227, 0.14236336302602126, 0.13334585463195753, 0.13148685127229875, 0.16010693768676212, 0.14827249752238292, 0.12537875386133732, 0.13996343789522653, 0.15618721009084183, 0.14308194298768598, 0.12565572065551572, 0.16416382294090787, 0.12102725272730878, 0.15186258234812194, 0.11845313663163846, 0.1371381318127285, 0.1444460308107884, 0.14149870606041565, 0.12147836476241079, 0.14336353007862318, 0.13015223073945192, 0.14343134477445169, 0.14372620927968494, 0.11510270396237209, 0.15572844214257955, 0.13744211155444913, 0.14850170800525142, 0.15293180424058142, 0.14070157655499105, 0.1541984766352866, 0.1366684600216604, 0.13587624294721937, 0.11526388579589063, 0.149029161759932, 0.14188139754628992, 0.14052586968984704, 0.17374781610176074, 0.155699982709309, 0.14095848969718072, 0.14039932888682585, 0.12154622913633595, 0.11514483682566516, 0.15129327159674422, 0.13260024291190603, 0.1477917658545152, 0.13725501560400108, 0.17391222231852652, 0.13358699895741738, 0.12303280962037541, 0.16628036307733238, 0.12147658769756359, 0.14311966593701977, 0.1476113304007078, 0.12478963514772516, 0.14661274214914224, 0.14050340456286925, 0.1435947413069879, 0.13315922878994163, 0.17001687448374847, 0.13493557930135314, 0.12027268085711552, 0.16151161982704132, 0.10983910083510932, 0.1352233697855893, 0.15608125796902717, 0.14116758885571398, 0.12064770443558666, 0.12057869445938094, 0.116113527256088, 0.12300090503018894, 0.1560255685793458, 0.11985514861584208, 0.13515511994409238, 0.1288527186296108, 0.14585408735783734, 0.12873029575008482, 0.13509806052115986, 0.14397420723151186, 0.1411492779949668, 0.14952372799269487, 0.1346989773384809, 0.14172234404095807, 0.14433992870701745, 0.11381417746312238, 0.16169674364451964, 0.14265326416311733, 0.1278435811202039, 0.1498659819374808, 0.16240046874457875, 0.13512650778894725, 0.13636049978790046, 0.15836280300910582, 0.15134991456240746, 0.13302795273087722, 0.1347110855490138]],
            'im_sum': [False, -56.868862802253915],
            'regn_sum': [False, 73.90656678751111],
            'npts_real': [True, 3251200],
            'profile': [False, 1.0035623595674308],
            'fit': [False, [0.9973981504909871, 14.040354685163061, 7.173426330542138]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [37.979806295018136, 36.942619899930506]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 70, 0, 0]), \
                (img+'.image', False, [40, 71, 0, 0]), \
                (img+'.image', True, [10, 40, 0, 0]), \
                (img+'.image', False, [9, 40, 0, 0])])

        # .image report
        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(image=img+'.mask')

        # test_standard_cube_briggsbwtaper.exp_mask_stats
        exp_mask_stats = {'npts': [True, 3251200],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'mask_pix': [False, 472],
            'mask_regns': [True, 1],
            'npts_real': [True, 3251200]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_mask_stats
        exp_pb_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20089669525623322],
            'im_rms': [False, 0.5782380936695958],
            'npts_0.2': [False, [2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997, 2997]],
            'npts_0.5': [False, [1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449, 1449]],
            'npts_real': [True, 3251200],
            'fit': [False, [1.0308176341918283, 46.61782691067407, 46.61233763883]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [40.00020058175598, 39.999935258038846]]}
                                                     

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(image=img+'.psf', fit_region = \
            'ellipse[[11.47648725deg, -73.25812003deg], [8.0291arcsec, 6.8080arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.19239011406898499],
            'min_val_pos': [True, [6, 13, 0, 6]],
            'im_rms': [False, 0.1368833033007673],
            'im_sum': [False, 6479.367987492627],
            'npts_real': [True, 3251200],
            'fit_0': [False, [1.0853206916188203, 8.11144818181605, 5.41751982967181]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [39.999047037790554, 39.992833306432374]],
            'fit_1': [False, [1.0853469585284308, 8.108788804262787, 5.415833699352728]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [39.99904889075988, 39.992845525585395]],
            'fit_2': [False, [1.0854080380567832, 8.10550076980935, 5.413890498307579]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [39.999051437972305, 39.99285529326793]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.9296272993087769],
            'max_val_pos': [True, [42, 43, 0, 256]],
            'min_val': [False, -0.6416710615158081],
            'min_val_pos': [True, [21, 27, 0, 163]],
            'im_rms': [False, 0.14099844382779275],
            'im_sum': [False, -108.81305603456369],
            'regn_sum': [False, 26.929447835078463],
            'npts_real': [True, 3251200]}


        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.3030758798122406],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, -0.03826450929045677],
            'min_val_pos': [True, [50, 42, 0, 254]],
            'im_rms': [False, 0.0002689667846165681],
            'im_sum': [False, 1.0062360279262066],
            'regn_sum': [False, 1.0445005372166634],
            'mask_non0': [True, 0],
            'npts_real': [True, 3251200]}


        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        # test_standard_cube_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 100.40980529785156],
            'max_val_pos': [True, [0, 0, 0, 152]],
            'min_val': [False, 100.40977478027344],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 100.40979754649751],
            'npts_real': [True, 508]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube_briggsbwtaper.exp_bmin_dict
            exp_bmin_dict = {"*0": 6.187705993652344, "*1": 6.187703609466553, "*10": 6.187617778778076, "*100": 6.187054634094238, "*101": 6.187032699584961, "*102": 6.187026500701904, "*103": 6.187020301818848, "*104": 6.187004566192627, "*105": 6.187004566192627, "*106": 6.186976432800293, "*107": 6.186976432800293, "*108": 6.186949253082275, "*109": 6.186949253082275, "*11": 6.187597751617432, "*110": 6.186949253082275, "*111": 6.186949253082275, "*112": 6.186949253082275, "*113": 6.186949253082275, "*114": 6.186949253082275, "*115": 6.186949253082275, "*116": 6.186949253082275, "*117": 6.186944007873535, "*118": 6.186914920806885, "*119": 6.186914920806885, "*12": 6.187572479248047, "*120": 6.186878204345703, "*121": 6.186878204345703, "*122": 6.186853408813477, "*123": 6.186824321746826, "*124": 6.186792850494385, "*125": 6.186830520629883, "*126": 6.186825275421143, "*127": 6.186825275421143, "*128": 6.186825275421143, "*129": 6.186817169189453, "*13": 6.187572479248047, "*130": 6.186817169189453, "*131": 6.186817169189453, "*132": 6.186817169189453, "*133": 6.186787128448486, "*134": 6.186787128448486, "*135": 6.186721324920654, "*136": 6.186721324920654, "*137": 6.186716079711914, "*138": 6.1866960525512695, "*139": 6.1866960525512695, "*14": 6.187572479248047, "*140": 6.1866774559021, "*141": 6.1866774559021, "*142": 6.1866774559021, "*143": 6.1866774559021, "*144": 6.1866774559021, "*145": 6.1866774559021, "*146": 6.186654567718506, "*147": 6.186654567718506, "*148": 6.186617374420166, "*149": 6.186596870422363, "*15": 6.187572479248047, "*150": 6.186596870422363, "*151": 6.186596870422363, "*152": 6.186582088470459, "*153": 6.186582088470459, "*154": 6.186563968658447, "*155": 6.186563968658447, "*156": 6.1865386962890625, "*157": 6.1865386962890625, "*158": 6.1865386962890625, "*159": 6.1865386962890625, "*16": 6.187555313110352, "*160": 6.1865386962890625, "*161": 6.186465263366699, "*162": 6.186465263366699, "*163": 6.186446666717529, "*164": 6.186429023742676, "*165": 6.186424255371094, "*166": 6.1864094734191895, "*167": 6.1864094734191895, "*168": 6.186404705047607, "*169": 6.186404705047607, "*17": 6.187524795532227, "*170": 6.186404705047607, "*171": 6.186404705047607, "*172": 6.18638801574707, "*173": 6.186382293701172, "*174": 6.186376094818115, "*175": 6.18633508682251, "*176": 6.186307907104492, "*177": 6.186295509338379, "*178": 6.186268329620361, "*179": 6.186268329620361, "*18": 6.1875128746032715, "*180": 6.186244010925293, "*181": 6.186244010925293, "*182": 6.186244010925293, "*183": 6.186238765716553, "*184": 6.186242580413818, "*185": 6.186242580413818, "*186": 6.186227321624756, "*187": 6.186227321624756, "*188": 6.186227321624756, "*189": 6.186227321624756, "*19": 6.187506675720215, "*190": 6.186227321624756, "*191": 6.186227321624756, "*192": 6.186222553253174, "*193": 6.186222553253174, "*194": 6.186222553253174, "*195": 6.186222076416016, "*196": 6.186222076416016, "*197": 6.186222076416016, "*198": 6.186221599578857, "*199": 6.186195373535156, "*2": 6.187680244445801, "*20": 6.187506675720215, "*200": 6.186195373535156, "*201": 6.186195373535156, "*202": 6.1861701011657715, "*203": 6.186167240142822, "*204": 6.186154365539551, "*205": 6.186154365539551, "*206": 6.186154365539551, "*207": 6.186154365539551, "*208": 6.186154365539551, "*209": 6.186154365539551, "*21": 6.187506675720215, "*210": 6.186148166656494, "*211": 6.186148166656494, "*212": 6.186148166656494, "*213": 6.186148166656494, "*214": 6.186148166656494, "*215": 6.186148166656494, "*216": 6.186131477355957, "*217": 6.186131477355957, "*218": 6.186120986938477, "*219": 6.186120986938477, "*22": 6.187507152557373, "*220": 6.186114311218262, "*221": 6.186114311218262, "*222": 6.186114311218262, "*223": 6.186114311218262, "*224": 6.186093807220459, "*225": 6.186089515686035, "*226": 6.186089515686035, "*227": 6.186079978942871, "*228": 6.186054706573486, "*229": 6.186054706573486, "*23": 6.187507152557373, "*230": 6.186054706573486, "*231": 6.186056137084961, "*232": 6.186056137084961, "*233": 6.186056613922119, "*234": 6.186056613922119, "*235": 6.1859893798828125, "*236": 6.1859893798828125, "*237": 6.185965061187744, "*238": 6.185965061187744, "*239": 6.185941696166992, "*24": 6.187501430511475, "*240": 6.1859283447265625, "*241": 6.1859283447265625, "*242": 6.1859283447265625, "*243": 6.1858954429626465, "*244": 6.185847282409668, "*245": 6.185847282409668, "*246": 6.185847282409668, "*247": 6.185847282409668, "*248": 6.185847282409668, "*249": 6.185847282409668, "*25": 6.187495231628418, "*250": 6.185821056365967, "*251": 6.185821056365967, "*252": 6.185774803161621, "*253": 6.185781478881836, "*254": 6.185781478881836, "*255": 6.1857476234436035, "*256": 6.18572998046875, "*257": 6.185718059539795, "*258": 6.185718059539795, "*259": 6.185718059539795, "*26": 6.187495231628418, "*260": 6.185718059539795, "*261": 6.185710430145264, "*262": 6.185710430145264, "*263": 6.185710430145264, "*264": 6.185710430145264, "*265": 6.185710430145264, "*266": 6.1856794357299805, "*267": 6.185654163360596, "*268": 6.185657978057861, "*269": 6.185657978057861, "*27": 6.187495231628418, "*270": 6.185628414154053, "*271": 6.185598850250244, "*272": 6.185598850250244, "*273": 6.18557071685791, "*274": 6.18557071685791, "*275": 6.18557071685791, "*276": 6.18557071685791, "*277": 6.185550212860107, "*278": 6.185550212860107, "*279": 6.185549259185791, "*28": 6.187490463256836, "*280": 6.185549259185791, "*281": 6.185549259185791, "*282": 6.185549259185791, "*283": 6.185484409332275, "*284": 6.185484409332275, "*285": 6.185484409332275, "*286": 6.185484409332275, "*287": 6.185484409332275, "*288": 6.185482978820801, "*289": 6.1854777336120605, "*29": 6.187490463256836, "*290": 6.1854777336120605, "*291": 6.1854634284973145, "*292": 6.185457706451416, "*293": 6.185457706451416, "*294": 6.185457706451416, "*295": 6.185457706451416, "*296": 6.1854448318481445, "*297": 6.185440540313721, "*298": 6.185404300689697, "*299": 6.185377597808838, "*3": 6.187680244445801, "*30": 6.187490463256836, "*300": 6.185345649719238, "*301": 6.185345649719238, "*302": 6.18532657623291, "*303": 6.185262680053711, "*304": 6.185262680053711, "*305": 6.18524169921875, "*306": 6.18524169921875, "*307": 6.18524169921875, "*308": 6.18524169921875, "*309": 6.18524169921875, "*31": 6.187490463256836, "*310": 6.185212135314941, "*311": 6.185202598571777, "*312": 6.185201168060303, "*313": 6.185201168060303, "*314": 6.185200214385986, "*315": 6.185200214385986, "*316": 6.185200214385986, "*317": 6.185200214385986, "*318": 6.185186862945557, "*319": 6.185186862945557, "*32": 6.187478065490723, "*320": 6.1851606369018555, "*321": 6.1851606369018555, "*322": 6.1851606369018555, "*323": 6.1851701736450195, "*324": 6.1851701736450195, "*325": 6.185135364532471, "*326": 6.185092449188232, "*327": 6.185092449188232, "*328": 6.185092449188232, "*329": 6.185092449188232, "*33": 6.187478065490723, "*330": 6.185092449188232, "*331": 6.185092449188232, "*332": 6.185092449188232, "*333": 6.185058116912842, "*334": 6.185058116912842, "*335": 6.185028553009033, "*336": 6.185028553009033, "*337": 6.185028553009033, "*338": 6.184974193572998, "*339": 6.184974193572998, "*34": 6.187478065490723, "*340": 6.1849365234375, "*341": 6.184927940368652, "*342": 6.184927940368652, "*343": 6.184874534606934, "*344": 6.184867858886719, "*345": 6.184861183166504, "*346": 6.184861183166504, "*347": 6.184861183166504, "*348": 6.184861183166504, "*349": 6.184861183166504, "*35": 6.187469005584717, "*350": 6.184861183166504, "*351": 6.184833526611328, "*352": 6.184816837310791, "*353": 6.184816837310791, "*354": 6.184807777404785, "*355": 6.184807777404785, "*356": 6.184795379638672, "*357": 6.184795379638672, "*358": 6.184795379638672, "*359": 6.184795379638672, "*36": 6.187469005584717, "*360": 6.184795379638672, "*361": 6.184802055358887, "*362": 6.184802055358887, "*363": 6.184802055358887, "*364": 6.184776782989502, "*365": 6.184757232666016, "*366": 6.1847357749938965, "*367": 6.1847357749938965, "*368": 6.1847357749938965, "*369": 6.184734344482422, "*37": 6.187469005584717, "*370": 6.184726715087891, "*371": 6.184726715087891, "*372": 6.184726715087891, "*373": 6.184726715087891, "*374": 6.184726715087891, "*375": 6.184720993041992, "*376": 6.184720993041992, "*377": 6.184720993041992, "*378": 6.184720993041992, "*379": 6.184720993041992, "*38": 6.1874613761901855, "*380": 6.184720993041992, "*381": 6.184720993041992, "*382": 6.184720993041992, "*383": 6.184698581695557, "*384": 6.184698581695557, "*385": 6.18469762802124, "*386": 6.184688568115234, "*387": 6.184693813323975, "*388": 6.184693813323975, "*389": 6.18468713760376, "*39": 6.187460899353027, "*390": 6.18468713760376, "*391": 6.1846604347229, "*392": 6.1846604347229, "*393": 6.1846604347229, "*394": 6.184654235839844, "*395": 6.184628486633301, "*396": 6.184628486633301, "*397": 6.184628486633301, "*398": 6.184628486633301, "*399": 6.184628486633301, "*4": 6.187680244445801, "*40": 6.187460899353027, "*400": 6.184612274169922, "*401": 6.184593200683594, "*402": 6.184593200683594, "*403": 6.184566974639893, "*404": 6.184552192687988, "*405": 6.184545993804932, "*406": 6.184545993804932, "*407": 6.184488296508789, "*408": 6.184488296508789, "*409": 6.184488296508789, "*41": 6.187460899353027, "*410": 6.184458255767822, "*411": 6.184458255767822, "*412": 6.184458255767822, "*413": 6.184458255767822, "*414": 6.184430122375488, "*415": 6.184410095214844, "*416": 6.184410095214844, "*417": 6.184410095214844, "*418": 6.184395790100098, "*419": 6.1843414306640625, "*42": 6.187466621398926, "*420": 6.184304714202881, "*421": 6.184304714202881, "*422": 6.184304714202881, "*423": 6.184304714202881, "*424": 6.184275150299072, "*425": 6.184275150299072, "*426": 6.184267520904541, "*427": 6.184267520904541, "*428": 6.184267520904541, "*429": 6.184267520904541, "*43": 6.187466621398926, "*430": 6.184267520904541, "*431": 6.184262275695801, "*432": 6.184262275695801, "*433": 6.184262275695801, "*434": 6.184262275695801, "*435": 6.184220314025879, "*436": 6.184220314025879, "*437": 6.184220314025879, "*438": 6.184220314025879, "*439": 6.184220314025879, "*44": 6.187450408935547, "*440": 6.184220314025879, "*441": 6.1842145919799805, "*442": 6.1842145919799805, "*443": 6.1842145919799805, "*444": 6.1842145919799805, "*445": 6.184215068817139, "*446": 6.184215068817139, "*447": 6.184215068817139, "*448": 6.184215068817139, "*449": 6.1842193603515625, "*45": 6.187427997589111, "*450": 6.18421745300293, "*451": 6.18421745300293, "*452": 6.184205532073975, "*453": 6.18418550491333, "*454": 6.18418550491333, "*455": 6.18418550491333, "*456": 6.1841840744018555, "*457": 6.1841840744018555, "*458": 6.18415641784668, "*459": 6.184157848358154, "*46": 6.187422752380371, "*460": 6.184157848358154, "*461": 6.184141635894775, "*462": 6.184139251708984, "*463": 6.184139251708984, "*464": 6.184139251708984, "*465": 6.184132099151611, "*466": 6.184132099151611, "*467": 6.184132099151611, "*468": 6.184132099151611, "*469": 6.184132099151611, "*47": 6.187422752380371, "*470": 6.184132099151611, "*471": 6.184110641479492, "*472": 6.184110641479492, "*473": 6.1841044425964355, "*474": 6.184077739715576, "*475": 6.184042930603027, "*476": 6.184006690979004, "*477": 6.184006690979004, "*478": 6.183976650238037, "*479": 6.183962345123291, "*48": 6.187422752380371, "*480": 6.183962345123291, "*481": 6.183962345123291, "*482": 6.183955669403076, "*483": 6.183955669403076, "*484": 6.183955669403076, "*485": 6.183930397033691, "*486": 6.183930397033691, "*487": 6.183930397033691, "*488": 6.183930397033691, "*489": 6.183930397033691, "*49": 6.187422752380371, "*490": 6.1839213371276855, "*491": 6.1839213371276855, "*492": 6.1839213371276855, "*493": 6.1839213371276855, "*494": 6.1839213371276855, "*495": 6.1839213371276855, "*496": 6.1839213371276855, "*497": 6.1839165687561035, "*498": 6.1839165687561035, "*499": 6.1839165687561035, "*5": 6.18768310546875, "*50": 6.18739652633667, "*500": 6.183862209320068, "*501": 6.183862209320068, "*502": 6.183862209320068, "*503": 6.183810234069824, "*504": 6.183810234069824, "*505": 6.1837921142578125, "*506": 6.1837921142578125, "*507": 6.1837921142578125, "*51": 6.18739652633667, "*52": 6.18739652633667, "*53": 6.18739652633667, "*54": 6.18739652633667, "*55": 6.18739652633667, "*56": 6.187367916107178, "*57": 6.1873626708984375, "*58": 6.187351703643799, "*59": 6.187351703643799, "*6": 6.18766975402832, "*60": 6.187342166900635, "*61": 6.187337398529053, "*62": 6.187328338623047, "*63": 6.187313556671143, "*64": 6.187313556671143, "*65": 6.187313556671143, "*66": 6.187314033508301, "*67": 6.187314033508301, "*68": 6.187314033508301, "*69": 6.187308311462402, "*7": 6.187633991241455, "*70": 6.187301158905029, "*71": 6.187301158905029, "*72": 6.187295436859131, "*73": 6.187263488769531, "*74": 6.1872639656066895, "*75": 6.1872639656066895, "*76": 6.187257766723633, "*77": 6.187231540679932, "*78": 6.187231540679932, "*79": 6.187231540679932, "*8": 6.187633991241455, "*80": 6.187226295471191, "*81": 6.187226295471191, "*82": 6.187189102172852, "*83": 6.187182903289795, "*84": 6.187182903289795, "*85": 6.187182903289795, "*86": 6.187180519104004, "*87": 6.187180519104004, "*88": 6.187180519104004, "*89": 6.187171936035156, "*9": 6.187617778778076, "*90": 6.187171936035156, "*91": 6.187171936035156, "*92": 6.187149524688721, "*93": 6.187122821807861, "*94": 6.187122821807861, "*95": 6.187107086181641, "*96": 6.1870927810668945, "*97": 6.187066555023193, "*98": 6.187066555023193, "*99": 6.187066555023193} 
            # test_standard_cube_briggsbwtaper.exp_bmaj_dict
            exp_bmaj_dict = {"*0": 8.891188621520996, "*1": 8.891138076782227, "*10": 8.891054153442383, "*100": 8.889816284179688, "*101": 8.889824867248535, "*102": 8.889791488647461, "*103": 8.88974380493164, "*104": 8.889734268188477, "*105": 8.889734268188477, "*106": 8.889701843261719, "*107": 8.889701843261719, "*108": 8.889681816101074, "*109": 8.889681816101074, "*11": 8.891039848327637, "*110": 8.889681816101074, "*111": 8.889681816101074, "*112": 8.889681816101074, "*113": 8.889681816101074, "*114": 8.889681816101074, "*115": 8.889681816101074, "*116": 8.889681816101074, "*117": 8.889642715454102, "*118": 8.889630317687988, "*119": 8.889630317687988, "*12": 8.89102554321289, "*120": 8.889619827270508, "*121": 8.889619827270508, "*122": 8.889603614807129, "*123": 8.889480590820312, "*124": 8.889419555664062, "*125": 8.889314651489258, "*126": 8.88927936553955, "*127": 8.88927936553955, "*128": 8.88927936553955, "*129": 8.88926887512207, "*13": 8.89102554321289, "*130": 8.88926887512207, "*131": 8.88926887512207, "*132": 8.88926887512207, "*133": 8.88925838470459, "*134": 8.88925838470459, "*135": 8.889309883117676, "*136": 8.889309883117676, "*137": 8.889289855957031, "*138": 8.889214515686035, "*139": 8.889214515686035, "*14": 8.89102554321289, "*140": 8.88923454284668, "*141": 8.88923454284668, "*142": 8.88923454284668, "*143": 8.88923454284668, "*144": 8.88923454284668, "*145": 8.88923454284668, "*146": 8.889192581176758, "*147": 8.889192581176758, "*148": 8.889154434204102, "*149": 8.889152526855469, "*15": 8.89102554321289, "*150": 8.889152526855469, "*151": 8.889152526855469, "*152": 8.889168739318848, "*153": 8.889168739318848, "*154": 8.889181137084961, "*155": 8.889181137084961, "*156": 8.889180183410645, "*157": 8.889180183410645, "*158": 8.889180183410645, "*159": 8.889180183410645, "*16": 8.891008377075195, "*160": 8.889180183410645, "*161": 8.889176368713379, "*162": 8.889176368713379, "*163": 8.889165878295898, "*164": 8.88913345336914, "*165": 8.889095306396484, "*166": 8.889080047607422, "*167": 8.889080047607422, "*168": 8.889060020446777, "*169": 8.889060020446777, "*17": 8.89089584350586, "*170": 8.889060020446777, "*171": 8.889060020446777, "*172": 8.889078140258789, "*173": 8.889034271240234, "*174": 8.88901138305664, "*175": 8.889034271240234, "*176": 8.889015197753906, "*177": 8.889002799987793, "*178": 8.888991355895996, "*179": 8.888991355895996, "*18": 8.89090347290039, "*180": 8.888973236083984, "*181": 8.888973236083984, "*182": 8.888973236083984, "*183": 8.88895034790039, "*184": 8.888908386230469, "*185": 8.888908386230469, "*186": 8.88890552520752, "*187": 8.88890552520752, "*188": 8.88890552520752, "*189": 8.88890552520752, "*19": 8.890859603881836, "*190": 8.88890552520752, "*191": 8.88890552520752, "*192": 8.888861656188965, "*193": 8.888861656188965, "*194": 8.888861656188965, "*195": 8.888860702514648, "*196": 8.888860702514648, "*197": 8.888860702514648, "*198": 8.888829231262207, "*199": 8.888810157775879, "*2": 8.891165733337402, "*20": 8.890859603881836, "*200": 8.888810157775879, "*201": 8.888810157775879, "*202": 8.88879108428955, "*203": 8.888731956481934, "*204": 8.888740539550781, "*205": 8.888740539550781, "*206": 8.888740539550781, "*207": 8.888740539550781, "*208": 8.888740539550781, "*209": 8.888740539550781, "*21": 8.890859603881836, "*210": 8.888670921325684, "*211": 8.888670921325684, "*212": 8.888670921325684, "*213": 8.888670921325684, "*214": 8.888670921325684, "*215": 8.888670921325684, "*216": 8.888643264770508, "*217": 8.888643264770508, "*218": 8.888569831848145, "*219": 8.888569831848145, "*22": 8.890859603881836, "*220": 8.888533592224121, "*221": 8.888533592224121, "*222": 8.888533592224121, "*223": 8.888533592224121, "*224": 8.888554573059082, "*225": 8.888519287109375, "*226": 8.888519287109375, "*227": 8.888473510742188, "*228": 8.88846206665039, "*229": 8.88846206665039, "*23": 8.890859603881836, "*230": 8.88846206665039, "*231": 8.888426780700684, "*232": 8.888426780700684, "*233": 8.888361930847168, "*234": 8.888361930847168, "*235": 8.888412475585938, "*236": 8.888412475585938, "*237": 8.888394355773926, "*238": 8.888394355773926, "*239": 8.888354301452637, "*24": 8.890816688537598, "*240": 8.888344764709473, "*241": 8.888344764709473, "*242": 8.888344764709473, "*243": 8.888361930847168, "*244": 8.888338088989258, "*245": 8.888338088989258, "*246": 8.888338088989258, "*247": 8.888338088989258, "*248": 8.888338088989258, "*249": 8.888338088989258, "*25": 8.890801429748535, "*250": 8.88831901550293, "*251": 8.88831901550293, "*252": 8.888285636901855, "*253": 8.888249397277832, "*254": 8.888249397277832, "*255": 8.888227462768555, "*256": 8.888239860534668, "*257": 8.888182640075684, "*258": 8.888182640075684, "*259": 8.888182640075684, "*26": 8.890801429748535, "*260": 8.888182640075684, "*261": 8.888158798217773, "*262": 8.888158798217773, "*263": 8.888158798217773, "*264": 8.888158798217773, "*265": 8.888158798217773, "*266": 8.888177871704102, "*267": 8.888092994689941, "*268": 8.888045310974121, "*269": 8.888045310974121, "*27": 8.890801429748535, "*270": 8.88806438446045, "*271": 8.888053894042969, "*272": 8.888053894042969, "*273": 8.888071060180664, "*274": 8.888071060180664, "*275": 8.888071060180664, "*276": 8.88802433013916, "*277": 8.88798999786377, "*278": 8.88798999786377, "*279": 8.887958526611328, "*28": 8.890766143798828, "*280": 8.887958526611328, "*281": 8.887958526611328, "*282": 8.887958526611328, "*283": 8.887968063354492, "*284": 8.887968063354492, "*285": 8.887968063354492, "*286": 8.887968063354492, "*287": 8.887968063354492, "*288": 8.8878755569458, "*289": 8.887828826904297, "*29": 8.890766143798828, "*290": 8.887828826904297, "*291": 8.887823104858398, "*292": 8.88769245147705, "*293": 8.88769245147705, "*294": 8.88769245147705, "*295": 8.88769245147705, "*296": 8.887701034545898, "*297": 8.88759708404541, "*298": 8.887588500976562, "*299": 8.887596130371094, "*3": 8.891165733337402, "*30": 8.890766143798828, "*300": 8.88758659362793, "*301": 8.88758659362793, "*302": 8.887445449829102, "*303": 8.887458801269531, "*304": 8.887458801269531, "*305": 8.88736629486084, "*306": 8.88736629486084, "*307": 8.88736629486084, "*308": 8.88736629486084, "*309": 8.88736629486084, "*31": 8.890766143798828, "*310": 8.887368202209473, "*311": 8.887341499328613, "*312": 8.887307167053223, "*313": 8.887307167053223, "*314": 8.88726806640625, "*315": 8.88726806640625, "*316": 8.88726806640625, "*317": 8.88726806640625, "*318": 8.887258529663086, "*319": 8.887258529663086, "*32": 8.890764236450195, "*320": 8.887239456176758, "*321": 8.887239456176758, "*322": 8.887239456176758, "*323": 8.887185096740723, "*324": 8.887185096740723, "*325": 8.887202262878418, "*326": 8.887165069580078, "*327": 8.887165069580078, "*328": 8.887165069580078, "*329": 8.887165069580078, "*33": 8.890764236450195, "*330": 8.887165069580078, "*331": 8.887165069580078, "*332": 8.887165069580078, "*333": 8.88715648651123, "*334": 8.88715648651123, "*335": 8.88715648651123, "*336": 8.88715648651123, "*337": 8.88715648651123, "*338": 8.887151718139648, "*339": 8.887151718139648, "*34": 8.890764236450195, "*340": 8.88717269897461, "*341": 8.887157440185547, "*342": 8.887157440185547, "*343": 8.887078285217285, "*344": 8.887033462524414, "*345": 8.886980056762695, "*346": 8.886980056762695, "*347": 8.886980056762695, "*348": 8.886980056762695, "*349": 8.886980056762695, "*35": 8.890774726867676, "*350": 8.886980056762695, "*351": 8.886980056762695, "*352": 8.886878967285156, "*353": 8.886878967285156, "*354": 8.886759757995605, "*355": 8.886759757995605, "*356": 8.886747360229492, "*357": 8.886747360229492, "*358": 8.886747360229492, "*359": 8.886747360229492, "*36": 8.890774726867676, "*360": 8.886747360229492, "*361": 8.886665344238281, "*362": 8.886665344238281, "*363": 8.886665344238281, "*364": 8.886652946472168, "*365": 8.886676788330078, "*366": 8.886693000793457, "*367": 8.886693000793457, "*368": 8.886693000793457, "*369": 8.8866548538208, "*37": 8.890774726867676, "*370": 8.886632919311523, "*371": 8.886632919311523, "*372": 8.886632919311523, "*373": 8.886632919311523, "*374": 8.886632919311523, "*375": 8.886564254760742, "*376": 8.886564254760742, "*377": 8.886564254760742, "*378": 8.886564254760742, "*379": 8.886564254760742, "*38": 8.890726089477539, "*380": 8.886564254760742, "*381": 8.886564254760742, "*382": 8.886564254760742, "*383": 8.886563301086426, "*384": 8.886563301086426, "*385": 8.886528015136719, "*386": 8.886497497558594, "*387": 8.886425018310547, "*388": 8.886425018310547, "*389": 8.886369705200195, "*39": 8.890726089477539, "*390": 8.886369705200195, "*391": 8.886350631713867, "*392": 8.886350631713867, "*393": 8.886350631713867, "*394": 8.886300086975098, "*395": 8.886270523071289, "*396": 8.886270523071289, "*397": 8.886270523071289, "*398": 8.886270523071289, "*399": 8.886270523071289, "*4": 8.891165733337402, "*40": 8.890726089477539, "*400": 8.886260986328125, "*401": 8.886275291442871, "*402": 8.886275291442871, "*403": 8.886244773864746, "*404": 8.886176109313965, "*405": 8.886120796203613, "*406": 8.886120796203613, "*407": 8.886150360107422, "*408": 8.886150360107422, "*409": 8.886150360107422, "*41": 8.890726089477539, "*410": 8.886139869689941, "*411": 8.886139869689941, "*412": 8.886139869689941, "*413": 8.886139869689941, "*414": 8.88612174987793, "*415": 8.886110305786133, "*416": 8.886110305786133, "*417": 8.886110305786133, "*418": 8.886038780212402, "*419": 8.885890007019043, "*42": 8.890661239624023, "*420": 8.885895729064941, "*421": 8.885895729064941, "*422": 8.885895729064941, "*423": 8.885895729064941, "*424": 8.885836601257324, "*425": 8.885836601257324, "*426": 8.885828971862793, "*427": 8.885828971862793, "*428": 8.885828971862793, "*429": 8.885828971862793, "*43": 8.890661239624023, "*430": 8.885828971862793, "*431": 8.885787010192871, "*432": 8.885787010192871, "*433": 8.885787010192871, "*434": 8.885787010192871, "*435": 8.885794639587402, "*436": 8.885794639587402, "*437": 8.885794639587402, "*438": 8.885794639587402, "*439": 8.885794639587402, "*44": 8.890636444091797, "*440": 8.885794639587402, "*441": 8.885746002197266, "*442": 8.885746002197266, "*443": 8.885746002197266, "*444": 8.885746002197266, "*445": 8.885746002197266, "*446": 8.885746002197266, "*447": 8.885746002197266, "*448": 8.885746002197266, "*449": 8.885686874389648, "*45": 8.890645027160645, "*450": 8.885637283325195, "*451": 8.885637283325195, "*452": 8.885584831237793, "*453": 8.88558292388916, "*454": 8.88558292388916, "*455": 8.88558292388916, "*456": 8.885557174682617, "*457": 8.885557174682617, "*458": 8.885503768920898, "*459": 8.88544750213623, "*46": 8.89060115814209, "*460": 8.88544750213623, "*461": 8.885454177856445, "*462": 8.885412216186523, "*463": 8.885412216186523, "*464": 8.885412216186523, "*465": 8.88538932800293, "*466": 8.88538932800293, "*467": 8.88538932800293, "*468": 8.88538932800293, "*469": 8.88538932800293, "*47": 8.89060115814209, "*470": 8.88538932800293, "*471": 8.885375022888184, "*472": 8.885375022888184, "*473": 8.88534164428711, "*474": 8.885329246520996, "*475": 8.885315895080566, "*476": 8.885306358337402, "*477": 8.885306358337402, "*478": 8.885294914245605, "*479": 8.885313034057617, "*48": 8.89060115814209, "*480": 8.885313034057617, "*481": 8.885313034057617, "*482": 8.885234832763672, "*483": 8.885234832763672, "*484": 8.885234832763672, "*485": 8.885233879089355, "*486": 8.885233879089355, "*487": 8.885233879089355, "*488": 8.885233879089355, "*489": 8.885233879089355, "*49": 8.89060115814209, "*490": 8.885217666625977, "*491": 8.885217666625977, "*492": 8.885217666625977, "*493": 8.885217666625977, "*494": 8.885217666625977, "*495": 8.885217666625977, "*496": 8.885217666625977, "*497": 8.885114669799805, "*498": 8.885114669799805, "*499": 8.885114669799805, "*5": 8.891115188598633, "*50": 8.890581130981445, "*500": 8.885082244873047, "*501": 8.885082244873047, "*502": 8.885082244873047, "*503": 8.885062217712402, "*504": 8.885062217712402, "*505": 8.88504409790039, "*506": 8.88504409790039, "*507": 8.88504409790039, "*51": 8.890581130981445, "*52": 8.890534400939941, "*53": 8.890534400939941, "*54": 8.890534400939941, "*55": 8.890534400939941, "*56": 8.890515327453613, "*57": 8.890485763549805, "*58": 8.890435218811035, "*59": 8.890436172485352, "*6": 8.89104175567627, "*60": 8.890392303466797, "*61": 8.890374183654785, "*62": 8.890360832214355, "*63": 8.890379905700684, "*64": 8.890379905700684, "*65": 8.890379905700684, "*66": 8.890379905700684, "*67": 8.890379905700684, "*68": 8.890379905700684, "*69": 8.890341758728027, "*7": 8.891036987304688, "*70": 8.89029598236084, "*71": 8.89029598236084, "*72": 8.890267372131348, "*73": 8.890254974365234, "*74": 8.8902006149292, "*75": 8.8902006149292, "*76": 8.890145301818848, "*77": 8.890096664428711, "*78": 8.890096664428711, "*79": 8.890096664428711, "*8": 8.891036987304688, "*80": 8.89009952545166, "*81": 8.89009952545166, "*82": 8.890077590942383, "*83": 8.890035629272461, "*84": 8.890035629272461, "*85": 8.890035629272461, "*86": 8.889941215515137, "*87": 8.889941215515137, "*88": 8.889941215515137, "*89": 8.889917373657227, "*9": 8.891054153442383, "*90": 8.889917373657227, "*91": 8.889917373657227, "*92": 8.88991641998291, "*93": 8.889906883239746, "*94": 8.889906883239746, "*95": 8.889881134033203, "*96": 8.88982105255127, "*97": 8.889833450317383, "*98": 8.889833450317383, "*99": 8.889833450317383}
            # test_standard_cube_briggsbwtaper.exp_pa_dict
            exp_pa_dict = {"*0": 67.88911437988281, "*1": 67.88862609863281, "*10": 67.88904571533203, "*100": 67.8851089477539, "*101": 67.88517761230469, "*102": 67.88480377197266, "*103": 67.8844985961914, "*104": 67.88472747802734, "*105": 67.88472747802734, "*106": 67.88450622558594, "*107": 67.88450622558594, "*108": 67.88478088378906, "*109": 67.88478088378906, "*11": 67.8893814086914, "*110": 67.88478088378906, "*111": 67.88478088378906, "*112": 67.88478088378906, "*113": 67.88478088378906, "*114": 67.88478088378906, "*115": 67.88478088378906, "*116": 67.88478088378906, "*117": 67.8846206665039, "*118": 67.88495635986328, "*119": 67.88495635986328, "*12": 67.88976287841797, "*120": 67.88533020019531, "*121": 67.88533020019531, "*122": 67.88558959960938, "*123": 67.88492584228516, "*124": 67.8849105834961, "*125": 67.88545989990234, "*126": 67.88529968261719, "*127": 67.88529968261719, "*128": 67.88529968261719, "*129": 67.88470458984375, "*13": 67.88976287841797, "*130": 67.88470458984375, "*131": 67.88470458984375, "*132": 67.88470458984375, "*133": 67.88504028320312, "*134": 67.88504028320312, "*135": 67.88516998291016, "*136": 67.88516998291016, "*137": 67.8848876953125, "*138": 67.8851318359375, "*139": 67.8851318359375, "*14": 67.88976287841797, "*140": 67.88496398925781, "*141": 67.88496398925781, "*142": 67.88496398925781, "*143": 67.88496398925781, "*144": 67.88496398925781, "*145": 67.88496398925781, "*146": 67.8846664428711, "*147": 67.8846664428711, "*148": 67.88465881347656, "*149": 67.88484954833984, "*15": 67.88976287841797, "*150": 67.88484954833984, "*151": 67.88484954833984, "*152": 67.88421630859375, "*153": 67.88421630859375, "*154": 67.88420104980469, "*155": 67.88420104980469, "*156": 67.88442993164062, "*157": 67.88442993164062, "*158": 67.88442993164062, "*159": 67.88442993164062, "*16": 67.89012145996094, "*160": 67.88442993164062, "*161": 67.88473510742188, "*162": 67.88473510742188, "*163": 67.88501739501953, "*164": 67.88580322265625, "*165": 67.88563537597656, "*166": 67.88595581054688, "*167": 67.88595581054688, "*168": 67.88567352294922, "*169": 67.88567352294922, "*17": 67.89022827148438, "*170": 67.88567352294922, "*171": 67.88567352294922, "*172": 67.88494873046875, "*173": 67.88475799560547, "*174": 67.88445281982422, "*175": 67.8838119506836, "*176": 67.88408660888672, "*177": 67.88433837890625, "*178": 67.88463592529297, "*179": 67.88463592529297, "*18": 67.89019775390625, "*180": 67.88488006591797, "*181": 67.88488006591797, "*182": 67.88488006591797, "*183": 67.88459777832031, "*184": 67.88463592529297, "*185": 67.88463592529297, "*186": 67.8847885131836, "*187": 67.8847885131836, "*188": 67.8847885131836, "*189": 67.8847885131836, "*19": 67.88999938964844, "*190": 67.8847885131836, "*191": 67.8847885131836, "*192": 67.88453674316406, "*193": 67.88453674316406, "*194": 67.88453674316406, "*195": 67.88453674316406, "*196": 67.88453674316406, "*197": 67.88453674316406, "*198": 67.88445281982422, "*199": 67.88471984863281, "*2": 67.88844299316406, "*20": 67.88999938964844, "*200": 67.88471984863281, "*201": 67.88471984863281, "*202": 67.88497161865234, "*203": 67.88440704345703, "*204": 67.88438415527344, "*205": 67.88438415527344, "*206": 67.88438415527344, "*207": 67.88438415527344, "*208": 67.88438415527344, "*209": 67.88438415527344, "*21": 67.88999938964844, "*210": 67.88417053222656, "*211": 67.88417053222656, "*212": 67.88417053222656, "*213": 67.88417053222656, "*214": 67.88417053222656, "*215": 67.88417053222656, "*216": 67.88472747802734, "*217": 67.88472747802734, "*218": 67.88542175292969, "*219": 67.88542175292969, "*22": 67.88999938964844, "*220": 67.88502502441406, "*221": 67.88502502441406, "*222": 67.88502502441406, "*223": 67.88502502441406, "*224": 67.8849105834961, "*225": 67.88475036621094, "*226": 67.88475036621094, "*227": 67.88414764404297, "*228": 67.88443756103516, "*229": 67.88443756103516, "*23": 67.88999938964844, "*230": 67.88443756103516, "*231": 67.88507843017578, "*232": 67.88507843017578, "*233": 67.88499450683594, "*234": 67.88499450683594, "*235": 67.8851547241211, "*236": 67.8851547241211, "*237": 67.88540649414062, "*238": 67.88540649414062, "*239": 67.88505554199219, "*24": 67.88982391357422, "*240": 67.88528442382812, "*241": 67.88528442382812, "*242": 67.88528442382812, "*243": 67.88533020019531, "*244": 67.88581848144531, "*245": 67.88581848144531, "*246": 67.88581848144531, "*247": 67.88581848144531, "*248": 67.88581848144531, "*249": 67.88581085205078, "*25": 67.88949584960938, "*250": 67.88623809814453, "*251": 67.88623809814453, "*252": 67.88670349121094, "*253": 67.88667297363281, "*254": 67.88667297363281, "*255": 67.88690948486328, "*256": 67.88687896728516, "*257": 67.88607788085938, "*258": 67.88607788085938, "*259": 67.88607788085938, "*26": 67.88949584960938, "*260": 67.88607788085938, "*261": 67.88653564453125, "*262": 67.88653564453125, "*263": 67.88653564453125, "*264": 67.88653564453125, "*265": 67.88653564453125, "*266": 67.88668060302734, "*267": 67.88717651367188, "*268": 67.88722229003906, "*269": 67.88722229003906, "*27": 67.88949584960938, "*270": 67.88722229003906, "*271": 67.88756561279297, "*272": 67.88756561279297, "*273": 67.8876724243164, "*274": 67.8876724243164, "*275": 67.8876724243164, "*276": 67.88760375976562, "*277": 67.88844299316406, "*278": 67.88844299316406, "*279": 67.88835906982422, "*28": 67.88933563232422, "*280": 67.88835906982422, "*281": 67.88835906982422, "*282": 67.88835906982422, "*283": 67.88874053955078, "*284": 67.88874053955078, "*285": 67.88874053955078, "*286": 67.88874053955078, "*287": 67.88874053955078, "*288": 67.8885269165039, "*289": 67.88822174072266, "*29": 67.88933563232422, "*290": 67.88822174072266, "*291": 67.88841247558594, "*292": 67.88805389404297, "*293": 67.88805389404297, "*294": 67.88805389404297, "*295": 67.88805389404297, "*296": 67.88802337646484, "*297": 67.88714599609375, "*298": 67.88751983642578, "*299": 67.88763427734375, "*3": 67.88844299316406, "*30": 67.88933563232422, "*300": 67.88797760009766, "*301": 67.88797760009766, "*302": 67.8873291015625, "*303": 67.88766479492188, "*304": 67.88766479492188, "*305": 67.887451171875, "*306": 67.887451171875, "*307": 67.887451171875, "*308": 67.887451171875, "*309": 67.887451171875, "*31": 67.88933563232422, "*310": 67.88771057128906, "*311": 67.88711547851562, "*312": 67.88702392578125, "*313": 67.88702392578125, "*314": 67.88690185546875, "*315": 67.88690185546875, "*316": 67.88690185546875, "*317": 67.88690185546875, "*318": 67.88713073730469, "*319": 67.88713073730469, "*32": 67.88945770263672, "*320": 67.88755798339844, "*321": 67.88755798339844, "*322": 67.88755798339844, "*323": 67.88748168945312, "*324": 67.88748168945312, "*325": 67.8875503540039, "*326": 67.88845825195312, "*327": 67.88845825195312, "*328": 67.88845825195312, "*329": 67.88845825195312, "*33": 67.88945770263672, "*330": 67.88845825195312, "*331": 67.88845825195312, "*332": 67.88845825195312, "*333": 67.88880920410156, "*334": 67.88880920410156, "*335": 67.88907623291016, "*336": 67.88907623291016, "*337": 67.88907623291016, "*338": 67.88957214355469, "*339": 67.88957214355469, "*34": 67.88945770263672, "*340": 67.88961029052734, "*341": 67.888916015625, "*342": 67.888916015625, "*343": 67.88896942138672, "*344": 67.88858032226562, "*345": 67.88825988769531, "*346": 67.88825988769531, "*347": 67.88825988769531, "*348": 67.88825988769531, "*349": 67.88825988769531, "*35": 67.88904571533203, "*350": 67.88825988769531, "*351": 67.88851165771484, "*352": 67.88846588134766, "*353": 67.88846588134766, "*354": 67.88778686523438, "*355": 67.88778686523438, "*356": 67.8880386352539, "*357": 67.8880386352539, "*358": 67.8880386352539, "*359": 67.8880386352539, "*36": 67.88904571533203, "*360": 67.8880386352539, "*361": 67.88784790039062, "*362": 67.88784790039062, "*363": 67.88784790039062, "*364": 67.88813781738281, "*365": 67.88797760009766, "*366": 67.88792419433594, "*367": 67.88792419433594, "*368": 67.88792419433594, "*369": 67.8878173828125, "*37": 67.88904571533203, "*370": 67.88822937011719, "*371": 67.88822937011719, "*372": 67.88822937011719, "*373": 67.88822937011719, "*374": 67.88822937011719, "*375": 67.88800811767578, "*376": 67.88800811767578, "*377": 67.88800811767578, "*378": 67.88800811767578, "*379": 67.88800811767578, "*38": 67.88861846923828, "*380": 67.88800811767578, "*381": 67.88800811767578, "*382": 67.88800811767578, "*383": 67.88822174072266, "*384": 67.88822174072266, "*385": 67.88812255859375, "*386": 67.88761138916016, "*387": 67.88765716552734, "*388": 67.88765716552734, "*389": 67.88729858398438, "*39": 67.88861846923828, "*390": 67.88729858398438, "*391": 67.88772583007812, "*392": 67.88772583007812, "*393": 67.88772583007812, "*394": 67.8874282836914, "*395": 67.88814544677734, "*396": 67.88814544677734, "*397": 67.88814544677734, "*398": 67.88814544677734, "*399": 67.88814544677734, "*4": 67.88844299316406, "*40": 67.88861846923828, "*400": 67.88836669921875, "*401": 67.8883285522461, "*402": 67.8883285522461, "*403": 67.88897705078125, "*404": 67.88809967041016, "*405": 67.88770294189453, "*406": 67.88770294189453, "*407": 67.887939453125, "*408": 67.887939453125, "*409": 67.887939453125, "*41": 67.88861846923828, "*410": 67.88827514648438, "*411": 67.88827514648438, "*412": 67.88827514648438, "*413": 67.88827514648438, "*414": 67.88873291015625, "*415": 67.88903045654297, "*416": 67.88903045654297, "*417": 67.88903045654297, "*418": 67.8884048461914, "*419": 67.88810729980469, "*42": 67.888671875, "*420": 67.88841247558594, "*421": 67.88841247558594, "*422": 67.88841247558594, "*423": 67.88841247558594, "*424": 67.88876342773438, "*425": 67.88876342773438, "*426": 67.88822937011719, "*427": 67.88822937011719, "*428": 67.88822937011719, "*429": 67.88822937011719, "*43": 67.888671875, "*430": 67.88822937011719, "*431": 67.88794708251953, "*432": 67.88794708251953, "*433": 67.88794708251953, "*434": 67.88794708251953, "*435": 67.88660430908203, "*436": 67.88660430908203, "*437": 67.88660430908203, "*438": 67.88660430908203, "*439": 67.88660430908203, "*44": 67.88928985595703, "*440": 67.88660430908203, "*441": 67.88628387451172, "*442": 67.88628387451172, "*443": 67.88628387451172, "*444": 67.88628387451172, "*445": 67.88628387451172, "*446": 67.88628387451172, "*447": 67.88628387451172, "*448": 67.88628387451172, "*449": 67.88634490966797, "*45": 67.88932037353516, "*450": 67.8858413696289, "*451": 67.8858413696289, "*452": 67.88516998291016, "*453": 67.88536834716797, "*454": 67.88536834716797, "*455": 67.88536834716797, "*456": 67.88578796386719, "*457": 67.88578796386719, "*458": 67.8861312866211, "*459": 67.88607788085938, "*46": 67.88905334472656, "*460": 67.88607788085938, "*461": 67.8861083984375, "*462": 67.88568115234375, "*463": 67.88568115234375, "*464": 67.88568115234375, "*465": 67.8852310180664, "*466": 67.8852310180664, "*467": 67.8852310180664, "*468": 67.8852310180664, "*469": 67.8852310180664, "*47": 67.88905334472656, "*470": 67.8852310180664, "*471": 67.88560485839844, "*472": 67.88560485839844, "*473": 67.8852310180664, "*474": 67.88553619384766, "*475": 67.88593292236328, "*476": 67.88629913330078, "*477": 67.88629913330078, "*478": 67.88662719726562, "*479": 67.88603210449219, "*48": 67.88905334472656, "*480": 67.88603210449219, "*481": 67.88603210449219, "*482": 67.88577270507812, "*483": 67.88577270507812, "*484": 67.88578033447266, "*485": 67.88600158691406, "*486": 67.88600158691406, "*487": 67.88600158691406, "*488": 67.88600158691406, "*489": 67.8860092163086, "*49": 67.88905334472656, "*490": 67.88530731201172, "*491": 67.88530731201172, "*492": 67.88530731201172, "*493": 67.88530731201172, "*494": 67.88530731201172, "*495": 67.88530731201172, "*496": 67.88530731201172, "*497": 67.88445281982422, "*498": 67.88445281982422, "*499": 67.88445281982422, "*5": 67.88935852050781, "*50": 67.88932037353516, "*500": 67.88520050048828, "*501": 67.88520050048828, "*502": 67.88520050048828, "*503": 67.8858413696289, "*504": 67.8858413696289, "*505": 67.88623046875, "*506": 67.88623046875, "*507": 67.88623046875, "*51": 67.88932037353516, "*52": 67.8892593383789, "*53": 67.8892593383789, "*54": 67.8892593383789, "*55": 67.8892593383789, "*56": 67.88983154296875, "*57": 67.88948822021484, "*58": 67.88876342773438, "*59": 67.88876342773438, "*6": 67.88946533203125, "*60": 67.88812255859375, "*61": 67.88846588134766, "*62": 67.88780212402344, "*63": 67.8871841430664, "*64": 67.8871841430664, "*65": 67.8871841430664, "*66": 67.8871841430664, "*67": 67.8871841430664, "*68": 67.8871841430664, "*69": 67.88702392578125, "*7": 67.88978576660156, "*70": 67.88661193847656, "*71": 67.88661193847656, "*72": 67.88624572753906, "*73": 67.88660430908203, "*74": 67.88655090332031, "*75": 67.88655090332031, "*76": 67.88615417480469, "*77": 67.88612365722656, "*78": 67.88612365722656, "*79": 67.88612365722656, "*8": 67.88978576660156, "*80": 67.8858871459961, "*81": 67.8858871459961, "*82": 67.88557434082031, "*83": 67.88539123535156, "*84": 67.88539123535156, "*85": 67.88539123535156, "*86": 67.8851547241211, "*87": 67.8851547241211, "*88": 67.8851547241211, "*89": 67.88460540771484, "*9": 67.88904571533203, "*90": 67.88460540771484, "*91": 67.88460540771484, "*92": 67.88481140136719, "*93": 67.88510131835938, "*94": 67.88510131835938, "*95": 67.88561248779297, "*96": 67.88472747802734, "*97": 67.8847427368164, "*98": 67.8847427368164, "*99": 67.8847427368164} 
            exp_pa_dict = {'*0': 67.88911437988281,'*1': 67.88862609863281,'*10': 67.88904571533203,'*100': 67.8851089477539,'*101': 67.88517761230469,'*102': 67.88480377197266,'*103': 67.8844985961914,'*104': 67.88472747802734,'*105': 67.88472747802734,'*106': 67.88450622558594,'*107': 67.88450622558594,'*108': 67.88478088378906,'*109': 67.88478088378906,'*11': 67.8893814086914,'*110': 67.88478088378906,'*111': 67.88478088378906,'*112': 67.88478088378906,'*113': 67.88478088378906,'*114': 67.88478088378906,'*115': 67.88478088378906,'*116': 67.88478088378906,'*117': 67.8846206665039,'*118': 67.88495635986328,'*119': 67.88495635986328,'*12': 67.88976287841797,'*120': 67.88533020019531,'*121': 67.88533020019531,'*122': 67.88558959960938,'*123': 67.88492584228516,'*124': 67.8849105834961,'*125': 67.88545989990234,'*126': 67.88529968261719,'*127': 67.88529968261719,'*128': 67.88529968261719,'*129': 67.88470458984375,'*13': 67.88976287841797,'*130': 67.88470458984375,'*131': 67.88470458984375,'*132': 67.88470458984375,'*133': 67.88504028320312,'*134': 67.88504028320312,'*135': 67.88516998291016,'*136': 67.88516998291016,'*137': 67.8848876953125,'*138': 67.8851318359375,'*139': 67.8851318359375,'*14': 67.88976287841797,'*140': 67.88496398925781,'*141': 67.88496398925781,'*142': 67.88496398925781,'*143': 67.88496398925781,'*144': 67.88496398925781,'*145': 67.88496398925781,'*146': 67.8846664428711,'*147': 67.8846664428711,'*148': 67.88465881347656,'*149': 67.88484954833984,'*15': 67.88976287841797,'*150': 67.88484954833984,'*151': 67.88484954833984,'*152': 67.88421630859375,'*153': 67.88421630859375,'*154': 67.88420104980469,'*155': 67.88420104980469,'*156': 67.88442993164062,'*157': 67.88442993164062,'*158': 67.88442993164062,'*159': 67.88442993164062,'*16': 67.89012145996094,'*160': 67.88442993164062,'*161': 67.88473510742188,'*162': 67.88473510742188,'*163': 67.88501739501953,'*164': 67.88580322265625,'*165': 67.88563537597656,'*166': 67.88595581054688,'*167': 67.88595581054688,'*168': 67.88567352294922,'*169': 67.88567352294922,'*17': 67.89022827148438,'*170': 67.88567352294922,'*171': 67.88567352294922,'*172': 67.88494873046875,'*173': 67.88475799560547,'*174': 67.88445281982422,'*175': 67.8838119506836,'*176': 67.88408660888672,'*177': 67.88433837890625,'*178': 67.88463592529297,'*179': 67.88463592529297,'*18': 67.89019775390625,'*180': 67.88488006591797,'*181': 67.88488006591797,'*182': 67.88488006591797,'*183': 67.88459777832031,'*184': 67.88463592529297,'*185': 67.88463592529297,'*186': 67.8847885131836,'*187': 67.8847885131836,'*188': 67.8847885131836,'*189': 67.8847885131836,'*19': 67.88999938964844,'*190': 67.8847885131836,'*191': 67.8847885131836,'*192': 67.88453674316406,'*193': 67.88453674316406,'*194': 67.88453674316406,'*195': 67.88453674316406,'*196': 67.88453674316406,'*197': 67.88453674316406,'*198': 67.88445281982422,'*199': 67.88471984863281,'*2': 67.88844299316406,'*20': 67.88999938964844,'*200': 67.88471984863281,'*201': 67.88471984863281,'*202': 67.88497161865234,'*203': 67.88440704345703,'*204': 67.88438415527344,'*205': 67.88438415527344,'*206': 67.88438415527344,'*207': 67.88438415527344,'*208': 67.88438415527344,'*209': 67.88438415527344,'*21': 67.88999938964844,'*210': 67.88417053222656,'*211': 67.88417053222656,'*212': 67.88417053222656,'*213': 67.88417053222656,'*214': 67.88417053222656,'*215': 67.88417053222656,'*216': 67.88472747802734,'*217': 67.88472747802734,'*218': 67.88542175292969,'*219': 67.88542175292969,'*22': 67.88999938964844,'*220': 67.88502502441406,'*221': 67.88502502441406,'*222': 67.88502502441406,'*223': 67.88502502441406,'*224': 67.8849105834961,'*225': 67.88475036621094,'*226': 67.88475036621094,'*227': 67.88414764404297,'*228': 67.88443756103516,'*229': 67.88443756103516,'*23': 67.88999938964844,'*230': 67.88443756103516,'*231': 67.88507843017578,'*232': 67.88507843017578,'*233': 67.88499450683594,'*234': 67.88499450683594,'*235': 67.8851547241211,'*236': 67.8851547241211,'*237': 67.88540649414062,'*238': 67.88540649414062,'*239': 67.88505554199219,'*24': 67.88982391357422,'*240': 67.88528442382812,'*241': 67.88528442382812,'*242': 67.88528442382812,'*243': 67.88533020019531,'*244': 67.88581848144531,'*245': 67.88581848144531,'*246': 67.88581848144531,'*247': 67.88581848144531,'*248': 67.88581848144531,'*249': 67.88581085205078,'*25': 67.88949584960938,'*250': 67.88623809814453,'*251': 67.88623809814453,'*252': 67.88670349121094,'*253': 67.88667297363281,'*254': 67.88667297363281,'*255': 67.88690948486328,'*256': 67.88687896728516,'*257': 67.88607788085938,'*258': 67.88607788085938,'*259': 67.88607788085938,'*26': 67.88949584960938,'*260': 67.88607788085938,'*261': 67.88653564453125,'*262': 67.88653564453125,'*263': 67.88653564453125,'*264': 67.88653564453125,'*265': 67.88653564453125,'*266': 67.88668060302734,'*267': 67.88717651367188,'*268': 67.88722229003906,'*269': 67.88722229003906,'*27': 67.88949584960938,'*270': 67.88722229003906,'*271': 67.88756561279297,'*272': 67.88756561279297,'*273': 67.8876724243164,'*274': 67.8876724243164,'*275': 67.8876724243164,'*276': 67.88760375976562,'*277': 67.88844299316406,'*278': 67.88844299316406,'*279': 67.88835906982422,'*28': 67.88933563232422,'*280': 67.88835906982422,'*281': 67.88835906982422,'*282': 67.88835906982422,'*283': 67.88874053955078,'*284': 67.88874053955078,'*285': 67.88874053955078,'*286': 67.88874053955078,'*287': 67.88874053955078,'*288': 67.8885269165039,'*289': 67.88822174072266,'*29': 67.88933563232422,'*290': 67.88822174072266,'*291': 67.88841247558594,'*292': 67.88805389404297,'*293': 67.88805389404297,'*294': 67.88805389404297,'*295': 67.88805389404297,'*296': 67.88802337646484,'*297': 67.88714599609375,'*298': 67.88751983642578,'*299': 67.88763427734375,'*3': 67.88844299316406,'*30': 67.88933563232422,'*300': 67.88797760009766,'*301': 67.88797760009766,'*302': 67.8873291015625,'*303': 67.88766479492188,'*304': 67.88766479492188,'*305': 67.887451171875,'*306': 67.887451171875,'*307': 67.887451171875,'*308': 67.887451171875,'*309': 67.887451171875,'*31': 67.88933563232422,'*310': 67.88771057128906,'*311': 67.88711547851562,'*312': 67.88702392578125,'*313': 67.88702392578125,'*314': 67.88690185546875,'*315': 67.88690185546875,'*316': 67.88690185546875,'*317': 67.88690185546875,'*318': 67.88713073730469,'*319': 67.88713073730469,'*32': 67.88945770263672,'*320': 67.88755798339844,'*321': 67.88755798339844,'*322': 67.88755798339844,'*323': 67.88748168945312,'*324': 67.88748168945312,'*325': 67.8875503540039,'*326': 67.88845825195312,'*327': 67.88845825195312,'*328': 67.88845825195312,'*329': 67.88845825195312,'*33': 67.88945770263672,'*330': 67.88845825195312,'*331': 67.88845825195312,'*332': 67.88845825195312,'*333': 67.88880920410156,'*334': 67.88880920410156,'*335': 67.88907623291016,'*336': 67.88907623291016,'*337': 67.88907623291016,'*338': 67.88957214355469,'*339': 67.88957214355469,'*34': 67.88945770263672,'*340': 67.88961029052734,'*341': 67.888916015625,'*342': 67.888916015625,'*343': 67.88896942138672,'*344': 67.88858032226562,'*345': 67.88825988769531,'*346': 67.88825988769531,'*347': 67.88825988769531,'*348': 67.88825988769531,'*349': 67.88825988769531,'*35': 67.88904571533203,'*350': 67.88825988769531,'*351': 67.88851165771484,'*352': 67.88846588134766,'*353': 67.88846588134766,'*354': 67.88778686523438,'*355': 67.88778686523438,'*356': 67.8880386352539,'*357': 67.8880386352539,'*358': 67.8880386352539,'*359': 67.8880386352539,'*36': 67.88904571533203,'*360': 67.8880386352539,'*361': 67.88784790039062,'*362': 67.88784790039062,'*363': 67.88784790039062,'*364': 67.88813781738281,'*365': 67.88797760009766,'*366': 67.88792419433594,'*367': 67.88792419433594,'*368': 67.88792419433594,'*369': 67.8878173828125,'*37': 67.88904571533203,'*370': 67.88822937011719,'*371': 67.88822937011719,'*372': 67.88822937011719,'*373': 67.88822937011719,'*374': 67.88822937011719,'*375': 67.88800811767578,'*376': 67.88800811767578,'*377': 67.88800811767578,'*378': 67.88800811767578,'*379': 67.88800811767578,'*38': 67.88861846923828,'*380': 67.88800811767578,'*381': 67.88800811767578,'*382': 67.88800811767578,'*383': 67.88822174072266,'*384': 67.88822174072266,'*385': 67.88812255859375,'*386': 67.88761138916016,'*387': 67.88765716552734,'*388': 67.88765716552734,'*389': 67.88729858398438,'*39': 67.88861846923828,'*390': 67.88729858398438,'*391': 67.88772583007812,'*392': 67.88772583007812,'*393': 67.88772583007812,'*394': 67.8874282836914,'*395': 67.88814544677734,'*396': 67.88814544677734,'*397': 67.88814544677734,'*398': 67.88814544677734,'*399': 67.88814544677734,'*4': 67.88844299316406,'*40': 67.88861846923828,'*400': 67.88836669921875,'*401': 67.8883285522461,'*402': 67.8883285522461,'*403': 67.88897705078125,'*404': 67.88809967041016,'*405': 67.88770294189453,'*406': 67.88770294189453,'*407': 67.887939453125,'*408': 67.887939453125,'*409': 67.887939453125,'*41': 67.88861846923828,'*410': 67.88827514648438,'*411': 67.88827514648438,'*412': 67.88827514648438,'*413': 67.88827514648438,'*414': 67.88873291015625,'*415': 67.88903045654297,'*416': 67.88903045654297,'*417': 67.88903045654297,'*418': 67.8884048461914,'*419': 67.88810729980469,'*42': 67.888671875,'*420': 67.88841247558594,'*421': 67.88841247558594,'*422': 67.88841247558594,'*423': 67.88841247558594,'*424': 67.88876342773438,'*425': 67.88876342773438,'*426': 67.88822937011719,'*427': 67.88822937011719,'*428': 67.88822937011719,'*429': 67.88822937011719,'*43': 67.888671875,'*430': 67.88822937011719,'*431': 67.88794708251953,'*432': 67.88794708251953,'*433': 67.88794708251953,'*434': 67.88794708251953,'*435': 67.88660430908203,'*436': 67.88660430908203,'*437': 67.88660430908203,'*438': 67.88660430908203,'*439': 67.88660430908203,'*44': 67.88928985595703,'*440': 67.88660430908203,'*441': 67.88628387451172,'*442': 67.88628387451172,'*443': 67.88628387451172,'*444': 67.88628387451172,'*445': 67.88628387451172,'*446': 67.88628387451172,'*447': 67.88628387451172,'*448': 67.88628387451172,'*449': 67.88634490966797,'*45': 67.88932037353516,'*450': 67.8858413696289,'*451': 67.8858413696289,'*452': 67.88516998291016,'*453': 67.88536834716797,'*454': 67.88536834716797,'*455': 67.88536834716797,'*456': 67.88578796386719,'*457': 67.88578796386719,'*458': 67.8861312866211,'*459': 67.88607788085938,'*46': 67.88905334472656,'*460': 67.88607788085938,'*461': 67.8861083984375,'*462': 67.88568115234375,'*463': 67.88568115234375,'*464': 67.88568115234375,'*465': 67.8852310180664,'*466': 67.8852310180664,'*467': 67.8852310180664,'*468': 67.8852310180664,'*469': 67.8852310180664,'*47': 67.88905334472656,'*470': 67.8852310180664,'*471': 67.88560485839844,'*472': 67.88560485839844,'*473': 67.8852310180664,'*474': 67.88553619384766,'*475': 67.88593292236328,'*476': 67.88629913330078,'*477': 67.88629913330078,'*478': 67.88662719726562,'*479': 67.88603210449219,'*48': 67.88905334472656,'*480': 67.88603210449219,'*481': 67.88603210449219,'*482': 67.88577270507812,'*483': 67.88577270507812,'*484': 67.88578033447266,'*485': 67.88600158691406,'*486': 67.88600158691406,'*487': 67.88600158691406,'*488': 67.88600158691406,'*489': 67.8860092163086,'*49': 67.88905334472656,'*490': 67.88530731201172,'*491': 67.88530731201172,'*492': 67.88530731201172,'*493': 67.88530731201172,'*494': 67.88530731201172,'*495': 67.88530731201172,'*496': 67.88530731201172,'*497': 67.88445281982422,'*498': 67.88445281982422,'*499': 67.88445281982422,'*5': 67.88935852050781,'*50': 67.88932037353516,'*500': 67.88520050048828,'*501': 67.88520050048828,'*502': 67.88520050048828,'*503': 67.8858413696289,'*504': 67.8858413696289,'*505': 67.88623046875,'*506': 67.88623046875,'*507': 67.88623046875,'*51': 67.88932037353516,'*52': 67.8892593383789,'*53': 67.8892593383789,'*54': 67.8892593383789,'*55': 67.8892593383789,'*56': 67.88983154296875,'*57': 67.88948822021484,'*58': 67.88876342773438,'*59': 67.88876342773438,'*6': 67.88946533203125,'*60': 67.88812255859375,'*61': 67.88846588134766,'*62': 67.88780212402344,'*63': 67.8871841430664,'*64': 67.8871841430664,'*65': 67.8871841430664,'*66': 67.8871841430664,'*67': 67.8871841430664,'*68': 67.8871841430664,'*69': 67.88702392578125,'*7': 67.88978576660156,'*70': 67.88661193847656,'*71': 67.88661193847656,'*72': 67.88624572753906,'*73': 67.88660430908203,'*74': 67.88655090332031,'*75': 67.88655090332031,'*76': 67.88615417480469,'*77': 67.88612365722656,'*78': 67.88612365722656,'*79': 67.88612365722656,'*8': 67.88978576660156,'*80': 67.8858871459961,'*81': 67.8858871459961,'*82': 67.88557434082031,'*83': 67.88539123535156,'*84': 67.88539123535156,'*85': 67.88539123535156,'*86': 67.8851547241211,'*87': 67.8851547241211,'*88': 67.8851547241211,'*89': 67.88460540771484,'*9': 67.88904571533203,'*90': 67.88460540771484,'*91': 67.88460540771484,'*92': 67.88481140136719,'*93': 67.88510131835938,'*94': 67.88510131835938,'*95': 67.88561248779297,'*96': 67.88472747802734,'*97': 67.8847427368164,'*98': 67.8847427368164,'*99': 67.8847427368164}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed=self.filter_report(report)
        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, 'test_standard_cube_briggsbwtaper', self.parallel)

        test_dict['test_standard_cube_briggsbwtaper']['report'] = report
        test_dict['test_standard_cube_briggsbwtaper']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(image=img+'.image', range_list=[0.3, 1.0])
        self.mom8_creator(image=img+'.residual', range_list=[0.3, 1.0])
        test_dict['test_standard_cube_briggsbwtaper']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_standard_cube_briggsbwtaper']['images'].append(img+'.image.profile.png')


        if savemetricdict:
            ### serialize ndarray in mask_stats_dcit
            mask_stats_mod_dict = copy.deepcopy(mask_stats_dict)
            mask_stats_mod_dict['mask'] = mask_stats_dict['mask'].tolist()

            #create a nested dictionary containing exp dictionaries to save
            savedict = {}
            #list of stats to save
            # im_stats, mask_stats, pb_stats, psf_stats,\
            # model_stats, resid_stats, sumwt_stats]
            savedict['im_stats_dict']=im_stats_dict
            savedict['mask_stats_dict']=mask_stats_mod_dict
            savedict['pb_stats_dict']=pb_stats_dict
            savedict['psf_stats_dict']=psf_stats_dict
            savedict['model_stats_dict']=model_stats_dict
            savedict['resid_stats_dict']=resid_stats_dict
            savedict['sumwt_stats_dict']=sumwt_stats_dict
            
            #if self.parallel: 
            # for serial get it from *.psf
            savedict['bmin_dict']=bmin_dict
            savedict['bmaj_dict']=bmaj_dict
            savedict['pa_dict']=pa_dict

            self.save_dict_to_file(test_name, savedict, test_name+'_cas13317mod_stats')


        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cube_briggsbwtaper
#-------------------------------------------------#

    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mfs(self):
        ''' Standard (single field) MFS imaging - central field of NGC5363 (field 2), spw 16 & 22 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
            savemodel='none', parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[209.03000552deg, 5.25503742deg], [16.2902arcsec, 10.3226arcsec], 90.00000000deg]')

        # test_standard_mfs_exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 17.8415584564209],
            'com_bmin': [False, 10.3130340576],
            'com_pa': [False, 86.4390563965],
            'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
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
            'fit': [False, [0.0368173095435, 17.888484296, 9.90872728645]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.84024497577988],
            'fit_pix': [False, [40.2022706573, 40.0784833662]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 74, 0, 0]), \
                (img+'.image', False, [40, 75, 0, 0]), \
                (img+'.image', True, [6, 40, 0, 0]), \
                (img+'.image', False, [5, 40, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_standard_mfs_exp_mask_stats
        exp_mask_stats = {'npts': [True, 6400],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 330],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[209.03005051deg, 5.25476861deg], [71.9366arcsec, 71.6106arcsec], 0.00000000deg]')

        # test_standard_mfs_exp_pb_stats
        exp_pb_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'im_rms': [False, 0.577629540243],
            'npts_0.2': [True, 3793],
            'npts_0.5': [True, 1813],
            'npts_real': [True, 6400],
            'fit': [False, [1.0467417343495562, 92.30725376920157, \
                       92.30671415384658]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.84024497577988],
            'fit_pix': [False, [39.99973335198128, 40.00036927599604]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[209.03018549deg, 5.25490301deg], [18.5483arcsec, 11.7743arcsec], 90.00000000deg]')

        # test_standard_mfs_exp_psf_stats
        exp_psf_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [1e-10, 15849921197.895538],
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
            'fit': [False, [1.097246906267534, 15.626704258596684, \
                        9.180460042245928]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.84024497577988],
            'fit_pix': [False, [40.01095621317507, 39.995429898147734]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[209.02978056deg, 5.25512703deg], [18.0644arcsec, 11.9355arcsec], 90.00000000deg]')

        # test_standard_mfs_exp_resid_stats
        exp_resid_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
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

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.resid', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse [[209.02974119deg, 5.25476301deg], [2.7621arcsec, 1.8750arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_mfs_exp_model_stats
        exp_model_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [1e-10, 15849921197.895538],
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
            'mask_non0': [True, 0],
            'npts_real': [True, 6400]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_mfs_exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 15849921197.895538],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed=self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict['test_standard_mfs']['self.parallel'] = self.parallel
        test_dict['test_standard_mfs']['report'] = report
        test_dict['test_standard_mfs']['images'] = []


        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.003, 0.04])
        self.mom8_creator(img+'.residual', range_list=[-0.003, 0.04])
        test_dict['test_standard_mfs']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_mfs
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mtmfs(self):
        ''' Single field mtmfs imaging - central field of NGC5363 (field 2), spw 16 & 22 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
            savemodel='none', parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image.tt0', fit_region = \
            'ellipse[[209.03000552deg, 5.25503742deg], [16.2902arcsec, 10.3226arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 17.8415584564209],
            'com_bmin': [False, 10.3130340576],
            'com_pa': [False, 86.4390563965],
            'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
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
            'fit': [False, [0.03820503393292664, 18.02503733453906, \
                            9.894877124019276]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.84024497577988],
            'fit_pix': [False, [40.2022706573, 40.0784833662]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [40, 74, 0, 0]), \
                (img+'.image.tt0', False, [40, 75, 0, 0]), \
                (img+'.image.tt0', True, [6, 40, 0, 0]), \
                (img+'.image.tt0', False, [5, 40, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image.tt0', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_standard_mtmfs_exp_mask_stats
        exp_mask_stats = {'npts': [True, 6400],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 330],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', fit_region = \
            'ellipse[[209.03005051deg, 5.25476861deg], [71.9366arcsec, 71.6106arcsec], 0.00000000deg]')

        # test_standard_mtmfs_exp_pb_stats
        exp_pb_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.200896695256],
            'im_rms': [False, 0.577629540243],
            'npts_0.2': [True, 3793],
            'npts_0.5': [True, 1813],
            'npts_real': [True, 6400],
            'fit': [False, [1.0467417343495562, 92.30725376920157, \
                       92.30671415384658]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.84024497577988],
            'fit_pix': [False, [39.99973335198128, 40.00036927599604]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb.tt0', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', fit_region = \
            'ellipse[[209.03018549deg, 5.25490301deg], [18.5483arcsec, 11.7743arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_psf_stats
        exp_psf_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [1e-10, 15849921197.895538],
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
            'fit': [False, [1.097246906267534, 15.626704258596684, \
                        9.180460042245928]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.84024497577988],
            'fit_pix': [False, [40.01095621317507, 39.995429898147734]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf.tt0', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            'ellipse[[209.02978056deg, 5.25512703deg], [18.0644arcsec, 11.9355arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_resid_stats
        exp_resid_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
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

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', fit_region = \
            'ellipse[[209.02974119deg, 5.25476301deg], [2.7621arcsec, 1.8750arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_mtmfs_exp_model_stats
        exp_model_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [1e-10, 15849921197.895538],
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
            'mask_non0': [True, 0],
            'npts_real': [True, 6400]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        # test_standard_mtmfs_exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 15849921197.895538],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0', epsilon=self.epsilon)


        # .image.tt1 report
        im1_stats_dict = self.image_stats(img+'.image.tt1', fit_region = \
            'ellipse[[209.03000552deg, 5.25503742deg], [16.2902arcsec, 10.3226arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_im1_stats
        exp_im1_stats = {'com_bmaj': [False, 17.8415584564209],
            'com_bmin': [False, 10.3130340576],
            'com_pa': [False, 86.4390563965],
            'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.072932086885],
            'max_val_pos': [True, [43, 39, 0, 0]],
            'min_val': [False, -0.0771512910724],
            'min_val_pos': [True, [40, 62, 0, 0]],
            'im_rms': [False, 0.0222513486701],
            'im_sum': [False, 3.60833236173],
            'regn_sum': [False, 1.71812647558],
            'npts_real': [True, 6400]}

        report9 = th.check_dict_vals(exp_im1_stats, im1_stats_dict, '.image.tt1', epsilon=self.epsilon)

        # .residual.tt1 report
        resid1_stats_dict = self.image_stats(img+'.residual.tt1', fit_region = \
            'ellipse[[209.02978056deg, 5.25512703deg], [18.0644arcsec, 11.9355arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_resid1_stats
        exp_resid1_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 3793.0],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.000205082003959],
            'max_val_pos': [True, [44, 47, 0, 0]],
            'min_val': [False, -0.000261118228082],
            'min_val_pos': [True, [39, 41, 0, 0]],
            'im_rms': [False, 7.21058902237e-05],
            'im_sum': [False, 0.000510152031676],
            'regn_sum': [False, -0.0077287665232],
            'npts_real': [True, 6400]}

        report10 = th.check_dict_vals(exp_resid1_stats, resid1_stats_dict, \
            '.residual.tt1', epsilon=self.epsilon)

        # .model.tt1 report
        model1_stats_dict = self.image_stats(img+'.model.tt1', fit_region = \
            'ellipse[[209.02974119deg, 5.25476301deg], [2.7621arcsec, 1.8750arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_mtmfs_exp_model1_stats
        exp_model1_stats = {'npts': [True, 6400],
            'npts_unmasked': [True, 6400.0],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0282121878117],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 0.000369540677759],
            'im_sum': [False, 0.0370472576469],
            'regn_sum': [False, 0.0370472576469],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400]}

        report11 = th.check_dict_vals(exp_model1_stats, model1_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt.tt1 report
        sumwt1_stats_dict = self.image_stats(img+'.sumwt.tt1')

        # test_standard_mtmfs_exp_sumwt1_stats
        exp_sumwt1_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 15849921197.895538],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, -103907.34375],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, -103907.34375],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 103907.345804],
            'npts_real': [True, 1]}

        report12 = th.check_dict_vals(exp_sumwt1_stats, sumwt1_stats_dict, \
            '.sumwt.tt1', epsilon=self.epsilon)
        
        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11 + report12

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict['test_standard_mtmfs']['self.parallel'] = self.parallel
        test_dict['test_standard_mtmfs']['report'] = report
        test_dict['test_standard_mtmfs']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image.tt0', range_list=[-0.005, 0.04])
        self.mom8_creator(img+'.residual.tt0', range_list=[-0.005, 0.04])
        test_dict['test_standard_mtmfs']['images'].extend( \
            (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_mtmfs
#-------------------------------------------------#
    @stats_dict(test_dict)
    def test_standard_cube_eph(self):
        ''' Single field multi-EB ephemeris cube imaging - field 21PGiacobini-Zinner, spw 20 '''

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
            parallel=self.parallel, verbose=True)

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
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        exp_im_stats = {'com_bmaj': [False, 4.49769604687],
            'com_bmin': [False, 3.3237527868],
            'com_pa': [False, 87.0964067383],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.64398193359],
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
            'regn_sum': [False, 254.645320565],
            'npts_real': [True, 6400000],
            'profile': [False, 2.73209990935],
            'fit': [False, [3.04538752387499, 5.974552107890284, \
                    5.6824086756315]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.5049652049504],
            'fit_pix': [False, [45.76223486395211, 41.08882263728372]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 6400000],
            'freq_bin': [1e-10, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'mask_pix': [False, 9362],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.64398193359],
            'start': [True, 3.544453e+11],
            'end': [True, 3.545672e+11],
            'start_delta': [False, 3.544453e+11],
            'end_delta': [False, 3.545672e+11],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20036059618],
            'im_rms': [False, 0.57684102],
            'npts_0.2': [False, [3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233]],
            'npts_0.5': [False, [1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549]],
            'npts_real': [True, 6400000],
            'fit': [False, [1.046847676114786, 28.075049566294457, \
                    28.075184571520158]],
            'fit_loc_chan': [True, 500],
            'fit_loc_freq': [1e-10, 354.5063079930342],
            'fit_pix': [False, [40.0, 40.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122071.64398193359],
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
            'fit_0': [False, [1.1041815481079804, 3.942786862512436, \
                      2.861422654767207]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.4453942426872],
            'fit_pix_0': [False, [39.99991697297065, 39.997612723971876]],
            'fit_1': [False, [1.1041897949047454, 3.943995482479015, \
                      2.858120721885849]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.5063079930342],
            'fit_pix_1': [False, [39.999914819825086, 39.997626281193952]],
            'fit_2': [False, [1.1041648382099558, 3.946029935006709, \
                      2.8548012968092817]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672217433812],
            'fit_pix_2': [False, [39.99991098336098, 39.99764268795339]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.64398193359],
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
            'regn_sum': [False, 50.9751291326],
            'npts_real': [True, 6400000]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122071.64398193359],
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
            'regn_sum': [False, 53.3284724131],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400000]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [1e-10, 122071.64398193359],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict['test_standard_cube_eph']['self.parallel'] = self.parallel
        test_dict['test_standard_cube_eph']['report'] = report
        test_dict['test_standard_cube_eph']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[0.0, 3.25])
        self.mom8_creator(img+'.residual', range_list=[0.0, 3.25])
        test_dict['test_standard_cube_eph']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_standard_cube_eph']['images'].append(img+'.image.profile.png')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of est_standard_cube_eph
#-------------------------------------------------#
    @stats_dict(test_dict)
    def test_standard_cube_eph_briggsbwtaper(self):
        ''' Single field multi-EB ephemeris cube imaging with briggsbwtaper - field 21PGiacobini-Zinner, spw 20 '''

        test_name = self._testMethodName
        file_name = 'standard_cube_eph_briggsbwtaper.iter'
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
            nchan=1000, start=1550, width=1, perchanweightdensity=True, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            usepointing=False, pblimit=0.2, deconvolver='hogbom', \
            restoration=False, restoringbeam='common', pbcor=False, \
            weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, savemodel='none', \
            parallel=self.parallel, verbose=True)

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
            nchan=1000, start=1550, width=1, perchanweightdensity=True, \
            gridder='standard', chanchunks=-1, mosweight=False, \
            usepointing=False, pblimit=0.2, deconvolver='hogbom', \
            restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggsbwtaper', robust=0.5, npixels=0, niter=30000, \
            threshold='0.274Jy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            calcres=False, calcpsf=False, savemodel='none', \
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 4.6153633847721665],
            'com_bmin': [False, 3.3352055145021593],
            'com_pa': [False, 84.78377068098479],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 3.278780937194824],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.36217406392097473],
            'min_val_pos': [True, [58, 39, 0, 491]],
            'im_rms': [False, 0.05609585143475306],
            'rms_per_chan': [False, [0.04980624428087379, 0.049181907069306884, 0.04801920643963151, 0.04583301914121372, 0.040720229967857464, 0.042013242218929284, 0.044496877201658624, 0.04355822791216707, 0.0449778600808559, 0.046856256495317716, 0.04525944482302899, 0.045430809986070304, 0.04299310281523898, 0.04186902570926661, 0.04043084679940776, 0.04545874813352345, 0.04364338087203972, 0.04315687368960783, 0.04018915741849983, 0.0423976686635374, 0.04566862276957951, 0.044980425608893815, 0.047986799828113234, 0.0454670082067846, 0.04426951895873982, 0.041307588444546905, 0.04018901097824569, 0.04913337295576702, 0.04483536246116845, 0.04414032475776276, 0.04752608973278683, 0.04657242307929253, 0.04793042158708105, 0.046760752326240204, 0.047214975532489274, 0.04839378091520341, 0.04825546043735031, 0.04816914656949764, 0.04600274306879198, 0.04130966699659313, 0.040011129111791784, 0.04165829187993401, 0.04486324130951055, 0.04220008804332373, 0.04022094813980892, 0.041294389596907316, 0.04025890515206982, 0.040833357228103564, 0.04426716647794076, 0.04599827748439303, 0.0424597022226665, 0.04192993720032838, 0.04438428153293445, 0.04897413600810136, 0.04517680416565816, 0.04296266432388383, 0.046479441601272604, 0.04726413407598469, 0.04410408431676381, 0.04134281936058372, 0.04093760864821317, 0.040132813240349764, 0.046391849904873765, 0.04897015487916892, 0.04734767804137239, 0.047270390127751416, 0.04736645605197829, 0.04495719984099806, 0.045496968493763904, 0.04500780177515373, 0.04258769147994997, 0.04486979757351609, 0.049459663707374674, 0.05043497955453825, 0.04956427350050683, 0.04500281929508872, 0.04515290050015542, 0.045714135571707136, 0.04469620538913332, 0.04835843054775528, 0.04715527152773822, 0.04448291905233272, 0.04244176276753743, 0.041084170210988875, 0.04404038926778707, 0.045362341517562056, 0.04583488732018053, 0.04469146348477125, 0.04663724370917155, 0.04566821915138833, 0.04245169387859445, 0.0412800342905669, 0.042486984976373106, 0.040813668303577175, 0.041808431005998015, 0.045121108305867255, 0.045385506017004876, 0.05118439007550736, 0.05320404309166031, 0.04548350508181054, 0.04118184049955647, 0.0411512750552513, 0.04493689285128879, 0.04870190853374057, 0.047284100835281756, 0.04631604803368088, 0.04754751972682775, 0.04675224540819776, 0.04416995141663751, 0.04642057015211928, 0.046011961665140474, 0.04032335026726956, 0.03657195169867859, 0.04040259153083771, 0.04909800863415589, 0.04646050373013411, 0.044957052801354105, 0.04207719445895267, 0.04606668497034549, 0.04524264553097416, 0.04544479089911701, 0.046936269184580486, 0.045706196311616404, 0.04353111520262883, 0.042313590377506476, 0.04482796667681791, 0.047613409807898445, 0.045693774297831384, 0.042940222292323305, 0.04441896756335694, 0.043263604081920426, 0.042026692593400854, 0.04544804217903177, 0.05011227353998639, 0.04816373413542914, 0.04138134768986196, 0.04355039695416439, 0.042681191552351436, 0.040408007255303545, 0.03992214286466301, 0.04173199051457398, 0.044210135133686816, 0.04294243138754665, 0.0422144351238895, 0.04091376046834093, 0.038239476687855954, 0.039478720169426565, 0.04305181582316579, 0.04406747870880562, 0.04055128740308246, 0.041638453836799344, 0.040455903262116934, 0.04659513423728624, 0.04834144721794595, 0.04643827529763154, 0.044026559129439216, 0.04197465211969988, 0.043893898715026365, 0.042919438548186375, 0.046964700377572784, 0.04891064676024664, 0.05084393677414183, 0.044932369275116006, 0.04313927419530443, 0.04786708903327418, 0.04307081434153513, 0.044417635039119246, 0.04894970938107093, 0.04378239901714948, 0.04453574987340776, 0.04825027081431791, 0.04296916121670574, 0.03815406804388665, 0.04303894848673802, 0.04527881738848557, 0.040376351107935436, 0.040056944612741806, 0.04345827760108247, 0.046334453557989194, 0.0488763264857944, 0.044697620657334944, 0.04075277212488179, 0.041123039367720715, 0.0401143358235041, 0.04485666839208719, 0.045301884409551886, 0.04028040261156905, 0.039716018847840914, 0.040527497522189926, 0.04632675276554114, 0.04700324442344388, 0.04592163974549857, 0.045171470052460215, 0.04012238355473271, 0.040455031735650404, 0.03999991868628664, 0.04194371888685938, 0.04311115896250561, 0.04408706358470665, 0.04347547496802239, 0.04203740289513717, 0.043263652478779994, 0.042166523855643914, 0.04352398978783646, 0.045461953083174536, 0.04419560398392808, 0.04421756495230233, 0.04022396265467646, 0.041524512027132705, 0.04565787052492574, 0.0500315042238843, 0.0445624797646394, 0.04319203811837168, 0.04887012919722226, 0.04883433827377812, 0.048631448744955674, 0.04475117265182139, 0.040963688594120674, 0.042319821523143196, 0.04462708353905811, 0.04313569364380933, 0.043846266999300984, 0.04669445958770959, 0.04548152335798154, 0.04430707277355682, 0.04422960089101078, 0.040764386740285416, 0.04336192842690154, 0.04229400254760735, 0.04442315057589513, 0.044068365596022445, 0.040604459568553845, 0.049625846821818864, 0.05021146088317837, 0.044999920343737065, 0.0457526104623044, 0.045032684477328874, 0.044740144700677296, 0.042239167734301315, 0.04509271724470693, 0.05006398482212349, 0.04656892447564126, 0.04572764977062891, 0.04628712624582241, 0.04464059279577941, 0.04698171445665529, 0.04953811270248742, 0.04718508503309449, 0.04196182626509631, 0.04692870733319241, 0.04463895881309594, 0.041878529431377064, 0.043762869106085484, 0.044432708185992355, 0.048273794856721196, 0.04795907402179907, 0.0422627028029166, 0.042199185655552494, 0.043924713614789335, 0.04444487724298314, 0.044563511785574224, 0.037535200343673004, 0.03936505040178308, 0.04510833595809811, 0.045553055990566235, 0.04478751011683643, 0.045386297064108946, 0.0432097877183372, 0.03985998576652175, 0.041022194738624444, 0.048654312884396596, 0.050985261529745696, 0.04542321727758847, 0.04160751407283635, 0.041281980936247506, 0.04271397441788915, 0.0418983775040208, 0.04044934706943166, 0.043006246333914894, 0.04609117451494283, 0.04546556823371546, 0.04221481597863437, 0.04041883750147706, 0.04770925430227707, 0.04823457263527067, 0.04726978377651058, 0.045106492474668525, 0.04225757845669231, 0.04678933310799391, 0.04787018414352443, 0.04108542457002423, 0.037331793915892274, 0.039019595089095034, 0.04144712138705571, 0.0426892476027771, 0.04202044133898991, 0.038951585411872006, 0.043542092485884316, 0.0468130417689795, 0.04526892805824094, 0.043783353777364344, 0.04642090076874412, 0.04639001862061178, 0.04123936755154233, 0.038114589767314314, 0.04190173307200588, 0.042634155404572376, 0.043529360933599895, 0.04540341847809748, 0.04387585155078013, 0.044515620896691095, 0.04372826970193465, 0.04230511425935403, 0.0462100510157415, 0.04335007140884004, 0.04384441184561817, 0.04757563232289588, 0.0448257469648976, 0.04211644847771746, 0.044154843109719705, 0.047906395764182076, 0.046817951635192694, 0.04685022216333965, 0.04620850864383186, 0.04109408941905257, 0.0415755776179549, 0.03986359408286068, 0.04270506604206936, 0.042762240419396066, 0.043714811486360076, 0.046760273258416134, 0.046310524876293786, 0.04186714820049629, 0.03942194578879492, 0.0408324088875011, 0.04066367375963188, 0.04098352805910296, 0.03845156034620269, 0.03758989617101229, 0.04226734224030958, 0.04446549313470781, 0.04323417229315278, 0.0433379273841091, 0.043231797417263775, 0.04471202345517597, 0.04404202140796481, 0.04559266883750324, 0.04511132921388279, 0.04783327174710853, 0.04785175235414532, 0.044080125813153424, 0.043563356500855535, 0.045228365893599326, 0.048119001847574415, 0.04516816145742154, 0.04375121280277562, 0.04270026510368359, 0.043416465093796645, 0.04435849706295363, 0.04714260630937985, 0.04587270438917867, 0.042613118958383364, 0.04166664124631433, 0.042658540605232505, 0.04242970607215766, 0.04569840307168379, 0.043884801262789944, 0.04368167828018671, 0.04616796550154052, 0.04557970819782034, 0.04596346940429089, 0.04744308632649514, 0.04581887153316369, 0.04791089473964647, 0.047974396989098764, 0.04295641770163905, 0.04037239375700781, 0.041736885828141494, 0.04629613476951637, 0.047360285637575744, 0.043680772307138925, 0.043248961808992986, 0.04684372667684942, 0.04453332754874046, 0.04410262135737886, 0.041115401967210745, 0.045490448195516564, 0.04941202988368553, 0.04561940652718864, 0.044565646356953946, 0.046990762276144335, 0.04188813097523235, 0.04043355352843307, 0.041730113847184454, 0.0447412248222876, 0.045212858836236616, 0.0446284276843788, 0.04229886045101561, 0.04339854974188936, 0.040571665304725446, 0.03838944654011296, 0.04088656460980091, 0.04180508809389838, 0.041185700819891265, 0.04327764421717837, 0.04554877162429928, 0.044296690188205735, 0.04548400765872825, 0.04807068118803301, 0.0462220295601339, 0.04434384783873343, 0.04207745772153928, 0.044530088703203935, 0.04599698242383248, 0.04792550783540083, 0.05023303339556688, 0.050942453049856035, 0.04888602392069365, 0.04357286831307993, 0.03750278441040552, 0.03846888223164225, 0.040712385191687246, 0.045182149891697714, 0.047869204815265534, 0.04323777491157295, 0.03920736450460925, 0.03959688656450024, 0.04539319679730105, 0.04790311517905508, 0.04269868015428697, 0.04613374136091106, 0.0425781691476122, 0.04249090023861118, 0.044113231973184466, 0.04299976873524666, 0.04582150130401724, 0.04891562008299229, 0.04492345452579839, 0.0409978218224756, 0.04243896965319094, 0.04604525804736929, 0.04396477249394169, 0.040130762446879485, 0.041822963928351366, 0.04215982196546785, 0.04164974805642088, 0.04536506461438564, 0.04777626706375024, 0.0431995035900221, 0.04354770707684929, 0.03895193340190124, 0.0389156484608813, 0.04238502096732911, 0.04424178492118985, 0.042052114010331175, 0.03972915428556974, 0.041546459434003466, 0.045020494797250295, 0.04718692080579099, 0.04451399601866224, 0.04366346105099904, 0.04365130489179828, 0.0395858634667152, 0.04396755835713188, 0.04461707938196565, 0.04604701621952156, 0.04783633452172411, 0.044758238154428094, 0.043063252861316674, 0.04201998361385465, 0.04073877380625209, 0.044402546913336734, 0.04222212984362886, 0.040332331390910005, 0.04861910094838351, 0.04380241485836759, 0.04215905929326605, 0.04338423844872642, 0.04403204676267041, 0.04497941798317015, 0.04578002115066703, 0.04470870953239883, 0.044184817977881034, 0.04177959133841808, 0.04387369068327449, 0.07305277156729907, 0.15544584137364367, 0.2702652737180338, 0.36394865385458536, 0.3865861515727367, 0.35011843468077647, 0.3041555482433153, 0.27664007072652724, 0.2564945596114512, 0.24787542897601644, 0.25251243580670746, 0.2656525510825005, 0.273718878392953, 0.293368866444176, 0.27651661697780644, 0.2101462791294832, 0.12929803220875385, 0.06723521727193423, 0.04378645674190316, 0.04314026444006957, 0.043367676477108194, 0.039959769938721454, 0.04370069001097637, 0.04821575247345042, 0.047693550789981685, 0.042152038442580135, 0.04157198737558815, 0.043403997121955015, 0.044507810107929244, 0.045794821826312865, 0.04955539940000539, 0.04664636113524865, 0.04570109643352515, 0.04591288430141476, 0.04234895999184957, 0.040750154354899266, 0.04444706083068934, 0.047600928243617315, 0.04711862728547932, 0.04572248410505566, 0.043822096253733525, 0.042020072499529604, 0.04670420063133669, 0.05206483244905589, 0.05212855409861097, 0.04915852922859097, 0.042879350068720554, 0.0416947473398087, 0.04533147123562846, 0.04941893827054203, 0.05149789903730183, 0.04840108966379967, 0.04546533811938858, 0.04446842686277916, 0.04355172925910375, 0.03807080320216515, 0.04204923252887683, 0.044675315898427624, 0.04489718537370469, 0.04480407160813274, 0.041420989758451564, 0.04171991898273599, 0.04194691197686678, 0.046366480009900436, 0.045797837338140586, 0.045539396058395715, 0.04525346407039292, 0.043591190311794946, 0.04349976562696881, 0.044149146210975036, 0.04273249108667336, 0.0395596337382392, 0.040493294537671386, 0.0453563375403294, 0.04435984275749569, 0.04192004894088999, 0.04246829130561042, 0.04461554331472894, 0.04174049759567807, 0.042534370275570166, 0.04169512036572455, 0.0442550106160829, 0.04487152187633878, 0.04674051544288151, 0.04241417773221762, 0.03951044660265742, 0.04335366868249594, 0.046743559350610746, 0.049836378625525694, 0.046185606692434045, 0.05068954907741317, 0.04973027133995008, 0.0391829648914183, 0.044119508227948384, 0.049653930129857304, 0.047451272723996823, 0.047695508241589746, 0.04749783571741011, 0.0453621862234796, 0.04594710933319943, 0.04355776792284987, 0.0407931770593412, 0.04035841817514514, 0.042302749889980715, 0.0429730130311627, 0.04466129565183151, 0.0445906243572373, 0.04097757673706596, 0.04156754572176273, 0.04273041960578711, 0.04813130436820179, 0.04207531116193894, 0.043589550286053796, 0.048831511214735224, 0.0430989582948261, 0.03856286877150434, 0.041245202425397726, 0.039567194657886676, 0.04054837100494838, 0.043805356444913464, 0.04645588227596993, 0.043426898199757945, 0.043083058854070454, 0.04332663286330242, 0.04277172448425333, 0.042336170052119156, 0.04181620443825819, 0.04574395530764865, 0.04527189771769063, 0.040782663128380116, 0.04088182461712406, 0.04191171199030378, 0.04384235896350345, 0.045395326132367256, 0.046818529164889884, 0.04638858233497556, 0.045871078289963174, 0.04344631259960327, 0.04193551291577763, 0.04180431133476504, 0.0373889774021823, 0.043491254415707954, 0.0505776020810951, 0.046934136196861424, 0.04628396447397184, 0.0439071579207898, 0.04408112598825481, 0.04345447636536801, 0.04011749585977608, 0.03905179924743059, 0.04192455414683994, 0.041773625875387306, 0.03971202674626663, 0.04147221174774904, 0.037676015964510694, 0.03743411557331315, 0.03968387227174169, 0.04168907741599824, 0.044602389235800476, 0.04522250935903086, 0.04515088409162237, 0.04675339295431293, 0.05177600964111318, 0.046269431292066655, 0.039923211834322035, 0.03996613798731274, 0.04426859312062878, 0.04370530624251765, 0.0489219749502112, 0.04783026017783438, 0.046416058275799345, 0.04530957136206231, 0.04421265544586097, 0.045168836978819994, 0.048309304212433596, 0.04837500109115457, 0.04506961145775217, 0.04112117099214982, 0.0426200258905851, 0.0434544251038104, 0.040377733070579486, 0.04264081766732213, 0.04354552534640634, 0.04783279481664696, 0.046905878521101736, 0.041332753980774406, 0.04098873033435018, 0.04682343103683209, 0.050759621938972384, 0.047107665961214476, 0.04245807805955911, 0.04087450375000468, 0.03831741518169265, 0.041921230773073026, 0.04689877060521487, 0.045988513937869256, 0.04795145527758773, 0.04511493075004005, 0.0403257594695132, 0.03796384782585744, 0.038451787457557215, 0.04435249935659175, 0.04359665587861431, 0.04419658898314747, 0.04216633482535092, 0.04442145448446196, 0.04444878902801034, 0.047233945779062605, 0.046777020540399644, 0.043673788010065995, 0.04019670505866605, 0.04232770491972007, 0.03752652573763034, 0.043486972021560245, 0.04654748847019314, 0.04318571445814521, 0.04389362194631678, 0.04153877008823151, 0.039093455486948946, 0.03940525581755273, 0.03958009551565283, 0.04271191651184301, 0.04106679153998524, 0.04087979891866519, 0.03875213964036597, 0.040902432104083876, 0.04205056723519903, 0.04432166442611053, 0.04756886170640058, 0.042596461904025794, 0.03978341644117022, 0.04343789460686302, 0.044419990244377115, 0.04315594594446142, 0.04650974276194317, 0.041985616855391664, 0.04414012991693783, 0.0413698436084358, 0.04043295125702463, 0.0411724489679433, 0.0449391352313776, 0.048601670089105355, 0.04709151565495982, 0.04198885837966443, 0.040394164303216316, 0.041183022670481524, 0.04497884529488561, 0.042733375614482344, 0.04263955184690891, 0.04587598190839218, 0.04748222235105568, 0.04522195012776497, 0.042534521572661216, 0.040610225980800045, 0.04349224646776112, 0.041254512035441684, 0.04210296158624721, 0.03887546899791237, 0.036808562495817106, 0.03912711179373199, 0.04928362259660566, 0.051548360721369986, 0.04272316663205788, 0.041620355450248094, 0.0481752738918194, 0.04687823832641236, 0.041563827594157304, 0.04101096181831323, 0.04299126782045179, 0.04610885210981244, 0.04392945000425368, 0.045040843502235375, 0.04741909054042633, 0.04633257589646319, 0.044040106534833576, 0.04007295666212543, 0.04144794392147843, 0.04211471299453388, 0.04273486862357951, 0.04583606325976122, 0.04457857056146514, 0.043169259177076064, 0.0415055680254373, 0.041587660676987916, 0.04469126758993004, 0.04574906484239683, 0.04161011710388962, 0.04384982503710152, 0.04616737105612526, 0.048353014951543, 0.04767659321729997, 0.04642393878703179, 0.0477388248362813, 0.047742871238186525, 0.045463015898933326, 0.04397636666246282, 0.04310492961261954, 0.04068791545197722, 0.04227921073297164, 0.041218937716768815, 0.0451556729966859, 0.04773556868704452, 0.04391115465060777, 0.04245659968616945, 0.04671939766349839, 0.04430279724961868, 0.04959418691923308, 0.052171974246664, 0.046174050462525026, 0.04372705821927368, 0.041821941949304574, 0.04352171765646556, 0.04667720830043348, 0.045785286838259776, 0.04529381420569239, 0.045264332256172914, 0.044862389147434435, 0.04389034552873384, 0.044589549523640686, 0.04109548487995436, 0.04162204076477356, 0.03984821115729712, 0.041939249507135797, 0.046242612670249934, 0.044646421091092944, 0.04365607347289076, 0.0421870981331197, 0.03993297114160635, 0.03946903193000705, 0.04029320566599106, 0.04552046168518562, 0.04813756924553471, 0.04705332993721689, 0.04841216473592258, 0.050197538140956205, 0.044296817001627097, 0.042485239459908695, 0.04391362929237573, 0.04132343869699482, 0.03994493615398885, 0.04494611359820636, 0.04537923226960366, 0.039830096078664126, 0.03852482278212683, 0.04071752740657917, 0.04187476389396255, 0.04247533417724338, 0.0429865340997313, 0.044228192126820456, 0.043945698523526, 0.03999219784680604, 0.042079374275514, 0.04488340639022528, 0.042060032208827286, 0.041530847120592934, 0.0488975188690411, 0.045881010926442534, 0.04287253550957093, 0.041926531339338095, 0.040231028971184406, 0.041257347868347054, 0.04501465599440565, 0.044992431432709236, 0.04200578222823483, 0.04190542554778625, 0.04019192260807578, 0.039140793287919806, 0.038744433317932904, 0.038148439569082336, 0.040968989845211114, 0.042164169786815156, 0.041270726307542444, 0.040222078452398835, 0.042256364473736575, 0.04032293388574375, 0.04101505604183243, 0.04191052762997484, 0.04528139294164947, 0.04777875066205739, 0.04837678048483792, 0.04972941651706581, 0.04640671671485831, 0.04113848247215838, 0.04058421572456904, 0.043038012396567064, 0.04171807008244231, 0.042132025635814224, 0.04539782448003977, 0.04781995564681897, 0.04853817461178558, 0.04859688563224961, 0.044902809524043595, 0.04189588138450321, 0.03995813525544304, 0.04074973490805897, 0.042752876538648094, 0.04357136448441821, 0.042756583672004314, 0.04222607115959047, 0.042901091292005415, 0.03998762820347158, 0.037664109896116175, 0.03757397694224933, 0.03971212078289213, 0.04424603781821552, 0.046232943395224, 0.04899749585053205, 0.04792630397887076, 0.042492823926643865, 0.041992887252723016, 0.042132046425240856, 0.03955067678224552, 0.04053676903254468, 0.042427559989564254, 0.04507214772981009, 0.042755877711259736, 0.04575481866837017, 0.04140646521993148, 0.03844441894472225, 0.03652163449465, 0.038052447790345945, 0.04368156350072061, 0.04261886474535546, 0.04451026536552225, 0.04444789906772739, 0.039927410085489494, 0.044273417121475825, 0.048012317715860585, 0.04172341237241333, 0.04078731341972721, 0.04307378570039496, 0.040218440461566876, 0.038857491538841306, 0.04190966270080123, 0.043585315489432945, 0.04032241257291049, 0.0449185454133472, 0.045616086466998156, 0.043874570286642475, 0.04876542469087236, 0.046319177207382725, 0.039917657544146785, 0.03886210248944117, 0.037495455623999516, 0.04327468446178233, 0.050035893725219036, 0.05073699844277334, 0.04752767717341416, 0.04393040003149914, 0.03702360612120744, 0.039809385494204375, 0.04092348485546371, 0.04155343225203003, 0.04282200096205279, 0.04119349856603349, 0.04149071297094778, 0.046441342718078336, 0.048362296929523416, 0.04891974439045019, 0.047445213050046564, 0.046323392621991505, 0.04421522337142135, 0.04437017903964814, 0.04709756459370677, 0.04546990259679002, 0.041376587315531774, 0.04263876974146591, 0.04427113354592288, 0.045323262391820225, 0.046248307276283346, 0.046844671815476524, 0.04849616244091082, 0.04319703638027593, 0.041192175821505285, 0.03747596736562667, 0.03811028618410946, 0.042800443034937505, 0.044180721374921884, 0.0422298003416534, 0.042464250469215634, 0.04250350336742838, 0.04729958460723437, 0.04952695203996929, 0.048979869530825955, 0.04277116592012446, 0.04015436011010988, 0.0446207539506111, 0.047771508304026745, 0.04195450094591805, 0.04485475514826325, 0.04732289226790657, 0.0485766334562288, 0.04405986312529342, 0.046542572386518745, 0.047094421010156935, 0.04666226192461151, 0.045428735940816534, 0.047162925646747324, 0.05200163240986284, 0.0478373283231855, 0.039233559611014414, 0.04303906667316332, 0.04225846037237894, 0.04586998644480703, 0.04706375833764583, 0.043672319069770905, 0.0461561805056443, 0.04539671333506246, 0.043145669144071884, 0.04143930371848491, 0.04098195247072255, 0.04498886325541145, 0.04489405486491206, 0.03762505349731428, 0.04451778288746828]],
            'im_sum': [False, 2550.144158477671],
            'regn_sum': [False, 266.1940816491842],
            'npts_real': [True, 6400000],
            'profile': [False, 2.765215599712529],
            'fit': [False, [3.1088767528803665, 6.085664910672221, 5.757693572227245]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.5049652049504],
            'fit_pix': [False, [45.67982128909444, 41.05347470202845]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_standard_cube_eph_briggsbwtaper.exp_mask_stats
        exp_mask_stats = {'npts': [True, 6400000],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'mask_pix': [False, 9311],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20036059617996216],
            'im_rms': [False, 0.5768410193579128],
            'npts_0.2': [False, [3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233]],
            'npts_0.5': [False, [1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549]],
            'npts_real': [True, 6400000],
            'fit': [False, [1.046847676114786, 28.075049566294457, 28.075049566280292]],
            'fit_loc_chan': [True, 500],
            'fit_loc_freq': [1e-10, 354.5063079930342],
            'fit_pix': [False, [40.0, 40.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.15693220496177673],
            'min_val_pos': [True, [34, 36, 0, 949]],
            'im_rms': [False, 0.08721424406426508],
            'im_sum': [False, -1836.9928123016398],
            'npts_real': [True, 6400000],
            'fit_0': [False, [1.095680038820762, 4.101209438181028, 2.9137339402226234]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.4453942426872],
            'fit_pix_0': [False, [40.00006531772031, 39.99740702688986]],
            'fit_1': [False, [1.0956946363854294, 4.100629995250517, 2.913432511720987]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.5063079930342],
            'fit_pix_1': [False, [40.00006539731349, 39.997408575244705]],
            'fit_2': [False, [1.0957370794438608, 4.1001549993337605, 2.9123273283427773]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672217433812],
            'fit_pix_2': [False, [40.000066216213064, 39.99741175467212]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 0.32706770300865173],
            'max_val_pos': [True, [59, 65, 0, 493]],
            'min_val': [False, -0.33076247572898865],
            'min_val_pos': [True, [32, 47, 0, 491]],
            'im_rms': [False, 0.04607551972676364],
            'im_sum': [False, 398.8007672918444],
            'regn_sum': [False, 51.87397705297917],
            'npts_real': [True, 6400000]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_eph_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.4025037288665771],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.07458142936229706],
            'min_val_pos': [True, [58, 39, 0, 491]],
            'im_rms': [False, 0.001708715350338477],
            'im_sum': [False, 53.728492775931954],
            'regn_sum': [False, 54.275748455896974],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400000]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_cube_eph_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1134.3048095703125],
            'max_val_pos': [True, [0, 0, 0, 70]],
            'min_val': [False, 1133.8365478515625],
            'min_val_pos': [True, [0, 0, 0, 666]],
            'im_rms': [False, 1134.0983575620767],
            'npts_real': [True, 1000]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict['test_standard_cube_eph_briggsbwtaper']['self.parallel'] = self.parallel
        test_dict['test_standard_cube_eph_briggsbwtaper']['report'] = report
        test_dict['test_standard_cube_eph_briggsbwtaper']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[0.0, 3.25])
        self.mom8_creator(img+'.residual', range_list=[0.0, 3.25])
        test_dict['test_standard_cube_eph_briggsbwtaper']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_standard_cube_eph_briggsbwtaper']['images'].append(img+'.image.profile.png')

        if savemetricdict:
            ### serialize ndarray in mask_stats_dcit
            mask_stats_mod_dict = copy.deepcopy(mask_stats_dict)
            mask_stats_mod_dict['mask'] = mask_stats_dict['mask'].tolist()
            #create a nested dictionary containing exp dictionaries to save
            savedict = {}
            #list of stats to save
            # im_stats, mask_stats, pb_stats, psf_stats,\
            # model_stats, resid_stats, sumwt_stats]
            savedict['im_stats_dict']=im_stats_dict
            savedict['mask_stats_dict']=mask_stats_mod_dict
            savedict['pb_stats_dict']=pb_stats_dict
            savedict['psf_stats_dict']=psf_stats_dict
            savedict['model_stats_dict']=model_stats_dict
            savedict['resid_stats_dict']=resid_stats_dict
            savedict['sumwt_stats_dict']=sumwt_stats_dict
            
            self.save_dict_to_file(test_name, savedict, test_name+'_cas13317mod_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cube_eph_briggsbwtaper
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mfs_eph(self):
        ''' Standard (single field) ephemeris mfs imaging - central field of Venus (field 2), spw 25 & 45 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel, verbose=True)

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
            savemodel='none', parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_im_stats = {'com_bmaj': [False, 0.875946879387],
            'com_bmin': [False, 0.673672378063],
            'com_pa': [False, 88.5368652344],
            'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 82944],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 82944]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[239.37090838deg, -16.96415647deg], [17.8437arcsec, 17.4772arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, 0.200061768293],
            'im_rms': [False, 0.569403001285],
            'npts_0.2': [1e-4, 47329],
            'npts_0.5': [1.5e-4, 22365],
            'npts_real': [True, 82944],
            'fit': [False, [1.0286609217550002, 22.907692916947163, \
                       22.90769291676479]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442221593894],
            'fit_pix': [False, [144.0, 144.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.37094084deg, -16.96415506deg], [1.1279arcsec, 0.7875arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 82944.0],
            'freq_bin': [1e-10, 16762501225.396851],
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
            'fit': [False, [0.9200466881631709, 0.9746655722260728, \
                        0.7626550313652652]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442221593894],
            'fit_pix': [False, [144.00051463175717, 144.00004766689185]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 82944.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [1e-4, 1.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict['test_standard_mfs_eph']['self.parallel'] = self.parallel
        test_dict['test_standard_mfs_eph']['report'] = report
        test_dict['test_standard_mfs_eph']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-1.05, 1.05])
        self.mom8_creator(img+'.residual', range_list=[-1.05, 1.05])
        test_dict['test_standard_mfs_eph']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_mfs_eph
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mtmfs_eph(self):
        ''' Standard (single field) ephemeris mtmfs imaging - central field of Venus (field 2), spw 25 & 45 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel, verbose=True)

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
            savemodel='none', parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mtmfs'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image.tt0', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_im_stats = {'com_bmaj': [False, 0.875946879387],
            'com_bmin': [False, 0.673672378063],
            'com_pa': [False, 88.5368652344],
            'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 82944],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 82944]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', fit_region = \
            'ellipse[[239.37090838deg, -16.96415647deg], [17.8437arcsec, 17.4772arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, 0.200061768293],
            'im_rms': [False, 0.569403001285],
            'npts_0.2': [1e-4, 47329],
            'npts_0.5': [1.5e-4, 22365],
            'npts_real': [True, 82944],
            'fit': [False, [1.0286609217550002, 22.907692916947163, \
                       22.90769291676479]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442221593894],
            'fit_pix': [False, [144.0, 144.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', fit_region = \
            'ellipse[[239.37094084deg, -16.96415506deg], [1.1279arcsec, 0.7875arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 82944.0],
            'freq_bin': [1e-10, 16762501225.396851],
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
            'fit': [False, [0.9200466881631709, 0.9746655722260728, \
                        0.7626550313652652]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442221593894],
            'fit_pix': [False, [144.00051463175717, 144.00004766689185]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 82944.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [1e-4, 1.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .image.tt1 report
        im1_stats_dict = self.image_stats(img+'.image.tt1', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_im1_stats = {'com_bmaj': [False, 0.875946879387],
            'com_bmin': [False, 0.673672378063],
            'com_pa': [False, 88.5368652344],
            'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 13.6926994324],
            'max_val_pos': [True, [213, 77, 0, 0]],
            'min_val': [False, -9.34282302856],
            'min_val_pos': [True, [190, 160, 0, 0]],
            'im_rms': [False, 2.9229102716],
            'im_sum': [False, 8106.40348176],
            'regn_sum': [False, -3094.95175185],
            'npts_real': [True, 82944]}

        report9 = th.check_dict_vals(exp_im1_stats, im1_stats_dict, '.image.tt1', epsilon=self.epsilon)

        # .residual.tt1 report
        resid1_stats_dict = self.image_stats(img+'.residual.tt1', \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid1_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 47329.0],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0125599103048],
            'max_val_pos': [True, [162, 91, 0, 0]],
            'min_val': [False, -0.0108835268766],
            'min_val_pos': [True, [150, 51, 0, 0]],
            'im_rms': [False, 0.00403524922747],
            'im_sum': [False, -4.25091994218],
            'regn_sum': [False, 24.7472725619],
            'npts_real': [True, 82944]}

        report10 = th.check_dict_vals(exp_resid1_stats, resid1_stats_dict, \
            '.residual.tt1', epsilon=self.epsilon)

        # .model.tt1 report
        model1_stats_dict = self.image_stats(img+'.model.tt1', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model1_stats = {'npts': [True, 82944],
            'npts_unmasked': [1e-4, 82944.0],
            'freq_bin': [1e-10, 16762501225.396851],
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

        report11 = th.check_dict_vals(exp_model1_stats, model1_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt.tt1 report
        sumwt1_stats_dict = self.image_stats(img+'.sumwt.tt1')

        exp_sumwt1_stats = {'npts': [True, 1],
            'npts_unmasked': [1e-4, 1.0],
            'freq_bin': [1e-10, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 192323.671875],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 192323.671875],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 192323.673842],
            'npts_real': [True, 1]}

        report12 = th.check_dict_vals(exp_sumwt1_stats, sumwt1_stats_dict, \
            '.sumwt.tt1', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11 + report12

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict['test_standard_mtmfs_eph']['self.parallel'] = self.parallel
        test_dict['test_standard_mtmfs_eph']['report'] = report
        test_dict['test_standard_mtmfs_eph']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image.tt0', range_list=[-1.05, 1.05])
        self.mom8_creator(img+'.residual.tt0', range_list=[-1.05, 1.05])
        test_dict['test_standard_mtmfs_eph']['images'].extend( \
            (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_mtmfs_eph
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cal(self):
        ''' Calibrator image - field J2258-2758, spw 22 '''

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
            parallel=False, verbose=True)

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
            parallel=False, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[344.52480945deg, -27.97253944deg], [12.1076arcsec, 6.2463arcsec], 90.00000000deg]')

        exp_im_stats = {'com_bmaj': [False, 9.626538276672363],
            'com_bmin': [False, 4.673975467681885],
            'com_pa': [False, -86.3871307373],
            'npts': [True, 8100],
            'npts_unmasked': [True, 5041.0],
            'freq_bin': [1e-10, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.40606927872],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, -0.0115236453712],
            'min_val_pos': [True, [50, 16, 0, 0]],
            'im_rms': [False, 0.2012872201556086],
            'im_sum': [False, 168.25755220036578],
            'regn_sum': [False, 167.478445882909],
            'npts_real': [True, 8100],
            'fit': [False, [2.40974849537, 9.60566422443952, 4.668418785483205]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 220.30076542192973],
            'fit_pix': [False, [45.000405766, 45.0014155577]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [45, 85, 0, 0]), \
                (img+'.image', False, [45, 86, 0, 0]), \
                (img+'.image', True, [5, 45, 0, 0]), \
                (img+'.image', False, [4, 45, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 8100],
            'freq_bin': [1e-10, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 400],
            'mask_regns': [True, 1],
            'npts_real': [True, 8100]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[344.52476094deg, -27.97251802deg], [34.7828arcsec, 34.7011arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 5041.0],
            'freq_bin': [1e-10, 125009872.91876221],
            'start': [True, 2.20301e+11],
            'end': [True, 2.20301e+11],
            'start_delta': [False, 2.20301e+11],
            'end_delta': [False, 2.20301e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [45, 45, 0, 0]],
            'min_val': [False, 0.200092822313],
            'im_rms': [False, 0.577170049921],
            'npts_0.2': [True, 5041],
            'npts_0.5': [True, 2409],
            'npts_real': [True, 8100],
            'fit': [False, [1.0468035426303963, 45.181424068122176, \
                       45.18134398951289]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 220.30076542192973],
            'fit_pix': [False, [45.000270927482546, 45.00030384048325]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[344.52484287deg, -27.97253611deg], [8.0802arcsec, 4.8086arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 8100.0],
            'freq_bin': [1e-10, 125009872.91876221],
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
            'fit': [False, [1.0640200932648511, 8.801094080240267, \
                        4.303338569406158]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 220.30076542192973],
            'fit_pix': [False, [44.99810399006913, 44.996587647973605]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[344.52480945deg, -27.97253944deg], [12.1076arcsec, 6.2463arcsec], 90.00000000deg]')

        exp_resid_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 5041.0],
            'freq_bin': [1e-10, 125009872.91876221],
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

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[344.52480945deg, -27.97253944deg], [12.1076arcsec, 6.2463arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 8100],
            'npts_unmasked': [True, 8100.0],
            'freq_bin': [1e-10, 125009872.91876221],
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
            'mask_non0': [True, 0],
            'npts_real': [True, 8100]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 125009872.91876221],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        failed = self.filter_report(report)

        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict['test_standard_cal']['self.parallel'] = self.parallel
        test_dict['test_standard_cal']['report'] = report
        test_dict['test_standard_cal']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.015, 2.5])
        self.mom8_creator(img+'.residual', range_list=[-0.015, 2.5])
        test_dict['test_standard_cal']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cal
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cal_eph(self):
        ''' Standard (single field) ephemeris calibrator imaging - central field of Venus (field 2), spw 25 & 45 '''

        file_name = 'standard_cal_eph.iter'
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
            fastnoise=False, savemodel='none', parallel=False, verbose=True)

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
            savemodel='none', parallel=False, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report (test_standard_cal_eph)
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_im_stats = {'com_bmaj': [False, 0.875946879387],
            'com_bmin': [False, 0.673672378063],
            'com_pa': [False, 88.5368652344],
            'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [False, 16762501225.396851],
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

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', self.epsilon)

        # .mask report (test_standard_cal_eph)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 82944],
            'freq_bin': [False, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 82944]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', self.epsilon)

        # .pb report (test_standard_cal_eph)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[239.37090838deg, -16.96415647deg], [17.8437arcsec, 17.4772arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [False, 16762501225.396851],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [144, 144, 0, 0]],
            'min_val': [False, 0.200061768293],
            'im_rms': [False, 0.569403001285],
            'npts_0.2': [True, 47329],
            'npts_0.5': [True, 22365],
            'npts_real': [True, 82944],
            'fit': [False, [1.0286609217550002, 22.907692916947163, \
                       22.90769291676479]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [False, 253.57442221593894],
            'fit_pix': [False, [144.0, 144.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', self.epsilon)

        # .psf report (test_standard_cal_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.37094084deg, -16.96415506deg], [1.1279arcsec, 0.7875arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 82944.0],
            'freq_bin': [False, 16762501225.396851],
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
            'fit': [False, [0.9200466881631709, 0.9746655722260728, \
                        0.7626550313652652]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [False, 253.57442221593894],
            'fit_pix': [False, [144.00051463175717, 144.00004766689185]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', self.epsilon)

        # .residual report (test_standard_cal_eph)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 47329.0],
            'freq_bin': [False, 16762501225.396851],
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

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', self.epsilon)

        # .model report (test_standard_cal_eph)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 82944],
            'npts_unmasked': [True, 82944.0],
            'freq_bin': [False, 16762501225.396851],
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

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', self.epsilon)

        # .sumwt report (test_standard_cal_eph)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [False, 16762501225.396851],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', self.epsilon)

        # report combination (test_standard_cal_eph)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict['test_standard_cal_eph']['self.parallel'] = self.parallel
        test_dict['test_standard_cal_eph']['report'] = report
        test_dict['test_standard_cal_eph']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-1.05, 1.05])
        self.mom8_creator(img+'.residual', range_list=[-1.05, 1.05])
        test_dict['test_standard_cal_eph']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cal_eph
###############################################
###############################################

class Test_mosaic(test_tclean_base):


    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_cube(self):
        ''' Mosaic cube imaging - field SMIDGE_NWCloud, spw 22 '''

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
            specmode='cube', nchan=508, start='220.2526743594GHz', \
            width='0.2441741MHz', outframe='LSRK', \
            perchanweightdensity=False, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
                start='220.2526743594GHz', width='0.2441741MHz', \
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
                parallel=True, verbose=True)

            # retrieve per-channel beam statistics
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(img+'.image')

            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526743594GHz', width='0.2441741MHz', \
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
                parallel=False, verbose=True)

        else:
            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'],scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526743594GHz', width='0.2441741MHz', \
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
                restoringbeam='common', parallel=False, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report (test_mosaic_cube)
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]', field_regions = \
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
            'npts_unmasked': [False, 3338068.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 71.2575168759],
            'npts_real': [True, 5925312],
            'rms_per_field': [False, [0.089055932035104174, 0.087593317854421537, 0.088754854977895078, 0.087673584891401882, 0.088859674046269987, 0.087793556001629414, 0.087858029812576927]],
            'profile': [False, 1.2198751],
            'fit': [False, [1.193476708202794, 8.742582360593987, \
                    8.04574047263242]],
            'fit_loc_chan': [True, 252],
            'fit_loc_freq': [1e-10, 220.31420623259388],
            'fit_pix': [False, [44.5001461142688, 38.26244952541851]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_cube)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 5925312],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'mask_pix': [False, 3929],
            'mask_regns': [True, 1],
            'npts_real': [True, 5925312]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        exp_pb_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3338068.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, 0.200100466609],
            'im_rms': [False, 0.615051728939],
            'npts_0.2': [False, [6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571]],
            'npts_0.5': [False, [3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594]],
            'npts_real': [True, 5925312],
            'fit': [False, [1.0810930483564398, 69.22608257076189, \
                    69.16658859812452]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [54.070053902647572, 54.013503538761256]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.1],
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
            # CAS-9386 update with Build 100 serial run
            #'im_sum': [False, 66.9960417559],
            'im_sum': [False, 63.68329861515883],
            'npts_real': [True, 5925312],
            'fit_0': [False, [1.1012023947984113, 7.881990108924573, \
                      5.249591565272608]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [53.9877555807218, 54.00478115211167]],
            'fit_1': [False, [1.1012785744622782, 7.879844016452006, \
                      5.242454867627971]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [53.987518610964671, 54.004592842286243]],
            'fit_2': [False, [1.1012668735944595, 7.881043827878078, \
                      5.23808759170775]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [53.987317744517256, 54.004239429673945]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        exp_resid_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3338068.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 19.2033108851],
            'npts_real': [True, 5925312]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 5.9057832174],
            'mask_non0': [True, 0],
            'npts_real': [True, 5925312]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.1],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube
        exp_wt_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.393758654594],
            'max_val_pos': [True, [54, 54, 0, 0]],
            # CAS-9386 update based on build 100 serial
            #'min_val': [False, 7.45326979086e-05],
            'min_val': [False, 7.177492079790682e-05],
            'im_rms': [False, 0.140904168376],
            'im_sum': [False, 506828.976527],
            'npts_0.2': [False, [6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571]],
            'npts_0.5': [False, [3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594]],
            'npts_real': [True, 5925312]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            exp_bmin_dict = {'*0': 6.098905563354492, '*1': 6.098905563354492, '*10': 6.098886489868164, '*100': 6.095863342285156, '*101': 6.095863342285156, '*102': 6.095863342285156, '*103': 6.095863342285156, '*104': 6.095863342285156, '*105': 6.095863342285156, '*106': 6.095863342285156, '*107': 6.095808982849121, '*108': 6.095808982849121, '*109': 6.095808982849121, '*11': 6.098886489868164, '*110': 6.095808982849121, '*111': 6.0957818031311035, '*112': 6.0957818031311035, '*113': 6.095760822296143, '*114': 6.095760822296143, '*115': 6.095760822296143, '*116': 6.095760822296143, '*117': 6.095751762390137, '*118': 6.095702648162842, '*119': 6.094822883605957, '*12': 6.09887170791626, '*120': 6.094706058502197, '*121': 6.094687461853027, '*122': 6.094687461853027, '*123': 6.094687461853027, '*124': 6.094687461853027, '*125': 6.094687461853027, '*126': 6.094687461853027, '*127': 6.09466552734375, '*128': 6.09466552734375, '*129': 6.094656944274902, '*13': 6.09887170791626, '*130': 6.094656944274902, '*131': 6.094616889953613, '*132': 6.094616889953613, '*133': 6.094616889953613, '*134': 6.094616889953613, '*135': 6.094616889953613, '*136': 6.094616889953613, '*137': 6.09455680847168, '*138': 6.09455680847168, '*139': 6.09455680847168, '*14': 6.096207141876221, '*140': 6.09455680847168, '*141': 6.094529628753662, '*142': 6.094529628753662, '*143': 6.094529628753662, '*144': 6.094529628753662, '*145': 6.094529628753662, '*146': 6.094518184661865, '*147': 6.094518184661865, '*148': 6.094518184661865, '*149': 6.094503402709961, '*15': 6.0961480140686035, '*150': 6.094503402709961, '*151': 6.094503402709961, '*152': 6.094494342803955, '*153': 6.094494342803955, '*154': 6.094484329223633, '*155': 6.094484329223633, '*156': 6.094484329223633, '*157': 6.094484329223633, '*158': 6.094449996948242, '*159': 6.094449996948242, '*16': 6.0961480140686035, '*160': 6.094449996948242, '*161': 6.094449996948242, '*162': 6.094449996948242, '*163': 6.094395637512207, '*164': 6.094395637512207, '*165': 6.094395637512207, '*166': 6.094395637512207, '*167': 6.094395637512207, '*168': 6.094388008117676, '*169': 6.094388008117676, '*17': 6.0961480140686035, '*170': 6.094380855560303, '*171': 6.094380855560303, '*172': 6.094343662261963, '*173': 6.094334602355957, '*174': 6.094306468963623, '*175': 6.09429931640625, '*176': 6.0942792892456055, '*177': 6.0942792892456055, '*178': 6.094268798828125, '*179': 6.094268798828125, '*18': 6.096136569976807, '*180': 6.094268798828125, '*181': 6.094228744506836, '*182': 6.094203948974609, '*183': 6.094203948974609, '*184': 6.094203948974609, '*185': 6.094203948974609, '*186': 6.094203948974609, '*187': 6.094203948974609, '*188': 6.094203948974609, '*189': 6.094205379486084, '*19': 6.096136569976807, '*190': 6.094205379486084, '*191': 6.094142436981201, '*192': 6.094142436981201, '*193': 6.094144344329834, '*194': 6.094126224517822, '*195': 6.094126224517822, '*196': 6.094126224517822, '*197': 6.094119071960449, '*198': 6.094119071960449, '*199': 6.094119071960449, '*2': 6.098905563354492, '*20': 6.096136569976807, '*200': 6.094120979309082, '*201': 6.094120979309082, '*202': 6.094120979309082, '*203': 6.094120979309082, '*204': 6.094120979309082, '*205': 6.094120979309082, '*206': 6.094120979309082, '*207': 6.094120979309082, '*208': 6.094120979309082, '*209': 6.094099044799805, '*21': 6.096136569976807, '*210': 6.094099044799805, '*211': 6.094067096710205, '*212': 6.09406042098999, '*213': 6.094041347503662, '*214': 6.094041347503662, '*215': 6.09397554397583, '*216': 6.09397554397583, '*217': 6.09397554397583, '*218': 6.09397554397583, '*219': 6.093908309936523, '*22': 6.0961174964904785, '*220': 6.093873977661133, '*221': 6.093847751617432, '*222': 6.093847751617432, '*223': 6.093847751617432, '*224': 6.093782424926758, '*225': 6.0906219482421875, '*226': 6.0906219482421875, '*227': 6.0906219482421875, '*228': 6.0906219482421875, '*229': 6.0906219482421875, '*23': 6.0961174964904785, '*230': 6.0906219482421875, '*231': 6.090587139129639, '*232': 6.090587139129639, '*233': 6.090587139129639, '*234': 6.090587139129639, '*235': 6.090561389923096, '*236': 6.090561389923096, '*237': 6.090516090393066, '*238': 6.090525150299072, '*239': 6.0905280113220215, '*24': 6.0961174964904785, '*240': 6.090687274932861, '*241': 6.090687274932861, '*242': 6.090687274932861, '*243': 6.090682029724121, '*244': 6.090682029724121, '*245': 6.090682029724121, '*246': 6.0906877517700195, '*247': 6.090632438659668, '*248': 6.090632438659668, '*249': 6.090632438659668, '*25': 6.0961174964904785, '*250': 6.090611457824707, '*251': 6.090611457824707, '*252': 6.090573787689209, '*253': 6.090564250946045, '*254': 6.090564250946045, '*255': 6.090564250946045, '*256': 6.090507984161377, '*257': 6.090507984161377, '*258': 6.090507984161377, '*259': 6.090509414672852, '*26': 6.0961174964904785, '*260': 6.090509414672852, '*261': 6.090509414672852, '*262': 6.090509414672852, '*263': 6.0904340744018555, '*264': 6.0904340744018555, '*265': 6.0904340744018555, '*266': 6.0904340744018555, '*267': 6.0904340744018555, '*268': 6.0904340744018555, '*269': 6.090435028076172, '*27': 6.096085071563721, '*270': 6.090435028076172, '*271': 6.090393543243408, '*272': 6.090388774871826, '*273': 6.090388774871826, '*274': 6.090388774871826, '*275': 6.090388774871826, '*276': 6.090373516082764, '*277': 6.090373516082764, '*278': 6.090373516082764, '*279': 6.090373516082764, '*28': 6.096085071563721, '*280': 6.090358734130859, '*281': 6.088774681091309, '*282': 6.0886993408203125, '*283': 6.088681221008301, '*284': 6.088681221008301, '*285': 6.088681221008301, '*286': 6.088681221008301, '*287': 6.088675022125244, '*288': 6.088675022125244, '*289': 6.088625907897949, '*29': 6.096073627471924, '*290': 6.088708400726318, '*291': 6.088708400726318, '*292': 6.088721752166748, '*293': 6.088721752166748, '*294': 6.088721752166748, '*295': 6.088707447052002, '*296': 6.088707447052002, '*297': 6.088707447052002, '*298': 6.088682174682617, '*299': 6.088682174682617, '*3': 6.098886489868164, '*30': 6.096073627471924, '*300': 6.088682174682617, '*301': 6.088682174682617, '*302': 6.088682174682617, '*303': 6.088648796081543, '*304': 6.088648796081543, '*305': 6.088648796081543, '*306': 6.088648796081543, '*307': 6.088648796081543, '*308': 6.088629722595215, '*309': 6.088621616363525, '*31': 6.09605598449707, '*310': 6.088621616363525, '*311': 6.088621616363525, '*312': 6.088621616363525, '*313': 6.088621616363525, '*314': 6.088621616363525, '*315': 6.088593482971191, '*316': 6.088555812835693, '*317': 6.088555812835693, '*318': 6.0885443687438965, '*319': 6.0885443687438965, '*32': 6.09605598449707, '*320': 6.0885443687438965, '*321': 6.0885443687438965, '*322': 6.0885443687438965, '*323': 6.0885443687438965, '*324': 6.0885443687438965, '*325': 6.0885443687438965, '*326': 6.0885443687438965, '*327': 6.0885443687438965, '*328': 6.0885443687438965, '*329': 6.0885443687438965, '*33': 6.09605598449707, '*330': 6.0885443687438965, '*331': 6.0885443687438965, '*332': 6.0885443687438965, '*333': 6.0885443687438965, '*334': 6.0885443687438965, '*335': 6.0885443687438965, '*336': 6.088516712188721, '*337': 6.089442729949951, '*338': 6.089435577392578, '*339': 6.089423656463623, '*34': 6.09605598449707, '*340': 6.089423656463623, '*341': 6.089423656463623, '*342': 6.089423656463623, '*343': 6.089423656463623, '*344': 6.089423656463623, '*345': 6.08940315246582, '*346': 6.08937931060791, '*347': 6.089338302612305, '*348': 6.089310646057129, '*349': 6.089306831359863, '*35': 6.096081733703613, '*350': 6.088985919952393, '*351': 6.088985919952393, '*352': 6.088985919952393, '*353': 6.088985919952393, '*354': 6.088986873626709, '*355': 6.088986873626709, '*356': 6.088986873626709, '*357': 6.0888848304748535, '*358': 6.0888848304748535, '*359': 6.088881969451904, '*36': 6.096024990081787, '*360': 6.088850021362305, '*361': 6.088778018951416, '*362': 6.088778018951416, '*363': 6.08873176574707, '*364': 6.08873176574707, '*365': 6.08873176574707, '*366': 6.08873176574707, '*367': 6.0887131690979, '*368': 6.0887131690979, '*369': 6.0887131690979, '*37': 6.096024990081787, '*370': 6.0887131690979, '*371': 6.088683128356934, '*372': 6.088683128356934, '*373': 6.0886335372924805, '*374': 6.088621139526367, '*375': 6.088621139526367, '*376': 6.088621139526367, '*377': 6.088621139526367, '*378': 6.088607311248779, '*379': 6.088596343994141, '*38': 6.096024990081787, '*380': 6.088578224182129, '*381': 6.088578224182129, '*382': 6.088578224182129, '*383': 6.088578224182129, '*384': 6.08854341506958, '*385': 6.08854341506958, '*386': 6.088532447814941, '*387': 6.088512897491455, '*388': 6.0884833335876465, '*389': 6.088479518890381, '*39': 6.096009731292725, '*390': 6.088479518890381, '*391': 6.088479518890381, '*392': 6.088479518890381, '*393': 6.088479518890381, '*394': 6.088479518890381, '*395': 6.0884809494018555, '*396': 6.0884809494018555, '*397': 6.088460922241211, '*398': 6.088454246520996, '*399': 6.088454246520996, '*4': 6.098886489868164, '*40': 6.096009731292725, '*400': 6.088454246520996, '*401': 6.088454246520996, '*402': 6.088454246520996, '*403': 6.088400840759277, '*404': 6.088386535644531, '*405': 6.088386535644531, '*406': 6.088373184204102, '*407': 6.088373184204102, '*408': 6.088373184204102, '*409': 6.088373184204102, '*41': 6.096009731292725, '*410': 6.088373184204102, '*411': 6.088373184204102, '*412': 6.088373184204102, '*413': 6.0883469581604, '*414': 6.0883469581604, '*415': 6.088336944580078, '*416': 6.0886454582214355, '*417': 6.0886454582214355, '*418': 6.0886454582214355, '*419': 6.0886454582214355, '*42': 6.095993995666504, '*420': 6.0886454582214355, '*421': 6.0886454582214355, '*422': 6.0886454582214355, '*423': 6.0886454582214355, '*424': 6.0886454582214355, '*425': 6.0886454582214355, '*426': 6.08853006362915, '*427': 6.0885114669799805, '*428': 6.088533878326416, '*429': 6.088533878326416, '*43': 6.095993995666504, '*430': 6.088533878326416, '*431': 6.088514804840088, '*432': 6.088514804840088, '*433': 6.088474273681641, '*434': 6.088474273681641, '*435': 6.088474273681641, '*436': 6.088474273681641, '*437': 6.088474273681641, '*438': 6.088474273681641, '*439': 6.088474273681641, '*44': 6.095993995666504, '*440': 6.088443756103516, '*441': 6.088443756103516, '*442': 6.088443756103516, '*443': 6.088443756103516, '*444': 6.088418483734131, '*445': 6.088418483734131, '*446': 6.08836030960083, '*447': 6.08573579788208, '*448': 6.085716247558594, '*449': 6.085716247558594, '*45': 6.095993995666504, '*450': 6.085716247558594, '*451': 6.085702896118164, '*452': 6.085702896118164, '*453': 6.085702896118164, '*454': 6.085669994354248, '*455': 6.085644245147705, '*456': 6.085601806640625, '*457': 6.085601806640625, '*458': 6.085601806640625, '*459': 6.085562229156494, '*46': 6.096076011657715, '*460': 6.085562229156494, '*461': 6.08554220199585, '*462': 6.08554220199585, '*463': 6.08554220199585, '*464': 6.08554220199585, '*465': 6.08554220199585, '*466': 6.08554220199585, '*467': 6.08554220199585, '*468': 6.08554220199585, '*469': 6.08554220199585, '*47': 6.096029758453369, '*470': 6.08554220199585, '*471': 6.08554220199585, '*472': 6.08554220199585, '*473': 6.08554220199585, '*474': 6.08554220199585, '*475': 6.08554220199585, '*476': 6.08554220199585, '*477': 6.0854997634887695, '*478': 6.0854997634887695, '*479': 6.085493087768555, '*48': 6.096001625061035, '*480': 6.085493087768555, '*481': 6.085493087768555, '*482': 6.085493087768555, '*483': 6.085474967956543, '*484': 6.085474967956543, '*485': 6.085445404052734, '*486': 6.085445404052734, '*487': 6.085445404052734, '*488': 6.085445404052734, '*489': 6.085446834564209, '*49': 6.09599494934082, '*490': 6.085446834564209, '*491': 6.085446834564209, '*492': 6.085446834564209, '*493': 6.085446834564209, '*494': 6.085422515869141, '*495': 6.085422515869141, '*496': 6.085402011871338, '*497': 6.085402011871338, '*498': 6.085402011871338, '*499': 6.085396766662598, '*5': 6.098886489868164, '*50': 6.09599494934082, '*500': 6.085396766662598, '*501': 6.085396766662598, '*502': 6.085396766662598, '*503': 6.085396766662598, '*504': 6.085379123687744, '*505': 6.085379123687744, '*506': 6.085379123687744, '*507': 6.085379123687744, '*51': 6.09599494934082, '*52': 6.09599494934082, '*53': 6.09599494934082, '*54': 6.09599494934082, '*55': 6.09599494934082, '*56': 6.09599494934082, '*57': 6.09599494934082, '*58': 6.095983505249023, '*59': 6.095983505249023, '*6': 6.098886489868164, '*60': 6.095983505249023, '*61': 6.095983505249023, '*62': 6.095969200134277, '*63': 6.095969200134277, '*64': 6.095969200134277, '*65': 6.095958709716797, '*66': 6.095958709716797, '*67': 6.095958709716797, '*68': 6.095958709716797, '*69': 6.095958709716797, '*7': 6.098886489868164, '*70': 6.095940113067627, '*71': 6.0959086418151855, '*72': 6.0959086418151855, '*73': 6.095920562744141, '*74': 6.095920562744141, '*75': 6.095920562744141, '*76': 6.095920562744141, '*77': 6.095907688140869, '*78': 6.09588623046875, '*79': 6.09588623046875, '*8': 6.098886489868164, '*80': 6.09588623046875, '*81': 6.09588623046875, '*82': 6.09588623046875, '*83': 6.09588623046875, '*84': 6.09588623046875, '*85': 6.09588623046875, '*86': 6.095888137817383, '*87': 6.095888137817383, '*88': 6.095888137817383, '*89': 6.095888137817383, '*9': 6.098886489868164, '*90': 6.095885753631592, '*91': 6.095885753631592, '*92': 6.095885753631592, '*93': 6.095885753631592, '*94': 6.095874309539795, '*95': 6.095874309539795, '*96': 6.095874309539795, '*97': 6.095863342285156, '*98': 6.095863342285156, '*99': 6.095863342285156}

            exp_bmaj_dict = {'*0': 8.788710594177246, '*1': 8.788710594177246, '*10': 8.788692474365234, '*100': 8.787450790405273, '*101': 8.787450790405273, '*102': 8.787450790405273, '*103': 8.787450790405273, '*104': 8.787450790405273, '*105': 8.787450790405273, '*106': 8.787450790405273, '*107': 8.78741455078125, '*108': 8.78741455078125, '*109': 8.78741455078125, '*11': 8.788692474365234, '*110': 8.78741455078125, '*111': 8.787422180175781, '*112': 8.787422180175781, '*113': 8.787412643432617, '*114': 8.787412643432617, '*115': 8.787412643432617, '*116': 8.787412643432617, '*117': 8.787357330322266, '*118': 8.787363052368164, '*119': 8.786246299743652, '*12': 8.788710594177246, '*120': 8.786258697509766, '*121': 8.786214828491211, '*122': 8.786214828491211, '*123': 8.786214828491211, '*124': 8.786214828491211, '*125': 8.786214828491211, '*126': 8.786214828491211, '*127': 8.786246299743652, '*128': 8.786246299743652, '*129': 8.78619384765625, '*13': 8.788710594177246, '*130': 8.78619384765625, '*131': 8.786169052124023, '*132': 8.786169052124023, '*133': 8.786169052124023, '*134': 8.786169052124023, '*135': 8.786169052124023, '*136': 8.786169052124023, '*137': 8.786203384399414, '*138': 8.786203384399414, '*139': 8.786203384399414, '*14': 8.788393020629883, '*140': 8.786203384399414, '*141': 8.7861909866333, '*142': 8.7861909866333, '*143': 8.7861909866333, '*144': 8.7861909866333, '*145': 8.7861909866333, '*146': 8.78610897064209, '*147': 8.78610897064209, '*148': 8.78610897064209, '*149': 8.786099433898926, '*15': 8.788331985473633, '*150': 8.786099433898926, '*151': 8.786099433898926, '*152': 8.786046981811523, '*153': 8.786046981811523, '*154': 8.786056518554688, '*155': 8.786056518554688, '*156': 8.786056518554688, '*157': 8.786056518554688, '*158': 8.786031723022461, '*159': 8.786031723022461, '*16': 8.788331985473633, '*160': 8.786031723022461, '*161': 8.786031723022461, '*162': 8.786031723022461, '*163': 8.78590202331543, '*164': 8.78590202331543, '*165': 8.78590202331543, '*166': 8.78590202331543, '*167': 8.78590202331543, '*168': 8.785841941833496, '*169': 8.785841941833496, '*17': 8.788331985473633, '*170': 8.78583812713623, '*171': 8.78583812713623, '*172': 8.785764694213867, '*173': 8.785714149475098, '*174': 8.78570556640625, '*175': 8.78567886352539, '*176': 8.785659790039062, '*177': 8.785659790039062, '*178': 8.785603523254395, '*179': 8.785603523254395, '*18': 8.788325309753418, '*180': 8.785603523254395, '*181': 8.785640716552734, '*182': 8.78563117980957, '*183': 8.78563117980957, '*184': 8.78563117980957, '*185': 8.78563117980957, '*186': 8.78563117980957, '*187': 8.78563117980957, '*188': 8.78563117980957, '*189': 8.785579681396484, '*19': 8.788325309753418, '*190': 8.785579681396484, '*191': 8.785578727722168, '*192': 8.785578727722168, '*193': 8.78551959991455, '*194': 8.785514831542969, '*195': 8.785514831542969, '*196': 8.785514831542969, '*197': 8.785470962524414, '*198': 8.785470962524414, '*199': 8.785470962524414, '*2': 8.788710594177246, '*20': 8.788325309753418, '*200': 8.785435676574707, '*201': 8.785435676574707, '*202': 8.785435676574707, '*203': 8.785435676574707, '*204': 8.785435676574707, '*205': 8.785435676574707, '*206': 8.785435676574707, '*207': 8.785435676574707, '*208': 8.785435676574707, '*209': 8.785418510437012, '*21': 8.788325309753418, '*210': 8.785418510437012, '*211': 8.785431861877441, '*212': 8.785388946533203, '*213': 8.785396575927734, '*214': 8.785396575927734, '*215': 8.785368919372559, '*216': 8.785368919372559, '*217': 8.785368919372559, '*218': 8.785368919372559, '*219': 8.785303115844727, '*22': 8.788300514221191, '*220': 8.785208702087402, '*221': 8.78520679473877, '*222': 8.78520679473877, '*223': 8.78520679473877, '*224': 8.785249710083008, '*225': 8.787298202514648, '*226': 8.787298202514648, '*227': 8.787298202514648, '*228': 8.787298202514648, '*229': 8.787298202514648, '*23': 8.788300514221191, '*230': 8.787298202514648, '*231': 8.787242889404297, '*232': 8.787242889404297, '*233': 8.787242889404297, '*234': 8.787242889404297, '*235': 8.787264823913574, '*236': 8.787264823913574, '*237': 8.787288665771484, '*238': 8.787248611450195, '*239': 8.78724193572998, '*24': 8.788300514221191, '*240': 8.786666870117188, '*241': 8.786666870117188, '*242': 8.786666870117188, '*243': 8.786641120910645, '*244': 8.786641120910645, '*245': 8.786641120910645, '*246': 8.786598205566406, '*247': 8.78667163848877, '*248': 8.78667163848877, '*249': 8.78667163848877, '*25': 8.788300514221191, '*250': 8.786662101745605, '*251': 8.786662101745605, '*252': 8.786628723144531, '*253': 8.786534309387207, '*254': 8.786534309387207, '*255': 8.786534309387207, '*256': 8.78652572631836, '*257': 8.78652572631836, '*258': 8.78652572631836, '*259': 8.78647232055664, '*26': 8.788300514221191, '*260': 8.78647232055664, '*261': 8.78647232055664, '*262': 8.78647232055664, '*263': 8.786478996276855, '*264': 8.786478996276855, '*265': 8.786478996276855, '*266': 8.786478996276855, '*267': 8.786478996276855, '*268': 8.786478996276855, '*269': 8.78643798828125, '*27': 8.788321495056152, '*270': 8.78643798828125, '*271': 8.786358833312988, '*272': 8.786345481872559, '*273': 8.786345481872559, '*274': 8.786345481872559, '*275': 8.786345481872559, '*276': 8.78632926940918, '*277': 8.78632926940918, '*278': 8.78632926940918, '*279': 8.78632926940918, '*28': 8.788321495056152, '*280': 8.786319732666016, '*281': 8.787357330322266, '*282': 8.78735637664795, '*283': 8.787260055541992, '*284': 8.787260055541992, '*285': 8.787260055541992, '*286': 8.787260055541992, '*287': 8.78719425201416, '*288': 8.78719425201416, '*289': 8.787174224853516, '*29': 8.788315773010254, '*290': 8.78727912902832, '*291': 8.78727912902832, '*292': 8.787213325500488, '*293': 8.787213325500488, '*294': 8.787213325500488, '*295': 8.787132263183594, '*296': 8.787132263183594, '*297': 8.787132263183594, '*298': 8.787094116210938, '*299': 8.787094116210938, '*3': 8.788692474365234, '*30': 8.788315773010254, '*300': 8.787094116210938, '*301': 8.787094116210938, '*302': 8.787094116210938, '*303': 8.787084579467773, '*304': 8.787084579467773, '*305': 8.787084579467773, '*306': 8.787084579467773, '*307': 8.787084579467773, '*308': 8.787069320678711, '*309': 8.787017822265625, '*31': 8.788280487060547, '*310': 8.787017822265625, '*311': 8.787017822265625, '*312': 8.787017822265625, '*313': 8.787017822265625, '*314': 8.787017822265625, '*315': 8.787042617797852, '*316': 8.786896705627441, '*317': 8.786896705627441, '*318': 8.78686237335205, '*319': 8.78686237335205, '*32': 8.788280487060547, '*320': 8.78686237335205, '*321': 8.78686237335205, '*322': 8.78686237335205, '*323': 8.78686237335205, '*324': 8.78686237335205, '*325': 8.78686237335205, '*326': 8.78686237335205, '*327': 8.78686237335205, '*328': 8.78686237335205, '*329': 8.78686237335205, '*33': 8.788280487060547, '*330': 8.78686237335205, '*331': 8.78686237335205, '*332': 8.78686237335205, '*333': 8.78686237335205, '*334': 8.78686237335205, '*335': 8.78686237335205, '*336': 8.786849021911621, '*337': 8.787076950073242, '*338': 8.787121772766113, '*339': 8.78708267211914, '*34': 8.788280487060547, '*340': 8.78708267211914, '*341': 8.78708267211914, '*342': 8.78708267211914, '*343': 8.78708267211914, '*344': 8.78708267211914, '*345': 8.787066459655762, '*346': 8.787047386169434, '*347': 8.78708267211914, '*348': 8.787079811096191, '*349': 8.78708553314209, '*35': 8.788140296936035, '*350': 8.787437438964844, '*351': 8.787437438964844, '*352': 8.787437438964844, '*353': 8.787437438964844, '*354': 8.787362098693848, '*355': 8.787362098693848, '*356': 8.787362098693848, '*357': 8.787331581115723, '*358': 8.787331581115723, '*359': 8.787214279174805, '*36': 8.788116455078125, '*360': 8.787230491638184, '*361': 8.787257194519043, '*362': 8.787257194519043, '*363': 8.787291526794434, '*364': 8.787291526794434, '*365': 8.787291526794434, '*366': 8.787291526794434, '*367': 8.787285804748535, '*368': 8.787285804748535, '*369': 8.787285804748535, '*37': 8.788116455078125, '*370': 8.787285804748535, '*371': 8.787185668945312, '*372': 8.787185668945312, '*373': 8.787177085876465, '*374': 8.787156105041504, '*375': 8.787156105041504, '*376': 8.787156105041504, '*377': 8.787156105041504, '*378': 8.787126541137695, '*379': 8.787050247192383, '*38': 8.788116455078125, '*380': 8.787004470825195, '*381': 8.787004470825195, '*382': 8.787004470825195, '*383': 8.787004470825195, '*384': 8.787025451660156, '*385': 8.787025451660156, '*386': 8.786970138549805, '*387': 8.78696060180664, '*388': 8.786956787109375, '*389': 8.786932945251465, '*39': 8.78808879852295, '*390': 8.786932945251465, '*391': 8.786932945251465, '*392': 8.786932945251465, '*393': 8.786932945251465, '*394': 8.786932945251465, '*395': 8.786917686462402, '*396': 8.786917686462402, '*397': 8.786901473999023, '*398': 8.786880493164062, '*399': 8.786880493164062, '*4': 8.788692474365234, '*40': 8.78808879852295, '*400': 8.786880493164062, '*401': 8.786880493164062, '*402': 8.786880493164062, '*403': 8.786869049072266, '*404': 8.78687858581543, '*405': 8.78687858581543, '*406': 8.786901473999023, '*407': 8.786901473999023, '*408': 8.786901473999023, '*409': 8.786901473999023, '*41': 8.78808879852295, '*410': 8.786901473999023, '*411': 8.786901473999023, '*412': 8.786901473999023, '*413': 8.786819458007812, '*414': 8.786819458007812, '*415': 8.786763191223145, '*416': 8.7879638671875, '*417': 8.7879638671875, '*418': 8.7879638671875, '*419': 8.7879638671875, '*42': 8.788002967834473, '*420': 8.7879638671875, '*421': 8.7879638671875, '*422': 8.7879638671875, '*423': 8.7879638671875, '*424': 8.7879638671875, '*425': 8.7879638671875, '*426': 8.788007736206055, '*427': 8.787994384765625, '*428': 8.787877082824707, '*429': 8.787877082824707, '*43': 8.788002967834473, '*430': 8.787877082824707, '*431': 8.787899017333984, '*432': 8.787899017333984, '*433': 8.787874221801758, '*434': 8.787874221801758, '*435': 8.787874221801758, '*436': 8.787874221801758, '*437': 8.787874221801758, '*438': 8.787874221801758, '*439': 8.787874221801758, '*44': 8.788002967834473, '*440': 8.787858009338379, '*441': 8.787858009338379, '*442': 8.787858009338379, '*443': 8.787858009338379, '*444': 8.787843704223633, '*445': 8.787843704223633, '*446': 8.787836074829102, '*447': 8.787530899047852, '*448': 8.787516593933105, '*449': 8.787516593933105, '*45': 8.788002967834473, '*450': 8.787516593933105, '*451': 8.7875337600708, '*452': 8.7875337600708, '*453': 8.7875337600708, '*454': 8.787500381469727, '*455': 8.787508964538574, '*456': 8.787493705749512, '*457': 8.787493705749512, '*458': 8.787493705749512, '*459': 8.787456512451172, '*46': 8.788106918334961, '*460': 8.787456512451172, '*461': 8.787471771240234, '*462': 8.787471771240234, '*463': 8.787471771240234, '*464': 8.787471771240234, '*465': 8.787471771240234, '*466': 8.787471771240234, '*467': 8.787471771240234, '*468': 8.787471771240234, '*469': 8.787471771240234, '*47': 8.788034439086914, '*470': 8.787471771240234, '*471': 8.787471771240234, '*472': 8.787471771240234, '*473': 8.787471771240234, '*474': 8.787471771240234, '*475': 8.787471771240234, '*476': 8.787471771240234, '*477': 8.7874116897583, '*478': 8.7874116897583, '*479': 8.787361145019531, '*48': 8.788006782531738, '*480': 8.787361145019531, '*481': 8.787361145019531, '*482': 8.787361145019531, '*483': 8.787348747253418, '*484': 8.787348747253418, '*485': 8.78724193572998, '*486': 8.78724193572998, '*487': 8.78724193572998, '*488': 8.78724193572998, '*489': 8.787201881408691, '*49': 8.787985801696777, '*490': 8.787201881408691, '*491': 8.787201881408691, '*492': 8.787201881408691, '*493': 8.787201881408691, '*494': 8.787139892578125, '*495': 8.787139892578125, '*496': 8.787134170532227, '*497': 8.787134170532227, '*498': 8.787134170532227, '*499': 8.787117004394531, '*5': 8.788692474365234, '*50': 8.787985801696777, '*500': 8.787117004394531, '*501': 8.787117004394531, '*502': 8.787117004394531, '*503': 8.787117004394531, '*504': 8.78710651397705, '*505': 8.78710651397705, '*506': 8.78710651397705, '*507': 8.78710651397705, '*51': 8.787985801696777, '*52': 8.787985801696777, '*53': 8.787985801696777, '*54': 8.787985801696777, '*55': 8.787985801696777, '*56': 8.787985801696777, '*57': 8.787985801696777, '*58': 8.787952423095703, '*59': 8.787952423095703, '*6': 8.788692474365234, '*60': 8.787952423095703, '*61': 8.787952423095703, '*62': 8.78796672821045, '*63': 8.78796672821045, '*64': 8.78796672821045, '*65': 8.78791332244873, '*66': 8.78791332244873, '*67': 8.78791332244873, '*68': 8.78791332244873, '*69': 8.78791332244873, '*7': 8.788692474365234, '*70': 8.787813186645508, '*71': 8.787795066833496, '*72': 8.787795066833496, '*73': 8.787667274475098, '*74': 8.787667274475098, '*75': 8.787667274475098, '*76': 8.787667274475098, '*77': 8.787690162658691, '*78': 8.7876558303833, '*79': 8.7876558303833, '*8': 8.788692474365234, '*80': 8.7876558303833, '*81': 8.7876558303833, '*82': 8.7876558303833, '*83': 8.7876558303833, '*84': 8.7876558303833, '*85': 8.7876558303833, '*86': 8.787596702575684, '*87': 8.787596702575684, '*88': 8.787596702575684, '*89': 8.787596702575684, '*9': 8.788692474365234, '*90': 8.787479400634766, '*91': 8.787479400634766, '*92': 8.787479400634766, '*93': 8.787479400634766, '*94': 8.787473678588867, '*95': 8.787473678588867, '*96': 8.787473678588867, '*97': 8.787450790405273, '*98': 8.787450790405273, '*99': 8.787450790405273}

            exp_pa_dict = {'*0': 64.9308090209961, '*1': 64.9308090209961, '*10': 64.9310531616211, '*100': 64.86714172363281, '*101': 64.86714172363281, '*102': 64.86714172363281, '*103': 64.86714172363281, '*104': 64.86714172363281, '*105': 64.86714172363281, '*106': 64.86714172363281, '*107': 64.86788940429688, '*108': 64.86788940429688, '*109': 64.86788940429688, '*11': 64.9310531616211, '*110': 64.86788940429688, '*111': 64.86804962158203, '*112': 64.86804962158203, '*113': 64.86830139160156, '*114': 64.86830139160156, '*115': 64.86830139160156, '*116': 64.86830139160156, '*117': 64.86802673339844, '*118': 64.86830139160156, '*119': 64.87887573242188, '*12': 64.93042755126953, '*120': 64.87953186035156, '*121': 64.88017272949219, '*122': 64.88017272949219, '*123': 64.88017272949219, '*124': 64.88017272949219, '*125': 64.88017272949219, '*126': 64.88017272949219, '*127': 64.8794174194336, '*128': 64.8794174194336, '*129': 64.87901306152344, '*13': 64.93042755126953, '*130': 64.87901306152344, '*131': 64.87952423095703, '*132': 64.87952423095703, '*133': 64.87952423095703, '*134': 64.87952423095703, '*135': 64.87952423095703, '*136': 64.87952423095703, '*137': 64.87982177734375, '*138': 64.87982177734375, '*139': 64.87982177734375, '*14': 64.8731689453125, '*140': 64.87982177734375, '*141': 64.880126953125, '*142': 64.880126953125, '*143': 64.880126953125, '*144': 64.880126953125, '*145': 64.880126953125, '*146': 64.87919616699219, '*147': 64.87919616699219, '*148': 64.87919616699219, '*149': 64.87939453125, '*15': 64.87181091308594, '*150': 64.87939453125, '*151': 64.87939453125, '*152': 64.8790283203125, '*153': 64.8790283203125, '*154': 64.8790283203125, '*155': 64.8790283203125, '*156': 64.8790283203125, '*157': 64.8790283203125, '*158': 64.87954711914062, '*159': 64.87954711914062, '*16': 64.87181091308594, '*160': 64.87954711914062, '*161': 64.87954711914062, '*162': 64.87954711914062, '*163': 64.87789916992188, '*164': 64.87789916992188, '*165': 64.87789916992188, '*166': 64.87789916992188, '*167': 64.87789916992188, '*168': 64.87753295898438, '*169': 64.87753295898438, '*17': 64.87181091308594, '*170': 64.87736511230469, '*171': 64.87736511230469, '*172': 64.87860107421875, '*173': 64.87818908691406, '*174': 64.87852478027344, '*175': 64.87816619873047, '*176': 64.8784408569336, '*177': 64.8784408569336, '*178': 64.87816619873047, '*179': 64.87816619873047, '*18': 64.87200164794922, '*180': 64.87816619873047, '*181': 64.87799835205078, '*182': 64.87834167480469, '*183': 64.87834167480469, '*184': 64.87834167480469, '*185': 64.87834167480469, '*186': 64.87834167480469, '*187': 64.87834167480469, '*188': 64.87834167480469, '*189': 64.87824249267578, '*19': 64.87200164794922, '*190': 64.87824249267578, '*191': 64.87889099121094, '*192': 64.87889099121094, '*193': 64.8786849975586, '*194': 64.87889099121094, '*195': 64.87889099121094, '*196': 64.87889099121094, '*197': 64.87833404541016, '*198': 64.87833404541016, '*199': 64.87833404541016, '*2': 64.9308090209961, '*20': 64.87200164794922, '*200': 64.87796020507812, '*201': 64.87796020507812, '*202': 64.87796020507812, '*203': 64.87796020507812, '*204': 64.87796020507812, '*205': 64.87796020507812, '*206': 64.87796020507812, '*207': 64.87796020507812, '*208': 64.87796020507812, '*209': 64.87831115722656, '*21': 64.87200164794922, '*210': 64.87831115722656, '*211': 64.87844848632812, '*212': 64.87785339355469, '*213': 64.8779525756836, '*214': 64.8779525756836, '*215': 64.87743377685547, '*216': 64.87743377685547, '*217': 64.87743377685547, '*218': 64.87743377685547, '*219': 64.87794494628906, '*22': 64.87059783935547, '*220': 64.8759536743164, '*221': 64.87623596191406, '*222': 64.87623596191406, '*223': 64.87623596191406, '*224': 64.87592315673828, '*225': 64.861083984375, '*226': 64.861083984375, '*227': 64.861083984375, '*228': 64.861083984375, '*229': 64.861083984375, '*23': 64.87059783935547, '*230': 64.861083984375, '*231': 64.86063385009766, '*232': 64.86063385009766, '*233': 64.86063385009766, '*234': 64.86063385009766, '*235': 64.860595703125, '*236': 64.860595703125, '*237': 64.86064147949219, '*238': 64.86036682128906, '*239': 64.86065673828125, '*24': 64.87059783935547, '*240': 64.8814697265625, '*241': 64.8814697265625, '*242': 64.8814697265625, '*243': 64.88191223144531, '*244': 64.88191223144531, '*245': 64.88191223144531, '*246': 64.88255310058594, '*247': 64.88206481933594, '*248': 64.88206481933594, '*249': 64.88206481933594, '*25': 64.87059783935547, '*250': 64.88232421875, '*251': 64.88232421875, '*252': 64.88300323486328, '*253': 64.88233184814453, '*254': 64.88233184814453, '*255': 64.88233184814453, '*256': 64.88296508789062, '*257': 64.88296508789062, '*258': 64.88296508789062, '*259': 64.88285064697266, '*26': 64.87059783935547, '*260': 64.88285064697266, '*261': 64.88285064697266, '*262': 64.88285064697266, '*263': 64.88339233398438, '*264': 64.88339233398438, '*265': 64.88339233398438, '*266': 64.88339233398438, '*267': 64.88339233398438, '*268': 64.88339233398438, '*269': 64.88331604003906, '*27': 64.87067413330078, '*270': 64.88331604003906, '*271': 64.88362884521484, '*272': 64.88343048095703, '*273': 64.88343048095703, '*274': 64.88343048095703, '*275': 64.88343048095703, '*276': 64.88369750976562, '*277': 64.88369750976562, '*278': 64.88369750976562, '*279': 64.88369750976562, '*28': 64.87067413330078, '*280': 64.88389587402344, '*281': 64.88187408447266, '*282': 64.88233184814453, '*283': 64.88140869140625, '*284': 64.88140869140625, '*285': 64.88140869140625, '*286': 64.88140869140625, '*287': 64.88103485107422, '*288': 64.88103485107422, '*289': 64.88016510009766, '*29': 64.870849609375, '*290': 64.87916564941406, '*291': 64.87916564941406, '*292': 64.87899017333984, '*293': 64.87899017333984, '*294': 64.87899017333984, '*295': 64.87834930419922, '*296': 64.87834930419922, '*297': 64.87834930419922, '*298': 64.87855529785156, '*299': 64.87855529785156, '*3': 64.9310531616211, '*30': 64.870849609375, '*300': 64.87855529785156, '*301': 64.87855529785156, '*302': 64.87855529785156, '*303': 64.87895202636719, '*304': 64.87895202636719, '*305': 64.87895202636719, '*306': 64.87895202636719, '*307': 64.87895202636719, '*308': 64.87928771972656, '*309': 64.87892150878906, '*31': 64.86991882324219, '*310': 64.87892150878906, '*311': 64.87892150878906, '*312': 64.87892150878906, '*313': 64.87892150878906, '*314': 64.87892150878906, '*315': 64.87873077392578, '*316': 64.87691497802734, '*317': 64.87691497802734, '*318': 64.87628173828125, '*319': 64.87628173828125, '*32': 64.86991882324219, '*320': 64.87628173828125, '*321': 64.87628173828125, '*322': 64.87628173828125, '*323': 64.87628173828125, '*324': 64.87628173828125, '*325': 64.87628173828125, '*326': 64.87628173828125, '*327': 64.87628173828125, '*328': 64.87628173828125, '*329': 64.87628173828125, '*33': 64.86991882324219, '*330': 64.87628173828125, '*331': 64.87628173828125, '*332': 64.87628173828125, '*333': 64.87628173828125, '*334': 64.87628173828125, '*335': 64.87628173828125, '*336': 64.87659454345703, '*337': 64.90101623535156, '*338': 64.90149688720703, '*339': 64.9006118774414, '*34': 64.86991882324219, '*340': 64.9006118774414, '*341': 64.9006118774414, '*342': 64.9006118774414, '*343': 64.9006118774414, '*344': 64.9006118774414, '*345': 64.90101623535156, '*346': 64.90128326416016, '*347': 64.9015884399414, '*348': 64.90193939208984, '*349': 64.90177154541016, '*35': 64.86958312988281, '*350': 64.88825988769531, '*351': 64.88825988769531, '*352': 64.88825988769531, '*353': 64.88825988769531, '*354': 64.88855743408203, '*355': 64.88855743408203, '*356': 64.88855743408203, '*357': 64.88982391357422, '*358': 64.88982391357422, '*359': 64.88941955566406, '*36': 64.87027740478516, '*360': 64.88955688476562, '*361': 64.89034271240234, '*362': 64.89034271240234, '*363': 64.89037322998047, '*364': 64.89037322998047, '*365': 64.89037322998047, '*366': 64.89037322998047, '*367': 64.89057159423828, '*368': 64.89057159423828, '*369': 64.89057159423828, '*37': 64.87027740478516, '*370': 64.89057159423828, '*371': 64.8902816772461, '*372': 64.8902816772461, '*373': 64.89068603515625, '*374': 64.89106750488281, '*375': 64.89106750488281, '*376': 64.89106750488281, '*377': 64.89106750488281, '*378': 64.8904800415039, '*379': 64.88965606689453, '*38': 64.87027740478516, '*380': 64.8896484375, '*381': 64.8896484375, '*382': 64.8896484375, '*383': 64.8896484375, '*384': 64.8897933959961, '*385': 64.8897933959961, '*386': 64.88924407958984, '*387': 64.8895034790039, '*388': 64.88985443115234, '*389': 64.89020538330078, '*39': 64.8707275390625, '*390': 64.89020538330078, '*391': 64.89020538330078, '*392': 64.89020538330078, '*393': 64.89020538330078, '*394': 64.89020538330078, '*395': 64.89041900634766, '*396': 64.89041900634766, '*397': 64.89082336425781, '*398': 64.89054107666016, '*399': 64.89054107666016, '*4': 64.9310531616211, '*40': 64.8707275390625, '*400': 64.89054107666016, '*401': 64.89054107666016, '*402': 64.89054107666016, '*403': 64.89093780517578, '*404': 64.89095306396484, '*405': 64.89095306396484, '*406': 64.89080047607422, '*407': 64.89080047607422, '*408': 64.89080047607422, '*409': 64.89080047607422, '*41': 64.8707275390625, '*410': 64.89080047607422, '*411': 64.89080047607422, '*412': 64.89080047607422, '*413': 64.8907470703125, '*414': 64.8907470703125, '*415': 64.89019012451172, '*416': 64.88300323486328, '*417': 64.88300323486328, '*418': 64.88300323486328, '*419': 64.88300323486328, '*42': 64.87004089355469, '*420': 64.88300323486328, '*421': 64.88300323486328, '*422': 64.88300323486328, '*423': 64.88300323486328, '*424': 64.88300323486328, '*425': 64.88300323486328, '*426': 64.8835678100586, '*427': 64.88374328613281, '*428': 64.88340759277344, '*429': 64.88340759277344, '*43': 64.87004089355469, '*430': 64.88340759277344, '*431': 64.88275909423828, '*432': 64.88275909423828, '*433': 64.88326263427734, '*434': 64.88326263427734, '*435': 64.88326263427734, '*436': 64.88326263427734, '*437': 64.88326263427734, '*438': 64.88326263427734, '*439': 64.88326263427734, '*44': 64.87004089355469, '*440': 64.8836669921875, '*441': 64.8836669921875, '*442': 64.8836669921875, '*443': 64.8836669921875, '*444': 64.88400268554688, '*445': 64.88400268554688, '*446': 64.88275146484375, '*447': 64.82654571533203, '*448': 64.82695007324219, '*449': 64.82695007324219, '*45': 64.87004089355469, '*450': 64.82695007324219, '*451': 64.82649230957031, '*452': 64.82649230957031, '*453': 64.82649230957031, '*454': 64.8271255493164, '*455': 64.82719421386719, '*456': 64.8277359008789, '*457': 64.8277359008789, '*458': 64.8277359008789, '*459': 64.82827758789062, '*46': 64.8690414428711, '*460': 64.82827758789062, '*461': 64.82827758789062, '*462': 64.82827758789062, '*463': 64.82827758789062, '*464': 64.82827758789062, '*465': 64.82827758789062, '*466': 64.82827758789062, '*467': 64.82827758789062, '*468': 64.82827758789062, '*469': 64.82827758789062, '*47': 64.86918640136719, '*470': 64.82827758789062, '*471': 64.82827758789062, '*472': 64.82827758789062, '*473': 64.82827758789062, '*474': 64.82827758789062, '*475': 64.82827758789062, '*476': 64.82827758789062, '*477': 64.82826232910156, '*478': 64.82826232910156, '*479': 64.82791137695312, '*48': 64.8697280883789, '*480': 64.82791137695312, '*481': 64.82791137695312, '*482': 64.82791137695312, '*483': 64.82807922363281, '*484': 64.82807922363281, '*485': 64.8299560546875, '*486': 64.8299560546875, '*487': 64.8299560546875, '*488': 64.8299560546875, '*489': 64.82986450195312, '*49': 64.86937713623047, '*490': 64.82986450195312, '*491': 64.82986450195312, '*492': 64.82986450195312, '*493': 64.82986450195312, '*494': 64.83010864257812, '*495': 64.83010864257812, '*496': 64.8304214477539, '*497': 64.8304214477539, '*498': 64.8304214477539, '*499': 64.83016967773438, '*5': 64.9310531616211, '*50': 64.86937713623047, '*500': 64.83016967773438, '*501': 64.83016967773438, '*502': 64.83016967773438, '*503': 64.83016967773438, '*504': 64.83045959472656, '*505': 64.83045959472656, '*506': 64.83045959472656, '*507': 64.83045959472656, '*51': 64.86937713623047, '*52': 64.86937713623047, '*53': 64.86937713623047, '*54': 64.86937713623047, '*55': 64.86937713623047, '*56': 64.86937713623047, '*57': 64.86937713623047, '*58': 64.86885070800781, '*59': 64.86885070800781, '*6': 64.9310531616211, '*60': 64.86885070800781, '*61': 64.86885070800781, '*62': 64.8682861328125, '*63': 64.8682861328125, '*64': 64.8682861328125, '*65': 64.86774444580078, '*66': 64.86774444580078, '*67': 64.86774444580078, '*68': 64.86774444580078, '*69': 64.86774444580078, '*7': 64.9310531616211, '*70': 64.86693572998047, '*71': 64.86734771728516, '*72': 64.86734771728516, '*73': 64.86736297607422, '*74': 64.86736297607422, '*75': 64.86736297607422, '*76': 64.86736297607422, '*77': 64.86717224121094, '*78': 64.86811065673828, '*79': 64.86811065673828, '*8': 64.9310531616211, '*80': 64.86811065673828, '*81': 64.86811065673828, '*82': 64.86811065673828, '*83': 64.86811065673828, '*84': 64.86811065673828, '*85': 64.86811065673828, '*86': 64.86791229248047, '*87': 64.86791229248047, '*88': 64.86791229248047, '*89': 64.86791229248047, '*9': 64.9310531616211, '*90': 64.86750793457031, '*91': 64.86750793457031, '*92': 64.86750793457031, '*93': 64.86750793457031, '*94': 64.86768341064453, '*95': 64.86768341064453, '*96': 64.86768341064453, '*97': 64.86714172363281, '*98': 64.86714172363281, '*99': 64.86714172363281}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed = self.filter_report(report)

        img = shutil._basename(img)
        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, 'test_mosaic_cube', self.parallel)

        test_dict['test_mosaic_cube']['report'] = report
        test_dict['test_mosaic_cube']['images'] = []

        self.mom8_creator(img+'.image', range_list=[0.15, 1.2])
        self.mom8_creator(img+'.residual', range_list=[0.15, 1.2])
        test_dict['test_mosaic_cube']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_mosaic_cube']['images'].append(img+'.image.profile.png')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_cube
#-------------------------------------------------#
    @stats_dict(test_dict)
    def test_mosaic_cube_briggsbwtaper(self):
        ''' Mosaic cube imaging wth briggsbwtaper - field SMIDGE_NWCloud, spw 22 '''

        test_name = self._testMethodName
        file_name = 'mosaic_cube_briggsbwtaper.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
            antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[108, 108], cell=['1.1arcsec'], \
            phasecenter='ICRS 00:45:54.3836 -073.15.29.413', stokes='I', \
            specmode='cube', nchan=508, start='220.2526743594GHz', \
            width='0.2441741MHz', outframe='LSRK', \
            perchanweightdensity=True, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, minbeamfrac=0.1, \
            growiterations=75, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
                start='220.2526743594GHz', width='0.2441741MHz', \
                outframe='LSRK', perchanweightdensity=True, \
                gridder='mosaic', chanchunks=-1, mosweight=True, \
                usepointing=False, pblimit=0.2, deconvolver='hogbom', \
                restoration=True, pbcor=True, weighting='briggsbwtaper', robust=0.5,\
                npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
                interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                savemodel='none', calcres=False, calcpsf=False, \
                parallel=True, verbose=True)

            # retrieve per-channel beam statistics
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(img+'.image')

            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'], scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526743594GHz', width='0.2441741MHz', \
                outframe='LSRK', perchanweightdensity=True, \
                gridder='mosaic', chanchunks=-1, mosweight=True, \
                usepointing=False, pblimit=0.2, deconvolver='hogbom', \
                restoration=True, restoringbeam='common', pbcor=True, \
                weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, \
                threshold='0.354Jy', nsigma=0.0, interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                savemodel='none', calcres=False, calcpsf=False, 
                parallel=False, verbose=True)

        else:
            tclean(vis=self.msfile, field='SMIDGE_NWCloud', spw=['0'], \
                antenna=['0,1,2,3,4,5,6,7,8'],scan=['8,12,16'], \
                intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
                imagename=file_name+'1', imsize=[108, 108], \
                cell=['1.1arcsec'], phasecenter='ICRS 00:45:54.3836'
                ' -073.15.29.413', stokes='I', specmode='cube', nchan=508, \
                start='220.2526743594GHz', width='0.2441741MHz', \
                outframe='LSRK', perchanweightdensity=True, \
                gridder='mosaic', chanchunks=-1, mosweight=True, \
                usepointing=False, pblimit=0.2, deconvolver='hogbom', \
                restoration=True, pbcor=True, weighting='briggsbwtaper', robust=0.5,\
                npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
                interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                savemodel='none', calcres=False, calcpsf=False, \
                restoringbeam='common', parallel=False, verbose=True)

            # retrieve per-channel beam statistics
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(img+'.psf')

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report (test_mosaic_cube_briggsbwtaper)
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]', field_regions = \
            ['circle[[00:45:54.383559, -73.15.29.41306], 22.45arcsec]',
             'circle[[00:45:49.435664, -73.15.35.13742], 22.45arcsec]',
             'circle[[00:45:53.057440, -73.15.50.79016], 22.45arcsec]',
             'circle[[00:45:50.762696, -73.15.13.76177], 22.45arcsec]',
             'circle[[00:45:58.006248, -73.15.45.06040], 22.45arcsec]',
             'circle[[00:45:55.708764, -73.15.08.03543], 22.45arcsec]',
             'circle[[00:45:59.330540, -73.15.23.68133], 22.45arcsec]'])

        # test_mosaic_cube_briggsbwtaper.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 8.79758467135],
            'com_bmin': [False, 6.10500627626],
            'com_pa': [False, 64.9303836873],
            'npts': [True, 5925312],
            'npts_unmasked': [False, 3338068.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 71.2575168759],
            'npts_real': [True, 5925312],
            'rms_per_field': [False, [0.089055932035104174, 0.087593317854421537, 0.088754854977895078, 0.087673584891401882, 0.088859674046269987, 0.087793556001629414, 0.087858029812576927]],
            'profile': [False, 1.2198751],
            'fit': [False, [1.193476708202794, 8.742582360593987, \
                    8.04574047263242]],
            'fit_loc_chan': [True, 252],
            'fit_loc_freq': [1e-10, 220.31420623259388],
            'fit_pix': [False, [44.5001461142688, 38.26244952541851]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_cube_briggsbwtaper)
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_mosaic_cube_briggsbwtaper.exp_mask_stats
        exp_mask_stats = {'npts': [True, 5925312],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'mask_pix': [False, 3929],
            'mask_regns': [True, 1],
            'npts_real': [True, 5925312]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube_briggsbwtaper)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3338068.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, 0.200100466609],
            'im_rms': [False, 0.615051728939],
            'npts_0.2': [False, [6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571]],
            'npts_0.5': [False, [3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594]],
            'npts_real': [True, 5925312],
            'fit': [False, [1.0810930483564398, 69.22608257076189, \
                    69.16658859812452]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [54.070053902647572, 54.013503538761256]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_briggsbwtaper)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.1],
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
            # CAS-9386 update with Build 100 serial run
            #'im_sum': [False, 66.9960417559],
            'im_sum': [False, 63.68329861515883],
            'npts_real': [True, 5925312],
            'fit_0': [False, [1.1012023947984113, 7.881990108924573, \
                      5.249591565272608]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [53.9877555807218, 54.00478115211167]],
            'fit_1': [False, [1.1012785744622782, 7.879844016452006, \
                      5.242454867627971]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [53.987518610964671, 54.004592842286243]],
            'fit_2': [False, [1.1012668735944595, 7.881043827878078, \
                      5.23808759170775]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [53.987317744517256, 54.004239429673945]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3338068.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 19.2033108851],
            'npts_real': [True, 5925312]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_briggsbwtaper)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.1],
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
            'regn_sum': [False, 5.9057832174],
            'mask_non0': [True, 0],
            'npts_real': [True, 5925312]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.1],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube_briggsbwtaper.exp_wt_stats
        exp_wt_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.1],
            'start': [True, 2.202527e+11],
            'end': [True, 2.203765e+11],
            'start_delta': [False, 2.202527e+11],
            'end_delta': [False, 2.203765e+11],
            'nchan': [True, 508],
            'max_val': [False, 0.393758654594],
            'max_val_pos': [True, [54, 54, 0, 0]],
            # CAS-9386 update based on build 100 serial
            #'min_val': [False, 7.45326979086e-05],
            'min_val': [False, 7.177492079790682e-05],
            'im_rms': [False, 0.140904168376],
            'im_sum': [False, 506828.976527],
            'npts_0.2': [False, [6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571, 6571]],
            'npts_0.5': [False, [3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3593, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594]],
            'npts_real': [True, 5925312]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            # test_mosaic_cube_briggsbwtaper.exp_bmin_dict
            exp_bmin_dict = {'*0': 6.098905563354492, '*1': 6.098905563354492, '*10': 6.098886489868164, '*100': 6.095863342285156, '*101': 6.095863342285156, '*102': 6.095863342285156, '*103': 6.095863342285156, '*104': 6.095863342285156, '*105': 6.095863342285156, '*106': 6.095863342285156, '*107': 6.095808982849121, '*108': 6.095808982849121, '*109': 6.095808982849121, '*11': 6.098886489868164, '*110': 6.095808982849121, '*111': 6.0957818031311035, '*112': 6.0957818031311035, '*113': 6.095760822296143, '*114': 6.095760822296143, '*115': 6.095760822296143, '*116': 6.095760822296143, '*117': 6.095751762390137, '*118': 6.095702648162842, '*119': 6.094822883605957, '*12': 6.09887170791626, '*120': 6.094706058502197, '*121': 6.094687461853027, '*122': 6.094687461853027, '*123': 6.094687461853027, '*124': 6.094687461853027, '*125': 6.094687461853027, '*126': 6.094687461853027, '*127': 6.09466552734375, '*128': 6.09466552734375, '*129': 6.094656944274902, '*13': 6.09887170791626, '*130': 6.094656944274902, '*131': 6.094616889953613, '*132': 6.094616889953613, '*133': 6.094616889953613, '*134': 6.094616889953613, '*135': 6.094616889953613, '*136': 6.094616889953613, '*137': 6.09455680847168, '*138': 6.09455680847168, '*139': 6.09455680847168, '*14': 6.096207141876221, '*140': 6.09455680847168, '*141': 6.094529628753662, '*142': 6.094529628753662, '*143': 6.094529628753662, '*144': 6.094529628753662, '*145': 6.094529628753662, '*146': 6.094518184661865, '*147': 6.094518184661865, '*148': 6.094518184661865, '*149': 6.094503402709961, '*15': 6.0961480140686035, '*150': 6.094503402709961, '*151': 6.094503402709961, '*152': 6.094494342803955, '*153': 6.094494342803955, '*154': 6.094484329223633, '*155': 6.094484329223633, '*156': 6.094484329223633, '*157': 6.094484329223633, '*158': 6.094449996948242, '*159': 6.094449996948242, '*16': 6.0961480140686035, '*160': 6.094449996948242, '*161': 6.094449996948242, '*162': 6.094449996948242, '*163': 6.094395637512207, '*164': 6.094395637512207, '*165': 6.094395637512207, '*166': 6.094395637512207, '*167': 6.094395637512207, '*168': 6.094388008117676, '*169': 6.094388008117676, '*17': 6.0961480140686035, '*170': 6.094380855560303, '*171': 6.094380855560303, '*172': 6.094343662261963, '*173': 6.094334602355957, '*174': 6.094306468963623, '*175': 6.09429931640625, '*176': 6.0942792892456055, '*177': 6.0942792892456055, '*178': 6.094268798828125, '*179': 6.094268798828125, '*18': 6.096136569976807, '*180': 6.094268798828125, '*181': 6.094228744506836, '*182': 6.094203948974609, '*183': 6.094203948974609, '*184': 6.094203948974609, '*185': 6.094203948974609, '*186': 6.094203948974609, '*187': 6.094203948974609, '*188': 6.094203948974609, '*189': 6.094205379486084, '*19': 6.096136569976807, '*190': 6.094205379486084, '*191': 6.094142436981201, '*192': 6.094142436981201, '*193': 6.094144344329834, '*194': 6.094126224517822, '*195': 6.094126224517822, '*196': 6.094126224517822, '*197': 6.094119071960449, '*198': 6.094119071960449, '*199': 6.094119071960449, '*2': 6.098905563354492, '*20': 6.096136569976807, '*200': 6.094120979309082, '*201': 6.094120979309082, '*202': 6.094120979309082, '*203': 6.094120979309082, '*204': 6.094120979309082, '*205': 6.094120979309082, '*206': 6.094120979309082, '*207': 6.094120979309082, '*208': 6.094120979309082, '*209': 6.094099044799805, '*21': 6.096136569976807, '*210': 6.094099044799805, '*211': 6.094067096710205, '*212': 6.09406042098999, '*213': 6.094041347503662, '*214': 6.094041347503662, '*215': 6.09397554397583, '*216': 6.09397554397583, '*217': 6.09397554397583, '*218': 6.09397554397583, '*219': 6.093908309936523, '*22': 6.0961174964904785, '*220': 6.093873977661133, '*221': 6.093847751617432, '*222': 6.093847751617432, '*223': 6.093847751617432, '*224': 6.093782424926758, '*225': 6.0906219482421875, '*226': 6.0906219482421875, '*227': 6.0906219482421875, '*228': 6.0906219482421875, '*229': 6.0906219482421875, '*23': 6.0961174964904785, '*230': 6.0906219482421875, '*231': 6.090587139129639, '*232': 6.090587139129639, '*233': 6.090587139129639, '*234': 6.090587139129639, '*235': 6.090561389923096, '*236': 6.090561389923096, '*237': 6.090516090393066, '*238': 6.090525150299072, '*239': 6.0905280113220215, '*24': 6.0961174964904785, '*240': 6.090687274932861, '*241': 6.090687274932861, '*242': 6.090687274932861, '*243': 6.090682029724121, '*244': 6.090682029724121, '*245': 6.090682029724121, '*246': 6.0906877517700195, '*247': 6.090632438659668, '*248': 6.090632438659668, '*249': 6.090632438659668, '*25': 6.0961174964904785, '*250': 6.090611457824707, '*251': 6.090611457824707, '*252': 6.090573787689209, '*253': 6.090564250946045, '*254': 6.090564250946045, '*255': 6.090564250946045, '*256': 6.090507984161377, '*257': 6.090507984161377, '*258': 6.090507984161377, '*259': 6.090509414672852, '*26': 6.0961174964904785, '*260': 6.090509414672852, '*261': 6.090509414672852, '*262': 6.090509414672852, '*263': 6.0904340744018555, '*264': 6.0904340744018555, '*265': 6.0904340744018555, '*266': 6.0904340744018555, '*267': 6.0904340744018555, '*268': 6.0904340744018555, '*269': 6.090435028076172, '*27': 6.096085071563721, '*270': 6.090435028076172, '*271': 6.090393543243408, '*272': 6.090388774871826, '*273': 6.090388774871826, '*274': 6.090388774871826, '*275': 6.090388774871826, '*276': 6.090373516082764, '*277': 6.090373516082764, '*278': 6.090373516082764, '*279': 6.090373516082764, '*28': 6.096085071563721, '*280': 6.090358734130859, '*281': 6.088774681091309, '*282': 6.0886993408203125, '*283': 6.088681221008301, '*284': 6.088681221008301, '*285': 6.088681221008301, '*286': 6.088681221008301, '*287': 6.088675022125244, '*288': 6.088675022125244, '*289': 6.088625907897949, '*29': 6.096073627471924, '*290': 6.088708400726318, '*291': 6.088708400726318, '*292': 6.088721752166748, '*293': 6.088721752166748, '*294': 6.088721752166748, '*295': 6.088707447052002, '*296': 6.088707447052002, '*297': 6.088707447052002, '*298': 6.088682174682617, '*299': 6.088682174682617, '*3': 6.098886489868164, '*30': 6.096073627471924, '*300': 6.088682174682617, '*301': 6.088682174682617, '*302': 6.088682174682617, '*303': 6.088648796081543, '*304': 6.088648796081543, '*305': 6.088648796081543, '*306': 6.088648796081543, '*307': 6.088648796081543, '*308': 6.088629722595215, '*309': 6.088621616363525, '*31': 6.09605598449707, '*310': 6.088621616363525, '*311': 6.088621616363525, '*312': 6.088621616363525, '*313': 6.088621616363525, '*314': 6.088621616363525, '*315': 6.088593482971191, '*316': 6.088555812835693, '*317': 6.088555812835693, '*318': 6.0885443687438965, '*319': 6.0885443687438965, '*32': 6.09605598449707, '*320': 6.0885443687438965, '*321': 6.0885443687438965, '*322': 6.0885443687438965, '*323': 6.0885443687438965, '*324': 6.0885443687438965, '*325': 6.0885443687438965, '*326': 6.0885443687438965, '*327': 6.0885443687438965, '*328': 6.0885443687438965, '*329': 6.0885443687438965, '*33': 6.09605598449707, '*330': 6.0885443687438965, '*331': 6.0885443687438965, '*332': 6.0885443687438965, '*333': 6.0885443687438965, '*334': 6.0885443687438965, '*335': 6.0885443687438965, '*336': 6.088516712188721, '*337': 6.089442729949951, '*338': 6.089435577392578, '*339': 6.089423656463623, '*34': 6.09605598449707, '*340': 6.089423656463623, '*341': 6.089423656463623, '*342': 6.089423656463623, '*343': 6.089423656463623, '*344': 6.089423656463623, '*345': 6.08940315246582, '*346': 6.08937931060791, '*347': 6.089338302612305, '*348': 6.089310646057129, '*349': 6.089306831359863, '*35': 6.096081733703613, '*350': 6.088985919952393, '*351': 6.088985919952393, '*352': 6.088985919952393, '*353': 6.088985919952393, '*354': 6.088986873626709, '*355': 6.088986873626709, '*356': 6.088986873626709, '*357': 6.0888848304748535, '*358': 6.0888848304748535, '*359': 6.088881969451904, '*36': 6.096024990081787, '*360': 6.088850021362305, '*361': 6.088778018951416, '*362': 6.088778018951416, '*363': 6.08873176574707, '*364': 6.08873176574707, '*365': 6.08873176574707, '*366': 6.08873176574707, '*367': 6.0887131690979, '*368': 6.0887131690979, '*369': 6.0887131690979, '*37': 6.096024990081787, '*370': 6.0887131690979, '*371': 6.088683128356934, '*372': 6.088683128356934, '*373': 6.0886335372924805, '*374': 6.088621139526367, '*375': 6.088621139526367, '*376': 6.088621139526367, '*377': 6.088621139526367, '*378': 6.088607311248779, '*379': 6.088596343994141, '*38': 6.096024990081787, '*380': 6.088578224182129, '*381': 6.088578224182129, '*382': 6.088578224182129, '*383': 6.088578224182129, '*384': 6.08854341506958, '*385': 6.08854341506958, '*386': 6.088532447814941, '*387': 6.088512897491455, '*388': 6.0884833335876465, '*389': 6.088479518890381, '*39': 6.096009731292725, '*390': 6.088479518890381, '*391': 6.088479518890381, '*392': 6.088479518890381, '*393': 6.088479518890381, '*394': 6.088479518890381, '*395': 6.0884809494018555, '*396': 6.0884809494018555, '*397': 6.088460922241211, '*398': 6.088454246520996, '*399': 6.088454246520996, '*4': 6.098886489868164, '*40': 6.096009731292725, '*400': 6.088454246520996, '*401': 6.088454246520996, '*402': 6.088454246520996, '*403': 6.088400840759277, '*404': 6.088386535644531, '*405': 6.088386535644531, '*406': 6.088373184204102, '*407': 6.088373184204102, '*408': 6.088373184204102, '*409': 6.088373184204102, '*41': 6.096009731292725, '*410': 6.088373184204102, '*411': 6.088373184204102, '*412': 6.088373184204102, '*413': 6.0883469581604, '*414': 6.0883469581604, '*415': 6.088336944580078, '*416': 6.0886454582214355, '*417': 6.0886454582214355, '*418': 6.0886454582214355, '*419': 6.0886454582214355, '*42': 6.095993995666504, '*420': 6.0886454582214355, '*421': 6.0886454582214355, '*422': 6.0886454582214355, '*423': 6.0886454582214355, '*424': 6.0886454582214355, '*425': 6.0886454582214355, '*426': 6.08853006362915, '*427': 6.0885114669799805, '*428': 6.088533878326416, '*429': 6.088533878326416, '*43': 6.095993995666504, '*430': 6.088533878326416, '*431': 6.088514804840088, '*432': 6.088514804840088, '*433': 6.088474273681641, '*434': 6.088474273681641, '*435': 6.088474273681641, '*436': 6.088474273681641, '*437': 6.088474273681641, '*438': 6.088474273681641, '*439': 6.088474273681641, '*44': 6.095993995666504, '*440': 6.088443756103516, '*441': 6.088443756103516, '*442': 6.088443756103516, '*443': 6.088443756103516, '*444': 6.088418483734131, '*445': 6.088418483734131, '*446': 6.08836030960083, '*447': 6.08573579788208, '*448': 6.085716247558594, '*449': 6.085716247558594, '*45': 6.095993995666504, '*450': 6.085716247558594, '*451': 6.085702896118164, '*452': 6.085702896118164, '*453': 6.085702896118164, '*454': 6.085669994354248, '*455': 6.085644245147705, '*456': 6.085601806640625, '*457': 6.085601806640625, '*458': 6.085601806640625, '*459': 6.085562229156494, '*46': 6.096076011657715, '*460': 6.085562229156494, '*461': 6.08554220199585, '*462': 6.08554220199585, '*463': 6.08554220199585, '*464': 6.08554220199585, '*465': 6.08554220199585, '*466': 6.08554220199585, '*467': 6.08554220199585, '*468': 6.08554220199585, '*469': 6.08554220199585, '*47': 6.096029758453369, '*470': 6.08554220199585, '*471': 6.08554220199585, '*472': 6.08554220199585, '*473': 6.08554220199585, '*474': 6.08554220199585, '*475': 6.08554220199585, '*476': 6.08554220199585, '*477': 6.0854997634887695, '*478': 6.0854997634887695, '*479': 6.085493087768555, '*48': 6.096001625061035, '*480': 6.085493087768555, '*481': 6.085493087768555, '*482': 6.085493087768555, '*483': 6.085474967956543, '*484': 6.085474967956543, '*485': 6.085445404052734, '*486': 6.085445404052734, '*487': 6.085445404052734, '*488': 6.085445404052734, '*489': 6.085446834564209, '*49': 6.09599494934082, '*490': 6.085446834564209, '*491': 6.085446834564209, '*492': 6.085446834564209, '*493': 6.085446834564209, '*494': 6.085422515869141, '*495': 6.085422515869141, '*496': 6.085402011871338, '*497': 6.085402011871338, '*498': 6.085402011871338, '*499': 6.085396766662598, '*5': 6.098886489868164, '*50': 6.09599494934082, '*500': 6.085396766662598, '*501': 6.085396766662598, '*502': 6.085396766662598, '*503': 6.085396766662598, '*504': 6.085379123687744, '*505': 6.085379123687744, '*506': 6.085379123687744, '*507': 6.085379123687744, '*51': 6.09599494934082, '*52': 6.09599494934082, '*53': 6.09599494934082, '*54': 6.09599494934082, '*55': 6.09599494934082, '*56': 6.09599494934082, '*57': 6.09599494934082, '*58': 6.095983505249023, '*59': 6.095983505249023, '*6': 6.098886489868164, '*60': 6.095983505249023, '*61': 6.095983505249023, '*62': 6.095969200134277, '*63': 6.095969200134277, '*64': 6.095969200134277, '*65': 6.095958709716797, '*66': 6.095958709716797, '*67': 6.095958709716797, '*68': 6.095958709716797, '*69': 6.095958709716797, '*7': 6.098886489868164, '*70': 6.095940113067627, '*71': 6.0959086418151855, '*72': 6.0959086418151855, '*73': 6.095920562744141, '*74': 6.095920562744141, '*75': 6.095920562744141, '*76': 6.095920562744141, '*77': 6.095907688140869, '*78': 6.09588623046875, '*79': 6.09588623046875, '*8': 6.098886489868164, '*80': 6.09588623046875, '*81': 6.09588623046875, '*82': 6.09588623046875, '*83': 6.09588623046875, '*84': 6.09588623046875, '*85': 6.09588623046875, '*86': 6.095888137817383, '*87': 6.095888137817383, '*88': 6.095888137817383, '*89': 6.095888137817383, '*9': 6.098886489868164, '*90': 6.095885753631592, '*91': 6.095885753631592, '*92': 6.095885753631592, '*93': 6.095885753631592, '*94': 6.095874309539795, '*95': 6.095874309539795, '*96': 6.095874309539795, '*97': 6.095863342285156, '*98': 6.095863342285156, '*99': 6.095863342285156}

            # test_mosaic_cube_briggsbwtaper.exp_bmaj_dict
            exp_bmaj_dict = {'*0': 8.788710594177246, '*1': 8.788710594177246, '*10': 8.788692474365234, '*100': 8.787450790405273, '*101': 8.787450790405273, '*102': 8.787450790405273, '*103': 8.787450790405273, '*104': 8.787450790405273, '*105': 8.787450790405273, '*106': 8.787450790405273, '*107': 8.78741455078125, '*108': 8.78741455078125, '*109': 8.78741455078125, '*11': 8.788692474365234, '*110': 8.78741455078125, '*111': 8.787422180175781, '*112': 8.787422180175781, '*113': 8.787412643432617, '*114': 8.787412643432617, '*115': 8.787412643432617, '*116': 8.787412643432617, '*117': 8.787357330322266, '*118': 8.787363052368164, '*119': 8.786246299743652, '*12': 8.788710594177246, '*120': 8.786258697509766, '*121': 8.786214828491211, '*122': 8.786214828491211, '*123': 8.786214828491211, '*124': 8.786214828491211, '*125': 8.786214828491211, '*126': 8.786214828491211, '*127': 8.786246299743652, '*128': 8.786246299743652, '*129': 8.78619384765625, '*13': 8.788710594177246, '*130': 8.78619384765625, '*131': 8.786169052124023, '*132': 8.786169052124023, '*133': 8.786169052124023, '*134': 8.786169052124023, '*135': 8.786169052124023, '*136': 8.786169052124023, '*137': 8.786203384399414, '*138': 8.786203384399414, '*139': 8.786203384399414, '*14': 8.788393020629883, '*140': 8.786203384399414, '*141': 8.7861909866333, '*142': 8.7861909866333, '*143': 8.7861909866333, '*144': 8.7861909866333, '*145': 8.7861909866333, '*146': 8.78610897064209, '*147': 8.78610897064209, '*148': 8.78610897064209, '*149': 8.786099433898926, '*15': 8.788331985473633, '*150': 8.786099433898926, '*151': 8.786099433898926, '*152': 8.786046981811523, '*153': 8.786046981811523, '*154': 8.786056518554688, '*155': 8.786056518554688, '*156': 8.786056518554688, '*157': 8.786056518554688, '*158': 8.786031723022461, '*159': 8.786031723022461, '*16': 8.788331985473633, '*160': 8.786031723022461, '*161': 8.786031723022461, '*162': 8.786031723022461, '*163': 8.78590202331543, '*164': 8.78590202331543, '*165': 8.78590202331543, '*166': 8.78590202331543, '*167': 8.78590202331543, '*168': 8.785841941833496, '*169': 8.785841941833496, '*17': 8.788331985473633, '*170': 8.78583812713623, '*171': 8.78583812713623, '*172': 8.785764694213867, '*173': 8.785714149475098, '*174': 8.78570556640625, '*175': 8.78567886352539, '*176': 8.785659790039062, '*177': 8.785659790039062, '*178': 8.785603523254395, '*179': 8.785603523254395, '*18': 8.788325309753418, '*180': 8.785603523254395, '*181': 8.785640716552734, '*182': 8.78563117980957, '*183': 8.78563117980957, '*184': 8.78563117980957, '*185': 8.78563117980957, '*186': 8.78563117980957, '*187': 8.78563117980957, '*188': 8.78563117980957, '*189': 8.785579681396484, '*19': 8.788325309753418, '*190': 8.785579681396484, '*191': 8.785578727722168, '*192': 8.785578727722168, '*193': 8.78551959991455, '*194': 8.785514831542969, '*195': 8.785514831542969, '*196': 8.785514831542969, '*197': 8.785470962524414, '*198': 8.785470962524414, '*199': 8.785470962524414, '*2': 8.788710594177246, '*20': 8.788325309753418, '*200': 8.785435676574707, '*201': 8.785435676574707, '*202': 8.785435676574707, '*203': 8.785435676574707, '*204': 8.785435676574707, '*205': 8.785435676574707, '*206': 8.785435676574707, '*207': 8.785435676574707, '*208': 8.785435676574707, '*209': 8.785418510437012, '*21': 8.788325309753418, '*210': 8.785418510437012, '*211': 8.785431861877441, '*212': 8.785388946533203, '*213': 8.785396575927734, '*214': 8.785396575927734, '*215': 8.785368919372559, '*216': 8.785368919372559, '*217': 8.785368919372559, '*218': 8.785368919372559, '*219': 8.785303115844727, '*22': 8.788300514221191, '*220': 8.785208702087402, '*221': 8.78520679473877, '*222': 8.78520679473877, '*223': 8.78520679473877, '*224': 8.785249710083008, '*225': 8.787298202514648, '*226': 8.787298202514648, '*227': 8.787298202514648, '*228': 8.787298202514648, '*229': 8.787298202514648, '*23': 8.788300514221191, '*230': 8.787298202514648, '*231': 8.787242889404297, '*232': 8.787242889404297, '*233': 8.787242889404297, '*234': 8.787242889404297, '*235': 8.787264823913574, '*236': 8.787264823913574, '*237': 8.787288665771484, '*238': 8.787248611450195, '*239': 8.78724193572998, '*24': 8.788300514221191, '*240': 8.786666870117188, '*241': 8.786666870117188, '*242': 8.786666870117188, '*243': 8.786641120910645, '*244': 8.786641120910645, '*245': 8.786641120910645, '*246': 8.786598205566406, '*247': 8.78667163848877, '*248': 8.78667163848877, '*249': 8.78667163848877, '*25': 8.788300514221191, '*250': 8.786662101745605, '*251': 8.786662101745605, '*252': 8.786628723144531, '*253': 8.786534309387207, '*254': 8.786534309387207, '*255': 8.786534309387207, '*256': 8.78652572631836, '*257': 8.78652572631836, '*258': 8.78652572631836, '*259': 8.78647232055664, '*26': 8.788300514221191, '*260': 8.78647232055664, '*261': 8.78647232055664, '*262': 8.78647232055664, '*263': 8.786478996276855, '*264': 8.786478996276855, '*265': 8.786478996276855, '*266': 8.786478996276855, '*267': 8.786478996276855, '*268': 8.786478996276855, '*269': 8.78643798828125, '*27': 8.788321495056152, '*270': 8.78643798828125, '*271': 8.786358833312988, '*272': 8.786345481872559, '*273': 8.786345481872559, '*274': 8.786345481872559, '*275': 8.786345481872559, '*276': 8.78632926940918, '*277': 8.78632926940918, '*278': 8.78632926940918, '*279': 8.78632926940918, '*28': 8.788321495056152, '*280': 8.786319732666016, '*281': 8.787357330322266, '*282': 8.78735637664795, '*283': 8.787260055541992, '*284': 8.787260055541992, '*285': 8.787260055541992, '*286': 8.787260055541992, '*287': 8.78719425201416, '*288': 8.78719425201416, '*289': 8.787174224853516, '*29': 8.788315773010254, '*290': 8.78727912902832, '*291': 8.78727912902832, '*292': 8.787213325500488, '*293': 8.787213325500488, '*294': 8.787213325500488, '*295': 8.787132263183594, '*296': 8.787132263183594, '*297': 8.787132263183594, '*298': 8.787094116210938, '*299': 8.787094116210938, '*3': 8.788692474365234, '*30': 8.788315773010254, '*300': 8.787094116210938, '*301': 8.787094116210938, '*302': 8.787094116210938, '*303': 8.787084579467773, '*304': 8.787084579467773, '*305': 8.787084579467773, '*306': 8.787084579467773, '*307': 8.787084579467773, '*308': 8.787069320678711, '*309': 8.787017822265625, '*31': 8.788280487060547, '*310': 8.787017822265625, '*311': 8.787017822265625, '*312': 8.787017822265625, '*313': 8.787017822265625, '*314': 8.787017822265625, '*315': 8.787042617797852, '*316': 8.786896705627441, '*317': 8.786896705627441, '*318': 8.78686237335205, '*319': 8.78686237335205, '*32': 8.788280487060547, '*320': 8.78686237335205, '*321': 8.78686237335205, '*322': 8.78686237335205, '*323': 8.78686237335205, '*324': 8.78686237335205, '*325': 8.78686237335205, '*326': 8.78686237335205, '*327': 8.78686237335205, '*328': 8.78686237335205, '*329': 8.78686237335205, '*33': 8.788280487060547, '*330': 8.78686237335205, '*331': 8.78686237335205, '*332': 8.78686237335205, '*333': 8.78686237335205, '*334': 8.78686237335205, '*335': 8.78686237335205, '*336': 8.786849021911621, '*337': 8.787076950073242, '*338': 8.787121772766113, '*339': 8.78708267211914, '*34': 8.788280487060547, '*340': 8.78708267211914, '*341': 8.78708267211914, '*342': 8.78708267211914, '*343': 8.78708267211914, '*344': 8.78708267211914, '*345': 8.787066459655762, '*346': 8.787047386169434, '*347': 8.78708267211914, '*348': 8.787079811096191, '*349': 8.78708553314209, '*35': 8.788140296936035, '*350': 8.787437438964844, '*351': 8.787437438964844, '*352': 8.787437438964844, '*353': 8.787437438964844, '*354': 8.787362098693848, '*355': 8.787362098693848, '*356': 8.787362098693848, '*357': 8.787331581115723, '*358': 8.787331581115723, '*359': 8.787214279174805, '*36': 8.788116455078125, '*360': 8.787230491638184, '*361': 8.787257194519043, '*362': 8.787257194519043, '*363': 8.787291526794434, '*364': 8.787291526794434, '*365': 8.787291526794434, '*366': 8.787291526794434, '*367': 8.787285804748535, '*368': 8.787285804748535, '*369': 8.787285804748535, '*37': 8.788116455078125, '*370': 8.787285804748535, '*371': 8.787185668945312, '*372': 8.787185668945312, '*373': 8.787177085876465, '*374': 8.787156105041504, '*375': 8.787156105041504, '*376': 8.787156105041504, '*377': 8.787156105041504, '*378': 8.787126541137695, '*379': 8.787050247192383, '*38': 8.788116455078125, '*380': 8.787004470825195, '*381': 8.787004470825195, '*382': 8.787004470825195, '*383': 8.787004470825195, '*384': 8.787025451660156, '*385': 8.787025451660156, '*386': 8.786970138549805, '*387': 8.78696060180664, '*388': 8.786956787109375, '*389': 8.786932945251465, '*39': 8.78808879852295, '*390': 8.786932945251465, '*391': 8.786932945251465, '*392': 8.786932945251465, '*393': 8.786932945251465, '*394': 8.786932945251465, '*395': 8.786917686462402, '*396': 8.786917686462402, '*397': 8.786901473999023, '*398': 8.786880493164062, '*399': 8.786880493164062, '*4': 8.788692474365234, '*40': 8.78808879852295, '*400': 8.786880493164062, '*401': 8.786880493164062, '*402': 8.786880493164062, '*403': 8.786869049072266, '*404': 8.78687858581543, '*405': 8.78687858581543, '*406': 8.786901473999023, '*407': 8.786901473999023, '*408': 8.786901473999023, '*409': 8.786901473999023, '*41': 8.78808879852295, '*410': 8.786901473999023, '*411': 8.786901473999023, '*412': 8.786901473999023, '*413': 8.786819458007812, '*414': 8.786819458007812, '*415': 8.786763191223145, '*416': 8.7879638671875, '*417': 8.7879638671875, '*418': 8.7879638671875, '*419': 8.7879638671875, '*42': 8.788002967834473, '*420': 8.7879638671875, '*421': 8.7879638671875, '*422': 8.7879638671875, '*423': 8.7879638671875, '*424': 8.7879638671875, '*425': 8.7879638671875, '*426': 8.788007736206055, '*427': 8.787994384765625, '*428': 8.787877082824707, '*429': 8.787877082824707, '*43': 8.788002967834473, '*430': 8.787877082824707, '*431': 8.787899017333984, '*432': 8.787899017333984, '*433': 8.787874221801758, '*434': 8.787874221801758, '*435': 8.787874221801758, '*436': 8.787874221801758, '*437': 8.787874221801758, '*438': 8.787874221801758, '*439': 8.787874221801758, '*44': 8.788002967834473, '*440': 8.787858009338379, '*441': 8.787858009338379, '*442': 8.787858009338379, '*443': 8.787858009338379, '*444': 8.787843704223633, '*445': 8.787843704223633, '*446': 8.787836074829102, '*447': 8.787530899047852, '*448': 8.787516593933105, '*449': 8.787516593933105, '*45': 8.788002967834473, '*450': 8.787516593933105, '*451': 8.7875337600708, '*452': 8.7875337600708, '*453': 8.7875337600708, '*454': 8.787500381469727, '*455': 8.787508964538574, '*456': 8.787493705749512, '*457': 8.787493705749512, '*458': 8.787493705749512, '*459': 8.787456512451172, '*46': 8.788106918334961, '*460': 8.787456512451172, '*461': 8.787471771240234, '*462': 8.787471771240234, '*463': 8.787471771240234, '*464': 8.787471771240234, '*465': 8.787471771240234, '*466': 8.787471771240234, '*467': 8.787471771240234, '*468': 8.787471771240234, '*469': 8.787471771240234, '*47': 8.788034439086914, '*470': 8.787471771240234, '*471': 8.787471771240234, '*472': 8.787471771240234, '*473': 8.787471771240234, '*474': 8.787471771240234, '*475': 8.787471771240234, '*476': 8.787471771240234, '*477': 8.7874116897583, '*478': 8.7874116897583, '*479': 8.787361145019531, '*48': 8.788006782531738, '*480': 8.787361145019531, '*481': 8.787361145019531, '*482': 8.787361145019531, '*483': 8.787348747253418, '*484': 8.787348747253418, '*485': 8.78724193572998, '*486': 8.78724193572998, '*487': 8.78724193572998, '*488': 8.78724193572998, '*489': 8.787201881408691, '*49': 8.787985801696777, '*490': 8.787201881408691, '*491': 8.787201881408691, '*492': 8.787201881408691, '*493': 8.787201881408691, '*494': 8.787139892578125, '*495': 8.787139892578125, '*496': 8.787134170532227, '*497': 8.787134170532227, '*498': 8.787134170532227, '*499': 8.787117004394531, '*5': 8.788692474365234, '*50': 8.787985801696777, '*500': 8.787117004394531, '*501': 8.787117004394531, '*502': 8.787117004394531, '*503': 8.787117004394531, '*504': 8.78710651397705, '*505': 8.78710651397705, '*506': 8.78710651397705, '*507': 8.78710651397705, '*51': 8.787985801696777, '*52': 8.787985801696777, '*53': 8.787985801696777, '*54': 8.787985801696777, '*55': 8.787985801696777, '*56': 8.787985801696777, '*57': 8.787985801696777, '*58': 8.787952423095703, '*59': 8.787952423095703, '*6': 8.788692474365234, '*60': 8.787952423095703, '*61': 8.787952423095703, '*62': 8.78796672821045, '*63': 8.78796672821045, '*64': 8.78796672821045, '*65': 8.78791332244873, '*66': 8.78791332244873, '*67': 8.78791332244873, '*68': 8.78791332244873, '*69': 8.78791332244873, '*7': 8.788692474365234, '*70': 8.787813186645508, '*71': 8.787795066833496, '*72': 8.787795066833496, '*73': 8.787667274475098, '*74': 8.787667274475098, '*75': 8.787667274475098, '*76': 8.787667274475098, '*77': 8.787690162658691, '*78': 8.7876558303833, '*79': 8.7876558303833, '*8': 8.788692474365234, '*80': 8.7876558303833, '*81': 8.7876558303833, '*82': 8.7876558303833, '*83': 8.7876558303833, '*84': 8.7876558303833, '*85': 8.7876558303833, '*86': 8.787596702575684, '*87': 8.787596702575684, '*88': 8.787596702575684, '*89': 8.787596702575684, '*9': 8.788692474365234, '*90': 8.787479400634766, '*91': 8.787479400634766, '*92': 8.787479400634766, '*93': 8.787479400634766, '*94': 8.787473678588867, '*95': 8.787473678588867, '*96': 8.787473678588867, '*97': 8.787450790405273, '*98': 8.787450790405273, '*99': 8.787450790405273}

            # test_mosaic_cube_briggsbwtaper.exp_pa_dict
            exp_pa_dict = {'*0': 64.9308090209961, '*1': 64.9308090209961, '*10': 64.9310531616211, '*100': 64.86714172363281, '*101': 64.86714172363281, '*102': 64.86714172363281, '*103': 64.86714172363281, '*104': 64.86714172363281, '*105': 64.86714172363281, '*106': 64.86714172363281, '*107': 64.86788940429688, '*108': 64.86788940429688, '*109': 64.86788940429688, '*11': 64.9310531616211, '*110': 64.86788940429688, '*111': 64.86804962158203, '*112': 64.86804962158203, '*113': 64.86830139160156, '*114': 64.86830139160156, '*115': 64.86830139160156, '*116': 64.86830139160156, '*117': 64.86802673339844, '*118': 64.86830139160156, '*119': 64.87887573242188, '*12': 64.93042755126953, '*120': 64.87953186035156, '*121': 64.88017272949219, '*122': 64.88017272949219, '*123': 64.88017272949219, '*124': 64.88017272949219, '*125': 64.88017272949219, '*126': 64.88017272949219, '*127': 64.8794174194336, '*128': 64.8794174194336, '*129': 64.87901306152344, '*13': 64.93042755126953, '*130': 64.87901306152344, '*131': 64.87952423095703, '*132': 64.87952423095703, '*133': 64.87952423095703, '*134': 64.87952423095703, '*135': 64.87952423095703, '*136': 64.87952423095703, '*137': 64.87982177734375, '*138': 64.87982177734375, '*139': 64.87982177734375, '*14': 64.8731689453125, '*140': 64.87982177734375, '*141': 64.880126953125, '*142': 64.880126953125, '*143': 64.880126953125, '*144': 64.880126953125, '*145': 64.880126953125, '*146': 64.87919616699219, '*147': 64.87919616699219, '*148': 64.87919616699219, '*149': 64.87939453125, '*15': 64.87181091308594, '*150': 64.87939453125, '*151': 64.87939453125, '*152': 64.8790283203125, '*153': 64.8790283203125, '*154': 64.8790283203125, '*155': 64.8790283203125, '*156': 64.8790283203125, '*157': 64.8790283203125, '*158': 64.87954711914062, '*159': 64.87954711914062, '*16': 64.87181091308594, '*160': 64.87954711914062, '*161': 64.87954711914062, '*162': 64.87954711914062, '*163': 64.87789916992188, '*164': 64.87789916992188, '*165': 64.87789916992188, '*166': 64.87789916992188, '*167': 64.87789916992188, '*168': 64.87753295898438, '*169': 64.87753295898438, '*17': 64.87181091308594, '*170': 64.87736511230469, '*171': 64.87736511230469, '*172': 64.87860107421875, '*173': 64.87818908691406, '*174': 64.87852478027344, '*175': 64.87816619873047, '*176': 64.8784408569336, '*177': 64.8784408569336, '*178': 64.87816619873047, '*179': 64.87816619873047, '*18': 64.87200164794922, '*180': 64.87816619873047, '*181': 64.87799835205078, '*182': 64.87834167480469, '*183': 64.87834167480469, '*184': 64.87834167480469, '*185': 64.87834167480469, '*186': 64.87834167480469, '*187': 64.87834167480469, '*188': 64.87834167480469, '*189': 64.87824249267578, '*19': 64.87200164794922, '*190': 64.87824249267578, '*191': 64.87889099121094, '*192': 64.87889099121094, '*193': 64.8786849975586, '*194': 64.87889099121094, '*195': 64.87889099121094, '*196': 64.87889099121094, '*197': 64.87833404541016, '*198': 64.87833404541016, '*199': 64.87833404541016, '*2': 64.9308090209961, '*20': 64.87200164794922, '*200': 64.87796020507812, '*201': 64.87796020507812, '*202': 64.87796020507812, '*203': 64.87796020507812, '*204': 64.87796020507812, '*205': 64.87796020507812, '*206': 64.87796020507812, '*207': 64.87796020507812, '*208': 64.87796020507812, '*209': 64.87831115722656, '*21': 64.87200164794922, '*210': 64.87831115722656, '*211': 64.87844848632812, '*212': 64.87785339355469, '*213': 64.8779525756836, '*214': 64.8779525756836, '*215': 64.87743377685547, '*216': 64.87743377685547, '*217': 64.87743377685547, '*218': 64.87743377685547, '*219': 64.87794494628906, '*22': 64.87059783935547, '*220': 64.8759536743164, '*221': 64.87623596191406, '*222': 64.87623596191406, '*223': 64.87623596191406, '*224': 64.87592315673828, '*225': 64.861083984375, '*226': 64.861083984375, '*227': 64.861083984375, '*228': 64.861083984375, '*229': 64.861083984375, '*23': 64.87059783935547, '*230': 64.861083984375, '*231': 64.86063385009766, '*232': 64.86063385009766, '*233': 64.86063385009766, '*234': 64.86063385009766, '*235': 64.860595703125, '*236': 64.860595703125, '*237': 64.86064147949219, '*238': 64.86036682128906, '*239': 64.86065673828125, '*24': 64.87059783935547, '*240': 64.8814697265625, '*241': 64.8814697265625, '*242': 64.8814697265625, '*243': 64.88191223144531, '*244': 64.88191223144531, '*245': 64.88191223144531, '*246': 64.88255310058594, '*247': 64.88206481933594, '*248': 64.88206481933594, '*249': 64.88206481933594, '*25': 64.87059783935547, '*250': 64.88232421875, '*251': 64.88232421875, '*252': 64.88300323486328, '*253': 64.88233184814453, '*254': 64.88233184814453, '*255': 64.88233184814453, '*256': 64.88296508789062, '*257': 64.88296508789062, '*258': 64.88296508789062, '*259': 64.88285064697266, '*26': 64.87059783935547, '*260': 64.88285064697266, '*261': 64.88285064697266, '*262': 64.88285064697266, '*263': 64.88339233398438, '*264': 64.88339233398438, '*265': 64.88339233398438, '*266': 64.88339233398438, '*267': 64.88339233398438, '*268': 64.88339233398438, '*269': 64.88331604003906, '*27': 64.87067413330078, '*270': 64.88331604003906, '*271': 64.88362884521484, '*272': 64.88343048095703, '*273': 64.88343048095703, '*274': 64.88343048095703, '*275': 64.88343048095703, '*276': 64.88369750976562, '*277': 64.88369750976562, '*278': 64.88369750976562, '*279': 64.88369750976562, '*28': 64.87067413330078, '*280': 64.88389587402344, '*281': 64.88187408447266, '*282': 64.88233184814453, '*283': 64.88140869140625, '*284': 64.88140869140625, '*285': 64.88140869140625, '*286': 64.88140869140625, '*287': 64.88103485107422, '*288': 64.88103485107422, '*289': 64.88016510009766, '*29': 64.870849609375, '*290': 64.87916564941406, '*291': 64.87916564941406, '*292': 64.87899017333984, '*293': 64.87899017333984, '*294': 64.87899017333984, '*295': 64.87834930419922, '*296': 64.87834930419922, '*297': 64.87834930419922, '*298': 64.87855529785156, '*299': 64.87855529785156, '*3': 64.9310531616211, '*30': 64.870849609375, '*300': 64.87855529785156, '*301': 64.87855529785156, '*302': 64.87855529785156, '*303': 64.87895202636719, '*304': 64.87895202636719, '*305': 64.87895202636719, '*306': 64.87895202636719, '*307': 64.87895202636719, '*308': 64.87928771972656, '*309': 64.87892150878906, '*31': 64.86991882324219, '*310': 64.87892150878906, '*311': 64.87892150878906, '*312': 64.87892150878906, '*313': 64.87892150878906, '*314': 64.87892150878906, '*315': 64.87873077392578, '*316': 64.87691497802734, '*317': 64.87691497802734, '*318': 64.87628173828125, '*319': 64.87628173828125, '*32': 64.86991882324219, '*320': 64.87628173828125, '*321': 64.87628173828125, '*322': 64.87628173828125, '*323': 64.87628173828125, '*324': 64.87628173828125, '*325': 64.87628173828125, '*326': 64.87628173828125, '*327': 64.87628173828125, '*328': 64.87628173828125, '*329': 64.87628173828125, '*33': 64.86991882324219, '*330': 64.87628173828125, '*331': 64.87628173828125, '*332': 64.87628173828125, '*333': 64.87628173828125, '*334': 64.87628173828125, '*335': 64.87628173828125, '*336': 64.87659454345703, '*337': 64.90101623535156, '*338': 64.90149688720703, '*339': 64.9006118774414, '*34': 64.86991882324219, '*340': 64.9006118774414, '*341': 64.9006118774414, '*342': 64.9006118774414, '*343': 64.9006118774414, '*344': 64.9006118774414, '*345': 64.90101623535156, '*346': 64.90128326416016, '*347': 64.9015884399414, '*348': 64.90193939208984, '*349': 64.90177154541016, '*35': 64.86958312988281, '*350': 64.88825988769531, '*351': 64.88825988769531, '*352': 64.88825988769531, '*353': 64.88825988769531, '*354': 64.88855743408203, '*355': 64.88855743408203, '*356': 64.88855743408203, '*357': 64.88982391357422, '*358': 64.88982391357422, '*359': 64.88941955566406, '*36': 64.87027740478516, '*360': 64.88955688476562, '*361': 64.89034271240234, '*362': 64.89034271240234, '*363': 64.89037322998047, '*364': 64.89037322998047, '*365': 64.89037322998047, '*366': 64.89037322998047, '*367': 64.89057159423828, '*368': 64.89057159423828, '*369': 64.89057159423828, '*37': 64.87027740478516, '*370': 64.89057159423828, '*371': 64.8902816772461, '*372': 64.8902816772461, '*373': 64.89068603515625, '*374': 64.89106750488281, '*375': 64.89106750488281, '*376': 64.89106750488281, '*377': 64.89106750488281, '*378': 64.8904800415039, '*379': 64.88965606689453, '*38': 64.87027740478516, '*380': 64.8896484375, '*381': 64.8896484375, '*382': 64.8896484375, '*383': 64.8896484375, '*384': 64.8897933959961, '*385': 64.8897933959961, '*386': 64.88924407958984, '*387': 64.8895034790039, '*388': 64.88985443115234, '*389': 64.89020538330078, '*39': 64.8707275390625, '*390': 64.89020538330078, '*391': 64.89020538330078, '*392': 64.89020538330078, '*393': 64.89020538330078, '*394': 64.89020538330078, '*395': 64.89041900634766, '*396': 64.89041900634766, '*397': 64.89082336425781, '*398': 64.89054107666016, '*399': 64.89054107666016, '*4': 64.9310531616211, '*40': 64.8707275390625, '*400': 64.89054107666016, '*401': 64.89054107666016, '*402': 64.89054107666016, '*403': 64.89093780517578, '*404': 64.89095306396484, '*405': 64.89095306396484, '*406': 64.89080047607422, '*407': 64.89080047607422, '*408': 64.89080047607422, '*409': 64.89080047607422, '*41': 64.8707275390625, '*410': 64.89080047607422, '*411': 64.89080047607422, '*412': 64.89080047607422, '*413': 64.8907470703125, '*414': 64.8907470703125, '*415': 64.89019012451172, '*416': 64.88300323486328, '*417': 64.88300323486328, '*418': 64.88300323486328, '*419': 64.88300323486328, '*42': 64.87004089355469, '*420': 64.88300323486328, '*421': 64.88300323486328, '*422': 64.88300323486328, '*423': 64.88300323486328, '*424': 64.88300323486328, '*425': 64.88300323486328, '*426': 64.8835678100586, '*427': 64.88374328613281, '*428': 64.88340759277344, '*429': 64.88340759277344, '*43': 64.87004089355469, '*430': 64.88340759277344, '*431': 64.88275909423828, '*432': 64.88275909423828, '*433': 64.88326263427734, '*434': 64.88326263427734, '*435': 64.88326263427734, '*436': 64.88326263427734, '*437': 64.88326263427734, '*438': 64.88326263427734, '*439': 64.88326263427734, '*44': 64.87004089355469, '*440': 64.8836669921875, '*441': 64.8836669921875, '*442': 64.8836669921875, '*443': 64.8836669921875, '*444': 64.88400268554688, '*445': 64.88400268554688, '*446': 64.88275146484375, '*447': 64.82654571533203, '*448': 64.82695007324219, '*449': 64.82695007324219, '*45': 64.87004089355469, '*450': 64.82695007324219, '*451': 64.82649230957031, '*452': 64.82649230957031, '*453': 64.82649230957031, '*454': 64.8271255493164, '*455': 64.82719421386719, '*456': 64.8277359008789, '*457': 64.8277359008789, '*458': 64.8277359008789, '*459': 64.82827758789062, '*46': 64.8690414428711, '*460': 64.82827758789062, '*461': 64.82827758789062, '*462': 64.82827758789062, '*463': 64.82827758789062, '*464': 64.82827758789062, '*465': 64.82827758789062, '*466': 64.82827758789062, '*467': 64.82827758789062, '*468': 64.82827758789062, '*469': 64.82827758789062, '*47': 64.86918640136719, '*470': 64.82827758789062, '*471': 64.82827758789062, '*472': 64.82827758789062, '*473': 64.82827758789062, '*474': 64.82827758789062, '*475': 64.82827758789062, '*476': 64.82827758789062, '*477': 64.82826232910156, '*478': 64.82826232910156, '*479': 64.82791137695312, '*48': 64.8697280883789, '*480': 64.82791137695312, '*481': 64.82791137695312, '*482': 64.82791137695312, '*483': 64.82807922363281, '*484': 64.82807922363281, '*485': 64.8299560546875, '*486': 64.8299560546875, '*487': 64.8299560546875, '*488': 64.8299560546875, '*489': 64.82986450195312, '*49': 64.86937713623047, '*490': 64.82986450195312, '*491': 64.82986450195312, '*492': 64.82986450195312, '*493': 64.82986450195312, '*494': 64.83010864257812, '*495': 64.83010864257812, '*496': 64.8304214477539, '*497': 64.8304214477539, '*498': 64.8304214477539, '*499': 64.83016967773438, '*5': 64.9310531616211, '*50': 64.86937713623047, '*500': 64.83016967773438, '*501': 64.83016967773438, '*502': 64.83016967773438, '*503': 64.83016967773438, '*504': 64.83045959472656, '*505': 64.83045959472656, '*506': 64.83045959472656, '*507': 64.83045959472656, '*51': 64.86937713623047, '*52': 64.86937713623047, '*53': 64.86937713623047, '*54': 64.86937713623047, '*55': 64.86937713623047, '*56': 64.86937713623047, '*57': 64.86937713623047, '*58': 64.86885070800781, '*59': 64.86885070800781, '*6': 64.9310531616211, '*60': 64.86885070800781, '*61': 64.86885070800781, '*62': 64.8682861328125, '*63': 64.8682861328125, '*64': 64.8682861328125, '*65': 64.86774444580078, '*66': 64.86774444580078, '*67': 64.86774444580078, '*68': 64.86774444580078, '*69': 64.86774444580078, '*7': 64.9310531616211, '*70': 64.86693572998047, '*71': 64.86734771728516, '*72': 64.86734771728516, '*73': 64.86736297607422, '*74': 64.86736297607422, '*75': 64.86736297607422, '*76': 64.86736297607422, '*77': 64.86717224121094, '*78': 64.86811065673828, '*79': 64.86811065673828, '*8': 64.9310531616211, '*80': 64.86811065673828, '*81': 64.86811065673828, '*82': 64.86811065673828, '*83': 64.86811065673828, '*84': 64.86811065673828, '*85': 64.86811065673828, '*86': 64.86791229248047, '*87': 64.86791229248047, '*88': 64.86791229248047, '*89': 64.86791229248047, '*9': 64.9310531616211, '*90': 64.86750793457031, '*91': 64.86750793457031, '*92': 64.86750793457031, '*93': 64.86750793457031, '*94': 64.86768341064453, '*95': 64.86768341064453, '*96': 64.86768341064453, '*97': 64.86714172363281, '*98': 64.86714172363281, '*99': 64.86714172363281}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed = self.filter_report(report)

        img = shutil._basename(img)
        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, 'test_mosaic_cube_briggsbwtaper', self.parallel)

        test_dict['test_mosaic_cube_briggsbwtaper']['report'] = report
        test_dict['test_mosaic_cube_briggsbwtaper']['images'] = []

        self.mom8_creator(img+'.image', range_list=[0.15, 1.2])
        self.mom8_creator(img+'.residual', range_list=[0.15, 1.2])
        test_dict['test_mosaic_cube_briggsbwtaper']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_mosaic_cube_briggsbwtaper']['images'].append(img+'.image.profile.png')

        if savemetricdict:
            ### serialize ndarray in mask_stats_dcit
            mask_stats_mod_dict = copy.deepcopy(mask_stats_dict)
            mask_stats_mod_dict['mask'] = mask_stats_dict['mask'].tolist()
            pb_stats_mod_dict = copy.deepcopy(pb_stats_dict)
            pb_stats_mod_dict['pb_mask_0.2'] = pb_stats_dict['pb_mask_0.2'].tolist()
            pb_stats_mod_dict['pb_mask_0.5'] = pb_stats_dict['pb_mask_0.5'].tolist()
            #create a nested dictionary containing exp dictionaries to save
            savedict = {}
            #list of stats to save
            # im_stats, mask_stats, pb_stats, psf_stats,\
            # model_stats, resid_stats, sumwt_stats, + wt_stats (mosaic)
            savedict['im_stats_dict']=im_stats_dict
            savedict['mask_stats_dict']=mask_stats_mod_dict
            savedict['pb_stats_dict']=pb_stats_mod_dict
            savedict['psf_stats_dict']=psf_stats_dict
            savedict['model_stats_dict']=model_stats_dict
            savedict['resid_stats_dict']=resid_stats_dict
            savedict['sumwt_stats_dict']=sumwt_stats_dict
            savedict['wt_stats_dict']=wt_stats_dict
            
            #if self.parallel:
            savedict['bmin_dict']=bmin_dict
            savedict['bmaj_dict']=bmaj_dict
            savedict['pa_dict']=pa_dict

            self.save_dict_to_file(test_name,savedict, test_name+'_cas13317_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_cube_briggsbwtaper
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_mfs(self):
        ''' Mosaic MFS imaging field NGC5363, spw 16 & 22 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))
        #self.save_dict_to_file(report0, 'myreport0')

        # .image report (test_mosaic_mfs)
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[209.02978056deg, 5.25471484deg], [18.8386arcsec, 9.9356arcsec], 90.00000000deg]', field_regions = \
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
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            'im_sum': [False, 1.4970276115098295],
            'regn_sum': [False, 1.58850853646],
            'npts_real': [True, 15876],
            'rms_per_field': [False, [0.0043183333329864697, 0.0035324567542234214, 0.0033643881411162453, 0.0034484378708067886, 0.0035013779149007367, 0.0033359473616749375, 0.0035589233954835845]],
            'fit': [False, [0.03522500582263719, 17.46093579518058, 
                       9.709830310449933]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.8402451422565],
            'fit_pix': [False, [62.9942562846151, 62.995885097033394]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [64, 114, 0, 0]), \
                      (img+'.image', False, [64, 115, 0, 0]), \
                      (img+'.image', True, [11, 60, 0, 0]), \
                      (img+'.image', False, [10, 60, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_mfs)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 15876],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 360],
            'mask_regns': [True, 1],
            'npts_real': [True, 15876]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_mfs)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[209.03003701deg, 5.25471148deg], [110.7528arcsec, 107.4584arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 15876],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 0.200080364943],
            'im_rms': [False, 0.604913988836],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 8454],
            #'npts_0.5': [True, 4497],
            'npts_0.2': [True, 8471],
            'npts_0.5': [True, 4500],
            'npts_real': [True, 15876],
            'fit': [False, [1.0693559655652305, 141.80580479462876, \
                       141.74549135472637]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.8402451422565],
            'fit_pix': [False, [62.975154097364715, 62.94725116661756]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_mfs)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[209.02988854deg, 5.25478653deg], [16.3870arcsec, 10.7097arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            # CAS-9386 update build100 serial
            #'im_sum': [False, 0.00266324029255],
            'im_sum': [False, 0.00074156666480],
            'npts_real': [True, 15876],
            'fit': [False, [1.088207720785799, 15.893701850875548, \
                        8.795192549423799]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.8402451422565],
            'fit_pix': [False, [62.98416058527938, 63.00086190688355]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, \
            '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_mfs)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[209.02978056deg, 5.25471484deg], [18.8386arcsec, 9.9356arcsec], 90.00000000deg]')

        exp_resid_stats = {'npts': [True, 15876],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            # CAS-9386 update build100 serial
            #'im_sum': [False, 0.134135401952],
            'im_sum': [False, 0.131760904704],
            'regn_sum': [False, 0.291504465893],
            'npts_real': [True, 15876]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_mfs)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[209.02978056deg, 5.25471484deg], [18.8386arcsec, 9.9356arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            'mask_non0': [True, 0],
            'npts_real': [True, 15876]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_mfs)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 15849925874.83342],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_mfs)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        exp_wt_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.426215469837],
            'max_val_pos': [True, [63, 63, 0, 0]],
            # CAS-9386 update build100 serial
            #'min_val': [False, 4.20771575591e-05],
            'min_val': [False, 4.26501137554e-05],
            'im_rms': [False, 0.144357241275],
            'im_sum': [False, 1347.6264633],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 8454],
            #'npts_0.5': [True, 4497],
            'npts_0.2': [True, 8471],
            'npts_0.5': [True, 4500],
            'npts_real': [True, 15876]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9
        failed = self.filter_report(report)

        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict['test_mosaic_mfs']['self.parallel'] = self.parallel
        test_dict['test_mosaic_mfs']['report'] = report
        test_dict['test_mosaic_mfs']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.002, 0.035])
        self.mom8_creator(img+'.residual', range_list=[-0.002, 0.035])
        test_dict['test_mosaic_mfs']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            #msg = report)
            msg = failed)

# End of  test_mosaic_mfs
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_mtmfs(self):
        ''' Mosaic mtmfs imaging - field NGC5363, spw 16 & 22 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel,
            verbose=True)

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
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mos_mtmfs'))

        # .image report (test_mosaic_mtmfs)
        im_stats_dict = self.image_stats(img+'.image.tt0', fit_region = \
            'ellipse[[209.02969058deg, 5.25475068deg], [16.7096arcsec, 10.9678arcsec], 90.00000000deg]', \
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
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            'im_sum': [False, 1.4601284610577565],
            'regn_sum': [False, 1.58228409861],
            'npts_real': [True, 15876],
            'rms_per_field': [False, [0.0043714750541124468, 0.003595955649075896, 0.003430812546940405, 0.0034721618952933994, 0.0035884838221589368, 0.0033594003986614097, 0.0036156061407187608]],
            'fit': [False, [0.03580520672441999, 17.19187684101627, \
                       9.68274896612347]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.8402451422565],
            'fit_pix': [False, [63.09049673358014, 62.94805812018937]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [64, 114, 0, 0]), \
                      (img+'.image.tt0', False, [64, 115, 0, 0]), \
                      (img+'.image.tt0', True, [11, 60, 0, 0]), \
                      (img+'.image.tt0', False, [10, 60, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image.tt0', epsilon=self.epsilon)

        # .mask report (test_mosaic_mtmfs)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 15876],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 360],
            'mask_regns': [True, 1],
            'npts_real': [True, 15876]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_mtmfs)
        pb_stats_dict = self.image_stats(img+'.pb.tt0', fit_region = \
            'ellipse[[209.03003701deg, 5.25478205deg], [109.2364arcsec, 107.1964arcsec], 0.00000000deg]')

        exp_pb_stats = {'npts': [True, 15876],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [63, 63, 0, 0]],
            'min_val': [False, 0.200080364943],
            'im_rms': [False, 0.604913988836],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 8454],
            #'npts_0.5': [True, 4497],
            'npts_0.2': [True, 8471],
            'npts_0.5': [True, 4500],
            'npts_real': [True, 15876],
            'fit': [False, [1.0693559655651996, 141.80580479464936, \
                       141.74549135470988]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.8402451422565],
            'fit_pix': [False, [62.975154097364715, 62.94725116661756]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb.tt0', epsilon=self.epsilon)

        # .psf report test_mosaic_mtmfs)
        psf_stats_dict = self.image_stats(img+'.psf.tt0', fit_region = \
            'ellipse[[209.02997853deg, 5.25475069deg], [16.3225arcsec, 9.9355arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            # CAS-9386 update build100 serial
            #'im_sum': [False, 0.00266324029255],
            'im_sum': [False, 0.00074156666480],
            'npts_real': [True, 15876],
            'fit': [False, [1.0781857293103545, 15.898196388608632, \
                        8.995969894587292]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 107.8402451422565],
            'fit_pix': [False, [62.991298508308404, 63.00339664380328]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf.tt0', epsilon=self.epsilon)

        # .residual report test_mosaic_mtmfs)
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            'ellipse[[209.02969058deg, 5.25475068deg], [16.7096arcsec, 10.9678arcsec], 90.00000000deg]')

        exp_resid_stats = {'npts': [True, 15876],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            # CAS-9386 update build100 serial
            #'im_sum': [False, 0.130538249029],
            'im_sum': [False, 0.128005618624],
            'regn_sum': [False, 0.283332834988],
            'npts_real': [True, 15876]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0', epsilon=self.epsilon)

        # .model report test_mosaic_mtmfs)
        model_stats_dict = self.image_stats(img+'.model.tt0', fit_region = \
            'ellipse[[209.02969058deg, 5.25475068deg], [16.7096arcsec, 10.9678arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
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
            'mask_non0': [True, 0],
            'npts_real': [True, 15876]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt report test_mosaic_mtmfs)
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 15849925874.83342],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0', epsilon=self.epsilon)

        # .weight report test_mosaic_mtmfs)
        wt_stats_dict = self.image_stats(img+'.weight.tt0', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        exp_wt_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.426215469837],
            'max_val_pos': [True, [63, 63, 0, 0]],
            # CAS-9386 update build100 serial
            #'min_val': [False, 4.20771575591e-05],
            'min_val': [False, 4.26501137554e-05],
            'im_rms': [False, 0.144357241275],
            'im_sum': [False, 1347.6264633],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 8454],
            #'npts_0.5': [True, 4497],
            'npts_0.2': [True, 8471],
            'npts_0.5': [True, 4500],
            'npts_real': [True, 15876]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, \
            '.weight.tt0', epsilon=self.epsilon)

        # .image.tt1 report test_mosaic_mtmfs)
        im1_stats_dict = self.image_stats(img+'.image.tt1', fit_region = \
            'ellipse[[209.02969058deg, 5.25475068deg], [16.7096arcsec, 10.9678arcsec], 90.00000000deg]', \
            field_regions = \
            ['circle[[13:56:07.210000, +05.15.17.20000], 45.9arcsec]',
             'circle[[13:56:06.355525, +05.15.59.74129], 45.9arcsec]',
             'circle[[13:56:04.316267, +05.15.27.41716], 45.9arcsec]',
             'circle[[13:56:09.249291, +05.15.49.52355], 45.9arcsec]',
             'circle[[13:56:05.170768, +05.14.44.87604], 45.9arcsec]',
             'circle[[13:56:10.103706, +05.15.06.98201], 45.9arcsec]',
             'circle[[13:56:08.064442, +05.14.34.65864], 45.9arcsec]'])

        exp_im1_stats = {'com_bmaj': [False, 17.6737785339],
            'com_bmin': [False, 10.060172081],
            'com_pa': [False, 86.6785964966],
            'npts': [True, 15876],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0410945117474],
            'max_val_pos': [True, [53, 72, 0, 0]],
            'min_val': [False, -0.0436826199293],
            'min_val_pos': [True, [70, 97, 0, 0]],
            'im_rms': [False, 0.0119215006547],
            # CAS-9386 update build100 serial
            #'im_sum': [False, -0.99675725757],
            #'regn_sum': [False, 0.144841228408],
            'im_sum': [False, -1.05556125437],
            'regn_sum': [False, 0.12572287418879569],
            'npts_real': [True, 15876],
            'rms_per_field': [False, [0.0118926721315, 0.0131530097001, 0.0123432407276, 0.0117928565232, 0.0110465636431, 0.0122420920176, 0.012233014507]]}

        report10 = th.check_dict_vals(exp_im1_stats, im1_stats_dict, '.image.tt1', epsilon=self.epsilon)

        # .residual.tt1 report test_mosaic_mtmfs)
        resid1_stats_dict = self.image_stats(img+'.residual.tt1', \
            'ellipse[[209.02969058deg, 5.25475068deg], [16.7096arcsec, 10.9678arcsec], 90.00000000deg]')

        exp_resid1_stats = {'npts': [True, 15876],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 8454.0],
            'npts_unmasked': [True, 8471.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.000123286881717],
            'max_val_pos': [True, [51, 71, 0, 0]],
            'min_val': [False, -0.000183079260751],
            'min_val_pos': [True, [61, 63, 0, 0]],
            'im_rms': [False, 4.57978644843e-05],
            # CAS-9386 update build100 serial
            #'im_sum': [False, -0.00662636467916],
            'im_sum': [False, -0.00669605520018],
            'regn_sum': [False, -0.0077922191209],
            'npts_real': [True, 15876]}

        report11 = th.check_dict_vals(exp_resid1_stats, resid1_stats_dict, \
            '.residual.tt1', epsilon=self.epsilon)

        # .model.tt1 report test_mosaic_mtmfs)
        model1_stats_dict = self.image_stats(img+'.model.tt1', fit_region = \
            'ellipse[[209.02969058deg, 5.25475068deg], [16.7096arcsec, 10.9678arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model1_stats = {'npts': [True, 15876],
            'npts_unmasked': [True, 15876.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.0],
            'max_val_pos': [True, [0, 0, 0, 0]],
            # CAS-9386 update build100 serial
            #'min_val': [False, -0.00433648051694],
            'min_val': [False, -0.00455870199949],
            'min_val_pos': [True, [63, 63, 0, 0]],
            # CAS-9386 update build100 serial
            #'im_rms': [False, 3.44165120392e-05],
            #'im_sum': [False, -0.00433648051694],
            #'regn_sum': [False, -0.00433648051694],
            'im_rms': [False, 3.61801745991e-05],
            'im_sum': [False, -0.00455870199949],
            'regn_sum': [False, -0.00455870199949],
            'mask_non0': [True, 0],
            'npts_real': [True, 15876]}

        report12 = th.check_dict_vals(exp_model1_stats, model1_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt.tt1 report test_mosaic_mtmfs)
        sumwt1_stats_dict = self.image_stats(img+'.sumwt.tt1')

        exp_sumwt1_stats = {'npts': [True, 1],
            'npts_unmasked': [True, 1.0],
            'freq_bin': [1e-10, 15849925874.83342],
            'start': [True, 1.0784e+11],
            'end': [True, 1.0784e+11],
            'start_delta': [False, 1.0784e+11],
            'end_delta': [False, 1.0784e+11],
            'nchan': [True, 1],
            'max_val': [False, -123654.40625],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, -123654.40625],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 123654.40631],
            'npts_real': [True, 1]}

        report13 = th.check_dict_vals(exp_sumwt1_stats, sumwt1_stats_dict, \
            '.sumwt.tt1', epsilon=self.epsilon)

        # report combination test_mosaic_mtmfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11 + report12 + report13

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict['test_mosaic_mtmfs']['self.parallel'] = self.parallel
        test_dict['test_mosaic_mtmfs']['report'] = report
        test_dict['test_mosaic_mtmfs']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image.tt0', range_list=[-0.003, 0.035])
        self.mom8_creator(img+'.residual.tt0', range_list=[-0.003, 0.035])
        test_dict['test_mosaic_mtmfs']['images'].extend( \
            (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_mtmfs

#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_cube_eph(self):
        ''' Mosaic ephemeris cube imaging - field Venus, spw 45 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel, verbose=True)

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
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))
  
        # test_mosaic_cube_eph
        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        # test_mosaic_cube_eph
        exp_im_stats = {'com_bmaj': [False, 0.9347811147697555],
            'com_bmin': [False, 0.7185930247591006],
            'com_pa': [False, -88.16793752979638],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 104998085.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.086364768445491791],
            'max_val_pos': [True, [318, 243, 0, 511]],
            'min_val': [False, -0.070862673223018646,],
            'min_val_pos': [True, [222, 212, 0, 674]],
            'im_rms': [False, 0.012281882598118963],
            'rms_per_chan': [False, [0.01566755765592744, 0.012499771976912814, 0.011123979147164407, 0.011561886220609966, 0.012292568630946268, 0.013503738082575103, 0.013895709096621548, 0.015371602081183504, 0.012571966053954607, 0.012110657886905454, 0.010681959860783033, 0.011448528933751815, 0.013363157615506101, 0.01177851095735433, 0.011819999220451467, 0.011063128145379053, 0.011211424009495315, 0.012675826502373685, 0.00964365231187898, 0.01230416428699256, 0.011816948425841614, 0.012524610132431753, 0.012091583554529339, 0.01415517292636185, 0.013038519870781123, 0.010927180475068366, 0.011348376416636784, 0.011452376144559037, 0.013880464215687604, 0.01271421858369905, 0.012365348078626935, 0.011134439502267753, 0.01173999969262898, 0.010275262798856601, 0.012024637616711497, 0.011827779313317115, 0.012005719677964755, 0.01126998085149768, 0.01263100142352357, 0.010919404120788623, 0.01230292923375664, 0.013824000423944424, 0.012030441249337255, 0.01471235121764837, 0.012690250219699167, 0.012264732434797237, 0.010338113974845756, 0.010588074612506274, 0.010658916275604544, 0.01192531675206317, 0.011963523090655131, 0.011739554348224612, 0.010864559965038227, 0.01135082917289156, 0.011495290472285927, 0.011805507617290682, 0.010933831622378246, 0.014576785073645402, 0.01297711633387727, 0.010222489107261269, 0.0110984415724654, 0.011692109085082973, 0.011605118458967478, 0.011686408270450569, 0.010719138244995948, 0.011490010832745277, 0.012016309876685964, 0.012905508417480243, 0.010424679511378012, 0.010476580776329382, 0.013555509132609708, 0.01168612606268359, 0.01089950574737809, 0.009744535896916217, 0.011655046195014002, 0.012674198608575759, 0.014310946665162668, 0.014293294641270153, 0.012638247297143402, 0.011225312863067707, 0.010533579687603354, 0.010619978615716235, 0.01130230533155928, 0.012504854882493843, 0.011302020261357171, 0.012161798863089938, 0.010608102908705218, 0.011122188031292404, 0.010132879137328434, 0.011293949972102117, 0.012010374392817949, 0.013191568018758376, 0.01233135069689494, 0.010769746462901904, 0.01153549017814873, 0.011979399621333483, 0.011111746944971652, 0.011843988636800718, 0.011876601129763641, 0.00985833760297452, 0.012805113173417818, 0.012001866506178401, 0.012972732329926176, 0.012947046238817135, 0.01169227386437645, 0.013660834492575508, 0.015240297428641192, 0.012045506099515298, 0.014349189506965473, 0.011753680556779432, 0.011011732485486841, 0.011048193549742145, 0.012814236743542255, 0.011674542313321304, 0.010770073501250303, 0.012405362003451131, 0.011262904273749614, 0.010016645329827439, 0.012059229625014212, 0.01128781720269892, 0.011844913248897927, 0.012245345012809291, 0.009952101453716862, 0.012237226805610778, 0.010830186561964461, 0.012038505819334888, 0.010469212178041205, 0.012712147168754424, 0.009858813880899741, 0.009618511420097878, 0.011251630406509621, 0.01094736018021904, 0.009938750246980809, 0.011806500650862979, 0.013746618610294636, 0.011890591602572501, 0.011285353902870112, 0.010746148131018936, 0.010780721254182541, 0.010400367901096838, 0.010548644745226116, 0.011701519719661485, 0.011515444299289019, 0.01137979009080972, 0.009897580126157974, 0.012679387512909252, 0.012393749358995263, 0.010126541811740334, 0.012078757530494299, 0.010422834662425521, 0.011491611179879445, 0.011350002754119041, 0.010786499538647526, 0.01241069099315559, 0.01287288681929542, 0.012544796642462093, 0.012486480230120337, 0.012027122098965512, 0.011436187647136224, 0.014507824634016918, 0.012781326055012877, 0.010971035372093395, 0.011222652982616295, 0.010381402164033299, 0.011789211848634065, 0.010467460759016958, 0.010519913231344411, 0.009804241543959213, 0.009826802489597598, 0.011168128286266319, 0.009898755834256175, 0.010638546697493598, 0.01075904603518445, 0.01186190395476391, 0.012698873200736769, 0.012275066865331018, 0.010120998522019762, 0.011856430473088262, 0.01187640156331042, 0.012858360146453808, 0.010458450314412523, 0.011686631685348224, 0.012324510594898344, 0.013706507117441952, 0.010507073125858515, 0.009727999262257602, 0.010749737212657888, 0.010496316622257333, 0.012266025310015635, 0.009968887987260507, 0.010656267754294232, 0.01137826346655607, 0.013691020253483116, 0.009712315969605423, 0.0116765832200688, 0.01204200794054857, 0.011085824144059132, 0.012319945895455337, 0.010766015294335025, 0.012434827540397717, 0.012730633406142734, 0.014922472720576237, 0.013340994079356561, 0.012184668987628627, 0.011860251479407752, 0.010611453465531486, 0.010927075341466757, 0.010737943994644793, 0.011284903658799174, 0.012101560213749602, 0.01253847688484235, 0.011820390980633899, 0.01141402994073023, 0.01060976878040442, 0.010983370981323954, 0.013427371951491312, 0.012374725410784092, 0.011400378106504076, 0.011527234422091191, 0.01340962759016125, 0.013570332253931302, 0.011804029668992844, 0.012197974733541946, 0.013983623168953468, 0.012725763762714332, 0.011354745750571254, 0.012145679916283425, 0.01126874546535086, 0.012261492463656863, 0.011178855897759667, 0.010264367951585555, 0.013345043592089396, 0.012409363114615287, 0.011495934975609993, 0.009277679225470462, 0.011720945646776957, 0.012561715975156802, 0.012325900140530158, 0.011173759164206036, 0.01444252189374611, 0.013235686426999993, 0.009545194033680014, 0.00989282465604266, 0.011317913505981915, 0.011172814496764794, 0.011097429969483154, 0.011431240017029451, 0.010546755637817182, 0.011214435474509551, 0.011763583664075946, 0.014689109193649632, 0.011251056734311256, 0.009776300805049154, 0.012860933158710076, 0.012509654758086465, 0.013307995464813777, 0.012069070120183812, 0.011422201031864795, 0.011512246169509426, 0.011243380202117025, 0.011621384685604818, 0.011008804779668683, 0.013029468912847373, 0.011484166233438071, 0.01154715840745958, 0.012030983875891486, 0.012644418892294014, 0.010840072063581964, 0.01217667528470696, 0.012348270662399291, 0.011295167588599378, 0.011131983820484052, 0.011525183914729311, 0.01245079346986009, 0.010668667918574357, 0.011110716151104136, 0.012536736404100099, 0.013297691638473148, 0.013088181650764244, 0.011822130327201136, 0.010774771031200035, 0.010775807119122434, 0.01105586988033813, 0.009420738313971375, 0.012013742695805617, 0.011326747652322514, 0.012096904239732881, 0.01150806457510923, 0.011675274227450103, 0.0135991855834598, 0.01044590308467654, 0.01066146998192241, 0.010712518291299183, 0.01193296299641928, 0.012947758061995076, 0.010130158406514036, 0.011855179848796187, 0.011440889779476115, 0.012458460329640899, 0.011076297729893039, 0.012814778866789587, 0.011811606075871006, 0.012867705357363268, 0.01621962105427537, 0.016465436698466954, 0.016201295070531984, 0.01058752407190847, 0.011716930684685505, 0.011008111220908223, 0.01309686038666478, 0.011249739924392245, 0.011298808361717275, 0.01286619495044334, 0.012457597012475598, 0.011357408247558936, 0.011635136666110633, 0.01147265786957545, 0.01411270947674405, 0.012138457217748517, 0.012013161207475307, 0.012399897890804487, 0.017261754698869394, 0.011751036356983229, 0.01220773352984818, 0.01169628539620397, 0.013572907681349047, 0.014279084224693019, 0.013350071117589353, 0.01119955478052971, 0.011678178134442125, 0.01143788295713284, 0.01338236323998995, 0.013148991760974126, 0.012261293827211488, 0.012230033568364952, 0.011982727185002299, 0.014628741472584738, 0.013936276355000195, 0.011606244163871274, 0.013217927391198507, 0.011877296698841999, 0.013751323956061321, 0.010502431019977986, 0.011623558011397342, 0.013225657353801752, 0.013687394057886515, 0.012047896602847821, 0.01247678361440421, 0.010956367674752512, 0.011928249343126845, 0.012456739592702704, 0.014412803794511948, 0.013091054766212506, 0.012391320075065916, 0.012041284843754036, 0.011175660150719864, 0.012796161329318703, 0.012376600001200107, 0.013213434092613835, 0.011391999754128477, 0.011735911173241135, 0.011338817659054366, 0.013685494524633305, 0.01405734652024605, 0.011873232352346962, 0.011492726398798213, 0.011357899840429751, 0.011924626623408293, 0.01292314975466112, 0.011568890237960018, 0.012450768675985815, 0.01615561678574021, 0.014389883701947508, 0.012020341777831386, 0.011691722848726802, 0.012565949873863497, 0.011694804496294954, 0.011846565903593035, 0.011247492738346963, 0.010809808599760048, 0.011330103597179153, 0.012046098499893751, 0.012792309821515609, 0.012055452071571462, 0.012529809308170823, 0.01166776384371269, 0.01271547225436619, 0.014155919212706156, 0.011590774944224877, 0.013272175555436434, 0.01230731481769422, 0.01134847189287814, 0.011467954772903752, 0.01184115575454566, 0.011988620387346051, 0.010501901957914332, 0.011616509715974813, 0.011471135705227897, 0.012451191281457752, 0.015319504006822212, 0.012649022813805498, 0.011475752166486771, 0.012390712505079823, 0.013092398166506612, 0.0108425779557803, 0.010850393792933026, 0.011858723037825333, 0.014367570168849372, 0.016151613939500507, 0.010857278233794036, 0.012348584717685072, 0.013369893472750983, 0.013520298224929771, 0.012490511233733568, 0.011838067139252433, 0.01379337132839034, 0.013615981503194808, 0.01213329972068078, 0.011926986383721444, 0.011989638160863852, 0.011459946102237622, 0.00998617556766575, 0.011909475646001892, 0.01144107317570622, 0.012796872801637271, 0.011812967842105799, 0.01226189488300729, 0.013800981007306753, 0.011403886599520798, 0.011873516666846865, 0.012409987676322875, 0.012188938970018453, 0.011353880208547378, 0.013058666569184637, 0.011338981854502068, 0.01409796019586988, 0.012231139093670572, 0.015107895235161905, 0.012326810087808738, 0.012946694516644243, 0.012447550801307723, 0.014909071523596254, 0.01397117897114995, 0.013954030858987044, 0.012642975889461814, 0.012387518693283873, 0.011610200546359207, 0.011526268526938685, 0.011758827798373622, 0.011908680699262349, 0.015003488394504251, 0.014659137370125483, 0.01303276571543743, 0.013752408592063675, 0.013708881194186243, 0.01147231577101178, 0.011546435459926428, 0.011541811993355754, 0.013315481649931753, 0.012924077797838931, 0.013861771686590075, 0.012908741618203564, 0.012185556034842455, 0.0123901699294395, 0.011089612618787055, 0.012565991045585109, 0.011172521028748187, 0.011930136297848918, 0.012217995900935783, 0.012567934765036304, 0.013402949630010139, 0.012677891119972322, 0.012831082643587562, 0.012077587392213665, 0.01192202179044923, 0.012437303467243326, 0.012739158128205348, 0.013229922671982922, 0.011707021012177004, 0.013242901791918983, 0.012429029653969133, 0.011897651409429517, 0.012333981656159804, 0.01652056933731621, 0.012803033285199423, 0.013291710033149702, 0.013915563269442936, 0.013293192094911472, 0.012809923485290182, 0.012329480430228985, 0.012327561250559955, 0.011635712112230329, 0.013309364248632287, 0.01298532449779078, 0.011141245068791102, 0.0123662111290416, 0.012952937497979091, 0.012652571307871136, 0.015633401429410065, 0.014029079448528728, 0.013085527098929505, 0.013671367699869673, 0.011449641495798529, 0.013870108540867539, 0.013267219241966837, 0.013009671745599366, 0.012061924878367703, 0.013344180204315093, 0.012472386713589599, 0.012819448624406097, 0.013468487050066785, 0.015061287987425221, 0.013101479251036413, 0.011695612693457496, 0.014215627129265261, 0.014030414081014604, 0.013144534677081394, 0.013769631498326561, 0.013665079857862566, 0.01278869825702419, 0.011758045653360501, 0.013379677594574436, 0.013176320438872265, 0.013726863311133137, 0.014093922857692474, 0.011978594418837533, 0.012688741742419463, 0.014213190054434053, 0.012778188415525355, 0.012893566603547115, 0.011031522915526809, 0.012099967208788644, 0.01197488835353761, 0.011886891227553936, 0.012806477665845199, 0.013276830507854871, 0.013241164133079063, 0.012154967597401761, 0.011697698448488806, 0.0138159512574303, 0.012144528154314117, 0.011982896302520097, 0.012088015079846063, 0.01258443418498413, 0.012417041503252886, 0.013316302302718812, 0.013253698599679495, 0.014875163789523499, 0.013611399587966033, 0.01270670566020336, 0.012428452848201328, 0.013258866017540937, 0.012735738693151316, 0.014906643094604965, 0.014207002883168778, 0.01351154568809692, 0.012428703834641097, 0.012531695617773475, 0.014339355016753071, 0.011835752014891681, 0.015198065031029119, 0.01430821328282015, 0.013257734215508951, 0.012994339259791049, 0.013865856665286771, 0.012158304848082375, 0.014198843802734912, 0.01419354016501019, 0.013678418586819103, 0.012056908737372238, 0.012088134020138442, 0.012708550023920906, 0.012056042831971507, 0.012944165533598666, 0.012998439619323971, 0.011738401817877483, 0.013534728428647752, 0.012933197018222002, 0.011888680423975263, 0.012833718856536632, 0.013026329306504825, 0.01389771764379911, 0.013170374270730664, 0.01394810209074441, 0.013726597650998646, 0.012620621703300312, 0.012195206631413907, 0.013644885542175174, 0.01339675796598244, 0.015953133266745596, 0.015572165751410377, 0.013918610648959726, 0.012539523963263447, 0.014863648164363799, 0.013557760935362215, 0.013767615205969776, 0.012387780847253898, 0.013111551187613596, 0.012178574316034968, 0.012804589291833082, 0.013764797217935443, 0.014682794411092211, 0.011923563140996867, 0.012457566642369073, 0.012577155733152602, 0.012268408936443286, 0.013433767866486363, 0.011254144886535153, 0.012239416835082412, 0.012719204139179307, 0.014153628502136058, 0.013714206504685713, 0.01555664321990409, 0.012482696827891392, 0.012798094219294787, 0.013725222959993684, 0.011297342405230785, 0.013866245776811888, 0.013948944319164524, 0.014752096346863636, 0.012392843752160694, 0.014416608533154353, 0.013416960727895713, 0.013014075062578258, 0.012457948814251782, 0.01214561019554649, 0.014573596703768801, 0.013355931500995009, 0.010935442745219525, 0.013474856369494768, 0.014093982151433494, 0.011282427815945854, 0.013022998498462857, 0.011932228038133832, 0.011193973084211214, 0.013372017260001274, 0.013994297538831645, 0.011742113041326014, 0.012987909251690535, 0.012569010525416004, 0.011966315747010852, 0.013958876937310783, 0.013110867016705616, 0.015176327530987377, 0.013271913419261386, 0.014815021510858065, 0.013019719016438232, 0.012593706022601611, 0.013636464785719227, 0.013024015052417843, 0.011867287358548854, 0.014432494354300884, 0.011513467360750092, 0.013751315885178503, 0.011944914326008424, 0.011579770223979212, 0.013104138863451702, 0.012572035675858769, 0.013719589991898679, 0.011799423913740791, 0.014222826331815612, 0.013130190236718657, 0.015687333785054065, 0.013977869745415655, 0.013856666637965972, 0.0119206019645679, 0.01340210329202167, 0.01236780694175228, 0.013548894032861069, 0.012186169719114176, 0.01048406598287907, 0.010655269577945959, 0.011891505227937497, 0.012232046724993523, 0.016338504257622802, 0.011951257476223706, 0.012787444872754723, 0.011710212821400005, 0.011695865895061741, 0.015838707904307876, 0.012939755683705046, 0.012477587036156625, 0.013330358565129482, 0.01082528855224014, 0.013773287643131382, 0.012570171555245143, 0.013162359160576995, 0.011218965550597317, 0.01342890868237659, 0.01202182099406552, 0.010697138162887927, 0.015232948941257044, 0.012711001004850505, 0.011506289930999608, 0.012369850055536856, 0.01233361532935207, 0.012360627932656498, 0.013326575850585057, 0.012934710498929827, 0.013248039695274097, 0.012226203594835677, 0.012818113815664448, 0.01341051851094283, 0.013541762993268074, 0.013056711886632514, 0.012115662158377625, 0.013357590759015063, 0.012032645498305006, 0.016207947768265084, 0.01389311652405909, 0.013847679128406589, 0.013886386730554553, 0.012067502182802831, 0.011753272427404053, 0.013279048057631762, 0.011214754296383541, 0.011683816996863522, 0.010609551860508486, 0.012357120345452005, 0.013050369382881564, 0.01174704009897493, 0.011743154945937362, 0.012531199758011487, 0.011798725028591634, 0.012419735337445137, 0.013975728513221264, 0.011781109647231187, 0.014598054146810256, 0.011754504071338916, 0.011305854302102692, 0.013318945027316854, 0.013042196249451354, 0.011052738741563284, 0.010558225013410463, 0.011020779774382213, 0.011340319330782552, 0.012033592610864789, 0.012455257536242884, 0.01110139391512171, 0.01046562816674095, 0.015266557792594495, 0.011918313536064278, 0.011360562910751727, 0.010365874108060247, 0.011737890229192168, 0.012793643279848948, 0.013256120317317155, 0.012383560012815365, 0.013301087661143126, 0.01531080427533124, 0.010224017979386747, 0.011257686276499953, 0.009965303526327884, 0.011481057685759922, 0.01274428412077991, 0.01237779048897638, 0.011054791139194602, 0.012765273103949742, 0.012472945953336118, 0.012608159036016397, 0.011204550231481033, 0.013905107233033231, 0.010773937005630393, 0.013310027013851727, 0.017420680889573912, 0.011673620547153839, 0.010614208861963346, 0.012253262631090077, 0.012560836225969919, 0.012560372974566972, 0.013315953743992456, 0.01250495576610427, 0.010300621736575054, 0.0126743669068127, 0.012497189500451297, 0.01058986672748917, 0.0140554647719132, 0.011586323266949704, 0.01051084444510618, 0.010723172131198219, 0.013670881790176151, 0.01294181100319149, 0.01329464895892603, 0.011257304758166548, 0.011406076857722271, 0.01192400994821658, 0.01090031653434306, 0.011516503202529208, 0.01164722421363161, 0.011278034900260045, 0.013101497734367626, 0.010571820551740189, 0.01258494483689583, 0.010690018512359867, 0.012288203401799393, 0.011394189183191222, 0.010661795126849621, 0.014600444436046287, 0.012545075612095091, 0.011433423838627848, 0.0108797962059277, 0.01179905039017321, 0.011201593386431653, 0.012398023234344002, 0.011585602038687984, 0.012105760854189395, 0.01012363568557298, 0.010608758997456596, 0.013055946950589595, 0.01115366731042303, 0.011356175896159434, 0.010664413285359914, 0.011095187626260742, 0.01055909008819056, 0.01199264964050027, 0.01198885985504728, 0.011856295143569895, 0.011865458906234921, 0.012803395900220201, 0.011576869098809376, 0.010189936253659717, 0.010195412935336356, 0.012280457502570626, 0.011418747542118402, 0.011054267711398783, 0.010663172535730193, 0.010788394178560413, 0.012825843175935922, 0.014830488452138243, 0.012284354036472657, 0.012694714212578187, 0.013724353996911465, 0.012851827129972288, 0.010592820880735418, 0.012302591864876384, 0.012065809850599015, 0.011956358889472574, 0.010872197413464683, 0.012323221424305216, 0.012123107304017921, 0.010815605182278264, 0.012638711922615793, 0.011655868729735695, 0.013318117884136317, 0.011374335257623935, 0.011336547412403876, 0.013538061509343441, 0.010884298085465788, 0.00983147858074236, 0.010395424381769653, 0.010546991873577608, 0.012797725536117008, 0.011287668516950705, 0.010972896319568352, 0.011817914808157692, 0.010995529479193617, 0.01009639486723553, 0.01141833937111661, 0.011690435639595306, 0.010061766319288705, 0.011971105561760871, 0.012700043131615801, 0.011318883390568383, 0.010909315349852766, 0.00955009744793359, 0.010498572864691814, 0.013751363636934978, 0.012414462975548525, 0.013914888435336918, 0.014227696895603312, 0.012244996318125041, 0.01179968071369893, 0.01269311294414022, 0.011707870106088975, 0.010958813990455807, 0.011651043993684336, 0.01030137511322718, 0.011034955595029074, 0.011010590344592711, 0.012071622334637981, 0.012156945636578436, 0.01207776327788143, 0.012467860244637925, 0.011202073492240526, 0.013522425017327757, 0.013101403096011953, 0.011576541797128533, 0.010161447675978133, 0.010808172743640786, 0.01053092595610692, 0.011264161352788734, 0.01233589290036475, 0.010814458779802489, 0.012333845855799757, 0.009642890958343297, 0.012173492497045947, 0.013243881943014759, 0.011496651178443193, 0.01015960484422759, 0.01181458373889043, 0.009783393881506904, 0.00995620066577112, 0.011522308540852155, 0.010744486398458035, 0.010880760248384722, 0.0103126533049451, 0.010249354156011579, 0.011065664705266258, 0.01142653024289266, 0.01119440735004512, 0.01097503503562323, 0.014023046074535458, 0.010432042532047749, 0.010194842729529698, 0.010817754756984273, 0.011134926202450978, 0.011429007186568798, 0.010583381706566801, 0.015263726359132578, 0.011921804488613408, 0.009753228291247933, 0.012196147705148373, 0.009316733907216289, 0.01146380450784002, 0.010523188985348219, 0.012950242966562036, 0.010905274766360632, 0.009936978616187102, 0.011966920815569286, 0.010385478642918732, 0.010089607608523066, 0.01317787233691099, 0.01482152961916169, 0.011166319039339316, 0.01158412070404005, 0.010686371837829486, 0.011619838402650143, 0.010062419209719961, 0.010294116467679117, 0.01054194871478958, 0.011768968106377857, 0.012003013751281379, 0.011635692550686211, 0.012427950542111656, 0.012436753739762036, 0.010597044719259695, 0.01047844931079727, 0.01300316381280212]],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 4241.3850895847481],
            'im_sum': [False, 4107.755953973777],
            'regn_sum': [False, 34.376076006740732],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.013727427598497976, 0.01265603739684268, 0.012647516154896329, 0.013363016802334838, 0.012374606023989609, 0.013129803882930785, 0.01267018256172788]],
            'profile': [False, 0.025674302295194625]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [240, 393, 0, 0]), \
                      # CAS-9386 update build100 serial
                      #(img+'.image', False, [240, 394, 0, 0]), \
                      (img+'.image', True, [240, 394, 0, 0]), \
                      (img+'.image', True, [49, 209, 0, 0]), \
                      (img+'.image', False, [48, 209, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_cube_eph)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 191116800],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            # CAS-9386 update build100 serial
            #'mask_pix': [False, 8519],
            'mask_pix': [False, 8194],
            'mask_regns': [True, 31],
            'npts_real': [True, 191116800]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 104998085.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.200000017881],
            'im_rms': [False, 0.638163448188254],
            'npts_0.2': [False, [110768, 110768, 110768, 110767, 110766, 110766, 110765, 110765, 110765, 110765, 110766, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110768, 110768, 110768, 110766, 110766, 110768, 110767, 110768, 110768, 110767, 110767, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110766, 110770, 110770, 110768, 110768, 110770, 110771, 110771, 110771, 110770, 110770, 110770, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110769, 110769, 110769, 110768, 110768, 110768, 110769, 110771, 110771, 110771, 110770, 110770, 110770, 110770, 110770, 110770, 110769, 110767, 110769, 110768, 110769, 110770, 110772, 110772, 110772, 110772, 110772, 110772, 110773, 110773, 110773, 110774, 110774, 110773, 110771, 110772, 110771, 110770, 110770, 110770, 110769, 110771, 110771, 110769, 110770, 110769, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110769, 110769, 110769, 110769, 110770, 110770, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110770, 110769, 110769, 110768, 110768, 110768, 110769, 110768, 110767, 110767, 110768, 110767, 110768, 110767, 110767, 110767, 110768, 110769, 110769, 110768, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110766, 110766, 110766, 110767, 110769, 110769, 110770, 110771, 110770, 110769, 110770, 110770, 110769, 110768, 110768, 110769, 110769, 110769, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110770, 110770, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110771, 110771, 110771, 110773, 110774, 110773, 110774, 110776, 110775, 110775, 110773, 110773, 110774, 110774, 110774, 110775, 110774, 110774, 110776, 110776, 110776, 110777, 110776, 110776, 110776, 110776, 110776, 110776, 110775, 110773, 110773, 110773, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110775, 110774, 110774, 110775, 110775, 110775, 110776, 110777, 110777, 110777, 110777, 110778, 110778, 110779, 110778, 110780, 110781, 110781, 110781, 110779, 110780, 110780, 110781, 110779, 110779, 110779, 110779, 110779, 110782, 110782, 110781, 110781, 110782, 110782, 110782, 110782, 110782, 110782, 110781, 110783, 110782, 110782, 110782, 110783, 110782, 110782, 110783, 110780, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110782, 110782, 110781, 110781, 110781, 110781, 110781, 110781, 110781, 110782, 110781, 110781, 110784, 110782, 110782, 110784, 110782, 110782, 110782, 110782, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110783, 110783, 110782, 110783, 110783, 110784, 110784, 110784, 110786, 110786, 110785, 110786, 110786, 110786, 110786, 110785, 110785, 110784, 110785, 110785, 110785, 110785, 110785, 110785, 110786, 110787, 110787, 110787, 110787, 110786, 110786, 110786, 110786, 110786, 110787, 110787, 110787, 110787, 110787, 110786, 110787, 110787, 110786, 110786, 110785, 110785, 110786, 110786, 110786, 110786, 110786, 110784, 110786, 110786, 110786, 110786, 110785, 110786, 110786, 110785, 110785, 110782, 110783, 110785, 110785, 110786, 110785, 110785, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110783, 110784, 110783, 110783, 110783, 110784, 110784, 110784, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110746, 110746, 110746, 110745, 110746, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110748, 110748, 110748, 110748, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110747, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110747, 110745, 110745, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110745, 110745, 110744, 110745, 110744, 110744, 110744, 110745, 110745, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110743, 110743, 110743, 110743, 110744, 110741, 110743, 110743, 110743, 110740, 110741, 110740, 110741, 110743, 110743, 110743, 110741, 110741, 110741, 110741, 110741, 110740, 110741, 110740, 110740, 110741, 110741, 110740, 110741, 110741, 110742, 110742, 110742, 110742, 110741, 110742, 110742, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110741, 110741, 110741, 110742, 110742, 110741, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110742, 110742, 110742, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110742, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110744, 110744, 110744, 110742, 110741, 110741, 110742, 110742, 110740, 110740, 110741, 110741, 110741, 110740, 110740, 110740, 110739, 110740, 110739, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110739, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110736, 110738, 110736, 110738, 110737, 110738, 110738, 110738, 110738, 110739, 110739, 110736, 110736, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110738, 110737, 110737, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110738, 110737, 110737, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110739, 110739, 110738, 110738, 110738, 110736, 110736, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110734, 110733, 110734, 110733, 110735, 110735, 110735, 110735, 110736, 110737, 110736, 110736, 110736, 110736, 110735, 110736, 110736, 110733, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110733, 110734, 110734, 110734, 110733, 110733, 110734, 110734, 110734, 110732, 110734, 110733, 110733, 110731, 110732, 110731, 110732, 110733, 110733, 110733, 110731, 110731, 110732, 110730, 110731, 110733, 110732, 110732]],
            'npts_0.5': [False, [63987, 63987, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63989, 63989, 63990, 63990, 63990, 63989, 63989, 63990, 63990, 63991, 63990, 63990, 63991, 63991, 63991, 63991, 63990, 63991, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63989, 63990, 63990, 63990, 63989, 63989, 63988, 63988, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63989, 63989, 63988, 63989, 63989, 63989, 63990, 63990, 63991, 63991, 63991, 63991, 63994, 63995, 63996, 63993, 63993, 63994, 63994, 63993, 63993, 63993, 63993, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63993, 63992, 63992, 63992, 63993, 63994, 63994, 63994, 63991, 63990, 63990, 63990, 63990, 63990, 63991, 63991, 63991, 63991, 63991, 63989, 63990, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63991, 63991, 63991, 63992, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63990, 63990, 63990, 63990, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63987, 63989, 63989, 63989, 63988, 63989, 63989, 63988, 63987, 63986, 63986, 63985, 63985, 63986, 63985, 63985, 63986, 63985, 63984, 63985, 63984, 63984, 63986, 63986, 63985, 63986, 63986, 63986, 63986, 63984, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63982, 63983, 63984, 63985, 63986, 63986, 63987, 63988, 63986, 63984, 63984, 63984, 63984, 63983, 63983, 63985, 63984, 63984, 63983, 63984, 63983, 63983, 63984, 63984, 63984, 63984, 63984, 63983, 63981, 63981, 63981, 63981, 63982, 63982, 63983, 63983, 63983, 63981, 63982, 63983, 63982, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63984, 63985, 63984, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63981, 63981, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63980, 63980, 63979, 63980, 63979, 63979, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63982, 63982, 63982, 63983, 63982, 63983, 63983, 63983, 63984, 63985, 63985, 63985, 63985, 63985, 63984, 63984, 63983, 63984, 63984, 63983, 63983, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63984, 63984, 63984, 63985, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63984, 63983, 63983, 63983, 63985, 63985, 63984, 63984, 63987, 63985, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63984, 63987, 63988, 63987, 63987, 63987, 63987, 63987, 63986, 63987, 63987, 63988, 63987, 63987, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63985, 63984, 63988, 63988, 63988, 63988, 63986, 63985, 63986, 63986, 63986, 63986, 63987, 63986, 63985, 63985, 63985, 63986, 63986, 63987, 63987, 63987, 63986, 63985, 63985, 63986, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63985, 63985, 63984, 63985, 63976, 63976, 63974, 63974, 63974, 63974, 63975, 63974, 63974, 63974, 63973, 63973, 63972, 63973, 63973, 63974, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63972, 63972, 63971, 63971, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63973, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63969, 63970, 63971, 63969, 63970, 63972, 63972, 63975, 63975, 63974, 63974, 63974, 63975, 63975, 63975, 63974, 63973, 63973, 63973, 63973, 63972, 63973, 63973, 63972, 63972, 63971, 63969, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63972, 63972, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63970, 63969, 63969, 63969, 63969, 63970, 63970, 63971, 63971, 63971, 63970, 63971, 63970, 63970, 63968, 63970, 63970, 63970, 63971, 63973, 63972, 63972, 63973, 63973, 63974, 63974, 63974, 63973, 63974, 63974, 63974, 63974, 63968, 63969, 63968, 63968, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63970, 63970, 63970, 63968, 63969, 63969, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63970, 63971, 63969, 63969, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63966, 63966, 63967, 63966, 63966, 63966, 63966, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63965, 63965, 63966, 63965, 63965, 63966, 63965, 63966, 63965, 63966, 63962, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63963, 63963, 63963, 63962, 63964, 63964, 63964, 63963, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63961, 63961, 63962, 63961, 63960, 63961, 63961, 63963, 63963, 63961, 63961, 63962, 63962, 63961, 63961, 63961, 63961, 63961, 63961, 63962, 63961, 63961, 63961, 63961, 63960, 63959, 63959, 63958, 63958, 63958, 63958, 63959, 63959, 63958, 63959, 63960, 63959, 63958, 63958, 63959, 63959, 63959, 63960, 63960, 63959, 63960, 63959, 63960, 63960, 63960, 63959, 63959, 63960, 63961, 63961, 63961, 63961, 63962, 63962, 63962, 63963, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63962, 63961, 63962, 63961, 63962, 63961, 63961, 63961, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63963, 63963, 63961, 63963, 63963, 63962, 63963, 63962, 63962, 63963, 63962, 63964, 63963, 63963, 63963, 63963, 63962, 63962, 63962, 63961, 63960, 63959, 63956, 63956, 63958, 63958, 63958, 63958, 63958, 63961, 63961, 63961, 63962, 63962, 63961, 63960, 63960, 63960, 63958, 63958, 63958, 63958, 63958, 63958, 63958, 63960, 63961, 63961, 63961, 63961, 63960, 63959, 63960, 63960, 63960, 63963, 63962, 63963, 63961, 63961, 63961, 63961, 63959, 63958, 63960, 63960, 63960, 63958, 63958, 63958, 63959, 63960, 63960, 63958, 63958, 63959, 63959, 63959, 63958, 63958, 63958, 63958, 63958, 63958, 63957, 63957, 63957, 63957, 63954, 63953, 63953, 63953, 63952, 63954, 63954, 63954, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63953, 63953, 63953, 63952, 63952, 63952, 63952, 63953, 63952, 63951, 63951, 63951, 63951, 63952, 63952, 63952, 63951, 63953, 63952, 63951, 63951, 63950, 63949, 63949, 63949, 63949, 63947, 63947, 63947, 63948, 63948, 63948, 63948, 63949, 63946, 63946, 63946, 63948, 63948, 63948]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1057529783407027, 36.994958712675974, \
                    36.71800149173757]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.8801035135706],
            'fit_pix': [False, [240.80157119155351, 209.98069221787847]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.044496592134237289],
            'min_val_pos': [True, [250, 204, 0, 13]],
            'im_rms': [False, 0.012118826675730165],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 203.0061296436524],
            'im_sum': [False, 193.48533807980857],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.8853663419051631, 1.0831964982105018, \
                      0.8486568935446293]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462000557837],
            'fit_pix_0': [False, [239.96619358406645, 209.9923811349359]],
            'fit_1': [False, [0.8855601515490675, 1.0833246653826742, \
                      0.8510940656955308]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.8801035135706],
            'fit_pix_1': [False, [239.96614158332466, 209.99207827032774]],
            'fit_2': [False, [0.8851418403402185, 1.0806791592573344, \
                      0.8488399641404538]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99558702156276],
            'fit_pix_2': [False, [239.96629893524417, 209.9920584854601]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 104998085.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.070553310215473175],
            'max_val_pos': [True, [304, 256, 0, 321]],
            'min_val': [False, -0.07086267322301865],
            'min_val_pos': [True, [222, 212, 0, 674]],
            'im_rms': [False, 0.012274356404197695],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 3864.4951676249366],
            'im_sum': [False, 3744.7025652566444],
            'regn_sum': [False, 37.53474615866071],
            'npts_real': [True, 191116800]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_eph)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.03622784838080406],
            'max_val_pos': [True, [315, 229, 0, 511]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            #CAS-9386 update build100 serial
            #'im_rms': [False, 2.4667777528231367e-05],
            #'im_sum': [False, 9.705396932782605],
            #'regn_sum': [False, 0.8115214343415573],
            'im_rms': [False, 2.413243705039696e-05],
            'im_sum': [False, 9.349006329197437], 
            'regn_sum': [False, 0.7349461655830964],
            'mask_non0': [True, 0],
            'npts_real': [True, 191116800]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, 
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube_eph)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [1e-10, 244151.1796875],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_cube_eph)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.3178490698337555],
            # CAS-9386 update build100 serial
            #'max_val_pos': [True, [240, 210, 0, 9]],
            #'min_val': [False, 6.582064816029742e-05],
            'max_val_pos': [True, [240, 209, 0, 944]],
            'min_val': [False, 7.376736175501719e-05],
            'im_rms': [False, 0.119822765111],
            'im_sum': [False, 13818688.5072],
            'npts_0.2': [False, [110768, 110768, 110768, 110767, 110766, 110766, 110765, 110765, 110765, 110765, 110766, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110768, 110768, 110768, 110766, 110766, 110768, 110767, 110768, 110768, 110767, 110767, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110766, 110770, 110770, 110768, 110768, 110770, 110771, 110771, 110771, 110770, 110770, 110770, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110769, 110769, 110769, 110768, 110768, 110768, 110769, 110771, 110771, 110771, 110770, 110770, 110770, 110770, 110770, 110770, 110769, 110767, 110769, 110768, 110769, 110770, 110772, 110772, 110772, 110772, 110772, 110772, 110773, 110773, 110773, 110774, 110774, 110773, 110771, 110772, 110771, 110770, 110770, 110770, 110769, 110771, 110771, 110769, 110770, 110769, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110769, 110769, 110769, 110769, 110770, 110770, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110770, 110769, 110769, 110768, 110768, 110768, 110769, 110768, 110767, 110767, 110768, 110767, 110768, 110767, 110767, 110767, 110768, 110769, 110769, 110768, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110766, 110766, 110766, 110767, 110769, 110769, 110770, 110771, 110770, 110769, 110770, 110770, 110769, 110768, 110768, 110769, 110769, 110769, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110770, 110770, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110771, 110771, 110771, 110773, 110774, 110773, 110774, 110776, 110775, 110775, 110773, 110773, 110774, 110774, 110774, 110775, 110774, 110774, 110776, 110776, 110776, 110777, 110776, 110776, 110776, 110776, 110776, 110776, 110775, 110773, 110773, 110773, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110775, 110774, 110774, 110775, 110775, 110775, 110776, 110777, 110777, 110777, 110777, 110778, 110778, 110779, 110778, 110780, 110781, 110781, 110781, 110779, 110780, 110780, 110781, 110779, 110779, 110779, 110779, 110779, 110782, 110782, 110781, 110781, 110782, 110782, 110782, 110782, 110782, 110782, 110781, 110783, 110782, 110782, 110782, 110783, 110782, 110782, 110783, 110780, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110782, 110782, 110781, 110781, 110781, 110781, 110781, 110781, 110781, 110782, 110781, 110781, 110784, 110782, 110782, 110784, 110782, 110782, 110782, 110782, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110783, 110783, 110782, 110783, 110783, 110784, 110784, 110784, 110786, 110786, 110785, 110786, 110786, 110786, 110786, 110785, 110785, 110784, 110785, 110785, 110785, 110785, 110785, 110785, 110786, 110787, 110787, 110787, 110787, 110786, 110786, 110786, 110786, 110786, 110787, 110787, 110787, 110787, 110787, 110786, 110787, 110787, 110786, 110786, 110785, 110785, 110786, 110786, 110786, 110786, 110786, 110784, 110786, 110786, 110786, 110786, 110785, 110786, 110786, 110785, 110785, 110782, 110783, 110785, 110785, 110786, 110785, 110785, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110783, 110784, 110783, 110783, 110783, 110784, 110784, 110784, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110746, 110746, 110746, 110745, 110746, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110748, 110748, 110748, 110748, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110747, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110747, 110745, 110745, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110745, 110745, 110744, 110745, 110744, 110744, 110744, 110745, 110745, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110743, 110743, 110743, 110743, 110744, 110741, 110743, 110743, 110743, 110740, 110741, 110740, 110741, 110743, 110743, 110743, 110741, 110741, 110741, 110741, 110741, 110740, 110741, 110740, 110740, 110741, 110741, 110740, 110741, 110741, 110742, 110742, 110742, 110742, 110741, 110742, 110742, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110741, 110741, 110741, 110742, 110742, 110741, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110742, 110742, 110742, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110742, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110744, 110744, 110744, 110742, 110741, 110741, 110742, 110742, 110740, 110740, 110741, 110741, 110741, 110740, 110740, 110740, 110739, 110740, 110739, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110739, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110736, 110738, 110736, 110738, 110737, 110738, 110738, 110738, 110738, 110739, 110739, 110736, 110736, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110738, 110737, 110737, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110738, 110737, 110737, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110739, 110739, 110738, 110738, 110738, 110736, 110736, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110734, 110733, 110734, 110733, 110735, 110735, 110735, 110735, 110736, 110737, 110736, 110736, 110736, 110736, 110735, 110736, 110736, 110733, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110733, 110734, 110734, 110734, 110733, 110733, 110734, 110734, 110734, 110732, 110734, 110733, 110733, 110731, 110732, 110731, 110732, 110733, 110733, 110733, 110731, 110731, 110732, 110730, 110731, 110733, 110732, 110732]],
            'npts_0.5': [False, [63987, 63987, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63989, 63989, 63990, 63990, 63990, 63989, 63989, 63990, 63990, 63991, 63990, 63990, 63991, 63991, 63991, 63991, 63990, 63991, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63989, 63990, 63990, 63990, 63989, 63989, 63988, 63988, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63989, 63989, 63988, 63989, 63989, 63989, 63990, 63990, 63991, 63991, 63991, 63991, 63994, 63995, 63996, 63993, 63993, 63994, 63994, 63993, 63993, 63993, 63993, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63993, 63992, 63992, 63992, 63993, 63994, 63994, 63994, 63991, 63990, 63990, 63990, 63990, 63990, 63991, 63991, 63991, 63991, 63991, 63989, 63990, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63991, 63991, 63991, 63992, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63990, 63990, 63990, 63990, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63987, 63989, 63989, 63989, 63988, 63989, 63989, 63988, 63987, 63986, 63986, 63985, 63985, 63986, 63985, 63985, 63986, 63985, 63984, 63985, 63984, 63984, 63986, 63986, 63985, 63986, 63986, 63986, 63986, 63984, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63982, 63983, 63984, 63985, 63986, 63986, 63987, 63988, 63986, 63984, 63984, 63984, 63984, 63983, 63983, 63985, 63984, 63984, 63983, 63984, 63983, 63983, 63984, 63984, 63984, 63984, 63984, 63983, 63981, 63981, 63981, 63981, 63982, 63982, 63983, 63983, 63983, 63981, 63982, 63983, 63982, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63984, 63985, 63984, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63981, 63981, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63980, 63980, 63979, 63980, 63979, 63979, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63982, 63982, 63982, 63983, 63982, 63983, 63983, 63983, 63984, 63985, 63985, 63985, 63985, 63985, 63984, 63984, 63983, 63984, 63984, 63983, 63983, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63984, 63984, 63984, 63985, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63984, 63983, 63983, 63983, 63985, 63985, 63984, 63984, 63987, 63985, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63984, 63987, 63988, 63987, 63987, 63987, 63987, 63987, 63986, 63987, 63987, 63988, 63987, 63987, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63985, 63984, 63988, 63988, 63988, 63988, 63986, 63985, 63986, 63986, 63986, 63986, 63987, 63986, 63985, 63985, 63985, 63986, 63986, 63987, 63987, 63987, 63986, 63985, 63985, 63986, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63985, 63985, 63984, 63985, 63976, 63976, 63974, 63974, 63974, 63974, 63975, 63974, 63974, 63974, 63973, 63973, 63972, 63973, 63973, 63974, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63972, 63972, 63971, 63971, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63973, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63969, 63970, 63971, 63969, 63970, 63972, 63972, 63975, 63975, 63974, 63974, 63974, 63975, 63975, 63975, 63974, 63973, 63973, 63973, 63973, 63972, 63973, 63973, 63972, 63972, 63971, 63969, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63972, 63972, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63970, 63969, 63969, 63969, 63969, 63970, 63970, 63971, 63971, 63971, 63970, 63971, 63970, 63970, 63968, 63970, 63970, 63970, 63971, 63973, 63972, 63972, 63973, 63973, 63974, 63974, 63974, 63973, 63974, 63974, 63974, 63974, 63968, 63969, 63968, 63968, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63970, 63970, 63970, 63968, 63969, 63969, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63970, 63971, 63969, 63969, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63966, 63966, 63967, 63966, 63966, 63966, 63966, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63965, 63965, 63966, 63965, 63965, 63966, 63965, 63966, 63965, 63966, 63962, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63963, 63963, 63963, 63962, 63964, 63964, 63964, 63963, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63961, 63961, 63962, 63961, 63960, 63961, 63961, 63963, 63963, 63961, 63961, 63962, 63962, 63961, 63961, 63961, 63961, 63961, 63961, 63962, 63961, 63961, 63961, 63961, 63960, 63959, 63959, 63958, 63958, 63958, 63958, 63959, 63959, 63958, 63959, 63960, 63959, 63958, 63958, 63959, 63959, 63959, 63960, 63960, 63959, 63960, 63959, 63960, 63960, 63960, 63959, 63959, 63960, 63961, 63961, 63961, 63961, 63962, 63962, 63962, 63963, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63962, 63961, 63962, 63961, 63962, 63961, 63961, 63961, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63963, 63963, 63961, 63963, 63963, 63962, 63963, 63962, 63962, 63963, 63962, 63964, 63963, 63963, 63963, 63963, 63962, 63962, 63962, 63961, 63960, 63959, 63956, 63956, 63958, 63958, 63958, 63958, 63958, 63961, 63961, 63961, 63962, 63962, 63961, 63960, 63960, 63960, 63958, 63958, 63958, 63958, 63958, 63958, 63958, 63960, 63961, 63961, 63961, 63961, 63960, 63959, 63960, 63960, 63960, 63963, 63962, 63963, 63961, 63961, 63961, 63961, 63959, 63958, 63960, 63960, 63960, 63958, 63958, 63958, 63959, 63960, 63960, 63958, 63958, 63959, 63959, 63959, 63958, 63958, 63958, 63958, 63958, 63958, 63957, 63957, 63957, 63957, 63954, 63953, 63953, 63953, 63952, 63954, 63954, 63954, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63953, 63953, 63953, 63952, 63952, 63952, 63952, 63953, 63952, 63951, 63951, 63951, 63951, 63952, 63952, 63952, 63951, 63953, 63952, 63951, 63951, 63950, 63949, 63949, 63949, 63949, 63947, 63947, 63947, 63948, 63948, 63948, 63948, 63949, 63946, 63946, 63946, 63948, 63948, 63948]],
            'npts_real': [True, 191116800]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict['test_mosaic_cube_eph']['self.parallel'] = self.parallel
        test_dict['test_mosaic_cube_eph']['report'] = report
        test_dict['test_mosaic_cube_eph']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.01, 0.1])
        self.mom8_creator(img+'.residual', range_list=[-0.01, 0.1])
        test_dict['test_mosaic_cube_eph']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_mosaic_cube_eph']['images'].append(img+'.image.profile.png')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of est_mosaic_cube_eph
#-------------------------------------------------#
    @stats_dict(test_dict)
    def test_mosaic_cube_eph_briggsbwtaper(self):
        ''' Mosaic ephemeris cube imaging with briggsbwtaper - field Venus, spw 45 '''

        test_name = self._testMethodName
        file_name = 'mosaic_cube_eph_briggsbwtaper.iter'
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
            perchanweightdensity=True, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=False, restoringbeam='common', \
            pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, \
            threshold='0.0mJy', nsigma=0.0, interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=2.0, noisethreshold=4.25, \
            lownoisethreshold=1.5, negativethreshold=15.0, minbeamfrac=0.3, \
            growiterations=50, dogrowprune=True, minpercentchange=1.0, \
            fastnoise=False, savemodel='none', parallel=self.parallel, verbose=True)

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
            perchanweightdensity=True, gridder='mosaic', chanchunks=-1, \
            mosweight=True, usepointing=False, pblimit=0.2, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggsbwtaper', robust=0.5, npixels=0, \
            niter=700000, threshold='0.0106Jy', nsigma=0.0, interactive=0, \
            usemask='auto-multithresh', sidelobethreshold=2.0, \
            noisethreshold=4.25, lownoisethreshold=1.5, \
            negativethreshold=15.0, minbeamfrac=0.3, growiterations=50, \
            dogrowprune=True, minpercentchange=1.0, fastnoise=False, \
            restart=True, savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))
  
        # test_mosaic_cube_eph_briggsbwtaper
        # .image report
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        # test_mosaic_cube_eph_briggsbwtaper.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 0.9347811147697555],
            'com_bmin': [False, 0.7185930247591006],
            'com_pa': [False, -88.16793752979638],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 104998085.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.086364768445491791],
            'max_val_pos': [True, [318, 243, 0, 511]],
            'min_val': [False, -0.070862673223018646,],
            'min_val_pos': [True, [222, 212, 0, 674]],
            'im_rms': [False, 0.012281882598118963],
            'rms_per_chan': [False, [0.01566755765592744, 0.012499771976912814, 0.011123979147164407, 0.011561886220609966, 0.012292568630946268, 0.013503738082575103, 0.013895709096621548, 0.015371602081183504, 0.012571966053954607, 0.012110657886905454, 0.010681959860783033, 0.011448528933751815, 0.013363157615506101, 0.01177851095735433, 0.011819999220451467, 0.011063128145379053, 0.011211424009495315, 0.012675826502373685, 0.00964365231187898, 0.01230416428699256, 0.011816948425841614, 0.012524610132431753, 0.012091583554529339, 0.01415517292636185, 0.013038519870781123, 0.010927180475068366, 0.011348376416636784, 0.011452376144559037, 0.013880464215687604, 0.01271421858369905, 0.012365348078626935, 0.011134439502267753, 0.01173999969262898, 0.010275262798856601, 0.012024637616711497, 0.011827779313317115, 0.012005719677964755, 0.01126998085149768, 0.01263100142352357, 0.010919404120788623, 0.01230292923375664, 0.013824000423944424, 0.012030441249337255, 0.01471235121764837, 0.012690250219699167, 0.012264732434797237, 0.010338113974845756, 0.010588074612506274, 0.010658916275604544, 0.01192531675206317, 0.011963523090655131, 0.011739554348224612, 0.010864559965038227, 0.01135082917289156, 0.011495290472285927, 0.011805507617290682, 0.010933831622378246, 0.014576785073645402, 0.01297711633387727, 0.010222489107261269, 0.0110984415724654, 0.011692109085082973, 0.011605118458967478, 0.011686408270450569, 0.010719138244995948, 0.011490010832745277, 0.012016309876685964, 0.012905508417480243, 0.010424679511378012, 0.010476580776329382, 0.013555509132609708, 0.01168612606268359, 0.01089950574737809, 0.009744535896916217, 0.011655046195014002, 0.012674198608575759, 0.014310946665162668, 0.014293294641270153, 0.012638247297143402, 0.011225312863067707, 0.010533579687603354, 0.010619978615716235, 0.01130230533155928, 0.012504854882493843, 0.011302020261357171, 0.012161798863089938, 0.010608102908705218, 0.011122188031292404, 0.010132879137328434, 0.011293949972102117, 0.012010374392817949, 0.013191568018758376, 0.01233135069689494, 0.010769746462901904, 0.01153549017814873, 0.011979399621333483, 0.011111746944971652, 0.011843988636800718, 0.011876601129763641, 0.00985833760297452, 0.012805113173417818, 0.012001866506178401, 0.012972732329926176, 0.012947046238817135, 0.01169227386437645, 0.013660834492575508, 0.015240297428641192, 0.012045506099515298, 0.014349189506965473, 0.011753680556779432, 0.011011732485486841, 0.011048193549742145, 0.012814236743542255, 0.011674542313321304, 0.010770073501250303, 0.012405362003451131, 0.011262904273749614, 0.010016645329827439, 0.012059229625014212, 0.01128781720269892, 0.011844913248897927, 0.012245345012809291, 0.009952101453716862, 0.012237226805610778, 0.010830186561964461, 0.012038505819334888, 0.010469212178041205, 0.012712147168754424, 0.009858813880899741, 0.009618511420097878, 0.011251630406509621, 0.01094736018021904, 0.009938750246980809, 0.011806500650862979, 0.013746618610294636, 0.011890591602572501, 0.011285353902870112, 0.010746148131018936, 0.010780721254182541, 0.010400367901096838, 0.010548644745226116, 0.011701519719661485, 0.011515444299289019, 0.01137979009080972, 0.009897580126157974, 0.012679387512909252, 0.012393749358995263, 0.010126541811740334, 0.012078757530494299, 0.010422834662425521, 0.011491611179879445, 0.011350002754119041, 0.010786499538647526, 0.01241069099315559, 0.01287288681929542, 0.012544796642462093, 0.012486480230120337, 0.012027122098965512, 0.011436187647136224, 0.014507824634016918, 0.012781326055012877, 0.010971035372093395, 0.011222652982616295, 0.010381402164033299, 0.011789211848634065, 0.010467460759016958, 0.010519913231344411, 0.009804241543959213, 0.009826802489597598, 0.011168128286266319, 0.009898755834256175, 0.010638546697493598, 0.01075904603518445, 0.01186190395476391, 0.012698873200736769, 0.012275066865331018, 0.010120998522019762, 0.011856430473088262, 0.01187640156331042, 0.012858360146453808, 0.010458450314412523, 0.011686631685348224, 0.012324510594898344, 0.013706507117441952, 0.010507073125858515, 0.009727999262257602, 0.010749737212657888, 0.010496316622257333, 0.012266025310015635, 0.009968887987260507, 0.010656267754294232, 0.01137826346655607, 0.013691020253483116, 0.009712315969605423, 0.0116765832200688, 0.01204200794054857, 0.011085824144059132, 0.012319945895455337, 0.010766015294335025, 0.012434827540397717, 0.012730633406142734, 0.014922472720576237, 0.013340994079356561, 0.012184668987628627, 0.011860251479407752, 0.010611453465531486, 0.010927075341466757, 0.010737943994644793, 0.011284903658799174, 0.012101560213749602, 0.01253847688484235, 0.011820390980633899, 0.01141402994073023, 0.01060976878040442, 0.010983370981323954, 0.013427371951491312, 0.012374725410784092, 0.011400378106504076, 0.011527234422091191, 0.01340962759016125, 0.013570332253931302, 0.011804029668992844, 0.012197974733541946, 0.013983623168953468, 0.012725763762714332, 0.011354745750571254, 0.012145679916283425, 0.01126874546535086, 0.012261492463656863, 0.011178855897759667, 0.010264367951585555, 0.013345043592089396, 0.012409363114615287, 0.011495934975609993, 0.009277679225470462, 0.011720945646776957, 0.012561715975156802, 0.012325900140530158, 0.011173759164206036, 0.01444252189374611, 0.013235686426999993, 0.009545194033680014, 0.00989282465604266, 0.011317913505981915, 0.011172814496764794, 0.011097429969483154, 0.011431240017029451, 0.010546755637817182, 0.011214435474509551, 0.011763583664075946, 0.014689109193649632, 0.011251056734311256, 0.009776300805049154, 0.012860933158710076, 0.012509654758086465, 0.013307995464813777, 0.012069070120183812, 0.011422201031864795, 0.011512246169509426, 0.011243380202117025, 0.011621384685604818, 0.011008804779668683, 0.013029468912847373, 0.011484166233438071, 0.01154715840745958, 0.012030983875891486, 0.012644418892294014, 0.010840072063581964, 0.01217667528470696, 0.012348270662399291, 0.011295167588599378, 0.011131983820484052, 0.011525183914729311, 0.01245079346986009, 0.010668667918574357, 0.011110716151104136, 0.012536736404100099, 0.013297691638473148, 0.013088181650764244, 0.011822130327201136, 0.010774771031200035, 0.010775807119122434, 0.01105586988033813, 0.009420738313971375, 0.012013742695805617, 0.011326747652322514, 0.012096904239732881, 0.01150806457510923, 0.011675274227450103, 0.0135991855834598, 0.01044590308467654, 0.01066146998192241, 0.010712518291299183, 0.01193296299641928, 0.012947758061995076, 0.010130158406514036, 0.011855179848796187, 0.011440889779476115, 0.012458460329640899, 0.011076297729893039, 0.012814778866789587, 0.011811606075871006, 0.012867705357363268, 0.01621962105427537, 0.016465436698466954, 0.016201295070531984, 0.01058752407190847, 0.011716930684685505, 0.011008111220908223, 0.01309686038666478, 0.011249739924392245, 0.011298808361717275, 0.01286619495044334, 0.012457597012475598, 0.011357408247558936, 0.011635136666110633, 0.01147265786957545, 0.01411270947674405, 0.012138457217748517, 0.012013161207475307, 0.012399897890804487, 0.017261754698869394, 0.011751036356983229, 0.01220773352984818, 0.01169628539620397, 0.013572907681349047, 0.014279084224693019, 0.013350071117589353, 0.01119955478052971, 0.011678178134442125, 0.01143788295713284, 0.01338236323998995, 0.013148991760974126, 0.012261293827211488, 0.012230033568364952, 0.011982727185002299, 0.014628741472584738, 0.013936276355000195, 0.011606244163871274, 0.013217927391198507, 0.011877296698841999, 0.013751323956061321, 0.010502431019977986, 0.011623558011397342, 0.013225657353801752, 0.013687394057886515, 0.012047896602847821, 0.01247678361440421, 0.010956367674752512, 0.011928249343126845, 0.012456739592702704, 0.014412803794511948, 0.013091054766212506, 0.012391320075065916, 0.012041284843754036, 0.011175660150719864, 0.012796161329318703, 0.012376600001200107, 0.013213434092613835, 0.011391999754128477, 0.011735911173241135, 0.011338817659054366, 0.013685494524633305, 0.01405734652024605, 0.011873232352346962, 0.011492726398798213, 0.011357899840429751, 0.011924626623408293, 0.01292314975466112, 0.011568890237960018, 0.012450768675985815, 0.01615561678574021, 0.014389883701947508, 0.012020341777831386, 0.011691722848726802, 0.012565949873863497, 0.011694804496294954, 0.011846565903593035, 0.011247492738346963, 0.010809808599760048, 0.011330103597179153, 0.012046098499893751, 0.012792309821515609, 0.012055452071571462, 0.012529809308170823, 0.01166776384371269, 0.01271547225436619, 0.014155919212706156, 0.011590774944224877, 0.013272175555436434, 0.01230731481769422, 0.01134847189287814, 0.011467954772903752, 0.01184115575454566, 0.011988620387346051, 0.010501901957914332, 0.011616509715974813, 0.011471135705227897, 0.012451191281457752, 0.015319504006822212, 0.012649022813805498, 0.011475752166486771, 0.012390712505079823, 0.013092398166506612, 0.0108425779557803, 0.010850393792933026, 0.011858723037825333, 0.014367570168849372, 0.016151613939500507, 0.010857278233794036, 0.012348584717685072, 0.013369893472750983, 0.013520298224929771, 0.012490511233733568, 0.011838067139252433, 0.01379337132839034, 0.013615981503194808, 0.01213329972068078, 0.011926986383721444, 0.011989638160863852, 0.011459946102237622, 0.00998617556766575, 0.011909475646001892, 0.01144107317570622, 0.012796872801637271, 0.011812967842105799, 0.01226189488300729, 0.013800981007306753, 0.011403886599520798, 0.011873516666846865, 0.012409987676322875, 0.012188938970018453, 0.011353880208547378, 0.013058666569184637, 0.011338981854502068, 0.01409796019586988, 0.012231139093670572, 0.015107895235161905, 0.012326810087808738, 0.012946694516644243, 0.012447550801307723, 0.014909071523596254, 0.01397117897114995, 0.013954030858987044, 0.012642975889461814, 0.012387518693283873, 0.011610200546359207, 0.011526268526938685, 0.011758827798373622, 0.011908680699262349, 0.015003488394504251, 0.014659137370125483, 0.01303276571543743, 0.013752408592063675, 0.013708881194186243, 0.01147231577101178, 0.011546435459926428, 0.011541811993355754, 0.013315481649931753, 0.012924077797838931, 0.013861771686590075, 0.012908741618203564, 0.012185556034842455, 0.0123901699294395, 0.011089612618787055, 0.012565991045585109, 0.011172521028748187, 0.011930136297848918, 0.012217995900935783, 0.012567934765036304, 0.013402949630010139, 0.012677891119972322, 0.012831082643587562, 0.012077587392213665, 0.01192202179044923, 0.012437303467243326, 0.012739158128205348, 0.013229922671982922, 0.011707021012177004, 0.013242901791918983, 0.012429029653969133, 0.011897651409429517, 0.012333981656159804, 0.01652056933731621, 0.012803033285199423, 0.013291710033149702, 0.013915563269442936, 0.013293192094911472, 0.012809923485290182, 0.012329480430228985, 0.012327561250559955, 0.011635712112230329, 0.013309364248632287, 0.01298532449779078, 0.011141245068791102, 0.0123662111290416, 0.012952937497979091, 0.012652571307871136, 0.015633401429410065, 0.014029079448528728, 0.013085527098929505, 0.013671367699869673, 0.011449641495798529, 0.013870108540867539, 0.013267219241966837, 0.013009671745599366, 0.012061924878367703, 0.013344180204315093, 0.012472386713589599, 0.012819448624406097, 0.013468487050066785, 0.015061287987425221, 0.013101479251036413, 0.011695612693457496, 0.014215627129265261, 0.014030414081014604, 0.013144534677081394, 0.013769631498326561, 0.013665079857862566, 0.01278869825702419, 0.011758045653360501, 0.013379677594574436, 0.013176320438872265, 0.013726863311133137, 0.014093922857692474, 0.011978594418837533, 0.012688741742419463, 0.014213190054434053, 0.012778188415525355, 0.012893566603547115, 0.011031522915526809, 0.012099967208788644, 0.01197488835353761, 0.011886891227553936, 0.012806477665845199, 0.013276830507854871, 0.013241164133079063, 0.012154967597401761, 0.011697698448488806, 0.0138159512574303, 0.012144528154314117, 0.011982896302520097, 0.012088015079846063, 0.01258443418498413, 0.012417041503252886, 0.013316302302718812, 0.013253698599679495, 0.014875163789523499, 0.013611399587966033, 0.01270670566020336, 0.012428452848201328, 0.013258866017540937, 0.012735738693151316, 0.014906643094604965, 0.014207002883168778, 0.01351154568809692, 0.012428703834641097, 0.012531695617773475, 0.014339355016753071, 0.011835752014891681, 0.015198065031029119, 0.01430821328282015, 0.013257734215508951, 0.012994339259791049, 0.013865856665286771, 0.012158304848082375, 0.014198843802734912, 0.01419354016501019, 0.013678418586819103, 0.012056908737372238, 0.012088134020138442, 0.012708550023920906, 0.012056042831971507, 0.012944165533598666, 0.012998439619323971, 0.011738401817877483, 0.013534728428647752, 0.012933197018222002, 0.011888680423975263, 0.012833718856536632, 0.013026329306504825, 0.01389771764379911, 0.013170374270730664, 0.01394810209074441, 0.013726597650998646, 0.012620621703300312, 0.012195206631413907, 0.013644885542175174, 0.01339675796598244, 0.015953133266745596, 0.015572165751410377, 0.013918610648959726, 0.012539523963263447, 0.014863648164363799, 0.013557760935362215, 0.013767615205969776, 0.012387780847253898, 0.013111551187613596, 0.012178574316034968, 0.012804589291833082, 0.013764797217935443, 0.014682794411092211, 0.011923563140996867, 0.012457566642369073, 0.012577155733152602, 0.012268408936443286, 0.013433767866486363, 0.011254144886535153, 0.012239416835082412, 0.012719204139179307, 0.014153628502136058, 0.013714206504685713, 0.01555664321990409, 0.012482696827891392, 0.012798094219294787, 0.013725222959993684, 0.011297342405230785, 0.013866245776811888, 0.013948944319164524, 0.014752096346863636, 0.012392843752160694, 0.014416608533154353, 0.013416960727895713, 0.013014075062578258, 0.012457948814251782, 0.01214561019554649, 0.014573596703768801, 0.013355931500995009, 0.010935442745219525, 0.013474856369494768, 0.014093982151433494, 0.011282427815945854, 0.013022998498462857, 0.011932228038133832, 0.011193973084211214, 0.013372017260001274, 0.013994297538831645, 0.011742113041326014, 0.012987909251690535, 0.012569010525416004, 0.011966315747010852, 0.013958876937310783, 0.013110867016705616, 0.015176327530987377, 0.013271913419261386, 0.014815021510858065, 0.013019719016438232, 0.012593706022601611, 0.013636464785719227, 0.013024015052417843, 0.011867287358548854, 0.014432494354300884, 0.011513467360750092, 0.013751315885178503, 0.011944914326008424, 0.011579770223979212, 0.013104138863451702, 0.012572035675858769, 0.013719589991898679, 0.011799423913740791, 0.014222826331815612, 0.013130190236718657, 0.015687333785054065, 0.013977869745415655, 0.013856666637965972, 0.0119206019645679, 0.01340210329202167, 0.01236780694175228, 0.013548894032861069, 0.012186169719114176, 0.01048406598287907, 0.010655269577945959, 0.011891505227937497, 0.012232046724993523, 0.016338504257622802, 0.011951257476223706, 0.012787444872754723, 0.011710212821400005, 0.011695865895061741, 0.015838707904307876, 0.012939755683705046, 0.012477587036156625, 0.013330358565129482, 0.01082528855224014, 0.013773287643131382, 0.012570171555245143, 0.013162359160576995, 0.011218965550597317, 0.01342890868237659, 0.01202182099406552, 0.010697138162887927, 0.015232948941257044, 0.012711001004850505, 0.011506289930999608, 0.012369850055536856, 0.01233361532935207, 0.012360627932656498, 0.013326575850585057, 0.012934710498929827, 0.013248039695274097, 0.012226203594835677, 0.012818113815664448, 0.01341051851094283, 0.013541762993268074, 0.013056711886632514, 0.012115662158377625, 0.013357590759015063, 0.012032645498305006, 0.016207947768265084, 0.01389311652405909, 0.013847679128406589, 0.013886386730554553, 0.012067502182802831, 0.011753272427404053, 0.013279048057631762, 0.011214754296383541, 0.011683816996863522, 0.010609551860508486, 0.012357120345452005, 0.013050369382881564, 0.01174704009897493, 0.011743154945937362, 0.012531199758011487, 0.011798725028591634, 0.012419735337445137, 0.013975728513221264, 0.011781109647231187, 0.014598054146810256, 0.011754504071338916, 0.011305854302102692, 0.013318945027316854, 0.013042196249451354, 0.011052738741563284, 0.010558225013410463, 0.011020779774382213, 0.011340319330782552, 0.012033592610864789, 0.012455257536242884, 0.01110139391512171, 0.01046562816674095, 0.015266557792594495, 0.011918313536064278, 0.011360562910751727, 0.010365874108060247, 0.011737890229192168, 0.012793643279848948, 0.013256120317317155, 0.012383560012815365, 0.013301087661143126, 0.01531080427533124, 0.010224017979386747, 0.011257686276499953, 0.009965303526327884, 0.011481057685759922, 0.01274428412077991, 0.01237779048897638, 0.011054791139194602, 0.012765273103949742, 0.012472945953336118, 0.012608159036016397, 0.011204550231481033, 0.013905107233033231, 0.010773937005630393, 0.013310027013851727, 0.017420680889573912, 0.011673620547153839, 0.010614208861963346, 0.012253262631090077, 0.012560836225969919, 0.012560372974566972, 0.013315953743992456, 0.01250495576610427, 0.010300621736575054, 0.0126743669068127, 0.012497189500451297, 0.01058986672748917, 0.0140554647719132, 0.011586323266949704, 0.01051084444510618, 0.010723172131198219, 0.013670881790176151, 0.01294181100319149, 0.01329464895892603, 0.011257304758166548, 0.011406076857722271, 0.01192400994821658, 0.01090031653434306, 0.011516503202529208, 0.01164722421363161, 0.011278034900260045, 0.013101497734367626, 0.010571820551740189, 0.01258494483689583, 0.010690018512359867, 0.012288203401799393, 0.011394189183191222, 0.010661795126849621, 0.014600444436046287, 0.012545075612095091, 0.011433423838627848, 0.0108797962059277, 0.01179905039017321, 0.011201593386431653, 0.012398023234344002, 0.011585602038687984, 0.012105760854189395, 0.01012363568557298, 0.010608758997456596, 0.013055946950589595, 0.01115366731042303, 0.011356175896159434, 0.010664413285359914, 0.011095187626260742, 0.01055909008819056, 0.01199264964050027, 0.01198885985504728, 0.011856295143569895, 0.011865458906234921, 0.012803395900220201, 0.011576869098809376, 0.010189936253659717, 0.010195412935336356, 0.012280457502570626, 0.011418747542118402, 0.011054267711398783, 0.010663172535730193, 0.010788394178560413, 0.012825843175935922, 0.014830488452138243, 0.012284354036472657, 0.012694714212578187, 0.013724353996911465, 0.012851827129972288, 0.010592820880735418, 0.012302591864876384, 0.012065809850599015, 0.011956358889472574, 0.010872197413464683, 0.012323221424305216, 0.012123107304017921, 0.010815605182278264, 0.012638711922615793, 0.011655868729735695, 0.013318117884136317, 0.011374335257623935, 0.011336547412403876, 0.013538061509343441, 0.010884298085465788, 0.00983147858074236, 0.010395424381769653, 0.010546991873577608, 0.012797725536117008, 0.011287668516950705, 0.010972896319568352, 0.011817914808157692, 0.010995529479193617, 0.01009639486723553, 0.01141833937111661, 0.011690435639595306, 0.010061766319288705, 0.011971105561760871, 0.012700043131615801, 0.011318883390568383, 0.010909315349852766, 0.00955009744793359, 0.010498572864691814, 0.013751363636934978, 0.012414462975548525, 0.013914888435336918, 0.014227696895603312, 0.012244996318125041, 0.01179968071369893, 0.01269311294414022, 0.011707870106088975, 0.010958813990455807, 0.011651043993684336, 0.01030137511322718, 0.011034955595029074, 0.011010590344592711, 0.012071622334637981, 0.012156945636578436, 0.01207776327788143, 0.012467860244637925, 0.011202073492240526, 0.013522425017327757, 0.013101403096011953, 0.011576541797128533, 0.010161447675978133, 0.010808172743640786, 0.01053092595610692, 0.011264161352788734, 0.01233589290036475, 0.010814458779802489, 0.012333845855799757, 0.009642890958343297, 0.012173492497045947, 0.013243881943014759, 0.011496651178443193, 0.01015960484422759, 0.01181458373889043, 0.009783393881506904, 0.00995620066577112, 0.011522308540852155, 0.010744486398458035, 0.010880760248384722, 0.0103126533049451, 0.010249354156011579, 0.011065664705266258, 0.01142653024289266, 0.01119440735004512, 0.01097503503562323, 0.014023046074535458, 0.010432042532047749, 0.010194842729529698, 0.010817754756984273, 0.011134926202450978, 0.011429007186568798, 0.010583381706566801, 0.015263726359132578, 0.011921804488613408, 0.009753228291247933, 0.012196147705148373, 0.009316733907216289, 0.01146380450784002, 0.010523188985348219, 0.012950242966562036, 0.010905274766360632, 0.009936978616187102, 0.011966920815569286, 0.010385478642918732, 0.010089607608523066, 0.01317787233691099, 0.01482152961916169, 0.011166319039339316, 0.01158412070404005, 0.010686371837829486, 0.011619838402650143, 0.010062419209719961, 0.010294116467679117, 0.01054194871478958, 0.011768968106377857, 0.012003013751281379, 0.011635692550686211, 0.012427950542111656, 0.012436753739762036, 0.010597044719259695, 0.01047844931079727, 0.01300316381280212]],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 4241.3850895847481],
            'im_sum': [False, 4107.755953973777],
            'regn_sum': [False, 34.376076006740732],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.013727427598497976, 0.01265603739684268, 0.012647516154896329, 0.013363016802334838, 0.012374606023989609, 0.013129803882930785, 0.01267018256172788]],
            'profile': [False, 0.025674302295194625]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [240, 393, 0, 0]), \
                      # CAS-9386 update build100 serial
                      #(img+'.image', False, [240, 394, 0, 0]), \
                      (img+'.image', True, [240, 394, 0, 0]), \
                      (img+'.image', True, [49, 209, 0, 0]), \
                      (img+'.image', False, [48, 209, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_cube_eph_briggsbwtaper)
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_mosaic_cube_eph_briggsbwtaper.exp_mask_stats
        exp_mask_stats = {'npts': [True, 191116800],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            # CAS-9386 update build100 serial
            #'mask_pix': [False, 8519],
            'mask_pix': [False, 8194],
            'mask_regns': [True, 31],
            'npts_real': [True, 191116800]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 104998085.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.200000017881],
            'im_rms': [False, 0.638163448188254],
            'npts_0.2': [False, [110768, 110768, 110768, 110767, 110766, 110766, 110765, 110765, 110765, 110765, 110766, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110768, 110768, 110768, 110766, 110766, 110768, 110767, 110768, 110768, 110767, 110767, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110766, 110770, 110770, 110768, 110768, 110770, 110771, 110771, 110771, 110770, 110770, 110770, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110769, 110769, 110769, 110768, 110768, 110768, 110769, 110771, 110771, 110771, 110770, 110770, 110770, 110770, 110770, 110770, 110769, 110767, 110769, 110768, 110769, 110770, 110772, 110772, 110772, 110772, 110772, 110772, 110773, 110773, 110773, 110774, 110774, 110773, 110771, 110772, 110771, 110770, 110770, 110770, 110769, 110771, 110771, 110769, 110770, 110769, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110769, 110769, 110769, 110769, 110770, 110770, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110770, 110769, 110769, 110768, 110768, 110768, 110769, 110768, 110767, 110767, 110768, 110767, 110768, 110767, 110767, 110767, 110768, 110769, 110769, 110768, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110766, 110766, 110766, 110767, 110769, 110769, 110770, 110771, 110770, 110769, 110770, 110770, 110769, 110768, 110768, 110769, 110769, 110769, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110770, 110770, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110771, 110771, 110771, 110773, 110774, 110773, 110774, 110776, 110775, 110775, 110773, 110773, 110774, 110774, 110774, 110775, 110774, 110774, 110776, 110776, 110776, 110777, 110776, 110776, 110776, 110776, 110776, 110776, 110775, 110773, 110773, 110773, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110775, 110774, 110774, 110775, 110775, 110775, 110776, 110777, 110777, 110777, 110777, 110778, 110778, 110779, 110778, 110780, 110781, 110781, 110781, 110779, 110780, 110780, 110781, 110779, 110779, 110779, 110779, 110779, 110782, 110782, 110781, 110781, 110782, 110782, 110782, 110782, 110782, 110782, 110781, 110783, 110782, 110782, 110782, 110783, 110782, 110782, 110783, 110780, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110782, 110782, 110781, 110781, 110781, 110781, 110781, 110781, 110781, 110782, 110781, 110781, 110784, 110782, 110782, 110784, 110782, 110782, 110782, 110782, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110783, 110783, 110782, 110783, 110783, 110784, 110784, 110784, 110786, 110786, 110785, 110786, 110786, 110786, 110786, 110785, 110785, 110784, 110785, 110785, 110785, 110785, 110785, 110785, 110786, 110787, 110787, 110787, 110787, 110786, 110786, 110786, 110786, 110786, 110787, 110787, 110787, 110787, 110787, 110786, 110787, 110787, 110786, 110786, 110785, 110785, 110786, 110786, 110786, 110786, 110786, 110784, 110786, 110786, 110786, 110786, 110785, 110786, 110786, 110785, 110785, 110782, 110783, 110785, 110785, 110786, 110785, 110785, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110783, 110784, 110783, 110783, 110783, 110784, 110784, 110784, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110746, 110746, 110746, 110745, 110746, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110748, 110748, 110748, 110748, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110747, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110747, 110745, 110745, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110745, 110745, 110744, 110745, 110744, 110744, 110744, 110745, 110745, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110743, 110743, 110743, 110743, 110744, 110741, 110743, 110743, 110743, 110740, 110741, 110740, 110741, 110743, 110743, 110743, 110741, 110741, 110741, 110741, 110741, 110740, 110741, 110740, 110740, 110741, 110741, 110740, 110741, 110741, 110742, 110742, 110742, 110742, 110741, 110742, 110742, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110741, 110741, 110741, 110742, 110742, 110741, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110742, 110742, 110742, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110742, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110744, 110744, 110744, 110742, 110741, 110741, 110742, 110742, 110740, 110740, 110741, 110741, 110741, 110740, 110740, 110740, 110739, 110740, 110739, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110739, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110736, 110738, 110736, 110738, 110737, 110738, 110738, 110738, 110738, 110739, 110739, 110736, 110736, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110738, 110737, 110737, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110738, 110737, 110737, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110739, 110739, 110738, 110738, 110738, 110736, 110736, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110734, 110733, 110734, 110733, 110735, 110735, 110735, 110735, 110736, 110737, 110736, 110736, 110736, 110736, 110735, 110736, 110736, 110733, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110733, 110734, 110734, 110734, 110733, 110733, 110734, 110734, 110734, 110732, 110734, 110733, 110733, 110731, 110732, 110731, 110732, 110733, 110733, 110733, 110731, 110731, 110732, 110730, 110731, 110733, 110732, 110732]],
            'npts_0.5': [False, [63987, 63987, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63989, 63989, 63990, 63990, 63990, 63989, 63989, 63990, 63990, 63991, 63990, 63990, 63991, 63991, 63991, 63991, 63990, 63991, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63989, 63990, 63990, 63990, 63989, 63989, 63988, 63988, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63989, 63989, 63988, 63989, 63989, 63989, 63990, 63990, 63991, 63991, 63991, 63991, 63994, 63995, 63996, 63993, 63993, 63994, 63994, 63993, 63993, 63993, 63993, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63993, 63992, 63992, 63992, 63993, 63994, 63994, 63994, 63991, 63990, 63990, 63990, 63990, 63990, 63991, 63991, 63991, 63991, 63991, 63989, 63990, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63991, 63991, 63991, 63992, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63990, 63990, 63990, 63990, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63987, 63989, 63989, 63989, 63988, 63989, 63989, 63988, 63987, 63986, 63986, 63985, 63985, 63986, 63985, 63985, 63986, 63985, 63984, 63985, 63984, 63984, 63986, 63986, 63985, 63986, 63986, 63986, 63986, 63984, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63982, 63983, 63984, 63985, 63986, 63986, 63987, 63988, 63986, 63984, 63984, 63984, 63984, 63983, 63983, 63985, 63984, 63984, 63983, 63984, 63983, 63983, 63984, 63984, 63984, 63984, 63984, 63983, 63981, 63981, 63981, 63981, 63982, 63982, 63983, 63983, 63983, 63981, 63982, 63983, 63982, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63984, 63985, 63984, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63981, 63981, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63980, 63980, 63979, 63980, 63979, 63979, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63982, 63982, 63982, 63983, 63982, 63983, 63983, 63983, 63984, 63985, 63985, 63985, 63985, 63985, 63984, 63984, 63983, 63984, 63984, 63983, 63983, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63984, 63984, 63984, 63985, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63984, 63983, 63983, 63983, 63985, 63985, 63984, 63984, 63987, 63985, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63984, 63987, 63988, 63987, 63987, 63987, 63987, 63987, 63986, 63987, 63987, 63988, 63987, 63987, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63985, 63984, 63988, 63988, 63988, 63988, 63986, 63985, 63986, 63986, 63986, 63986, 63987, 63986, 63985, 63985, 63985, 63986, 63986, 63987, 63987, 63987, 63986, 63985, 63985, 63986, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63985, 63985, 63984, 63985, 63976, 63976, 63974, 63974, 63974, 63974, 63975, 63974, 63974, 63974, 63973, 63973, 63972, 63973, 63973, 63974, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63972, 63972, 63971, 63971, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63973, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63969, 63970, 63971, 63969, 63970, 63972, 63972, 63975, 63975, 63974, 63974, 63974, 63975, 63975, 63975, 63974, 63973, 63973, 63973, 63973, 63972, 63973, 63973, 63972, 63972, 63971, 63969, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63972, 63972, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63970, 63969, 63969, 63969, 63969, 63970, 63970, 63971, 63971, 63971, 63970, 63971, 63970, 63970, 63968, 63970, 63970, 63970, 63971, 63973, 63972, 63972, 63973, 63973, 63974, 63974, 63974, 63973, 63974, 63974, 63974, 63974, 63968, 63969, 63968, 63968, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63970, 63970, 63970, 63968, 63969, 63969, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63970, 63971, 63969, 63969, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63966, 63966, 63967, 63966, 63966, 63966, 63966, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63965, 63965, 63966, 63965, 63965, 63966, 63965, 63966, 63965, 63966, 63962, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63963, 63963, 63963, 63962, 63964, 63964, 63964, 63963, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63961, 63961, 63962, 63961, 63960, 63961, 63961, 63963, 63963, 63961, 63961, 63962, 63962, 63961, 63961, 63961, 63961, 63961, 63961, 63962, 63961, 63961, 63961, 63961, 63960, 63959, 63959, 63958, 63958, 63958, 63958, 63959, 63959, 63958, 63959, 63960, 63959, 63958, 63958, 63959, 63959, 63959, 63960, 63960, 63959, 63960, 63959, 63960, 63960, 63960, 63959, 63959, 63960, 63961, 63961, 63961, 63961, 63962, 63962, 63962, 63963, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63962, 63961, 63962, 63961, 63962, 63961, 63961, 63961, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63963, 63963, 63961, 63963, 63963, 63962, 63963, 63962, 63962, 63963, 63962, 63964, 63963, 63963, 63963, 63963, 63962, 63962, 63962, 63961, 63960, 63959, 63956, 63956, 63958, 63958, 63958, 63958, 63958, 63961, 63961, 63961, 63962, 63962, 63961, 63960, 63960, 63960, 63958, 63958, 63958, 63958, 63958, 63958, 63958, 63960, 63961, 63961, 63961, 63961, 63960, 63959, 63960, 63960, 63960, 63963, 63962, 63963, 63961, 63961, 63961, 63961, 63959, 63958, 63960, 63960, 63960, 63958, 63958, 63958, 63959, 63960, 63960, 63958, 63958, 63959, 63959, 63959, 63958, 63958, 63958, 63958, 63958, 63958, 63957, 63957, 63957, 63957, 63954, 63953, 63953, 63953, 63952, 63954, 63954, 63954, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63953, 63953, 63953, 63952, 63952, 63952, 63952, 63953, 63952, 63951, 63951, 63951, 63951, 63952, 63952, 63952, 63951, 63953, 63952, 63951, 63951, 63950, 63949, 63949, 63949, 63949, 63947, 63947, 63947, 63948, 63948, 63948, 63948, 63949, 63946, 63946, 63946, 63948, 63948, 63948]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1057529783407027, 36.994958712675974, \
                    36.71800149173757]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.8801035135706],
            'fit_pix': [False, [240.80157119155351, 209.98069221787847]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.044496592134237289],
            'min_val_pos': [True, [250, 204, 0, 13]],
            'im_rms': [False, 0.012118826675730165],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 203.0061296436524],
            'im_sum': [False, 193.48533807980857],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.8853663419051631, 1.0831964982105018, \
                      0.8486568935446293]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462000557837],
            'fit_pix_0': [False, [239.96619358406645, 209.9923811349359]],
            'fit_1': [False, [0.8855601515490675, 1.0833246653826742, \
                      0.8510940656955308]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.8801035135706],
            'fit_pix_1': [False, [239.96614158332466, 209.99207827032774]],
            'fit_2': [False, [0.8851418403402185, 1.0806791592573344, \
                      0.8488399641404538]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99558702156276],
            'fit_pix_2': [False, [239.96629893524417, 209.9920584854601]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 104998085.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.070553310215473175],
            'max_val_pos': [True, [304, 256, 0, 321]],
            'min_val': [False, -0.07086267322301865],
            'min_val_pos': [True, [222, 212, 0, 674]],
            'im_rms': [False, 0.012274356404197695],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 3864.4951676249366],
            'im_sum': [False, 3744.7025652566444],
            'regn_sum': [False, 37.53474615866071],
            'npts_real': [True, 191116800]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_eph_briggsbwtaper)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_eph_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.03622784838080406],
            'max_val_pos': [True, [315, 229, 0, 511]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            #CAS-9386 update build100 serial
            #'im_rms': [False, 2.4667777528231367e-05],
            #'im_sum': [False, 9.705396932782605],
            #'regn_sum': [False, 0.8115214343415573],
            'im_rms': [False, 2.413243705039696e-05],
            'im_sum': [False, 9.349006329197437], 
            'regn_sum': [False, 0.7349461655830964],
            'mask_non0': [True, 0],
            'npts_real': [True, 191116800]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, 
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube_eph_briggsbwtaper)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_eph_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [1e-10, 244151.1796875],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_cube_eph_briggsbwtaper)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        # test_mosaic_cube_eph_briggsbwtaper.exp_wt_stats
        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 2.617644e+11],
            'end': [True, 2.619956e+11],
            'start_delta': [False, 2.617644e+11],
            'end_delta': [False, 2.619956e+11],
            'nchan': [True, 948],
            'max_val': [False, 0.3178490698337555],
            # CAS-9386 update build100 serial
            #'max_val_pos': [True, [240, 210, 0, 9]],
            #'min_val': [False, 6.582064816029742e-05],
            'max_val_pos': [True, [240, 209, 0, 944]],
            'min_val': [False, 7.376736175501719e-05],
            'im_rms': [False, 0.119822765111],
            'im_sum': [False, 13818688.5072],
            'npts_0.2': [False, [110768, 110768, 110768, 110767, 110766, 110766, 110765, 110765, 110765, 110765, 110766, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110768, 110768, 110768, 110766, 110766, 110768, 110767, 110768, 110768, 110767, 110767, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110768, 110766, 110770, 110770, 110768, 110768, 110770, 110771, 110771, 110771, 110770, 110770, 110770, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110768, 110769, 110769, 110769, 110768, 110768, 110768, 110769, 110771, 110771, 110771, 110770, 110770, 110770, 110770, 110770, 110770, 110769, 110767, 110769, 110768, 110769, 110770, 110772, 110772, 110772, 110772, 110772, 110772, 110773, 110773, 110773, 110774, 110774, 110773, 110771, 110772, 110771, 110770, 110770, 110770, 110769, 110771, 110771, 110769, 110770, 110769, 110768, 110769, 110769, 110769, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110772, 110771, 110771, 110769, 110769, 110769, 110769, 110770, 110770, 110769, 110769, 110769, 110769, 110770, 110770, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110770, 110769, 110769, 110768, 110768, 110768, 110769, 110768, 110767, 110767, 110768, 110767, 110768, 110767, 110767, 110767, 110768, 110769, 110769, 110768, 110768, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110767, 110766, 110766, 110766, 110767, 110769, 110769, 110770, 110771, 110770, 110769, 110770, 110770, 110769, 110768, 110768, 110769, 110769, 110769, 110770, 110770, 110770, 110771, 110771, 110770, 110770, 110771, 110771, 110770, 110771, 110771, 110771, 110771, 110770, 110770, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110771, 110772, 110771, 110771, 110771, 110773, 110774, 110773, 110774, 110776, 110775, 110775, 110773, 110773, 110774, 110774, 110774, 110775, 110774, 110774, 110776, 110776, 110776, 110777, 110776, 110776, 110776, 110776, 110776, 110776, 110775, 110773, 110773, 110773, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110774, 110775, 110774, 110774, 110775, 110775, 110775, 110776, 110777, 110777, 110777, 110777, 110778, 110778, 110779, 110778, 110780, 110781, 110781, 110781, 110779, 110780, 110780, 110781, 110779, 110779, 110779, 110779, 110779, 110782, 110782, 110781, 110781, 110782, 110782, 110782, 110782, 110782, 110782, 110781, 110783, 110782, 110782, 110782, 110783, 110782, 110782, 110783, 110780, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110782, 110782, 110781, 110781, 110781, 110781, 110781, 110781, 110781, 110782, 110781, 110781, 110784, 110782, 110782, 110784, 110782, 110782, 110782, 110782, 110781, 110782, 110781, 110782, 110782, 110782, 110781, 110782, 110782, 110783, 110783, 110782, 110783, 110783, 110784, 110784, 110784, 110786, 110786, 110785, 110786, 110786, 110786, 110786, 110785, 110785, 110784, 110785, 110785, 110785, 110785, 110785, 110785, 110786, 110787, 110787, 110787, 110787, 110786, 110786, 110786, 110786, 110786, 110787, 110787, 110787, 110787, 110787, 110786, 110787, 110787, 110786, 110786, 110785, 110785, 110786, 110786, 110786, 110786, 110786, 110784, 110786, 110786, 110786, 110786, 110785, 110786, 110786, 110785, 110785, 110782, 110783, 110785, 110785, 110786, 110785, 110785, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110784, 110783, 110784, 110783, 110783, 110783, 110784, 110784, 110784, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110746, 110746, 110746, 110745, 110746, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110748, 110748, 110748, 110748, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110747, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110748, 110747, 110745, 110745, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110747, 110746, 110747, 110745, 110745, 110744, 110745, 110744, 110744, 110744, 110745, 110745, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110744, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110745, 110743, 110743, 110743, 110743, 110744, 110741, 110743, 110743, 110743, 110740, 110741, 110740, 110741, 110743, 110743, 110743, 110741, 110741, 110741, 110741, 110741, 110740, 110741, 110740, 110740, 110741, 110741, 110740, 110741, 110741, 110742, 110742, 110742, 110742, 110741, 110742, 110742, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110741, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110742, 110741, 110741, 110741, 110742, 110742, 110741, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110742, 110742, 110742, 110742, 110743, 110742, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110742, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110743, 110743, 110744, 110743, 110743, 110743, 110743, 110743, 110744, 110744, 110744, 110742, 110741, 110741, 110742, 110742, 110740, 110740, 110741, 110741, 110741, 110740, 110740, 110740, 110739, 110740, 110739, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110739, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110740, 110739, 110739, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110736, 110738, 110736, 110738, 110737, 110738, 110738, 110738, 110738, 110739, 110739, 110736, 110736, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110739, 110739, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110740, 110740, 110738, 110737, 110737, 110739, 110739, 110740, 110740, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110740, 110739, 110739, 110739, 110739, 110739, 110740, 110738, 110737, 110737, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110739, 110738, 110739, 110739, 110739, 110738, 110738, 110738, 110738, 110738, 110739, 110739, 110738, 110739, 110739, 110739, 110740, 110740, 110740, 110740, 110739, 110739, 110738, 110738, 110738, 110736, 110736, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110734, 110733, 110734, 110733, 110735, 110735, 110735, 110735, 110736, 110737, 110736, 110736, 110736, 110736, 110735, 110736, 110736, 110733, 110734, 110734, 110734, 110734, 110734, 110733, 110734, 110733, 110733, 110734, 110734, 110734, 110733, 110733, 110734, 110734, 110734, 110732, 110734, 110733, 110733, 110731, 110732, 110731, 110732, 110733, 110733, 110733, 110731, 110731, 110732, 110730, 110731, 110733, 110732, 110732]],
            'npts_0.5': [False, [63987, 63987, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63989, 63989, 63990, 63990, 63990, 63989, 63989, 63990, 63990, 63991, 63990, 63990, 63991, 63991, 63991, 63991, 63990, 63991, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63990, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63990, 63989, 63990, 63990, 63990, 63989, 63989, 63988, 63988, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63990, 63989, 63989, 63988, 63989, 63989, 63989, 63990, 63990, 63991, 63991, 63991, 63991, 63994, 63995, 63996, 63993, 63993, 63994, 63994, 63993, 63993, 63993, 63993, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63994, 63993, 63992, 63992, 63992, 63993, 63994, 63994, 63994, 63991, 63990, 63990, 63990, 63990, 63990, 63991, 63991, 63991, 63991, 63991, 63989, 63990, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63991, 63991, 63991, 63992, 63992, 63992, 63992, 63992, 63992, 63992, 63991, 63990, 63990, 63990, 63990, 63989, 63988, 63988, 63988, 63988, 63988, 63988, 63987, 63989, 63989, 63989, 63988, 63989, 63989, 63988, 63987, 63986, 63986, 63985, 63985, 63986, 63985, 63985, 63986, 63985, 63984, 63985, 63984, 63984, 63986, 63986, 63985, 63986, 63986, 63986, 63986, 63984, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63982, 63983, 63984, 63985, 63986, 63986, 63987, 63988, 63986, 63984, 63984, 63984, 63984, 63983, 63983, 63985, 63984, 63984, 63983, 63984, 63983, 63983, 63984, 63984, 63984, 63984, 63984, 63983, 63981, 63981, 63981, 63981, 63982, 63982, 63983, 63983, 63983, 63981, 63982, 63983, 63982, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63984, 63985, 63984, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63981, 63981, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63980, 63980, 63979, 63980, 63979, 63979, 63980, 63980, 63980, 63980, 63980, 63981, 63981, 63982, 63982, 63982, 63983, 63982, 63983, 63983, 63983, 63984, 63985, 63985, 63985, 63985, 63985, 63984, 63984, 63983, 63984, 63984, 63983, 63983, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63983, 63984, 63984, 63984, 63985, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63986, 63985, 63984, 63984, 63984, 63984, 63984, 63984, 63983, 63984, 63983, 63983, 63983, 63985, 63985, 63984, 63984, 63987, 63985, 63985, 63986, 63985, 63985, 63985, 63985, 63984, 63984, 63987, 63988, 63987, 63987, 63987, 63987, 63987, 63986, 63987, 63987, 63988, 63987, 63987, 63986, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63989, 63985, 63984, 63988, 63988, 63988, 63988, 63986, 63985, 63986, 63986, 63986, 63986, 63987, 63986, 63985, 63985, 63985, 63986, 63986, 63987, 63987, 63987, 63986, 63985, 63985, 63986, 63985, 63986, 63985, 63986, 63985, 63985, 63985, 63985, 63985, 63985, 63984, 63985, 63985, 63984, 63985, 63976, 63976, 63974, 63974, 63974, 63974, 63975, 63974, 63974, 63974, 63973, 63973, 63972, 63973, 63973, 63974, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63972, 63972, 63971, 63971, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63973, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63969, 63970, 63971, 63969, 63970, 63972, 63972, 63975, 63975, 63974, 63974, 63974, 63975, 63975, 63975, 63974, 63973, 63973, 63973, 63973, 63972, 63973, 63973, 63972, 63972, 63971, 63969, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63972, 63972, 63973, 63972, 63972, 63972, 63972, 63972, 63972, 63971, 63970, 63970, 63971, 63971, 63971, 63971, 63970, 63969, 63969, 63969, 63969, 63970, 63970, 63971, 63971, 63971, 63970, 63971, 63970, 63970, 63968, 63970, 63970, 63970, 63971, 63973, 63972, 63972, 63973, 63973, 63974, 63974, 63974, 63973, 63974, 63974, 63974, 63974, 63968, 63969, 63968, 63968, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63970, 63970, 63970, 63968, 63969, 63969, 63970, 63971, 63971, 63971, 63971, 63971, 63971, 63971, 63972, 63970, 63971, 63969, 63969, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63967, 63966, 63966, 63967, 63966, 63966, 63966, 63966, 63968, 63968, 63968, 63968, 63968, 63968, 63967, 63965, 63965, 63966, 63965, 63965, 63966, 63965, 63966, 63965, 63966, 63962, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63963, 63963, 63963, 63962, 63964, 63964, 63964, 63963, 63963, 63962, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63961, 63961, 63962, 63961, 63960, 63961, 63961, 63963, 63963, 63961, 63961, 63962, 63962, 63961, 63961, 63961, 63961, 63961, 63961, 63962, 63961, 63961, 63961, 63961, 63960, 63959, 63959, 63958, 63958, 63958, 63958, 63959, 63959, 63958, 63959, 63960, 63959, 63958, 63958, 63959, 63959, 63959, 63960, 63960, 63959, 63960, 63959, 63960, 63960, 63960, 63959, 63959, 63960, 63961, 63961, 63961, 63961, 63962, 63962, 63962, 63963, 63963, 63962, 63962, 63962, 63962, 63962, 63962, 63963, 63963, 63963, 63962, 63961, 63962, 63961, 63962, 63961, 63961, 63961, 63962, 63962, 63962, 63962, 63962, 63962, 63961, 63963, 63963, 63961, 63963, 63963, 63962, 63963, 63962, 63962, 63963, 63962, 63964, 63963, 63963, 63963, 63963, 63962, 63962, 63962, 63961, 63960, 63959, 63956, 63956, 63958, 63958, 63958, 63958, 63958, 63961, 63961, 63961, 63962, 63962, 63961, 63960, 63960, 63960, 63958, 63958, 63958, 63958, 63958, 63958, 63958, 63960, 63961, 63961, 63961, 63961, 63960, 63959, 63960, 63960, 63960, 63963, 63962, 63963, 63961, 63961, 63961, 63961, 63959, 63958, 63960, 63960, 63960, 63958, 63958, 63958, 63959, 63960, 63960, 63958, 63958, 63959, 63959, 63959, 63958, 63958, 63958, 63958, 63958, 63958, 63957, 63957, 63957, 63957, 63954, 63953, 63953, 63953, 63952, 63954, 63954, 63954, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63955, 63953, 63953, 63953, 63952, 63952, 63952, 63952, 63953, 63952, 63951, 63951, 63951, 63951, 63952, 63952, 63952, 63951, 63953, 63952, 63951, 63951, 63950, 63949, 63949, 63949, 63949, 63947, 63947, 63947, 63948, 63948, 63948, 63948, 63949, 63946, 63946, 63946, 63948, 63948, 63948]],
            'npts_real': [True, 191116800]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict['test_mosaic_cube_eph_briggsbwtaper']['self.parallel'] = self.parallel
        test_dict['test_mosaic_cube_eph_briggsbwtaper']['report'] = report
        test_dict['test_mosaic_cube_eph_briggsbwtaper']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.01, 0.1])
        self.mom8_creator(img+'.residual', range_list=[-0.01, 0.1])
        test_dict['test_mosaic_cube_eph_briggsbwtaper']['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict['test_mosaic_cube_eph_briggsbwtaper']['images'].append(img+'.image.profile.png')

        if savemetricdict:
            ### serialize ndarray in mask_stats_dcit
            mask_stats_mod_dict = copy.deepcopy(mask_stats_dict)
            mask_stats_mod_dict['mask'] = mask_stats_dict['mask'].tolist()
            #create a nested dictionary containing exp dictionaries to save
            savedict = {}
            #list of stats to save
            # im_stats, mask_stats, pb_stats, psf_stats,\
            # model_stats, resid_stats, sumwt_stats, + wt_stats (mosaic)
            savedict['im_stats_dict']=im_stats_dict
            savedict['mask_stats_dict']=mask_stats_mod_dict
            savedict['pb_stats_dict']=pb_stats_dict
            savedict['psf_stats_dict']=psf_stats_dict
            savedict['model_stats_dict']=model_stats_dict
            savedict['resid_stats_dict']=resid_stats_dict
            savedict['sumwt_stats_dict']=sumwt_stats_dict
            savedict['wt_stats_dict']=wt_stats_dict
            
            self.save_dict_to_file(test_name,savedict, test_name+'_cas13317_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_cube_eph_briggsbwtaper
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_mfs_eph(self):
        ''' Mosaic ephemeris mfs imaging - field Venus, spw 25 & 45 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel, verbose=True)

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
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report (test_mosaic_mfs_eph)
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        exp_im_stats = {'com_bmaj': [False, 0.914250373840332],
            'com_bmin': [False, 0.708504319190979],
            'com_pa': [False, -89.3270263671875],
            'npts': [True, 201600],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.057927370071411],
            'max_val_pos': [True, [320, 219, 0, 0]],
            'min_val': [False, -2.188678741455078],
            'min_val_pos': [True, [319, 159, 0, 0]],
            'im_rms': [False, 0.676557465791],
            'im_sum': [False, 5498.32523989],
            # CAS-9386 update build100 serial
            #'regn_sum': [False, 1404.0369635481302],
            'regn_sum': [False, 1381.8595210092262],
            'npts_real': [True, 201600],
            'rms_per_field': [False, [0.956495428743, 0.812393759263, 0.850127531746, 0.857228981072, 0.826569686485, 0.855403274215, 0.791376139627]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [239, 396, 0, 0]), \
                      # CAS-9386 update build100 serial
                      #(img+'.image', False, [239, 397, 0, 0]), \
                      (img+'.image', True, [239, 397, 0, 0]), \
                      (img+'.image', True, [47, 210, 0, 0]), \
                      (img+'.image', False, [46, 210, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_mfs_eph)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 201600],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 201600]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_mfs_eph)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[239.37095571deg, -16.96411290deg], [28.1142arcsec, 27.7734arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 201600],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000052452087402],
            'im_rms': [False, 0.6312907916397755],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 113589],
            #'npts_0.5': [True, 64574],
            'npts_0.2': [1e-4, 114034],
            'npts_0.5': [1e-4, 64662],
            'npts_real': [True, 201600],
            'fit': [False, [1.0977556311256869, 37.34956230416832, \
                    36.99775156676905]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442262646273],
            'fit_pix': [False, [240.86482317132828, 210.08148532276593]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_mfs_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36978024deg, -16.96392002deg], [1.1516arcsec, 0.9492arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.05114944279193878],
            'min_val_pos': [True, [250, 204, 0, 0]],
            'im_rms': [False, 0.011101956781551448],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 0.17176903218875061],
            'im_sum': [False, 0.16332336102429368],
            'npts_real': [True, 201600],
            'fit': [False, [0.8980171961947707, 1.0457922779210116, \
                    0.8221985921765811]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442262646273],
            'fit_pix': [False, [239.99681799687036, 209.99815643954449]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_mfs_eph)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 201600],
            # CAS-9386  update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.057927370071411],
            'max_val_pos': [True, [320, 219, 0, 0]],
            'min_val': [False, -2.188678741455078],
            'min_val_pos': [True, [319, 159, 0, 0]],
            'im_rms': [False, 0.6765198782434925],
            'im_sum': [False, 5498.612235217561],
            # CAS-9386 update build100 serial
            #'regn_sum': [False, 1404.0369635481302],
            'regn_sum': [False, 1381.8595210092262],
            'npts_real': [True, 201600]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_mfs_eph)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
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

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_mfs_eph)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [1e-4, 1.0],
            'freq_bin': [1e-10, 16762504556.453674],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_mfs_eph)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        exp_wt_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.3335382044315338],
            'max_val_pos': [True, [240, 210, 0, 0]],
            # CAS-9386 update build100 serial
            #'min_val': [False, 9.392295760335401e-05],
            'min_val': [False, 9.849096386460587e-05],
            'im_rms': [False, 0.12506836881],
            'im_sum': [False, 15366.9703442],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 113589],
            #'npts_0.5': [True, 64574],
            'npts_0.2': [1e-4, 114034],
            'npts_0.5': [1e-4, 64662],
            'npts_real': [True, 201600]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict['test_mosaic_mfs_eph']['self.parallel'] = self.parallel
        test_dict['test_mosaic_mfs_eph']['report'] = report
        test_dict['test_mosaic_mfs_eph']['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-2.2, 2.1])
        self.mom8_creator(img+'.residual', range_list=[-2.2, 2.1])
        test_dict['test_mosaic_mfs_eph']['images'].extend( \
                (img+'.image.moment8.png',img+'.residual.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_mfs_eph
#-------------------------------------------------#
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_mtmfs_eph(self):
        ''' Mosaic ephemeris mtmfs imaging - field Venus, spw 25 & 45 '''

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
            fastnoise=False, savemodel='none', parallel=self.parallel, verbose=True)

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
            parallel=self.parallel, verbose=True)

        report0 = th.checkall(imgexist = self.image_list(img, 'mos_mtmfs'))

        # .image report (test_mosaic_mtmfs_eph)
        im_stats_dict = self.image_stats(img+'.image.tt0', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]', \
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
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.054377794265747],
            'max_val_pos': [True, [320, 218, 0, 0]],
            'min_val': [False, -2.195145845413208],
            'min_val_pos': [True, [319, 159, 0, 0]],
            'im_rms': [False, 0.6725999217209672],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 5402.754884014255],
            #'regn_sum': [False, 1406.6249240019824],
            'im_sum': [False, 5458.724023531031],
            'regn_sum': [False, 1383.755136466585],
            'npts_real': [True, 201600],
            'rms_per_field': [False, [0.948501720585, 0.806080633033, 0.843574363826, 0.852078553708, 0.821227271327, 0.850437711204, 0.784834956688]]}

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [239, 396, 0, 0]), \
                      # CAS-9386 update build100 serial
                      #(img+'.image.tt0', False, [239, 397, 0, 0]), \
                      (img+'.image.tt0', True, [239, 397, 0, 0]), \
                      (img+'.image.tt0', True, [47, 210, 0, 0]), \
                      (img+'.image.tt0', False, [46, 210, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image.tt0', epsilon=self.epsilon)

        # .mask report (test_mosaic_mtmfs_eph)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = {'npts': [True, 201600],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'mask_pix': [True, 0],
            'mask_regns': [True, 0],
            'npts_real': [True, 201600]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_mtmfs_eph)
        pb_stats_dict = self.image_stats(img+'.pb.tt0', fit_region = \
            'ellipse[[239.37095571deg, -16.96411290deg], [28.1142arcsec, 27.7734arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 201600],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000052452087402],
            'im_rms': [False, 0.6312907916397755],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 113589],
            #'npts_0.5': [True, 64574],
            'npts_0.2': [1e-4, 114034],
            'npts_0.5': [1e-4, 64662],
            'npts_real': [True, 201600],
            'fit': [False, [1.0977556311256869, 37.34956230416832, \
                    36.99775156676905]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442262646273],
            'fit_pix': [False, [240.86482317132828, 210.08148532276593]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb.tt0', epsilon=self.epsilon)

        # .psf report (test_mosaic_mtmfs_eph)
        psf_stats_dict = self.image_stats(img+'.psf.tt0', fit_region = \
            'ellipse[[239.36978024deg, -16.96392002deg], [1.1516arcsec, 0.9492arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.05114944279193878],
            'min_val_pos': [True, [250, 204, 0, 0]],
            'im_rms': [False, 0.011101956781551448],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 0.17176903218875061],
            'im_sum': [False, 0.16332336102429368],
            'npts_real': [True, 201600],
            'fit': [False, [0.8980212570855989, 1.0458854777504984, \
                        0.8222593788495552]],
            'fit_loc_chan': [True, 0],
            'fit_loc_freq': [1e-10, 253.57442262646273],
            'fit_pix': [False, [239.96621779301014, 209.99390876796625]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf.tt0', epsilon=self.epsilon)

        # .residual report (test_mosaic_mtmfs_eph)
        resid_stats_dict = self.image_stats(img+'.residual.tt0', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 201600],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 2.057927370071411],
            'max_val_pos': [True, [320, 219, 0, 0]],
            'min_val': [False, -2.188678741455078],
            'min_val_pos': [True, [319, 159, 0, 0]],
            'im_rms': [False, 0.6765198782434925],
            'im_sum': [False, 5498.612235217561],
            # CAS-9386 update build100 serial
            #'regn_sum': [False, 1404.0369635481302],
            'regn_sum': [False, 1381.8595210092262],
            'npts_real': [True, 201600]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0', epsilon=self.epsilon)

        # .model report (test_mosaic_mtmfs_eph)
        model_stats_dict = self.image_stats(img+'.model.tt0', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
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

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_mtmfs_eph)
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = {'npts': [True, 1],
            'npts_unmasked': [1e-4, 1.0],
            'freq_bin': [1e-10, 16762504556.453674],
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

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0', epsilon=self.epsilon)

        # .weight report (test_mosaic_mtmfs_eph)
        wt_stats_dict = self.image_stats(img+'.weight.tt0', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        exp_wt_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.3335382044315338],
            'max_val_pos': [True, [240, 210, 0, 0]],
            # CAS-9386 update build100 serial
            #'min_val': [False, 9.392295760335401e-05],
            'min_val': [False, 9.849096386460587e-05],
            'im_rms': [False, 0.12506831619046738],
            'im_sum': [False, 15366.970359400737],
            # CAS-9386 update build100 serial
            #'npts_0.2': [True, 113589],
            #'npts_0.5': [True, 64574],
            'npts_0.2': [1e-4, 114034],
            'npts_0.5': [1e-4, 64662],
            'npts_real': [True, 201600]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, \
            '.weight.tt0', epsilon=self.epsilon)

        # .image.tt1 report (test_mosaic_mtmfs_eph)
        im1_stats_dict = self.image_stats(img+'.image.tt1', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]', \
            field_regions = \
            ['circle[[15:57:28.454567, -16.57.49.11051], 11.5arcsec]',
             'circle[[15:57:28.112222, -16.57.59.87434], 11.5arcsec]',
             'circle[[15:57:29.051302, -16.58.00.17973], 11.5arcsec]',
             'circle[[15:57:27.877217, -16.57.49.98258], 11.5arcsec]',
             'circle[[15:57:29.755349, -16.57.50.59334], 11.5arcsec]',
             'circle[[15:57:28.581274, -16.57.40.39634], 11.5arcsec]',
             'circle[[15:57:29.520326, -16.57.40.70171], 11.5arcsec]'])

        exp_im1_stats = {'com_bmaj': [False, 0.914250373840332],
            'com_bmin': [False, 0.708504319190979],
            'com_pa': [False, -89.3270263671875],
            'npts': [True, 201600],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 14.496827125549316],
            'max_val_pos': [True, [202, 238, 0, 0]],
            'min_val': [False, -15.082664489746094],
            'min_val_pos': [True, [225, 203, 0, 0]],
            'im_rms': [False, 3.6352690257406723],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 14670.480392508209],
            #'regn_sum': [False, -437.08825725317],
            'im_sum': [False, 13760.151578558609],
            'regn_sum': [False, -320.53889639303],
            'npts_real': [True, 201600],
            'rms_per_field': [False, [4.8419718773, 3.97862920107, 3.92391811, 4.1641374813, 3.58102697509, 3.96398521308, 3.53341315536]]}

        report10 = th.check_dict_vals(exp_im1_stats, im1_stats_dict, '.image.tt1', epsilon=self.epsilon)

        # .residual.tt1 report (test_mosaic_mtmfs_eph)
        resid1_stats_dict = self.image_stats(img+'.residual.tt1', \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]')

        exp_resid1_stats = {'npts': [True, 201600],
            # CAS-9386 update build100 serial
            #'npts_unmasked': [True, 113589.0],
            'npts_unmasked': [1e-4, 114034.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 0.018855653703212738],
            'max_val_pos': [True, [302, 155, 0, 0]],
            'min_val': [False, -0.023053178563714027],
            'min_val_pos': [True, [225, 203, 0, 0]],
            'im_rms': [False, 0.006006708967196801],
            # CAS-9386 update build100 serial
            #'im_sum': [False, 47.65366425671894],
            'im_sum': [False, 46.97670805919702],
            'regn_sum': [False, 7.861296703074515],
            'npts_real': [True, 201600]}

        report11 = th.check_dict_vals(exp_resid1_stats, resid1_stats_dict, \
            '.residual.tt1', epsilon=self.epsilon)

        # .model.tt1 report (test_mosaic_mtmfs_eph)
        model1_stats_dict = self.image_stats(img+'.model.tt1', fit_region = \
            'ellipse[[239.37089670deg, -16.96420698deg], [13.2095arcsec, 13.1423arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model1_stats = {'npts': [True, 201600],
            'npts_unmasked': [1e-4, 201600.0],
            'freq_bin': [1e-10, 16762504556.453674],
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

        report12 = th.check_dict_vals(exp_model1_stats, model1_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt.tt1 report (test_mosaic_mtmfs_eph)
        sumwt1_stats_dict = self.image_stats(img+'.sumwt.tt1')

        exp_sumwt1_stats = {'npts': [True, 1],
            'npts_unmasked': [1e-4, 1.0],
            'freq_bin': [1e-10, 16762504556.453674],
            'start': [True, 2.53574e+11],
            'end': [True, 2.53574e+11],
            'start_delta': [False, 2.53574e+11],
            'end_delta': [False, 2.53574e+11],
            'nchan': [True, 1],
            'max_val': [False, 213949.859375],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 213949.859375],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 213949.85580738305],
            'npts_real': [True, 1]}

        report13 = th.check_dict_vals(exp_sumwt1_stats, sumwt1_stats_dict, \
            '.sumwt.tt1', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11 + report12 + report13

        failed = self.filter_report(report)

        add_to_dict(self, output=test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        img = shutil._basename(img)
        test_dict['test_mosaic_mtmfs_eph']['self.parallel'] = self.parallel
        test_dict['test_mosaic_mtmfs_eph']['report'] = report
        test_dict['test_mosaic_mtmfs_eph']['images'] = []

        self.mom8_creator(img+'.image.tt0', range_list=[-2.2, 2.1])
        self.mom8_creator(img+'.residual.tt0', range_list=[-2.2, 2.1])
        test_dict['test_mosaic_mtmfs_eph']['images'].extend( \
            (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_mtmfs_eph


def suite():
     return [Test_standard, Test_mosaic]

# Main #
if __name__ == '__main__':
    unittest.main()

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

Test list - 21 total
1a.  Single field(SF) cube with perchanweightdensity=False(pcwdF), weighting=briggs - E2E6.1.00034.S
1b.  SF cube with pcwdT, weighting=briggs - E2E6.1.00034.S
1c.  SF cube with pcwdT, weighting=briggsbwtaper - E2E6.1.00034.S
2.   SF MFS - E2E6.1.00020.S
3.   SF mtmfs - E2E6.1.00020.S
4a.  SF ephemeris cube (multi-EB) with pcwdF+briggs - 2017.1.00750.T
4b.  SF ephemeris cube (multi-EB) with pcwdT+briggs - 2017.1.00750.T
4c.  SF ephemeris cube (multi-EB) with pcwdT+briggsbwtaper - 2017.1.00750.T
5.   SF ephemeris MFS - 2018.1.00879.S
6.   SF ephemeris mtmfs - 2018.1.00879.S
7.   SF Calibrator - E2E6.1.00034.S
8.   SF ephemeris Calibrator - 2018.1.00879.S  
9a.  Mosaic cube with pcwdF, briggs- E2E6.1.00034.S
9b.  Mosaic cube with pcwdT+brigs- E2E6.1.00034.S
9c.  Mosaic cube with pcwdT+briggsbwtaper- E2E6.1.00034.S
10.  Mosaic MFS - E2E6.1.00020.S
11.  Mosaic mtmfs - E2E6.1.00020.S
12a. Mosaic ephemeris cube with pcwdF- 2018.1.00879.S
12b. Mosaic ephemeris cube with pcwdT+briggs - 2018.1.00879.S
12c. Mosaic ephemeris cube with pcwdT+briggsbwtaper - 2018.1.00879.S
13.  Mosaic ephemeris MFS - 2018.1.00879.S
14.  Mosaic ephemeris mtmfs - 2018.1.00879.S

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

    def remove_prefix(self,string, prefix):
        ''' Remove a specified prefix string from string '''
        return string[string.startswith(prefix) and len(prefix):]  

##############################################
# TESTS
##############################################
test_dict = {}
class Test_standard(test_tclean_base):


    #Test 1a
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
            gridder='standard',  mosweight=False, \
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
                gridder='standard',  mosweight=False, \
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
                 mosweight=False, usepointing=False, \
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
                gridder='standard',  mosweight=False, \
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
    # Test 1b
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cube_pcwdT(self):
        ''' Standard (single field) cube imaging with pcwdT and briggs - central field of SMIDGE_NWCloud (field 3), spw 22 '''
        
        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
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
            gridder='standard', mosweight=False, \
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
                outframe='LSRK', perchanweightdensity=True, \
                usepointing=False, pblimit=0.2, nsigma=0.0, \
                gridder='standard', mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggs', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
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
                mosweight=False, usepointing=False, \
                pblimit=0.2, deconvolver='hogbom', restoration=True, \
                restoringbeam='common', pbcor=True, weighting='briggs', \
                robust=0.5, npixels=0, niter=0, threshold='0.354Jy', \
                nsigma=0.0, interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
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
                gridder='standard', mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggs', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
                minpercentchange=1.0, fastnoise=False, restart=True, \
                calcres=False, calcpsf=False, savemodel='none', \
                restoringbeam='common', parallel=False, verbose=True)

            # retrieve per-channel beam statistics 
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(image=img+'.psf')


        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report(test_standard_cube_pcwdT)
        im_stats_dict = self.image_stats(image=img+'.image', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 8.814602965342319],
            'com_bmin': [False, 6.206190071773739],
            'com_pa': [False, 66.71102881167695],
            'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0721781253814697],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, -0.6544958353042603],
            'min_val_pos': [True, [35, 39, 0, 249]],
            'im_rms': [False, 0.1417189969049196],
            'rms_per_chan': [False, [0.12524489491701085, 0.14325395222222756, 0.12584678102392538, 0.11145207324867318, 0.14800700454273893, 0.12384229098179399, 0.13250800971797416, 0.1090860883231357, 0.1531594196918675, 0.15835686519019695, 0.11148824643078865, 0.1402621423451732, 0.16825655657840805, 0.1324912139739508, 0.1576496247292607, 0.13677337518999708, 0.10676780725128465, 0.12578482004349698, 0.14054293457152536, 0.17541239977654075, 0.13238700370985992, 0.15362042708391535, 0.11417798107270065, 0.11107475613548828, 0.16902137589960362, 0.13582918235373295, 0.13295711576394656, 0.11946748181001618, 0.1411508312437078, 0.14131722026096524, 0.14404028109252656, 0.11966484223548988, 0.12628873297234783, 0.1216324847351966, 0.15573166407340042, 0.15638903025707565, 0.1584263825445491, 0.14220007944895904, 0.14366597604460787, 0.13051159896010986, 0.13145100820788863, 0.15131927349323315, 0.14277176576019135, 0.13688680749695017, 0.13974237035008388, 0.14198244695259946, 0.13009522695285364, 0.1343284552987026, 0.16895692538942175, 0.14005971920558666, 0.12083455408368526, 0.15749305555121293, 0.13711632732757945, 0.1283089177013551, 0.11983756684724207, 0.13163394650171878, 0.1521539060704143, 0.17515168166004855, 0.10835609122124179, 0.13889193335697847, 0.15123779456299932, 0.13873116826399493, 0.12069415502471983, 0.12402433877274839, 0.12549093291979233, 0.1392944953403145, 0.12718710993893223, 0.1465783709621025, 0.13422666853729287, 0.15796141542839223, 0.12711458229382847, 0.13461425687392203, 0.14954169375981388, 0.14080049697795835, 0.17275285554611003, 0.16474890876979764, 0.11963948380326994, 0.14776285264096514, 0.15106303405101967, 0.1204216585796539, 0.12339519904195956, 0.1645693969933035, 0.11217981199770478, 0.18008993584141478, 0.13496304338851975, 0.1432310718111576, 0.11447824071414166, 0.1342803808466825, 0.1585707369393876, 0.13374434141248762, 0.15015631733109366, 0.1317426470629964, 0.13739313067512826, 0.143413646163012, 0.12130786382117466, 0.1437681823348491, 0.12350172538926665, 0.16439956572028933, 0.14533593594406713, 0.12947613236332578, 0.1552915685320751, 0.12061872222208996, 0.14016806515007463, 0.1500214372138822, 0.1324421188236096, 0.1281523950745231, 0.15807100724699694, 0.14288762079648842, 0.15240431990199288, 0.1263952011902941, 0.1436807888624474, 0.12911318445893907, 0.12546759708434568, 0.1369038010904738, 0.13074792425516618, 0.13080655867597857, 0.11585277990785434, 0.15869508823496187, 0.1252950901665166, 0.15124752246701284, 0.13277597096694072, 0.14202227406378076, 0.1373221943926423, 0.14551525168140142, 0.12318705602452809, 0.15605531947598278, 0.142758926144495, 0.14552935351330806, 0.1291163429736027, 0.1222646805789003, 0.1450321552688544, 0.12178048305667369, 0.15675174838925315, 0.15361290460694835, 0.15692902228521985, 0.12245020068894151, 0.15963914072329885, 0.13490156638614706, 0.13594528591514057, 0.1489070817421063, 0.11642417722148743, 0.15093366607685868, 0.14466410641565836, 0.11465472475108811, 0.1603493243327231, 0.12455325290954068, 0.14199395848651591, 0.1623054512365012, 0.13040747008900846, 0.13582532354110208, 0.142021865826451, 0.10177138787940485, 0.14793468250259353, 0.1228018607033082, 0.16461575804695594, 0.14717368983455842, 0.148571614271376, 0.12806978764193133, 0.11603319136864312, 0.14791560010077368, 0.1185335464669558, 0.15396244771281936, 0.15385278642408748, 0.1732531991163082, 0.1555277548641792, 0.13772162751577746, 0.12123100282920028, 0.1269501209460228, 0.1284220158831818, 0.1411998721712427, 0.13587329957916597, 0.1293140094247406, 0.1404720858391, 0.16702366856389905, 0.12590433329668535, 0.13303074023832714, 0.18525730595573806, 0.14339782000654405, 0.15138724050655633, 0.14470290019259652, 0.1294433369229332, 0.1389458931281672, 0.13715052619480636, 0.1496455808347652, 0.10399460120245078, 0.12342442083747604, 0.1736035481448696, 0.12498233680740921, 0.1306613744068757, 0.15224491920658553, 0.11805367064811335, 0.12241730310443885, 0.14831800881058274, 0.16215939628538384, 0.12528440628292406, 0.13975796253832268, 0.1323004775538386, 0.15517188470426754, 0.16838325920077513, 0.1558285836935377, 0.14847492703998635, 0.16462899376418544, 0.13379110624729804, 0.15225896189555624, 0.1623114620050739, 0.12538544744850075, 0.13764786207647164, 0.12579172650774645, 0.12675095012739554, 0.1366135960365169, 0.15400018537671403, 0.15031806815026125, 0.14383076852055457, 0.14560401320442373, 0.1493957680186803, 0.12289816887471859, 0.1193376057647752, 0.12094510527330993, 0.1677510608703905, 0.15739855505862038, 0.14721015132193335, 0.11927101951517409, 0.13022462919569572, 0.11819170503194255, 0.13230923822936705, 0.15905747529151806, 0.18894547354377714, 0.12401463304219389, 0.12776232926845077, 0.136837390021789, 0.15969418558451334, 0.17298036859665913, 0.13360450829021434, 0.16809124367325454, 0.1431560244209695, 0.12784652836349178, 0.12069839694040896, 0.14797138610376223, 0.15840827327295529, 0.1786308247420827, 0.1367436221598973, 0.12495177192951143, 0.1628345900372722, 0.16139259937585423, 0.1288794288431209, 0.15428105300475353, 0.12654434943101284, 0.16490778876470474, 0.15694746882511199, 0.175988382534047, 0.16654059769301435, 0.16597237835285358, 0.20116918809524773, 0.20776580804911968, 0.18805825309242327, 0.2001763533633269, 0.16770189475163647, 0.16870919905647316, 0.1371351140834906, 0.12226661931450258, 0.13166713312105557, 0.15051684370939142, 0.13594438735308018, 0.12633219093605338, 0.15302153542699828, 0.16619299620756908, 0.15813515343915474, 0.14992029427456371, 0.11827856756338502, 0.13056252917289174, 0.1377697816185522, 0.13318526651907872, 0.13861090344615612, 0.14446664257263822, 0.14257712478602858, 0.15069908697496895, 0.1317034505734666, 0.14332645635038485, 0.13161528193302177, 0.15386088567766357, 0.1103966808850666, 0.1282658305217911, 0.1465277454816055, 0.12503054748666037, 0.13920212749475475, 0.11981213441975987, 0.1497498046163771, 0.11762412635314032, 0.11868325512198513, 0.13498260009077395, 0.1606264724384225, 0.15324369234189578, 0.1217086143256347, 0.1345566650056483, 0.14200625130750413, 0.13262518438188683, 0.16860787698371996, 0.14984080671930644, 0.11433075097020196, 0.10468965950025874, 0.11561260289751259, 0.1477010998904688, 0.146841171291768, 0.14124201226849853, 0.1255019666728538, 0.13597647427776116, 0.12265192563217224, 0.12347313111136181, 0.14190079009883896, 0.1625179379730895, 0.1568594342682481, 0.13396158803746328, 0.17330148152334945, 0.1689302060076563, 0.1539497875208728, 0.12476260614495098, 0.14189093277471782, 0.1376764567752186, 0.13121165216057507, 0.14483788718555038, 0.13488933911524514, 0.1528910640773625, 0.13535884600424447, 0.13894617633431885, 0.1293400819114347, 0.1337463738482851, 0.13856345437170964, 0.12854938736571653, 0.13487554005511215, 0.13149058000126856, 0.15666333476099914, 0.1284294397035068, 0.14122650758408795, 0.15305113906204645, 0.14151054701461038, 0.1614694553737268, 0.12847797015580203, 0.11819993769653907, 0.1576747779304388, 0.14281131620150772, 0.14927273778925843, 0.16179411692755996, 0.13396773098252637, 0.11774099764971954, 0.16002741575308263, 0.1383176260962914, 0.14918886641762139, 0.12008167913815912, 0.10739974761002842, 0.13847768103325622, 0.1361455058067519, 0.12429948629013535, 0.1480970320530088, 0.12908330720560632, 0.14668124670238084, 0.1530531512452235, 0.14713094169930319, 0.12647940007478897, 0.13300520548088177, 0.13168906912335726, 0.15982329481262783, 0.13656717627267317, 0.1591167308028954, 0.14558688960656205, 0.13168874567646155, 0.12827758095921585, 0.11630201772826844, 0.12222397077390755, 0.14506263131140554, 0.13507049225439993, 0.15759786976860962, 0.12584233220097027, 0.12835926122822208, 0.1586608813429134, 0.14648391709056546, 0.1251338002861393, 0.14110561676930122, 0.13838164065975408, 0.13814159522991348, 0.1336540925243465, 0.1571504036759359, 0.16545947425629867, 0.13035676779440047, 0.13346428372641092, 0.14002601377825621, 0.1360984019694493, 0.16340403640542722, 0.14920181811100805, 0.14053213351285324, 0.1577812551882318, 0.15684482263932834, 0.13524020447021376, 0.13909520615926244, 0.13789578243675427, 0.1810337835386548, 0.11663413531353857, 0.15411744267891908, 0.11988002055030943, 0.1327600310012873, 0.16425118892677, 0.1149242098465087, 0.13702970434184217, 0.12178018201903118, 0.1671827710574231, 0.12769824560019175, 0.1449758248779768, 0.13610894431786952, 0.16933639436481315, 0.13019925181846595, 0.12471874652541953, 0.16889987791222305, 0.15335521703062963, 0.14791937513521608, 0.14159610061902936, 0.1435335719227469, 0.13278867064532876, 0.13211983845851852, 0.161264026219656, 0.14807010612203764, 0.13007758805736488, 0.1420334313484663, 0.15511995111891716, 0.1414986878017996, 0.12905390438717543, 0.16408391858608887, 0.12180556296245272, 0.15237051557711043, 0.11678925841124789, 0.13509354877668883, 0.14569738371952598, 0.14578504775018494, 0.12239757032741803, 0.14391176727621435, 0.13045959310499075, 0.1424408640255801, 0.14529699008281616, 0.11430192893389196, 0.15487233091301458, 0.13874598153950946, 0.1493024132383649, 0.15436460961794887, 0.1418913711885544, 0.15407620539089198, 0.1386836170284381, 0.1350322028984825, 0.11685926518256017, 0.14941042889022763, 0.14151723831364962, 0.14178768525925706, 0.17919246325969468, 0.15390441593799337, 0.13870450330007966, 0.13993197694767937, 0.12127753792324049, 0.11573355240419174, 0.15152402470033696, 0.13525599071453207, 0.1491050868019604, 0.14213446563773416, 0.17295173376752984, 0.1322176759212858, 0.12189924401912204, 0.16398592369696519, 0.12092522928420052, 0.14449204459146553, 0.15141980732716634, 0.12532636289355095, 0.1466645778982128, 0.14063715850339065, 0.14215528693894702, 0.13173894200724207, 0.17359441947936782, 0.13917205856186796, 0.1197739466405848, 0.1590686304898883, 0.1118377938561468, 0.13457973924991823, 0.15557421951268213, 0.14105351810985578, 0.12107166255775609, 0.12027666370834235, 0.118907502961726, 0.1215842983448503, 0.15759232940772627, 0.1201644098483269, 0.13504715858058158, 0.12719096282920364, 0.14409922333300987, 0.12827148632374744, 0.1347399543186439, 0.14213989908570096, 0.14168386196133823, 0.15024249806811477, 0.13395288793441512, 0.14067095175856845, 0.14652859369257204, 0.11509059255131812, 0.16272606364316666, 0.14439702733405144, 0.1300006769831858, 0.15036020935641703, 0.16554836581318094, 0.13650957739690425, 0.1388647373335627, 0.1580022405022512, 0.14879162756998998, 0.1340061587735869, 0.13509079081931988]],
            'im_sum': [False, -25.784847760371548],
            'regn_sum': [False, 73.32121007144451],
            'npts_real': [True, 3251200],
            'profile': [False, 0.9989975377104944],
            'fit': [False, [0.9974838931599738, 14.02714300108566, 7.1313286173316905]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [38.04662126023306, 36.97689772381557]]}


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

        # test_standard_cube_pcwdT.exp_mask_stats
        exp_mask_stats = {'npts': [True, 3251200],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'mask_pix': [False, 925],
            'mask_regns': [True, 2],
            'npts_real': [True, 3251200]}


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_mask_stats
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

        # test_standard_cube_pcwdT.exp_psf_stats
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
            'min_val': [False, -0.2037811130285263],
            'min_val_pos': [True, [7, 13, 0, 473]],
            'im_rms': [False, 0.13685870644179274],
            'im_sum': [False, 6929.015781511987],
            'npts_real': [True, 3251200],
            'fit_0': [False, [1.0878730820498594, 7.966064940466852, 5.426217238869401]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [39.99905004600111, 39.9928745515336]],
            'fit_1': [False, [1.0879006481338078, 7.963335715996354, 5.42448001757671]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [39.99905224428489, 39.992886816553664]],
            'fit_2': [False, [1.0879618127708734, 7.960153326499446, 5.4225612867089]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [39.99905497352608, 39.992896143439694]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)


        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_resid_stats
        exp_resid_stats = {'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.8373207449913025],
            'max_val_pos': [True, [41, 43, 0, 255]],
            'min_val': [False, -0.6544958353042603],
            'min_val_pos': [True, [35, 39, 0, 249]],
            'im_rms': [False, 0.14157333436045208],
            'im_sum': [False, -110.88880667924738],
            'regn_sum': [False, 33.618323663948104],
            'npts_real': [True, 3251200]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_pcwdT.exp_model_stats
        exp_model_stats = {'npts': [True, 3251200],
            'npts_unmasked': [True, 3251200.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.5203386545181274],
            'max_val_pos': [True, [42, 43, 0, 256]],
            'min_val': [False, -0.03937779739499092],
            'min_val_pos': [True, [30, 38, 0, 256]],
            'im_rms': [False, 0.0004073724174074555],
            'im_sum': [False, 1.6612802222371101],
            'regn_sum': [False, 1.6987815089523792],
            'mask_non0': [True, 0],
            'npts_real': [True, 3251200]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        # test_standard_cube_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 90.36437225341797],
            'max_val_pos': [True, [0, 0, 0, 138]],
            'min_val': [False, 90.36434936523438],
            'min_val_pos': [True, [0, 0, 0, 419]],
            'im_rms': [False, 90.36436077605586],
            'npts_real': [True, 508]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube_pcwdT.exp_bmin_dict
            exp_bmin_dict = {'*0': 6.202598571777344,'*1': 6.2025957107543945,'*10': 6.202511787414551,'*100': 6.201929569244385,'*101': 6.201907634735107,'*102': 6.201900959014893,'*103': 6.201895236968994,'*104': 6.201879501342773,'*105': 6.201879501342773,'*106': 6.2018513679504395,'*107': 6.2018513679504395,'*108': 6.201825141906738,'*109': 6.201825141906738,'*11': 6.202492713928223,'*110': 6.201825141906738,'*111': 6.201825141906738,'*112': 6.201825141906738,'*113': 6.201825141906738,'*114': 6.201825141906738,'*115': 6.201825141906738,'*116': 6.201825141906738,'*117': 6.201816558837891,'*118': 6.20178747177124,'*119': 6.20178747177124,'*12': 6.202468395233154,'*120': 6.201752185821533,'*121': 6.201752185821533,'*122': 6.201727390289307,'*123': 6.20169734954834,'*124': 6.201666831970215,'*125': 6.201706409454346,'*126': 6.201698303222656,'*127': 6.201698303222656,'*128': 6.201698303222656,'*129': 6.201690196990967,'*13': 6.202468395233154,'*130': 6.201690196990967,'*131': 6.201690196990967,'*132': 6.201690196990967,'*133': 6.201659679412842,'*134': 6.201659679412842,'*135': 6.201589584350586,'*136': 6.201589584350586,'*137': 6.2015838623046875,'*138': 6.201565265655518,'*139': 6.201565265655518,'*14': 6.202468395233154,'*140': 6.201546669006348,'*141': 6.201546669006348,'*142': 6.201546669006348,'*143': 6.201546669006348,'*144': 6.201546669006348,'*145': 6.201546669006348,'*146': 6.201524257659912,'*147': 6.201524257659912,'*148': 6.201486587524414,'*149': 6.201467037200928,'*15': 6.202468395233154,'*150': 6.201467037200928,'*151': 6.201467037200928,'*152': 6.201451778411865,'*153': 6.201451778411865,'*154': 6.201433181762695,'*155': 6.201433181762695,'*156': 6.2014079093933105,'*157': 6.2014079093933105,'*158': 6.2014079093933105,'*159': 6.2014079093933105,'*16': 6.202452182769775,'*160': 6.2014079093933105,'*161': 6.201333045959473,'*162': 6.201333045959473,'*163': 6.201314926147461,'*164': 6.201298713684082,'*165': 6.201289653778076,'*166': 6.2012763023376465,'*167': 6.2012763023376465,'*168': 6.201270580291748,'*169': 6.201270580291748,'*17': 6.202420711517334,'*170': 6.201270580291748,'*171': 6.201270580291748,'*172': 6.201253890991211,'*173': 6.201243877410889,'*174': 6.20123815536499,'*175': 6.20119571685791,'*176': 6.201169490814209,'*177': 6.201157569885254,'*178': 6.201131820678711,'*179': 6.201131820678711,'*18': 6.202409267425537,'*180': 6.201108455657959,'*181': 6.201108455657959,'*182': 6.201108455657959,'*183': 6.2011027336120605,'*184': 6.201106548309326,'*185': 6.201106548309326,'*186': 6.201091766357422,'*187': 6.201091766357422,'*188': 6.201091766357422,'*189': 6.201091766357422,'*19': 6.202399730682373,'*190': 6.201091766357422,'*191': 6.201091766357422,'*192': 6.201086521148682,'*193': 6.201086521148682,'*194': 6.201086521148682,'*195': 6.201086521148682,'*196': 6.201086521148682,'*197': 6.201086521148682,'*198': 6.201085567474365,'*199': 6.2010602951049805,'*2': 6.202572345733643,'*20': 6.202399730682373,'*200': 6.2010602951049805,'*201': 6.2010602951049805,'*202': 6.201035499572754,'*203': 6.2010321617126465,'*204': 6.201019287109375,'*205': 6.201019287109375,'*206': 6.201019287109375,'*207': 6.201019287109375,'*208': 6.201019287109375,'*209': 6.201019287109375,'*21': 6.202399730682373,'*210': 6.201013565063477,'*211': 6.201013565063477,'*212': 6.201013565063477,'*213': 6.201013565063477,'*214': 6.201013565063477,'*215': 6.201013565063477,'*216': 6.200997352600098,'*217': 6.200997352600098,'*218': 6.200987339019775,'*219': 6.200987339019775,'*22': 6.202399730682373,'*220': 6.200980186462402,'*221': 6.200980186462402,'*222': 6.200980186462402,'*223': 6.200980186462402,'*224': 6.200960159301758,'*225': 6.200952053070068,'*226': 6.200952053070068,'*227': 6.200942516326904,'*228': 6.200918197631836,'*229': 6.200918197631836,'*23': 6.202399730682373,'*230': 6.200918197631836,'*231': 6.200920104980469,'*232': 6.200920104980469,'*233': 6.200920581817627,'*234': 6.200920581817627,'*235': 6.200849533081055,'*236': 6.200849533081055,'*237': 6.2008256912231445,'*238': 6.2008256912231445,'*239': 6.200802326202393,'*24': 6.202390670776367,'*240': 6.200789928436279,'*241': 6.200789928436279,'*242': 6.200789928436279,'*243': 6.200757026672363,'*244': 6.200710296630859,'*245': 6.200710296630859,'*246': 6.200710296630859,'*247': 6.200710296630859,'*248': 6.200710296630859,'*249': 6.200710296630859,'*25': 6.2023844718933105,'*250': 6.200685024261475,'*251': 6.200685024261475,'*252': 6.200636386871338,'*253': 6.200643062591553,'*254': 6.200643062591553,'*255': 6.20060920715332,'*256': 6.200591564178467,'*257': 6.200578689575195,'*258': 6.200578689575195,'*259': 6.200578689575195,'*26': 6.2023844718933105,'*260': 6.200578689575195,'*261': 6.2005720138549805,'*262': 6.2005720138549805,'*263': 6.2005720138549805,'*264': 6.2005720138549805,'*265': 6.2005720138549805,'*266': 6.200540542602539,'*267': 6.200516223907471,'*268': 6.200520038604736,'*269': 6.200520038604736,'*27': 6.2023844718933105,'*270': 6.200491905212402,'*271': 6.200461387634277,'*272': 6.200461387634277,'*273': 6.200433254241943,'*274': 6.200433254241943,'*275': 6.200433254241943,'*276': 6.20043420791626,'*277': 6.200414180755615,'*278': 6.200414180755615,'*279': 6.200412750244141,'*28': 6.202376365661621,'*280': 6.200412750244141,'*281': 6.200412750244141,'*282': 6.200412750244141,'*283': 6.200348854064941,'*284': 6.200348854064941,'*285': 6.200348854064941,'*286': 6.200348854064941,'*287': 6.200348854064941,'*288': 6.200347423553467,'*289': 6.200340747833252,'*29': 6.202376365661621,'*290': 6.200340747833252,'*291': 6.200327396392822,'*292': 6.200321197509766,'*293': 6.200321197509766,'*294': 6.200321197509766,'*295': 6.200321197509766,'*296': 6.200308322906494,'*297': 6.200302600860596,'*298': 6.200267314910889,'*299': 6.2002410888671875,'*3': 6.202572345733643,'*30': 6.202376365661621,'*300': 6.200209617614746,'*301': 6.200209617614746,'*302': 6.20019006729126,'*303': 6.200123310089111,'*304': 6.200123310089111,'*305': 6.200102806091309,'*306': 6.200102806091309,'*307': 6.200102806091309,'*308': 6.200102806091309,'*309': 6.200102806091309,'*31': 6.202376365661621,'*310': 6.200073719024658,'*311': 6.200063705444336,'*312': 6.200062274932861,'*313': 6.200062274932861,'*314': 6.200060844421387,'*315': 6.200060844421387,'*316': 6.200060844421387,'*317': 6.200060844421387,'*318': 6.200047969818115,'*319': 6.200047969818115,'*32': 6.202364444732666,'*320': 6.2000226974487305,'*321': 6.2000226974487305,'*322': 6.2000226974487305,'*323': 6.200031757354736,'*324': 6.200031757354736,'*325': 6.199997425079346,'*326': 6.199954986572266,'*327': 6.199954986572266,'*328': 6.199954986572266,'*329': 6.199954986572266,'*33': 6.202364444732666,'*330': 6.199954986572266,'*331': 6.199954986572266,'*332': 6.199954986572266,'*333': 6.199921607971191,'*334': 6.199921607971191,'*335': 6.199893474578857,'*336': 6.199893474578857,'*337': 6.199893474578857,'*338': 6.199840068817139,'*339': 6.199840068817139,'*34': 6.202364444732666,'*340': 6.199802875518799,'*341': 6.199792385101318,'*342': 6.199792385101318,'*343': 6.199737548828125,'*344': 6.19673490524292,'*345': 6.196728229522705,'*346': 6.196728229522705,'*347': 6.196728229522705,'*348': 6.196728229522705,'*349': 6.196728229522705,'*35': 6.202354431152344,'*350': 6.196728229522705,'*351': 6.196700572967529,'*352': 6.196684837341309,'*353': 6.196684837341309,'*354': 6.1966753005981445,'*355': 6.1966753005981445,'*356': 6.196662902832031,'*357': 6.196662902832031,'*358': 6.196662902832031,'*359': 6.196662902832031,'*36': 6.202354431152344,'*360': 6.196662902832031,'*361': 6.196669578552246,'*362': 6.196669578552246,'*363': 6.196669578552246,'*364': 6.196645259857178,'*365': 6.196625709533691,'*366': 6.1966047286987305,'*367': 6.1966047286987305,'*368': 6.1966047286987305,'*369': 6.196602821350098,'*37': 6.202354431152344,'*370': 6.196596145629883,'*371': 6.196596145629883,'*372': 6.196596145629883,'*373': 6.196596145629883,'*374': 6.196596145629883,'*375': 6.196590423583984,'*376': 6.196590423583984,'*377': 6.196590423583984,'*378': 6.196590423583984,'*379': 6.196590423583984,'*38': 6.202346324920654,'*380': 6.196590423583984,'*381': 6.196590423583984,'*382': 6.196590423583984,'*383': 6.196568965911865,'*384': 6.196568965911865,'*385': 6.196567058563232,'*386': 6.196557521820068,'*387': 6.19656229019165,'*388': 6.19656229019165,'*389': 6.196555137634277,'*39': 6.202346324920654,'*390': 6.196555137634277,'*391': 6.196529865264893,'*392': 6.196529865264893,'*393': 6.196529865264893,'*394': 6.196523189544678,'*395': 6.196498394012451,'*396': 6.196498394012451,'*397': 6.196498394012451,'*398': 6.196498394012451,'*399': 6.196498394012451,'*4': 6.202572345733643,'*40': 6.202346324920654,'*400': 6.1964826583862305,'*401': 6.1964640617370605,'*402': 6.1964640617370605,'*403': 6.196438789367676,'*404': 6.196423053741455,'*405': 6.196415901184082,'*406': 6.196415901184082,'*407': 6.196358680725098,'*408': 6.196358680725098,'*409': 6.196358680725098,'*41': 6.202346324920654,'*410': 6.196329116821289,'*411': 6.196329116821289,'*412': 6.196329116821289,'*413': 6.196329116821289,'*414': 6.19630241394043,'*415': 6.196282863616943,'*416': 6.196282863616943,'*417': 6.196282863616943,'*418': 6.196267604827881,'*419': 6.196212291717529,'*42': 6.202351093292236,'*420': 6.196177005767822,'*421': 6.196177005767822,'*422': 6.196177005767822,'*423': 6.196177005767822,'*424': 6.196147918701172,'*425': 6.196147918701172,'*426': 6.196139812469482,'*427': 6.196139812469482,'*428': 6.196139812469482,'*429': 6.196139812469482,'*43': 6.202351093292236,'*430': 6.196139812469482,'*431': 6.196134090423584,'*432': 6.196134090423584,'*433': 6.196134090423584,'*434': 6.196134090423584,'*435': 6.196091175079346,'*436': 6.196091175079346,'*437': 6.196091175079346,'*438': 6.196091175079346,'*439': 6.196091175079346,'*44': 6.202335357666016,'*440': 6.196091175079346,'*441': 6.196084976196289,'*442': 6.196084976196289,'*443': 6.196084976196289,'*444': 6.196084976196289,'*445': 6.196084976196289,'*446': 6.196084976196289,'*447': 6.196084976196289,'*448': 6.196084976196289,'*449': 6.196089267730713,'*45': 6.202314376831055,'*450': 6.1960859298706055,'*451': 6.1960859298706055,'*452': 6.196073532104492,'*453': 6.196053981781006,'*454': 6.196053981781006,'*455': 6.196053981781006,'*456': 6.196053504943848,'*457': 6.196053504943848,'*458': 6.196025848388672,'*459': 6.19602632522583,'*46': 6.2023091316223145,'*460': 6.19602632522583,'*461': 6.196011066436768,'*462': 6.196008682250977,'*463': 6.196008682250977,'*464': 6.196008682250977,'*465': 6.196001052856445,'*466': 6.196001052856445,'*467': 6.196001052856445,'*468': 6.196001052856445,'*469': 6.196001052856445,'*47': 6.2023091316223145,'*470': 6.196001052856445,'*471': 6.195979595184326,'*472': 6.195979595184326,'*473': 6.1959733963012695,'*474': 6.195947170257568,'*475': 6.195913314819336,'*476': 6.195878028869629,'*477': 6.195878028869629,'*478': 6.1958489418029785,'*479': 6.195834159851074,'*48': 6.2023091316223145,'*480': 6.195834159851074,'*481': 6.195834159851074,'*482': 6.195827007293701,'*483': 6.195827007293701,'*484': 6.195827007293701,'*485': 6.195802211761475,'*486': 6.195802211761475,'*487': 6.195802211761475,'*488': 6.195802211761475,'*489': 6.195802211761475,'*49': 6.2023091316223145,'*490': 6.195792198181152,'*491': 6.195792198181152,'*492': 6.195792198181152,'*493': 6.195792198181152,'*494': 6.195792198181152,'*495': 6.195792198181152,'*496': 6.195792198181152,'*497': 6.195785999298096,'*498': 6.195785999298096,'*499': 6.195785999298096,'*5': 6.202576637268066,'*50': 6.2022833824157715,'*500': 6.195734024047852,'*501': 6.195734024047852,'*502': 6.195734024047852,'*503': 6.195683479309082,'*504': 6.195683479309082,'*505': 6.195666790008545,'*506': 6.195666790008545,'*507': 6.195666790008545,'*51': 6.2022833824157715,'*52': 6.2022833824157715,'*53': 6.2022833824157715,'*54': 6.2022833824157715,'*55': 6.2022833824157715,'*56': 6.202255725860596,'*57': 6.202250003814697,'*58': 6.2022385597229,'*59': 6.2022385597229,'*6': 6.202563762664795,'*60': 6.20222806930542,'*61': 6.20222282409668,'*62': 6.202212333679199,'*63': 6.202197074890137,'*64': 6.202197074890137,'*65': 6.202197074890137,'*66': 6.202197074890137,'*67': 6.202197074890137,'*68': 6.202197074890137,'*69': 6.202188968658447,'*7': 6.202528953552246,'*70': 6.202180862426758,'*71': 6.202180862426758,'*72': 6.202174663543701,'*73': 6.202142715454102,'*74': 6.202143669128418,'*75': 6.202143669128418,'*76': 6.202136516571045,'*77': 6.202112197875977,'*78': 6.202112197875977,'*79': 6.202112197875977,'*8': 6.202528953552246,'*80': 6.202106475830078,'*81': 6.202106475830078,'*82': 6.20206880569458,'*83': 6.202059268951416,'*84': 6.202059268951416,'*85': 6.202059268951416,'*86': 6.202056884765625,'*87': 6.202056884765625,'*88': 6.202056884765625,'*89': 6.202047824859619,'*9': 6.202511787414551,'*90': 6.202047824859619,'*91': 6.202047824859619,'*92': 6.202025413513184,'*93': 6.202000617980957,'*94': 6.202000617980957,'*95': 6.2019853591918945,'*96': 6.201965808868408,'*97': 6.201940536499023,'*98': 6.201940536499023,'*99': 6.201940536499023}
            # test_standard_cube_pcwdT.exp_bmaj_dict
            exp_bmaj_dict = {'*0': 8.796586990356445,'*1': 8.796539306640625,'*10': 8.796457290649414,'*100': 8.795126914978027,'*101': 8.795134544372559,'*102': 8.79510498046875,'*103': 8.795060157775879,'*104': 8.795051574707031,'*105': 8.795051574707031,'*106': 8.795022010803223,'*107': 8.795022010803223,'*108': 8.795001983642578,'*109': 8.795001983642578,'*11': 8.796442031860352,'*110': 8.795001983642578,'*111': 8.795001983642578,'*112': 8.795001983642578,'*113': 8.795001983642578,'*114': 8.795001983642578,'*115': 8.795001983642578,'*116': 8.795001983642578,'*117': 8.794939994812012,'*118': 8.794926643371582,'*119': 8.794926643371582,'*12': 8.796428680419922,'*120': 8.794916152954102,'*121': 8.794916152954102,'*122': 8.79489803314209,'*123': 8.794774055480957,'*124': 8.79471492767334,'*125': 8.794610023498535,'*126': 8.79455280303955,'*127': 8.79455280303955,'*128': 8.79455280303955,'*129': 8.794547080993652,'*13': 8.796428680419922,'*130': 8.794547080993652,'*131': 8.794547080993652,'*132': 8.794547080993652,'*133': 8.794533729553223,'*134': 8.794533729553223,'*135': 8.794589042663574,'*136': 8.794589042663574,'*137': 8.794570922851562,'*138': 8.794499397277832,'*139': 8.794499397277832,'*14': 8.796428680419922,'*140': 8.794520378112793,'*141': 8.794520378112793,'*142': 8.794520378112793,'*143': 8.794520378112793,'*144': 8.794520378112793,'*145': 8.794520378112793,'*146': 8.794479370117188,'*147': 8.794479370117188,'*148': 8.794443130493164,'*149': 8.794440269470215,'*15': 8.796428680419922,'*150': 8.794440269470215,'*151': 8.794440269470215,'*152': 8.794458389282227,'*153': 8.794458389282227,'*154': 8.79447078704834,'*155': 8.79447078704834,'*156': 8.79446792602539,'*157': 8.79446792602539,'*158': 8.79446792602539,'*159': 8.79446792602539,'*16': 8.796411514282227,'*160': 8.79446792602539,'*161': 8.794461250305176,'*162': 8.794461250305176,'*163': 8.794451713562012,'*164': 8.794419288635254,'*165': 8.794356346130371,'*166': 8.794341087341309,'*167': 8.794341087341309,'*168': 8.794323921203613,'*169': 8.794323921203613,'*17': 8.796299934387207,'*170': 8.794323921203613,'*171': 8.794323921203613,'*172': 8.794342994689941,'*173': 8.794272422790527,'*174': 8.794251441955566,'*175': 8.794276237487793,'*176': 8.794256210327148,'*177': 8.794244766235352,'*178': 8.794232368469238,'*179': 8.794232368469238,'*18': 8.796308517456055,'*180': 8.794214248657227,'*181': 8.794214248657227,'*182': 8.794214248657227,'*183': 8.794193267822266,'*184': 8.79415512084961,'*185': 8.79415512084961,'*186': 8.794151306152344,'*187': 8.794151306152344,'*188': 8.794151306152344,'*189': 8.794151306152344,'*19': 8.796236038208008,'*190': 8.794151306152344,'*191': 8.794151306152344,'*192': 8.794109344482422,'*193': 8.794109344482422,'*194': 8.794109344482422,'*195': 8.794109344482422,'*196': 8.794109344482422,'*197': 8.794109344482422,'*198': 8.794079780578613,'*199': 8.794061660766602,'*2': 8.7965669631958,'*20': 8.796236038208008,'*200': 8.794061660766602,'*201': 8.794061660766602,'*202': 8.794042587280273,'*203': 8.793986320495605,'*204': 8.793993949890137,'*205': 8.793993949890137,'*206': 8.793993949890137,'*207': 8.793993949890137,'*208': 8.793993949890137,'*209': 8.793993949890137,'*21': 8.796236038208008,'*210': 8.793926239013672,'*211': 8.793926239013672,'*212': 8.793926239013672,'*213': 8.793926239013672,'*214': 8.793926239013672,'*215': 8.793926239013672,'*216': 8.793899536132812,'*217': 8.793899536132812,'*218': 8.793828964233398,'*219': 8.793828964233398,'*22': 8.796236038208008,'*220': 8.793795585632324,'*221': 8.793795585632324,'*222': 8.793795585632324,'*223': 8.793795585632324,'*224': 8.793817520141602,'*225': 8.793760299682617,'*226': 8.793760299682617,'*227': 8.793719291687012,'*228': 8.793706893920898,'*229': 8.793706893920898,'*23': 8.796236038208008,'*230': 8.793706893920898,'*231': 8.793670654296875,'*232': 8.793670654296875,'*233': 8.793610572814941,'*234': 8.793610572814941,'*235': 8.79366397857666,'*236': 8.79366397857666,'*237': 8.793645858764648,'*238': 8.793645858764648,'*239': 8.793609619140625,'*24': 8.796170234680176,'*240': 8.793600082397461,'*241': 8.793600082397461,'*242': 8.793600082397461,'*243': 8.793617248535156,'*244': 8.793591499328613,'*245': 8.793591499328613,'*246': 8.793591499328613,'*247': 8.793591499328613,'*248': 8.793591499328613,'*249': 8.793591499328613,'*25': 8.79615592956543,'*250': 8.793573379516602,'*251': 8.793573379516602,'*252': 8.793537139892578,'*253': 8.793503761291504,'*254': 8.793503761291504,'*255': 8.79348373413086,'*256': 8.793496131896973,'*257': 8.793444633483887,'*258': 8.793444633483887,'*259': 8.793444633483887,'*26': 8.79615592956543,'*260': 8.793444633483887,'*261': 8.793420791625977,'*262': 8.793420791625977,'*263': 8.793420791625977,'*264': 8.793420791625977,'*265': 8.793420791625977,'*266': 8.793438911437988,'*267': 8.793355941772461,'*268': 8.79331111907959,'*269': 8.79331111907959,'*27': 8.79615592956543,'*270': 8.793331146240234,'*271': 8.793317794799805,'*272': 8.793317794799805,'*273': 8.793335914611816,'*274': 8.793335914611816,'*275': 8.793335914611816,'*276': 8.793292999267578,'*277': 8.793257713317871,'*278': 8.793257713317871,'*279': 8.793228149414062,'*28': 8.796098709106445,'*280': 8.793228149414062,'*281': 8.793228149414062,'*282': 8.793228149414062,'*283': 8.793237686157227,'*284': 8.793237686157227,'*285': 8.793237686157227,'*286': 8.793237686157227,'*287': 8.793237686157227,'*288': 8.7931489944458,'*289': 8.79310417175293,'*29': 8.796098709106445,'*290': 8.79310417175293,'*291': 8.793097496032715,'*292': 8.792970657348633,'*293': 8.792970657348633,'*294': 8.792970657348633,'*295': 8.792970657348633,'*296': 8.79297924041748,'*297': 8.79287338256836,'*298': 8.792862892150879,'*299': 8.79287052154541,'*3': 8.7965669631958,'*30': 8.796098709106445,'*300': 8.792859077453613,'*301': 8.792859077453613,'*302': 8.792716026306152,'*303': 8.792729377746582,'*304': 8.792729377746582,'*305': 8.792640686035156,'*306': 8.792640686035156,'*307': 8.792640686035156,'*308': 8.792640686035156,'*309': 8.792640686035156,'*31': 8.796098709106445,'*310': 8.792640686035156,'*311': 8.792617797851562,'*312': 8.792586326599121,'*313': 8.792586326599121,'*314': 8.792549133300781,'*315': 8.792549133300781,'*316': 8.792549133300781,'*317': 8.792549133300781,'*318': 8.792539596557617,'*319': 8.792539596557617,'*32': 8.796095848083496,'*320': 8.792521476745605,'*321': 8.792521476745605,'*322': 8.792521476745605,'*323': 8.792470932006836,'*324': 8.792470932006836,'*325': 8.792487144470215,'*326': 8.792447090148926,'*327': 8.792447090148926,'*328': 8.792447090148926,'*329': 8.792447090148926,'*33': 8.796095848083496,'*330': 8.792447090148926,'*331': 8.792447090148926,'*332': 8.792447090148926,'*333': 8.792437553405762,'*334': 8.792437553405762,'*335': 8.792438507080078,'*336': 8.792438507080078,'*337': 8.792438507080078,'*338': 8.792431831359863,'*339': 8.792431831359863,'*34': 8.796095848083496,'*340': 8.792451858520508,'*341': 8.792439460754395,'*342': 8.792439460754395,'*343': 8.792364120483398,'*344': 8.80567455291748,'*345': 8.805624961853027,'*346': 8.805624961853027,'*347': 8.805624961853027,'*348': 8.805624961853027,'*349': 8.805624961853027,'*35': 8.796107292175293,'*350': 8.805624961853027,'*351': 8.805624008178711,'*352': 8.805524826049805,'*353': 8.805524826049805,'*354': 8.805410385131836,'*355': 8.805410385131836,'*356': 8.805398941040039,'*357': 8.805398941040039,'*358': 8.805398941040039,'*359': 8.805398941040039,'*36': 8.796107292175293,'*360': 8.805398941040039,'*361': 8.80532169342041,'*362': 8.80532169342041,'*363': 8.80532169342041,'*364': 8.80530834197998,'*365': 8.805331230163574,'*366': 8.805347442626953,'*367': 8.805347442626953,'*368': 8.805347442626953,'*369': 8.805312156677246,'*37': 8.796107292175293,'*370': 8.805290222167969,'*371': 8.805290222167969,'*372': 8.805290222167969,'*373': 8.805290222167969,'*374': 8.805290222167969,'*375': 8.80522346496582,'*376': 8.80522346496582,'*377': 8.80522346496582,'*378': 8.80522346496582,'*379': 8.80522346496582,'*38': 8.796062469482422,'*380': 8.80522346496582,'*381': 8.80522346496582,'*382': 8.80522346496582,'*383': 8.805221557617188,'*384': 8.805221557617188,'*385': 8.80518913269043,'*386': 8.805163383483887,'*387': 8.805095672607422,'*388': 8.805095672607422,'*389': 8.80504322052002,'*39': 8.796062469482422,'*390': 8.80504322052002,'*391': 8.805025100708008,'*392': 8.805025100708008,'*393': 8.805025100708008,'*394': 8.804977416992188,'*395': 8.804946899414062,'*396': 8.804946899414062,'*397': 8.804946899414062,'*398': 8.804946899414062,'*399': 8.804946899414062,'*4': 8.7965669631958,'*40': 8.796062469482422,'*400': 8.804937362670898,'*401': 8.804951667785645,'*402': 8.804951667785645,'*403': 8.804920196533203,'*404': 8.804859161376953,'*405': 8.804807662963867,'*406': 8.804807662963867,'*407': 8.80483627319336,'*408': 8.80483627319336,'*409': 8.80483627319336,'*41': 8.796062469482422,'*410': 8.804824829101562,'*411': 8.804824829101562,'*412': 8.804824829101562,'*413': 8.804824829101562,'*414': 8.804805755615234,'*415': 8.804793357849121,'*416': 8.804793357849121,'*417': 8.804793357849121,'*418': 8.804725646972656,'*419': 8.804572105407715,'*42': 8.796002388000488,'*420': 8.80457878112793,'*421': 8.80457878112793,'*422': 8.80457878112793,'*423': 8.80457878112793,'*424': 8.804522514343262,'*425': 8.804522514343262,'*426': 8.804516792297363,'*427': 8.804516792297363,'*428': 8.804516792297363,'*429': 8.804516792297363,'*43': 8.796002388000488,'*430': 8.804516792297363,'*431': 8.80447769165039,'*432': 8.80447769165039,'*433': 8.80447769165039,'*434': 8.80447769165039,'*435': 8.804488182067871,'*436': 8.804488182067871,'*437': 8.804488182067871,'*438': 8.804488182067871,'*439': 8.804488182067871,'*44': 8.795976638793945,'*440': 8.804488182067871,'*441': 8.804443359375,'*442': 8.804443359375,'*443': 8.804443359375,'*444': 8.804443359375,'*445': 8.804443359375,'*446': 8.804443359375,'*447': 8.804443359375,'*448': 8.804443359375,'*449': 8.804388999938965,'*45': 8.795985221862793,'*450': 8.804343223571777,'*451': 8.804343223571777,'*452': 8.804296493530273,'*453': 8.804293632507324,'*454': 8.804293632507324,'*455': 8.804293632507324,'*456': 8.804268836975098,'*457': 8.804268836975098,'*458': 8.804217338562012,'*459': 8.804165840148926,'*46': 8.795942306518555,'*460': 8.804165840148926,'*461': 8.804170608520508,'*462': 8.804133415222168,'*463': 8.804133415222168,'*464': 8.804133415222168,'*465': 8.804113388061523,'*466': 8.804113388061523,'*467': 8.804113388061523,'*468': 8.804113388061523,'*469': 8.804113388061523,'*47': 8.795942306518555,'*470': 8.804113388061523,'*471': 8.804097175598145,'*472': 8.804097175598145,'*473': 8.804068565368652,'*474': 8.804055213928223,'*475': 8.804040908813477,'*476': 8.804031372070312,'*477': 8.804031372070312,'*478': 8.8040189743042,'*479': 8.804037094116211,'*48': 8.795942306518555,'*480': 8.804037094116211,'*481': 8.804037094116211,'*482': 8.803960800170898,'*483': 8.803960800170898,'*484': 8.803960800170898,'*485': 8.803958892822266,'*486': 8.803958892822266,'*487': 8.803958892822266,'*488': 8.803958892822266,'*489': 8.803958892822266,'*49': 8.795942306518555,'*490': 8.803946495056152,'*491': 8.803946495056152,'*492': 8.803946495056152,'*493': 8.803946495056152,'*494': 8.803946495056152,'*495': 8.803946495056152,'*496': 8.803946495056152,'*497': 8.80384349822998,'*498': 8.80384349822998,'*499': 8.80384349822998,'*5': 8.796513557434082,'*50': 8.795923233032227,'*500': 8.803811073303223,'*501': 8.803811073303223,'*502': 8.803811073303223,'*503': 8.803791046142578,'*504': 8.803791046142578,'*505': 8.803772926330566,'*506': 8.803772926330566,'*507': 8.803772926330566,'*51': 8.795923233032227,'*52': 8.795879364013672,'*53': 8.795879364013672,'*54': 8.795879364013672,'*55': 8.795879364013672,'*56': 8.795859336853027,'*57': 8.795833587646484,'*58': 8.795787811279297,'*59': 8.795787811279297,'*6': 8.796442985534668,'*60': 8.795748710632324,'*61': 8.795729637145996,'*62': 8.795721054077148,'*63': 8.795741081237793,'*64': 8.795741081237793,'*65': 8.795741081237793,'*66': 8.795741081237793,'*67': 8.795741081237793,'*68': 8.795741081237793,'*69': 8.79567813873291,'*7': 8.796438217163086,'*70': 8.795636177062988,'*71': 8.795636177062988,'*72': 8.795610427856445,'*73': 8.795598030090332,'*74': 8.795546531677246,'*75': 8.795546531677246,'*76': 8.79549503326416,'*77': 8.795449256896973,'*78': 8.795449256896973,'*79': 8.795449256896973,'*8': 8.796438217163086,'*80': 8.795453071594238,'*81': 8.795453071594238,'*82': 8.79543399810791,'*83': 8.795365333557129,'*84': 8.795365333557129,'*85': 8.795365333557129,'*86': 8.795273780822754,'*87': 8.795273780822754,'*88': 8.795273780822754,'*89': 8.79525375366211,'*9': 8.796457290649414,'*90': 8.79525375366211,'*91': 8.79525375366211,'*92': 8.79525089263916,'*93': 8.79524040222168,'*94': 8.79524040222168,'*95': 8.79521656036377,'*96': 8.795132637023926,'*97': 8.795145034790039,'*98': 8.795145034790039,'*99': 8.795145034790039}
            # test_standard_cube_pcwdT.exp_pa_dict
            exp_pa_dict = {'*0': 66.75178527832031,'*1': 66.75125122070312,'*10': 66.75174713134766,'*100': 66.74695587158203,'*101': 66.7470474243164,'*102': 66.74665832519531,'*103': 66.746337890625,'*104': 66.74657440185547,'*105': 66.74657440185547,'*106': 66.74636840820312,'*107': 66.74636840820312,'*108': 66.74665832519531,'*109': 66.74665832519531,'*11': 66.75209045410156,'*110': 66.74665832519531,'*111': 66.74665832519531,'*112': 66.74665832519531,'*113': 66.74665832519531,'*114': 66.74665832519531,'*115': 66.74665832519531,'*116': 66.74665832519531,'*117': 66.74635314941406,'*118': 66.7467269897461,'*119': 66.7467269897461,'*12': 66.75248718261719,'*120': 66.74711608886719,'*121': 66.74711608886719,'*122': 66.74738311767578,'*123': 66.7466049194336,'*124': 66.74657440185547,'*125': 66.74715423583984,'*126': 66.74685668945312,'*127': 66.74685668945312,'*128': 66.74685668945312,'*129': 66.74626159667969,'*13': 66.75248718261719,'*130': 66.74626159667969,'*131': 66.74626159667969,'*132': 66.74626159667969,'*133': 66.74662780761719,'*134': 66.74662780761719,'*135': 66.74685668945312,'*136': 66.74685668945312,'*137': 66.74657440185547,'*138': 66.7468032836914,'*139': 66.7468032836914,'*14': 66.75248718261719,'*140': 66.74665832519531,'*141': 66.74665832519531,'*142': 66.74665832519531,'*143': 66.74665832519531,'*144': 66.74665832519531,'*145': 66.74665832519531,'*146': 66.74634552001953,'*147': 66.74634552001953,'*148': 66.7463607788086,'*149': 66.74656677246094,'*15': 66.75248718261719,'*150': 66.74656677246094,'*151': 66.74656677246094,'*152': 66.74593353271484,'*153': 66.74593353271484,'*154': 66.74592590332031,'*155': 66.74592590332031,'*156': 66.74618530273438,'*157': 66.74618530273438,'*158': 66.74618530273438,'*159': 66.74618530273438,'*16': 66.75284576416016,'*160': 66.74618530273438,'*161': 66.74657440185547,'*162': 66.74657440185547,'*163': 66.74687194824219,'*164': 66.7476806640625,'*165': 66.74737548828125,'*166': 66.7477035522461,'*167': 66.7477035522461,'*168': 66.74742126464844,'*169': 66.74742126464844,'*17': 66.75294494628906,'*170': 66.74742126464844,'*171': 66.74742126464844,'*172': 66.7467269897461,'*173': 66.74636840820312,'*174': 66.74606323242188,'*175': 66.74544525146484,'*176': 66.74573516845703,'*177': 66.74598693847656,'*178': 66.74629211425781,'*179': 66.74629211425781,'*18': 66.75291442871094,'*180': 66.74655151367188,'*181': 66.74655151367188,'*182': 66.74655151367188,'*183': 66.74626159667969,'*184': 66.74629974365234,'*185': 66.74629974365234,'*186': 66.7464599609375,'*187': 66.7464599609375,'*188': 66.7464599609375,'*189': 66.7464599609375,'*19': 66.7525634765625,'*190': 66.7464599609375,'*191': 66.7464599609375,'*192': 66.74620056152344,'*193': 66.74620056152344,'*194': 66.74620056152344,'*195': 66.74620056152344,'*196': 66.74620056152344,'*197': 66.74620056152344,'*198': 66.74610137939453,'*199': 66.74637603759766,'*2': 66.75110626220703,'*20': 66.7525634765625,'*200': 66.74637603759766,'*201': 66.74637603759766,'*202': 66.74664306640625,'*203': 66.74602508544922,'*204': 66.74602508544922,'*205': 66.74602508544922,'*206': 66.74602508544922,'*207': 66.74602508544922,'*208': 66.74602508544922,'*209': 66.74602508544922,'*21': 66.7525634765625,'*210': 66.74577331542969,'*211': 66.74577331542969,'*212': 66.74577331542969,'*213': 66.74577331542969,'*214': 66.74577331542969,'*215': 66.74577331542969,'*216': 66.746337890625,'*217': 66.746337890625,'*218': 66.74706268310547,'*219': 66.74706268310547,'*22': 66.7525634765625,'*220': 66.74665069580078,'*221': 66.74665069580078,'*222': 66.74665069580078,'*223': 66.74665069580078,'*224': 66.7465591430664,'*225': 66.74626922607422,'*226': 66.74626922607422,'*227': 66.74566650390625,'*228': 66.74595642089844,'*229': 66.74595642089844,'*23': 66.7525634765625,'*230': 66.74595642089844,'*231': 66.7466049194336,'*232': 66.7466049194336,'*233': 66.74650573730469,'*234': 66.74650573730469,'*235': 66.74677276611328,'*236': 66.74677276611328,'*237': 66.74703216552734,'*238': 66.74703216552734,'*239': 66.74669647216797,'*24': 66.75222778320312,'*240': 66.74691772460938,'*241': 66.74691772460938,'*242': 66.74691772460938,'*243': 66.74700164794922,'*244': 66.74751281738281,'*245': 66.74751281738281,'*246': 66.74751281738281,'*247': 66.74751281738281,'*248': 66.74751281738281,'*249': 66.74751281738281,'*25': 66.75190734863281,'*250': 66.7479476928711,'*251': 66.7479476928711,'*252': 66.74848175048828,'*253': 66.74842834472656,'*254': 66.74842834472656,'*255': 66.74869537353516,'*256': 66.7486801147461,'*257': 66.74786376953125,'*258': 66.74786376953125,'*259': 66.74786376953125,'*26': 66.75190734863281,'*260': 66.74786376953125,'*261': 66.74832153320312,'*262': 66.74832153320312,'*263': 66.74832153320312,'*264': 66.74832153320312,'*265': 66.74832153320312,'*266': 66.74850463867188,'*267': 66.74901580810547,'*268': 66.7490463256836,'*269': 66.7490463256836,'*27': 66.75190734863281,'*270': 66.74906921386719,'*271': 66.74942779541016,'*272': 66.74942779541016,'*273': 66.74958038330078,'*274': 66.74958038330078,'*275': 66.74958038330078,'*276': 66.74950408935547,'*277': 66.75039672851562,'*278': 66.75039672851562,'*279': 66.75029754638672,'*28': 66.75161743164062,'*280': 66.75029754638672,'*281': 66.75029754638672,'*282': 66.75029754638672,'*283': 66.75074768066406,'*284': 66.75074768066406,'*285': 66.75074768066406,'*286': 66.75074768066406,'*287': 66.75074768066406,'*288': 66.75048828125,'*289': 66.75017547607422,'*29': 66.75161743164062,'*290': 66.75017547607422,'*291': 66.75037384033203,'*292': 66.74996948242188,'*293': 66.74996948242188,'*294': 66.74996948242188,'*295': 66.74996948242188,'*296': 66.74993896484375,'*297': 66.74894714355469,'*298': 66.74934387207031,'*299': 66.74947357177734,'*3': 66.75110626220703,'*30': 66.75161743164062,'*300': 66.74984741210938,'*301': 66.74984741210938,'*302': 66.74908447265625,'*303': 66.74947357177734,'*304': 66.74947357177734,'*305': 66.7492446899414,'*306': 66.7492446899414,'*307': 66.7492446899414,'*308': 66.7492446899414,'*309': 66.7492446899414,'*31': 66.75161743164062,'*310': 66.74951934814453,'*311': 66.7489242553711,'*312': 66.74882507324219,'*313': 66.74882507324219,'*314': 66.74870300292969,'*315': 66.74870300292969,'*316': 66.74870300292969,'*317': 66.74870300292969,'*318': 66.74893188476562,'*319': 66.74893188476562,'*32': 66.75174713134766,'*320': 66.7493667602539,'*321': 66.7493667602539,'*322': 66.7493667602539,'*323': 66.74928283691406,'*324': 66.74928283691406,'*325': 66.74938201904297,'*326': 66.75033569335938,'*327': 66.75033569335938,'*328': 66.75033569335938,'*329': 66.75033569335938,'*33': 66.75174713134766,'*330': 66.75033569335938,'*331': 66.75033569335938,'*332': 66.75033569335938,'*333': 66.7507095336914,'*334': 66.7507095336914,'*335': 66.7509994506836,'*336': 66.7509994506836,'*337': 66.7509994506836,'*338': 66.75153350830078,'*339': 66.75153350830078,'*34': 66.75174713134766,'*340': 66.7516098022461,'*341': 66.75091552734375,'*342': 66.75091552734375,'*343': 66.75101470947266,'*344': 66.69735717773438,'*345': 66.697021484375,'*346': 66.697021484375,'*347': 66.697021484375,'*348': 66.697021484375,'*349': 66.697021484375,'*35': 66.75133514404297,'*350': 66.697021484375,'*351': 66.69729614257812,'*352': 66.69722747802734,'*353': 66.69722747802734,'*354': 66.69647979736328,'*355': 66.69647979736328,'*356': 66.69673156738281,'*357': 66.69673156738281,'*358': 66.69673156738281,'*359': 66.69673156738281,'*36': 66.75133514404297,'*360': 66.69673156738281,'*361': 66.696533203125,'*362': 66.696533203125,'*363': 66.696533203125,'*364': 66.69682312011719,'*365': 66.69669342041016,'*366': 66.69666290283203,'*367': 66.69666290283203,'*368': 66.69666290283203,'*369': 66.69654846191406,'*37': 66.75133514404297,'*370': 66.69695281982422,'*371': 66.69695281982422,'*372': 66.69695281982422,'*373': 66.69695281982422,'*374': 66.69695281982422,'*375': 66.69670867919922,'*376': 66.69670867919922,'*377': 66.69670867919922,'*378': 66.69670867919922,'*379': 66.69670867919922,'*38': 66.75089263916016,'*380': 66.69670867919922,'*381': 66.69670867919922,'*382': 66.69670867919922,'*383': 66.69692993164062,'*384': 66.69692993164062,'*385': 66.69683074951172,'*386': 66.69631958007812,'*387': 66.69635009765625,'*388': 66.69635009765625,'*389': 66.69596862792969,'*39': 66.75089263916016,'*390': 66.69596862792969,'*391': 66.69639587402344,'*392': 66.69639587402344,'*393': 66.69639587402344,'*394': 66.69609069824219,'*395': 66.69683074951172,'*396': 66.69683074951172,'*397': 66.69683074951172,'*398': 66.69683074951172,'*399': 66.69683074951172,'*4': 66.75110626220703,'*40': 66.75089263916016,'*400': 66.69706726074219,'*401': 66.69705200195312,'*402': 66.69705200195312,'*403': 66.69770812988281,'*404': 66.69681549072266,'*405': 66.6963882446289,'*406': 66.6963882446289,'*407': 66.69669342041016,'*408': 66.69669342041016,'*409': 66.69669342041016,'*41': 66.75089263916016,'*410': 66.6970443725586,'*411': 66.6970443725586,'*412': 66.6970443725586,'*413': 66.6970443725586,'*414': 66.69752502441406,'*415': 66.69783020019531,'*416': 66.69783020019531,'*417': 66.69783020019531,'*418': 66.69718933105469,'*419': 66.6968002319336,'*42': 66.75093078613281,'*420': 66.69713592529297,'*421': 66.69713592529297,'*422': 66.69713592529297,'*423': 66.69713592529297,'*424': 66.69749450683594,'*425': 66.69749450683594,'*426': 66.69696807861328,'*427': 66.69696807861328,'*428': 66.69696807861328,'*429': 66.69696807861328,'*43': 66.75093078613281,'*430': 66.69696807861328,'*431': 66.6966552734375,'*432': 66.6966552734375,'*433': 66.6966552734375,'*434': 66.6966552734375,'*435': 66.69532775878906,'*436': 66.69532775878906,'*437': 66.69532775878906,'*438': 66.69532775878906,'*439': 66.69532775878906,'*44': 66.75159454345703,'*440': 66.69532775878906,'*441': 66.69498443603516,'*442': 66.69498443603516,'*443': 66.69498443603516,'*444': 66.69498443603516,'*445': 66.69498443603516,'*446': 66.69498443603516,'*447': 66.69498443603516,'*448': 66.69498443603516,'*449': 66.69502258300781,'*45': 66.75163269042969,'*450': 66.69448852539062,'*451': 66.69448852539062,'*452': 66.69380950927734,'*453': 66.69401550292969,'*454': 66.69401550292969,'*455': 66.69401550292969,'*456': 66.6944351196289,'*457': 66.6944351196289,'*458': 66.69477844238281,'*459': 66.69471740722656,'*46': 66.7513427734375,'*460': 66.69471740722656,'*461': 66.69476318359375,'*462': 66.69429779052734,'*463': 66.69429779052734,'*464': 66.69429779052734,'*465': 66.69385528564453,'*466': 66.69385528564453,'*467': 66.69385528564453,'*468': 66.69385528564453,'*469': 66.69385528564453,'*47': 66.7513427734375,'*470': 66.69385528564453,'*471': 66.6942367553711,'*472': 66.6942367553711,'*473': 66.69386291503906,'*474': 66.69416809082031,'*475': 66.69458770751953,'*476': 66.69498443603516,'*477': 66.69498443603516,'*478': 66.69532012939453,'*479': 66.69474792480469,'*48': 66.7513427734375,'*480': 66.69474792480469,'*481': 66.69474792480469,'*482': 66.6944580078125,'*483': 66.6944580078125,'*484': 66.6944580078125,'*485': 66.69470977783203,'*486': 66.69470977783203,'*487': 66.69470977783203,'*488': 66.69470977783203,'*489': 66.69470977783203,'*49': 66.7513427734375,'*490': 66.69402313232422,'*491': 66.69402313232422,'*492': 66.69402313232422,'*493': 66.69402313232422,'*494': 66.69402313232422,'*495': 66.69402313232422,'*496': 66.69402313232422,'*497': 66.69302368164062,'*498': 66.69302368164062,'*499': 66.69302368164062,'*5': 66.75205993652344,'*50': 66.75162506103516,'*500': 66.69379425048828,'*501': 66.69379425048828,'*502': 66.69379425048828,'*503': 66.69446563720703,'*504': 66.69446563720703,'*505': 66.69486236572266,'*506': 66.69486236572266,'*507': 66.69486236572266,'*51': 66.75162506103516,'*52': 66.75154876708984,'*53': 66.75154876708984,'*54': 66.75154876708984,'*55': 66.75154876708984,'*56': 66.75214385986328,'*57': 66.75179290771484,'*58': 66.75106048583984,'*59': 66.75106048583984,'*6': 66.75213623046875,'*60': 66.75039672851562,'*61': 66.75074768066406,'*62': 66.75008392333984,'*63': 66.74947357177734,'*64': 66.74947357177734,'*65': 66.74947357177734,'*66': 66.74947357177734,'*67': 66.74947357177734,'*68': 66.74947357177734,'*69': 66.7491683959961,'*7': 66.75247955322266,'*70': 66.74874114990234,'*71': 66.74874114990234,'*72': 66.74837493896484,'*73': 66.74876403808594,'*74': 66.74869537353516,'*75': 66.74869537353516,'*76': 66.74827575683594,'*77': 66.74824523925781,'*78': 66.74824523925781,'*79': 66.74824523925781,'*8': 66.75247955322266,'*80': 66.74801635742188,'*81': 66.74801635742188,'*82': 66.74771118164062,'*83': 66.74738311767578,'*84': 66.74738311767578,'*85': 66.74738311767578,'*86': 66.74710845947266,'*87': 66.74710845947266,'*88': 66.74710845947266,'*89': 66.7465591430664,'*9': 66.75174713134766,'*90': 66.7465591430664,'*91': 66.7465591430664,'*92': 66.74678039550781,'*93': 66.7470703125,'*94': 66.7470703125,'*95': 66.74758911132812,'*96': 66.74655151367188,'*97': 66.74658203125,'*98': 66.74658203125,'*99': 66.74658203125}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed=self.filter_report(report)
        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, test_name, self.parallel)

        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(image=img+'.image', range_list=[0.3, 1.0])
        self.mom8_creator(image=img+'.residual', range_list=[0.3, 1.0])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')


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

# End of test_standard_cube_pcwdT
#-------------------------------------------------#
    # Test 1c
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cube_briggsbwtaper(self):
        ''' Standard (single field) cube imaging with briggsbwtaper - central field of SMIDGE_NWCloud (field 3), spw 22 '''
        
        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name,'test_')+'.iter'
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
            gridder='standard', mosweight=False, \
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
                gridder='standard', mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggsbwtaper', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
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
                mosweight=False, usepointing=False, \
                pblimit=0.2, deconvolver='hogbom', restoration=True, \
                restoringbeam='common', pbcor=True, weighting='briggsbwtaper', \
                robust=0.5, npixels=0, niter=0, threshold='0.354Jy', \
                nsigma=0.0, interactive=0, usemask='auto-multithresh', \
                sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
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
                gridder='standard', mosweight=False, \
                deconvolver='hogbom', restoration=True, pbcor=True, \
                weighting='briggsbwtaper', robust=0.5, npixels=0, niter=20000, \
                threshold='0.354Jy', interactive=0, usemask='auto'
                '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
                lownoisethreshold=2.0, negativethreshold=0.0, \
                minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
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
        exp_im_stats = {'com_bmaj': [False, 8.814554388244703],
            'com_bmin': [False, 6.202822301517538],
            'com_pa': [False, 66.69898795336522],
            'npts': [True, 3251200],
            'npts_unmasked': [False, 1522476.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0720406770706177],
            'max_val_pos': [True, [38, 36, 0, 254]],
            'min_val': [False, -0.6544804573059082],
            'min_val_pos': [True, [35, 39, 0, 249]],
            'im_rms': [False, 0.141714622981388],
            'rms_per_chan': [False, [0.12523990320953438, 0.14324924776024664, 0.12583950705148642, 0.11144521686953844, 0.14801911642978904, 0.12384414332196574, 0.13250378569694637, 0.10908336752138008, 0.15314765802689306, 0.1583436246441215, 0.11148486298886229, 0.14026536455604707, 0.16826333021708917, 0.13248216110406938, 0.15764612972434358, 0.13676761900823206, 0.10676431428877683, 0.12577388209100313, 0.14053759687358927, 0.17540610946945773, 0.1323837493527306, 0.15362016243097526, 0.11417316739635917, 0.11108196045813304, 0.16902560122983376, 0.13581396797082926, 0.1329436902855199, 0.11946738611307549, 0.14114087952870244, 0.14131605584009532, 0.14403875481300618, 0.11966224764187047, 0.1262892338131058, 0.12162446249753633, 0.15572429272778068, 0.1563711726373196, 0.15842446875968766, 0.1422062409898373, 0.1436638074808527, 0.13050965429929662, 0.13145123246119894, 0.15131515970975923, 0.14277396397877132, 0.1368859509936139, 0.1397316365701372, 0.1419720108419799, 0.13008409616201266, 0.13432683037237741, 0.16895022395199208, 0.14005692770520914, 0.12083566316690438, 0.1574960269395061, 0.1371135587512026, 0.1282966203055203, 0.11983225024629629, 0.13162758833729024, 0.1521450831821477, 0.17514796833920834, 0.10835354506445961, 0.13888487134945685, 0.1512265169967286, 0.1387244947181035, 0.12068125217468258, 0.12401348905968908, 0.12549184914108094, 0.13929440451026806, 0.12718524060890687, 0.14656878589787323, 0.1342193073823361, 0.15795224105299954, 0.12711010726110766, 0.13460207863316664, 0.14954177655117828, 0.14078922615688697, 0.17274352235994472, 0.16473779895523946, 0.11964488117889586, 0.1477568465758766, 0.15106122991701645, 0.12041895256198092, 0.12339454930558275, 0.16455778722050338, 0.11217696010101033, 0.18007400852146138, 0.13496695697534009, 0.14322603270185877, 0.11447308069071326, 0.1342779703244271, 0.15857195668142554, 0.1337513000940909, 0.1501517081622709, 0.13173220833360696, 0.13739247672443639, 0.14340668611668503, 0.12130531495372678, 0.14376342970557113, 0.1234912487725557, 0.1643869855701406, 0.14533114035000744, 0.12946880599173138, 0.15529524522203336, 0.12060794943842192, 0.14016285763408454, 0.15001192106344322, 0.1324329725847986, 0.12815597883111274, 0.15807066170120912, 0.14287343148088255, 0.15240396951008667, 0.12639275765574423, 0.14367929483550576, 0.1291083213332945, 0.12546232492096743, 0.13689215625626247, 0.13072954130611233, 0.13080680664054334, 0.11585040601884666, 0.15869033302859498, 0.12529404503652164, 0.15124305786275177, 0.13278132470748547, 0.14201522467073688, 0.13731687999210507, 0.1455160423705179, 0.12319121701562512, 0.15604767815479806, 0.14275263129456994, 0.14552597524837496, 0.12911201234548528, 0.12226225041379932, 0.14501221704501832, 0.12177996237313185, 0.15674756268839654, 0.15361337637920233, 0.15692418011617992, 0.12244514355098816, 0.15963945884914252, 0.13489339606836004, 0.13594362380091116, 0.14890052037969656, 0.11642527534791237, 0.15093293515633302, 0.14465342744892626, 0.11465862560559961, 0.160344398253983, 0.1245492611681814, 0.14199230371775157, 0.16229598256409908, 0.1303991415961471, 0.13582195028757157, 0.1420181278589793, 0.1017718380146455, 0.14792638312962866, 0.12280005351122836, 0.16460742508028833, 0.14717036550027546, 0.1485620744674475, 0.1280633312712909, 0.11602697765572954, 0.14791979960358906, 0.11853896408385961, 0.15395768056232112, 0.15383920951833827, 0.1732559063233046, 0.15552404361299874, 0.13771394778202495, 0.12123406510349648, 0.12695475263626937, 0.12842019700275953, 0.14119657578976594, 0.1358780027646155, 0.12930560092217513, 0.1404636702740045, 0.16702551828252538, 0.12589301838537403, 0.13302962329839318, 0.1852553086260136, 0.14338775061585476, 0.15138144537920686, 0.14469997182275277, 0.1294434963058186, 0.1389366231001781, 0.13714987100981843, 0.14964147606229003, 0.10399362890984527, 0.12341902258265508, 0.17360259200378875, 0.12498048031851769, 0.13065253281046377, 0.15224287270299725, 0.11804679628385578, 0.12241212406338949, 0.1483057604174402, 0.16214203232937022, 0.12528077037312707, 0.1397542648426262, 0.1322894600585039, 0.15516535922705219, 0.1683904466350801, 0.1558204781965461, 0.14846961328462063, 0.1646352907520968, 0.1337893106316923, 0.15225693182492928, 0.1623153101061988, 0.12537442473220844, 0.13764109418856924, 0.12578984056816664, 0.12674193835547487, 0.13661251291589657, 0.15399305188750886, 0.15031409544998142, 0.14383478935233593, 0.1455994215141421, 0.14939241957980293, 0.12289512904076387, 0.11933505363640073, 0.12093625713206506, 0.16774153557579063, 0.1573899933879673, 0.14721360899669808, 0.11926559262029535, 0.13021917335734925, 0.11819527600748192, 0.13230611011121265, 0.1590549212878262, 0.1889329202104889, 0.12401029402661594, 0.1277568344234904, 0.1368304026602026, 0.15968950052498893, 0.1729749248861455, 0.13360188047951554, 0.16808069770590067, 0.1431546391916005, 0.12783913900569127, 0.12068914811832127, 0.14794931992452193, 0.1584167790342406, 0.17862021952702084, 0.13674458438737186, 0.12495084893731379, 0.16282523925293957, 0.16138956122251336, 0.12888068804090588, 0.1542726922656351, 0.12653738267343279, 0.16489299411876693, 0.1569414467050826, 0.17598195544610923, 0.16653026514675875, 0.1659623044302847, 0.20115863656408775, 0.2077561987431515, 0.18802930230962378, 0.2001614410566781, 0.16767913455512257, 0.1687059108089183, 0.1371334607147074, 0.12226683295624655, 0.13166531585677832, 0.15051445470841793, 0.1359436664467732, 0.12632949023785098, 0.15301217999670708, 0.16619153901133457, 0.1581315524051311, 0.14993009238936136, 0.11827620559034412, 0.1305592374364079, 0.13777012863086804, 0.133182750736316, 0.1386089633027505, 0.14446350825720872, 0.14257403881381536, 0.15069064315613953, 0.13169317408300663, 0.14332278745942761, 0.13161273194190543, 0.15385506121807455, 0.11039828814828209, 0.12825756396321258, 0.14651762713634642, 0.12502819315155014, 0.13919164444370638, 0.11980808335953722, 0.14974049985671262, 0.11761516324645621, 0.11868470746734361, 0.13497453229650766, 0.16061827923328664, 0.15323567662060397, 0.12170167258533864, 0.1345547368132984, 0.14199950670555758, 0.13262376952254404, 0.16860258204366008, 0.14983660771171844, 0.11432993665606445, 0.10469686872030824, 0.11560205025414177, 0.1477026701946162, 0.14682798177864456, 0.14124411247434374, 0.12549436995580673, 0.13597484136362314, 0.1226443655782383, 0.12346657364310346, 0.14189038160652306, 0.16251787156110653, 0.15685849097353452, 0.13396097440784885, 0.17330263245641195, 0.16892315331393676, 0.15394018740030385, 0.12476394806751057, 0.1418872421195087, 0.13766400384141936, 0.13120391786442592, 0.14483681782405977, 0.13488309691920633, 0.15289690168979655, 0.13535832124401406, 0.13893481759699391, 0.12933859899264505, 0.1337385471001014, 0.1385609402851234, 0.12855694649194774, 0.13487289342730532, 0.1314879769034824, 0.15665437806908408, 0.12843219016026014, 0.14121932211111263, 0.15304227621161443, 0.14151261318899291, 0.16147436692227715, 0.12847114724378497, 0.11819794804965005, 0.15766687384251313, 0.1428098659865244, 0.1492665677401512, 0.16179241587473858, 0.13395707889595074, 0.11773878709128975, 0.16002199761568775, 0.1383204336230241, 0.14917819583463224, 0.12008321097679758, 0.1074060244648183, 0.13847709256527688, 0.13614390897153583, 0.12429744533612494, 0.14808673357709526, 0.12908553412849066, 0.1466775991227559, 0.15304338014946395, 0.14713159025186043, 0.12647691746085746, 0.1330103694469739, 0.13168316682008346, 0.15981577136331093, 0.1365577963382396, 0.15911118653169384, 0.14558238093646186, 0.1316812380063819, 0.12827171299841258, 0.1163013354431386, 0.12222041849683479, 0.1450612224624337, 0.13507007440044178, 0.1576027985034685, 0.12583835501269164, 0.12835162052844085, 0.1586501534733811, 0.14647146558664084, 0.12513500511868775, 0.1411048740444981, 0.13837722655241508, 0.13813318110873773, 0.1336471756974836, 0.1571461743077331, 0.16545636525496607, 0.13035278026445915, 0.13345355591347133, 0.14001400874558187, 0.1360898550799681, 0.16339444833840214, 0.14920905726172168, 0.1405336494284344, 0.15777812849965847, 0.15683966068730068, 0.13524264828698487, 0.13909648301796249, 0.13788995266747353, 0.18102666433970682, 0.1166307220292982, 0.15411615624663488, 0.11987764878583637, 0.13274984904523748, 0.16424289189597435, 0.11491788800991166, 0.13702550119433862, 0.12178506502249363, 0.16718163474852588, 0.1276879151404633, 0.1449795119126693, 0.13610686163435434, 0.1693382159203699, 0.130194030114716, 0.12472450970875938, 0.168900612656342, 0.15334396781135223, 0.1479222964111586, 0.1415849131114789, 0.14352600585338088, 0.13278140992612253, 0.132113855227346, 0.16125801095514644, 0.14807108781294695, 0.13007143629896492, 0.14203126668313262, 0.15511790688775692, 0.14149841060536295, 0.129044198647649, 0.16407838283989334, 0.1218092004422616, 0.1523737943499214, 0.11679351071579344, 0.13509227343596528, 0.1456962589461304, 0.1457619443800225, 0.12238625675859113, 0.1439047561282197, 0.13045843226723303, 0.14243637508190748, 0.1452917327769675, 0.11430486349741108, 0.1548706653614608, 0.13873713667686777, 0.1492975120624088, 0.15436076434012666, 0.14188774628035392, 0.154076630625794, 0.13867511128871202, 0.13502483548286764, 0.11686534807193642, 0.14940289012783628, 0.14150929621693897, 0.1417824035288388, 0.17918265557676552, 0.15389824235009916, 0.13870231729553045, 0.1399238896746827, 0.12127927150314834, 0.1157243646229968, 0.15152906766016336, 0.13524754961670793, 0.14909866156411875, 0.14212914142416852, 0.17294457805032637, 0.13221269092826315, 0.12189744256869132, 0.16399462773356113, 0.12091681498289851, 0.1444892723916532, 0.15141165054936703, 0.12531871850158352, 0.14665362911307916, 0.14062938671781208, 0.14214302989482006, 0.13173382976626075, 0.17359083958420152, 0.13916968008553207, 0.11976894103745454, 0.15906462695211024, 0.11184180348360881, 0.1345781030878482, 0.15557142577239674, 0.1410529965682585, 0.12106350352635133, 0.12027248175949688, 0.11890366667896228, 0.12158878002892544, 0.15759183317476963, 0.12015862680345264, 0.13505065445074577, 0.1271856793406001, 0.14408675739331528, 0.12827537422816912, 0.1347408597552298, 0.1421241641642072, 0.1416804336744793, 0.1502300827550599, 0.13395241846349626, 0.14066736926849258, 0.14652468743039632, 0.11509657653458907, 0.16272224671565466, 0.14439549171278654, 0.13000149240471043, 0.15036059536004276, 0.16554505010366696, 0.13650148067415863, 0.1388575775665499, 0.1579981350659498, 0.14879474704282886, 0.13400413195998784, 0.13508382353307088]],
            'im_sum': [False, -25.70624720742208],
            'regn_sum': [False, 73.29206447303295],
            'npts_real': [True, 3251200],
            'profile': [False, 0.9988724790368498],
            'fit': [False, [0.9971958047380823, 14.026386572107223, 7.130826791545065]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [38.04655656544654, 36.97683757400734]]}

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
            'mask_pix': [False, 925],
            'mask_regns': [True, 2],
            'npts_real': [True, 3251200]}


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_pb_stats
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
            'min_val': [False, -0.20378327369689941],
            'min_val_pos': [True, [7, 13, 0, 473]],
            'im_rms': [False, 0.13685788807110835],
            'im_sum': [False, 6927.265729894329],
            'npts_real': [True, 3251200],
            'fit_0': [False, [1.0878921738729006, 7.965647165978641, 5.42560698574255]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [39.9990504588436, 39.99287718926547]],
            'fit_1': [False, [1.087919792860718, 7.96291787231842, 5.423869803559424]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [39.99905265666077, 39.99288945679226]],
            'fit_2': [False, [1.0879808688789583, 7.959735590112945, 5.421951218801309]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [39.99905538610993, 39.992898784789176]]}

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
            'max_val': [False, 0.8372206091880798],
            'max_val_pos': [True, [41, 43, 0, 255]],
            'min_val': [False, -0.6544804573059082],
            'min_val_pos': [True, [35, 39, 0, 249]],
            'im_rms': [False, 0.14156906077438056],
            'im_sum': [False, -110.75728893821375],
            'regn_sum': [False, 33.60975038865581],
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
            'max_val': [False, 0.5203002691268921],
            'max_val_pos': [True, [42, 43, 0, 256]],
            'min_val': [False, -0.03937331587076187],
            'min_val_pos': [True, [30, 38, 0, 256]],
            'im_rms': [False, 0.0004073412175002982],
            'im_sum': [False, 1.661157950758934],
            'regn_sum': [False, 1.6986552141606808],
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
            'max_val': [False, 90.41231536865234],
            'max_val_pos': [True, [0, 0, 0, 140]],
            'min_val': [False, 90.41229248046875],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 90.4123080137647],
            'npts_real': [True, 508]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube_briggsbwtaper.exp_bmin_dict
            exp_bmin_dict = {'*0': 6.202030658721924,'*1': 6.202027320861816,'*10': 6.2019429206848145,'*100': 6.201361179351807,'*101': 6.201338768005371,'*102': 6.2013325691223145,'*103': 6.2013258934021,'*104': 6.201310634613037,'*105': 6.201310634613037,'*106': 6.201282978057861,'*107': 6.201282978057861,'*108': 6.201256275177002,'*109': 6.201256275177002,'*11': 6.201923370361328,'*110': 6.201256275177002,'*111': 6.201256275177002,'*112': 6.201256275177002,'*113': 6.201256275177002,'*114': 6.201256275177002,'*115': 6.201256275177002,'*116': 6.201256275177002,'*117': 6.201247692108154,'*118': 6.201218605041504,'*119': 6.201218605041504,'*12': 6.20189905166626,'*120': 6.201183319091797,'*121': 6.201183319091797,'*122': 6.20115852355957,'*123': 6.20112943649292,'*124': 6.201098918914795,'*125': 6.201138496398926,'*126': 6.201130390167236,'*127': 6.201130390167236,'*128': 6.201130390167236,'*129': 6.2011213302612305,'*13': 6.20189905166626,'*130': 6.2011213302612305,'*131': 6.2011213302612305,'*132': 6.2011213302612305,'*133': 6.201091766357422,'*134': 6.201091766357422,'*135': 6.201021671295166,'*136': 6.201021671295166,'*137': 6.201016426086426,'*138': 6.200996398925781,'*139': 6.200996398925781,'*14': 6.20189905166626,'*140': 6.2009782791137695,'*141': 6.2009782791137695,'*142': 6.2009782791137695,'*143': 6.2009782791137695,'*144': 6.2009782791137695,'*145': 6.2009782791137695,'*146': 6.200955867767334,'*147': 6.200955867767334,'*148': 6.200918674468994,'*149': 6.20089864730835,'*15': 6.20189905166626,'*150': 6.20089864730835,'*151': 6.20089864730835,'*152': 6.200883388519287,'*153': 6.200883388519287,'*154': 6.200864791870117,'*155': 6.200864791870117,'*156': 6.200839996337891,'*157': 6.200839996337891,'*158': 6.200839996337891,'*159': 6.200839996337891,'*16': 6.201883316040039,'*160': 6.200839996337891,'*161': 6.200765132904053,'*162': 6.200765132904053,'*163': 6.200747013092041,'*164': 6.200730323791504,'*165': 6.200721263885498,'*166': 6.200707912445068,'*167': 6.200707912445068,'*168': 6.200702667236328,'*169': 6.200702667236328,'*17': 6.201852798461914,'*170': 6.200702667236328,'*171': 6.200702667236328,'*172': 6.200685501098633,'*173': 6.200675964355469,'*174': 6.200669765472412,'*175': 6.20062780380249,'*176': 6.200601577758789,'*177': 6.200589179992676,'*178': 6.200563430786133,'*179': 6.200563430786133,'*18': 6.201840877532959,'*180': 6.200540065765381,'*181': 6.200540065765381,'*182': 6.200540065765381,'*183': 6.200534820556641,'*184': 6.200538158416748,'*185': 6.200538158416748,'*186': 6.200523376464844,'*187': 6.200523376464844,'*188': 6.200523376464844,'*189': 6.200523376464844,'*19': 6.201831340789795,'*190': 6.200523376464844,'*191': 6.200523376464844,'*192': 6.2005181312561035,'*193': 6.2005181312561035,'*194': 6.2005181312561035,'*195': 6.2005181312561035,'*196': 6.2005181312561035,'*197': 6.2005181312561035,'*198': 6.200517177581787,'*199': 6.200491905212402,'*2': 6.202004432678223,'*20': 6.201831340789795,'*200': 6.200491905212402,'*201': 6.200491905212402,'*202': 6.200467109680176,'*203': 6.200464248657227,'*204': 6.200451374053955,'*205': 6.200451374053955,'*206': 6.200451374053955,'*207': 6.200451374053955,'*208': 6.200451374053955,'*209': 6.200451374053955,'*21': 6.201831340789795,'*210': 6.200445175170898,'*211': 6.200445175170898,'*212': 6.200445175170898,'*213': 6.200445175170898,'*214': 6.200445175170898,'*215': 6.200445175170898,'*216': 6.200429439544678,'*217': 6.200429439544678,'*218': 6.2004194259643555,'*219': 6.200418949127197,'*22': 6.201831340789795,'*220': 6.200412273406982,'*221': 6.200412273406982,'*222': 6.200412273406982,'*223': 6.200412273406982,'*224': 6.20039176940918,'*225': 6.20038366317749,'*226': 6.20038366317749,'*227': 6.200374126434326,'*228': 6.200349807739258,'*229': 6.200349807739258,'*23': 6.201831340789795,'*230': 6.200349807739258,'*231': 6.200352191925049,'*232': 6.200352191925049,'*233': 6.200352668762207,'*234': 6.200352668762207,'*235': 6.200281143188477,'*236': 6.200281143188477,'*237': 6.200257778167725,'*238': 6.200257778167725,'*239': 6.200234413146973,'*24': 6.201822280883789,'*240': 6.200222015380859,'*241': 6.200222015380859,'*242': 6.200222015380859,'*243': 6.200189113616943,'*244': 6.2001423835754395,'*245': 6.2001423835754395,'*246': 6.2001423835754395,'*247': 6.2001423835754395,'*248': 6.2001423835754395,'*249': 6.2001423835754395,'*25': 6.201816082000732,'*250': 6.200117111206055,'*251': 6.200117111206055,'*252': 6.20006799697876,'*253': 6.200074672698975,'*254': 6.200074672698975,'*255': 6.200040817260742,'*256': 6.200023651123047,'*257': 6.200010776519775,'*258': 6.200010776519775,'*259': 6.200010776519775,'*26': 6.201816082000732,'*260': 6.200010776519775,'*261': 6.2000041007995605,'*262': 6.2000041007995605,'*263': 6.2000041007995605,'*264': 6.2000041007995605,'*265': 6.2000041007995605,'*266': 6.199972629547119,'*267': 6.199947834014893,'*268': 6.199952125549316,'*269': 6.199952125549316,'*27': 6.201816082000732,'*270': 6.199923038482666,'*271': 6.199893474578857,'*272': 6.199893474578857,'*273': 6.199865341186523,'*274': 6.199865341186523,'*275': 6.199865341186523,'*276': 6.199865341186523,'*277': 6.199844837188721,'*278': 6.199844837188721,'*279': 6.1998443603515625,'*28': 6.201807975769043,'*280': 6.1998443603515625,'*281': 6.1998443603515625,'*282': 6.1998443603515625,'*283': 6.199779987335205,'*284': 6.199779987335205,'*285': 6.199779987335205,'*286': 6.199779987335205,'*287': 6.199779987335205,'*288': 6.199778079986572,'*289': 6.199771881103516,'*29': 6.201807975769043,'*290': 6.199771881103516,'*291': 6.199758052825928,'*292': 6.199752330780029,'*293': 6.199752330780029,'*294': 6.199752330780029,'*295': 6.199752330780029,'*296': 6.199739456176758,'*297': 6.199734210968018,'*298': 6.199698448181152,'*299': 6.199672222137451,'*3': 6.202004432678223,'*30': 6.201807975769043,'*300': 6.199641227722168,'*301': 6.199641227722168,'*302': 6.196625232696533,'*303': 6.196559429168701,'*304': 6.196559429168701,'*305': 6.196538925170898,'*306': 6.196538925170898,'*307': 6.196538925170898,'*308': 6.196538925170898,'*309': 6.196538925170898,'*31': 6.201807975769043,'*310': 6.196509838104248,'*311': 6.196499824523926,'*312': 6.196498394012451,'*313': 6.196498394012451,'*314': 6.196496963500977,'*315': 6.196496963500977,'*316': 6.196496963500977,'*317': 6.196496963500977,'*318': 6.196484088897705,'*319': 6.196484088897705,'*32': 6.20179557800293,'*320': 6.19645881652832,'*321': 6.19645881652832,'*322': 6.19645881652832,'*323': 6.196468353271484,'*324': 6.196468353271484,'*325': 6.1964335441589355,'*326': 6.19639253616333,'*327': 6.19639253616333,'*328': 6.19639253616333,'*329': 6.19639253616333,'*33': 6.20179557800293,'*330': 6.19639253616333,'*331': 6.19639253616333,'*332': 6.19639253616333,'*333': 6.196359157562256,'*334': 6.196359157562256,'*335': 6.1963300704956055,'*336': 6.1963300704956055,'*337': 6.1963300704956055,'*338': 6.196277618408203,'*339': 6.196277141571045,'*34': 6.20179557800293,'*340': 6.196240425109863,'*341': 6.196230888366699,'*342': 6.196230888366699,'*343': 6.196174144744873,'*344': 6.1961669921875,'*345': 6.196159839630127,'*346': 6.196159839630127,'*347': 6.196159839630127,'*348': 6.196159839630127,'*349': 6.196159839630127,'*35': 6.201786041259766,'*350': 6.196159839630127,'*351': 6.196132659912109,'*352': 6.1961164474487305,'*353': 6.1961164474487305,'*354': 6.196106910705566,'*355': 6.196106910705566,'*356': 6.196094989776611,'*357': 6.196094989776611,'*358': 6.196094989776611,'*359': 6.196094989776611,'*36': 6.201786041259766,'*360': 6.196094989776611,'*361': 6.196101188659668,'*362': 6.196101188659668,'*363': 6.196101188659668,'*364': 6.1960768699646,'*365': 6.196057319641113,'*366': 6.1960368156433105,'*367': 6.1960368156433105,'*368': 6.1960368156433105,'*369': 6.196034908294678,'*37': 6.201786041259766,'*370': 6.196028232574463,'*371': 6.196028232574463,'*372': 6.196028232574463,'*373': 6.196028232574463,'*374': 6.196028232574463,'*375': 6.196022033691406,'*376': 6.196022033691406,'*377': 6.196022033691406,'*378': 6.196022033691406,'*379': 6.196022033691406,'*38': 6.201777935028076,'*380': 6.196022033691406,'*381': 6.196022033691406,'*382': 6.196022033691406,'*383': 6.196000576019287,'*384': 6.196000576019287,'*385': 6.1959991455078125,'*386': 6.19598913192749,'*387': 6.1959943771362305,'*388': 6.1959943771362305,'*389': 6.195987224578857,'*39': 6.201777935028076,'*390': 6.195987224578857,'*391': 6.1959614753723145,'*392': 6.1959614753723145,'*393': 6.1959614753723145,'*394': 6.195955276489258,'*395': 6.195930004119873,'*396': 6.195930004119873,'*397': 6.195930004119873,'*398': 6.195930004119873,'*399': 6.195930004119873,'*4': 6.202004432678223,'*40': 6.201777935028076,'*400': 6.1959147453308105,'*401': 6.195895671844482,'*402': 6.195895671844482,'*403': 6.195870399475098,'*404': 6.195854663848877,'*405': 6.195847988128662,'*406': 6.195847988128662,'*407': 6.1957902908325195,'*408': 6.1957902908325195,'*409': 6.1957902908325195,'*41': 6.201777935028076,'*410': 6.195760726928711,'*411': 6.195760726928711,'*412': 6.195760726928711,'*413': 6.195760726928711,'*414': 6.195734024047852,'*415': 6.195714950561523,'*416': 6.195714950561523,'*417': 6.195714950561523,'*418': 6.195699214935303,'*419': 6.195644855499268,'*42': 6.201782703399658,'*420': 6.195609092712402,'*421': 6.195609092712402,'*422': 6.195609092712402,'*423': 6.195609092712402,'*424': 6.195580005645752,'*425': 6.195580005645752,'*426': 6.195571422576904,'*427': 6.195571422576904,'*428': 6.195571422576904,'*429': 6.195571422576904,'*43': 6.201782703399658,'*430': 6.195571422576904,'*431': 6.195566177368164,'*432': 6.195566177368164,'*433': 6.195566177368164,'*434': 6.195566177368164,'*435': 6.195522785186768,'*436': 6.195522785186768,'*437': 6.195523262023926,'*438': 6.195523262023926,'*439': 6.195523262023926,'*44': 6.2017669677734375,'*440': 6.195523262023926,'*441': 6.195517063140869,'*442': 6.195517063140869,'*443': 6.195517063140869,'*444': 6.195517063140869,'*445': 6.195517063140869,'*446': 6.195517063140869,'*447': 6.195517063140869,'*448': 6.195517063140869,'*449': 6.195521354675293,'*45': 6.201745510101318,'*450': 6.195518493652344,'*451': 6.195518493652344,'*452': 6.195505619049072,'*453': 6.1954851150512695,'*454': 6.1954851150512695,'*455': 6.1954851150512695,'*456': 6.1954851150512695,'*457': 6.1954851150512695,'*458': 6.195457935333252,'*459': 6.19545841217041,'*46': 6.201740264892578,'*460': 6.195458889007568,'*461': 6.195443630218506,'*462': 6.195440292358398,'*463': 6.195440292358398,'*464': 6.195440292358398,'*465': 6.195432662963867,'*466': 6.195432662963867,'*467': 6.195432662963867,'*468': 6.195432662963867,'*469': 6.195433139801025,'*47': 6.201740264892578,'*470': 6.195433139801025,'*471': 6.195411682128906,'*472': 6.195411682128906,'*473': 6.195405006408691,'*474': 6.19537878036499,'*475': 6.195344924926758,'*476': 6.195309638977051,'*477': 6.195309638977051,'*478': 6.1952805519104,'*479': 6.195265769958496,'*48': 6.201740264892578,'*480': 6.195265769958496,'*481': 6.195265769958496,'*482': 6.195259094238281,'*483': 6.195259094238281,'*484': 6.195259094238281,'*485': 6.195234298706055,'*486': 6.195234298706055,'*487': 6.195234298706055,'*488': 6.195234298706055,'*489': 6.195234298706055,'*49': 6.201740264892578,'*490': 6.195224285125732,'*491': 6.195224285125732,'*492': 6.195224285125732,'*493': 6.195224285125732,'*494': 6.195224285125732,'*495': 6.195224285125732,'*496': 6.195223808288574,'*497': 6.195218563079834,'*498': 6.195218563079834,'*499': 6.195218563079834,'*5': 6.20200777053833,'*50': 6.201714515686035,'*500': 6.195166110992432,'*501': 6.195166110992432,'*502': 6.195166110992432,'*503': 6.195115566253662,'*504': 6.195115566253662,'*505': 6.195098400115967,'*506': 6.195098400115967,'*507': 6.195098400115967,'*51': 6.201714515686035,'*52': 6.201714992523193,'*53': 6.201714992523193,'*54': 6.201714992523193,'*55': 6.201714992523193,'*56': 6.201687812805176,'*57': 6.201681613922119,'*58': 6.201670169830322,'*59': 6.201670169830322,'*6': 6.201994895935059,'*60': 6.201659202575684,'*61': 6.20165491104126,'*62': 6.201644420623779,'*63': 6.201629161834717,'*64': 6.201629161834717,'*65': 6.201629161834717,'*66': 6.201629161834717,'*67': 6.201629161834717,'*68': 6.201629161834717,'*69': 6.201620578765869,'*7': 6.201960563659668,'*70': 6.201612949371338,'*71': 6.201612949371338,'*72': 6.201606273651123,'*73': 6.201574802398682,'*74': 6.201575756072998,'*75': 6.201575756072998,'*76': 6.201569080352783,'*77': 6.20154333114624,'*78': 6.20154333114624,'*79': 6.20154333114624,'*8': 6.201960563659668,'*80': 6.201538562774658,'*81': 6.201538562774658,'*82': 6.201499938964844,'*83': 6.201490879058838,'*84': 6.201490879058838,'*85': 6.201490879058838,'*86': 6.201488494873047,'*87': 6.201488494873047,'*88': 6.201488494873047,'*89': 6.201478958129883,'*9': 6.2019429206848145,'*90': 6.201478958129883,'*91': 6.201478958129883,'*92': 6.201456546783447,'*93': 6.201431751251221,'*94': 6.201431751251221,'*95': 6.201416492462158,'*96': 6.20139741897583,'*97': 6.201372146606445,'*98': 6.201372146606445,'*99': 6.201372146606445}
            # test_standard_cube_briggsbwtaper.exp_bmaj_dict
            exp_bmaj_dict = {'*0': 8.796247482299805,'*1': 8.796197891235352,'*10': 8.79611587524414,'*100': 8.794785499572754,'*101': 8.794793128967285,'*102': 8.79476261138916,'*103': 8.794717788696289,'*104': 8.794709205627441,'*105': 8.794709205627441,'*106': 8.79468059539795,'*107': 8.79468059539795,'*108': 8.794660568237305,'*109': 8.794660568237305,'*11': 8.796100616455078,'*110': 8.794660568237305,'*111': 8.794660568237305,'*112': 8.794660568237305,'*113': 8.794660568237305,'*114': 8.794660568237305,'*115': 8.794660568237305,'*116': 8.794660568237305,'*117': 8.794598579406738,'*118': 8.794584274291992,'*119': 8.794584274291992,'*12': 8.796087265014648,'*120': 8.794573783874512,'*121': 8.794573783874512,'*122': 8.7945556640625,'*123': 8.79443359375,'*124': 8.794374465942383,'*125': 8.794269561767578,'*126': 8.794212341308594,'*127': 8.794212341308594,'*128': 8.794212341308594,'*129': 8.794205665588379,'*13': 8.796087265014648,'*130': 8.794205665588379,'*131': 8.794205665588379,'*132': 8.794205665588379,'*133': 8.794194221496582,'*134': 8.794194221496582,'*135': 8.794248580932617,'*136': 8.794248580932617,'*137': 8.794230461120605,'*138': 8.794158935546875,'*139': 8.794158935546875,'*14': 8.796087265014648,'*140': 8.79417896270752,'*141': 8.79417896270752,'*142': 8.79417896270752,'*143': 8.79417896270752,'*144': 8.79417896270752,'*145': 8.79417896270752,'*146': 8.79413890838623,'*147': 8.79413890838623,'*148': 8.79410171508789,'*149': 8.794098854064941,'*15': 8.796087265014648,'*150': 8.794098854064941,'*151': 8.794098854064941,'*152': 8.79411792755127,'*153': 8.79411792755127,'*154': 8.794129371643066,'*155': 8.794129371643066,'*156': 8.794127464294434,'*157': 8.794127464294434,'*158': 8.794127464294434,'*159': 8.794127464294434,'*16': 8.796070098876953,'*160': 8.794127464294434,'*161': 8.794120788574219,'*162': 8.794120788574219,'*163': 8.794110298156738,'*164': 8.79407787322998,'*165': 8.794014930725098,'*166': 8.793999671936035,'*167': 8.793999671936035,'*168': 8.79398250579834,'*169': 8.79398250579834,'*17': 8.795960426330566,'*170': 8.79398250579834,'*171': 8.79398250579834,'*172': 8.794002532958984,'*173': 8.79393196105957,'*174': 8.79391098022461,'*175': 8.79393482208252,'*176': 8.793915748596191,'*177': 8.793904304504395,'*178': 8.793890953063965,'*179': 8.793890953063965,'*18': 8.795968055725098,'*180': 8.79387378692627,'*181': 8.79387378692627,'*182': 8.79387378692627,'*183': 8.793851852416992,'*184': 8.793813705444336,'*185': 8.793813705444336,'*186': 8.793810844421387,'*187': 8.793810844421387,'*188': 8.793810844421387,'*189': 8.793810844421387,'*19': 8.795896530151367,'*190': 8.793810844421387,'*191': 8.793810844421387,'*192': 8.793768882751465,'*193': 8.793768882751465,'*194': 8.793768882751465,'*195': 8.793768882751465,'*196': 8.793768882751465,'*197': 8.793768882751465,'*198': 8.793739318847656,'*199': 8.793720245361328,'*2': 8.796226501464844,'*20': 8.795896530151367,'*200': 8.793720245361328,'*201': 8.793720245361328,'*202': 8.793702125549316,'*203': 8.793644905090332,'*204': 8.79365348815918,'*205': 8.79365348815918,'*206': 8.79365348815918,'*207': 8.79365348815918,'*208': 8.79365348815918,'*209': 8.79365348815918,'*21': 8.795896530151367,'*210': 8.793585777282715,'*211': 8.793585777282715,'*212': 8.793585777282715,'*213': 8.793585777282715,'*214': 8.793585777282715,'*215': 8.793585777282715,'*216': 8.793559074401855,'*217': 8.793559074401855,'*218': 8.793488502502441,'*219': 8.793488502502441,'*22': 8.795896530151367,'*220': 8.793455123901367,'*221': 8.793455123901367,'*222': 8.793455123901367,'*223': 8.793455123901367,'*224': 8.793476104736328,'*225': 8.79341983795166,'*226': 8.79341983795166,'*227': 8.793378829956055,'*228': 8.793365478515625,'*229': 8.793366432189941,'*23': 8.795896530151367,'*230': 8.793366432189941,'*231': 8.793330192565918,'*232': 8.793330192565918,'*233': 8.793269157409668,'*234': 8.793269157409668,'*235': 8.793323516845703,'*236': 8.793323516845703,'*237': 8.793305397033691,'*238': 8.793305397033691,'*239': 8.793269157409668,'*24': 8.795827865600586,'*240': 8.793258666992188,'*241': 8.793258666992188,'*242': 8.793258666992188,'*243': 8.793275833129883,'*244': 8.793251037597656,'*245': 8.793251037597656,'*246': 8.793251037597656,'*247': 8.793251037597656,'*248': 8.793251037597656,'*249': 8.793251037597656,'*25': 8.795814514160156,'*250': 8.793231964111328,'*251': 8.793231964111328,'*252': 8.793196678161621,'*253': 8.793163299560547,'*254': 8.793163299560547,'*255': 8.793144226074219,'*256': 8.793155670166016,'*257': 8.79310417175293,'*258': 8.79310417175293,'*259': 8.79310417175293,'*26': 8.795814514160156,'*260': 8.79310417175293,'*261': 8.79308032989502,'*262': 8.79308032989502,'*263': 8.79308032989502,'*264': 8.79308032989502,'*265': 8.79308032989502,'*266': 8.793098449707031,'*267': 8.793015480041504,'*268': 8.792970657348633,'*269': 8.792970657348633,'*27': 8.795814514160156,'*270': 8.792988777160645,'*271': 8.792977333068848,'*272': 8.792977333068848,'*273': 8.79299545288086,'*274': 8.79299545288086,'*275': 8.79299545288086,'*276': 8.792950630187988,'*277': 8.792915344238281,'*278': 8.792915344238281,'*279': 8.792886734008789,'*28': 8.795757293701172,'*280': 8.792886734008789,'*281': 8.792886734008789,'*282': 8.792886734008789,'*283': 8.792896270751953,'*284': 8.792896270751953,'*285': 8.792896270751953,'*286': 8.792896270751953,'*287': 8.792896270751953,'*288': 8.792807579040527,'*289': 8.792762756347656,'*29': 8.795757293701172,'*290': 8.792762756347656,'*291': 8.792756080627441,'*292': 8.79262924194336,'*293': 8.79262924194336,'*294': 8.79262924194336,'*295': 8.79262924194336,'*296': 8.792637825012207,'*297': 8.792531967163086,'*298': 8.792521476745605,'*299': 8.792529106140137,'*3': 8.796226501464844,'*30': 8.795757293701172,'*300': 8.792518615722656,'*301': 8.792518615722656,'*302': 8.805731773376465,'*303': 8.805745124816895,'*304': 8.805745124816895,'*305': 8.80565357208252,'*306': 8.80565357208252,'*307': 8.80565357208252,'*308': 8.805654525756836,'*309': 8.805654525756836,'*31': 8.795757293701172,'*310': 8.805654525756836,'*311': 8.805631637573242,'*312': 8.805599212646484,'*313': 8.805599212646484,'*314': 8.805563926696777,'*315': 8.805563926696777,'*316': 8.805563926696777,'*317': 8.805563926696777,'*318': 8.805554389953613,'*319': 8.805554389953613,'*32': 8.795754432678223,'*320': 8.805535316467285,'*321': 8.805535316467285,'*322': 8.805535316467285,'*323': 8.8054838180542,'*324': 8.8054838180542,'*325': 8.805500030517578,'*326': 8.805460929870605,'*327': 8.805460929870605,'*328': 8.805460929870605,'*329': 8.805460929870605,'*33': 8.795754432678223,'*330': 8.805460929870605,'*331': 8.805460929870605,'*332': 8.805460929870605,'*333': 8.805451393127441,'*334': 8.805451393127441,'*335': 8.805450439453125,'*336': 8.805450439453125,'*337': 8.805450439453125,'*338': 8.805444717407227,'*339': 8.80544376373291,'*34': 8.795754432678223,'*340': 8.805462837219238,'*341': 8.805451393127441,'*342': 8.805451393127441,'*343': 8.805376052856445,'*344': 8.805334091186523,'*345': 8.805285453796387,'*346': 8.805285453796387,'*347': 8.805285453796387,'*348': 8.805285453796387,'*349': 8.805285453796387,'*35': 8.79576587677002,'*350': 8.805285453796387,'*351': 8.805283546447754,'*352': 8.805183410644531,'*353': 8.805184364318848,'*354': 8.805069923400879,'*355': 8.805069923400879,'*356': 8.805058479309082,'*357': 8.805058479309082,'*358': 8.805058479309082,'*359': 8.805058479309082,'*36': 8.79576587677002,'*360': 8.805058479309082,'*361': 8.804980278015137,'*362': 8.804980278015137,'*363': 8.804981231689453,'*364': 8.804967880249023,'*365': 8.804991722106934,'*366': 8.805006980895996,'*367': 8.805006980895996,'*368': 8.805006980895996,'*369': 8.804971694946289,'*37': 8.79576587677002,'*370': 8.804949760437012,'*371': 8.804949760437012,'*372': 8.804949760437012,'*373': 8.804949760437012,'*374': 8.804949760437012,'*375': 8.80488395690918,'*376': 8.80488395690918,'*377': 8.80488395690918,'*378': 8.80488395690918,'*379': 8.80488395690918,'*38': 8.795722007751465,'*380': 8.80488395690918,'*381': 8.80488395690918,'*382': 8.80488395690918,'*383': 8.80488109588623,'*384': 8.80488109588623,'*385': 8.804848670959473,'*386': 8.80482292175293,'*387': 8.804755210876465,'*388': 8.804755210876465,'*389': 8.804702758789062,'*39': 8.795722007751465,'*390': 8.804702758789062,'*391': 8.80468463897705,'*392': 8.80468463897705,'*393': 8.80468463897705,'*394': 8.80463695526123,'*395': 8.804607391357422,'*396': 8.804607391357422,'*397': 8.804607391357422,'*398': 8.804607391357422,'*399': 8.804607391357422,'*4': 8.796226501464844,'*40': 8.795722007751465,'*400': 8.804597854614258,'*401': 8.804611206054688,'*402': 8.804611206054688,'*403': 8.804579734802246,'*404': 8.804519653320312,'*405': 8.804468154907227,'*406': 8.804468154907227,'*407': 8.804494857788086,'*408': 8.804494857788086,'*409': 8.804494857788086,'*41': 8.795722007751465,'*410': 8.804483413696289,'*411': 8.804483413696289,'*412': 8.804483413696289,'*413': 8.804483413696289,'*414': 8.804465293884277,'*415': 8.80445384979248,'*416': 8.80445384979248,'*417': 8.80445384979248,'*418': 8.804384231567383,'*419': 8.804231643676758,'*42': 8.795660972595215,'*420': 8.804238319396973,'*421': 8.804238319396973,'*422': 8.804238319396973,'*423': 8.804238319396973,'*424': 8.804182052612305,'*425': 8.804182052612305,'*426': 8.804177284240723,'*427': 8.804177284240723,'*428': 8.804177284240723,'*429': 8.804177284240723,'*43': 8.795660972595215,'*430': 8.804177284240723,'*431': 8.80413818359375,'*432': 8.80413818359375,'*433': 8.80413818359375,'*434': 8.80413818359375,'*435': 8.804149627685547,'*436': 8.804149627685547,'*437': 8.804149627685547,'*438': 8.804149627685547,'*439': 8.804149627685547,'*44': 8.795635223388672,'*440': 8.804149627685547,'*441': 8.80410385131836,'*442': 8.80410385131836,'*443': 8.80410385131836,'*444': 8.80410385131836,'*445': 8.80410385131836,'*446': 8.80410385131836,'*447': 8.80410385131836,'*448': 8.80410385131836,'*449': 8.804048538208008,'*45': 8.79564380645752,'*450': 8.80400276184082,'*451': 8.80400276184082,'*452': 8.803956031799316,'*453': 8.803953170776367,'*454': 8.803953170776367,'*455': 8.803953170776367,'*456': 8.803927421569824,'*457': 8.803927421569824,'*458': 8.803877830505371,'*459': 8.803825378417969,'*46': 8.795600891113281,'*460': 8.803825378417969,'*461': 8.803831100463867,'*462': 8.803791999816895,'*463': 8.803791999816895,'*464': 8.803791999816895,'*465': 8.803772926330566,'*466': 8.803772926330566,'*467': 8.803772926330566,'*468': 8.803772926330566,'*469': 8.803772926330566,'*47': 8.795600891113281,'*470': 8.803772926330566,'*471': 8.803756713867188,'*472': 8.803756713867188,'*473': 8.803728103637695,'*474': 8.803715705871582,'*475': 8.80370044708252,'*476': 8.803690910339355,'*477': 8.803690910339355,'*478': 8.803678512573242,'*479': 8.803696632385254,'*48': 8.795600891113281,'*480': 8.803696632385254,'*481': 8.803696632385254,'*482': 8.803621292114258,'*483': 8.803621292114258,'*484': 8.803621292114258,'*485': 8.803618431091309,'*486': 8.803618431091309,'*487': 8.803618431091309,'*488': 8.803618431091309,'*489': 8.803618431091309,'*49': 8.795600891113281,'*490': 8.803606986999512,'*491': 8.803606986999512,'*492': 8.803606986999512,'*493': 8.803606986999512,'*494': 8.803606986999512,'*495': 8.803606986999512,'*496': 8.803606033325195,'*497': 8.80350399017334,'*498': 8.80350399017334,'*499': 8.80350399017334,'*5': 8.796172142028809,'*50': 8.79558277130127,'*500': 8.803470611572266,'*501': 8.803470611572266,'*502': 8.803470611572266,'*503': 8.803450584411621,'*504': 8.803450584411621,'*505': 8.80343246459961,'*506': 8.80343246459961,'*507': 8.80343246459961,'*51': 8.79558277130127,'*52': 8.795537948608398,'*53': 8.795537948608398,'*54': 8.795537948608398,'*55': 8.795537948608398,'*56': 8.795517921447754,'*57': 8.795492172241211,'*58': 8.795446395874023,'*59': 8.795446395874023,'*6': 8.796101570129395,'*60': 8.795408248901367,'*61': 8.795389175415039,'*62': 8.795380592346191,'*63': 8.795400619506836,'*64': 8.795400619506836,'*65': 8.795400619506836,'*66': 8.795400619506836,'*67': 8.795400619506836,'*68': 8.795400619506836,'*69': 8.795337677001953,'*7': 8.796096801757812,'*70': 8.795294761657715,'*71': 8.795294761657715,'*72': 8.795269966125488,'*73': 8.795256614685059,'*74': 8.795206069946289,'*75': 8.795206069946289,'*76': 8.795154571533203,'*77': 8.795106887817383,'*78': 8.795106887817383,'*79': 8.795106887817383,'*8': 8.796096801757812,'*80': 8.795111656188965,'*81': 8.795111656188965,'*82': 8.79509162902832,'*83': 8.795023918151855,'*84': 8.795023918151855,'*85': 8.795023918151855,'*86': 8.79493236541748,'*87': 8.79493236541748,'*88': 8.79493236541748,'*89': 8.79491138458252,'*9': 8.79611587524414,'*90': 8.79491138458252,'*91': 8.79491138458252,'*92': 8.794909477233887,'*93': 8.794898986816406,'*94': 8.794898986816406,'*95': 8.794875144958496,'*96': 8.794791221618652,'*97': 8.794803619384766,'*98': 8.794803619384766,'*99': 8.794803619384766}
            # test_standard_cube_briggsbwtaper.exp_pa_dict
            exp_pa_dict = {'*0': 66.7549057006836,'*1': 66.75436401367188,'*10': 66.75486755371094,'*100': 66.75007629394531,'*101': 66.75016021728516,'*102': 66.74978637695312,'*103': 66.74946594238281,'*104': 66.74970245361328,'*105': 66.74970245361328,'*106': 66.74949645996094,'*107': 66.74949645996094,'*108': 66.7497787475586,'*109': 66.7497787475586,'*11': 66.75521087646484,'*110': 66.7497787475586,'*111': 66.7497787475586,'*112': 66.7497787475586,'*113': 66.7497787475586,'*114': 66.7497787475586,'*115': 66.7497787475586,'*116': 66.7497787475586,'*117': 66.74947357177734,'*118': 66.74983978271484,'*119': 66.74983978271484,'*12': 66.75560760498047,'*120': 66.750244140625,'*121': 66.750244140625,'*122': 66.75050354003906,'*123': 66.74972534179688,'*124': 66.74969482421875,'*125': 66.75027465820312,'*126': 66.74998474121094,'*127': 66.74998474121094,'*128': 66.74998474121094,'*129': 66.7493896484375,'*13': 66.75560760498047,'*130': 66.7493896484375,'*131': 66.7493896484375,'*132': 66.7493896484375,'*133': 66.74974822998047,'*134': 66.74974822998047,'*135': 66.7499771118164,'*136': 66.7499771118164,'*137': 66.74969482421875,'*138': 66.74993133544922,'*139': 66.74993133544922,'*14': 66.75560760498047,'*140': 66.74977111816406,'*141': 66.74977111816406,'*142': 66.74977111816406,'*143': 66.74977111816406,'*144': 66.74977111816406,'*145': 66.74977111816406,'*146': 66.74946594238281,'*147': 66.74946594238281,'*148': 66.74948120117188,'*149': 66.74968719482422,'*15': 66.75560760498047,'*150': 66.74968719482422,'*151': 66.74968719482422,'*152': 66.74905395507812,'*153': 66.74905395507812,'*154': 66.74905395507812,'*155': 66.74905395507812,'*156': 66.74929809570312,'*157': 66.74929809570312,'*158': 66.74929809570312,'*159': 66.74929809570312,'*16': 66.75597381591797,'*160': 66.74929809570312,'*161': 66.74969482421875,'*162': 66.74969482421875,'*163': 66.75,'*164': 66.75080108642578,'*165': 66.75050354003906,'*166': 66.75082397460938,'*167': 66.75082397460938,'*168': 66.75054168701172,'*169': 66.75054168701172,'*17': 66.75606536865234,'*170': 66.75054168701172,'*171': 66.75054168701172,'*172': 66.74984741210938,'*173': 66.7494888305664,'*174': 66.74918365478516,'*175': 66.74856567382812,'*176': 66.74885559082031,'*177': 66.74911499023438,'*178': 66.7494125366211,'*179': 66.7494125366211,'*18': 66.75603485107422,'*180': 66.74967193603516,'*181': 66.74967193603516,'*182': 66.74967193603516,'*183': 66.74938201904297,'*184': 66.74942016601562,'*185': 66.74942016601562,'*186': 66.74957275390625,'*187': 66.74957275390625,'*188': 66.74957275390625,'*189': 66.74957275390625,'*19': 66.75568389892578,'*190': 66.74957275390625,'*191': 66.74957275390625,'*192': 66.74932098388672,'*193': 66.74932098388672,'*194': 66.74932098388672,'*195': 66.74932098388672,'*196': 66.74932098388672,'*197': 66.74932098388672,'*198': 66.74922943115234,'*199': 66.74949645996094,'*2': 66.75421905517578,'*20': 66.75568389892578,'*200': 66.74949645996094,'*201': 66.74949645996094,'*202': 66.74976348876953,'*203': 66.74916076660156,'*204': 66.7491455078125,'*205': 66.7491455078125,'*206': 66.7491455078125,'*207': 66.7491455078125,'*208': 66.7491455078125,'*209': 66.7491455078125,'*21': 66.75568389892578,'*210': 66.74889373779297,'*211': 66.7489013671875,'*212': 66.7489013671875,'*213': 66.7489013671875,'*214': 66.7489013671875,'*215': 66.7489013671875,'*216': 66.74946594238281,'*217': 66.74946594238281,'*218': 66.75018310546875,'*219': 66.75018310546875,'*22': 66.75568389892578,'*220': 66.74977111816406,'*221': 66.74977111816406,'*222': 66.74977111816406,'*223': 66.74977111816406,'*224': 66.74968719482422,'*225': 66.74939727783203,'*226': 66.74939727783203,'*227': 66.74878692626953,'*228': 66.74908447265625,'*229': 66.74908447265625,'*23': 66.75568389892578,'*230': 66.74908447265625,'*231': 66.74972534179688,'*232': 66.74972534179688,'*233': 66.7496337890625,'*234': 66.7496337890625,'*235': 66.7499008178711,'*236': 66.7499008178711,'*237': 66.75015258789062,'*238': 66.75015258789062,'*239': 66.74981689453125,'*24': 66.7553482055664,'*240': 66.75004577636719,'*241': 66.75004577636719,'*242': 66.75004577636719,'*243': 66.7501220703125,'*244': 66.75064086914062,'*245': 66.75064086914062,'*246': 66.75064086914062,'*247': 66.75064086914062,'*248': 66.75064086914062,'*249': 66.75064086914062,'*25': 66.7550277709961,'*250': 66.75106811523438,'*251': 66.75106811523438,'*252': 66.75160217285156,'*253': 66.75155639648438,'*254': 66.75155639648438,'*255': 66.75181579589844,'*256': 66.75181579589844,'*257': 66.75099182128906,'*258': 66.75099182128906,'*259': 66.75099182128906,'*26': 66.7550277709961,'*260': 66.75099182128906,'*261': 66.7514419555664,'*262': 66.7514419555664,'*263': 66.7514419555664,'*264': 66.7514419555664,'*265': 66.7514419555664,'*266': 66.75162506103516,'*267': 66.75213623046875,'*268': 66.75216674804688,'*269': 66.75216674804688,'*27': 66.7550277709961,'*270': 66.752197265625,'*271': 66.75255584716797,'*272': 66.75255584716797,'*273': 66.75270080566406,'*274': 66.75270080566406,'*275': 66.7527084350586,'*276': 66.75263214111328,'*277': 66.7535171508789,'*278': 66.7535171508789,'*279': 66.75341796875,'*28': 66.75473022460938,'*280': 66.75341796875,'*281': 66.75341796875,'*282': 66.75341796875,'*283': 66.75386810302734,'*284': 66.75386810302734,'*285': 66.75386810302734,'*286': 66.75386810302734,'*287': 66.75386810302734,'*288': 66.75361633300781,'*289': 66.75330352783203,'*29': 66.75473022460938,'*290': 66.75330352783203,'*291': 66.75349426269531,'*292': 66.75308990478516,'*293': 66.75308990478516,'*294': 66.75308990478516,'*295': 66.75308990478516,'*296': 66.75306701660156,'*297': 66.75206756591797,'*298': 66.7524642944336,'*299': 66.75260162353516,'*3': 66.75421905517578,'*30': 66.75473022460938,'*300': 66.75296783447266,'*301': 66.75296783447266,'*302': 66.69892883300781,'*303': 66.6993179321289,'*304': 66.6993179321289,'*305': 66.69908905029297,'*306': 66.69908905029297,'*307': 66.69908905029297,'*308': 66.69908905029297,'*309': 66.69908905029297,'*31': 66.75473022460938,'*310': 66.6993637084961,'*311': 66.69876861572266,'*312': 66.69866943359375,'*313': 66.69866943359375,'*314': 66.69854736328125,'*315': 66.69854736328125,'*316': 66.69854736328125,'*317': 66.69854736328125,'*318': 66.69877624511719,'*319': 66.69877624511719,'*32': 66.75486755371094,'*320': 66.69921112060547,'*321': 66.69921112060547,'*322': 66.69921112060547,'*323': 66.69913482666016,'*324': 66.69913482666016,'*325': 66.6992416381836,'*326': 66.70018005371094,'*327': 66.70018005371094,'*328': 66.70018005371094,'*329': 66.70018005371094,'*33': 66.75486755371094,'*330': 66.70018005371094,'*331': 66.70018005371094,'*332': 66.70018005371094,'*333': 66.70055389404297,'*334': 66.70055389404297,'*335': 66.70084381103516,'*336': 66.70084381103516,'*337': 66.70084381103516,'*338': 66.70137786865234,'*339': 66.70137786865234,'*34': 66.75486755371094,'*340': 66.70146179199219,'*341': 66.70076751708984,'*342': 66.70076751708984,'*343': 66.70085906982422,'*344': 66.7004623413086,'*345': 66.70011901855469,'*346': 66.70011901855469,'*347': 66.70011901855469,'*348': 66.70011901855469,'*349': 66.70011901855469,'*35': 66.75445556640625,'*350': 66.70011901855469,'*351': 66.70039367675781,'*352': 66.70033264160156,'*353': 66.70033264160156,'*354': 66.6995849609375,'*355': 66.6995849609375,'*356': 66.69983673095703,'*357': 66.69984436035156,'*358': 66.69984436035156,'*359': 66.69984436035156,'*36': 66.75445556640625,'*360': 66.69984436035156,'*361': 66.69963836669922,'*362': 66.69963836669922,'*363': 66.69963836669922,'*364': 66.6999282836914,'*365': 66.69979858398438,'*366': 66.69976806640625,'*367': 66.69976806640625,'*368': 66.69976806640625,'*369': 66.69964599609375,'*37': 66.75445556640625,'*370': 66.70005798339844,'*371': 66.70005798339844,'*372': 66.70005798339844,'*373': 66.70005798339844,'*374': 66.70005798339844,'*375': 66.69981384277344,'*376': 66.69981384277344,'*377': 66.69981384277344,'*378': 66.69981384277344,'*379': 66.69981384277344,'*38': 66.75402069091797,'*380': 66.69981384277344,'*381': 66.69981384277344,'*382': 66.69981384277344,'*383': 66.70003509521484,'*384': 66.70003509521484,'*385': 66.69993591308594,'*386': 66.69942474365234,'*387': 66.69945526123047,'*388': 66.69945526123047,'*389': 66.6990737915039,'*39': 66.75402069091797,'*390': 66.6990737915039,'*391': 66.69950103759766,'*392': 66.69950103759766,'*393': 66.69950103759766,'*394': 66.69918823242188,'*395': 66.69993591308594,'*396': 66.69993591308594,'*397': 66.69993591308594,'*398': 66.69993591308594,'*399': 66.69993591308594,'*4': 66.75421905517578,'*40': 66.75402069091797,'*400': 66.7001724243164,'*401': 66.70014953613281,'*402': 66.70014953613281,'*403': 66.7008056640625,'*404': 66.69992065429688,'*405': 66.69949340820312,'*406': 66.69949340820312,'*407': 66.6998062133789,'*408': 66.6998062133789,'*409': 66.6998062133789,'*41': 66.75402069091797,'*410': 66.70014953613281,'*411': 66.70014953613281,'*412': 66.70014953613281,'*413': 66.70014953613281,'*414': 66.70063018798828,'*415': 66.70093536376953,'*416': 66.70093536376953,'*417': 66.70093536376953,'*418': 66.70028686523438,'*419': 66.69990539550781,'*42': 66.75405883789062,'*420': 66.70024108886719,'*421': 66.70024108886719,'*422': 66.70024108886719,'*423': 66.70024108886719,'*424': 66.70059204101562,'*425': 66.70059204101562,'*426': 66.70006561279297,'*427': 66.70006561279297,'*428': 66.70006561279297,'*429': 66.70006561279297,'*43': 66.75405883789062,'*430': 66.70006561279297,'*431': 66.69976043701172,'*432': 66.69976043701172,'*433': 66.69976043701172,'*434': 66.69976043701172,'*435': 66.69843292236328,'*436': 66.69843292236328,'*437': 66.69843292236328,'*438': 66.69843292236328,'*439': 66.69843292236328,'*44': 66.75470733642578,'*440': 66.69843292236328,'*441': 66.69808959960938,'*442': 66.69808959960938,'*443': 66.69808959960938,'*444': 66.69808959960938,'*445': 66.69808959960938,'*446': 66.69808959960938,'*447': 66.69808959960938,'*448': 66.69808959960938,'*449': 66.6981201171875,'*45': 66.75475311279297,'*450': 66.69759368896484,'*451': 66.69759368896484,'*452': 66.69691467285156,'*453': 66.6971206665039,'*454': 66.6971206665039,'*455': 66.6971206665039,'*456': 66.69754028320312,'*457': 66.69754028320312,'*458': 66.69789123535156,'*459': 66.69782257080078,'*46': 66.75447082519531,'*460': 66.69782257080078,'*461': 66.69786834716797,'*462': 66.6974105834961,'*463': 66.6974105834961,'*464': 66.6974105834961,'*465': 66.69696044921875,'*466': 66.69696044921875,'*467': 66.69696044921875,'*468': 66.69696044921875,'*469': 66.69696044921875,'*47': 66.75447082519531,'*470': 66.69696044921875,'*471': 66.69734191894531,'*472': 66.69734191894531,'*473': 66.69696807861328,'*474': 66.69727325439453,'*475': 66.69769287109375,'*476': 66.69808959960938,'*477': 66.69808959960938,'*478': 66.69841766357422,'*479': 66.6978530883789,'*48': 66.75447082519531,'*480': 66.6978530883789,'*481': 66.6978530883789,'*482': 66.69756317138672,'*483': 66.69756317138672,'*484': 66.69756317138672,'*485': 66.69781494140625,'*486': 66.69781494140625,'*487': 66.69781494140625,'*488': 66.69781494140625,'*489': 66.69781494140625,'*49': 66.75447082519531,'*490': 66.6971206665039,'*491': 66.6971206665039,'*492': 66.6971206665039,'*493': 66.6971206665039,'*494': 66.6971206665039,'*495': 66.6971206665039,'*496': 66.6971206665039,'*497': 66.69612121582031,'*498': 66.69612121582031,'*499': 66.69612121582031,'*5': 66.75518035888672,'*50': 66.75474548339844,'*500': 66.6968994140625,'*501': 66.6968994140625,'*502': 66.6968994140625,'*503': 66.69757080078125,'*504': 66.69757080078125,'*505': 66.69796752929688,'*506': 66.69796752929688,'*507': 66.69796752929688,'*51': 66.75474548339844,'*52': 66.75466918945312,'*53': 66.75466918945312,'*54': 66.75466918945312,'*55': 66.75466918945312,'*56': 66.75526428222656,'*57': 66.75492095947266,'*58': 66.75418090820312,'*59': 66.75418090820312,'*6': 66.75525665283203,'*60': 66.7535171508789,'*61': 66.75386810302734,'*62': 66.75320434570312,'*63': 66.75259399414062,'*64': 66.75259399414062,'*65': 66.75259399414062,'*66': 66.75259399414062,'*67': 66.75259399414062,'*68': 66.75259399414062,'*69': 66.75228881835938,'*7': 66.75559997558594,'*70': 66.75186920166016,'*71': 66.75186920166016,'*72': 66.75150299072266,'*73': 66.75187683105469,'*74': 66.75181579589844,'*75': 66.75181579589844,'*76': 66.75139617919922,'*77': 66.7513656616211,'*78': 66.7513656616211,'*79': 66.7513656616211,'*8': 66.75559997558594,'*80': 66.75113677978516,'*81': 66.75113677978516,'*82': 66.75083923339844,'*83': 66.75049591064453,'*84': 66.75049591064453,'*85': 66.75049591064453,'*86': 66.75023651123047,'*87': 66.75023651123047,'*88': 66.75023651123047,'*89': 66.74967956542969,'*9': 66.75486755371094,'*90': 66.74967956542969,'*91': 66.74967956542969,'*92': 66.7499008178711,'*93': 66.75019836425781,'*94': 66.75019836425781,'*95': 66.7507095336914,'*96': 66.74967193603516,'*97': 66.74970245361328,'*98': 66.74970245361328,'*99': 66.74970245361328}


            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed=self.filter_report(report)
        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, test_name, self.parallel)

        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(image=img+'.image', range_list=[0.3, 1.0])
        self.mom8_creator(image=img+'.residual', range_list=[0.3, 1.0])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')


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
    # Test 2
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
            perchanweightdensity=False, gridder='standard', \
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
            perchanweightdensity=False, gridder='standard', \
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
    # Test 3
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
            perchanweightdensity=False, gridder='standard', \
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
            perchanweightdensity=False, gridder='standard', \
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
    # Test 4a
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
            gridder='standard',  mosweight=False, \
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
            gridder='standard',  mosweight=False, \
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
    # Test 4b
    @stats_dict(test_dict)
    def test_standard_cube_eph_pcwdT(self):
        ''' Single field multi-EB ephemeris cube imaging with pcwdT - field 21PGiacobini-Zinner, spw 20 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
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
            gridder='standard', mosweight=False, \
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
            nchan=1000, start=1550, width=1, perchanweightdensity=True, \
            gridder='standard', mosweight=False, \
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

        # test_standard_cube_eph_pcwdT.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 4.555337133782812],
            'com_bmin': [False, 3.366170346656177],
            'com_pa': [False, 84.29991156839789],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 3.3081090450286865],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.37456169724464417],
            'min_val_pos': [True, [32, 47, 0, 491]],
            'im_rms': [False, 0.056959688262383505],
            'rms_per_chan': [False, [0.05117105242006511, 0.04910691216790692, 0.047426932102982156, 0.04594767095362991, 0.041504941815894876, 0.042348547290818675, 0.04466488810812934, 0.04397778192383839, 0.04495588560268184, 0.046872868313580854, 0.04510881392600156, 0.045145926959446615, 0.04341463331588715, 0.04305482906444181, 0.04212200609655667, 0.04635639230513644, 0.044076785702013815, 0.04311240874610752, 0.040855518536251786, 0.04423487535406727, 0.04787651698491268, 0.04677416291092365, 0.04848685078860121, 0.04517068043777972, 0.04537073358397038, 0.04403804170919719, 0.040989529182403055, 0.0484494984650129, 0.04396335844489689, 0.04387496716514192, 0.047806992944944865, 0.04666399854985885, 0.04777174072652818, 0.04694321756320211, 0.04743551806011033, 0.04825807780742468, 0.047894916312616186, 0.048787605009964444, 0.046723568372120036, 0.041415014001875354, 0.04030057620231381, 0.04260151752797399, 0.04547708367823241, 0.04291947350526876, 0.040855695067557145, 0.0413084679617386, 0.040466562185853844, 0.04260971134211316, 0.045372865950864566, 0.046321368931081214, 0.04299562858119001, 0.042574993927384494, 0.04546885694262151, 0.05022916945596966, 0.046300435845711654, 0.04329172625923692, 0.04765233578300124, 0.04854790061062798, 0.04506885480649343, 0.04191620395550035, 0.04192980176533679, 0.04074019127993974, 0.04694706139903479, 0.049882967045744506, 0.04766895462125716, 0.04709465837436243, 0.047402657285253794, 0.04497870462090571, 0.0458231779936337, 0.04543864517471242, 0.04268241742399205, 0.0457298390292671, 0.05067584463903453, 0.05136541554257022, 0.04948533602791158, 0.04533129542533506, 0.04585308303515753, 0.04620793156157702, 0.04468226554877215, 0.048261744257699926, 0.047887355767261996, 0.044980406688003016, 0.04249627628376679, 0.04134999933646306, 0.04463547839790661, 0.04557342151112838, 0.04549918062582748, 0.04434763217313727, 0.04654122959669387, 0.045892138858735415, 0.043338310860292356, 0.04243060373610809, 0.04343286908840588, 0.041763305635470434, 0.04273259972215914, 0.045599551555544765, 0.045885998511596, 0.05178195338382946, 0.054374353803415515, 0.045646703279901925, 0.0422111016560355, 0.04169373200503522, 0.045028273186962044, 0.04857388661025543, 0.04719437748312007, 0.046572139189858926, 0.048347379470069715, 0.04757372570672316, 0.04504940661965529, 0.04917427524889455, 0.047447815864952336, 0.04012891537298516, 0.03679587439411312, 0.040374357313197606, 0.049041228250044686, 0.04670787205248779, 0.046289341158637955, 0.044571417053124306, 0.04890685929840605, 0.047634872436643425, 0.04588840940858705, 0.04640668204687149, 0.0449550644882811, 0.043537829280317145, 0.043625133823568994, 0.04506651595944801, 0.047697259834801484, 0.04636071162880931, 0.043861149846464005, 0.04555436979659354, 0.04460252331685553, 0.04458384027228497, 0.04858149242247695, 0.0513023316282301, 0.04844088996855067, 0.042354599333352075, 0.0461534126629484, 0.04358293665038448, 0.04049576757405801, 0.04047208371559128, 0.042466880839238275, 0.045152784462721134, 0.04426044590420293, 0.043123110799153964, 0.040737362894599095, 0.03813129905647845, 0.03972801471509725, 0.04293960647277185, 0.04407495256829299, 0.040561528857132, 0.043393394160625674, 0.041996152594041554, 0.04685644378822259, 0.04856280580805144, 0.046809178126973265, 0.04424396363307654, 0.041648334111766744, 0.04457339576947193, 0.0438704706372499, 0.04722660822955548, 0.048975317721380494, 0.05099751955984494, 0.045196632152572254, 0.04375949513010004, 0.04823866984913788, 0.04283943546528775, 0.04450218575354305, 0.04959442237813569, 0.044598450570579345, 0.045383444545267694, 0.04889086313647761, 0.04345621542802586, 0.03854385804903462, 0.04374060954710464, 0.046377868479817466, 0.041197667049519125, 0.04095096726615447, 0.044333026742799854, 0.04623895604407687, 0.04826185971132163, 0.0441627163404013, 0.04063182071999084, 0.041338900071159104, 0.040203181836910226, 0.04492267621141167, 0.045984188296184385, 0.04040230185341752, 0.03978185849843574, 0.04079943094380277, 0.047007096003839545, 0.04812590599970763, 0.04717030754481706, 0.0459307097216564, 0.04066339741489514, 0.041312227425623446, 0.040711529180376425, 0.04268242277140402, 0.04357952610678519, 0.04460514148809379, 0.04360334642450869, 0.04217011514508262, 0.04334356732292601, 0.0426601554327807, 0.04532069637427389, 0.04823676862487624, 0.045667511984581695, 0.044231867004805796, 0.04067698213716437, 0.042301026031533816, 0.046032650730078124, 0.0502632040836994, 0.045268956715898714, 0.04422710022390402, 0.050048924885368135, 0.04936238741438512, 0.04884961098281242, 0.04541478364228696, 0.042051393694077135, 0.04285437875092125, 0.045111370413742455, 0.04316333020062166, 0.04409392902623922, 0.047355949425623246, 0.045271543984281286, 0.04452050015459909, 0.04423801381696368, 0.04146997568706274, 0.044123659698225896, 0.042347738743070705, 0.044512751199966796, 0.04493190395041527, 0.042063564085400074, 0.04958463923320106, 0.049813459387438405, 0.0447534592641763, 0.04581001339076328, 0.045164188326244664, 0.04534191229420978, 0.0433577860044516, 0.046267592288139314, 0.05048425119409419, 0.04697048251643655, 0.04576237334347273, 0.045984144391259675, 0.04465997077028123, 0.04660665694779917, 0.04997508324166489, 0.04804601469145342, 0.04233935040503736, 0.046429209443933056, 0.0445208876654168, 0.04208585164139094, 0.044578260742284166, 0.045000201891995925, 0.04792656158563008, 0.047775413462989726, 0.04208273790613025, 0.041706847522972545, 0.04393964273681822, 0.04530998058979156, 0.04550218025463763, 0.03832524689358054, 0.04100502849706543, 0.04641374035734766, 0.04584973453276269, 0.04580159413050548, 0.046784541550137705, 0.04385214018801729, 0.04029276159689517, 0.040686039077996033, 0.04873627362006497, 0.05130127206345877, 0.045431482631544916, 0.04133250893654304, 0.0416386491868112, 0.04294266273331856, 0.04196162033028605, 0.0406649633039519, 0.04384021080545183, 0.047661755902046916, 0.04600532631585428, 0.04300169756660305, 0.04085291913602358, 0.04910985147085584, 0.04929081694290266, 0.04729309365118642, 0.04517590609210008, 0.04202213915403456, 0.04668505375514549, 0.04786878245167835, 0.04107916567933025, 0.03850802748365839, 0.04081987633979023, 0.042737366425823375, 0.043413895415187864, 0.04218703348381408, 0.03945270295711564, 0.04391428585894295, 0.047884010334107316, 0.04757222467291381, 0.04503266366115083, 0.04676719857057471, 0.04640028386165369, 0.04131335118237111, 0.0379253864981139, 0.04235123672668414, 0.043777418434369164, 0.04464069748203905, 0.04602984443274384, 0.0444235419039592, 0.0450731542177212, 0.04431506616919288, 0.04272779749146518, 0.04643692167421165, 0.04334976640637093, 0.04394655541351319, 0.047525258197149894, 0.04533764381144537, 0.04301682840740343, 0.044820413491902196, 0.04887846499311324, 0.04738121939684068, 0.047036076038968866, 0.0462849801038192, 0.041589895501756007, 0.04228381380296282, 0.040439798848707256, 0.04305937788157469, 0.043036144192028045, 0.04359666524754318, 0.047310139494077794, 0.047141173978784714, 0.04232036166029365, 0.04018438616230611, 0.04091041466820697, 0.04072898120800978, 0.041001015435087475, 0.038665856486825634, 0.037520577564654635, 0.0429466986273869, 0.04491811616903804, 0.04377956449186469, 0.04430656624472989, 0.044093035435001766, 0.04554197088329801, 0.045109936393839145, 0.04568168763293486, 0.04503904702996461, 0.048361819488318664, 0.04856231299386956, 0.044173851807740076, 0.043969312716185184, 0.04535774012258389, 0.04843241189336424, 0.04556991912059686, 0.04393936551427311, 0.042317539877358666, 0.04310168135235653, 0.04492402972541233, 0.047781370875450535, 0.04634689797069972, 0.04277704277954059, 0.042155675898404715, 0.043618610312103095, 0.042994782109459215, 0.046109042761003624, 0.04421444544575973, 0.0443468062917926, 0.04703752185381213, 0.04614662815255516, 0.04669195603661659, 0.0492082392726853, 0.0481163703837761, 0.050114354281338286, 0.049509732277893806, 0.043243896196099786, 0.041194536878337254, 0.04253857886594465, 0.0465506824662963, 0.04846761677643726, 0.044805828619211154, 0.04325724380684927, 0.04750366389782178, 0.0455562513451488, 0.044414208903052955, 0.04113213846794192, 0.0452527281825917, 0.04910609472034695, 0.045546230887711384, 0.04453497847058054, 0.04753522680779338, 0.04247062535062559, 0.04024724928038783, 0.041260329742860334, 0.044950808714978696, 0.046379001145607134, 0.04560023251082345, 0.04314569614206234, 0.04559141853523529, 0.04284225993304479, 0.03925095201871566, 0.04109047420288943, 0.0419338531923751, 0.041420046002953924, 0.04353655525388458, 0.046007448171279865, 0.04456476993235748, 0.045373239181726706, 0.04777332114716908, 0.04626115342707098, 0.0447112807041298, 0.04220418749341254, 0.04470197717838656, 0.04577310954667083, 0.047693222072940036, 0.04994428513762421, 0.0510868257750944, 0.049513846122354994, 0.0444376271074247, 0.03895732343521757, 0.04001880584455591, 0.04215951055590354, 0.045868391946666816, 0.04862617192679427, 0.04419161646429283, 0.04106155962097399, 0.04147744096233234, 0.04712807458262805, 0.04877882022292148, 0.042825287679554944, 0.04563594789027461, 0.04211661673706236, 0.04251196968188058, 0.04405432710511507, 0.043812118750595444, 0.04697643064765009, 0.04897610661944044, 0.04444426226071319, 0.04083415775407893, 0.042801833061298776, 0.04648722245792107, 0.04483431277413052, 0.040179658539282376, 0.0418300044882904, 0.04246092750101033, 0.04317485130898904, 0.04685253181731342, 0.04789146269762151, 0.04277139233419926, 0.043603695738828, 0.03942086334021424, 0.039977581662169005, 0.04295899016130693, 0.04443937022995942, 0.04236624630097043, 0.04028808981296091, 0.0416731889841848, 0.04594942894110228, 0.048431672208516255, 0.04487167533068471, 0.04473761380737271, 0.044089297350886285, 0.04056077569652711, 0.04496300034741248, 0.04458917833391016, 0.046203896841951336, 0.04902372388096177, 0.045789832956717604, 0.04339741875046452, 0.04288864632051049, 0.041506872482375355, 0.04539014131826737, 0.04392604424499593, 0.04091754425001383, 0.04873362918434114, 0.04491556290440024, 0.042799345317106705, 0.04400514957656507, 0.044786587108495236, 0.0458768089804505, 0.04654466119290923, 0.04573782763147599, 0.04479591158613801, 0.04187599907014865, 0.04463356969404915, 0.07383305384898482, 0.15751382430218738, 0.27472115638004646, 0.368783044191283, 0.3912284847815476, 0.35360357348210525, 0.3119389046642328, 0.2830165653068165, 0.2650162506083908, 0.2531648859079913, 0.2593019773964714, 0.27107519659604734, 0.27968795434300076, 0.30047562343367507, 0.2823867490849573, 0.21481516341963103, 0.133144150238149, 0.0688307938988062, 0.04473130930565228, 0.04385127435822522, 0.044553736421294565, 0.04127664968207463, 0.04454622987207048, 0.04884031127046417, 0.048122871573130485, 0.042082200348678094, 0.042039352432849765, 0.04436835561225282, 0.045091165600214096, 0.04636928307044969, 0.05039329910350356, 0.046666733355361215, 0.045976521769172364, 0.04644611424242668, 0.042497365614117745, 0.04085586699527161, 0.044718871140089406, 0.049191848527691415, 0.04831230615580526, 0.04604736950731744, 0.0444485065966959, 0.04296083827915726, 0.04834935939501341, 0.05344790413704256, 0.05264127537186831, 0.04917340154250918, 0.04295018624969584, 0.04189545817965424, 0.04593404060246553, 0.04943685140881612, 0.051237239860025714, 0.04800403498095999, 0.04563253866917184, 0.04524124207474611, 0.04458011630345187, 0.03991521062851137, 0.04435924813543686, 0.04610314800106169, 0.04521929955241773, 0.04516156699302046, 0.04336122436996295, 0.04342669135028366, 0.0423213128071148, 0.046103989438112114, 0.04548510824793136, 0.04573135088272092, 0.045558630524557386, 0.0436007703230872, 0.04399841505484506, 0.04395345016107138, 0.04249579019867823, 0.03958433567052239, 0.040140458456784926, 0.04490212011778133, 0.04466717634753759, 0.043296476918217794, 0.0438984591376328, 0.04554898299571972, 0.04185759632961884, 0.04311776740719089, 0.04231304973898216, 0.04462883216372842, 0.04544288016945238, 0.04799980861802739, 0.04352392205837489, 0.04058689724493553, 0.043855783107433675, 0.047089921264400994, 0.050339723780790116, 0.04715675872483396, 0.050799211436584515, 0.04969153772503435, 0.040255819132993885, 0.046138190353369087, 0.051121964219688634, 0.0481747261575124, 0.04747408206498616, 0.047631044595418816, 0.04571777860139955, 0.046718688281015525, 0.044688513071596395, 0.04139068773816812, 0.040585791964882094, 0.04248485300973531, 0.044442388589751795, 0.04599719804788843, 0.044881004867125054, 0.04203409121686328, 0.04247716343413772, 0.042716448129993735, 0.0487247629943006, 0.04283717177511844, 0.04385926210360578, 0.04915414938870767, 0.043670920673914114, 0.03897746419631116, 0.0416799890357524, 0.03986292922438967, 0.0409603023864511, 0.044296318390389375, 0.047826194910518996, 0.045084927372510485, 0.044773663643888295, 0.04459979856154452, 0.04324230427029563, 0.04312015855238737, 0.042447979854178974, 0.04536531831457631, 0.04457949965728459, 0.040977602515928205, 0.041266003136808775, 0.042298016596590056, 0.043935383206010585, 0.04514509818460419, 0.04634437233613818, 0.046146532533406434, 0.04699243961223198, 0.04532946642124093, 0.042776208747042556, 0.04156547220215017, 0.03778673194603044, 0.04331178066018882, 0.05005967063107201, 0.046732290718825564, 0.04599330837289554, 0.04409241571393664, 0.04604737900642473, 0.0450948906101456, 0.040620631434856, 0.039849014521260474, 0.04265038675166294, 0.041418824248551735, 0.03986374956982926, 0.04243198746181498, 0.03860709515664676, 0.037923485276092336, 0.04009129923018547, 0.042289994033405505, 0.04511397983866875, 0.046050672985150226, 0.04704286838612696, 0.04755649394255143, 0.051720559853125965, 0.04641884351606359, 0.04105030742776407, 0.04175176225163367, 0.04519128392111783, 0.04437434984697224, 0.04952423834310004, 0.048496370640049406, 0.04687715662471721, 0.04578936141706773, 0.04563481201293892, 0.04598734688571377, 0.048303042155078836, 0.04855980247639367, 0.046344091867390864, 0.04224802294392569, 0.04280104610765892, 0.04311968991869175, 0.03973833726243323, 0.04252639199315207, 0.04396495355421774, 0.04867315541821018, 0.04771258346442159, 0.04145585174086321, 0.04124460746100822, 0.04753754721706368, 0.05078308687374586, 0.04774925528269905, 0.0442549672299769, 0.041549414324750386, 0.03819716495335711, 0.04182845046500639, 0.04741954458919868, 0.04633725652409185, 0.04821793525615074, 0.04593968716153478, 0.040952521979733494, 0.03812852920639058, 0.03994385264217905, 0.04542514811247429, 0.04396578671188969, 0.045120592979894694, 0.043434386671745685, 0.04486421236374555, 0.04511074565844008, 0.04802499233485515, 0.04724204220949786, 0.044951730334179665, 0.04112152694771015, 0.043078488066104724, 0.0385162410576129, 0.043293086976470004, 0.0471132771108006, 0.04410307782726108, 0.04387401019599393, 0.04180951017279642, 0.04004578767269875, 0.040480313607543346, 0.04076801407252453, 0.04352460306385267, 0.04178515272241207, 0.04217311932068119, 0.03969030678515222, 0.04132126848635888, 0.042157798921646436, 0.044255642418621906, 0.04770981347817974, 0.04311719627768798, 0.03982909176990479, 0.04398239494038796, 0.045403657731088995, 0.043541927113499225, 0.04718488479279637, 0.04304443128437063, 0.04512819219304071, 0.04250099490188445, 0.04140277290677278, 0.04210850202424541, 0.04590567675129869, 0.04998384034128927, 0.048631225898609146, 0.04282383230527812, 0.04068610757469241, 0.04175022986564506, 0.04561173491918591, 0.043037181542401326, 0.042735729420677646, 0.04711015550138621, 0.04876725773901864, 0.045947931047572055, 0.04292538928315598, 0.042094957050042514, 0.04506976761344216, 0.0420390276373629, 0.04384401856736436, 0.0407438404286659, 0.037309382244192654, 0.04017186285174903, 0.04996072374539154, 0.051546742792534433, 0.042900168545030495, 0.04231307447086904, 0.04841202017534444, 0.046803479085083054, 0.0422392870468687, 0.041837921239609804, 0.04387364302666453, 0.04771511882584397, 0.04489245415874867, 0.045566105556081246, 0.0479484693358972, 0.0466099433972263, 0.044332331907810434, 0.04023773080761353, 0.04149386785933833, 0.04169059303916562, 0.04320477272567272, 0.04649729163008696, 0.044685998781404786, 0.04334826052416294, 0.041075292453660865, 0.0409086435943719, 0.04438677754807575, 0.04623259599303403, 0.04262395865811467, 0.04421842895646172, 0.046397035979818126, 0.04815302890960915, 0.04718995215078028, 0.04671884022048119, 0.04880515874114139, 0.049024376791472216, 0.047156568441048236, 0.04496106537168301, 0.04360819573711831, 0.04151226151966121, 0.04291079956707027, 0.04129661590010412, 0.04542311368537912, 0.04831368525423336, 0.04534866465794432, 0.04342639028240374, 0.04739395365963675, 0.0443999663306812, 0.04954808861642869, 0.05264055219523586, 0.04679345298278846, 0.044364933642740155, 0.04265013981604563, 0.04550962607805169, 0.0479896422051963, 0.04550963923161413, 0.04724227547655426, 0.04986736947859932, 0.04803139668261429, 0.044552796614394444, 0.0443755754252993, 0.04063766096222834, 0.04144587333223653, 0.0398394486203249, 0.042480381984739295, 0.047019434032308675, 0.04546161733596401, 0.04450221601501025, 0.04210006715334784, 0.04005226407244173, 0.03955949729867939, 0.04071413038101421, 0.04535234698255664, 0.047714640377782024, 0.04723143924564565, 0.0489996468945443, 0.0507639034504261, 0.04463070405136204, 0.04343480122322345, 0.04402213712737706, 0.04222225935054869, 0.04086928432127433, 0.04457500088051577, 0.045325125972177915, 0.04021537227984362, 0.03903167796065588, 0.04113412321386732, 0.04241378029802823, 0.0436563105005204, 0.04408515905176545, 0.04426686939765035, 0.044623459224155934, 0.04061947506992713, 0.04294247726242183, 0.04590406146179182, 0.04381319115688296, 0.04238119157479209, 0.0494737738204695, 0.046978391853747796, 0.044047842018357736, 0.04253331101940259, 0.04071700511056543, 0.041635819303615755, 0.045563108402065934, 0.04569843807581536, 0.04312239571980319, 0.04273464950121801, 0.040612618695548684, 0.039769212246394035, 0.03932998633452101, 0.03885884200096123, 0.04190581341238799, 0.04273499713849331, 0.04164896183985993, 0.04070007031063211, 0.0420272316357413, 0.040230273993299724, 0.04197888149995172, 0.04272821747535997, 0.045700136612551075, 0.048276502336969146, 0.04881748037066385, 0.05010541036583158, 0.046733237833424014, 0.04117830790285544, 0.04072704524999467, 0.04357418494661676, 0.04176725979945469, 0.04231750930504891, 0.04604853491665204, 0.048257271608920696, 0.049035286283240595, 0.04945044080814003, 0.04551648487109603, 0.04243225442695252, 0.04039790435173626, 0.04119839463206468, 0.044516603090115694, 0.04573853728607062, 0.04346624722433721, 0.042715819226247724, 0.043341412414272466, 0.040886030133109694, 0.037861121369597656, 0.038449945701281824, 0.040609789948322074, 0.044563882834827734, 0.04638592010123272, 0.049017386518984865, 0.04777675130642749, 0.042238944842525175, 0.042340543158101526, 0.04249511716555553, 0.039916689768369894, 0.04145340605624698, 0.04351090665114707, 0.0453638255368915, 0.04343789144783731, 0.04606007450910063, 0.04152050332847474, 0.039293060177288716, 0.03722550195491536, 0.03824560577383229, 0.043921156697001966, 0.04261266738229992, 0.04414562100055506, 0.044960650567102405, 0.04057196995268912, 0.04481940872175942, 0.04945682385555079, 0.0429563280561257, 0.04120694478319821, 0.04332187458329594, 0.04083072941761571, 0.039886773825731435, 0.042963070735647146, 0.043932456301492416, 0.040544992287834286, 0.045073176950790605, 0.04631541767733007, 0.04439695360787865, 0.04807195554531096, 0.04584552845209921, 0.04004618109919624, 0.039551773849527495, 0.03794331000939169, 0.0435501070977116, 0.05039333427607532, 0.050699015971540806, 0.046966925996745464, 0.04380061638722383, 0.037397083511204125, 0.04032121419597859, 0.04153244007680435, 0.041968878264015005, 0.043188078701954576, 0.041566792909484236, 0.04173550070837873, 0.04686472330266697, 0.049266541534457545, 0.05048503132281127, 0.048218814614690034, 0.04698858208160473, 0.04441195152096434, 0.044496322214965835, 0.04647433807270051, 0.04555297755048278, 0.04159091108652965, 0.0440746536149738, 0.04594032549158145, 0.04670074914681843, 0.04683434360347864, 0.047537842473615734, 0.04905478939125279, 0.04315813455821698, 0.041128316651189975, 0.03724682040194883, 0.03879446979130044, 0.04331714600968578, 0.044345941310747816, 0.04213540362224531, 0.0425772211095165, 0.04323900487303508, 0.048100405392711945, 0.049985649271879234, 0.05052176954811673, 0.04367419759415172, 0.04096814524905894, 0.04563879798658588, 0.04901742905437848, 0.042731503179514886, 0.045659182274236, 0.04810819737106573, 0.04984569553629141, 0.04495550807203978, 0.04682187356898137, 0.04680512719708847, 0.046997339447880024, 0.046837688209724636, 0.048303264827979614, 0.05305189543377626, 0.04892989883788012, 0.039621168341751285, 0.043704957441342986, 0.04271354049908887, 0.04602263214439308, 0.047798471604359884, 0.045067635939786695, 0.047011509534419856, 0.04606975188074629, 0.04357823887179839, 0.042849397283256103, 0.042539196309035666, 0.04537907193640049, 0.04540241491597833, 0.0374899227030091, 0.044552780114536276]],
            'im_sum': [False, 2529.879352938107],
            'regn_sum': [False, 270.0249538719654],
            'npts_real': [True, 6400000],
            'profile': [False, 2.8040067211040336],
            'fit': [False, [3.1292950788543084, 6.121292471264506, 5.778675891349287]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.5049652049504],
            'fit_pix': [False, [45.635362794917455, 41.05818466058917]]}



        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_standard_cube_eph_pcwdT.exp_mask_stats
        exp_mask_stats = {'npts': [True, 6400000],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'mask_pix': [False, 9553],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_pb_stats
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

        # test_standard_cube_eph_pcwdT.exp_psf_stats
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
            'min_val': [False, -0.16455982625484467],
            'min_val_pos': [True, [34, 36, 0, 949]],
            'im_rms': [False, 0.08653504136921199],
            'im_sum': [False, 2857.8352289416644],
            'npts_real': [True, 6400000],
            'fit_0': [False, [1.095566785329151, 4.020892595865358, 2.9706317101525186]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.4453942426872],
            'fit_pix_0': [False, [40.00007105295384, 39.99750303438477]],
            'fit_1': [False, [1.0955778815758896, 4.020331875996908, 2.9703856624539005]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.5063079930342],
            'fit_pix_1': [False, [40.00007131349202, 39.99750438784483]],
            'fit_2': [False, [1.0956198671703976, 4.019978907059439, 2.969213090005495]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672217433812],
            'fit_pix_2': [False, [40.00007226788366, 39.99750678565682]]}


        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_resid_stats
        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 0.32578590512275696],
            'max_val_pos': [True, [20, 53, 0, 497]],
            'min_val': [False, -0.3510037958621979],
            'min_val_pos': [True, [32, 47, 0, 490]],
            'im_rms': [False, 0.04672249030226864],
            'im_sum': [False, 331.0043877220577],
            'regn_sum': [False, 51.16192187555134],
            'npts_real': [True, 6400000]}


        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_eph_pcwdT.exp_model_stats
        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.4211243391036987],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.09897369146347046],
            'min_val_pos': [True, [58, 39, 0, 493]],
            'im_rms': [False, 0.001760109582469154],
            'im_sum': [False, 55.12738844938576],
            'regn_sum': [False, 55.89786731824279],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400000]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_cube_eph_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1011.9810180664062],
            'max_val_pos': [True, [0, 0, 0, 857]],
            'min_val': [False, 1011.7631225585938],
            'min_val_pos': [True, [0, 0, 0, 475]],
            'im_rms': [False, 1011.8444616762994],
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

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[0.0, 3.25])
        self.mom8_creator(img+'.residual', range_list=[0.0, 3.25])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')

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

# End of test_standard_cube_eph_pcwdT
#-------------------------------------------------#
    # Test 4c
    @stats_dict(test_dict)
    def test_standard_cube_eph_briggsbwtaper(self):
        ''' Single field multi-EB ephemeris cube imaging with briggsbwtaper - field 21PGiacobini-Zinner, spw 20 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
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
            gridder='standard', mosweight=False, \
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
            gridder='standard', mosweight=False, \
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

        # test_standard_cube_eph_griggsbwtaper.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 4.555308268539054],
            'com_bmin': [False, 3.365998340500264],
            'com_pa': [False, 84.3021528064608],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122071.6439819336],
            'start': [True, 354445300000.0],
            'end': [True, 354567200000.0],
            'start_delta': [False, 354445300000.0],
            'end_delta': [False, 354567200000.0],
            'nchan': [True, 1000],
            'max_val': [False, 3.3076226711273193],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.37107840180397034],
            'min_val_pos': [True, [32, 47, 0, 491]],
            'im_rms': [False, 0.056953113383683536],
            'rms_per_chan': [False, [0.051167301929209416, 0.04910436885589881, 0.047426032870225364, 0.045946628051608096, 0.04150331116183864, 0.04234693717979622, 0.044663430812569856, 0.043975939214317375, 0.044954768263404686, 0.04687210851723484, 0.04510794530454117, 0.045144645019892946, 0.0434134335857173, 0.04305283981314241, 0.04211918242030221, 0.04635357418539099, 0.044074957871829246, 0.04311158306164848, 0.0408540834566592, 0.04423331468629409, 0.04787649581968905, 0.04677294689350967, 0.04848520754401066, 0.04516929134299152, 0.04536894846353064, 0.044035800523968396, 0.04098901435500981, 0.04844887693236758, 0.043962842880926674, 0.04387303288953025, 0.04780452139728523, 0.046661958475228006, 0.0477697101819565, 0.046940860796687064, 0.04743399138507472, 0.04825726825255441, 0.04789540670908134, 0.04878631142893934, 0.04672167116570699, 0.04141340836398444, 0.04029918536841034, 0.04260048332283974, 0.045475921739909934, 0.0429176695825316, 0.04085400717296316, 0.041307406256973336, 0.040466319575463046, 0.0426085735663968, 0.04537123976891103, 0.04631957526377112, 0.04299460623427574, 0.04257424125694552, 0.04546736789215301, 0.05022678418418266, 0.04629758484695974, 0.04328924551119241, 0.04765008616691762, 0.04854532697994634, 0.045066336735905765, 0.04191367592048695, 0.04192701088077219, 0.0407387960314003, 0.046945789363474, 0.04988122591851493, 0.047667710953245124, 0.04709489520385384, 0.04740250932056808, 0.0449781491298397, 0.045822647172837086, 0.04543818744363286, 0.04268119754407437, 0.04572685844427544, 0.050673702265752574, 0.05136403630064312, 0.0494839872897862, 0.045330709491184444, 0.04585223345023971, 0.04620617921870042, 0.04468015259889807, 0.04825989687208062, 0.04788559152872278, 0.044978708346781, 0.04249494868487235, 0.04134911962046679, 0.04463444405566262, 0.04557250758980735, 0.045499310801702025, 0.04434777790568684, 0.04654028760660557, 0.04589049853243483, 0.0433363038018721, 0.04242912925827068, 0.04343235419189811, 0.041761834316465796, 0.042730589319263095, 0.04559773609819441, 0.04588525418434693, 0.051779926474691554, 0.05437147848196889, 0.045644787614140224, 0.042209476545486055, 0.04169334542072149, 0.045027644188701636, 0.04857259833436285, 0.04719430134362024, 0.04657152785039264, 0.04834640355158802, 0.0475728190845166, 0.04504857260942262, 0.0491726190582117, 0.04744733272307072, 0.04012916808671003, 0.0367953655153898, 0.04037278805805138, 0.049038930898401945, 0.04670610471281036, 0.04628753577486229, 0.0445697407364478, 0.04890409695612506, 0.047632272141199535, 0.045886098435698366, 0.04640537262021118, 0.044955482856975856, 0.04353760796440149, 0.043622731879179394, 0.04506446170113727, 0.047695162056200405, 0.0463583645840669, 0.04385891796783386, 0.045552416177023644, 0.044600898903087896, 0.04458178214587721, 0.04857922867613965, 0.051301407484537265, 0.048439653475603105, 0.04235267893020762, 0.0461515349984245, 0.04358174893943753, 0.04049563415153134, 0.04047251488121909, 0.04246601301169996, 0.04515141323177437, 0.044257975850365575, 0.04312022110492146, 0.04073565293862463, 0.03812972922927284, 0.03972641637397884, 0.042937653810839894, 0.044072517240549156, 0.04056004606449945, 0.04339114395572016, 0.04199345120787946, 0.046854068843551124, 0.04856033781955519, 0.04680725030563083, 0.044242286761731746, 0.04164677260005263, 0.04457305197892608, 0.04387041573045802, 0.04722490304590815, 0.048974080215254207, 0.05099662599062068, 0.045195284045144175, 0.04375784822438968, 0.04823754087986354, 0.04283899495148726, 0.044500276092780476, 0.04959292156396535, 0.044597222630010196, 0.04538169334911277, 0.04888858245058824, 0.043453615305275924, 0.03854136780722329, 0.04373778547001929, 0.04637572612557875, 0.041195799934323, 0.040948683695996344, 0.04433086238295178, 0.04623740630074435, 0.04826054226953814, 0.04416156316861346, 0.040630413837077585, 0.04133739303791001, 0.04020163911182845, 0.04492182537116222, 0.0459837626521909, 0.040402241183903384, 0.0397815636700942, 0.04079827130110012, 0.04700630350258279, 0.04812567095249552, 0.047168449414278966, 0.04592877243654917, 0.04066233989118301, 0.041310551604324876, 0.0407108642911298, 0.04268145609913962, 0.043577196277068776, 0.04460330589891903, 0.043602214087116974, 0.04217003617777811, 0.04334290704724179, 0.04265880019494553, 0.04531858174058441, 0.048233501829010206, 0.0456650860083565, 0.04423112185311475, 0.04067563146453821, 0.0422991073375842, 0.04603082211930534, 0.0502620094872356, 0.04526694235243964, 0.044225058731084656, 0.050047395612416314, 0.04936085135871181, 0.04884823788380334, 0.04541376471664849, 0.042050457882486195, 0.04285394065441368, 0.045110310645639376, 0.04316208097301357, 0.04409213188630146, 0.04735472591579041, 0.04527046707370761, 0.044519148528667, 0.04423654779127111, 0.041469065705295925, 0.04412285846920254, 0.04234706653246038, 0.04451157206127082, 0.04492930398065287, 0.04206095741692099, 0.04958340621914414, 0.0498125918624894, 0.04475192531685877, 0.0458082292578735, 0.04516196328102129, 0.04533862153185725, 0.04335497745635367, 0.04626568431847512, 0.05048203605154181, 0.04696784784160007, 0.045760147324596705, 0.04598224882780486, 0.04465940378839121, 0.0466056083488391, 0.04997359122425284, 0.04804564406120846, 0.0423395980895589, 0.046428721977784526, 0.044519148761021185, 0.042083923178046115, 0.044576206491888315, 0.044997894498215035, 0.047925676028339144, 0.047775183842988025, 0.042081907236965733, 0.04170554707682394, 0.0439375401391505, 0.04530758722177607, 0.045500514545616884, 0.03832423963373504, 0.04100305029255071, 0.0464108731092276, 0.04584781852553112, 0.04580075552392769, 0.046784024034035636, 0.04385144124950881, 0.040291637289974355, 0.04068479186568497, 0.04873437479021023, 0.0512990647688585, 0.04543030266967547, 0.04133204268948601, 0.04163725373577503, 0.04294082289552146, 0.04195960948919333, 0.040662850405114796, 0.04383753741472004, 0.04765904713502353, 0.046002717014045465, 0.042999085865481866, 0.040851193677844316, 0.04910814524399428, 0.04928897629086158, 0.04729198593194097, 0.045175214948102636, 0.04202232492272581, 0.04668553085860465, 0.04786883737180495, 0.04107864276692618, 0.03850689954435268, 0.0408180477745673, 0.04273490951725388, 0.043412046301446246, 0.04218561026582016, 0.039450384489229495, 0.04391295305651671, 0.04788200402974977, 0.047570122198778965, 0.045030752178062086, 0.04676598933305676, 0.04639921662502428, 0.0413126803082296, 0.03792353301984149, 0.04234946152482861, 0.04377659037687563, 0.04463905809107571, 0.046028394986651115, 0.04442150008306717, 0.04507016964246473, 0.044313368869529135, 0.04272617956673372, 0.04643436264956652, 0.0433483593238604, 0.0439459127868522, 0.04752564014226432, 0.04533714229095887, 0.04301552627613123, 0.04481942837028821, 0.04887696695891762, 0.04738038171012137, 0.04703491709237416, 0.046283442598835454, 0.04158802224800722, 0.042281258005511996, 0.040437156805376614, 0.04305831351448433, 0.0430359736773737, 0.043595355789154744, 0.04730794463740114, 0.047138681394011735, 0.04231909826098159, 0.04018340196891829, 0.04090947507944172, 0.040728197021056324, 0.0409998747472632, 0.038665020789864404, 0.03751947473990552, 0.04294535063590318, 0.044916627427408914, 0.043777779092230816, 0.04430459857468188, 0.044091128385108945, 0.045540223780373756, 0.045108372977417814, 0.0456802901859014, 0.04503692694773148, 0.0483601875123851, 0.048560886497349944, 0.044172332056563936, 0.04396835619286222, 0.045357209400154114, 0.04843110825482586, 0.04556834288345021, 0.0439377071301504, 0.04231507304030345, 0.04309999165533348, 0.044921854827832286, 0.04777768426717635, 0.04634493370777653, 0.04277546264150503, 0.042154339932426316, 0.043617158706355255, 0.04299355718010844, 0.04610841782349669, 0.04421316503306352, 0.044345074763503356, 0.04703535473860608, 0.04614490874454537, 0.04668973561926332, 0.049206004494300536, 0.04811433318636512, 0.05011208602086733, 0.04950746973015089, 0.04324283897073358, 0.04119341247395139, 0.04253682798287216, 0.04654792845755986, 0.04846386476585031, 0.04480342757601164, 0.043256207038701595, 0.04750182647445544, 0.045553518548771224, 0.044411711154185546, 0.04113084346680847, 0.04525186481492797, 0.04910354221721132, 0.04554475988014388, 0.044535017657930935, 0.047534605903032136, 0.04246905191070787, 0.040245272521727986, 0.04125900941914659, 0.044948539251804115, 0.046376536600511276, 0.04559854547193087, 0.0431446867633259, 0.04558912321932169, 0.04284075160272995, 0.039250748805109766, 0.04109010492909606, 0.04193348635332057, 0.04141936794040166, 0.043534707556986074, 0.046005594809248866, 0.04456385259607724, 0.04537272398838617, 0.04777073758693133, 0.04625883056101605, 0.04471021935014107, 0.042203331426204145, 0.04469996382735254, 0.04577058107415554, 0.04769077407550813, 0.04994193402226542, 0.05108495036290946, 0.049512946311032424, 0.04443695411291232, 0.03895613680452749, 0.040017483183998975, 0.04215722435861593, 0.045865895605531205, 0.04862345863563123, 0.044190105685044216, 0.04106039270347619, 0.0414757829997012, 0.04712598239736813, 0.04877746817365123, 0.04282456338885783, 0.045634940524389114, 0.04211511498305421, 0.0425098479427093, 0.04405180003312188, 0.04381038626272056, 0.046975412244943976, 0.048975030371218704, 0.04444320269641191, 0.04083252379251179, 0.04279983894288616, 0.04648628966994779, 0.04483416022622179, 0.040179073884111946, 0.04182822845082981, 0.04245922561578296, 0.04317254609259317, 0.046851132593229726, 0.04789146799301545, 0.04277056970988489, 0.043601795578819616, 0.03941929009322431, 0.03997619008848362, 0.042958415909501585, 0.044439106466319094, 0.04236428378265602, 0.04028584925330022, 0.041672731843963286, 0.04594884521390053, 0.04843000678740922, 0.04487089877775389, 0.04473670753821486, 0.0440873433437084, 0.04055832856840577, 0.04496184555127177, 0.04458752976189551, 0.04620181440634824, 0.04902146190715618, 0.04578649157903686, 0.043394252037612636, 0.04288568627339057, 0.041504510591013205, 0.0453884448537182, 0.043924780722496705, 0.040916075844982615, 0.0487318144606915, 0.04491361754618209, 0.042799373417346964, 0.04400459570738337, 0.0447841958806134, 0.045873648621299876, 0.04654138553580639, 0.045734498902470076, 0.04479355388398826, 0.04187423488021391, 0.04463081332358924, 0.07382726420640791, 0.15750018670556137, 0.27469763970001776, 0.36875164866294413, 0.39111510028261853, 0.3532378527482003, 0.31180255965861425, 0.2829911733694037, 0.26499202307947617, 0.25314250276309064, 0.25927964804019227, 0.27105178016381376, 0.27966324998915093, 0.300446661451936, 0.28236162365968176, 0.21473008805460964, 0.13316231085407948, 0.06882529881515627, 0.04472945530847343, 0.043849207449791684, 0.04455102653491761, 0.0412747833449623, 0.04454495295272398, 0.04883827015793973, 0.0481207638516202, 0.04208085602856363, 0.04203691120974326, 0.04436606808988699, 0.04508968540731526, 0.046368784414902116, 0.050391761638093634, 0.04666551678301536, 0.04597546653682626, 0.04644448425549627, 0.0424961808520706, 0.04085466533306353, 0.044718021727328916, 0.0491894412608543, 0.0483102666251601, 0.046047883997327295, 0.04444716289490372, 0.04295805578453906, 0.04834681878467435, 0.05344481946638864, 0.052639213000270295, 0.049172302488352186, 0.042948590256469694, 0.04189301658108049, 0.045931527837792425, 0.049435184524373454, 0.051235930728614834, 0.048003964531681764, 0.04563284373742887, 0.045240128333993046, 0.04457818143480976, 0.03991292200106255, 0.044357593680924436, 0.04610233912101146, 0.04521840250579686, 0.0451604229069975, 0.043359496113145865, 0.04342516765304109, 0.042320945108052455, 0.04610464003322387, 0.045484702672244345, 0.04573066738291446, 0.04555745813319026, 0.043598864214186184, 0.04399613749357084, 0.0439522974935923, 0.04249466053097157, 0.039582871810896425, 0.04013969689870331, 0.04490130182157931, 0.04466558999488214, 0.043294361034486556, 0.04389709612437497, 0.04554753046215814, 0.04185663951270954, 0.04311695453627529, 0.04231207209229626, 0.04462770606708744, 0.045441501108417465, 0.04799829383496495, 0.04352218885678298, 0.04058494143324349, 0.04385412133285266, 0.047089181489416015, 0.05033822470115026, 0.04715436905690689, 0.05079709779105723, 0.04968907830525051, 0.040253823081395336, 0.04613634979788893, 0.051120464671912454, 0.04817401985191957, 0.047472633344022436, 0.04762866725280526, 0.04571593743504577, 0.0467171106484833, 0.0446870697725789, 0.04138929531255531, 0.04058412003133212, 0.04248304665297155, 0.04443983645090936, 0.04599481878794733, 0.04487895123125822, 0.042032011935095315, 0.042474745056319124, 0.04271522096979415, 0.04872322857321779, 0.04283549760149551, 0.04385796133203299, 0.049153272352852155, 0.04366933855894109, 0.03897548549238558, 0.04167898947444148, 0.0398617745202848, 0.04095868993221706, 0.044294652076103086, 0.04782379512262603, 0.045082894604297645, 0.04477207047875748, 0.044598371185995495, 0.04324115497153855, 0.04311894533623701, 0.042446771578127075, 0.045364302328767606, 0.04457843730735818, 0.04097751201319639, 0.04126604143893337, 0.042298147326806014, 0.04393512476763845, 0.04514396345649417, 0.046343925561560854, 0.04614592036802372, 0.04699069625101227, 0.045327711849206025, 0.04277538574603792, 0.04156495226165961, 0.037785075009890805, 0.04331097025660082, 0.05005921891901261, 0.046731739076414196, 0.04599332496043983, 0.04409199436569674, 0.046045077813329414, 0.045092927473862975, 0.04061986970715164, 0.039848069932693395, 0.042649158644332955, 0.041417904358686766, 0.039861830745551524, 0.042430091706017845, 0.03860582494343624, 0.03792189977506492, 0.04008875199389502, 0.042288318878351695, 0.0451123779732812, 0.046048272710029356, 0.04703972355472392, 0.04755373436865446, 0.05171915330504005, 0.04641792662867825, 0.041048386034591884, 0.04174967346316177, 0.045189799578819204, 0.04437253834653995, 0.04952322819596821, 0.048495220856919215, 0.04687564943710363, 0.045787312534946015, 0.045632461143513055, 0.045985460324898635, 0.04830367673346383, 0.048559423115293376, 0.04634203605010507, 0.042246323554976796, 0.0427993213662828, 0.043118202113991214, 0.0397376581035835, 0.042524906159124125, 0.04396301667115443, 0.04867039684003009, 0.04771072334118596, 0.041454266227342454, 0.04124213206255866, 0.047534682376393966, 0.05078110508940085, 0.04774850181271163, 0.04425430133123059, 0.04154804589053509, 0.038195909231608906, 0.04182649568456345, 0.0474171674946069, 0.04633589844148835, 0.04821749084563349, 0.04593936456256249, 0.04095114998039674, 0.03812693906280689, 0.039941669264537245, 0.04542310406977286, 0.04396482029255902, 0.045118704371473625, 0.043432105984383636, 0.04486230669861454, 0.04510733096996943, 0.04802194422979877, 0.04724024712376392, 0.04495031542359843, 0.0411194520421141, 0.04307633567133907, 0.03851478577185336, 0.043292083853173235, 0.0471114697385979, 0.04410141360899882, 0.04387259467792761, 0.041807039642388265, 0.040043913393272734, 0.04047875080697718, 0.04076540227519291, 0.04352287758531234, 0.041783690696520596, 0.04217269633418065, 0.039689257855054616, 0.04132008608671728, 0.04215640756748706, 0.044254987305214086, 0.047708929166049574, 0.04311592985947795, 0.03982766185296869, 0.04398050285108967, 0.04540210913541616, 0.04354137754695668, 0.04718340520711346, 0.04304256380716064, 0.04512527662191849, 0.04249792496569365, 0.04140133703418367, 0.042107551747422305, 0.04590343220607745, 0.04998206285564875, 0.048629512438683395, 0.04282212524888276, 0.040684721061896245, 0.04174877236176018, 0.04560911999910148, 0.04303492240365023, 0.042734335765768586, 0.04710858231610932, 0.04876592174621365, 0.045945345260402094, 0.042922587054203855, 0.042092430936217035, 0.045067319751666814, 0.04203697983369511, 0.043842431307806455, 0.04074253852244252, 0.03730837043776631, 0.040169986103056096, 0.049958670279766665, 0.051545257540286375, 0.04289900693250782, 0.04231154449535931, 0.048410145694794515, 0.046801763550440134, 0.042237999555196394, 0.041835741015615366, 0.04387174112759666, 0.047713634320126326, 0.044891065593481695, 0.04556497635112717, 0.04794837468728872, 0.04660960032050682, 0.04433117268479811, 0.04023655866424818, 0.04149363768522796, 0.041690792221083656, 0.043203198634303226, 0.046495052241677354, 0.044684251637310594, 0.04334641768609086, 0.04107400430020838, 0.04090743230875665, 0.0443848934695111, 0.04623037541052707, 0.042623150393662705, 0.044217609300187824, 0.04639531848757333, 0.0481520420900174, 0.04718865908878824, 0.04671711156852026, 0.048803227155624906, 0.0490226704377737, 0.0471549241226826, 0.044960104565566684, 0.04360673873038872, 0.041511010429551014, 0.04290965674085523, 0.04129455212589241, 0.04542054094606748, 0.04831137602239049, 0.045346424303656534, 0.04342491848394353, 0.04739216255774936, 0.04439836302730523, 0.04954618021641574, 0.052638705253746816, 0.046792080923896216, 0.04436289139391307, 0.04264844718808553, 0.04550752875794458, 0.04798740860215397, 0.04550803434517415, 0.04724035344319523, 0.04986446736834466, 0.04802973240934296, 0.04455164665564684, 0.044373904837360305, 0.04063588252774808, 0.04144441773364576, 0.03983918746791325, 0.0424804947993889, 0.047018660064111754, 0.04545985815260141, 0.04449947428600084, 0.042099323709878234, 0.04005113353150754, 0.03955788234617767, 0.04071352738557676, 0.045351231046503876, 0.0477126740245017, 0.047229559305357244, 0.04899802107204078, 0.05076153031732595, 0.044629172520261814, 0.04343444733107019, 0.04402179756596247, 0.042221052187320245, 0.04086826498911407, 0.044574422991234826, 0.045323910920803845, 0.040212857615484855, 0.03902977090448637, 0.041132616185642436, 0.04241213087564171, 0.0436542791640779, 0.044083832754888856, 0.04426605151363138, 0.04462255609528038, 0.04061864373614638, 0.0429400707866506, 0.045902528145711774, 0.04381111093869836, 0.04237787542530216, 0.04947183084580828, 0.04697737452930587, 0.044046360463075064, 0.04253220777488787, 0.0407159615018944, 0.04163414451257969, 0.04556149223589182, 0.04569673671244251, 0.043120357929141824, 0.042732979592021124, 0.04061139816703361, 0.03976743918766305, 0.03932856020564626, 0.03885791526589002, 0.04190413329291068, 0.04273352200708448, 0.041647761184739875, 0.040698575125602175, 0.04202658249920415, 0.040229396143461105, 0.04197693015296109, 0.042726557253184416, 0.045699936204365556, 0.04827638244111895, 0.04881570564070077, 0.050102315462853854, 0.04673081519263201, 0.04117653038068325, 0.04072458762791453, 0.043572012575425915, 0.04176665811737533, 0.042316232239207885, 0.046046812137092975, 0.048255447920317794, 0.049033926199807065, 0.04944897226097929, 0.04551393001028015, 0.04242912882554298, 0.04039516413932767, 0.04119694127065856, 0.04451483711829591, 0.04573624747480866, 0.04346421188185337, 0.0427146500856587, 0.04333979841547466, 0.04088438813057513, 0.037859715493757606, 0.038447521143876257, 0.04060762075251992, 0.044562393517217046, 0.046384746321203284, 0.04901567610019312, 0.04777436870239017, 0.042236712292383655, 0.04233795891850927, 0.04249361885394317, 0.039915150412748936, 0.04145143738683157, 0.04350834565541829, 0.045361798682876965, 0.043436424429605755, 0.04605866200153005, 0.04151852810440533, 0.03929083888673877, 0.037223529111601275, 0.03824526783616125, 0.04392140238427932, 0.042612605077156605, 0.04414573067805473, 0.044959265198177216, 0.040570065306508044, 0.04481740267380546, 0.04945561582910571, 0.042954874676428806, 0.041205599739528, 0.043320293153747566, 0.040829074065675694, 0.039884886191473706, 0.04296107292815479, 0.0439308364966333, 0.04054339991629799, 0.04507336256491878, 0.04631568968599266, 0.04439691394795106, 0.04807152200967308, 0.04584385024794673, 0.040044848315188133, 0.03955006413399268, 0.037941584039717896, 0.04354878936943471, 0.05039147049732044, 0.050696086791550014, 0.0469648802481131, 0.04379968031691941, 0.03739583573052905, 0.04031862176647399, 0.04153045398845755, 0.04196750051214071, 0.04318772069700794, 0.041566329133091676, 0.04173469363676586, 0.04686300679110433, 0.0492638134657306, 0.05048186071435193, 0.04821661116257725, 0.0469875352797971, 0.04441032301751122, 0.04449445694351213, 0.04647388594109532, 0.045551552679701506, 0.04158841558460514, 0.044072162346218684, 0.045938123693354344, 0.046697695398172874, 0.046832196059486415, 0.047536791756643645, 0.04905393549702959, 0.04315731085213599, 0.04112724707541983, 0.037245873106063504, 0.03879307057309298, 0.04331530492284225, 0.044344349928532464, 0.04213498535625566, 0.04257757660907315, 0.04323826012079497, 0.048099498781260235, 0.04998518526074056, 0.050520818786158266, 0.04367307746640965, 0.04096701587362801, 0.04563795807946692, 0.04901560326215511, 0.04272935212145646, 0.045657315713128775, 0.048105396916226284, 0.04984215845249612, 0.044952426588283766, 0.04682018163855276, 0.04680326381718331, 0.0469957976367488, 0.04683588375680137, 0.048301950232442, 0.0530498103022339, 0.048927032849742355, 0.03961895410218039, 0.043703846226113424, 0.042712896545018655, 0.04601984596954175, 0.04779591730200789, 0.04506612775799487, 0.047010494757446036, 0.04606848855700907, 0.04357608433849806, 0.0428466969667353, 0.042537402226638224, 0.04537809148259306, 0.04540093417047663, 0.037489112262732466, 0.04455117686583346]],
            'im_sum': [False, 2528.328722114691],
            'regn_sum': [False, 269.9288174510002],
            'npts_real': [True, 6400000],
            'profile': [False, 2.8036859283441817],
            'fit': [False, [3.129263783629313, 6.117903684384051, 5.779428322041025]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.5049652049504],
            'fit_pix': [False, [45.635862765960375, 41.06037101750577]]}



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
            'mask_pix': [False, 9548],
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
            'min_val': [False, -0.1645863801240921],
            'min_val_pos': [True, [34, 36, 0, 949]],
            'im_rms': [False, 0.08653598892456672],
            'im_sum': [False, 2855.507289054503],
            'npts_real': [True, 6400000],
            'fit_0': [False, [1.0955776433750248, 4.020838907836922, 2.970440524540877]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.4453942426872],
            'fit_pix_0': [False, [40.00007098273814, 39.997503447408945]],
            'fit_1': [False, [1.0955887833717917, 4.020278215143762, 2.9701943170887772]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.5063079930342],
            'fit_pix_1': [False, [40.00007124298777, 39.99750480108641]],
            'fit_2': [False, [1.095630862046276, 4.019925393715577, 2.96902170708751]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672217433812],
            'fit_pix_2': [False, [40.000072197166006, 39.99750720003138]]}



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
            'max_val': [False, 0.32576698064804077],
            'max_val_pos': [True, [20, 53, 0, 497]],
            'min_val': [False, -0.33937618136405945],
            'min_val_pos': [True, [32, 47, 0, 491]],
            'im_rms': [False, 0.046717514466101825],
            'im_sum': [False, 330.98844356753955],
            'regn_sum': [False, 51.15742180077359],
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
            'max_val': [False, 1.4210411310195923],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.09896747767925262],
            'min_val_pos': [True, [58, 39, 0, 493]],
            'im_rms': [False, 0.0017599677695868156],
            'im_sum': [False, 55.0920758638531],
            'regn_sum': [False, 55.892402071505785],
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
            'max_val': [False, 1012.2938232421875],
            'max_val_pos': [True, [0, 0, 0, 857]],
            'min_val': [False, 1012.0758666992188],
            'min_val_pos': [True, [0, 0, 0, 475]],
            'im_rms': [False, 1012.157220494919],
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

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[0.0, 3.25])
        self.mom8_creator(img+'.residual', range_list=[0.0, 3.25])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')

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
    # Test 5
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
            perchanweightdensity=False, gridder='standard', \
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
            perchanweightdensity=False, gridder='standard', \
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
    # Test 6
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
            gridder='standard', mosweight=False, \
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
            gridder='standard', mosweight=False, \
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
    # Test 7
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
            perchanweightdensity=False, gridder='standard',  \
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
            perchanweightdensity=False, gridder='standard',  \
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
    # Test 8
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
            perchanweightdensity=False, gridder='standard',  \
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
            perchanweightdensity=False, gridder='standard',  \
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

    # Test 9a
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
            perchanweightdensity=False, gridder='mosaic',  \
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
                gridder='mosaic',  mosweight=True, \
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
                gridder='mosaic',  mosweight=True, \
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
                gridder='mosaic',  mosweight=True, \
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
    # Test 9b
    @stats_dict(test_dict)
    def test_mosaic_cube_pcwdT(self):
        ''' Mosaic cube imaging wth pcwdT - field SMIDGE_NWCloud, spw 22 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name,'test_')+'.iter'
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
            perchanweightdensity=True, gridder='mosaic', \
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
                outframe='LSRK', perchanweightdensity=True, \
                gridder='mosaic', mosweight=True, \
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
                outframe='LSRK', perchanweightdensity=True, \
                gridder='mosaic', mosweight=True, \
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
                outframe='LSRK', perchanweightdensity=True, \
                gridder='mosaic', mosweight=True, \
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

            # retrieve per-channel beam statistics
            bmin_dict, bmaj_dict, pa_dict = \
                self.cube_beam_stats(img+'.psf')

        report0 = th.checkall(imgexist = self.image_list(img, 'mosaic'))

        # .image report (test_mosaic_cube_pcwdT)
        im_stats_dict = self.image_stats(img+'.image', fit_region = \
            'ellipse[[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]', field_regions = \
            ['circle[[00:45:54.383559, -73.15.29.41306], 22.45arcsec]',
             'circle[[00:45:49.435664, -73.15.35.13742], 22.45arcsec]',
             'circle[[00:45:53.057440, -73.15.50.79016], 22.45arcsec]',
             'circle[[00:45:50.762696, -73.15.13.76177], 22.45arcsec]',
             'circle[[00:45:58.006248, -73.15.45.06040], 22.45arcsec]',
             'circle[[00:45:55.708764, -73.15.08.03543], 22.45arcsec]',
             'circle[[00:45:59.330540, -73.15.23.68133], 22.45arcsec]'])

        # test_mosaic_cube_pcwdT.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 8.800156671735735],
            'com_bmin': [False, 6.100528967518075],
            'com_pa': [False, 67.46076115284168],
            'npts': [True, 5925312],
            'npts_unmasked': [False, 3352579.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.2677836418151855],
            'max_val_pos': [True, [45, 38, 0, 252]],
            'min_val': [False, -0.42875778675079346],
            'min_val_pos': [True, [48, 32, 0, 272]],
            'im_rms': [False, 0.08758318024493557],
            'rms_per_chan': [False, [0.08454784229054828, 0.08081074403400043, 0.08214515986480901, 0.0864265424214896, 0.08559577979993536, 0.08267240188831336, 0.080761658995135, 0.08864799515742947, 0.08511980018155442, 0.08710887468182987, 0.07455842794282884, 0.08679442967899591, 0.07672218555279199, 0.08339666396843515, 0.08342785056774209, 0.09508920756603946, 0.0782094630317009, 0.08229799600553456, 0.07865856445846503, 0.09174566220454843, 0.08281923825141245, 0.08850692298599323, 0.08720173656505815, 0.07983263245354182, 0.08901072023943249, 0.08676957195739378, 0.09228845736653081, 0.07952345773483506, 0.08360633847507558, 0.08953635419478691, 0.07687670056583094, 0.08547223871127813, 0.0722237566671663, 0.08122845114951444, 0.08472592112526439, 0.08869096329346654, 0.07975221748089455, 0.08101562173616539, 0.09005674685271724, 0.09064666056412814, 0.08432542487793476, 0.08442355846577512, 0.08225473424350162, 0.08660407244837623, 0.08419525514858209, 0.09092145500811731, 0.07700112553811718, 0.09263511564881605, 0.08996190950262276, 0.08805014739965264, 0.08334781274096709, 0.07995398373475042, 0.08361156044820355, 0.08139675186347044, 0.09640674977612612, 0.07344656015150823, 0.08411367354319836, 0.09181260778377988, 0.08404551028650238, 0.07488449786630803, 0.08726772698002463, 0.07592742634137034, 0.07869515124436191, 0.07713395840673742, 0.0889530068070014, 0.07886056200862707, 0.08433334766591734, 0.08015425283696251, 0.08083995495566454, 0.08232677388380731, 0.09437255030284213, 0.08143245640464152, 0.07951883717941302, 0.08567072593060555, 0.08200455489528632, 0.0754633650915257, 0.085767144139085, 0.08390698936543563, 0.08166615185276183, 0.08294388019809562, 0.08175161359532591, 0.09310134967850235, 0.08956426553496326, 0.08863772275001713, 0.08768742060060232, 0.08730596808772133, 0.09921418996758216, 0.08894272329661541, 0.08525290298222582, 0.09105389161552142, 0.08294466834746969, 0.08279888726356258, 0.0897444459110899, 0.09149609612906127, 0.09736826247293821, 0.0874580064771212, 0.08377256309157666, 0.08827127166827303, 0.08310655247718046, 0.08040840966622335, 0.08279817886158096, 0.08359300310750502, 0.07499876992627746, 0.08207994631296497, 0.08305464424011334, 0.0884947503049501, 0.0925505858545925, 0.08342969809216969, 0.09154929930536125, 0.08641008468572106, 0.09579222473022453, 0.09355982016700347, 0.08179196573941395, 0.08321297682807202, 0.08513449068981359, 0.08097172513845138, 0.08524563237921018, 0.08573522506850681, 0.08499235100075715, 0.08746568480259742, 0.08375366834698306, 0.08454591366146298, 0.0792188660098284, 0.08994735558535912, 0.09777852980854092, 0.09415849395298291, 0.09046545845055788, 0.0842601855414381, 0.09085331779460411, 0.08948416551204, 0.09252892609732252, 0.07885104429804914, 0.08507358084323587, 0.08924576919673373, 0.07632048372672504, 0.08816338783038676, 0.08914822508330522, 0.08521591226020223, 0.08768830495208184, 0.0846907567560094, 0.08245086433336567, 0.09359281278102287, 0.08795141538978407, 0.09133110724322849, 0.07839522749787323, 0.08256068993054838, 0.08444040010558822, 0.08955676262394105, 0.0744924857089921, 0.08285930527824735, 0.07537422704404396, 0.08946040581991314, 0.09225166169714621, 0.08883166824425388, 0.08658120437505346, 0.08078069041521534, 0.08287579235391287, 0.08308256941295859, 0.08279883646854855, 0.08296190701446955, 0.07872964853638084, 0.08409447686471307, 0.08881653564859565, 0.08627522825028244, 0.08466202617969501, 0.08689665906604596, 0.08378386027946005, 0.08593924740670453, 0.08088009767099892, 0.09030660703560082, 0.08574306737360128, 0.09585536044874188, 0.09097544289439698, 0.08739969646413384, 0.08146351183345106, 0.09036429531001434, 0.09594228565098932, 0.0754104244037596, 0.08586515429925357, 0.08548469955576445, 0.08854036826326296, 0.08314347783373752, 0.07538708501494322, 0.08424708401279762, 0.09200924971549221, 0.09298268742698543, 0.08405581086804013, 0.08391297043078556, 0.0874382738242522, 0.0865092157360675, 0.08393760101913277, 0.08423272275147145, 0.08736110502958849, 0.0838603278385726, 0.08713571076049065, 0.08914556563084271, 0.08286535381253005, 0.08941163018320565, 0.08795375788031498, 0.09857995877392249, 0.09038544328988601, 0.08630730158556425, 0.08074711150002968, 0.08833797601773176, 0.09339815243237862, 0.07916245668160474, 0.08722511491702026, 0.08445118836252768, 0.08818398062593365, 0.08991223252324108, 0.08958872745528162, 0.08958386695427573, 0.0859888541494378, 0.08244239643793537, 0.0880315249878136, 0.08511800689012519, 0.08403594070780207, 0.08095531225260892, 0.08964352885887882, 0.08400125252601336, 0.07870701425109543, 0.08068388039662826, 0.08996025461045817, 0.09060555025414195, 0.08710450066403135, 0.08165280951098304, 0.09597457340970073, 0.08848613818411262, 0.0804047532371268, 0.08883036023456667, 0.09002443864984841, 0.09590048065642524, 0.09201433452134399, 0.0784344219766899, 0.08993276472738429, 0.09043764691641651, 0.0950959374542557, 0.09045912506173952, 0.09470292143658518, 0.08309636587828109, 0.08295102264940497, 0.08542679476962445, 0.07967965077132584, 0.09250105979150762, 0.0852278939641754, 0.09389034506206444, 0.08732400651392125, 0.08736604365132138, 0.09405141874915855, 0.10004870526095205, 0.10723503893069163, 0.1257996589564165, 0.1436823034257318, 0.13810223983811945, 0.12066239867371778, 0.12437047875377277, 0.1072957271110957, 0.09041193150401249, 0.09101063129053057, 0.08324128092771393, 0.08942784831142817, 0.08803176186511635, 0.08754740699401242, 0.08812373385693446, 0.08266503446000094, 0.0932862597870741, 0.09509115889831377, 0.09341024923745748, 0.09716746558147825, 0.09371462532304208, 0.08493757979633329, 0.0840823563745246, 0.09225289907267915, 0.09721690331208403, 0.09102532563372895, 0.09191209341084346, 0.08222408950279445, 0.09549165468015078, 0.0828588141432136, 0.08621846741135158, 0.08102441646313648, 0.08528151685496348, 0.08971436016606618, 0.08495831892979763, 0.0751464983385108, 0.0825741797203384, 0.09301501446208484, 0.07820718000462748, 0.07851985599497896, 0.08026537559234427, 0.0853785181946889, 0.0841522267503986, 0.0902730950775931, 0.09320813048341908, 0.07921760517505672, 0.09120259125276349, 0.08662028051854927, 0.09319862578946961, 0.07917380885021928, 0.09081458871341828, 0.08413009427107505, 0.08075357399784142, 0.09038395373240243, 0.08884888986855069, 0.08580177025802656, 0.08358774933196578, 0.08616012189921222, 0.08907622261398132, 0.09093535732495094, 0.08410001063840633, 0.08836017031066742, 0.08462812202527752, 0.08166111457676953, 0.08952117220794778, 0.09485868771459792, 0.09079675620681926, 0.08093656378261929, 0.08593405117838876, 0.08922950899879277, 0.09227057153668085, 0.09401743722118289, 0.0928745289013179, 0.09934734172329057, 0.08540478186283113, 0.08535868985033361, 0.07921735506645071, 0.0811612544987567, 0.09309982858993383, 0.07275034569377031, 0.0834881145891589, 0.09304592981457883, 0.09545763433561774, 0.09168474196164955, 0.09653939179768366, 0.0779519101763873, 0.09053680873231403, 0.08730008263451008, 0.08039772102057097, 0.08585683588093039, 0.08556787357829983, 0.0947936311691715, 0.08246248891415837, 0.08792661871208168, 0.08015270038923986, 0.08256943387674712, 0.08927097573624211, 0.09207091351089977, 0.0896824859830704, 0.09234003709632913, 0.09463882111960684, 0.09137188816666182, 0.07969475057664313, 0.09414973400567359, 0.07531585580072055, 0.09009555546093831, 0.08235003628877575, 0.0911614507907243, 0.08482314118710578, 0.09094087612722215, 0.09298430761879853, 0.089550761769171, 0.0930529724204851, 0.09083286908816636, 0.08872869703344648, 0.09002174330250973, 0.08841686476596793, 0.0918771091483413, 0.08646475892780502, 0.09077999744273917, 0.0902324205734482, 0.09349418892452889, 0.08810103522129416, 0.08291325070003897, 0.07852868976600987, 0.08775973676251207, 0.07902683588925725, 0.09598239811752228, 0.07948624599332443, 0.08409186031395341, 0.0914668197499112, 0.0903906178514704, 0.0842839858889817, 0.08398783251971939, 0.08441879571202276, 0.08183365803853736, 0.08893715844011056, 0.07759090752716989, 0.08972723424366359, 0.08265048295523349, 0.092849420659467, 0.08179774845604361, 0.08685230205561319, 0.07989560032838225, 0.07796989239570334, 0.09311637874695278, 0.08696577375107635, 0.08951796642941898, 0.09463860562188808, 0.0839413707332889, 0.09297955537921224, 0.09367506665625282, 0.09089445993342953, 0.08581530413616678, 0.0902762324619946, 0.08640085730063643, 0.08823854852055034, 0.08981444309815341, 0.08565107832544971, 0.09193775702658576, 0.09260805786371104, 0.0921952197320759, 0.08452376050579682, 0.09328092342673133, 0.08562765496362403, 0.09032549920735691, 0.08798132852578897, 0.08905222099300897, 0.0867767542165491, 0.09194697859135229, 0.08932513698442741, 0.08889527402694361, 0.0840919921506925, 0.0836052486963255, 0.08727126227483044, 0.09209250380191883, 0.0844537628792832, 0.08297890500129225, 0.08178118773838311, 0.08459516026822875, 0.08053761474597218, 0.08381233835627229, 0.0839457466542763, 0.07949608290558277, 0.0904509018877568, 0.08200764321569995, 0.08050139462870623, 0.08720947157941367, 0.08591229399475535, 0.09018745443845397, 0.07849635318354746, 0.09387222104195914, 0.08209109556354426, 0.09124711008141882, 0.0905665380196555, 0.08085067430598548, 0.09119973417058405, 0.09155593055648911, 0.1010284928807239, 0.08970058474320176, 0.08232882192986131, 0.08398569353110392, 0.0877754526883052, 0.09208999746494177, 0.09472528023494405, 0.09034240395103484, 0.08518087490637068, 0.08350695964551678, 0.09209472883427205, 0.09036402807277293, 0.08296471302270847, 0.0882127589998411, 0.09355184833535637, 0.10061460079849864, 0.08536691235438341, 0.08399933918812515, 0.08712106652079493, 0.08387947393796223, 0.08737767801657396, 0.08931061517222441, 0.08712614309972459, 0.08601226423648409, 0.08938036440801743, 0.08869488330549792, 0.09355910147394612, 0.08626686052481321, 0.08448212598521213, 0.08859708838016166, 0.0971088804234492, 0.09103845664628896, 0.0849462566840279, 0.08758063749841716, 0.09067273697131494, 0.0839629156095452, 0.09617762805426901, 0.072522464162689, 0.0880036422297596, 0.0865640015743217, 0.09185579403530092, 0.09340925578879293, 0.07504254317822588, 0.09081850337601476, 0.09661397296295061, 0.09051672239952115, 0.10366343434446765, 0.08147205490776653, 0.09872900868720667, 0.07845316742855823, 0.09131002219444805, 0.08882995083036258, 0.0921941802749485, 0.09639031609999545, 0.08889271297115287, 0.08985346497743775, 0.09422387248199929, 0.08380520336843732, 0.09670680610402677, 0.08310685183674356, 0.08782985429831294]],
            'im_sum': [False, 148.72369868747853],
            'regn_sum': [False, 75.70919078961015],
            'npts_real': [True, 5925312],
            'rms_per_field': [False, [0.08869172867021574, 0.08749690204094095, 0.08899766355856122, 0.08697053128120366, 0.08862283721987575, 0.08720280700494934, 0.08755366033276048]],
            'profile': [False, 1.2530369225508722],
            'fit': [False, [1.2785385831665967, 8.996787448271887, 7.751896012285196]],
            'fit_loc_chan': [True, 252],
            'fit_loc_freq': [1e-10, 220.31420623259388],
            'fit_pix': [False, [44.60565437474867, 38.21219207866642]]}


        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_cube_pcwdT)
        mask_stats_dict = self.image_stats(img+'.mask')

        # test_mosaic_cube_pcwdT.exp_mask_stats
        exp_mask_stats = {'npts': [True, 5925312],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'mask_pix': [False, 4105],
            'mask_regns': [True, 1],
            'npts_real': [True, 5925312]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube_pcwdT)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        # test_mosaic_cube_pcwdT.exp_pb_stats
        exp_pb_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3352579.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, 0.20000170171260834],
            'im_rms': [False, 0.6139604304840641],
            'npts_0.2': [False, [6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600]],
            'npts_0.5': [False, [3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3598, 3598, 3598, 3598, 3598, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599]],
            'npts_real': [True, 5925312],
            'fit': [False, [1.079756729983638, 69.32945924752613, 69.27414584027905]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [54.058144958203236, 53.9848903021958]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_pcwdT)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        # test_mosaic_cube_pcwdT.exp_psf_stats
        exp_psf_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, -0.16422684490680695],
            'min_val_pos': [True, [63, 52, 0, 381]],
            'im_rms': [False, 0.06139080806163987],
            'im_sum': [False, 27.159552334707996],
            'npts_real': [True, 5925312],
            'fit_0': [False, [1.0996072179514504, 7.80778199842682, 5.287148667879593]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [53.99103608954199, 54.00525892337048]],
            'fit_1': [False, [1.0996013414468444, 7.809207265062828, 5.287348602363145]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [53.990709312925155, 54.00474442800999]],
            'fit_2': [False, [1.0995499572167458, 7.8001792765571105, 5.293446765535396]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [53.991293887361174, 54.00416063682335]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_pcwdT)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        # test_mosaic_cube_pcwdT.exp_resid_stats
        exp_resid_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3352579.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.45122435688972473],
            'max_val_pos': [True, [48, 60, 0, 249]],
            'min_val': [False, -0.42875778675079346],
            'min_val_pos': [True, [48, 32, 0, 272]],
            'im_rms': [False, 0.08724054766799626],
            'im_sum': [False, -139.18306628067614],
            'regn_sum': [False, 20.81485258124303],
            'npts_real': [True, 5925312]}


        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_pcwdT)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_pcwdT.exp_model_stats
        exp_model_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.5603721737861633],
            'max_val_pos': [True, [45, 39, 0, 251]],
            'min_val': [False, -0.08046908676624298],
            'min_val_pos': [True, [46, 57, 0, 255]],
            'im_rms': [False, 0.0005318959551228511],
            'im_sum': [False, 5.726837918162346],
            'regn_sum': [False, 5.726837918162346],
            'mask_non0': [True, 0],
            'npts_real': [True, 5925312]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 120.96144104003906],
            'max_val_pos': [True, [0, 0, 0, 336]],
            'min_val': [False, 120.7746353149414],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 120.86384620156764],
            'npts_real': [True, 508]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube_pcwdT.exp_wt_stats
        exp_wt_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.39721736311912537],
            'max_val_pos': [True, [54, 54, 0, 381]],
            'min_val': [False, 7.159297092584893e-05],
            'im_rms': [False, 0.14206311600056568],
            'im_sum': [False, 511866.22841321985],
            'npts_0.2': [False, [6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600]],
            'npts_0.5': [False, [3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3598, 3598, 3598, 3598, 3598, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599]],
            'npts_real': [True, 5925312]}
        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            # test_mosaic_cube_pcwdT.exp_bmin_dict
            exp_bmin_dict = {'*0': 6.0930256843566895,'*1': 6.0930256843566895,'*10': 6.093120098114014,'*100': 6.093204021453857,'*101': 6.093204021453857,'*102': 6.093204021453857,'*103': 6.093204021453857,'*104': 6.093204021453857,'*105': 6.093204021453857,'*106': 6.093204021453857,'*107': 6.093108177185059,'*108': 6.093108177185059,'*109': 6.093108177185059,'*11': 6.093120098114014,'*110': 6.093108177185059,'*111': 6.093087673187256,'*112': 6.091174125671387,'*113': 6.091147422790527,'*114': 6.091147422790527,'*115': 6.091147422790527,'*116': 6.091147422790527,'*117': 6.0911359786987305,'*118': 6.0911030769348145,'*119': 6.0911030769348145,'*12': 6.093110084533691,'*120': 6.090994834899902,'*121': 6.090970516204834,'*122': 6.090970516204834,'*123': 6.090970516204834,'*124': 6.090970516204834,'*125': 6.090970516204834,'*126': 6.090970516204834,'*127': 6.090961933135986,'*128': 6.090961933135986,'*129': 6.090953350067139,'*13': 6.093110084533691,'*130': 6.090953350067139,'*131': 6.090926170349121,'*132': 6.090926170349121,'*133': 6.090926170349121,'*134': 6.090926170349121,'*135': 6.090926170349121,'*136': 6.090926170349121,'*137': 6.090806484222412,'*138': 6.090806484222412,'*139': 6.090806484222412,'*14': 6.093110084533691,'*140': 6.090806484222412,'*141': 6.090785026550293,'*142': 6.090785026550293,'*143': 6.090785026550293,'*144': 6.090785026550293,'*145': 6.090785026550293,'*146': 6.090779781341553,'*147': 6.090779781341553,'*148': 6.090779781341553,'*149': 6.090763568878174,'*15': 6.093103885650635,'*150': 6.090763568878174,'*151': 6.090763568878174,'*152': 6.090755462646484,'*153': 6.090755462646484,'*154': 6.090744495391846,'*155': 6.090744495391846,'*156': 6.090744495391846,'*157': 6.090744495391846,'*158': 6.0907111167907715,'*159': 6.0907111167907715,'*16': 6.093103885650635,'*160': 6.0907111167907715,'*161': 6.0907111167907715,'*162': 6.0907111167907715,'*163': 6.090665340423584,'*164': 6.090665340423584,'*165': 6.090665340423584,'*166': 6.090665340423584,'*167': 6.090665340423584,'*168': 6.090653419494629,'*169': 6.090653419494629,'*17': 6.093103885650635,'*170': 6.090653419494629,'*171': 6.090653419494629,'*172': 6.090628623962402,'*173': 6.090611457824707,'*174': 6.090585231781006,'*175': 6.090578556060791,'*176': 6.0897040367126465,'*177': 6.090185642242432,'*178': 6.090174198150635,'*179': 6.090174198150635,'*18': 6.0930891036987305,'*180': 6.090174198150635,'*181': 6.090117454528809,'*182': 6.090091228485107,'*183': 6.090091228485107,'*184': 6.090091228485107,'*185': 6.090091228485107,'*186': 6.090091228485107,'*187': 6.090091228485107,'*188': 6.090091228485107,'*189': 6.09009313583374,'*19': 6.0930891036987305,'*190': 6.09009313583374,'*191': 6.09002685546875,'*192': 6.09002685546875,'*193': 6.090029716491699,'*194': 6.090010166168213,'*195': 6.090010166168213,'*196': 6.090010166168213,'*197': 6.090006351470947,'*198': 6.090006351470947,'*199': 6.090006351470947,'*2': 6.0930256843566895,'*20': 6.0930891036987305,'*200': 6.090010643005371,'*201': 6.090010643005371,'*202': 6.090010643005371,'*203': 6.090010643005371,'*204': 6.083273887634277,'*205': 6.083273887634277,'*206': 6.085849761962891,'*207': 6.085849761962891,'*208': 6.086237907409668,'*209': 6.086216926574707,'*21': 6.0930891036987305,'*210': 6.086216926574707,'*211': 6.086202144622803,'*212': 6.08619499206543,'*213': 6.0861496925354,'*214': 6.0861496925354,'*215': 6.086057662963867,'*216': 6.086057662963867,'*217': 6.086057662963867,'*218': 6.086057662963867,'*219': 6.085986137390137,'*22': 6.093074321746826,'*220': 6.085943222045898,'*221': 6.085077285766602,'*222': 6.085077285766602,'*223': 6.085077285766602,'*224': 6.085077285766602,'*225': 6.086266994476318,'*226': 6.086266994476318,'*227': 6.086266994476318,'*228': 6.086266994476318,'*229': 6.086266994476318,'*23': 6.093074321746826,'*230': 6.086266994476318,'*231': 6.0862321853637695,'*232': 6.0862321853637695,'*233': 6.0862321853637695,'*234': 6.0862321853637695,'*235': 6.086203098297119,'*236': 6.086203098297119,'*237': 6.0861616134643555,'*238': 6.08617639541626,'*239': 6.08617639541626,'*24': 6.093074321746826,'*240': 6.086090087890625,'*241': 6.086090087890625,'*242': 6.086090087890625,'*243': 6.086079120635986,'*244': 6.086079120635986,'*245': 6.086079120635986,'*246': 6.086080551147461,'*247': 6.08603048324585,'*248': 6.08603048324585,'*249': 6.08603048324585,'*25': 6.093074321746826,'*250': 6.086004257202148,'*251': 6.086004257202148,'*252': 6.0859761238098145,'*253': 6.085971355438232,'*254': 6.085971355438232,'*255': 6.085971355438232,'*256': 6.08592414855957,'*257': 6.08592414855957,'*258': 6.08592414855957,'*259': 6.085926055908203,'*26': 6.093074321746826,'*260': 6.087209701538086,'*261': 6.087515830993652,'*262': 6.087515830993652,'*263': 6.0874738693237305,'*264': 6.0874738693237305,'*265': 6.0874738693237305,'*266': 6.0874738693237305,'*267': 6.0874738693237305,'*268': 6.0874738693237305,'*269': 6.087475299835205,'*27': 6.093019485473633,'*270': 6.087475299835205,'*271': 6.087457180023193,'*272': 6.0874528884887695,'*273': 6.0874528884887695,'*274': 6.0874528884887695,'*275': 6.0874528884887695,'*276': 6.087437152862549,'*277': 6.087437152862549,'*278': 6.087437152862549,'*279': 6.087437152862549,'*28': 6.093019485473633,'*280': 6.0874199867248535,'*281': 6.0874199867248535,'*282': 6.087393283843994,'*283': 6.087385177612305,'*284': 6.087385177612305,'*285': 6.087385177612305,'*286': 6.087225914001465,'*287': 6.087218761444092,'*288': 6.087218761444092,'*289': 6.0871663093566895,'*29': 6.093006610870361,'*290': 6.086967468261719,'*291': 6.086967468261719,'*292': 6.086978435516357,'*293': 6.086978435516357,'*294': 6.086978435516357,'*295': 6.086970329284668,'*296': 6.086970329284668,'*297': 6.086970329284668,'*298': 6.086947441101074,'*299': 6.086947441101074,'*3': 6.0929975509643555,'*30': 6.093006610870361,'*300': 6.086947441101074,'*301': 6.086947441101074,'*302': 6.086947441101074,'*303': 6.086921691894531,'*304': 6.086921691894531,'*305': 6.086921691894531,'*306': 6.086921691894531,'*307': 6.086921691894531,'*308': 6.0869035720825195,'*309': 6.086895942687988,'*31': 6.092994213104248,'*310': 6.086895942687988,'*311': 6.086895942687988,'*312': 6.086895942687988,'*313': 6.086895942687988,'*314': 6.086895942687988,'*315': 6.0868377685546875,'*316': 6.08682107925415,'*317': 6.08682107925415,'*318': 6.0868120193481445,'*319': 6.0868120193481445,'*32': 6.093119144439697,'*320': 6.0868120193481445,'*321': 6.0868120193481445,'*322': 6.0868120193481445,'*323': 6.0868120193481445,'*324': 6.0868120193481445,'*325': 6.0868120193481445,'*326': 6.0868120193481445,'*327': 6.0868120193481445,'*328': 6.0868120193481445,'*329': 6.0868120193481445,'*33': 6.093119144439697,'*330': 6.0868120193481445,'*331': 6.086823463439941,'*332': 6.086823463439941,'*333': 6.086823463439941,'*334': 6.086823463439941,'*335': 6.086823463439941,'*336': 6.086922645568848,'*337': 6.086922645568848,'*338': 6.086894989013672,'*339': 6.086889743804932,'*34': 6.093119144439697,'*340': 6.086889743804932,'*341': 6.086889743804932,'*342': 6.086889743804932,'*343': 6.086889743804932,'*344': 6.086954116821289,'*345': 6.086935997009277,'*346': 6.086904525756836,'*347': 6.086871147155762,'*348': 6.086841583251953,'*349': 6.086841583251953,'*35': 6.093148231506348,'*350': 6.08683443069458,'*351': 6.08683443069458,'*352': 6.08683443069458,'*353': 6.08683443069458,'*354': 6.086830139160156,'*355': 6.086830139160156,'*356': 6.086830139160156,'*357': 6.086794853210449,'*358': 6.086794853210449,'*359': 6.086792469024658,'*36': 6.093117713928223,'*360': 6.0867600440979,'*361': 6.086698055267334,'*362': 6.086698055267334,'*363': 6.086672782897949,'*364': 6.086672782897949,'*365': 6.086672782897949,'*366': 6.086672782897949,'*367': 6.086653232574463,'*368': 6.086653232574463,'*369': 6.092501640319824,'*37': 6.093286037445068,'*370': 6.092501640319824,'*371': 6.092469692230225,'*372': 6.092469692230225,'*373': 6.0918426513671875,'*374': 6.091516017913818,'*375': 6.091516017913818,'*376': 6.091516017913818,'*377': 6.091516017913818,'*378': 6.091509819030762,'*379': 6.09150505065918,'*38': 6.093286037445068,'*380': 6.091473579406738,'*381': 6.091257095336914,'*382': 6.0920186042785645,'*383': 6.0920186042785645,'*384': 6.091982364654541,'*385': 6.091982364654541,'*386': 6.091970920562744,'*387': 6.091944217681885,'*388': 6.091914653778076,'*389': 6.091904163360596,'*39': 6.093268394470215,'*390': 6.091904163360596,'*391': 6.091904163360596,'*392': 6.091904163360596,'*393': 6.091904163360596,'*394': 6.091904163360596,'*395': 6.0919036865234375,'*396': 6.0919036865234375,'*397': 6.091885089874268,'*398': 6.091878414154053,'*399': 6.091878414154053,'*4': 6.0929975509643555,'*40': 6.093268394470215,'*400': 6.091878414154053,'*401': 6.091878414154053,'*402': 6.091878414154053,'*403': 6.091833114624023,'*404': 6.0918145179748535,'*405': 6.0918145179748535,'*406': 6.091785907745361,'*407': 6.091785907745361,'*408': 6.091785907745361,'*409': 6.091785907745361,'*41': 6.093268394470215,'*410': 6.091785907745361,'*411': 6.091785907745361,'*412': 6.091785907745361,'*413': 6.091761112213135,'*414': 6.091761112213135,'*415': 6.09174919128418,'*416': 6.09174919128418,'*417': 6.09174919128418,'*418': 6.09174919128418,'*419': 6.09174919128418,'*42': 6.093251705169678,'*420': 6.09174919128418,'*421': 6.09174919128418,'*422': 6.09174919128418,'*423': 6.09174919128418,'*424': 6.09174919128418,'*425': 6.09174919128418,'*426': 6.091667652130127,'*427': 6.0916428565979,'*428': 6.091669082641602,'*429': 6.091669082641602,'*43': 6.093251705169678,'*430': 6.091669082641602,'*431': 6.0916523933410645,'*432': 6.0916523933410645,'*433': 6.091625690460205,'*434': 6.091625690460205,'*435': 6.091625690460205,'*436': 6.091625690460205,'*437': 6.091625690460205,'*438': 6.091625690460205,'*439': 6.091625690460205,'*44': 6.093251705169678,'*440': 6.091555118560791,'*441': 6.091555118560791,'*442': 6.091555118560791,'*443': 6.091555118560791,'*444': 6.091536998748779,'*445': 6.091536998748779,'*446': 6.091536998748779,'*447': 6.091536998748779,'*448': 6.091518402099609,'*449': 6.091518402099609,'*45': 6.093251705169678,'*450': 6.091518402099609,'*451': 6.091506481170654,'*452': 6.091506481170654,'*453': 6.091506481170654,'*454': 6.091459274291992,'*455': 6.091434478759766,'*456': 6.091385364532471,'*457': 6.091385364532471,'*458': 6.091385364532471,'*459': 6.091338157653809,'*46': 6.093251705169678,'*460': 6.091338157653809,'*461': 6.091400146484375,'*462': 6.091400146484375,'*463': 6.091400146484375,'*464': 6.091400146484375,'*465': 6.091400146484375,'*466': 6.091400146484375,'*467': 6.091400146484375,'*468': 6.091400146484375,'*469': 6.091400146484375,'*47': 6.093205451965332,'*470': 6.091400146484375,'*471': 6.091400146484375,'*472': 6.091400146484375,'*473': 6.091400146484375,'*474': 6.091400146484375,'*475': 6.091400146484375,'*476': 6.091400146484375,'*477': 6.091357707977295,'*478': 6.091365337371826,'*479': 6.091360569000244,'*48': 6.09317684173584,'*480': 6.091360569000244,'*481': 6.091360569000244,'*482': 6.091360569000244,'*483': 6.091335296630859,'*484': 6.091335296630859,'*485': 6.091316223144531,'*486': 6.091316223144531,'*487': 6.091316223144531,'*488': 6.091316223144531,'*489': 6.091318130493164,'*49': 6.093170642852783,'*490': 6.091318130493164,'*491': 6.091318130493164,'*492': 6.091318130493164,'*493': 6.091318130493164,'*494': 6.091293811798096,'*495': 6.091293811798096,'*496': 6.091269493103027,'*497': 6.091269493103027,'*498': 6.091269493103027,'*499': 6.0912652015686035,'*5': 6.0929975509643555,'*50': 6.093170642852783,'*500': 6.0912652015686035,'*501': 6.0912652015686035,'*502': 6.0912652015686035,'*503': 6.0912652015686035,'*504': 6.091235637664795,'*505': 6.091235637664795,'*506': 6.091235637664795,'*507': 6.091235637664795,'*51': 6.093170642852783,'*52': 6.093170642852783,'*53': 6.093170642852783,'*54': 6.093170642852783,'*55': 6.093170642852783,'*56': 6.093170642852783,'*57': 6.093170642852783,'*58': 6.0931549072265625,'*59': 6.0931549072265625,'*6': 6.0931525230407715,'*60': 6.0931549072265625,'*61': 6.0931549072265625,'*62': 6.093143939971924,'*63': 6.093143939971924,'*64': 6.093143939971924,'*65': 6.0931315422058105,'*66': 6.0931315422058105,'*67': 6.0931315422058105,'*68': 6.0931315422058105,'*69': 6.0931315422058105,'*7': 6.0931525230407715,'*70': 6.093120098114014,'*71': 6.093102931976318,'*72': 6.093102931976318,'*73': 6.093114852905273,'*74': 6.093114852905273,'*75': 6.093114852905273,'*76': 6.093114852905273,'*77': 6.09310245513916,'*78': 6.093089580535889,'*79': 6.093089580535889,'*8': 6.0931525230407715,'*80': 6.093089580535889,'*81': 6.093079566955566,'*82': 6.093079566955566,'*83': 6.093079566955566,'*84': 6.093079566955566,'*85': 6.093079566955566,'*86': 6.093082427978516,'*87': 6.093417167663574,'*88': 6.093417167663574,'*89': 6.093417167663574,'*9': 6.0931525230407715,'*90': 6.093415260314941,'*91': 6.093415260314941,'*92': 6.093415260314941,'*93': 6.093415260314941,'*94': 6.093403339385986,'*95': 6.093403339385986,'*96': 6.093403339385986,'*97': 6.093393325805664,'*98': 6.093393325805664,'*99': 6.093204021453857}

            # test_mosaic_cube_pcwdT.exp_bmaj_dict
            exp_bmaj_dict = {'*0': 8.785683631896973,'*1': 8.785683631896973,'*10': 8.785573959350586,'*100': 8.786060333251953,'*101': 8.786060333251953,'*102': 8.786060333251953,'*103': 8.786060333251953,'*104': 8.786060333251953,'*105': 8.786060333251953,'*106': 8.786060333251953,'*107': 8.786009788513184,'*108': 8.786009788513184,'*109': 8.786009788513184,'*11': 8.785573959350586,'*110': 8.786009788513184,'*111': 8.786016464233398,'*112': 8.789488792419434,'*113': 8.789477348327637,'*114': 8.789477348327637,'*115': 8.789477348327637,'*116': 8.789477348327637,'*117': 8.78940486907959,'*118': 8.789399147033691,'*119': 8.789399147033691,'*12': 8.78558349609375,'*120': 8.789487838745117,'*121': 8.789450645446777,'*122': 8.789450645446777,'*123': 8.789450645446777,'*124': 8.789450645446777,'*125': 8.789450645446777,'*126': 8.789450645446777,'*127': 8.78946304321289,'*128': 8.78946304321289,'*129': 8.789395332336426,'*13': 8.78558349609375,'*130': 8.789395332336426,'*131': 8.78938102722168,'*132': 8.78938102722168,'*133': 8.78938102722168,'*134': 8.78938102722168,'*135': 8.78938102722168,'*136': 8.78938102722168,'*137': 8.789447784423828,'*138': 8.789447784423828,'*139': 8.789447784423828,'*14': 8.78558349609375,'*140': 8.789447784423828,'*141': 8.789440155029297,'*142': 8.789440155029297,'*143': 8.789440155029297,'*144': 8.789440155029297,'*145': 8.789440155029297,'*146': 8.78938102722168,'*147': 8.78938102722168,'*148': 8.78938102722168,'*149': 8.789372444152832,'*15': 8.78553581237793,'*150': 8.789372444152832,'*151': 8.789372444152832,'*152': 8.78931713104248,'*153': 8.78931713104248,'*154': 8.789324760437012,'*155': 8.789324760437012,'*156': 8.789324760437012,'*157': 8.789324760437012,'*158': 8.789305686950684,'*159': 8.789305686950684,'*16': 8.78553581237793,'*160': 8.789305686950684,'*161': 8.789305686950684,'*162': 8.789305686950684,'*163': 8.789173126220703,'*164': 8.789173126220703,'*165': 8.789173126220703,'*166': 8.789173126220703,'*167': 8.789173126220703,'*168': 8.789041519165039,'*169': 8.789041519165039,'*17': 8.78553581237793,'*170': 8.789041519165039,'*171': 8.789041519165039,'*172': 8.789009094238281,'*173': 8.788891792297363,'*174': 8.788886070251465,'*175': 8.788853645324707,'*176': 8.789477348327637,'*177': 8.788763999938965,'*178': 8.788691520690918,'*179': 8.788691520690918,'*18': 8.785530090332031,'*180': 8.788691520690918,'*181': 8.788743019104004,'*182': 8.788736343383789,'*183': 8.788736343383789,'*184': 8.788736343383789,'*185': 8.788736343383789,'*186': 8.788736343383789,'*187': 8.788736343383789,'*188': 8.788736343383789,'*189': 8.788628578186035,'*19': 8.785530090332031,'*190': 8.788628578186035,'*191': 8.788630485534668,'*192': 8.788630485534668,'*193': 8.788564682006836,'*194': 8.788558959960938,'*195': 8.788558959960938,'*196': 8.788558959960938,'*197': 8.788519859313965,'*198': 8.788519859313965,'*199': 8.788519859313965,'*2': 8.785683631896973,'*20': 8.785530090332031,'*200': 8.788466453552246,'*201': 8.788466453552246,'*202': 8.788466453552246,'*203': 8.788466453552246,'*204': 8.791274070739746,'*205': 8.791274070739746,'*206': 8.789750099182129,'*207': 8.789750099182129,'*208': 8.790533065795898,'*209': 8.790518760681152,'*21': 8.785530090332031,'*210': 8.790518760681152,'*211': 8.790525436401367,'*212': 8.790457725524902,'*213': 8.790474891662598,'*214': 8.790474891662598,'*215': 8.790464401245117,'*216': 8.790464401245117,'*217': 8.790464401245117,'*218': 8.790464401245117,'*219': 8.790465354919434,'*22': 8.785497665405273,'*220': 8.790297508239746,'*221': 8.791502952575684,'*222': 8.791502952575684,'*223': 8.791502952575684,'*224': 8.791502952575684,'*225': 8.788743019104004,'*226': 8.788743019104004,'*227': 8.788743019104004,'*228': 8.788743019104004,'*229': 8.788743019104004,'*23': 8.785497665405273,'*230': 8.788743019104004,'*231': 8.78868579864502,'*232': 8.78868579864502,'*233': 8.78868579864502,'*234': 8.78868579864502,'*235': 8.78870964050293,'*236': 8.78870964050293,'*237': 8.788731575012207,'*238': 8.788661003112793,'*239': 8.788661003112793,'*24': 8.785497665405273,'*240': 8.788551330566406,'*241': 8.788551330566406,'*242': 8.788551330566406,'*243': 8.788514137268066,'*244': 8.788514137268066,'*245': 8.788514137268066,'*246': 8.788484573364258,'*247': 8.788549423217773,'*248': 8.788549423217773,'*249': 8.788549423217773,'*25': 8.785497665405273,'*250': 8.788538932800293,'*251': 8.788538932800293,'*252': 8.788518905639648,'*253': 8.788459777832031,'*254': 8.788459777832031,'*255': 8.788459777832031,'*256': 8.78845500946045,'*257': 8.78845500946045,'*258': 8.78845500946045,'*259': 8.788349151611328,'*26': 8.785497665405273,'*260': 8.787582397460938,'*261': 8.787413597106934,'*262': 8.787413597106934,'*263': 8.787412643432617,'*264': 8.787412643432617,'*265': 8.787412643432617,'*266': 8.787412643432617,'*267': 8.787412643432617,'*268': 8.787412643432617,'*269': 8.787325859069824,'*27': 8.785531997680664,'*270': 8.787325859069824,'*271': 8.787242889404297,'*272': 8.78722858428955,'*273': 8.78722858428955,'*274': 8.78722858428955,'*275': 8.78722858428955,'*276': 8.787214279174805,'*277': 8.787214279174805,'*278': 8.787214279174805,'*279': 8.787214279174805,'*28': 8.785531997680664,'*280': 8.787205696105957,'*281': 8.787205696105957,'*282': 8.787191390991211,'*283': 8.787128448486328,'*284': 8.787128448486328,'*285': 8.787128448486328,'*286': 8.787383079528809,'*287': 8.78729248046875,'*288': 8.78729248046875,'*289': 8.787189483642578,'*29': 8.785526275634766,'*290': 8.784912109375,'*291': 8.784912109375,'*292': 8.78486156463623,'*293': 8.78486156463623,'*294': 8.78486156463623,'*295': 8.784805297851562,'*296': 8.784805297851562,'*297': 8.784805297851562,'*298': 8.78474235534668,'*299': 8.78474235534668,'*3': 8.785660743713379,'*30': 8.785526275634766,'*300': 8.78474235534668,'*301': 8.78474235534668,'*302': 8.78474235534668,'*303': 8.784736633300781,'*304': 8.784736633300781,'*305': 8.784736633300781,'*306': 8.784736633300781,'*307': 8.784736633300781,'*308': 8.784725189208984,'*309': 8.784669876098633,'*31': 8.785492897033691,'*310': 8.784669876098633,'*311': 8.784669876098633,'*312': 8.784669876098633,'*313': 8.784669876098633,'*314': 8.784669876098633,'*315': 8.784721374511719,'*316': 8.784637451171875,'*317': 8.784637451171875,'*318': 8.784601211547852,'*319': 8.784601211547852,'*32': 8.785429000854492,'*320': 8.784601211547852,'*321': 8.784601211547852,'*322': 8.784601211547852,'*323': 8.784601211547852,'*324': 8.784601211547852,'*325': 8.784601211547852,'*326': 8.784601211547852,'*327': 8.784601211547852,'*328': 8.784601211547852,'*329': 8.784601211547852,'*33': 8.785429000854492,'*330': 8.784601211547852,'*331': 8.784573554992676,'*332': 8.784573554992676,'*333': 8.784573554992676,'*334': 8.784573554992676,'*335': 8.784573554992676,'*336': 8.784466743469238,'*337': 8.784466743469238,'*338': 8.78449821472168,'*339': 8.784478187561035,'*34': 8.785429000854492,'*340': 8.784478187561035,'*341': 8.784478187561035,'*342': 8.784478187561035,'*343': 8.784478187561035,'*344': 8.784396171569824,'*345': 8.784385681152344,'*346': 8.78436279296875,'*347': 8.784391403198242,'*348': 8.784390449523926,'*349': 8.784390449523926,'*35': 8.785266876220703,'*350': 8.784299850463867,'*351': 8.784299850463867,'*352': 8.784299850463867,'*353': 8.784299850463867,'*354': 8.784197807312012,'*355': 8.784197807312012,'*356': 8.784197807312012,'*357': 8.78419017791748,'*358': 8.78419017791748,'*359': 8.784122467041016,'*36': 8.785256385803223,'*360': 8.784140586853027,'*361': 8.784163475036621,'*362': 8.784163475036621,'*363': 8.784181594848633,'*364': 8.784181594848633,'*365': 8.784181594848633,'*366': 8.784181594848633,'*367': 8.784175872802734,'*368': 8.784175872802734,'*369': 8.780726432800293,'*37': 8.785148620605469,'*370': 8.780726432800293,'*371': 8.780641555786133,'*372': 8.780641555786133,'*373': 8.778660774230957,'*374': 8.778858184814453,'*375': 8.778858184814453,'*376': 8.778858184814453,'*377': 8.778858184814453,'*378': 8.778841972351074,'*379': 8.778786659240723,'*38': 8.785148620605469,'*380': 8.778735160827637,'*381': 8.77867317199707,'*382': 8.78063678741455,'*383': 8.78063678741455,'*384': 8.780657768249512,'*385': 8.780657768249512,'*386': 8.780573844909668,'*387': 8.780562400817871,'*388': 8.780559539794922,'*389': 8.780515670776367,'*39': 8.785120964050293,'*390': 8.780515670776367,'*391': 8.780515670776367,'*392': 8.780515670776367,'*393': 8.780515670776367,'*394': 8.780515670776367,'*395': 8.780497550964355,'*396': 8.780497550964355,'*397': 8.780487060546875,'*398': 8.780458450317383,'*399': 8.780458450317383,'*4': 8.785660743713379,'*40': 8.785120964050293,'*400': 8.780458450317383,'*401': 8.780458450317383,'*402': 8.780458450317383,'*403': 8.780440330505371,'*404': 8.780451774597168,'*405': 8.780451774597168,'*406': 8.780499458312988,'*407': 8.780499458312988,'*408': 8.780499458312988,'*409': 8.780499458312988,'*41': 8.785120964050293,'*410': 8.780499458312988,'*411': 8.780499458312988,'*412': 8.780499458312988,'*413': 8.78039836883545,'*414': 8.78039836883545,'*415': 8.780316352844238,'*416': 8.780316352844238,'*417': 8.780316352844238,'*418': 8.780316352844238,'*419': 8.780316352844238,'*42': 8.785016059875488,'*420': 8.780316352844238,'*421': 8.780316352844238,'*422': 8.780316352844238,'*423': 8.780316352844238,'*424': 8.780316352844238,'*425': 8.780316352844238,'*426': 8.780336380004883,'*427': 8.780320167541504,'*428': 8.78018569946289,'*429': 8.78018569946289,'*43': 8.785016059875488,'*430': 8.78018569946289,'*431': 8.780202865600586,'*432': 8.780202865600586,'*433': 8.780187606811523,'*434': 8.780187606811523,'*435': 8.780187606811523,'*436': 8.780187606811523,'*437': 8.780187606811523,'*438': 8.780187606811523,'*439': 8.780187606811523,'*44': 8.785016059875488,'*440': 8.780158042907715,'*441': 8.780158042907715,'*442': 8.780158042907715,'*443': 8.780158042907715,'*444': 8.780149459838867,'*445': 8.780149459838867,'*446': 8.780149459838867,'*447': 8.780149459838867,'*448': 8.78013801574707,'*449': 8.78013801574707,'*45': 8.785016059875488,'*450': 8.78013801574707,'*451': 8.780152320861816,'*452': 8.780152320861816,'*453': 8.780152320861816,'*454': 8.780113220214844,'*455': 8.780120849609375,'*456': 8.780107498168945,'*457': 8.780107498168945,'*458': 8.780107498168945,'*459': 8.780070304870605,'*46': 8.785016059875488,'*460': 8.780070304870605,'*461': 8.780322074890137,'*462': 8.780322074890137,'*463': 8.780322074890137,'*464': 8.780322074890137,'*465': 8.780322074890137,'*466': 8.780322074890137,'*467': 8.780322074890137,'*468': 8.780322074890137,'*469': 8.780322074890137,'*47': 8.784873962402344,'*470': 8.780322074890137,'*471': 8.780322074890137,'*472': 8.780322074890137,'*473': 8.780322074890137,'*474': 8.780322074890137,'*475': 8.780322074890137,'*476': 8.780322074890137,'*477': 8.780200004577637,'*478': 8.780244827270508,'*479': 8.780193328857422,'*48': 8.784850120544434,'*480': 8.780193328857422,'*481': 8.780193328857422,'*482': 8.780193328857422,'*483': 8.780176162719727,'*484': 8.780176162719727,'*485': 8.780128479003906,'*486': 8.780128479003906,'*487': 8.780128479003906,'*488': 8.780128479003906,'*489': 8.78003978729248,'*49': 8.78482437133789,'*490': 8.78003978729248,'*491': 8.78003978729248,'*492': 8.78003978729248,'*493': 8.78003978729248,'*494': 8.779926300048828,'*495': 8.779926300048828,'*496': 8.77992057800293,'*497': 8.77992057800293,'*498': 8.77992057800293,'*499': 8.779901504516602,'*5': 8.785660743713379,'*50': 8.78482437133789,'*500': 8.779901504516602,'*501': 8.779901504516602,'*502': 8.779901504516602,'*503': 8.779901504516602,'*504': 8.779887199401855,'*505': 8.779887199401855,'*506': 8.779887199401855,'*507': 8.779887199401855,'*51': 8.78482437133789,'*52': 8.78482437133789,'*53': 8.78482437133789,'*54': 8.78482437133789,'*55': 8.78482437133789,'*56': 8.78482437133789,'*57': 8.78482437133789,'*58': 8.784768104553223,'*59': 8.784768104553223,'*6': 8.785565376281738,'*60': 8.784768104553223,'*61': 8.784768104553223,'*62': 8.784777641296387,'*63': 8.784777641296387,'*64': 8.784777641296387,'*65': 8.784693717956543,'*66': 8.784693717956543,'*67': 8.784693717956543,'*68': 8.784693717956543,'*69': 8.784693717956543,'*7': 8.785565376281738,'*70': 8.784612655639648,'*71': 8.784605026245117,'*72': 8.784605026245117,'*73': 8.784469604492188,'*74': 8.784469604492188,'*75': 8.784469604492188,'*76': 8.784469604492188,'*77': 8.784489631652832,'*78': 8.784473419189453,'*79': 8.784473419189453,'*8': 8.785565376281738,'*80': 8.784473419189453,'*81': 8.784463882446289,'*82': 8.784463882446289,'*83': 8.784463882446289,'*84': 8.784463882446289,'*85': 8.784463882446289,'*86': 8.784398078918457,'*87': 8.78575611114502,'*88': 8.78575611114502,'*89': 8.78575611114502,'*9': 8.785565376281738,'*90': 8.785689353942871,'*91': 8.785689353942871,'*92': 8.785689353942871,'*93': 8.785689353942871,'*94': 8.785684585571289,'*95': 8.785684585571289,'*96': 8.785684585571289,'*97': 8.785658836364746,'*98': 8.785658836364746,'*99': 8.786060333251953}


            # test_mosaic_cube_pcwdT.exp_pa_dict
            exp_pa_dict = {'*0': 67.42774200439453,'*1': 67.42774200439453,'*10': 67.42865753173828,'*100': 67.46749114990234,'*101': 67.46749114990234,'*102': 67.46749114990234,'*103': 67.46749114990234,'*104': 67.46749114990234,'*105': 67.46749114990234,'*106': 67.46749114990234,'*107': 67.4686050415039,'*108': 67.4686050415039,'*109': 67.4686050415039,'*11': 67.42865753173828,'*110': 67.4686050415039,'*111': 67.46869659423828,'*112': 67.40909576416016,'*113': 67.40939331054688,'*114': 67.40939331054688,'*115': 67.40939331054688,'*116': 67.40939331054688,'*117': 67.40909576416016,'*118': 67.40943145751953,'*119': 67.40943145751953,'*12': 67.42815399169922,'*120': 67.40927124023438,'*121': 67.40985870361328,'*122': 67.40985870361328,'*123': 67.40985870361328,'*124': 67.40985870361328,'*125': 67.40985870361328,'*126': 67.40985870361328,'*127': 67.40949249267578,'*128': 67.40949249267578,'*129': 67.40902709960938,'*13': 67.42815399169922,'*130': 67.40902709960938,'*131': 67.40931701660156,'*132': 67.40931701660156,'*133': 67.40931701660156,'*134': 67.40931701660156,'*135': 67.40931701660156,'*136': 67.40931701660156,'*137': 67.40974426269531,'*138': 67.40974426269531,'*139': 67.40974426269531,'*14': 67.42815399169922,'*140': 67.40974426269531,'*141': 67.40997314453125,'*142': 67.40997314453125,'*143': 67.40997314453125,'*144': 67.40997314453125,'*145': 67.40997314453125,'*146': 67.40936279296875,'*147': 67.40936279296875,'*148': 67.40936279296875,'*149': 67.4095687866211,'*15': 67.42796325683594,'*150': 67.4095687866211,'*151': 67.4095687866211,'*152': 67.40923309326172,'*153': 67.40923309326172,'*154': 67.40921020507812,'*155': 67.40921020507812,'*156': 67.40921020507812,'*157': 67.40921020507812,'*158': 67.40969848632812,'*159': 67.40969848632812,'*16': 67.42796325683594,'*160': 67.40969848632812,'*161': 67.40969848632812,'*162': 67.40969848632812,'*163': 67.40887451171875,'*164': 67.40887451171875,'*165': 67.40887451171875,'*166': 67.40887451171875,'*167': 67.40887451171875,'*168': 67.40816497802734,'*169': 67.40816497802734,'*17': 67.42796325683594,'*170': 67.40816497802734,'*171': 67.40816497802734,'*172': 67.40882873535156,'*173': 67.40799713134766,'*174': 67.40827941894531,'*175': 67.40790557861328,'*176': 67.4001235961914,'*177': 67.39026641845703,'*178': 67.38996124267578,'*179': 67.38996124267578,'*18': 67.42818450927734,'*180': 67.38996124267578,'*181': 67.3896255493164,'*182': 67.38994598388672,'*183': 67.38994598388672,'*184': 67.38994598388672,'*185': 67.38994598388672,'*186': 67.38994598388672,'*187': 67.38994598388672,'*188': 67.38994598388672,'*189': 67.38983917236328,'*19': 67.42818450927734,'*190': 67.38983917236328,'*191': 67.39041137695312,'*192': 67.39041137695312,'*193': 67.39024353027344,'*194': 67.39043426513672,'*195': 67.39043426513672,'*196': 67.39043426513672,'*197': 67.38996124267578,'*198': 67.38996124267578,'*199': 67.38996124267578,'*2': 67.42774200439453,'*20': 67.42818450927734,'*200': 67.38945770263672,'*201': 67.38945770263672,'*202': 67.38945770263672,'*203': 67.38945770263672,'*204': 67.45548248291016,'*205': 67.45548248291016,'*206': 67.47216033935547,'*207': 67.47216033935547,'*208': 67.47383117675781,'*209': 67.47413635253906,'*21': 67.42818450927734,'*210': 67.47413635253906,'*211': 67.47418212890625,'*212': 67.47332763671875,'*213': 67.4735107421875,'*214': 67.4735107421875,'*215': 67.47306823730469,'*216': 67.47306823730469,'*217': 67.47306823730469,'*218': 67.47306823730469,'*219': 67.47328186035156,'*22': 67.4267807006836,'*220': 67.47045135498047,'*221': 67.43890380859375,'*222': 67.43890380859375,'*223': 67.43890380859375,'*224': 67.43890380859375,'*225': 67.4367446899414,'*226': 67.4367446899414,'*227': 67.4367446899414,'*228': 67.4367446899414,'*229': 67.4367446899414,'*23': 67.4267807006836,'*230': 67.4367446899414,'*231': 67.43634796142578,'*232': 67.43634796142578,'*233': 67.43634796142578,'*234': 67.43634796142578,'*235': 67.4362564086914,'*236': 67.4362564086914,'*237': 67.43624877929688,'*238': 67.43585968017578,'*239': 67.43585968017578,'*24': 67.4267807006836,'*240': 67.43636322021484,'*241': 67.43636322021484,'*242': 67.43636322021484,'*243': 67.43710327148438,'*244': 67.43710327148438,'*245': 67.43710327148438,'*246': 67.43762969970703,'*247': 67.43688201904297,'*248': 67.43688201904297,'*249': 67.43688201904297,'*25': 67.4267807006836,'*250': 67.43717956542969,'*251': 67.43717956542969,'*252': 67.43766021728516,'*253': 67.43732452392578,'*254': 67.43732452392578,'*255': 67.43732452392578,'*256': 67.43777465820312,'*257': 67.43777465820312,'*258': 67.43777465820312,'*259': 67.4376449584961,'*26': 67.4267807006836,'*260': 67.44227600097656,'*261': 67.44729614257812,'*262': 67.44729614257812,'*263': 67.4476547241211,'*264': 67.4476547241211,'*265': 67.4476547241211,'*266': 67.4476547241211,'*267': 67.4476547241211,'*268': 67.4476547241211,'*269': 67.44754791259766,'*27': 67.42683410644531,'*270': 67.44754791259766,'*271': 67.44758605957031,'*272': 67.44739532470703,'*273': 67.44739532470703,'*274': 67.44739532470703,'*275': 67.44739532470703,'*276': 67.44766235351562,'*277': 67.44766235351562,'*278': 67.44766235351562,'*279': 67.44766235351562,'*28': 67.42683410644531,'*280': 67.4478759765625,'*281': 67.4478759765625,'*282': 67.44815826416016,'*283': 67.4476547241211,'*284': 67.4476547241211,'*285': 67.4476547241211,'*286': 67.44598388671875,'*287': 67.44554901123047,'*288': 67.44554901123047,'*289': 67.44438934326172,'*29': 67.42701721191406,'*290': 67.47785186767578,'*291': 67.47785186767578,'*292': 67.4777603149414,'*293': 67.4777603149414,'*294': 67.4777603149414,'*295': 67.47737121582031,'*296': 67.47737121582031,'*297': 67.47737121582031,'*298': 67.47763061523438,'*299': 67.47763061523438,'*3': 67.42806243896484,'*30': 67.42701721191406,'*300': 67.47763061523438,'*301': 67.47763061523438,'*302': 67.47763061523438,'*303': 67.47791290283203,'*304': 67.47791290283203,'*305': 67.47791290283203,'*306': 67.47791290283203,'*307': 67.47791290283203,'*308': 67.47820281982422,'*309': 67.47786712646484,'*31': 67.42622375488281,'*310': 67.47786712646484,'*311': 67.47786712646484,'*312': 67.47786712646484,'*313': 67.47786712646484,'*314': 67.47786712646484,'*315': 67.47734832763672,'*316': 67.47640228271484,'*317': 67.47640228271484,'*318': 67.47581481933594,'*319': 67.47581481933594,'*32': 67.42555236816406,'*320': 67.47581481933594,'*321': 67.47581481933594,'*322': 67.47581481933594,'*323': 67.47581481933594,'*324': 67.47581481933594,'*325': 67.47581481933594,'*326': 67.47581481933594,'*327': 67.47581481933594,'*328': 67.47581481933594,'*329': 67.47581481933594,'*33': 67.42555236816406,'*330': 67.47581481933594,'*331': 67.47685241699219,'*332': 67.47685241699219,'*333': 67.47685241699219,'*334': 67.47685241699219,'*335': 67.47685241699219,'*336': 67.47913360595703,'*337': 67.47913360595703,'*338': 67.4788818359375,'*339': 67.4780502319336,'*34': 67.42555236816406,'*340': 67.4780502319336,'*341': 67.4780502319336,'*342': 67.4780502319336,'*343': 67.4780502319336,'*344': 67.48170471191406,'*345': 67.4820556640625,'*346': 67.48235321044922,'*347': 67.4825439453125,'*348': 67.48287200927734,'*349': 67.48287200927734,'*35': 67.42527770996094,'*350': 67.48243713378906,'*351': 67.48243713378906,'*352': 67.48243713378906,'*353': 67.48243713378906,'*354': 67.48257446289062,'*355': 67.48257446289062,'*356': 67.48257446289062,'*357': 67.48297119140625,'*358': 67.48297119140625,'*359': 67.48279571533203,'*36': 67.42559814453125,'*360': 67.48289489746094,'*361': 67.48348999023438,'*362': 67.48348999023438,'*363': 67.48347473144531,'*364': 67.48347473144531,'*365': 67.48347473144531,'*366': 67.48347473144531,'*367': 67.48367309570312,'*368': 67.48367309570312,'*369': 67.52080535888672,'*37': 67.42813110351562,'*370': 67.52080535888672,'*371': 67.52066802978516,'*372': 67.52066802978516,'*373': 67.53384399414062,'*374': 67.53177642822266,'*375': 67.53177642822266,'*376': 67.53177642822266,'*377': 67.53177642822266,'*378': 67.53146362304688,'*379': 67.53092193603516,'*38': 67.42813110351562,'*380': 67.5309829711914,'*381': 67.52033996582031,'*382': 67.5829086303711,'*383': 67.5829086303711,'*384': 67.58301544189453,'*385': 67.58301544189453,'*386': 67.5822525024414,'*387': 67.58255004882812,'*388': 67.58284759521484,'*389': 67.58366394042969,'*39': 67.42864990234375,'*390': 67.58366394042969,'*391': 67.58366394042969,'*392': 67.58366394042969,'*393': 67.58366394042969,'*394': 67.58366394042969,'*395': 67.58395385742188,'*396': 67.58395385742188,'*397': 67.58430480957031,'*398': 67.58397674560547,'*399': 67.58397674560547,'*4': 67.42806243896484,'*40': 67.42864990234375,'*400': 67.58397674560547,'*401': 67.58397674560547,'*402': 67.58397674560547,'*403': 67.58432006835938,'*404': 67.58432006835938,'*405': 67.58432006835938,'*406': 67.5838623046875,'*407': 67.5838623046875,'*408': 67.5838623046875,'*409': 67.5838623046875,'*41': 67.42864990234375,'*410': 67.5838623046875,'*411': 67.5838623046875,'*412': 67.5838623046875,'*413': 67.5837173461914,'*414': 67.5837173461914,'*415': 67.58311462402344,'*416': 67.58311462402344,'*417': 67.58311462402344,'*418': 67.58311462402344,'*419': 67.58311462402344,'*42': 67.42789459228516,'*420': 67.58311462402344,'*421': 67.58311462402344,'*422': 67.58311462402344,'*423': 67.58311462402344,'*424': 67.58311462402344,'*425': 67.58311462402344,'*426': 67.58353424072266,'*427': 67.5837173461914,'*428': 67.58344268798828,'*429': 67.58344268798828,'*43': 67.42789459228516,'*430': 67.58344268798828,'*431': 67.58275604248047,'*432': 67.58275604248047,'*433': 67.58304595947266,'*434': 67.58304595947266,'*435': 67.58304595947266,'*436': 67.58304595947266,'*437': 67.58304595947266,'*438': 67.58304595947266,'*439': 67.58304595947266,'*44': 67.42789459228516,'*440': 67.5838851928711,'*441': 67.5838851928711,'*442': 67.5838851928711,'*443': 67.5838851928711,'*444': 67.5840835571289,'*445': 67.5840835571289,'*446': 67.5840835571289,'*447': 67.5840835571289,'*448': 67.58444213867188,'*449': 67.58444213867188,'*45': 67.42789459228516,'*450': 67.58444213867188,'*451': 67.58393096923828,'*452': 67.58393096923828,'*453': 67.58393096923828,'*454': 67.58477783203125,'*455': 67.5848159790039,'*456': 67.58538055419922,'*457': 67.58538055419922,'*458': 67.58538055419922,'*459': 67.58595275878906,'*46': 67.42789459228516,'*460': 67.58595275878906,'*461': 67.5853042602539,'*462': 67.5853042602539,'*463': 67.5853042602539,'*464': 67.5853042602539,'*465': 67.5853042602539,'*466': 67.5853042602539,'*467': 67.5853042602539,'*468': 67.5853042602539,'*469': 67.5853042602539,'*47': 67.42756652832031,'*470': 67.5853042602539,'*471': 67.5853042602539,'*472': 67.5853042602539,'*473': 67.5853042602539,'*474': 67.5853042602539,'*475': 67.5853042602539,'*476': 67.5853042602539,'*477': 67.58475494384766,'*478': 67.58609771728516,'*479': 67.58577728271484,'*48': 67.4280776977539,'*480': 67.58577728271484,'*481': 67.58577728271484,'*482': 67.58577728271484,'*483': 67.5859603881836,'*484': 67.5859603881836,'*485': 67.58699798583984,'*486': 67.58699798583984,'*487': 67.58699798583984,'*488': 67.58699798583984,'*489': 67.5868911743164,'*49': 67.42770385742188,'*490': 67.5868911743164,'*491': 67.5868911743164,'*492': 67.5868911743164,'*493': 67.5868911743164,'*494': 67.58709716796875,'*495': 67.58709716796875,'*496': 67.58740997314453,'*497': 67.58740997314453,'*498': 67.58740997314453,'*499': 67.58714294433594,'*5': 67.42806243896484,'*50': 67.42770385742188,'*500': 67.58714294433594,'*501': 67.58714294433594,'*502': 67.58714294433594,'*503': 67.58714294433594,'*504': 67.58756256103516,'*505': 67.58756256103516,'*506': 67.58756256103516,'*507': 67.58756256103516,'*51': 67.42770385742188,'*52': 67.42770385742188,'*53': 67.42770385742188,'*54': 67.42770385742188,'*55': 67.42770385742188,'*56': 67.42770385742188,'*57': 67.42770385742188,'*58': 67.42686462402344,'*59': 67.42686462402344,'*6': 67.42777252197266,'*60': 67.42686462402344,'*61': 67.42686462402344,'*62': 67.42632293701172,'*63': 67.42632293701172,'*64': 67.42632293701172,'*65': 67.42557525634766,'*66': 67.42557525634766,'*67': 67.42557525634766,'*68': 67.42557525634766,'*69': 67.42557525634766,'*7': 67.42777252197266,'*70': 67.42498016357422,'*71': 67.42517852783203,'*72': 67.42517852783203,'*73': 67.42532348632812,'*74': 67.42532348632812,'*75': 67.42532348632812,'*76': 67.42532348632812,'*77': 67.42509460449219,'*78': 67.42562866210938,'*79': 67.42562866210938,'*8': 67.42777252197266,'*80': 67.42562866210938,'*81': 67.4272689819336,'*82': 67.4272689819336,'*83': 67.4272689819336,'*84': 67.4272689819336,'*85': 67.4272689819336,'*86': 67.4271011352539,'*87': 67.46137237548828,'*88': 67.46137237548828,'*89': 67.46137237548828,'*9': 67.42777252197266,'*90': 67.46118927001953,'*91': 67.46118927001953,'*92': 67.46118927001953,'*93': 67.46118927001953,'*94': 67.46138000488281,'*95': 67.46138000488281,'*96': 67.46138000488281,'*97': 67.46085357666016,'*98': 67.46085357666016,'*99': 67.46749114990234}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed = self.filter_report(report)

        img = shutil._basename(img)
        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, test_name, self.parallel)

        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        self.mom8_creator(img+'.image', range_list=[0.15, 1.2])
        self.mom8_creator(img+'.residual', range_list=[0.15, 1.2])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')

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

# End of test_mosaic_cube_pcwdT
#-------------------------------------------------#
    # Test 9c
    @stats_dict(test_dict)
    def test_mosaic_cube_briggsbwtaper(self):
        ''' Mosaic cube imaging wth briggsbwtaper - field SMIDGE_NWCloud, spw 22 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
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
            perchanweightdensity=True, gridder='mosaic', \
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
                gridder='mosaic', mosweight=True, \
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
                gridder='mosaic', mosweight=True, \
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
                gridder='mosaic', mosweight=True, \
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
        exp_im_stats = {'com_bmaj': [False, 8.799594251773845],
            'com_bmin': [False, 6.099815949869395],
            'com_pa': [False, 67.46480038838807],
            'npts': [True, 5925312],
            'npts_unmasked': [False, 3352579.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.2676681280136108],
            'max_val_pos': [True, [45, 38, 0, 252]],
            'min_val': [False, -0.42875850200653076],
            'min_val_pos': [True, [48, 32, 0, 272]],
            'im_rms': [False, 0.08757929422825417],
            'rms_per_chan': [False, [0.08454448908911984, 0.08080263648768907, 0.08214279393909761, 0.08642318640888869, 0.08559389777148955, 0.0826678120056367, 0.08075963993970686, 0.0886406129800311, 0.08511467368448941, 0.08710165235214407, 0.07455555693559782, 0.08679355150897247, 0.07672297599313567, 0.0833949173922621, 0.08342601525186094, 0.09507914872509343, 0.07820515271768179, 0.08229263254437022, 0.07865683185812504, 0.09174598374745607, 0.08281488919913275, 0.08850180977098479, 0.08719597000576199, 0.0798279186018233, 0.08900869427637728, 0.08676195660808413, 0.09228273094904177, 0.07952046370618707, 0.08360526100714491, 0.08952729021335463, 0.0768737383586313, 0.08546788200832263, 0.07222108446959373, 0.08123054705441912, 0.08472955843645151, 0.08868644789356141, 0.07975308601037645, 0.0810112252815736, 0.09005101054216544, 0.09064277514742965, 0.08432143850019164, 0.08441895928657504, 0.08225354794825805, 0.08659678111034046, 0.0841907859596131, 0.09091523885740316, 0.07700200845663062, 0.09263466958661558, 0.08995362713140286, 0.08804900898004732, 0.08334611304149844, 0.079945163521021, 0.08361064805760578, 0.08139605046045667, 0.09639965651079935, 0.07344614343865323, 0.08410777206856682, 0.0918076547949407, 0.08404037292008604, 0.07488390254410975, 0.08726593961357626, 0.07592264028606217, 0.07869036731155865, 0.07712995508538403, 0.088946727901769, 0.0788589961405431, 0.08433288075890634, 0.08015469391373116, 0.08083529533844117, 0.08232100440295041, 0.09436694422580345, 0.08143081001362033, 0.07951910549408184, 0.08566876607196534, 0.08200417174011049, 0.07546036679889183, 0.08576592178824521, 0.08390201783746216, 0.08166300713054572, 0.08294476974708974, 0.08175034453695106, 0.09309839710925705, 0.08956282324776302, 0.08863806729732439, 0.08768177555631491, 0.08730355096033299, 0.09920658448951851, 0.08893986473865599, 0.08525145324493652, 0.0910542109223899, 0.08294171360146108, 0.08279749095443037, 0.08973901618353275, 0.0914915048775252, 0.09736128797765987, 0.08744920343331888, 0.08376733832441756, 0.08826914937611756, 0.08310085489638293, 0.08040603596032231, 0.0827948727015347, 0.08359258244124981, 0.07499660799933862, 0.0820795524519901, 0.08305418170311316, 0.08848395429721194, 0.09254540361584482, 0.0834250850720901, 0.09154432446675109, 0.08640487720343262, 0.0957902478621614, 0.09355943319944192, 0.08179046125426691, 0.0832043530817441, 0.08513172983407978, 0.08096771376062059, 0.08524244719511807, 0.08572843807117855, 0.08498689152548969, 0.0874652896630942, 0.08374668608074906, 0.0845395998019829, 0.07921815161199633, 0.08994484169183976, 0.09777083346759391, 0.0941529115071451, 0.09045874195400665, 0.08425352444354596, 0.09084473386135752, 0.08948140388989684, 0.09252015455406842, 0.07885057346611088, 0.08507495740462206, 0.08924207416636859, 0.07631696775562914, 0.08815350286139056, 0.08914844758149096, 0.08521386911186729, 0.08768673766479529, 0.08468883620070247, 0.08244560709624789, 0.09358733735181854, 0.08794678247456619, 0.09132743936261446, 0.07839269531347082, 0.08255793725800023, 0.08443809350504065, 0.08955331247016357, 0.07448797173292698, 0.08285343106495359, 0.0753733515398988, 0.08945893324209181, 0.09224320541897765, 0.08882548684208846, 0.08657529808469516, 0.08077973458973389, 0.08287353607837047, 0.08308038941701508, 0.082797098512978, 0.08296035329644438, 0.07872480428865179, 0.08408744829552543, 0.08880896807459428, 0.0862690974471767, 0.08466232473038289, 0.08689289347694795, 0.08377787525221031, 0.08593553889147006, 0.08088236733837953, 0.09030298769329086, 0.08573942112746871, 0.09584468678628222, 0.09097002040190538, 0.08739387947505368, 0.08145790245342659, 0.09035904742640996, 0.09593563107517897, 0.07540521888048991, 0.08585748167972262, 0.08547553945642043, 0.08853935168419275, 0.08314534550351242, 0.07538158044700248, 0.08424487424148978, 0.09200544898572927, 0.09297870174866711, 0.08405201082243019, 0.08390768117419219, 0.08743294278196162, 0.08651000707887543, 0.08393807100599077, 0.08422835676910953, 0.08735397302755507, 0.08385091781592317, 0.0871260702735577, 0.08914487593626652, 0.08285721054292167, 0.08940838209689977, 0.08795480019506897, 0.09857963204977876, 0.090381777862029, 0.08630513827079805, 0.08074908272369553, 0.08833487032474728, 0.09339171375548354, 0.07915710363427544, 0.08722057086151329, 0.08445367854086823, 0.08817932886122357, 0.08991099021186998, 0.08958524620572378, 0.08958065804161094, 0.08598668010091035, 0.08244158781235865, 0.08803063482643884, 0.08511356109050992, 0.08403270223876957, 0.08094977549656493, 0.08964053650671833, 0.08400091878326434, 0.07870514009908415, 0.08068398524361658, 0.08996085176724132, 0.0906021854937618, 0.08709621234440584, 0.08164750820156866, 0.09597221098008116, 0.0884837020770368, 0.08040656710659008, 0.08882638693901852, 0.09001765250906213, 0.09589572158449743, 0.09200839055429873, 0.07843005418350324, 0.0899325220701504, 0.09043548373137923, 0.09509161758073698, 0.09045806423953247, 0.09470251868256363, 0.08309363749899958, 0.08294602851910425, 0.08542798116911876, 0.07967048535468631, 0.09249791046662145, 0.0852274428680305, 0.09388758188399719, 0.08732525084838393, 0.08735877068070008, 0.0940410796179019, 0.1000470438532653, 0.10722457094262983, 0.12578767243058564, 0.1436667744488629, 0.13808358370364404, 0.12065097931507894, 0.12435471666070047, 0.1072904741792376, 0.09040416625117817, 0.09100495304305524, 0.08324076431723941, 0.08942139276736703, 0.08803060876415686, 0.08754649851127198, 0.08811845372824391, 0.08266397025271952, 0.09328648083839139, 0.09508555582661059, 0.09340697992731502, 0.09716083332735184, 0.09370927423681424, 0.08493054052247156, 0.0840808386930263, 0.09225158111501996, 0.09720692354754973, 0.09102389876407037, 0.09190655670798578, 0.08222058377571392, 0.09548728919353866, 0.08285457770644784, 0.08621580960396537, 0.08102250602684535, 0.08528004311241288, 0.08971109756328768, 0.08495170225977651, 0.0751408536037243, 0.082568100994246, 0.0930147171842407, 0.07820490996722416, 0.07851334579592707, 0.08025647416536694, 0.08537280404694918, 0.08414771876833929, 0.09026605164985937, 0.0932069792922012, 0.0792168447753108, 0.09119788607776576, 0.08661683696956987, 0.09319673196141545, 0.0791714604288018, 0.09080747338359901, 0.08412818503997524, 0.08075005612938946, 0.0903801707941507, 0.08884241789906402, 0.08580036997699231, 0.08359031747632817, 0.0861522060004678, 0.08907133367469561, 0.09093398326326332, 0.08409855882413776, 0.08835671448802153, 0.08462314115551772, 0.08165947017924235, 0.08951802119060676, 0.09484951993173886, 0.09079195300054037, 0.08093162327815479, 0.08592684885631356, 0.08922288526960626, 0.09226875974759953, 0.09400908617559704, 0.09288180132440255, 0.09933750932216401, 0.0854001354312934, 0.0853565421252504, 0.07921483568167087, 0.08115639994221235, 0.09309796131865276, 0.07274957185490932, 0.08348594331936901, 0.09304174677243612, 0.09545682026043886, 0.09168006307285233, 0.09653391294935647, 0.07795155023726419, 0.09053736111647923, 0.08729776662190053, 0.08039301741487047, 0.08585192113925974, 0.08556364369233134, 0.0947909880946024, 0.08246086455041851, 0.08792522763195232, 0.08014983717531496, 0.0825721384853619, 0.08926401880520694, 0.09206493775579697, 0.08967875484420444, 0.09233535900079426, 0.09463409101792417, 0.09136512087468808, 0.07968756583775934, 0.09414676966118163, 0.07531204428349275, 0.0900887514881628, 0.082350929720856, 0.09116057896110935, 0.08482555349635658, 0.09093505541401749, 0.09298421351843844, 0.08954880952216236, 0.09304818858251004, 0.09082859863387693, 0.08872671972251833, 0.09001536037810437, 0.08841056982585542, 0.09187630751484305, 0.08645673093346762, 0.09077234000341755, 0.0902324335743101, 0.09348825660114995, 0.08809817603438537, 0.08291850818331703, 0.07852343515479626, 0.08775670723907779, 0.07902160982413586, 0.09597807127521123, 0.07948026146106765, 0.08408574454405071, 0.09146042908184006, 0.0903829727899694, 0.08427870951568463, 0.08398533770750455, 0.08441905701930158, 0.08182850905372979, 0.08894100255371722, 0.07759031614716395, 0.08972510934934153, 0.08264479878145543, 0.09284448327533745, 0.08179855322516799, 0.08685106852108863, 0.07989272358171175, 0.07797016408294571, 0.09311475127614742, 0.08695678074949989, 0.08951173230590632, 0.09463464230622358, 0.0839422544679647, 0.0929779941399182, 0.09367276852206438, 0.09088690550241038, 0.08581284101179938, 0.09027131901774704, 0.08639636853197506, 0.08823440350999151, 0.08981140679650465, 0.08565009649292102, 0.09193698820112021, 0.09260481419788986, 0.09219162694349527, 0.08452138517403869, 0.09327975194383413, 0.08561876574638069, 0.09031930480610689, 0.08797619409942842, 0.08904158410395703, 0.08677449954241682, 0.09193836862557542, 0.08932336591380748, 0.08888734947455079, 0.08408822267354282, 0.08360185223175279, 0.08726903776947037, 0.0920928944052632, 0.08445285983786797, 0.08297436231859105, 0.08177558616249073, 0.08458920929055264, 0.08053289681362301, 0.08381024158839587, 0.08394247455720552, 0.07949335809162283, 0.09044280027855939, 0.08200577315972736, 0.08049760190845602, 0.08720597796672637, 0.08590545317806317, 0.09018328426885515, 0.07849384920136365, 0.09387003917122505, 0.0820839948336462, 0.09124627988014393, 0.09056621113177035, 0.08084378476065637, 0.09119682201299913, 0.09154924297404958, 0.10102238519951766, 0.08969739162264642, 0.08232612120804222, 0.08397760055038085, 0.0877776654297572, 0.09208767585452475, 0.09472256425923893, 0.09034008482053962, 0.08517679713173801, 0.08349990233298345, 0.0920914456782035, 0.09036263603380006, 0.08296210480531124, 0.0882103922157589, 0.09354976887291509, 0.10061026666508853, 0.08536083083675522, 0.08399716153322566, 0.08711958320544728, 0.08387366165386896, 0.08737460383578502, 0.08930710422940288, 0.08712188701373794, 0.0860077767706973, 0.08937967502709068, 0.08869269499997347, 0.09355441125716389, 0.08626325555932915, 0.08447958942334301, 0.08859132123644435, 0.09710644965882617, 0.09103282000673704, 0.08494051335588168, 0.08757694104631766, 0.09066729195468025, 0.08396210340146089, 0.09616904666129485, 0.07252059042743, 0.08800188254352609, 0.08656231360176422, 0.09185176600070785, 0.09340475912384315, 0.07503944262782218, 0.09081409994500066, 0.09660314498027173, 0.09050787961551574, 0.10365284142111499, 0.08146582714790225, 0.09872441731201544, 0.07845170891629007, 0.0913113298273851, 0.08882789944712782, 0.09219245061695863, 0.09638213812465103, 0.0888878320253745, 0.08984413821767516, 0.09422194635651554, 0.0838001137555367, 0.09670389499597991, 0.0831012788573468, 0.08782080309377456]],
            'im_sum': [False, 148.66257001338366],
            'regn_sum': [False, 75.69317591190338],
            'npts_real': [True, 5925312],
            'rms_per_field': [False, [0.08868748848086015, 0.08749214376728211, 0.08899292220049493, 0.08696635882457056, 0.08861907837656693, 0.0871993844302931, 0.0875503409774584]],
            'profile': [False, 1.2529497439026103],
            'fit': [False, [1.278427614020036, 8.99650305866635, 7.750912802084286]],
            'fit_loc_chan': [True, 252],
            'fit_loc_freq': [1e-10, 220.31420623259388],
            'fit_pix': [False, [44.6056573431561, 38.21224371975814]]}

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
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'mask_pix': [False, 4105],
            'mask_regns': [True, 1],
            'npts_real': [True, 5925312]}


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube_briggsbwtaper)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3352579.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, 0.20000091195106506],
            'im_rms': [False, 0.6139604370757835],
            'npts_0.2': [False, [6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600]],
            'npts_0.5': [False, [3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3598, 3598, 3598, 3598, 3598, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599]],
            'npts_real': [True, 5925312],
            'fit': [False, [1.0797565979383854, 69.3294373980139, 69.27416917859814]],
            'fit_loc_chan': [True, 254],
            'fit_loc_freq': [1e-10, 220.31469458079383],
            'fit_pix': [False, [54.058064547261154, 53.98492340521667]]}


        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_briggsbwtaper)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [54, 54, 0, 0]],
            'min_val': [False, -0.1642415076494217],
            'min_val_pos': [True, [63, 52, 0, 381]],
            'im_rms': [False, 0.06138798977040539],
            'im_sum': [False, 27.14456590548383],
            'npts_real': [True, 5925312],
            'fit_0': [False, [1.0996293058564999, 7.80716186452208, 5.286416613043531]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 220.2529185335],
            'fit_pix_0': [False, [53.991038646146926, 54.00525608246625]],
            'fit_1': [False, [1.0996235203307996, 7.808582641190577, 5.2866183210543465]],
            'fit_loc_chan_1': [True, 254],
            'fit_loc_freq_1': [1e-10, 220.31469458079383],
            'fit_pix_1': [False, [53.990712163511134, 54.0047420190993]],
            'fit_2': [False, [1.0995721334697781, 7.799558939086169, 5.292710802755232]],
            'fit_loc_chan_2': [True, 507],
            'fit_loc_freq_2': [1e-10, 220.37647062808767],
            'fit_pix_2': [False, [53.99129632179105, 54.004158696033315]]}


        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 5925312],
            'npts_unmasked': [False, 3352579.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.4512392282485962],
            'max_val_pos': [True, [48, 60, 0, 249]],
            'min_val': [False, -0.42875850200653076],
            'min_val_pos': [True, [48, 32, 0, 272]],
            'im_rms': [False, 0.08723675828952442],
            'im_sum': [False, -139.17177731977606],
            'regn_sum': [False, 20.808523556333967],
            'npts_real': [True, 5925312]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_briggsbwtaper)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.5603349208831787],
            'max_val_pos': [True, [45, 39, 0, 251]],
            'min_val': [False, -0.08046812564134598],
            'min_val_pos': [True, [46, 57, 0, 255]],
            'im_rms': [False, 0.0005318594401396416],
            'im_sum': [False, 5.72643318772316],
            'regn_sum': [False, 5.72643318772316],
            'mask_non0': [True, 0],
            'npts_real': [True, 5925312]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 508],
            'npts_unmasked': [True, 508.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 121.04640197753906],
            'max_val_pos': [True, [0, 0, 0, 336]],
            'min_val': [False, 120.85958862304688],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 120.94878783552005],
            'npts_real': [True, 508]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube_briggsbwtaper.exp_wt_stats
        exp_wt_stats = {'npts': [True, 5925312],
            'npts_unmasked': [True, 5925312.0],
            'freq_bin': [1e-10, 244174.09997558594],
            'start': [True, 220252700000.0],
            'end': [True, 220376500000.0],
            'start_delta': [False, 220252700000.0],
            'end_delta': [False, 220376500000.0],
            'nchan': [True, 508],
            'max_val': [False, 0.39721715450286865],
            'max_val_pos': [True, [54, 54, 0, 381]],
            'min_val': [False, 7.157817162806168e-05],
            'im_rms': [False, 0.1420631138560751],
            'im_sum': [False, 511866.2260002265],
            'npts_0.2': [False, [6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6601, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6598, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600]],
            'npts_0.5': [False, [3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3598, 3598, 3598, 3598, 3598, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3597, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3594, 3594, 3594, 3594, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3595, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3596, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599, 3599]],
            'npts_real': [True, 5925312]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            # test_mosaic_cube_briggsbwtaper.exp_bmin_dict
            exp_bmin_dict = {'*0': 6.092315673828125,'*1': 6.092315673828125,'*10': 6.092411994934082,'*100': 6.092495441436768,'*101': 6.092495441436768,'*102': 6.092495441436768,'*103': 6.092495441436768,'*104': 6.092495441436768,'*105': 6.092495441436768,'*106': 6.092495441436768,'*107': 6.0923991203308105,'*108': 6.0923991203308105,'*109': 6.0923991203308105,'*11': 6.092411518096924,'*110': 6.0923991203308105,'*111': 6.092378616333008,'*112': 6.090467929840088,'*113': 6.0904412269592285,'*114': 6.0904412269592285,'*115': 6.0904412269592285,'*116': 6.0904412269592285,'*117': 6.09043025970459,'*118': 6.090396404266357,'*119': 6.090396404266357,'*12': 6.092401504516602,'*120': 6.0902886390686035,'*121': 6.090263843536377,'*122': 6.090263843536377,'*123': 6.090263843536377,'*124': 6.090263843536377,'*125': 6.090263843536377,'*126': 6.090263843536377,'*127': 6.090255260467529,'*128': 6.090255260467529,'*129': 6.090246677398682,'*13': 6.092401504516602,'*130': 6.090246677398682,'*131': 6.0902204513549805,'*132': 6.0902204513549805,'*133': 6.0902204513549805,'*134': 6.0902204513549805,'*135': 6.0902204513549805,'*136': 6.0902204513549805,'*137': 6.0901007652282715,'*138': 6.0901007652282715,'*139': 6.0901007652282715,'*14': 6.092401504516602,'*140': 6.0901007652282715,'*141': 6.090078830718994,'*142': 6.090078830718994,'*143': 6.090078830718994,'*144': 6.090078830718994,'*145': 6.090078830718994,'*146': 6.090074062347412,'*147': 6.090073585510254,'*148': 6.090073585510254,'*149': 6.090056896209717,'*15': 6.092396259307861,'*150': 6.090056896209717,'*151': 6.090056896209717,'*152': 6.090048789978027,'*153': 6.090048789978027,'*154': 6.090037822723389,'*155': 6.090037822723389,'*156': 6.090038299560547,'*157': 6.090038299560547,'*158': 6.090004920959473,'*159': 6.090004920959473,'*16': 6.092396259307861,'*160': 6.090004920959473,'*161': 6.090004920959473,'*162': 6.090005397796631,'*163': 6.089959144592285,'*164': 6.089959144592285,'*165': 6.089959144592285,'*166': 6.089959144592285,'*167': 6.089959144592285,'*168': 6.08994722366333,'*169': 6.08994722366333,'*17': 6.092396259307861,'*170': 6.08994722366333,'*171': 6.08994722366333,'*172': 6.089921951293945,'*173': 6.089905261993408,'*174': 6.089879035949707,'*175': 6.089872360229492,'*176': 6.088998317718506,'*177': 6.089479446411133,'*178': 6.089467525482178,'*179': 6.089467525482178,'*18': 6.092380523681641,'*180': 6.089467525482178,'*181': 6.089411735534668,'*182': 6.089385509490967,'*183': 6.089385509490967,'*184': 6.089385509490967,'*185': 6.089385509490967,'*186': 6.089385509490967,'*187': 6.089385509490967,'*188': 6.089385509490967,'*189': 6.089386940002441,'*19': 6.092380523681641,'*190': 6.089386940002441,'*191': 6.089321136474609,'*192': 6.089321136474609,'*193': 6.0893235206604,'*194': 6.089304447174072,'*195': 6.089304447174072,'*196': 6.089304447174072,'*197': 6.089300155639648,'*198': 6.089300155639648,'*199': 6.089300155639648,'*2': 6.092315673828125,'*20': 6.092380523681641,'*200': 6.0893049240112305,'*201': 6.0893049240112305,'*202': 6.0893049240112305,'*203': 6.0893049240112305,'*204': 6.0912861824035645,'*205': 6.0912861824035645,'*206': 6.085134983062744,'*207': 6.085134983062744,'*208': 6.0855231285095215,'*209': 6.0855021476745605,'*21': 6.092380523681641,'*210': 6.0855021476745605,'*211': 6.0854878425598145,'*212': 6.085480690002441,'*213': 6.085435390472412,'*214': 6.085435390472412,'*215': 6.085343837738037,'*216': 6.085343837738037,'*217': 6.085343837738037,'*218': 6.085343837738037,'*219': 6.08527135848999,'*22': 6.092365741729736,'*220': 6.085229396820068,'*221': 6.08436393737793,'*222': 6.08436393737793,'*223': 6.08436393737793,'*224': 6.08436393737793,'*225': 6.0855536460876465,'*226': 6.0855536460876465,'*227': 6.0855536460876465,'*228': 6.0855536460876465,'*229': 6.0855536460876465,'*23': 6.092365741729736,'*230': 6.0855536460876465,'*231': 6.0855183601379395,'*232': 6.0855183601379395,'*233': 6.085517883300781,'*234': 6.085517883300781,'*235': 6.085488319396973,'*236': 6.085488319396973,'*237': 6.085447311401367,'*238': 6.0854620933532715,'*239': 6.0854620933532715,'*24': 6.092365741729736,'*240': 6.085376262664795,'*241': 6.085376262664795,'*242': 6.085376262664795,'*243': 6.085365295410156,'*244': 6.085365295410156,'*245': 6.085365295410156,'*246': 6.085367202758789,'*247': 6.0853166580200195,'*248': 6.0853166580200195,'*249': 6.0853166580200195,'*25': 6.092365741729736,'*250': 6.08528995513916,'*251': 6.08528995513916,'*252': 6.085261821746826,'*253': 6.085257530212402,'*254': 6.085257530212402,'*255': 6.085257530212402,'*256': 6.08521032333374,'*257': 6.08521032333374,'*258': 6.08521032333374,'*259': 6.085212230682373,'*26': 6.092365741729736,'*260': 6.086493492126465,'*261': 6.086802005767822,'*262': 6.086802005767822,'*263': 6.086759090423584,'*264': 6.086759090423584,'*265': 6.086759090423584,'*266': 6.086759090423584,'*267': 6.086759090423584,'*268': 6.086759090423584,'*269': 6.086760520935059,'*27': 6.092310905456543,'*270': 6.086760520935059,'*271': 6.086742877960205,'*272': 6.086738109588623,'*273': 6.086738109588623,'*274': 6.086738109588623,'*275': 6.086738109588623,'*276': 6.086722373962402,'*277': 6.086722373962402,'*278': 6.086722373962402,'*279': 6.086722373962402,'*28': 6.092310905456543,'*280': 6.086705207824707,'*281': 6.086705207824707,'*282': 6.086678981781006,'*283': 6.086670398712158,'*284': 6.086670398712158,'*285': 6.086670398712158,'*286': 6.086511135101318,'*287': 6.086503982543945,'*288': 6.086503982543945,'*289': 6.086452007293701,'*29': 6.09229850769043,'*290': 6.0862531661987305,'*291': 6.0862531661987305,'*292': 6.086263656616211,'*293': 6.086263656616211,'*294': 6.086263656616211,'*295': 6.08625602722168,'*296': 6.08625602722168,'*297': 6.08625602722168,'*298': 6.086233139038086,'*299': 6.086233139038086,'*3': 6.092287540435791,'*30': 6.09229850769043,'*300': 6.086233139038086,'*301': 6.086233139038086,'*302': 6.086233139038086,'*303': 6.086206912994385,'*304': 6.086206912994385,'*305': 6.086206912994385,'*306': 6.086206912994385,'*307': 6.086206912994385,'*308': 6.086188793182373,'*309': 6.086181163787842,'*31': 6.092285633087158,'*310': 6.086181163787842,'*311': 6.086181163787842,'*312': 6.086181163787842,'*313': 6.086181163787842,'*314': 6.086181163787842,'*315': 6.086123466491699,'*316': 6.086106777191162,'*317': 6.086106777191162,'*318': 6.086097717285156,'*319': 6.086097717285156,'*32': 6.092410564422607,'*320': 6.086097717285156,'*321': 6.086097717285156,'*322': 6.086097717285156,'*323': 6.086097717285156,'*324': 6.086097717285156,'*325': 6.086097717285156,'*326': 6.086097717285156,'*327': 6.086097717285156,'*328': 6.086097717285156,'*329': 6.086097717285156,'*33': 6.092410564422607,'*330': 6.086097717285156,'*331': 6.086109161376953,'*332': 6.086109161376953,'*333': 6.086109161376953,'*334': 6.086109161376953,'*335': 6.086109161376953,'*336': 6.086208343505859,'*337': 6.086208343505859,'*338': 6.086180686950684,'*339': 6.086175441741943,'*34': 6.092410564422607,'*340': 6.086174964904785,'*341': 6.086174964904785,'*342': 6.086174964904785,'*343': 6.086174964904785,'*344': 6.086239814758301,'*345': 6.086220741271973,'*346': 6.0861897468566895,'*347': 6.086155891418457,'*348': 6.086126804351807,'*349': 6.086126804351807,'*35': 6.092439651489258,'*350': 6.086119651794434,'*351': 6.086119651794434,'*352': 6.086119651794434,'*353': 6.086119651794434,'*354': 6.086115837097168,'*355': 6.086115837097168,'*356': 6.086115837097168,'*357': 6.086080074310303,'*358': 6.086080074310303,'*359': 6.08607816696167,'*36': 6.092409133911133,'*360': 6.086044788360596,'*361': 6.0859832763671875,'*362': 6.0859832763671875,'*363': 6.0859575271606445,'*364': 6.0859575271606445,'*365': 6.0859575271606445,'*366': 6.0859575271606445,'*367': 6.085938930511475,'*368': 6.085938930511475,'*369': 6.091782093048096,'*37': 6.09257698059082,'*370': 6.091782093048096,'*371': 6.091750621795654,'*372': 6.091750621795654,'*373': 6.091124057769775,'*374': 6.090796947479248,'*375': 6.090796947479248,'*376': 6.090796947479248,'*377': 6.090796947479248,'*378': 6.090790748596191,'*379': 6.090785980224609,'*38': 6.09257698059082,'*380': 6.090754985809326,'*381': 6.090539455413818,'*382': 6.091299057006836,'*383': 6.091299057006836,'*384': 6.091263294219971,'*385': 6.091263294219971,'*386': 6.091251373291016,'*387': 6.0912251472473145,'*388': 6.091195583343506,'*389': 6.091184616088867,'*39': 6.092559337615967,'*390': 6.091184616088867,'*391': 6.091184616088867,'*392': 6.091184616088867,'*393': 6.091184616088867,'*394': 6.091184616088867,'*395': 6.091185092926025,'*396': 6.091185092926025,'*397': 6.091165542602539,'*398': 6.091159343719482,'*399': 6.091159343719482,'*4': 6.092287540435791,'*40': 6.092559337615967,'*400': 6.091159343719482,'*401': 6.091159343719482,'*402': 6.091159343719482,'*403': 6.091114521026611,'*404': 6.091094970703125,'*405': 6.091094970703125,'*406': 6.091066837310791,'*407': 6.091066837310791,'*408': 6.091066837310791,'*409': 6.091066837310791,'*41': 6.092559337615967,'*410': 6.091066837310791,'*411': 6.091066837310791,'*412': 6.091066837310791,'*413': 6.0910420417785645,'*414': 6.0910420417785645,'*415': 6.091029644012451,'*416': 6.091029644012451,'*417': 6.091029644012451,'*418': 6.091029644012451,'*419': 6.091029644012451,'*42': 6.09254264831543,'*420': 6.091029644012451,'*421': 6.091029644012451,'*422': 6.091029644012451,'*423': 6.091030120849609,'*424': 6.091030120849609,'*425': 6.091030120849609,'*426': 6.090948104858398,'*427': 6.09092378616333,'*428': 6.090950012207031,'*429': 6.090950012207031,'*43': 6.09254264831543,'*430': 6.090950012207031,'*431': 6.090933322906494,'*432': 6.090933322906494,'*433': 6.090906143188477,'*434': 6.090906143188477,'*435': 6.090906143188477,'*436': 6.090906143188477,'*437': 6.090906143188477,'*438': 6.090906143188477,'*439': 6.090906143188477,'*44': 6.09254264831543,'*440': 6.0908355712890625,'*441': 6.0908355712890625,'*442': 6.0908355712890625,'*443': 6.0908355712890625,'*444': 6.090817928314209,'*445': 6.090817928314209,'*446': 6.090817928314209,'*447': 6.090817928314209,'*448': 6.090799331665039,'*449': 6.090799331665039,'*45': 6.09254264831543,'*450': 6.090799331665039,'*451': 6.090786933898926,'*452': 6.090786933898926,'*453': 6.090786933898926,'*454': 6.090739727020264,'*455': 6.090715408325195,'*456': 6.090665817260742,'*457': 6.090665817260742,'*458': 6.090665817260742,'*459': 6.090618133544922,'*46': 6.09254264831543,'*460': 6.090618133544922,'*461': 6.090680122375488,'*462': 6.0906805992126465,'*463': 6.090680122375488,'*464': 6.090680122375488,'*465': 6.090680122375488,'*466': 6.0906805992126465,'*467': 6.0906805992126465,'*468': 6.0906805992126465,'*469': 6.0906805992126465,'*47': 6.0924973487854,'*470': 6.0906805992126465,'*471': 6.0906805992126465,'*472': 6.0906805992126465,'*473': 6.0906805992126465,'*474': 6.0906805992126465,'*475': 6.0906805992126465,'*476': 6.0906805992126465,'*477': 6.090638160705566,'*478': 6.090646266937256,'*479': 6.090641021728516,'*48': 6.09246826171875,'*480': 6.090641021728516,'*481': 6.090641021728516,'*482': 6.090641021728516,'*483': 6.090616226196289,'*484': 6.090616226196289,'*485': 6.090596675872803,'*486': 6.090596675872803,'*487': 6.090596675872803,'*488': 6.090596675872803,'*489': 6.0905985832214355,'*49': 6.092462062835693,'*490': 6.0905985832214355,'*491': 6.0905985832214355,'*492': 6.0905985832214355,'*493': 6.0905985832214355,'*494': 6.090574264526367,'*495': 6.090574264526367,'*496': 6.090549945831299,'*497': 6.090549945831299,'*498': 6.090549945831299,'*499': 6.090546131134033,'*5': 6.092287540435791,'*50': 6.092462062835693,'*500': 6.090546131134033,'*501': 6.090546131134033,'*502': 6.090546131134033,'*503': 6.090545654296875,'*504': 6.090516567230225,'*505': 6.090516567230225,'*506': 6.090516567230225,'*507': 6.090516567230225,'*51': 6.092462062835693,'*52': 6.092462062835693,'*53': 6.092462062835693,'*54': 6.092462062835693,'*55': 6.092462062835693,'*56': 6.092462062835693,'*57': 6.092462062835693,'*58': 6.092446804046631,'*59': 6.092446804046631,'*6': 6.09244441986084,'*60': 6.092446804046631,'*61': 6.092446804046631,'*62': 6.092434883117676,'*63': 6.092434883117676,'*64': 6.092435359954834,'*65': 6.092422962188721,'*66': 6.092422962188721,'*67': 6.092422962188721,'*68': 6.092422962188721,'*69': 6.092422962188721,'*7': 6.09244441986084,'*70': 6.092411041259766,'*71': 6.0923943519592285,'*72': 6.0923943519592285,'*73': 6.092405796051025,'*74': 6.092405796051025,'*75': 6.092405796051025,'*76': 6.092405796051025,'*77': 6.09239387512207,'*78': 6.092380523681641,'*79': 6.092380523681641,'*8': 6.09244441986084,'*80': 6.092380523681641,'*81': 6.092370986938477,'*82': 6.092370986938477,'*83': 6.092370986938477,'*84': 6.092370986938477,'*85': 6.092370986938477,'*86': 6.092373847961426,'*87': 6.092708110809326,'*88': 6.092708110809326,'*89': 6.092708110809326,'*9': 6.09244441986084,'*90': 6.092706203460693,'*91': 6.092706203460693,'*92': 6.092706203460693,'*93': 6.092706203460693,'*94': 6.09269380569458,'*95': 6.09269380569458,'*96': 6.09269380569458,'*97': 6.092684268951416,'*98': 6.092684268951416,'*99': 6.092495441436768}

            # test_mosaic_cube_briggsbwtaper.exp_bmaj_dict
            exp_bmaj_dict = {'*0': 8.785123825073242,'*1': 8.785123825073242,'*10': 8.785015106201172,'*100': 8.785499572753906,'*101': 8.785499572753906,'*102': 8.785499572753906,'*103': 8.785499572753906,'*104': 8.785499572753906,'*105': 8.785499572753906,'*106': 8.785499572753906,'*107': 8.785449028015137,'*108': 8.785449028015137,'*109': 8.785449028015137,'*11': 8.785014152526855,'*110': 8.785449028015137,'*111': 8.785455703735352,'*112': 8.788925170898438,'*113': 8.788914680480957,'*114': 8.788914680480957,'*115': 8.788914680480957,'*116': 8.788914680480957,'*117': 8.788841247558594,'*118': 8.788835525512695,'*119': 8.788835525512695,'*12': 8.78502368927002,'*120': 8.788923263549805,'*121': 8.788887023925781,'*122': 8.788887023925781,'*123': 8.788887023925781,'*124': 8.788887023925781,'*125': 8.788887023925781,'*126': 8.788887023925781,'*127': 8.788899421691895,'*128': 8.788899421691895,'*129': 8.78883171081543,'*13': 8.78502368927002,'*130': 8.78883171081543,'*131': 8.788817405700684,'*132': 8.788817405700684,'*133': 8.788817405700684,'*134': 8.788817405700684,'*135': 8.788817405700684,'*136': 8.788817405700684,'*137': 8.788885116577148,'*138': 8.788885116577148,'*139': 8.788885116577148,'*14': 8.78502368927002,'*140': 8.788885116577148,'*141': 8.788875579833984,'*142': 8.788875579833984,'*143': 8.788875579833984,'*144': 8.788875579833984,'*145': 8.788875579833984,'*146': 8.788817405700684,'*147': 8.788816452026367,'*148': 8.788816452026367,'*149': 8.788808822631836,'*15': 8.784976959228516,'*150': 8.788808822631836,'*151': 8.788808822631836,'*152': 8.788753509521484,'*153': 8.788753509521484,'*154': 8.788761138916016,'*155': 8.788761138916016,'*156': 8.788761138916016,'*157': 8.788761138916016,'*158': 8.788742065429688,'*159': 8.788742065429688,'*16': 8.784976959228516,'*160': 8.788742065429688,'*161': 8.788742065429688,'*162': 8.788742065429688,'*163': 8.78860855102539,'*164': 8.78860855102539,'*165': 8.78860855102539,'*166': 8.78860855102539,'*167': 8.78860855102539,'*168': 8.788476943969727,'*169': 8.788476943969727,'*17': 8.784976959228516,'*170': 8.788476943969727,'*171': 8.788476943969727,'*172': 8.788444519042969,'*173': 8.788328170776367,'*174': 8.788322448730469,'*175': 8.788290023803711,'*176': 8.788912773132324,'*177': 8.788199424743652,'*178': 8.788125991821289,'*179': 8.788125991821289,'*18': 8.7849702835083,'*180': 8.788125991821289,'*181': 8.788178443908691,'*182': 8.788171768188477,'*183': 8.788171768188477,'*184': 8.788171768188477,'*185': 8.788171768188477,'*186': 8.788171768188477,'*187': 8.788171768188477,'*188': 8.788171768188477,'*189': 8.788064002990723,'*19': 8.7849702835083,'*190': 8.788064002990723,'*191': 8.788066864013672,'*192': 8.788066864013672,'*193': 8.788000106811523,'*194': 8.787995338439941,'*195': 8.787995338439941,'*196': 8.787995338439941,'*197': 8.787955284118652,'*198': 8.787955284118652,'*199': 8.787955284118652,'*2': 8.785123825073242,'*20': 8.7849702835083,'*200': 8.787901878356934,'*201': 8.787901878356934,'*202': 8.787901878356934,'*203': 8.787901878356934,'*204': 8.787531852722168,'*205': 8.787531852722168,'*206': 8.789189338684082,'*207': 8.789189338684082,'*208': 8.789972305297852,'*209': 8.789958000183105,'*21': 8.7849702835083,'*210': 8.789958000183105,'*211': 8.789963722229004,'*212': 8.789897918701172,'*213': 8.78991413116455,'*214': 8.78991413116455,'*215': 8.789904594421387,'*216': 8.789904594421387,'*217': 8.789904594421387,'*218': 8.789904594421387,'*219': 8.789904594421387,'*22': 8.784937858581543,'*220': 8.789737701416016,'*221': 8.790940284729004,'*222': 8.790940284729004,'*223': 8.79094123840332,'*224': 8.79094123840332,'*225': 8.788185119628906,'*226': 8.788185119628906,'*227': 8.788185119628906,'*228': 8.788185119628906,'*229': 8.788185119628906,'*23': 8.784937858581543,'*230': 8.788185119628906,'*231': 8.788126945495605,'*232': 8.788126945495605,'*233': 8.788125991821289,'*234': 8.788125991821289,'*235': 8.7881498336792,'*236': 8.7881498336792,'*237': 8.788171768188477,'*238': 8.788103103637695,'*239': 8.788103103637695,'*24': 8.784937858581543,'*240': 8.787992477416992,'*241': 8.787992477416992,'*242': 8.787992477416992,'*243': 8.787955284118652,'*244': 8.787955284118652,'*245': 8.787955284118652,'*246': 8.78792667388916,'*247': 8.78799057006836,'*248': 8.78799057006836,'*249': 8.78799057006836,'*25': 8.784937858581543,'*250': 8.787979125976562,'*251': 8.787979125976562,'*252': 8.787959098815918,'*253': 8.787900924682617,'*254': 8.787900924682617,'*255': 8.787900924682617,'*256': 8.787896156311035,'*257': 8.787896156311035,'*258': 8.787896156311035,'*259': 8.787790298461914,'*26': 8.784937858581543,'*260': 8.787022590637207,'*261': 8.78685474395752,'*262': 8.78685474395752,'*263': 8.786852836608887,'*264': 8.786852836608887,'*265': 8.786852836608887,'*266': 8.786852836608887,'*267': 8.786852836608887,'*268': 8.786852836608887,'*269': 8.786766052246094,'*27': 8.784971237182617,'*270': 8.786766052246094,'*271': 8.786684036254883,'*272': 8.786669731140137,'*273': 8.786669731140137,'*274': 8.786669731140137,'*275': 8.786669731140137,'*276': 8.786654472351074,'*277': 8.786654472351074,'*278': 8.786654472351074,'*279': 8.786654472351074,'*28': 8.784972190856934,'*280': 8.786645889282227,'*281': 8.786645889282227,'*282': 8.786632537841797,'*283': 8.786568641662598,'*284': 8.786568641662598,'*285': 8.786568641662598,'*286': 8.786824226379395,'*287': 8.786733627319336,'*288': 8.786733627319336,'*289': 8.786630630493164,'*29': 8.784966468811035,'*290': 8.784354209899902,'*291': 8.784354209899902,'*292': 8.784302711486816,'*293': 8.784302711486816,'*294': 8.784302711486816,'*295': 8.784248352050781,'*296': 8.784248352050781,'*297': 8.784248352050781,'*298': 8.784184455871582,'*299': 8.784184455871582,'*3': 8.785101890563965,'*30': 8.784966468811035,'*300': 8.784184455871582,'*301': 8.784184455871582,'*302': 8.784184455871582,'*303': 8.784178733825684,'*304': 8.784178733825684,'*305': 8.784178733825684,'*306': 8.784178733825684,'*307': 8.784178733825684,'*308': 8.784167289733887,'*309': 8.784111976623535,'*31': 8.784932136535645,'*310': 8.784111976623535,'*311': 8.784111976623535,'*312': 8.784111976623535,'*313': 8.784111976623535,'*314': 8.784111976623535,'*315': 8.784163475036621,'*316': 8.784079551696777,'*317': 8.784079551696777,'*318': 8.784043312072754,'*319': 8.784043312072754,'*32': 8.784869194030762,'*320': 8.784043312072754,'*321': 8.784043312072754,'*322': 8.784043312072754,'*323': 8.784043312072754,'*324': 8.784043312072754,'*325': 8.784043312072754,'*326': 8.784043312072754,'*327': 8.784043312072754,'*328': 8.784043312072754,'*329': 8.784043312072754,'*33': 8.784869194030762,'*330': 8.784043312072754,'*331': 8.784015655517578,'*332': 8.784015655517578,'*333': 8.784015655517578,'*334': 8.784015655517578,'*335': 8.784015655517578,'*336': 8.78390884399414,'*337': 8.78390884399414,'*338': 8.783940315246582,'*339': 8.783921241760254,'*34': 8.784869194030762,'*340': 8.783920288085938,'*341': 8.783920288085938,'*342': 8.783920288085938,'*343': 8.783920288085938,'*344': 8.783839225769043,'*345': 8.783827781677246,'*346': 8.783804893493652,'*347': 8.783832550048828,'*348': 8.783832550048828,'*349': 8.783832550048828,'*35': 8.78470516204834,'*350': 8.78374195098877,'*351': 8.78374195098877,'*352': 8.78374195098877,'*353': 8.78374195098877,'*354': 8.783639907836914,'*355': 8.783639907836914,'*356': 8.783639907836914,'*357': 8.783632278442383,'*358': 8.783632278442383,'*359': 8.783564567565918,'*36': 8.784695625305176,'*360': 8.783581733703613,'*361': 8.783604621887207,'*362': 8.783604621887207,'*363': 8.783622741699219,'*364': 8.783622741699219,'*365': 8.783622741699219,'*366': 8.783622741699219,'*367': 8.783618927001953,'*368': 8.783618927001953,'*369': 8.780170440673828,'*37': 8.784588813781738,'*370': 8.780170440673828,'*371': 8.780085563659668,'*372': 8.780085563659668,'*373': 8.778106689453125,'*374': 8.778304100036621,'*375': 8.778304100036621,'*376': 8.778304100036621,'*377': 8.778304100036621,'*378': 8.77828598022461,'*379': 8.778231620788574,'*38': 8.784588813781738,'*380': 8.778181076049805,'*381': 8.778120040893555,'*382': 8.780081748962402,'*383': 8.780081748962402,'*384': 8.780102729797363,'*385': 8.780102729797363,'*386': 8.78001880645752,'*387': 8.780008316040039,'*388': 8.78000545501709,'*389': 8.779960632324219,'*39': 8.784561157226562,'*390': 8.779960632324219,'*391': 8.779960632324219,'*392': 8.779960632324219,'*393': 8.779960632324219,'*394': 8.779960632324219,'*395': 8.779943466186523,'*396': 8.779943466186523,'*397': 8.779932022094727,'*398': 8.77990436553955,'*399': 8.77990436553955,'*4': 8.785101890563965,'*40': 8.784561157226562,'*400': 8.77990436553955,'*401': 8.77990436553955,'*402': 8.77990436553955,'*403': 8.779886245727539,'*404': 8.77989673614502,'*405': 8.77989673614502,'*406': 8.779945373535156,'*407': 8.779945373535156,'*408': 8.779945373535156,'*409': 8.779945373535156,'*41': 8.784561157226562,'*410': 8.779945373535156,'*411': 8.779945373535156,'*412': 8.779945373535156,'*413': 8.7798433303833,'*414': 8.7798433303833,'*415': 8.77976131439209,'*416': 8.77976131439209,'*417': 8.77976131439209,'*418': 8.77976131439209,'*419': 8.77976131439209,'*42': 8.784456253051758,'*420': 8.77976131439209,'*421': 8.77976131439209,'*422': 8.77976131439209,'*423': 8.77976131439209,'*424': 8.77976131439209,'*425': 8.77976131439209,'*426': 8.77978229522705,'*427': 8.779766082763672,'*428': 8.779631614685059,'*429': 8.779631614685059,'*43': 8.784456253051758,'*430': 8.779631614685059,'*431': 8.779647827148438,'*432': 8.779647827148438,'*433': 8.779633522033691,'*434': 8.779633522033691,'*435': 8.779633522033691,'*436': 8.779633522033691,'*437': 8.779633522033691,'*438': 8.779633522033691,'*439': 8.779633522033691,'*44': 8.784456253051758,'*440': 8.779603004455566,'*441': 8.779603004455566,'*442': 8.779603004455566,'*443': 8.779603004455566,'*444': 8.779594421386719,'*445': 8.779594421386719,'*446': 8.779594421386719,'*447': 8.779594421386719,'*448': 8.779583930969238,'*449': 8.779583930969238,'*45': 8.784456253051758,'*450': 8.779583930969238,'*451': 8.779597282409668,'*452': 8.779597282409668,'*453': 8.779597282409668,'*454': 8.779558181762695,'*455': 8.779566764831543,'*456': 8.779552459716797,'*457': 8.779552459716797,'*458': 8.779552459716797,'*459': 8.77951431274414,'*46': 8.784456253051758,'*460': 8.77951431274414,'*461': 8.779766082763672,'*462': 8.779766082763672,'*463': 8.779766082763672,'*464': 8.779766082763672,'*465': 8.779766082763672,'*466': 8.779766082763672,'*467': 8.779766082763672,'*468': 8.779766082763672,'*469': 8.779766082763672,'*47': 8.78431510925293,'*470': 8.779766082763672,'*471': 8.779766082763672,'*472': 8.779766082763672,'*473': 8.779766082763672,'*474': 8.779766082763672,'*475': 8.779766082763672,'*476': 8.779766082763672,'*477': 8.779644966125488,'*478': 8.77968978881836,'*479': 8.779638290405273,'*48': 8.78429126739502,'*480': 8.779638290405273,'*481': 8.779638290405273,'*482': 8.779638290405273,'*483': 8.779622077941895,'*484': 8.779622077941895,'*485': 8.779573440551758,'*486': 8.779573440551758,'*487': 8.779573440551758,'*488': 8.779573440551758,'*489': 8.779485702514648,'*49': 8.784265518188477,'*490': 8.779485702514648,'*491': 8.779485702514648,'*492': 8.779485702514648,'*493': 8.779485702514648,'*494': 8.779372215270996,'*495': 8.779372215270996,'*496': 8.779366493225098,'*497': 8.779366493225098,'*498': 8.779366493225098,'*499': 8.779346466064453,'*5': 8.785101890563965,'*50': 8.784265518188477,'*500': 8.779346466064453,'*501': 8.779346466064453,'*502': 8.779346466064453,'*503': 8.779345512390137,'*504': 8.779332160949707,'*505': 8.779332160949707,'*506': 8.779332160949707,'*507': 8.779332160949707,'*51': 8.784265518188477,'*52': 8.784265518188477,'*53': 8.784265518188477,'*54': 8.784265518188477,'*55': 8.784265518188477,'*56': 8.784265518188477,'*57': 8.784265518188477,'*58': 8.784209251403809,'*59': 8.784209251403809,'*6': 8.785005569458008,'*60': 8.784209251403809,'*61': 8.784209251403809,'*62': 8.784217834472656,'*63': 8.784217834472656,'*64': 8.784217834472656,'*65': 8.784133911132812,'*66': 8.784133911132812,'*67': 8.784133911132812,'*68': 8.784133911132812,'*69': 8.784133911132812,'*7': 8.785005569458008,'*70': 8.784052848815918,'*71': 8.78404426574707,'*72': 8.78404426574707,'*73': 8.783910751342773,'*74': 8.783910751342773,'*75': 8.783910751342773,'*76': 8.783910751342773,'*77': 8.783930778503418,'*78': 8.783914566040039,'*79': 8.783914566040039,'*8': 8.785005569458008,'*80': 8.783914566040039,'*81': 8.783904075622559,'*82': 8.783904075622559,'*83': 8.783904075622559,'*84': 8.783904075622559,'*85': 8.783904075622559,'*86': 8.783838272094727,'*87': 8.785196304321289,'*88': 8.785196304321289,'*89': 8.785196304321289,'*9': 8.785005569458008,'*90': 8.78512954711914,'*91': 8.78512954711914,'*92': 8.78512954711914,'*93': 8.78512954711914,'*94': 8.785123825073242,'*95': 8.785123825073242,'*96': 8.785123825073242,'*97': 8.785099029541016,'*98': 8.785099029541016,'*99': 8.785499572753906}

            # test_mosaic_cube_briggsbwtaper.exp_pa_dict
            exp_pa_dict = {'*0': 67.43175506591797,'*1': 67.43175506591797,'*10': 67.43268585205078,'*100': 67.47151184082031,'*101': 67.47151184082031,'*102': 67.47151184082031,'*103': 67.47151184082031,'*104': 67.47151184082031,'*105': 67.47151184082031,'*106': 67.47151184082031,'*107': 67.47262573242188,'*108': 67.47262573242188,'*109': 67.47262573242188,'*11': 67.43269348144531,'*110': 67.47262573242188,'*111': 67.47271728515625,'*112': 67.41316223144531,'*113': 67.41345977783203,'*114': 67.41345977783203,'*115': 67.41345977783203,'*116': 67.41345977783203,'*117': 67.41315460205078,'*118': 67.41349029541016,'*119': 67.41349029541016,'*12': 67.43218231201172,'*120': 67.41333770751953,'*121': 67.41392517089844,'*122': 67.4139175415039,'*123': 67.4139175415039,'*124': 67.4139175415039,'*125': 67.4139175415039,'*126': 67.4139175415039,'*127': 67.4135513305664,'*128': 67.4135513305664,'*129': 67.41309356689453,'*13': 67.43218231201172,'*130': 67.41309356689453,'*131': 67.41338348388672,'*132': 67.41338348388672,'*133': 67.41338348388672,'*134': 67.41338348388672,'*135': 67.41338348388672,'*136': 67.41338348388672,'*137': 67.41381072998047,'*138': 67.41381072998047,'*139': 67.41381072998047,'*14': 67.43218231201172,'*140': 67.41381072998047,'*141': 67.41403198242188,'*142': 67.41403198242188,'*143': 67.41403198242188,'*144': 67.41403198242188,'*145': 67.41403198242188,'*146': 67.41342163085938,'*147': 67.41342163085938,'*148': 67.41342163085938,'*149': 67.41363525390625,'*15': 67.43199157714844,'*150': 67.41363525390625,'*151': 67.41363525390625,'*152': 67.41329956054688,'*153': 67.41329956054688,'*154': 67.41326904296875,'*155': 67.41326904296875,'*156': 67.41326904296875,'*157': 67.41326904296875,'*158': 67.41375732421875,'*159': 67.41375732421875,'*16': 67.43199157714844,'*160': 67.41375732421875,'*161': 67.41375732421875,'*162': 67.41375732421875,'*163': 67.4129409790039,'*164': 67.4129409790039,'*165': 67.4129409790039,'*166': 67.4129409790039,'*167': 67.4129409790039,'*168': 67.4122314453125,'*169': 67.4122314453125,'*17': 67.43199157714844,'*170': 67.4122314453125,'*171': 67.4122314453125,'*172': 67.41289520263672,'*173': 67.41206359863281,'*174': 67.41234588623047,'*175': 67.4119644165039,'*176': 67.40419006347656,'*177': 67.39434814453125,'*178': 67.39404296875,'*179': 67.39404296875,'*18': 67.43222045898438,'*180': 67.39404296875,'*181': 67.39370727539062,'*182': 67.39402770996094,'*183': 67.39402770996094,'*184': 67.39402770996094,'*185': 67.39402770996094,'*186': 67.39402770996094,'*187': 67.39402770996094,'*188': 67.39402770996094,'*189': 67.39391326904297,'*19': 67.43222045898438,'*190': 67.39391326904297,'*191': 67.39449310302734,'*192': 67.39449310302734,'*193': 67.39432525634766,'*194': 67.39451599121094,'*195': 67.39451599121094,'*196': 67.39451599121094,'*197': 67.39404296875,'*198': 67.39404296875,'*199': 67.39404296875,'*2': 67.43175506591797,'*20': 67.43222045898438,'*200': 67.39353942871094,'*201': 67.39353942871094,'*202': 67.39353942871094,'*203': 67.39353942871094,'*204': 67.33964538574219,'*205': 67.33964538574219,'*206': 67.4762954711914,'*207': 67.4762954711914,'*208': 67.47795867919922,'*209': 67.47826385498047,'*21': 67.43222045898438,'*210': 67.47826385498047,'*211': 67.47830963134766,'*212': 67.47745513916016,'*213': 67.4776382446289,'*214': 67.4776382446289,'*215': 67.47718811035156,'*216': 67.47718811035156,'*217': 67.47718811035156,'*218': 67.47718811035156,'*219': 67.47740936279297,'*22': 67.43081665039062,'*220': 67.47457122802734,'*221': 67.44305419921875,'*222': 67.44305419921875,'*223': 67.44305419921875,'*224': 67.44305419921875,'*225': 67.44090270996094,'*226': 67.44090270996094,'*227': 67.44090270996094,'*228': 67.44090270996094,'*229': 67.44090270996094,'*23': 67.43081665039062,'*230': 67.44090270996094,'*231': 67.44050598144531,'*232': 67.44050598144531,'*233': 67.44050598144531,'*234': 67.44050598144531,'*235': 67.44041442871094,'*236': 67.44041442871094,'*237': 67.4404067993164,'*238': 67.44001770019531,'*239': 67.44001770019531,'*24': 67.43081665039062,'*240': 67.44052124023438,'*241': 67.44051361083984,'*242': 67.44051361083984,'*243': 67.4412612915039,'*244': 67.4412612915039,'*245': 67.4412612915039,'*246': 67.44178771972656,'*247': 67.4410400390625,'*248': 67.4410400390625,'*249': 67.4410400390625,'*25': 67.43081665039062,'*250': 67.44133758544922,'*251': 67.44133758544922,'*252': 67.44181060791016,'*253': 67.44148254394531,'*254': 67.44148254394531,'*255': 67.44148254394531,'*256': 67.44193267822266,'*257': 67.44193267822266,'*258': 67.44193267822266,'*259': 67.44180297851562,'*26': 67.43081665039062,'*260': 67.4464340209961,'*261': 67.45145416259766,'*262': 67.45145416259766,'*263': 67.45181274414062,'*264': 67.45181274414062,'*265': 67.45181274414062,'*266': 67.45181274414062,'*267': 67.45181274414062,'*268': 67.45181274414062,'*269': 67.45170593261719,'*27': 67.43086242675781,'*270': 67.45170593261719,'*271': 67.45174407958984,'*272': 67.45155334472656,'*273': 67.45155334472656,'*274': 67.45155334472656,'*275': 67.45155334472656,'*276': 67.45182037353516,'*277': 67.45182037353516,'*278': 67.45182037353516,'*279': 67.45182037353516,'*28': 67.43086242675781,'*280': 67.4520263671875,'*281': 67.4520263671875,'*282': 67.45232391357422,'*283': 67.45181274414062,'*284': 67.45181274414062,'*285': 67.45181274414062,'*286': 67.45014190673828,'*287': 67.44970703125,'*288': 67.44970703125,'*289': 67.44854736328125,'*29': 67.4310531616211,'*290': 67.48199462890625,'*291': 67.48199462890625,'*292': 67.48190307617188,'*293': 67.48190307617188,'*294': 67.48190307617188,'*295': 67.48151397705078,'*296': 67.48151397705078,'*297': 67.48151397705078,'*298': 67.48177337646484,'*299': 67.48177337646484,'*3': 67.43208312988281,'*30': 67.4310531616211,'*300': 67.48177337646484,'*301': 67.48177337646484,'*302': 67.48177337646484,'*303': 67.4820556640625,'*304': 67.4820556640625,'*305': 67.4820556640625,'*306': 67.4820556640625,'*307': 67.4820556640625,'*308': 67.48234558105469,'*309': 67.48200988769531,'*31': 67.43025970458984,'*310': 67.48200988769531,'*311': 67.48200988769531,'*312': 67.48200988769531,'*313': 67.48200988769531,'*314': 67.48200988769531,'*315': 67.48149871826172,'*316': 67.48054504394531,'*317': 67.48054504394531,'*318': 67.47996520996094,'*319': 67.47996520996094,'*32': 67.42958068847656,'*320': 67.47996520996094,'*321': 67.47996520996094,'*322': 67.47996520996094,'*323': 67.47996520996094,'*324': 67.47996520996094,'*325': 67.47996520996094,'*326': 67.47996520996094,'*327': 67.47996520996094,'*328': 67.47996520996094,'*329': 67.47996520996094,'*33': 67.42958068847656,'*330': 67.47996520996094,'*331': 67.48100280761719,'*332': 67.48100280761719,'*333': 67.48100280761719,'*334': 67.48100280761719,'*335': 67.48100280761719,'*336': 67.4832763671875,'*337': 67.4832763671875,'*338': 67.48302459716797,'*339': 67.48219299316406,'*34': 67.42958068847656,'*340': 67.48219299316406,'*341': 67.48219299316406,'*342': 67.48219299316406,'*343': 67.48219299316406,'*344': 67.48584747314453,'*345': 67.48619842529297,'*346': 67.48649597167969,'*347': 67.48668670654297,'*348': 67.48701477050781,'*349': 67.48701477050781,'*35': 67.42930603027344,'*350': 67.48657989501953,'*351': 67.48657989501953,'*352': 67.48657989501953,'*353': 67.48657989501953,'*354': 67.48670959472656,'*355': 67.48670959472656,'*356': 67.48670959472656,'*357': 67.48712158203125,'*358': 67.48712158203125,'*359': 67.4869384765625,'*36': 67.42963409423828,'*360': 67.4870376586914,'*361': 67.48763275146484,'*362': 67.48763275146484,'*363': 67.48761749267578,'*364': 67.48761749267578,'*365': 67.48761749267578,'*366': 67.48761749267578,'*367': 67.48780822753906,'*368': 67.48780822753906,'*369': 67.52494812011719,'*37': 67.43215942382812,'*370': 67.52494812011719,'*371': 67.5248031616211,'*372': 67.5248031616211,'*373': 67.53797149658203,'*374': 67.5359115600586,'*375': 67.5359115600586,'*376': 67.5359115600586,'*377': 67.5359115600586,'*378': 67.53559112548828,'*379': 67.53504943847656,'*38': 67.43215942382812,'*380': 67.53511810302734,'*381': 67.52452087402344,'*382': 67.5870361328125,'*383': 67.5870361328125,'*384': 67.58714294433594,'*385': 67.58714294433594,'*386': 67.58638000488281,'*387': 67.58667755126953,'*388': 67.58697509765625,'*389': 67.5877914428711,'*39': 67.43267822265625,'*390': 67.5877914428711,'*391': 67.5877914428711,'*392': 67.5877914428711,'*393': 67.5877914428711,'*394': 67.5877914428711,'*395': 67.58808135986328,'*396': 67.58808135986328,'*397': 67.58843231201172,'*398': 67.58810424804688,'*399': 67.58810424804688,'*4': 67.43208312988281,'*40': 67.43267822265625,'*400': 67.58810424804688,'*401': 67.58810424804688,'*402': 67.58810424804688,'*403': 67.58843994140625,'*404': 67.58844757080078,'*405': 67.58844757080078,'*406': 67.5879898071289,'*407': 67.5879898071289,'*408': 67.5879898071289,'*409': 67.5879898071289,'*41': 67.43267822265625,'*410': 67.5879898071289,'*411': 67.5879898071289,'*412': 67.5879898071289,'*413': 67.58783721923828,'*414': 67.58783721923828,'*415': 67.58724212646484,'*416': 67.58724212646484,'*417': 67.58724212646484,'*418': 67.58724212646484,'*419': 67.58724212646484,'*42': 67.43192291259766,'*420': 67.58724212646484,'*421': 67.58724212646484,'*422': 67.58724212646484,'*423': 67.58724212646484,'*424': 67.58724212646484,'*425': 67.58724212646484,'*426': 67.58766174316406,'*427': 67.58784484863281,'*428': 67.58757019042969,'*429': 67.58757019042969,'*43': 67.43192291259766,'*430': 67.58757019042969,'*431': 67.58688354492188,'*432': 67.58688354492188,'*433': 67.58717346191406,'*434': 67.58717346191406,'*435': 67.58717346191406,'*436': 67.58717346191406,'*437': 67.58717346191406,'*438': 67.58717346191406,'*439': 67.58716583251953,'*44': 67.43192291259766,'*440': 67.58800506591797,'*441': 67.58800506591797,'*442': 67.58800506591797,'*443': 67.58800506591797,'*444': 67.58821105957031,'*445': 67.58821105957031,'*446': 67.58821105957031,'*447': 67.58821105957031,'*448': 67.58856964111328,'*449': 67.58856964111328,'*45': 67.43192291259766,'*450': 67.58856964111328,'*451': 67.58806610107422,'*452': 67.58806610107422,'*453': 67.58806610107422,'*454': 67.58890533447266,'*455': 67.58894348144531,'*456': 67.58950805664062,'*457': 67.58950805664062,'*458': 67.58950805664062,'*459': 67.59008026123047,'*46': 67.43192291259766,'*460': 67.59008026123047,'*461': 67.58943176269531,'*462': 67.58943176269531,'*463': 67.58943176269531,'*464': 67.58943176269531,'*465': 67.58943176269531,'*466': 67.58943176269531,'*467': 67.58943176269531,'*468': 67.58943176269531,'*469': 67.58943176269531,'*47': 67.43159484863281,'*470': 67.58943176269531,'*471': 67.58943176269531,'*472': 67.58943176269531,'*473': 67.58943176269531,'*474': 67.58943176269531,'*475': 67.58943176269531,'*476': 67.58943176269531,'*477': 67.58888244628906,'*478': 67.5902328491211,'*479': 67.58990478515625,'*48': 67.4321060180664,'*480': 67.58990478515625,'*481': 67.58990478515625,'*482': 67.58990478515625,'*483': 67.590087890625,'*484': 67.590087890625,'*485': 67.59111785888672,'*486': 67.59111785888672,'*487': 67.59111785888672,'*488': 67.59111785888672,'*489': 67.59101104736328,'*49': 67.43173217773438,'*490': 67.59101104736328,'*491': 67.59101104736328,'*492': 67.59101104736328,'*493': 67.59101104736328,'*494': 67.59122467041016,'*495': 67.59122467041016,'*496': 67.59153747558594,'*497': 67.59153747558594,'*498': 67.59153747558594,'*499': 67.59127044677734,'*5': 67.43208312988281,'*50': 67.43173217773438,'*500': 67.59127044677734,'*501': 67.59127044677734,'*502': 67.59127044677734,'*503': 67.59126281738281,'*504': 67.59169006347656,'*505': 67.59169006347656,'*506': 67.59169006347656,'*507': 67.59169006347656,'*51': 67.43173217773438,'*52': 67.43173217773438,'*53': 67.43173217773438,'*54': 67.43173217773438,'*55': 67.43173217773438,'*56': 67.43173217773438,'*57': 67.43173217773438,'*58': 67.43089294433594,'*59': 67.43089294433594,'*6': 67.43180084228516,'*60': 67.43089294433594,'*61': 67.43089294433594,'*62': 67.43035888671875,'*63': 67.43035888671875,'*64': 67.43035888671875,'*65': 67.42960357666016,'*66': 67.42960357666016,'*67': 67.42960357666016,'*68': 67.42960357666016,'*69': 67.42960357666016,'*7': 67.43180084228516,'*70': 67.42900848388672,'*71': 67.42920684814453,'*72': 67.42920684814453,'*73': 67.42935180664062,'*74': 67.42935180664062,'*75': 67.42935180664062,'*76': 67.42935180664062,'*77': 67.42913055419922,'*78': 67.42965698242188,'*79': 67.42965698242188,'*8': 67.43180084228516,'*80': 67.42965698242188,'*81': 67.4312973022461,'*82': 67.4312973022461,'*83': 67.4312973022461,'*84': 67.4312973022461,'*85': 67.4312973022461,'*86': 67.4311294555664,'*87': 67.46540832519531,'*88': 67.46540832519531,'*89': 67.46540832519531,'*9': 67.43180084228516,'*90': 67.46522521972656,'*91': 67.46522521972656,'*92': 67.46522521972656,'*93': 67.46522521972656,'*94': 67.46541595458984,'*95': 67.46541595458984,'*96': 67.46541595458984,'*97': 67.46488952636719,'*98': 67.46488952636719,'*99': 67.47151184082031}

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed = self.filter_report(report)

        img = shutil._basename(img)
        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, test_name, self.parallel)

        test_dict['test_mosaic_cube_briggsbwtaper']['report'] = report
        test_dict['test_mosaic_cube_briggsbwtaper']['images'] = []

        self.mom8_creator(img+'.image', range_list=[0.15, 1.2])
        self.mom8_creator(img+'.residual', range_list=[0.15, 1.2])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')

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
    # Test 10
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
            perchanweightdensity=False, gridder='mosaic', \
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
            perchanweightdensity=False, gridder='mosaic', \
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
    # Test 11
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
            mosweight=True, usepointing=False, pblimit=0.2, \
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
            mosweight=True, usepointing=False, pblimit=0.2, \
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
    # Test 12a
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_mosaic_cube_eph(self):
        ''' Mosaic ephemeris cube imaging with pcwdF - field Venus, spw 45 '''

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
            perchanweightdensity=False, gridder='mosaic', \
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
            perchanweightdensity=False, gridder='mosaic', \
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
    # Test 12b
    @stats_dict(test_dict)
    def test_mosaic_cube_eph_pcwdT(self):
        ''' Mosaic ephemeris cube imaging with briggsbwtaper - field Venus, spw 45 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name,'test_')+'.iter'
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
            perchanweightdensity=True, gridder='mosaic', \
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
            perchanweightdensity=True, gridder='mosaic', \
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
  
        # test_mosaic_cube_eph_pcwdT
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

        # test_mosaic_cube_eph_pcwdT.exp_im_stats
        exp_im_stats = {'com_bmaj': [False, 0.9379540349590176],
            'com_bmin': [False, 0.7149098943227834],
            'com_pa': [False, -87.11030518848933],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 105415437.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.07764589786529541],
            'max_val_pos': [True, [319, 243, 0, 511]],
            'min_val': [False, -0.06537871807813644],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.011317758055898897],
            'rms_per_chan': [False, [0.013952824312672517, 0.011223327205098879, 0.01036411751622865, 0.010832455981697539, 0.011001718655577645, 0.012215811443519072, 0.01268480761423732, 0.014111828714418134, 0.011716815359172444, 0.011494545651590103, 0.010381387448426992, 0.010197685290497619, 0.012687468807117559, 0.01071566037747204, 0.010218405129633897, 0.00997905869753466, 0.010974153867643105, 0.01261961582579337, 0.00879838605324991, 0.011328126920983988, 0.011807105313732922, 0.011958532751411672, 0.011534909428539419, 0.012489514123597719, 0.011416277677934286, 0.010497340184665157, 0.010539470035923301, 0.010747755863646354, 0.013022671037398481, 0.012127856200895518, 0.010947504464902483, 0.010407856457587449, 0.010568310387002766, 0.009859366202430282, 0.011292797381614357, 0.01120602665517635, 0.010635943981574782, 0.011152124807944974, 0.010417536838069542, 0.010687742549113984, 0.010990171604217426, 0.011667347614963519, 0.011030008829057606, 0.012589538770680378, 0.011462920168516305, 0.010922228554987546, 0.010224620257577808, 0.009836982192368069, 0.010066212613859342, 0.010424307010180016, 0.011265021728864967, 0.011719394603471465, 0.010837379810025081, 0.010923920792423776, 0.010705927682567265, 0.011008248519198283, 0.010624861394264669, 0.014305773074673213, 0.011563415984269916, 0.00921208809813708, 0.009548894999835755, 0.011729367638065251, 0.010125635330541376, 0.01054867912022443, 0.01023352325427546, 0.010544042070910543, 0.010949270811889384, 0.011275444519638194, 0.01053108631785352, 0.00942392918950115, 0.012413345214870029, 0.0116069301291608, 0.010644934011768756, 0.009121797788205231, 0.01050384156097305, 0.012143875755542615, 0.012868463358860945, 0.012654364968631033, 0.012074341110803839, 0.011425164932901773, 0.010293880308579896, 0.010190408607867149, 0.010682109658478366, 0.011430958393335971, 0.010841714240402563, 0.010685178581039297, 0.010402286258579127, 0.010642738017406754, 0.009551526339592026, 0.010796150367443079, 0.010999266850767804, 0.011784822818137843, 0.010922162098344578, 0.009787003856290962, 0.010157967087551642, 0.01145915671680839, 0.011455026543021505, 0.011205240954884553, 0.010950659068542716, 0.009024171824664445, 0.011119777883099412, 0.011918294092352273, 0.012747590677585624, 0.011702895908100587, 0.011998524384729503, 0.012203151702688007, 0.013307091886199078, 0.01134223162713411, 0.012406902058416826, 0.01152484148381434, 0.010517874737052211, 0.010414271625402196, 0.011418733463761599, 0.010575999198963586, 0.010173487718504019, 0.01153432065542299, 0.010565830441642953, 0.010003823824159212, 0.011685277412959596, 0.011090626232779417, 0.00970521281304701, 0.011360437051448121, 0.009552842307196717, 0.011354149037183731, 0.01084545448744486, 0.011335700055106887, 0.010113526848132907, 0.012465823279149564, 0.010335045591247425, 0.009142404099591486, 0.010651906140592684, 0.010301550304724244, 0.009623243517698745, 0.011228152468291248, 0.012205988383916795, 0.010759913680451167, 0.010853463808646668, 0.010695302949308328, 0.009965690165150249, 0.010525641901303801, 0.011079524713215988, 0.010280705676998585, 0.009786192000947028, 0.01011409277808253, 0.008911847199352785, 0.01133437124074135, 0.011101970298572178, 0.009449435170700791, 0.01131953192839757, 0.009682073918958293, 0.010482697924882508, 0.010999641467440938, 0.010521701952727874, 0.010362848663538125, 0.011472211213910531, 0.010742157712604564, 0.011146542364040838, 0.010935968125652323, 0.010261561965420156, 0.013520077135145376, 0.013332670912511555, 0.010761449304700943, 0.009621794437748521, 0.00966268414159511, 0.010284528051079834, 0.009838118735407256, 0.009985420577107941, 0.009071079911384868, 0.009605618751526216, 0.010505332374731467, 0.009306295484045681, 0.010387683455025595, 0.010125081587307539, 0.010623781172863259, 0.011935171419953665, 0.011440649818733916, 0.0095628741923247, 0.011775691372410291, 0.012010921447161744, 0.011737959326522185, 0.009986030652084645, 0.011396031863618553, 0.010837932479341542, 0.012155113664181217, 0.010473041065082817, 0.0101562548787553, 0.01132426204266206, 0.009881383975162548, 0.012006245806477684, 0.010238011264374085, 0.010685260555692001, 0.010168555240532078, 0.011809352644626493, 0.009944294845793433, 0.010525309188935192, 0.011730911535539991, 0.010055350768496927, 0.01157771693079674, 0.0101092477719043, 0.01129598250113666, 0.011920580739780565, 0.013782891848275932, 0.012445845761524712, 0.01092555721450914, 0.010140171543239549, 0.009293499389782895, 0.009962632340998355, 0.010220152391767715, 0.011239342830666245, 0.01084923441182909, 0.011495887865788071, 0.011757150790091308, 0.011109284744976775, 0.010793689523644512, 0.00984614186034825, 0.011962597539343091, 0.01165554957026697, 0.010143279064414751, 0.010574676746458447, 0.012418837256078848, 0.012533292537981756, 0.011012268934608039, 0.01108488519052533, 0.012726957339984891, 0.01180816085984399, 0.009868731545120971, 0.011452510582904031, 0.011125189803021945, 0.011483996681784193, 0.01037763371474214, 0.009916189112952474, 0.012562106474178993, 0.01150132769338838, 0.011000935336306202, 0.009177713198635083, 0.011857228118415093, 0.0115243439704556, 0.011073291962102393, 0.010684357065035448, 0.013312686912497545, 0.0117057077536855, 0.009170557306868811, 0.009445814531668774, 0.010816518706640793, 0.01146311204550163, 0.010052530436927857, 0.010415957724708793, 0.010443502922100157, 0.011212573663772945, 0.011249967595220981, 0.01253815613397136, 0.010319388830770914, 0.00928336545678846, 0.01149773523068167, 0.011520549091695398, 0.01223860069815009, 0.011868495320598037, 0.010203453014934, 0.010648109456837898, 0.010297098687063192, 0.010997194695480187, 0.010954763627491572, 0.012343371454359708, 0.010724721644302446, 0.01087517214501339, 0.011523989807723513, 0.012020740855099772, 0.009987132648173734, 0.011454034648731902, 0.012182733321934905, 0.010735649017090267, 0.010266666269737136, 0.010617699696872837, 0.011268380153436237, 0.010401768633923796, 0.010587938980920227, 0.011853474102889755, 0.01246460337617053, 0.012291269167647285, 0.011004360353158896, 0.010019182378432284, 0.010669655273471983, 0.010971142299029952, 0.008613537942785899, 0.01096145572560422, 0.011433160429406735, 0.011587217978381347, 0.01135552283119047, 0.011563672604530095, 0.012426018246877425, 0.010060963580743598, 0.010034211953962284, 0.010882484927333376, 0.010844828477259247, 0.011535616128216611, 0.009552619997280907, 0.010473875541195332, 0.010280134790309776, 0.01119423208131282, 0.010053917198381878, 0.012532771224598636, 0.011140105633137394, 0.011675770834832597, 0.014912054355573212, 0.01501051856998091, 0.014954332295666546, 0.010070430201772824, 0.011387762084968063, 0.01047345543327577, 0.012957616190149278, 0.010495115064887253, 0.010700756794468037, 0.011838486763406554, 0.012196800790136634, 0.011747988623876115, 0.010992003552792309, 0.011105210237915096, 0.013622452020116603, 0.011466159835297515, 0.011239226908512047, 0.011142274150407838, 0.014798838352234478, 0.011209748251748206, 0.01102200221778611, 0.010642754263368269, 0.012562160567806655, 0.013473057156237487, 0.011580135812968389, 0.010120074705931355, 0.011060698889996202, 0.010564693813279097, 0.012423393895583324, 0.0110296900574634, 0.011570593470325161, 0.011454071654125873, 0.010201686153124615, 0.013323394333397343, 0.012065384547996508, 0.010767297608160451, 0.013084673972560803, 0.010352944971730187, 0.012728744853484195, 0.00943615650080395, 0.010533203849201138, 0.011378164080512588, 0.012932503497605972, 0.010260831384463404, 0.011570263239871856, 0.009935504019897646, 0.010467526254547955, 0.011543680068866033, 0.01229990299372627, 0.011998336188785624, 0.01115934471241602, 0.011097450948225791, 0.010372803652759107, 0.012166490875411965, 0.011562045107424375, 0.012095885417249222, 0.011516736594259596, 0.010816394082419818, 0.01120303739243391, 0.012964269825197204, 0.013377500424450907, 0.010378092687060046, 0.010898983212619187, 0.010409836720495879, 0.011374006736210113, 0.012511649126982553, 0.01075244774789806, 0.011450923685264533, 0.012850182321387377, 0.01194086379476733, 0.010344985051262981, 0.011885587909198012, 0.011188050135460786, 0.010625832387479208, 0.011386118588732312, 0.01116281235553508, 0.009543583181773263, 0.010363575982559813, 0.010307364722482126, 0.012092911799169059, 0.010384394573573588, 0.011647625332981968, 0.010492211196921324, 0.011901348744574314, 0.012810961511228012, 0.010333351672581046, 0.012749655798956039, 0.01195351169484717, 0.011000007833491171, 0.010464063716056963, 0.011205948981761976, 0.011025417429002165, 0.0097992426108751, 0.009815036362815664, 0.0110297022277036, 0.011357887772969604, 0.01468585782824629, 0.011134265942934042, 0.011137536981130378, 0.011093592530829855, 0.012205842461159331, 0.010969472567821026, 0.00969020343949966, 0.010773392423826316, 0.013231136488352287, 0.014702163715618968, 0.010814244748936684, 0.011382508673510788, 0.01129185030333294, 0.012449566983247166, 0.011463551180731892, 0.011121432575804089, 0.012900263235359623, 0.012400935694124515, 0.010964874501608952, 0.011398556449387863, 0.010858640784838862, 0.010514327355180828, 0.008914628867965069, 0.010274499271704157, 0.010192334388266899, 0.011189152871131482, 0.011188946913221053, 0.011271733595041047, 0.012217987979284313, 0.010802688399001139, 0.010907874329411992, 0.011275507048271077, 0.011043166354049899, 0.011305268157042778, 0.011299533345062476, 0.010401949169261219, 0.013047844348972945, 0.011282719715366263, 0.013209160272599692, 0.011810797616627147, 0.011595215632218065, 0.011413269867523071, 0.01365052459812259, 0.012883364752765226, 0.012025977955336994, 0.010501096913984027, 0.011192618075677481, 0.011203649623950469, 0.011034526432307436, 0.010641353768992378, 0.01141431476803484, 0.013100434709351618, 0.0132338630268423, 0.011533110713007843, 0.012416137624019915, 0.011579819861835027, 0.010401665748566746, 0.010659537404489949, 0.009944888037186038, 0.012224336022255065, 0.011123894707755462, 0.012565107101289547, 0.011219020869422289, 0.011134595738613015, 0.011374971202626281, 0.009909127716957617, 0.011525711053936363, 0.010540163407190151, 0.010664024561240474, 0.011416735928609714, 0.011125994671102373, 0.012293742718393012, 0.010649615183432362, 0.01127886840705332, 0.010656485468072808, 0.011005638434076419, 0.01142766714380127, 0.011183125156206529, 0.012644247687324821, 0.01087342596802822, 0.010923673836102306, 0.011377783644584957, 0.011095605485000949, 0.011239334909628433, 0.015128129992789618, 0.01123769143338418, 0.011950181234987364, 0.011943356050185013, 0.011943032808770274, 0.011435624764665307, 0.011693797052649172, 0.011863035240482633, 0.010802821549491798, 0.011530343914101697, 0.012220500299541063, 0.010222076675163698, 0.011545056691706133, 0.011541475089686585, 0.011365000012224534, 0.013129178661257903, 0.012109432485834133, 0.011705998603036942, 0.01332801222667836, 0.010421579746041753, 0.012332299119259166, 0.011695921851666293, 0.011861382037600104, 0.010733953868731591, 0.012443838237761275, 0.011749537257781906, 0.01127972820046856, 0.012417150758772474, 0.012876061405771487, 0.01229329706234427, 0.010463342639651282, 0.012452727831893464, 0.012780593608205014, 0.01235331635034423, 0.013470487854621271, 0.012988951448184115, 0.01199799885516839, 0.010079337518303772, 0.01287819088183971, 0.012074009169895558, 0.012024851076463856, 0.012357613757578743, 0.00986387398146561, 0.011156155647000541, 0.012986487007673282, 0.011167787621532713, 0.011692299307396932, 0.01077714725559941, 0.011058888782912312, 0.010765026318804142, 0.011166103216687135, 0.011894043626894535, 0.011942375440198388, 0.012839833776826164, 0.010718647520225305, 0.010835141114292672, 0.011851477300355531, 0.011742386487314387, 0.011066190044379092, 0.011247274314614874, 0.01174393228440886, 0.011267708906861922, 0.01148375256025399, 0.010722110221756747, 0.012309630916545324, 0.011671125387304273, 0.011420823180310237, 0.011403318068867855, 0.01130301052757781, 0.011638983004257196, 0.01338697924283547, 0.012343605832537987, 0.01205347572658705, 0.01098976767153549, 0.011437737818735963, 0.012124140759087246, 0.010459835531461544, 0.012123752899503785, 0.01336596520486306, 0.012123825918338877, 0.011378014223326372, 0.012291964825954542, 0.010722933804027448, 0.012184264858052407, 0.013810457001453188, 0.012096782298768602, 0.011228295508276611, 0.011085783824281143, 0.011051370469863544, 0.011083587614460051, 0.011909818273358692, 0.01263710350914078, 0.010442496793236982, 0.0114785080445132, 0.01149669773287715, 0.010908081354985178, 0.012178331609156838, 0.01215530340263708, 0.013399386053498358, 0.011565088133306681, 0.012284147517492032, 0.013220015139400316, 0.011110081108814026, 0.010522513543356152, 0.01210904367528662, 0.012046601276756446, 0.013750160295284754, 0.014029966993857815, 0.012708091908788645, 0.012132924215589775, 0.013660640733947167, 0.012448956372369237, 0.012152947682804217, 0.0110001504648406, 0.012302657614990867, 0.011025537510597195, 0.011086762934324377, 0.013494205946309241, 0.013143644557605886, 0.010516688032448693, 0.011260213564550594, 0.011306432683006412, 0.011216093747808395, 0.012084504685375964, 0.010130958442808402, 0.010836743616453778, 0.011195852477542742, 0.013314054355506207, 0.012468563281596363, 0.014531971439973198, 0.011234239975507543, 0.010679711206529123, 0.012910813782816036, 0.010771735485643015, 0.012959759597824645, 0.012811319545115869, 0.012840670288131933, 0.010900885481767994, 0.012187926874856899, 0.012086530391646822, 0.012544838111611417, 0.0125511334619372, 0.011126023800626897, 0.01291878296294955, 0.01201853675003584, 0.010459174616133314, 0.012834965996958124, 0.013413063733023372, 0.010708581807417157, 0.012463970379745917, 0.011024775587372719, 0.010075416119290878, 0.013061769383076216, 0.011447959812710486, 0.0106457178615456, 0.011723198358367916, 0.011838764905609225, 0.0116146498615186, 0.012524575628961108, 0.011349320835877977, 0.012961014566638878, 0.011929326921873174, 0.013290113130593873, 0.011110509338503912, 0.010329372233588555, 0.01255146380043951, 0.011815514998235147, 0.010548381071779899, 0.012762420175408174, 0.010751422768256316, 0.012323676128727633, 0.011094263613742716, 0.011311146546259658, 0.012124538487397554, 0.01148651691264351, 0.012116351428395063, 0.01134095895600492, 0.011880782764075419, 0.012147392302590055, 0.013804766040526675, 0.01277474875861892, 0.013228435205014353, 0.01110706693517133, 0.011651182186623341, 0.011157059065348567, 0.012195062239501939, 0.011224814709552586, 0.009603476714627315, 0.010112963894491367, 0.011131771188546889, 0.011436980318400864, 0.01450323953068241, 0.01097846080917362, 0.011605739673614645, 0.011074826006397639, 0.010675750111599052, 0.014711848007813342, 0.011261796207339771, 0.010719518887139111, 0.011810741112981166, 0.00997974269900419, 0.012314084085449388, 0.011685768447345297, 0.011233492639835585, 0.010423496093417918, 0.01212588767332324, 0.011017966854857396, 0.009797916651975491, 0.012367729331982075, 0.010452487308931509, 0.010371428083309906, 0.011128038848013659, 0.012207250252079914, 0.011933037402319524, 0.012215410651608355, 0.011607850541806762, 0.011545594804284295, 0.0118246051623973, 0.011635044046742185, 0.011901504103437679, 0.012935392158118486, 0.011649547056165003, 0.01129820772751571, 0.012138028008517791, 0.010514438523041464, 0.01358604905070241, 0.012481254863205438, 0.012252864397004844, 0.012965824070529413, 0.010426321827602118, 0.010717190217717682, 0.01157053277768032, 0.010242083737897588, 0.011123672082680695, 0.009874965858109803, 0.011138289835489487, 0.01176703814402652, 0.010198826088840618, 0.011789964682078666, 0.012404506523331115, 0.010611664332278843, 0.010378647634358374, 0.012979899001323627, 0.01096991669040994, 0.012791486459182216, 0.010771479489720024, 0.01079691883277352, 0.012046325162459867, 0.011836892692105807, 0.009933998781391123, 0.009722272914188677, 0.009866507344817373, 0.010814948632624244, 0.010675753208637, 0.012103219880598483, 0.010378855925225614, 0.010253829841538224, 0.013177281332760788, 0.011304135583421627, 0.010433311845928465, 0.009796279511695063, 0.01058018738446356, 0.011725972168137094, 0.011634554694044017, 0.010231738214832366, 0.011495808796132446, 0.01396779327996553, 0.01009026723361304, 0.010138283146139059, 0.009628673332888369, 0.010821703457303149, 0.012659600852274402, 0.011311575614730379, 0.010723850715594379, 0.011173622109809794, 0.011527914715418012, 0.011215266541811642, 0.01026689089944772, 0.014046297620689903, 0.011215505058165284, 0.012481376045881789, 0.01485446702804911, 0.01066108597373479, 0.010370778790703018, 0.012429412193305458, 0.012094586932484313, 0.012250363157907257, 0.011946716370198753, 0.010812668134338467, 0.01023598998276974, 0.01191249921891237, 0.011315496232350286, 0.010471754950069773, 0.011600932703482125, 0.010622930565577882, 0.010489782344523966, 0.010311017337175318, 0.012770340802573773, 0.011208804070517598, 0.012131109175132496, 0.010987464573892632, 0.011967223067430375, 0.011970948546115074, 0.011224393276117833, 0.011006568453344507, 0.01198576711889043, 0.01218658424605089, 0.012637704324893652, 0.009600503750844019, 0.010801160346116454, 0.00986866077578768, 0.010735656701440619, 0.010530327636750247, 0.009359327067672016, 0.01287830912945531, 0.010471850866332072, 0.010738959243559904, 0.010003666141775466, 0.01102915499122303, 0.010508010642465081, 0.010849987424999586, 0.010325713630890573, 0.011768643067424927, 0.010118159300559228, 0.01054600505828154, 0.012355373624261817, 0.011111346290705504, 0.011015459743564824, 0.009212744798395595, 0.011329466972907454, 0.009655640336682782, 0.010885639339990122, 0.011076087510704572, 0.011321483162672081, 0.011195051208515832, 0.012761975176240175, 0.01090027558960384, 0.009643628031204924, 0.009989849548270598, 0.011922936664847835, 0.010042803323676255, 0.010282380723108475, 0.00957402975318643, 0.010210843357696415, 0.012248462756202227, 0.013610912133505047, 0.011213983948505417, 0.011134348693513587, 0.012525038632003346, 0.011547419928411509, 0.010860328311161491, 0.011121518000162374, 0.0114137826461746, 0.01163458238392168, 0.010858717313946737, 0.011258034393661432, 0.010722296694329496, 0.01021839900575297, 0.012117431323588852, 0.010650364986139109, 0.012338125022154102, 0.01050816240341804, 0.010403647891587866, 0.011192171638835253, 0.00954567872611519, 0.009236045075260093, 0.009913649734802738, 0.010395303720862954, 0.011042575215920595, 0.010106110083623323, 0.010735068963138965, 0.011263695579921096, 0.011043841547598967, 0.009963117778748613, 0.010707973988856399, 0.010408083495578311, 0.009490454103057131, 0.01080648002610078, 0.012241135842630306, 0.010620553871851712, 0.009862373504115455, 0.009331036350230526, 0.009649550023746711, 0.011940791594072702, 0.011410389489454003, 0.012055541222527388, 0.012336868718636282, 0.010980677652329796, 0.010995474229251482, 0.012038621560508668, 0.011186941028217571, 0.009656071111081561, 0.009974345764455729, 0.009576076393054064, 0.010365676817708228, 0.010106976402648933, 0.011488823365170852, 0.01173293740078898, 0.01095432425152728, 0.011540914266038059, 0.010748436996980188, 0.01286752120964017, 0.012045687568002155, 0.01122847507695563, 0.00930542482549422, 0.0101291326260736, 0.009783844436450386, 0.009866853724654226, 0.010935013957735277, 0.010620708435289816, 0.011057657522831872, 0.009016748530654777, 0.01148003092753894, 0.012271850467752607, 0.010356972212958098, 0.010259579491849207, 0.010476616844666924, 0.009689993790738574, 0.009283575526453475, 0.010493392793737066, 0.01012370211001683, 0.00974804550840488, 0.010047508267840052, 0.009483008337857575, 0.010785311767897343, 0.010650841724220426, 0.010556073780271546, 0.01078959550765961, 0.013063889646218468, 0.01047977460831886, 0.009779918444738836, 0.010642369627587678, 0.01070798534835405, 0.010929944392939194, 0.010786460264211813, 0.013797355651233521, 0.011274012386486659, 0.009732247084459638, 0.010718452764614661, 0.009396182169069069, 0.011851837191195699, 0.010382810525665151, 0.01249163497170802, 0.0108684304316671, 0.009148068447102541, 0.011080172069328603, 0.010075780316782579, 0.009213029057359812, 0.012548189084418225, 0.012889302962971397, 0.010492908601949125, 0.010171151000591627, 0.009709201419742765, 0.010464665989395647, 0.010290300853461984, 0.0098475369890548, 0.010022678790603472, 0.01074246829964346, 0.011195775376208254, 0.010650533392173407, 0.011980640706155829, 0.01240593001230758, 0.010247278904611095, 0.009945054982318895, 0.011287081080246428]],
            'im_sum': [False, 9614.635165501642],
            'regn_sum': [False, 26.24731379647983],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.012549503583922393, 0.01153146849781589, 0.011820639611066072, 0.01187125876570391, 0.011437741462926757, 0.012342459120999507, 0.011678152191266191]],
            'profile': [False, 0.022861428161519218]}

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

        # test_mosaic_cube_eph_pcwdT.exp_mask_stats
        exp_mask_stats = {'npts': [True, 191116800],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'mask_pix': [False, 9396],
            'mask_regns': [True, 36],
            'npts_real': [True, 191116800]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_pcwdT.exp_pb_stats
        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415437.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000001788139343],
            'im_rms': [False, 0.6370611558093784],
            'npts_0.2': [False, [111207, 111211, 111208, 111208, 111206, 111204, 111203, 111206, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111207, 111206, 111206, 111207, 111206, 111205, 111205, 111206, 111205, 111205, 111207, 111206, 111206, 111206, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111204, 111206, 111204, 111206, 111206, 111206, 111206, 111208, 111209, 111207, 111207, 111208, 111207, 111207, 111207, 111207, 111208, 111207, 111207, 111208, 111209, 111209, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111206, 111206, 111205, 111207, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111209, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111208, 111207, 111205, 111205, 111206, 111206, 111205, 111204, 111206, 111206, 111205, 111205, 111206, 111206, 111205, 111206, 111204, 111206, 111206, 111204, 111203, 111204, 111204, 111203, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111203, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111200, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111197, 111198, 111199, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111206, 111205, 111206, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111206, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111203, 111203, 111204, 111202, 111202, 111203, 111202, 111202, 111200, 111203, 111201, 111200, 111200, 111201, 111201, 111201, 111196, 111198, 111196, 111196, 111196, 111197, 111196, 111196, 111196, 111196, 111197, 111197, 111198, 111197, 111199, 111199, 111199, 111198, 111200, 111202, 111200, 111200, 111200, 111201, 111201, 111200, 111199, 111200, 111200, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111200, 111198, 111198, 111198, 111201, 111199, 111199, 111199, 111199, 111199, 111200, 111201, 111200, 111201, 111199, 111199, 111201, 111200, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111203, 111201, 111204, 111204, 111202, 111199, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111202, 111203, 111201, 111199, 111201, 111202, 111203, 111202, 111203, 111205, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111203, 111202, 111201, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111201, 111202, 111198, 111198, 111201, 111199, 111200, 111199, 111200, 111197, 111195, 111194, 111194, 111194, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111199, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111203, 111202, 111203, 111203, 111203, 111204, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111195, 111199, 111199, 111197, 111197, 111196, 111196, 111196, 111199, 111200, 111200, 111198, 111197, 111195, 111196, 111194, 111194, 111196, 111197, 111197, 111196, 111196, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111193, 111194, 111192, 111192, 111192, 111192, 111193, 111190, 111190, 111190, 111191, 111188, 111188, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111184, 111186, 111185, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111197, 111198, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111198, 111196, 111195, 111194, 111194, 111194, 111194, 111193, 111191, 111190, 111192, 111191, 111190, 111188, 111188, 111188, 111190, 111188, 111188, 111186, 111186, 111185, 111185, 111184, 111184, 111185, 111183, 111180, 111179, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111181, 111180, 111179, 111179, 111181, 111185, 111184, 111184, 111185, 111185, 111185, 111185, 111185, 111187, 111187, 111188, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111186, 111188, 111188, 111186, 111185, 111191, 111188, 111187, 111188, 111188, 111188, 111190, 111189, 111192, 111191, 111192, 111190, 111189, 111188, 111189, 111192, 111191, 111192, 111192, 111190, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111194, 111193, 111192, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111197, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111193, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111194, 111195, 111196, 111196, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111196, 111198, 111198, 111197, 111198, 111200, 111198, 111197, 111201, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111196, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111191, 111191, 111194, 111194, 111194, 111193, 111192, 111192, 111190, 111191, 111192, 111190, 111190, 111190, 111192, 111192, 111192, 111192, 111193, 111193, 111193, 111193, 111192, 111193, 111192, 111192, 111197, 111197, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111201, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111198, 111198, 111197, 111197, 111197, 111199, 111199, 111197, 111197, 111198, 111195, 111195, 111194, 111194, 111194, 111193, 111193, 111191, 111191, 111190, 111192, 111192, 111192, 111192, 111192, 111193, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111193, 111193, 111194]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64059, 64060, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64060, 64059, 64059, 64061, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64061, 64061, 64065, 64065, 64063, 64064, 64065, 64063, 64063, 64064, 64064, 64063, 64065, 64066, 64065, 64065, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64065, 64066, 64065, 64065, 64065, 64067, 64067, 64068, 64069, 64069, 64069, 64069, 64070, 64069, 64069, 64069, 64069, 64070, 64070, 64069, 64069, 64068, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64069, 64069, 64069, 64069, 64069, 64068, 64066, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64068, 64066, 64067, 64067, 64067, 64065, 64066, 64065, 64066, 64066, 64067, 64067, 64067, 64066, 64066, 64066, 64064, 64063, 64064, 64065, 64064, 64063, 64059, 64060, 64061, 64064, 64063, 64061, 64058, 64060, 64058, 64059, 64058, 64058, 64059, 64059, 64060, 64059, 64059, 64058, 64058, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64067, 64065, 64067, 64065, 64064, 64064, 64064, 64064, 64064, 64063, 64065, 64066, 64065, 64068, 64068, 64067, 64067, 64066, 64067, 64067, 64067, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64068, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64060, 64060, 64060, 64060, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64058, 64055, 64055, 64056, 64055, 64054, 64054, 64053, 64055, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64060, 64058, 64059, 64058, 64059, 64059, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64056, 64057, 64056, 64056, 64057, 64059, 64058, 64057, 64058, 64057, 64058, 64058, 64057, 64058, 64059, 64056, 64058, 64058, 64058, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64064, 64063, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64060, 64060, 64060, 64062, 64061, 64060, 64061, 64061, 64057, 64058, 64057, 64057, 64057, 64056, 64056, 64064, 64062, 64060, 64061, 64062, 64063, 64062, 64060, 64059, 64059, 64060, 64061, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64063, 64065, 64064, 64064, 64064, 64062, 64063, 64062, 64063, 64059, 64058, 64060, 64058, 64058, 64059, 64059, 64057, 64056, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64061, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64062, 64063, 64064, 64062, 64061, 64063, 64062, 64061, 64062, 64061, 64059, 64060, 64061, 64062, 64059, 64058, 64059, 64058, 64059, 64058, 64058, 64058, 64060, 64058, 64057, 64056, 64055, 64055, 64054, 64055, 64054, 64055, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64052, 64053, 64053, 64053, 64054, 64053, 64053, 64054, 64054, 64054, 64052, 64050, 64050, 64049, 64048, 64048, 64048, 64049, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64058, 64057, 64057, 64058, 64059, 64057, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64050, 64051, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64047, 64046, 64047, 64047, 64043, 64043, 64044, 64044, 64044, 64044, 64044, 64042, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64041, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64043, 64041, 64041, 64042, 64044, 64043, 64043, 64045, 64045, 64044, 64046, 64046, 64048, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64047, 64047, 64051, 64048, 64050, 64049, 64050, 64050, 64051, 64051, 64052, 64051, 64051, 64047, 64048, 64048, 64050, 64050, 64049, 64049, 64052, 64048, 64051, 64051, 64052, 64052, 64049, 64049, 64049, 64049, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64057, 64056, 64056, 64056, 64056, 64055, 64054, 64054, 64055, 64057, 64059, 64059, 64059, 64058, 64057, 64057, 64056, 64056, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64060, 64061, 64061, 64061, 64058, 64057, 64059, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64064, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64062, 64063, 64062, 64060, 64060, 64060, 64059, 64059, 64057, 64055, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64061, 64060, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64064, 64063, 64064, 64062, 64062, 64062, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64062, 64060, 64061, 64061, 64060, 64059, 64058, 64057, 64056, 64055, 64055, 64053, 64052, 64052, 64054, 64053, 64052, 64053, 64052, 64050, 64049, 64049, 64050, 64048, 64049, 64048, 64049, 64050, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64049, 64048, 64048, 64060, 64059, 64059, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64063, 64062, 64061, 64061, 64062, 64062, 64060, 64063, 64065, 64065, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64058, 64059, 64060, 64059, 64059, 64058, 64059, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64062, 64062, 64061, 64061, 64061, 64060, 64060, 64059, 64059, 64061, 64059, 64058, 64058, 64058, 64057, 64056, 64055, 64053, 64053, 64054, 64053, 64050, 64050, 64050, 64052, 64053, 64052, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64054, 64056, 64055, 64055]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1042808716193668, 37.06017594715826, 36.758467524557425]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.8801035135706],
            'fit_pix': [False, [240.7427698068262, 210.14394033791376]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_pcwdT.exp_psf_stats
        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.05497631058096886],
            'min_val_pos': [True, [250, 204, 0, 947]],
            'im_rms': [False, 0.012125671143964905],
            'im_sum': [False, 75.7544672754139],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.8864100738559655, 1.0840855127823963, 0.8443649865524704]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462000557837],
            'fit_pix_0': [False, [239.99836406712743, 209.9973158998747]],
            'fit_1': [False, [0.8870506845789847, 1.083263480243368, 0.8441946239816467]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.8801035135706],
            'fit_pix_1': [False, [239.9985183221882, 209.99724649489312]],
            'fit_2': [False, [0.8876668105211563, 1.083930025052625, 0.8429083762110405]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99558702156276],
            'fit_pix_2': [False, [239.99845384938143, 209.9971631888005]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        # test_mosaic_cube_eph_pcwdT.exp_resid_stats
        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415437.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.06226184964179993],
            'max_val_pos': [True, [320, 214, 0, 617]],
            'min_val': [False, -0.06537871807813644],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.011310927582337987],
            'im_sum': [False, 9265.72619977831],
            'regn_sum': [False, -14.81155890946161],
            'npts_real': [True, 191116800]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_eph_pcwdT)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_eph_pcwdT.exp_model_stats
        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.032684747129678726],
            'max_val_pos': [True, [318, 238, 0, 589]],
            'min_val': [False, -0.0034886628855019808],
            'min_val_pos': [True, [323, 245, 0, 589]],
            'im_rms': [False, 2.2748885556543053e-05],
            'im_sum': [False, 9.005414494196884],
            'regn_sum': [False, 1.1523443398764357],
            'mask_non0': [True, 0],
            'npts_real': [True, 191116800]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, 
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube_eph_briggsbwtaper)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_eph_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 45346.8984375],
            'max_val_pos': [True, [0, 0, 0, 283]],
            'min_val': [False, 45300.625],
            'min_val_pos': [True, [0, 0, 0, 856]],
            'im_rms': [False, 45327.50481915811],
            'npts_real': [True, 948]}

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_cube_eph_pcwdT)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        # test_mosaic_cube_eph_pcwdT.exp_wt_stats
        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.3204350769519806],
            'max_val_pos': [True, [240, 210, 0, 595]],
            'min_val': [False, 7.413938146783039e-05],
            'im_rms': [False, 0.12064656547461543],
            'im_sum': [False, 13934226.29788425],
            'npts_0.2': [False, [111207, 111211, 111208, 111208, 111206, 111204, 111203, 111206, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111207, 111206, 111206, 111207, 111206, 111205, 111205, 111206, 111205, 111205, 111207, 111206, 111206, 111206, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111204, 111206, 111204, 111206, 111206, 111206, 111206, 111208, 111209, 111207, 111207, 111208, 111207, 111207, 111207, 111207, 111208, 111207, 111207, 111208, 111209, 111209, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111206, 111206, 111205, 111207, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111209, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111208, 111207, 111205, 111205, 111206, 111206, 111205, 111204, 111206, 111206, 111205, 111205, 111206, 111206, 111205, 111206, 111204, 111206, 111206, 111204, 111203, 111204, 111204, 111203, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111203, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111200, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111197, 111198, 111199, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111206, 111205, 111206, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111206, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111203, 111203, 111204, 111202, 111202, 111203, 111202, 111202, 111200, 111203, 111201, 111200, 111200, 111201, 111201, 111201, 111196, 111198, 111196, 111196, 111196, 111197, 111196, 111196, 111196, 111196, 111197, 111197, 111198, 111197, 111199, 111199, 111199, 111198, 111200, 111202, 111200, 111200, 111200, 111201, 111201, 111200, 111199, 111200, 111200, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111200, 111198, 111198, 111198, 111201, 111199, 111199, 111199, 111199, 111199, 111200, 111201, 111200, 111201, 111199, 111199, 111201, 111200, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111203, 111201, 111204, 111204, 111202, 111199, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111202, 111203, 111201, 111199, 111201, 111202, 111203, 111202, 111203, 111205, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111203, 111202, 111201, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111201, 111202, 111198, 111198, 111201, 111199, 111200, 111199, 111200, 111197, 111195, 111194, 111194, 111194, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111199, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111203, 111202, 111203, 111203, 111203, 111204, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111195, 111199, 111199, 111197, 111197, 111196, 111196, 111196, 111199, 111200, 111200, 111198, 111197, 111195, 111196, 111194, 111194, 111196, 111197, 111197, 111196, 111196, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111193, 111194, 111192, 111192, 111192, 111192, 111193, 111190, 111190, 111190, 111191, 111188, 111188, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111184, 111186, 111185, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111197, 111198, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111198, 111196, 111195, 111194, 111194, 111194, 111194, 111193, 111191, 111190, 111192, 111191, 111190, 111188, 111188, 111188, 111190, 111188, 111188, 111186, 111186, 111185, 111185, 111184, 111184, 111185, 111183, 111180, 111179, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111181, 111180, 111179, 111179, 111181, 111185, 111184, 111184, 111185, 111185, 111185, 111185, 111185, 111187, 111187, 111188, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111186, 111188, 111188, 111186, 111185, 111191, 111188, 111187, 111188, 111188, 111188, 111190, 111189, 111192, 111191, 111192, 111190, 111189, 111188, 111189, 111192, 111191, 111192, 111192, 111190, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111194, 111193, 111192, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111197, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111193, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111194, 111195, 111196, 111196, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111196, 111198, 111198, 111197, 111198, 111200, 111198, 111197, 111201, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111196, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111191, 111191, 111194, 111194, 111194, 111193, 111192, 111192, 111190, 111191, 111192, 111190, 111190, 111190, 111192, 111192, 111192, 111192, 111193, 111193, 111193, 111193, 111192, 111193, 111192, 111192, 111197, 111197, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111201, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111198, 111198, 111197, 111197, 111197, 111199, 111199, 111197, 111197, 111198, 111195, 111195, 111194, 111194, 111194, 111193, 111193, 111191, 111191, 111190, 111192, 111192, 111192, 111192, 111192, 111193, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111193, 111193, 111194]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64059, 64060, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64060, 64059, 64059, 64061, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64061, 64061, 64065, 64065, 64063, 64064, 64065, 64063, 64063, 64064, 64064, 64063, 64065, 64066, 64065, 64065, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64065, 64066, 64065, 64065, 64065, 64067, 64067, 64068, 64069, 64069, 64069, 64069, 64070, 64069, 64069, 64069, 64069, 64070, 64070, 64069, 64069, 64068, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64069, 64069, 64069, 64069, 64069, 64068, 64066, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64068, 64066, 64067, 64067, 64067, 64065, 64066, 64065, 64066, 64066, 64067, 64067, 64067, 64066, 64066, 64066, 64064, 64063, 64064, 64065, 64064, 64063, 64059, 64060, 64061, 64064, 64063, 64061, 64058, 64060, 64058, 64059, 64058, 64058, 64059, 64059, 64060, 64059, 64059, 64058, 64058, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64067, 64065, 64067, 64065, 64064, 64064, 64064, 64064, 64064, 64063, 64065, 64066, 64065, 64068, 64068, 64067, 64067, 64066, 64067, 64067, 64067, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64068, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64060, 64060, 64060, 64060, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64058, 64055, 64055, 64056, 64055, 64054, 64054, 64053, 64055, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64060, 64058, 64059, 64058, 64059, 64059, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64056, 64057, 64056, 64056, 64057, 64059, 64058, 64057, 64058, 64057, 64058, 64058, 64057, 64058, 64059, 64056, 64058, 64058, 64058, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64064, 64063, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64060, 64060, 64060, 64062, 64061, 64060, 64061, 64061, 64057, 64058, 64057, 64057, 64057, 64056, 64056, 64064, 64062, 64060, 64061, 64062, 64063, 64062, 64060, 64059, 64059, 64060, 64061, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64063, 64065, 64064, 64064, 64064, 64062, 64063, 64062, 64063, 64059, 64058, 64060, 64058, 64058, 64059, 64059, 64057, 64056, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64061, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64062, 64063, 64064, 64062, 64061, 64063, 64062, 64061, 64062, 64061, 64059, 64060, 64061, 64062, 64059, 64058, 64059, 64058, 64059, 64058, 64058, 64058, 64060, 64058, 64057, 64056, 64055, 64055, 64054, 64055, 64054, 64055, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64052, 64053, 64053, 64053, 64054, 64053, 64053, 64054, 64054, 64054, 64052, 64050, 64050, 64049, 64048, 64048, 64048, 64049, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64058, 64057, 64057, 64058, 64059, 64057, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64050, 64051, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64047, 64046, 64047, 64047, 64043, 64043, 64044, 64044, 64044, 64044, 64044, 64042, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64041, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64043, 64041, 64041, 64042, 64044, 64043, 64043, 64045, 64045, 64044, 64046, 64046, 64048, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64047, 64047, 64051, 64048, 64050, 64049, 64050, 64050, 64051, 64051, 64052, 64051, 64051, 64047, 64048, 64048, 64050, 64050, 64049, 64049, 64052, 64048, 64051, 64051, 64052, 64052, 64049, 64049, 64049, 64049, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64057, 64056, 64056, 64056, 64056, 64055, 64054, 64054, 64055, 64057, 64059, 64059, 64059, 64058, 64057, 64057, 64056, 64056, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64060, 64061, 64061, 64061, 64058, 64057, 64059, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64064, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64062, 64063, 64062, 64060, 64060, 64060, 64059, 64059, 64057, 64055, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64061, 64060, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64064, 64063, 64064, 64062, 64062, 64062, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64062, 64060, 64061, 64061, 64060, 64059, 64058, 64057, 64056, 64055, 64055, 64053, 64052, 64052, 64054, 64053, 64052, 64053, 64052, 64050, 64049, 64049, 64050, 64048, 64049, 64048, 64049, 64050, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64049, 64048, 64048, 64060, 64059, 64059, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64063, 64062, 64061, 64061, 64062, 64062, 64060, 64063, 64065, 64065, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64058, 64059, 64060, 64059, 64059, 64058, 64059, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64062, 64062, 64061, 64061, 64061, 64060, 64060, 64059, 64059, 64061, 64059, 64058, 64058, 64058, 64057, 64056, 64055, 64053, 64053, 64054, 64053, 64050, 64050, 64050, 64052, 64053, 64052, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64054, 64056, 64055, 64055]],
            'npts_real': [True, 191116800]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.01, 0.1])
        self.mom8_creator(img+'.residual', range_list=[-0.01, 0.1])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')

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
            
            self.save_dict_to_file(test_name,savedict, test_name+'_cas13317_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_cube_eph_pcwdT
#-------------------------------------------------#
    # Test 12c
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
            perchanweightdensity=True, gridder='mosaic', \
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
            perchanweightdensity=True, gridder='mosaic', \
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
        exp_im_stats = {'com_bmaj': [False, 0.9369738912503739],
            'com_bmin': [False, 0.7143172172695156],
            'com_pa': [False, -87.11471410651349],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 105415358.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.07762322574853897],
            'max_val_pos': [True, [319, 243, 0, 511]],
            'min_val': [False, -0.06524113565683365],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.01129179395717224],
            'rms_per_chan': [False, [0.01391552112260973, 0.011198396463193943, 0.010343933322022065, 0.010809416285053133, 0.010977990340729648, 0.012186306949614938, 0.012651767973330064, 0.014073053968191699, 0.011689241433724188, 0.01146726238232732, 0.010360057141555091, 0.010177265492376415, 0.012655764393299554, 0.010692999372625406, 0.010197460774929456, 0.009958039676528503, 0.010949263290257555, 0.012588522262073435, 0.008784171059800679, 0.011302257900173541, 0.0117787085568739, 0.01192900798109664, 0.011508006983753056, 0.012458716972351689, 0.01139023253347075, 0.010475022530977396, 0.010516796959728576, 0.01072399948798447, 0.012989196225412126, 0.012096989230959278, 0.010923440141567086, 0.010386804123588342, 0.01054613376091111, 0.009839170596410032, 0.011266597787145678, 0.011180208551680509, 0.01061306963262033, 0.011125305217235185, 0.010395052684136414, 0.010664896804158538, 0.010965564849589136, 0.011639596741143059, 0.011004326514788818, 0.01255767419671084, 0.011436760851161416, 0.010897629748286747, 0.010203862167219361, 0.009817739300640465, 0.01004740112200826, 0.010403404890970375, 0.011238606729570218, 0.011691532734109457, 0.010813975948494021, 0.010901475750564334, 0.010684322260151173, 0.01098304233220916, 0.010602280961452704, 0.01426648022699637, 0.01153735161583526, 0.009197121765880426, 0.009532115142635137, 0.011701320409887268, 0.01010561813988295, 0.010527492412262647, 0.01021247662894456, 0.010521304934798295, 0.010924949338486114, 0.011249584560179968, 0.010509240125843207, 0.00940799622619301, 0.012382289657307588, 0.011578648394620573, 0.010622373274679698, 0.009106163607450332, 0.010480814441017831, 0.012113439637458012, 0.012834856656003209, 0.012621576871742016, 0.012044942054955344, 0.011398526259345078, 0.010271918744797829, 0.010169992217735176, 0.010659667929461268, 0.011403352581790618, 0.010817467875967105, 0.010663181365935233, 0.010380907066466214, 0.010619922261589206, 0.009534057888651906, 0.010772919409697972, 0.010974885520884053, 0.011756044941656614, 0.010897897666976807, 0.009768199103077167, 0.010137322276836702, 0.011432570776799024, 0.011428269156318269, 0.011180581634084693, 0.01092732283810307, 0.009008971699236616, 0.011094094936623975, 0.011890030915220102, 0.012714469605527105, 0.011675270335864477, 0.01196888720235862, 0.012173302500335862, 0.01327120828641395, 0.011316360171815119, 0.012375625356419135, 0.011497782822133675, 0.010496183974832365, 0.010392383778172887, 0.011393160298779735, 0.010553412955603539, 0.010154148554427634, 0.0115074205310838, 0.010543375263794471, 0.009985402640483542, 0.011658219075259, 0.011066125474239217, 0.009687966637730104, 0.011334237941646202, 0.009534164684758761, 0.011328703723824355, 0.010821960402560202, 0.011310329641981127, 0.010093181571386394, 0.012434103784634004, 0.010312998652122903, 0.009126457106059851, 0.01062905143364573, 0.01027977788667286, 0.009604693031018581, 0.011203392952889574, 0.012175711122706682, 0.010736953769353097, 0.010829947520616402, 0.010670447051806092, 0.009945016686363722, 0.010503182953480837, 0.011054036426172278, 0.010259633675392177, 0.00976760023591767, 0.010093999970156078, 0.00889770849006833, 0.011307380889253623, 0.011076985873300735, 0.009431870550311201, 0.011292088007541667, 0.009663344112211005, 0.010460607340389555, 0.01097586126451534, 0.010500037171990353, 0.010341243533789152, 0.011445130488160834, 0.010719315777534339, 0.011120547173124168, 0.010911344933068655, 0.010240379864532464, 0.013484076272462099, 0.013296263039374088, 0.010737342604532131, 0.009603428692977508, 0.009645065481359157, 0.010262817299865709, 0.009819371619852675, 0.009965277405298345, 0.009055282338189717, 0.009586652931226633, 0.010481815501477035, 0.009289546552318636, 0.010365885329723873, 0.010103872846563756, 0.010600748972317442, 0.011905581081075839, 0.011413682167288993, 0.009545122943576003, 0.011747835122206797, 0.01198145117697354, 0.011710160960439841, 0.00996712452031176, 0.011369363212434041, 0.010813927221427926, 0.01212462234223, 0.010450404701814532, 0.010135846517441377, 0.011298001485090444, 0.009862074781774489, 0.011976019214079757, 0.010216024445230754, 0.010661861621396488, 0.010147551507939613, 0.011781167797481399, 0.009924771392013632, 0.01050320279901618, 0.011703743552674098, 0.010034410148838243, 0.011550480778495178, 0.010088713175294192, 0.011269806784639496, 0.011892306373653863, 0.01374597798653151, 0.012414091139083077, 0.010900898782020309, 0.010119232795671475, 0.009276780295921822, 0.009942327956983245, 0.010198851570012282, 0.011213533676259935, 0.010825870531098576, 0.011468561946841285, 0.011728910827083015, 0.011083476565502052, 0.010769001794555193, 0.00982655263594942, 0.011934504307180797, 0.01162826573193241, 0.010121880905910088, 0.010551722361175409, 0.01238726547867744, 0.012501632251717163, 0.010988762835532098, 0.011060427328795802, 0.012693781198233474, 0.011779446083693159, 0.009850364085612965, 0.011425461383082688, 0.011100265889725534, 0.011456678941830934, 0.01035581823666254, 0.0098960855434141, 0.012529067537023679, 0.011474107741700589, 0.010975887847209211, 0.009161356358172393, 0.011828144957072836, 0.011497648468475449, 0.011047302377789139, 0.010660621204848892, 0.013277136336152068, 0.011677221227331206, 0.00915460351794319, 0.009428985545839296, 0.010793182299579958, 0.011437180114632982, 0.010033333250168707, 0.010394632197880487, 0.010421144212822829, 0.011186910007496819, 0.011223423265429692, 0.012506807628565065, 0.010297749456457839, 0.009265900904257296, 0.011470258369577373, 0.011493480177761395, 0.012208647069534506, 0.01183949088063995, 0.01018258397518846, 0.010625984301686257, 0.010276584340383418, 0.010972108797526341, 0.01093000887346839, 0.012313947342454397, 0.010701440456944976, 0.010852022258113521, 0.011496175948316239, 0.011990847009793163, 0.009967159743546759, 0.011427261876768895, 0.012152300301912632, 0.010712436653209744, 0.010245123936610178, 0.010594556038532476, 0.011242300600707147, 0.010379633896139934, 0.010565109637854697, 0.011824473876155239, 0.01243244015341194, 0.012260167743711492, 0.01097903900057018, 0.009999147852798937, 0.010645926602965684, 0.010946820950604987, 0.008600435137409589, 0.010937189305399417, 0.011406495656560668, 0.01155984146154604, 0.011330002250252961, 0.011536676786056781, 0.012394571683428546, 0.01004088254138235, 0.01001458928029528, 0.010858140608079417, 0.010821618327605946, 0.011509482602924406, 0.009535222101165907, 0.010452064215699325, 0.010258114147910448, 0.011168867988531922, 0.01003341509689923, 0.01250002742156024, 0.01111406120168579, 0.011648474020214605, 0.01487042232672349, 0.014968014717286642, 0.014912044497928274, 0.010049657877993548, 0.011361474028869693, 0.010450736551600396, 0.012924759020456636, 0.01047269570677801, 0.01067798988986376, 0.011809299981466669, 0.012166834715303531, 0.011719974739243695, 0.010967355236127034, 0.011079358828134425, 0.013585751153939314, 0.0114397119232406, 0.011214231123430446, 0.011116776119994298, 0.014757481176270312, 0.011182988316413787, 0.010998045561777966, 0.010621592983351061, 0.012531498207540829, 0.013437953285084379, 0.011553247205728687, 0.010100374383503592, 0.011034995779197981, 0.010542581447883408, 0.012392044636652022, 0.011005896093877896, 0.011544313668244037, 0.011426996448672216, 0.01017993209564668, 0.01328749272057435, 0.012035981878614838, 0.010745168882430557, 0.013050322705029963, 0.010331443251391938, 0.012696532233140927, 0.009418567376640163, 0.010511268392686278, 0.011351815096869578, 0.01289801865864792, 0.010239046200121244, 0.011542222722144758, 0.009916383321385823, 0.010445314180810928, 0.011516844334473678, 0.0122694012799528, 0.011969362435971848, 0.011133288239148493, 0.011072187534401893, 0.01035150935637823, 0.01213657632753204, 0.01153504451961029, 0.01206594845317023, 0.011489966598334045, 0.010792915903095317, 0.011176893202432649, 0.01293138086791276, 0.01334243892339474, 0.010357050981237601, 0.010875617372884103, 0.010387946759901795, 0.011347179280819408, 0.012479472003997849, 0.01072950080237677, 0.011425176052804456, 0.01281799477616515, 0.011911859707582108, 0.010323188173236155, 0.011856571505421215, 0.011161691257014172, 0.010603779289994588, 0.01136142401900738, 0.011137892980180576, 0.009525740041064542, 0.010342320922698734, 0.01028672088730379, 0.012062763141598485, 0.010362860573444634, 0.011620182719982221, 0.010470234994187609, 0.011872522956010085, 0.01277907342972206, 0.010311842285642912, 0.012716217505088138, 0.011924208520129712, 0.01097481488044642, 0.010442252786586392, 0.011180528488909934, 0.01100030448990524, 0.009780102752006609, 0.009795227428078238, 0.01100552952939575, 0.011332424564401583, 0.014644063149767118, 0.011109062085278315, 0.011112287636348766, 0.011069664864464828, 0.012176442520039597, 0.010945180941767626, 0.009672923699368362, 0.010749971974514117, 0.01319640364213632, 0.014661636144175913, 0.010790524657865279, 0.011356739883718044, 0.011266133663662881, 0.012418745299879037, 0.011437575300394201, 0.011096409645036404, 0.012867706551998223, 0.012370757008851577, 0.010941169842382232, 0.011373191324738175, 0.010835402813871211, 0.010491924181292928, 0.008900259733189487, 0.010253540271952811, 0.010171644574121463, 0.01116430613954575, 0.011163407581305058, 0.011246658266871415, 0.012187951194384952, 0.010778404251978922, 0.01088427578140813, 0.011249911095120071, 0.011019141443474307, 0.01128020288188004, 0.011274481328144825, 0.010381529823489485, 0.01301497001284998, 0.01125642552039001, 0.013174756453017652, 0.011782450350197368, 0.011568772247716163, 0.011386671017979017, 0.013614480462608339, 0.012849886796338464, 0.01199692991905488, 0.010479203883010876, 0.011168245509840169, 0.011178316266981311, 0.01101013426925998, 0.010618226765243757, 0.011387461958685167, 0.013066281012805096, 0.01319855746012196, 0.011507316908488817, 0.012385309728060873, 0.011553268132355455, 0.010380201579075023, 0.010636976573218973, 0.009924749863242715, 0.012193647269262703, 0.01109868956111707, 0.012533886772384676, 0.011193706690753877, 0.01110849183005928, 0.011349176593576371, 0.009889609548710212, 0.01149872408324248, 0.010518655426434612, 0.010641531162000239, 0.011389799962659259, 0.011101021570044861, 0.012263541559989195, 0.010626222792129767, 0.011252242526876156, 0.010633394314831434, 0.010980607111152394, 0.011400695142225982, 0.011158377934131441, 0.012611511645868498, 0.01085047761891118, 0.010900271998382305, 0.011351638084714576, 0.011071194310012461, 0.011213832289207382, 0.015086062387120281, 0.011211797978263572, 0.011921695025218042, 0.011915360620387767, 0.011914439829529201, 0.011409733559730389, 0.011666742260575048, 0.011835896769492817, 0.010780372469583814, 0.011503611078264333, 0.012190906385723243, 0.010201794147592338, 0.011517199364141128, 0.011514250604873063, 0.011339320692666615, 0.013095772374127812, 0.012081215691636865, 0.01167870579906459, 0.013294500517003673, 0.010400274814914977, 0.012302996571727112, 0.011669734254711028, 0.011832658897866846, 0.010709610374694799, 0.012412932645115966, 0.011720702277626416, 0.011254746375650496, 0.012386139882724133, 0.01284507733588335, 0.012262346804734515, 0.010441843530507814, 0.012422117761770048, 0.01274881111696537, 0.012323210070529432, 0.013435331326741674, 0.012955939321448747, 0.011968721333003933, 0.010058823074007657, 0.012844731995447112, 0.012045054525003868, 0.01199587371291059, 0.012326286254450865, 0.009845084541710505, 0.011130486230741149, 0.012953538699511463, 0.011142113041188992, 0.011663717515592531, 0.010752740412162459, 0.011033546521962988, 0.010741806287994945, 0.011141711789612142, 0.011865394144902692, 0.011913938660192419, 0.012807388831490104, 0.010695983604251759, 0.010810973764566826, 0.01182388923107644, 0.011714101763328919, 0.011041126967030695, 0.011221645851562188, 0.01171570825069557, 0.011242986010453378, 0.011459259269154002, 0.010699052660184531, 0.012279268591655275, 0.01164298312122537, 0.011394119989051502, 0.01137792631891344, 0.011276947065139485, 0.01161109049111181, 0.01335231811318581, 0.012313972587404244, 0.012024389762422474, 0.010965879733210396, 0.011411890858713658, 0.012094763290635621, 0.010438183876828143, 0.012095057132373304, 0.013330003889317623, 0.012094958131995215, 0.011352345418891883, 0.012261559985070144, 0.010700129208378149, 0.01215578545310331, 0.013774233972922453, 0.012067107629934012, 0.011203174435945918, 0.011061969907609348, 0.011026618754079799, 0.011058650394077385, 0.011880966688155435, 0.012604234415151903, 0.010420076185615548, 0.011451647086075439, 0.011470257535051005, 0.010883601868666343, 0.01214935188633573, 0.01212620167771195, 0.01336416247709094, 0.01153798001556891, 0.012254619076193345, 0.013186080040769965, 0.011085265031964976, 0.010500618653592102, 0.012079916033932613, 0.012017841018753125, 0.013713873022785324, 0.013992960766276588, 0.012675763979996233, 0.012102177685112476, 0.013624649887855637, 0.012418253361231615, 0.012123056751937224, 0.010975262316582106, 0.012272153153510644, 0.011000964937257366, 0.011062129810947976, 0.01345843293334032, 0.013109262487988783, 0.010494003048493612, 0.011234712712262757, 0.011280596585901026, 0.011190662177615414, 0.012055008614906649, 0.010110022635582907, 0.010811958265203661, 0.011169763751830625, 0.013278508937467402, 0.012437322824421107, 0.014491960891414907, 0.01120831700747021, 0.010656437015920527, 0.012877356315977279, 0.010748434515587118, 0.012927102203282866, 0.012778635536341261, 0.012808319019161086, 0.010876451624590111, 0.012157599741825153, 0.01205640295753723, 0.012513737098713382, 0.012520130494081913, 0.01110127112907588, 0.012886360786114296, 0.011988596766017507, 0.010437173803573407, 0.012803178648571712, 0.013377319242311165, 0.010685703739657305, 0.012432317566056129, 0.011000685198200185, 0.010055378497488549, 0.013027217681346057, 0.011422103031368392, 0.010622424680790455, 0.01169517805919308, 0.011810925235764293, 0.011588023636216473, 0.012492483972232392, 0.01132464438396137, 0.012928041051797828, 0.011900866811310975, 0.013256253808160971, 0.011086625407640278, 0.010308502435335239, 0.012519620801843111, 0.011787486214853881, 0.010527404643714457, 0.01273101512153473, 0.010728110176605589, 0.012292865273394517, 0.011069875310313134, 0.011285102421407943, 0.012094626597776052, 0.011460696951244067, 0.012086398321361385, 0.01131602486284606, 0.011852138663998724, 0.012118574898453754, 0.013767773044714413, 0.012741061274580039, 0.013193541903863182, 0.011082459956534571, 0.011624220206139659, 0.011132162589114335, 0.012165508004319155, 0.011199868940797848, 0.009586378954027428, 0.010093136115865813, 0.011106932997820878, 0.011409897569441527, 0.014462850416500325, 0.010953485877208824, 0.011579041416334493, 0.011050698137006331, 0.010653133964167415, 0.014671508128760398, 0.011235644198452702, 0.010696159396036273, 0.011782375670566583, 0.009959957725922149, 0.012285114555421079, 0.011659234959911296, 0.011208617712754422, 0.010401120704625842, 0.012096745151680556, 0.01099494021556944, 0.009779837696181291, 0.012337531038455963, 0.010432047593846402, 0.010350095799656199, 0.01110397339839623, 0.012177411456587902, 0.01190496934778857, 0.01218601942111184, 0.01158032671626807, 0.011518971864208614, 0.011797897071130036, 0.011608606776569911, 0.011873994357711463, 0.012903089038410604, 0.011622742730973336, 0.0112721919482202, 0.012108449893431513, 0.010493702734973845, 0.013550570106305808, 0.012450754997505944, 0.01222264668068609, 0.0129317596138413, 0.010405203844920692, 0.010695101167482273, 0.011543484384249916, 0.01022213518020711, 0.011099459033023559, 0.00985564685260886, 0.01111296828173446, 0.01173923857546939, 0.010177841112496413, 0.011762594517028034, 0.012373724431677453, 0.010588725346733637, 0.010357376982988307, 0.01294663424453868, 0.010945446837699895, 0.012758181073638927, 0.010748098078665311, 0.010773550669874107, 0.01201689084364728, 0.011808714563029022, 0.009914440363372912, 0.009703219040000328, 0.009847603705416534, 0.010792168007821878, 0.010653653942806983, 0.012073601165067813, 0.010357284549865996, 0.010232505091602128, 0.013142644645674294, 0.011278097171992855, 0.010411772261604783, 0.009777926197370342, 0.010557876279618223, 0.01169776144855825, 0.01160698098869019, 0.010210840038921982, 0.011469961085604443, 0.013930172817395479, 0.010069506171264434, 0.010118403384111126, 0.00961120792435384, 0.010798192219377698, 0.01262696290898053, 0.011285766134244209, 0.01070043562065915, 0.011148171375914366, 0.01150132292664256, 0.011189540154382352, 0.010245285697351563, 0.014007603768924102, 0.011189645224871548, 0.012450067494127652, 0.014813454696038016, 0.010639868733716803, 0.010349839851105776, 0.012398776176217726, 0.012064597725173936, 0.01221994148940866, 0.011918970688819471, 0.010788843575465884, 0.010214967570739784, 0.011884441723651303, 0.011289296668807729, 0.010449051508741557, 0.011574329340125528, 0.010601289251650813, 0.010467688654061653, 0.010291968004713617, 0.012738053792305299, 0.011183438510395133, 0.012101326119942243, 0.01096368915023177, 0.01193806556534271, 0.011942403256438574, 0.011198287869649505, 0.010982761307659312, 0.011956896477592087, 0.012156717112390957, 0.012606505380360473, 0.009583606759183622, 0.010778435869455778, 0.00985013407324276, 0.010713836730758916, 0.01050919182877181, 0.009342819398409497, 0.01284523552408765, 0.010450381738152126, 0.010715817875274779, 0.00998505206154893, 0.011005264278406174, 0.010486670898006971, 0.010826178132963934, 0.01030504554026967, 0.011740874445339115, 0.010097661740348898, 0.010523899268981635, 0.01232497575976974, 0.011086387531309395, 0.010990482361696302, 0.009196767574539225, 0.011302748583633301, 0.009638381259110458, 0.010861814937452556, 0.011051924169079148, 0.011297121125404867, 0.011169905131133823, 0.012729343696629055, 0.01087665294554214, 0.009626357939898914, 0.009971019834659372, 0.011894472248965338, 0.010023041478220104, 0.010261091638596119, 0.00955730957538121, 0.010190965576011415, 0.012218669613733765, 0.013575782990596717, 0.011188780622430809, 0.01110923614440874, 0.012494384714873932, 0.01152161537362583, 0.01083610045980062, 0.011096755976328976, 0.011387984037470376, 0.011607497135092855, 0.010835448824275963, 0.011232855926961728, 0.010698634522633874, 0.01019811150140011, 0.012087869747125357, 0.01062800860232684, 0.012307010149479863, 0.010485314534694324, 0.010381638669462541, 0.011166642078924342, 0.009528517750261178, 0.009219966414349275, 0.009893937188609321, 0.01037325447479396, 0.011016877336031768, 0.010085811579706203, 0.010711228645828337, 0.011237730503854194, 0.01101878291301042, 0.009943060178147398, 0.010684603344342377, 0.010387585891597923, 0.009473761066755416, 0.010783345186202265, 0.012210889317997518, 0.010597302290037582, 0.009843258562160785, 0.009314400088976679, 0.009631041792153976, 0.011910715198681514, 0.011383564051389633, 0.012026537406552584, 0.012305143453734336, 0.010954884977935724, 0.010971173833196618, 0.01200996178031455, 0.011160723162695748, 0.009638285741310223, 0.00995479100027356, 0.00955886284852247, 0.010343874819538829, 0.010087788517892689, 0.011461551563193121, 0.011703353671492537, 0.010928971266745646, 0.011513332029277533, 0.010723860104163681, 0.01283378749693058, 0.012015340643280964, 0.011202781686609702, 0.009289346189445841, 0.010108304948268568, 0.009764801630256848, 0.009847497368910435, 0.010910573376823786, 0.010598551180669397, 0.011032076992544953, 0.009002280390016937, 0.011453624267679484, 0.012240622287238693, 0.010334991383745943, 0.010238032677970348, 0.0104552764406377, 0.009672457297922708, 0.009266758732895683, 0.010472260196509089, 0.010104350538404649, 0.009729403609383227, 0.010026901420582702, 0.009465560765216433, 0.010760793522209149, 0.010627689844371919, 0.010532986396577826, 0.01076598517468708, 0.013031900542068056, 0.010457689636671751, 0.009761123990214407, 0.010619408183739359, 0.010684152310898594, 0.010905519898812412, 0.01076248435002148, 0.013760013008656527, 0.011247740809603227, 0.009713306059121022, 0.01069577712199662, 0.00937853821358317, 0.011822502521068113, 0.010360099088483211, 0.012459330525831714, 0.010844430331161241, 0.00913156187088497, 0.011054266876339694, 0.010055114697346778, 0.009197520362908738, 0.012516623411902498, 0.012855716131619998, 0.010469981149541133, 0.0101500166006628, 0.009690469253029847, 0.010442118528290653, 0.010269280672524347, 0.009827210980388081, 0.010001809924929807, 0.01071932054204389, 0.011169647425339812, 0.01062816738366297, 0.011952248179244092, 0.012375994994089994, 0.010226028064128328, 0.009925983558943821, 0.011261621413049205]],
            'im_sum': [False, 9577.760838691873],
            'regn_sum': [False, 26.161059313285023],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.012518767238255986, 0.011504943829556804, 0.011792719193721971, 0.011843068051059186, 0.01141142197484636, 0.012311975745280002, 0.011650528959617571]],
            'profile': [False, 0.022807379768021347]}

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
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'mask_pix': [False, 9366],
            'mask_regns': [True, 36],
            'npts_real': [True, 191116800]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415358.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000003278255463],
            'im_rms': [False, 0.6370609915857086],
            'npts_0.2': [False, [111206, 111210, 111208, 111207, 111205, 111204, 111201, 111204, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111208, 111206, 111206, 111207, 111205, 111205, 111205, 111206, 111204, 111205, 111206, 111205, 111205, 111205, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111203, 111206, 111204, 111206, 111206, 111206, 111206, 111207, 111208, 111207, 111207, 111207, 111206, 111207, 111207, 111206, 111207, 111207, 111207, 111207, 111209, 111208, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111207, 111206, 111205, 111207, 111205, 111206, 111206, 111206, 111205, 111205, 111207, 111208, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111205, 111206, 111206, 111206, 111205, 111204, 111205, 111205, 111205, 111205, 111206, 111205, 111205, 111206, 111204, 111206, 111206, 111203, 111203, 111204, 111204, 111204, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111203, 111202, 111202, 111202, 111202, 111203, 111203, 111201, 111202, 111202, 111203, 111203, 111202, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111196, 111198, 111198, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111207, 111205, 111205, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111206, 111204, 111205, 111206, 111206, 111205, 111206, 111206, 111205, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111206, 111203, 111203, 111203, 111201, 111202, 111203, 111201, 111203, 111200, 111202, 111201, 111200, 111200, 111201, 111201, 111201, 111197, 111198, 111197, 111196, 111197, 111196, 111196, 111196, 111196, 111196, 111197, 111197, 111197, 111197, 111198, 111199, 111198, 111197, 111200, 111202, 111200, 111200, 111200, 111200, 111200, 111200, 111199, 111200, 111199, 111197, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111200, 111197, 111198, 111198, 111201, 111199, 111199, 111199, 111198, 111199, 111199, 111200, 111200, 111201, 111198, 111199, 111201, 111199, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111203, 111201, 111203, 111204, 111203, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111202, 111203, 111201, 111199, 111201, 111201, 111203, 111202, 111202, 111204, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111203, 111201, 111201, 111202, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111200, 111201, 111199, 111199, 111201, 111199, 111199, 111200, 111200, 111197, 111195, 111194, 111194, 111195, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111198, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111196, 111198, 111199, 111197, 111196, 111197, 111196, 111196, 111199, 111200, 111199, 111197, 111197, 111195, 111195, 111194, 111194, 111196, 111198, 111196, 111196, 111195, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111192, 111193, 111192, 111192, 111192, 111191, 111193, 111190, 111190, 111190, 111191, 111188, 111189, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111185, 111186, 111184, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111198, 111196, 111195, 111194, 111194, 111194, 111193, 111193, 111191, 111190, 111191, 111191, 111191, 111188, 111188, 111187, 111190, 111188, 111188, 111187, 111186, 111185, 111185, 111184, 111184, 111184, 111184, 111181, 111181, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111179, 111181, 111181, 111180, 111178, 111181, 111181, 111185, 111184, 111184, 111185, 111185, 111184, 111184, 111185, 111188, 111187, 111188, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111186, 111188, 111187, 111186, 111186, 111191, 111186, 111187, 111188, 111188, 111188, 111190, 111190, 111192, 111191, 111192, 111190, 111190, 111188, 111189, 111191, 111191, 111192, 111192, 111189, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111193, 111192, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111196, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111193, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111193, 111195, 111195, 111195, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111197, 111198, 111196, 111198, 111198, 111198, 111198, 111200, 111198, 111197, 111201, 111200, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111197, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111197, 111199, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111192, 111191, 111194, 111194, 111193, 111193, 111191, 111192, 111190, 111190, 111192, 111190, 111190, 111190, 111192, 111192, 111191, 111192, 111192, 111193, 111193, 111192, 111192, 111192, 111191, 111191, 111197, 111196, 111197, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111202, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111199, 111197, 111197, 111198, 111195, 111195, 111194, 111193, 111194, 111193, 111192, 111191, 111191, 111190, 111191, 111192, 111192, 111193, 111192, 111192, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111194, 111193, 111193]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64058, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64060, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64062, 64061, 64065, 64065, 64063, 64062, 64064, 64063, 64061, 64063, 64064, 64063, 64065, 64065, 64065, 64064, 64064, 64064, 64064, 64064, 64064, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64064, 64066, 64064, 64065, 64066, 64067, 64067, 64068, 64068, 64069, 64070, 64069, 64070, 64069, 64070, 64069, 64069, 64070, 64070, 64069, 64069, 64067, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64067, 64069, 64069, 64069, 64068, 64068, 64065, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64066, 64066, 64068, 64066, 64065, 64065, 64065, 64065, 64066, 64066, 64067, 64067, 64066, 64066, 64066, 64066, 64065, 64063, 64064, 64065, 64064, 64062, 64059, 64060, 64061, 64064, 64063, 64060, 64058, 64060, 64058, 64059, 64058, 64058, 64058, 64059, 64059, 64059, 64059, 64058, 64057, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64066, 64065, 64067, 64065, 64064, 64064, 64064, 64063, 64064, 64063, 64064, 64066, 64064, 64068, 64068, 64066, 64066, 64066, 64067, 64067, 64067, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64067, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64059, 64060, 64060, 64059, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64057, 64055, 64056, 64056, 64055, 64054, 64054, 64053, 64054, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64055, 64056, 64055, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64058, 64058, 64057, 64058, 64059, 64056, 64057, 64058, 64057, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64059, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64064, 64062, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64061, 64060, 64060, 64061, 64061, 64060, 64060, 64061, 64057, 64058, 64057, 64057, 64056, 64055, 64056, 64064, 64062, 64059, 64061, 64062, 64062, 64061, 64060, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64061, 64062, 64062, 64062, 64063, 64064, 64064, 64064, 64063, 64061, 64062, 64062, 64063, 64059, 64058, 64059, 64058, 64058, 64059, 64059, 64057, 64056, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64054, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64060, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64063, 64063, 64064, 64062, 64062, 64062, 64062, 64061, 64061, 64061, 64059, 64060, 64061, 64060, 64059, 64057, 64058, 64058, 64059, 64058, 64058, 64058, 64060, 64058, 64057, 64055, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64054, 64053, 64051, 64050, 64050, 64050, 64048, 64048, 64048, 64050, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64058, 64057, 64057, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64054, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64050, 64051, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64047, 64046, 64046, 64046, 64044, 64044, 64043, 64044, 64044, 64044, 64044, 64042, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64042, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64043, 64040, 64041, 64042, 64044, 64044, 64043, 64044, 64045, 64044, 64045, 64044, 64048, 64048, 64047, 64047, 64048, 64048, 64046, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64047, 64047, 64051, 64049, 64049, 64049, 64051, 64050, 64051, 64051, 64052, 64051, 64051, 64047, 64048, 64047, 64049, 64050, 64049, 64049, 64051, 64049, 64051, 64051, 64051, 64051, 64049, 64049, 64048, 64048, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64055, 64055, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64056, 64056, 64055, 64054, 64055, 64055, 64057, 64058, 64058, 64058, 64058, 64057, 64057, 64056, 64057, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64059, 64060, 64061, 64061, 64058, 64058, 64058, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64063, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64064, 64064, 64062, 64063, 64062, 64060, 64060, 64060, 64059, 64058, 64057, 64054, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64060, 64060, 64063, 64064, 64062, 64063, 64064, 64063, 64063, 64062, 64064, 64063, 64063, 64062, 64062, 64062, 64063, 64063, 64061, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64061, 64061, 64061, 64061, 64060, 64059, 64058, 64057, 64056, 64055, 64055, 64053, 64052, 64052, 64054, 64052, 64052, 64053, 64052, 64050, 64049, 64048, 64049, 64048, 64049, 64048, 64049, 64049, 64050, 64050, 64049, 64049, 64049, 64050, 64049, 64048, 64048, 64048, 64059, 64059, 64059, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64063, 64062, 64061, 64061, 64062, 64062, 64060, 64064, 64064, 64065, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64060, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64061, 64059, 64061, 64060, 64060, 64061, 64061, 64062, 64062, 64060, 64061, 64061, 64060, 64060, 64059, 64059, 64061, 64058, 64058, 64058, 64058, 64057, 64056, 64055, 64054, 64053, 64054, 64053, 64050, 64050, 64051, 64052, 64053, 64053, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64054, 64056, 64056, 64055]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1042804076994885, 37.06023635817038, 36.75839002675014]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.8801035135706],
            'fit_pix': [False, [240.74284851156537, 210.14459456183218]]}


        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.055055033415555954],
            'min_val_pos': [True, [250, 204, 0, 947]],
            'im_rms': [False, 0.012116211767164618],
            'im_sum': [False, 75.47834240821368],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.8863266542103199, 1.0827418437842484, 0.8432045654114863]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462000557837],
            'fit_pix_0': [False, [239.9983688269633, 209.9973327835061]],
            'fit_1': [False, [0.886964262720369, 1.0819246776219098, 0.8430340347754101]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.8801035135706],
            'fit_pix_1': [False, [239.99852220850903, 209.99726355342318]],
            'fit_2': [False, [0.8875771221733474, 1.082587354736208, 0.8417507627592157]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99558702156276],
            'fit_pix_2': [False, [239.99845820502563, 209.9971807078754]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415358.0],
            'freq_bin': [1e-10, 244151.1796875],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.06217697635293007],
            'max_val_pos': [True, [320, 214, 0, 617]],
            'min_val': [False, -0.06524113565683365],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.01128499026981346],
            'im_sum': [False, 9230.282666749963],
            'regn_sum': [False, -14.764186072577559],
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
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.03260615095496178],
            'max_val_pos': [True, [318, 238, 0, 589]],
            'min_val': [False, -0.003487145993858576],
            'min_val_pos': [True, [323, 245, 0, 589]],
            'im_rms': [False, 2.2776357806895966e-05],
            'im_sum': [False, 8.985301391105168],
            'regn_sum': [False, 1.1452979089226574],
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
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 45549.24609375],
            'max_val_pos': [True, [0, 0, 0, 283]],
            'min_val': [False, 45503.1953125],
            'min_val_pos': [True, [0, 0, 0, 856]],
            'im_rms': [False, 45529.93669806667],
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
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.3204345107078552],
            'max_val_pos': [True, [240, 210, 0, 595]],
            'min_val': [False, 7.415285654133186e-05],
            'im_rms': [False, 0.12064659005311969],
            'im_sum': [False, 13934226.28923168],
            'npts_0.2': [False, [111206, 111210, 111208, 111207, 111205, 111204, 111201, 111204, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111208, 111206, 111206, 111207, 111205, 111205, 111205, 111206, 111204, 111205, 111206, 111205, 111205, 111205, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111203, 111206, 111204, 111206, 111206, 111206, 111206, 111207, 111208, 111207, 111207, 111207, 111206, 111207, 111207, 111206, 111207, 111207, 111207, 111207, 111209, 111208, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111207, 111206, 111205, 111207, 111205, 111206, 111206, 111206, 111205, 111205, 111207, 111208, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111205, 111206, 111206, 111206, 111205, 111204, 111205, 111205, 111205, 111205, 111206, 111205, 111205, 111206, 111204, 111206, 111206, 111203, 111203, 111204, 111204, 111204, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111203, 111202, 111202, 111202, 111202, 111203, 111203, 111201, 111202, 111202, 111203, 111203, 111202, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111196, 111198, 111198, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111207, 111205, 111205, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111206, 111204, 111205, 111206, 111206, 111205, 111206, 111206, 111205, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111206, 111203, 111203, 111203, 111201, 111202, 111203, 111201, 111203, 111200, 111202, 111201, 111200, 111200, 111201, 111201, 111201, 111197, 111198, 111197, 111196, 111197, 111196, 111196, 111196, 111196, 111196, 111197, 111197, 111197, 111197, 111198, 111199, 111198, 111197, 111200, 111202, 111200, 111200, 111200, 111200, 111200, 111200, 111199, 111200, 111199, 111197, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111200, 111197, 111198, 111198, 111201, 111199, 111199, 111199, 111198, 111199, 111199, 111200, 111200, 111201, 111198, 111199, 111201, 111199, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111203, 111201, 111203, 111204, 111203, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111202, 111203, 111201, 111199, 111201, 111201, 111203, 111202, 111202, 111204, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111203, 111201, 111201, 111202, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111200, 111201, 111199, 111199, 111201, 111199, 111199, 111200, 111200, 111197, 111195, 111194, 111194, 111195, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111198, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111196, 111198, 111199, 111197, 111196, 111197, 111196, 111196, 111199, 111200, 111199, 111197, 111197, 111195, 111195, 111194, 111194, 111196, 111198, 111196, 111196, 111195, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111192, 111193, 111192, 111192, 111192, 111191, 111193, 111190, 111190, 111190, 111191, 111188, 111189, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111185, 111186, 111184, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111198, 111196, 111195, 111194, 111194, 111194, 111193, 111193, 111191, 111190, 111191, 111191, 111191, 111188, 111188, 111187, 111190, 111188, 111188, 111187, 111186, 111185, 111185, 111184, 111184, 111184, 111184, 111181, 111181, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111179, 111181, 111181, 111180, 111178, 111181, 111181, 111185, 111184, 111184, 111185, 111185, 111184, 111184, 111185, 111188, 111187, 111188, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111186, 111188, 111187, 111186, 111186, 111191, 111186, 111187, 111188, 111188, 111188, 111190, 111190, 111192, 111191, 111192, 111190, 111190, 111188, 111189, 111191, 111191, 111192, 111192, 111189, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111193, 111192, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111196, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111193, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111193, 111195, 111195, 111195, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111197, 111198, 111196, 111198, 111198, 111198, 111198, 111200, 111198, 111197, 111201, 111200, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111197, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111197, 111199, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111192, 111191, 111194, 111194, 111193, 111193, 111191, 111192, 111190, 111190, 111192, 111190, 111190, 111190, 111192, 111192, 111191, 111192, 111192, 111193, 111193, 111192, 111192, 111192, 111191, 111191, 111197, 111196, 111197, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111202, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111199, 111197, 111197, 111198, 111195, 111195, 111194, 111193, 111194, 111193, 111192, 111191, 111191, 111190, 111191, 111192, 111192, 111193, 111192, 111192, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111194, 111193, 111193]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64058, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64060, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64062, 64061, 64065, 64065, 64063, 64062, 64064, 64063, 64061, 64063, 64064, 64063, 64065, 64065, 64065, 64064, 64064, 64064, 64064, 64064, 64064, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64064, 64066, 64064, 64065, 64066, 64067, 64067, 64068, 64068, 64069, 64070, 64069, 64070, 64069, 64070, 64069, 64069, 64070, 64070, 64069, 64069, 64067, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64067, 64069, 64069, 64069, 64068, 64068, 64065, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64066, 64066, 64068, 64066, 64065, 64065, 64065, 64065, 64066, 64066, 64067, 64067, 64066, 64066, 64066, 64066, 64065, 64063, 64064, 64065, 64064, 64062, 64059, 64060, 64061, 64064, 64063, 64060, 64058, 64060, 64058, 64059, 64058, 64058, 64058, 64059, 64059, 64059, 64059, 64058, 64057, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64066, 64065, 64067, 64065, 64064, 64064, 64064, 64063, 64064, 64063, 64064, 64066, 64064, 64068, 64068, 64066, 64066, 64066, 64067, 64067, 64067, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64067, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64059, 64060, 64060, 64059, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64057, 64055, 64056, 64056, 64055, 64054, 64054, 64053, 64054, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64055, 64056, 64055, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64058, 64058, 64057, 64058, 64059, 64056, 64057, 64058, 64057, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64059, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64064, 64062, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64061, 64060, 64060, 64061, 64061, 64060, 64060, 64061, 64057, 64058, 64057, 64057, 64056, 64055, 64056, 64064, 64062, 64059, 64061, 64062, 64062, 64061, 64060, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64061, 64062, 64062, 64062, 64063, 64064, 64064, 64064, 64063, 64061, 64062, 64062, 64063, 64059, 64058, 64059, 64058, 64058, 64059, 64059, 64057, 64056, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64054, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64060, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64063, 64063, 64064, 64062, 64062, 64062, 64062, 64061, 64061, 64061, 64059, 64060, 64061, 64060, 64059, 64057, 64058, 64058, 64059, 64058, 64058, 64058, 64060, 64058, 64057, 64055, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64054, 64053, 64051, 64050, 64050, 64050, 64048, 64048, 64048, 64050, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64058, 64057, 64057, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64054, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64050, 64051, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64047, 64046, 64046, 64046, 64044, 64044, 64043, 64044, 64044, 64044, 64044, 64042, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64042, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64043, 64040, 64041, 64042, 64044, 64044, 64043, 64044, 64045, 64044, 64045, 64044, 64048, 64048, 64047, 64047, 64048, 64048, 64046, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64047, 64047, 64051, 64049, 64049, 64049, 64051, 64050, 64051, 64051, 64052, 64051, 64051, 64047, 64048, 64047, 64049, 64050, 64049, 64049, 64051, 64049, 64051, 64051, 64051, 64051, 64049, 64049, 64048, 64048, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64055, 64055, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64056, 64056, 64055, 64054, 64055, 64055, 64057, 64058, 64058, 64058, 64058, 64057, 64057, 64056, 64057, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64059, 64060, 64061, 64061, 64058, 64058, 64058, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64063, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64064, 64064, 64062, 64063, 64062, 64060, 64060, 64060, 64059, 64058, 64057, 64054, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64060, 64060, 64063, 64064, 64062, 64063, 64064, 64063, 64063, 64062, 64064, 64063, 64063, 64062, 64062, 64062, 64063, 64063, 64061, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64061, 64061, 64061, 64061, 64060, 64059, 64058, 64057, 64056, 64055, 64055, 64053, 64052, 64052, 64054, 64052, 64052, 64053, 64052, 64050, 64049, 64048, 64049, 64048, 64049, 64048, 64049, 64049, 64050, 64050, 64049, 64049, 64049, 64050, 64049, 64048, 64048, 64048, 64059, 64059, 64059, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64063, 64062, 64061, 64061, 64062, 64062, 64060, 64064, 64064, 64065, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64060, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64061, 64059, 64061, 64060, 64060, 64061, 64061, 64062, 64062, 64060, 64061, 64061, 64060, 64060, 64059, 64059, 64061, 64058, 64058, 64058, 64058, 64057, 64056, 64055, 64054, 64053, 64054, 64053, 64050, 64050, 64051, 64052, 64053, 64053, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64054, 64056, 64056, 64055]],
            'npts_real': [True, 191116800]}

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            ["2017.1.00750.T_tclean_exe1.ms", \
             "2017.1.00750.T_tclean_exe2.ms"])

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.01, 0.1])
        self.mom8_creator(img+'.residual', range_list=[-0.01, 0.1])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

        test_dict[test_name]['images'].append(img+'.image.profile.png')

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
            
            self.save_dict_to_file(test_name,savedict, test_name+'_cas13317_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_mosaic_cube_eph_briggsbwtaper
#-------------------------------------------------#
    # Test 13
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
            perchanweightdensity=False, gridder='mosaic',  \
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
            perchanweightdensity=False, gridder='mosaic',  \
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
    # Test 14
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
            gridder='mosaic',  mosweight=True, \
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
            mosweight=True, usepointing=False, pblimit=0.2, \
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

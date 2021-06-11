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
#from casatestutils import almastktestutils
import sys
sys.path.append('/export/home/murasame/casa/casa6/casatests/stakeholder/')
import almastktestutils


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
savemetricdict=False

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
        self.expdict_jsonfile = data_path+'almastk_exp_dicts.json'

    def tearDown(self):
        #print("TEST_DICT=",test_dict)
        generate_weblog("tclean_ALMA_pipeline",test_dict)
        print("Closing ia tool")
        self._myia.done()
        """ don't delete it all """
#        self.delData()

    def getExpdicts(self, testname):
        ''' read the json file containung exp_dicts (fiducial metric values)
        for a specific test '''
        self.exp_dicts=almastktestutils.read_testcase_expdicts(self.expdict_jsonfile, testname)
         
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

     
    def save_dict_to_file(self, topkey, indict, outfilename, appendversion=True, outformat='JSON'):
        """ function that will save input Python dictionaries to a JSON file (default)
            or pickle file. topkey is will be added as a top key for output (nested) dictionary 
            and indict is stored under the key.  
            Create a separate file with outfilename if appendversion=True casa version (based on
            casatasks version) will be appended to the output file name.
        """
        try: 
            import casatasks as __casatasks
            casaversion = __casatasks.version_string()
            del __casatasks
        except:
            casaversion = ''
            
        if casaversion !='':
            casaversion = '_' + casaversion
        if type(indict) != dict:
            print("indict is not a dict. Saved file may not be in correct format") 
        nestedDict={}
        nestedDict[topkey]=indict 
        print("Saving %s dictionaries", len(indict))
        if outformat == 'pickle':
            # writing to pickle: note if writing this way (without protocol=2) 
            # in casa6 and read in casa5 it will fail 
            with open(outfilename+casaversion+'.pickle', 'wb') as outf:
                pickle.dump(nestedDict, outf)
        elif outformat== 'JSON':
            with open(outfilename+casaversion+'.json', 'w') as outf:
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
            if 'taskcall' in output[testname] and len(output[testname]['taskcall'])==3: 
                if parallel:
                    # 0,1,2th in the list are used pop last one
                    output[testname]['taskcall'].pop()
                else:
                    output[testname]['taskcall'].pop(1)
                    #output[testname]['taskcall'].pop(1)
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

        test_name = self._testMethodName 
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')
        self.getExpdicts(test_name)

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
        tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
            '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
            outframe='LSRK', perchanweightdensity=False, \
            usepointing=False, pblimit=0.2, nsigma=0.0, \
            gridder='standard',  mosweight=False, \
            deconvolver='hogbom', restoringbeam='common', restoration=True, pbcor=True, \
            weighting='briggs', robust=0.5, npixels=0, niter=20000, \
            threshold='0.354Jy', interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            calcres=False, calcpsf=False, savemodel='none', \
            parallel=self.parallel, verbose=True)

        # retrieve per-channel beam statistics
        bmin_dict, bmaj_dict, pa_dict = \
            self.cube_beam_stats(image=img+'.psf')

        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report(test_standard_cube)
        im_stats_dict = self.image_stats(image=img+'.image', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube.exp_im_stats
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube.exp_mask_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(image=img+'.psf', fit_region = \
            'ellipse[[11.47648725deg, -73.25812003deg], [8.0291arcsec, 6.8080arcsec], 90.00000000deg]')

        # test_standard_cube.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        # test_standard_cube.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube.exp_bmin_dict
            exp_bmin_dict = self.exp_dicts['exp_bmin_dict']
            # test_standard_cube.exp_bmaj_dict
            exp_bmaj_dict = self.exp_dicts['exp_bmaj_dict']
            # test_standard_cube.exp_pa_dict
            exp_pa_dict = self.exp_dicts['exp_pa_dict']


            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed=self.filter_report(report)
        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, 'test_standard_cube', self.parallel)

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
            
            savedict['bmin_dict']=bmin_dict
            savedict['bmaj_dict']=bmaj_dict
            savedict['pa_dict']=pa_dict

            self.save_dict_to_file(test_name,savedict, test_name+'_cur_stats')

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
        self.getExpdicts(test_name)

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
        tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
            '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
            outframe='LSRK', perchanweightdensity=True, \
            usepointing=False, pblimit=0.2, nsigma=0.0, \
            gridder='standard', mosweight=False, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggs', robust=0.5, npixels=0, niter=20000, \
            threshold='0.354Jy', interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            calcres=False, calcpsf=False, savemodel='none', \
            parallel=self.parallel, verbose=True)

        # retrieve per-channel beam statistics 
        bmin_dict, bmaj_dict, pa_dict = \
            self.cube_beam_stats(image=img+'.psf')


        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report(test_standard_cube_pcwdT)
        im_stats_dict = self.image_stats(image=img+'.image', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_im_stats
        exp_im_stats = self.exp_dicts['exp_im_stats']


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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_mask_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(image=img+'.psf', fit_region = \
            'ellipse[[11.47648725deg, -73.25812003deg], [8.0291arcsec, 6.8080arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']


        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_pcwdT.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_pcwdT.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        # test_standard_cube_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube_pcwdT.exp_bmin_dict
            exp_bmin_dict = self.exp_dicts['exp_bmin_dict']
            # test_standard_cube_pcwdT.exp_bmaj_dict
            exp_bmaj_dict = self.exp_dicts['exp_bmaj_dict']
            # test_standard_cube_pcwdT.exp_pa_dict
            exp_pa_dict = self.exp_dicts['exp_pa_dict']

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

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')


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
        self.getExpdicts(test_name)

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
        tclean(vis=self.msfile, imagename=file_name+'1', field='1', \
            spw=['0'], imsize=[80, 80], antenna=['0,1,2,3,4,5,6,7,8'], \
            scan=['8,12,16'], intent='OBSERVE_TARGET#ON_SOURCE', \
            datacolumn='data', cell=['1.1arcsec'], phasecenter='ICRS '
            '00:45:54.3836 -073.15.29.413', stokes='I', specmode='cube', \
            nchan=508, start='220.2526743594GHz', width='0.2441741MHz',\
            outframe='LSRK', perchanweightdensity=True, \
            usepointing=False, pblimit=0.2, nsigma=0.0, \
            gridder='standard', mosweight=False, \
            deconvolver='hogbom', restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggsbwtaper', robust=0.5, npixels=0, niter=20000, \
            threshold='0.354Jy', interactive=0, usemask='auto'
            '-multithresh', sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.08, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            calcres=False, calcpsf=False, savemodel='none', \
            parallel=self.parallel, verbose=True)

        # retrieve per-channel beam statistics 
        bmin_dict, bmaj_dict, pa_dict = \
            self.cube_beam_stats(image=img+'.psf')


        report0 = th.checkall(imgexist = self.image_list(img, 'standard'))

        # .image report(test_standard_cube_briggsbwtaper)
        im_stats_dict = self.image_stats(image=img+'.image', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_im_stats
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(image=img+'.pb', fit_region = \
            'ellipse[[11.47659846deg, -73.25817055deg], [23.1086arcsec, 23.0957arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(image=img+'.psf', fit_region = \
            'ellipse[[11.47648725deg, -73.25812003deg], [8.0291arcsec, 6.8080arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(image=img+'.residual', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]')

        # test_standard_cube_briggsbwtaper.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(image=img+'.model', fit_region = \
            'ellipse[[11.47881897deg, -73.25881015deg], [9.0414arcsec, 8.4854arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_briggsbwtaper.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(image=img+'.sumwt')

        # test_standard_cube_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        if self.parallel:
            # test_standard_cube_briggsbwtaper.exp_bmin_dict
            exp_bmin_dict = {'*0': 6.202030658721924,'*1': 6.202027320861816,'*10': 6.2019429206848145,'*100': 6.201361179351807,'*101': 6.201338768005371,'*102': 6.2013325691223145,'*103': 6.2013258934021,'*104': 6.201310634613037,'*105': 6.201310634613037,'*106': 6.201282978057861,'*107': 6.201282978057861,'*108': 6.201256275177002,'*109': 6.201256275177002,'*11': 6.201923370361328,'*110': 6.201256275177002,'*111': 6.201256275177002,'*112': 6.201256275177002,'*113': 6.201256275177002,'*114': 6.201256275177002,'*115': 6.201256275177002,'*116': 6.201256275177002,'*117': 6.201247692108154,'*118': 6.201218605041504,'*119': 6.201218605041504,'*12': 6.20189905166626,'*120': 6.201183319091797,'*121': 6.201183319091797,'*122': 6.20115852355957,'*123': 6.20112943649292,'*124': 6.201098918914795,'*125': 6.201138496398926,'*126': 6.201130390167236,'*127': 6.201130390167236,'*128': 6.201130390167236,'*129': 6.2011213302612305,'*13': 6.20189905166626,'*130': 6.2011213302612305,'*131': 6.2011213302612305,'*132': 6.2011213302612305,'*133': 6.201091766357422,'*134': 6.201091766357422,'*135': 6.201021671295166,'*136': 6.201021671295166,'*137': 6.201016426086426,'*138': 6.200996398925781,'*139': 6.200996398925781,'*14': 6.20189905166626,'*140': 6.2009782791137695,'*141': 6.2009782791137695,'*142': 6.2009782791137695,'*143': 6.2009782791137695,'*144': 6.2009782791137695,'*145': 6.2009782791137695,'*146': 6.200955867767334,'*147': 6.200955867767334,'*148': 6.200918674468994,'*149': 6.20089864730835,'*15': 6.20189905166626,'*150': 6.20089864730835,'*151': 6.20089864730835,'*152': 6.200883388519287,'*153': 6.200883388519287,'*154': 6.200864791870117,'*155': 6.200864791870117,'*156': 6.200839996337891,'*157': 6.200839996337891,'*158': 6.200839996337891,'*159': 6.200839996337891,'*16': 6.201883316040039,'*160': 6.200839996337891,'*161': 6.200765132904053,'*162': 6.200765132904053,'*163': 6.200747013092041,'*164': 6.200730323791504,'*165': 6.200721263885498,'*166': 6.200707912445068,'*167': 6.200707912445068,'*168': 6.200702667236328,'*169': 6.200702667236328,'*17': 6.201852798461914,'*170': 6.200702667236328,'*171': 6.200702667236328,'*172': 6.200685501098633,'*173': 6.200675964355469,'*174': 6.200669765472412,'*175': 6.20062780380249,'*176': 6.200601577758789,'*177': 6.200589179992676,'*178': 6.200563430786133,'*179': 6.200563430786133,'*18': 6.201840877532959,'*180': 6.200540065765381,'*181': 6.200540065765381,'*182': 6.200540065765381,'*183': 6.200534820556641,'*184': 6.200538158416748,'*185': 6.200538158416748,'*186': 6.200523376464844,'*187': 6.200523376464844,'*188': 6.200523376464844,'*189': 6.200523376464844,'*19': 6.201831340789795,'*190': 6.200523376464844,'*191': 6.200523376464844,'*192': 6.2005181312561035,'*193': 6.2005181312561035,'*194': 6.2005181312561035,'*195': 6.2005181312561035,'*196': 6.2005181312561035,'*197': 6.2005181312561035,'*198': 6.200517177581787,'*199': 6.200491905212402,'*2': 6.202004432678223,'*20': 6.201831340789795,'*200': 6.200491905212402,'*201': 6.200491905212402,'*202': 6.200467109680176,'*203': 6.200464248657227,'*204': 6.200451374053955,'*205': 6.200451374053955,'*206': 6.200451374053955,'*207': 6.200451374053955,'*208': 6.200451374053955,'*209': 6.200451374053955,'*21': 6.201831340789795,'*210': 6.200445175170898,'*211': 6.200445175170898,'*212': 6.200445175170898,'*213': 6.200445175170898,'*214': 6.200445175170898,'*215': 6.200445175170898,'*216': 6.200429439544678,'*217': 6.200429439544678,'*218': 6.2004194259643555,'*219': 6.200418949127197,'*22': 6.201831340789795,'*220': 6.200412273406982,'*221': 6.200412273406982,'*222': 6.200412273406982,'*223': 6.200412273406982,'*224': 6.20039176940918,'*225': 6.20038366317749,'*226': 6.20038366317749,'*227': 6.200374126434326,'*228': 6.200349807739258,'*229': 6.200349807739258,'*23': 6.201831340789795,'*230': 6.200349807739258,'*231': 6.200352191925049,'*232': 6.200352191925049,'*233': 6.200352668762207,'*234': 6.200352668762207,'*235': 6.200281143188477,'*236': 6.200281143188477,'*237': 6.200257778167725,'*238': 6.200257778167725,'*239': 6.200234413146973,'*24': 6.201822280883789,'*240': 6.200222015380859,'*241': 6.200222015380859,'*242': 6.200222015380859,'*243': 6.200189113616943,'*244': 6.2001423835754395,'*245': 6.2001423835754395,'*246': 6.2001423835754395,'*247': 6.2001423835754395,'*248': 6.2001423835754395,'*249': 6.2001423835754395,'*25': 6.201816082000732,'*250': 6.200117111206055,'*251': 6.200117111206055,'*252': 6.20006799697876,'*253': 6.200074672698975,'*254': 6.200074672698975,'*255': 6.200040817260742,'*256': 6.200023651123047,'*257': 6.200010776519775,'*258': 6.200010776519775,'*259': 6.200010776519775,'*26': 6.201816082000732,'*260': 6.200010776519775,'*261': 6.2000041007995605,'*262': 6.2000041007995605,'*263': 6.2000041007995605,'*264': 6.2000041007995605,'*265': 6.2000041007995605,'*266': 6.199972629547119,'*267': 6.199947834014893,'*268': 6.199952125549316,'*269': 6.199952125549316,'*27': 6.201816082000732,'*270': 6.199923038482666,'*271': 6.199893474578857,'*272': 6.199893474578857,'*273': 6.199865341186523,'*274': 6.199865341186523,'*275': 6.199865341186523,'*276': 6.199865341186523,'*277': 6.199844837188721,'*278': 6.199844837188721,'*279': 6.1998443603515625,'*28': 6.201807975769043,'*280': 6.1998443603515625,'*281': 6.1998443603515625,'*282': 6.1998443603515625,'*283': 6.199779987335205,'*284': 6.199779987335205,'*285': 6.199779987335205,'*286': 6.199779987335205,'*287': 6.199779987335205,'*288': 6.199778079986572,'*289': 6.199771881103516,'*29': 6.201807975769043,'*290': 6.199771881103516,'*291': 6.199758052825928,'*292': 6.199752330780029,'*293': 6.199752330780029,'*294': 6.199752330780029,'*295': 6.199752330780029,'*296': 6.199739456176758,'*297': 6.199734210968018,'*298': 6.199698448181152,'*299': 6.199672222137451,'*3': 6.202004432678223,'*30': 6.201807975769043,'*300': 6.199641227722168,'*301': 6.199641227722168,'*302': 6.196625232696533,'*303': 6.196559429168701,'*304': 6.196559429168701,'*305': 6.196538925170898,'*306': 6.196538925170898,'*307': 6.196538925170898,'*308': 6.196538925170898,'*309': 6.196538925170898,'*31': 6.201807975769043,'*310': 6.196509838104248,'*311': 6.196499824523926,'*312': 6.196498394012451,'*313': 6.196498394012451,'*314': 6.196496963500977,'*315': 6.196496963500977,'*316': 6.196496963500977,'*317': 6.196496963500977,'*318': 6.196484088897705,'*319': 6.196484088897705,'*32': 6.20179557800293,'*320': 6.19645881652832,'*321': 6.19645881652832,'*322': 6.19645881652832,'*323': 6.196468353271484,'*324': 6.196468353271484,'*325': 6.1964335441589355,'*326': 6.19639253616333,'*327': 6.19639253616333,'*328': 6.19639253616333,'*329': 6.19639253616333,'*33': 6.20179557800293,'*330': 6.19639253616333,'*331': 6.19639253616333,'*332': 6.19639253616333,'*333': 6.196359157562256,'*334': 6.196359157562256,'*335': 6.1963300704956055,'*336': 6.1963300704956055,'*337': 6.1963300704956055,'*338': 6.196277618408203,'*339': 6.196277141571045,'*34': 6.20179557800293,'*340': 6.196240425109863,'*341': 6.196230888366699,'*342': 6.196230888366699,'*343': 6.196174144744873,'*344': 6.1961669921875,'*345': 6.196159839630127,'*346': 6.196159839630127,'*347': 6.196159839630127,'*348': 6.196159839630127,'*349': 6.196159839630127,'*35': 6.201786041259766,'*350': 6.196159839630127,'*351': 6.196132659912109,'*352': 6.1961164474487305,'*353': 6.1961164474487305,'*354': 6.196106910705566,'*355': 6.196106910705566,'*356': 6.196094989776611,'*357': 6.196094989776611,'*358': 6.196094989776611,'*359': 6.196094989776611,'*36': 6.201786041259766,'*360': 6.196094989776611,'*361': 6.196101188659668,'*362': 6.196101188659668,'*363': 6.196101188659668,'*364': 6.1960768699646,'*365': 6.196057319641113,'*366': 6.1960368156433105,'*367': 6.1960368156433105,'*368': 6.1960368156433105,'*369': 6.196034908294678,'*37': 6.201786041259766,'*370': 6.196028232574463,'*371': 6.196028232574463,'*372': 6.196028232574463,'*373': 6.196028232574463,'*374': 6.196028232574463,'*375': 6.196022033691406,'*376': 6.196022033691406,'*377': 6.196022033691406,'*378': 6.196022033691406,'*379': 6.196022033691406,'*38': 6.201777935028076,'*380': 6.196022033691406,'*381': 6.196022033691406,'*382': 6.196022033691406,'*383': 6.196000576019287,'*384': 6.196000576019287,'*385': 6.1959991455078125,'*386': 6.19598913192749,'*387': 6.1959943771362305,'*388': 6.1959943771362305,'*389': 6.195987224578857,'*39': 6.201777935028076,'*390': 6.195987224578857,'*391': 6.1959614753723145,'*392': 6.1959614753723145,'*393': 6.1959614753723145,'*394': 6.195955276489258,'*395': 6.195930004119873,'*396': 6.195930004119873,'*397': 6.195930004119873,'*398': 6.195930004119873,'*399': 6.195930004119873,'*4': 6.202004432678223,'*40': 6.201777935028076,'*400': 6.1959147453308105,'*401': 6.195895671844482,'*402': 6.195895671844482,'*403': 6.195870399475098,'*404': 6.195854663848877,'*405': 6.195847988128662,'*406': 6.195847988128662,'*407': 6.1957902908325195,'*408': 6.1957902908325195,'*409': 6.1957902908325195,'*41': 6.201777935028076,'*410': 6.195760726928711,'*411': 6.195760726928711,'*412': 6.195760726928711,'*413': 6.195760726928711,'*414': 6.195734024047852,'*415': 6.195714950561523,'*416': 6.195714950561523,'*417': 6.195714950561523,'*418': 6.195699214935303,'*419': 6.195644855499268,'*42': 6.201782703399658,'*420': 6.195609092712402,'*421': 6.195609092712402,'*422': 6.195609092712402,'*423': 6.195609092712402,'*424': 6.195580005645752,'*425': 6.195580005645752,'*426': 6.195571422576904,'*427': 6.195571422576904,'*428': 6.195571422576904,'*429': 6.195571422576904,'*43': 6.201782703399658,'*430': 6.195571422576904,'*431': 6.195566177368164,'*432': 6.195566177368164,'*433': 6.195566177368164,'*434': 6.195566177368164,'*435': 6.195522785186768,'*436': 6.195522785186768,'*437': 6.195523262023926,'*438': 6.195523262023926,'*439': 6.195523262023926,'*44': 6.2017669677734375,'*440': 6.195523262023926,'*441': 6.195517063140869,'*442': 6.195517063140869,'*443': 6.195517063140869,'*444': 6.195517063140869,'*445': 6.195517063140869,'*446': 6.195517063140869,'*447': 6.195517063140869,'*448': 6.195517063140869,'*449': 6.195521354675293,'*45': 6.201745510101318,'*450': 6.195518493652344,'*451': 6.195518493652344,'*452': 6.195505619049072,'*453': 6.1954851150512695,'*454': 6.1954851150512695,'*455': 6.1954851150512695,'*456': 6.1954851150512695,'*457': 6.1954851150512695,'*458': 6.195457935333252,'*459': 6.19545841217041,'*46': 6.201740264892578,'*460': 6.195458889007568,'*461': 6.195443630218506,'*462': 6.195440292358398,'*463': 6.195440292358398,'*464': 6.195440292358398,'*465': 6.195432662963867,'*466': 6.195432662963867,'*467': 6.195432662963867,'*468': 6.195432662963867,'*469': 6.195433139801025,'*47': 6.201740264892578,'*470': 6.195433139801025,'*471': 6.195411682128906,'*472': 6.195411682128906,'*473': 6.195405006408691,'*474': 6.19537878036499,'*475': 6.195344924926758,'*476': 6.195309638977051,'*477': 6.195309638977051,'*478': 6.1952805519104,'*479': 6.195265769958496,'*48': 6.201740264892578,'*480': 6.195265769958496,'*481': 6.195265769958496,'*482': 6.195259094238281,'*483': 6.195259094238281,'*484': 6.195259094238281,'*485': 6.195234298706055,'*486': 6.195234298706055,'*487': 6.195234298706055,'*488': 6.195234298706055,'*489': 6.195234298706055,'*49': 6.201740264892578,'*490': 6.195224285125732,'*491': 6.195224285125732,'*492': 6.195224285125732,'*493': 6.195224285125732,'*494': 6.195224285125732,'*495': 6.195224285125732,'*496': 6.195223808288574,'*497': 6.195218563079834,'*498': 6.195218563079834,'*499': 6.195218563079834,'*5': 6.20200777053833,'*50': 6.201714515686035,'*500': 6.195166110992432,'*501': 6.195166110992432,'*502': 6.195166110992432,'*503': 6.195115566253662,'*504': 6.195115566253662,'*505': 6.195098400115967,'*506': 6.195098400115967,'*507': 6.195098400115967,'*51': 6.201714515686035,'*52': 6.201714992523193,'*53': 6.201714992523193,'*54': 6.201714992523193,'*55': 6.201714992523193,'*56': 6.201687812805176,'*57': 6.201681613922119,'*58': 6.201670169830322,'*59': 6.201670169830322,'*6': 6.201994895935059,'*60': 6.201659202575684,'*61': 6.20165491104126,'*62': 6.201644420623779,'*63': 6.201629161834717,'*64': 6.201629161834717,'*65': 6.201629161834717,'*66': 6.201629161834717,'*67': 6.201629161834717,'*68': 6.201629161834717,'*69': 6.201620578765869,'*7': 6.201960563659668,'*70': 6.201612949371338,'*71': 6.201612949371338,'*72': 6.201606273651123,'*73': 6.201574802398682,'*74': 6.201575756072998,'*75': 6.201575756072998,'*76': 6.201569080352783,'*77': 6.20154333114624,'*78': 6.20154333114624,'*79': 6.20154333114624,'*8': 6.201960563659668,'*80': 6.201538562774658,'*81': 6.201538562774658,'*82': 6.201499938964844,'*83': 6.201490879058838,'*84': 6.201490879058838,'*85': 6.201490879058838,'*86': 6.201488494873047,'*87': 6.201488494873047,'*88': 6.201488494873047,'*89': 6.201478958129883,'*9': 6.2019429206848145,'*90': 6.201478958129883,'*91': 6.201478958129883,'*92': 6.201456546783447,'*93': 6.201431751251221,'*94': 6.201431751251221,'*95': 6.201416492462158,'*96': 6.20139741897583,'*97': 6.201372146606445,'*98': 6.201372146606445,'*99': 6.201372146606445}
            exp_bmin_dict = self.exp_dicts['exp_bmin_dict']
            # test_standard_cube_briggsbwtaper.exp_bmaj_dict
            exp_bmaj_dict = self.exp_dicts['exp_bmaj_dict']
            # test_standard_cube_briggsbwtaper.exp_pa_dict
            exp_pa_dict = self.exp_dicts['exp_pa_dict']


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

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')


        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cube_briggsbwtaper
#-------------------------------------------------#
    # Test 2
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mfs(self):
        ''' Standard (single field) MFS imaging - central field of NGC5363 (field 2), spw 16 & 22 '''

        test_name = self._testMethodName 
        file_name = 'standard_mfs.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00020.S_tclean.ms')
        self.getExpdicts(test_name)

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
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[209.03005051deg, 5.25476861deg], [71.9366arcsec, 71.6106arcsec], 0.00000000deg]')

        # test_standard_mfs_exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[209.03018549deg, 5.25490301deg], [18.5483arcsec, 11.7743arcsec], 90.00000000deg]')

        # test_standard_mfs_exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[209.02978056deg, 5.25512703deg], [18.0644arcsec, 11.9355arcsec], 90.00000000deg]')

        # test_standard_mfs_exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.resid', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse [[209.02974119deg, 5.25476301deg], [2.7621arcsec, 1.8750arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_mfs_exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_mfs_exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed=self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "E2E6.1.00020.S_tclean.ms")

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []


        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.003, 0.04])
        self.mom8_creator(img+'.residual', range_list=[-0.003, 0.04])
        test_dict[test_name]['images'].extend( \
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

        test_name = self._testMethodName
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
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', fit_region = \
            'ellipse[[209.03005051deg, 5.25476861deg], [71.9366arcsec, 71.6106arcsec], 0.00000000deg]')

        # test_standard_mtmfs_exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb.tt0', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', fit_region = \
            'ellipse[[209.03018549deg, 5.25490301deg], [18.5483arcsec, 11.7743arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf.tt0', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            'ellipse[[209.02978056deg, 5.25512703deg], [18.0644arcsec, 11.9355arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual.tt0', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', fit_region = \
            'ellipse[[209.02974119deg, 5.25476301deg], [2.7621arcsec, 1.8750arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_mtmfs_exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        # test_standard_mtmfs_exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt.tt0', epsilon=self.epsilon)


        # .image.tt1 report
        im1_stats_dict = self.image_stats(img+'.image.tt1', fit_region = \
            'ellipse[[209.03000552deg, 5.25503742deg], [16.2902arcsec, 10.3226arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_im1_stats
        exp_im1_stats = self.exp_dicts['exp_im1_stats']
 
        report9 = th.check_dict_vals(exp_im1_stats, im1_stats_dict, '.image.tt1', epsilon=self.epsilon)

        # .residual.tt1 report
        resid1_stats_dict = self.image_stats(img+'.residual.tt1', fit_region = \
            'ellipse[[209.02978056deg, 5.25512703deg], [18.0644arcsec, 11.9355arcsec], 90.00000000deg]')

        # test_standard_mtmfs_exp_resid1_stats
        exp_resid1_stats =  self.exp_dicts['exp_resid1_stats']

        report10 = th.check_dict_vals(exp_resid1_stats, resid1_stats_dict, \
            '.residual.tt1', epsilon=self.epsilon)

        # .model.tt1 report
        model1_stats_dict = self.image_stats(img+'.model.tt1', fit_region = \
            'ellipse[[209.02974119deg, 5.25476301deg], [2.7621arcsec, 1.8750arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_mtmfs_exp_model1_stats
        exp_model1_stats = self.exp_dicts['exp_model1_stats']

        report11 = th.check_dict_vals(exp_model1_stats, model1_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt.tt1 report
        sumwt1_stats_dict = self.image_stats(img+'.sumwt.tt1')

        # test_standard_mtmfs_exp_sumwt1_stats
        exp_sumwt1_stats = self.exp_dicts['exp_sumwt1_stats']

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

       
        test_name = self._testMethodName
        file_name = 'standard_cube_eph.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData([data_path+'2017.1.00750.T_tclean_exe1.ms', 
            data_path+'2017.1.00750.T_tclean_exe2.ms'])
        self.getExpdicts(test_name)

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='21PGiacobini-Zinner', spw=['0', '0'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11', '0,1,2,3,4,5,6,7,8,9'], \
            scan=['7,11,15,19,23','8,12,16,20,24'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[80, 80], cell=['0.66arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='cubesource', \
            #nchan=1000, start=1550, width=1, perchanweightdensity=False, \
            nchan=1000, start='354.4452721710GHz', width='0.1221004MHz', perchanweightdensity=False, \
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
            #nchan=1000, start=1550, width=1, perchanweightdensity=False, \
            nchan=1000, start='354.4452721710GHz', width='0.1221004MHz', perchanweightdensity=False, \
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

        exp_im_stats = self.exp_dicts['exp_im_stats']


        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [40, 72, 0, 0]), \
                (img+'.image', False, [40, 73, 0, 0]), \
                (img+'.image', True, [8, 40, 0, 0]), \
                (img+'.image', False, [7, 40, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

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

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')



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
        self.getExpdicts(test_name)

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='21PGiacobini-Zinner', spw=['0', '0'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11', '0,1,2,3,4,5,6,7,8,9'], \
            scan=['7,11,15,19,23','8,12,16,20,24'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[80, 80], cell=['0.66arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='cubesource', \
            nchan=1000, start='354.4452721710GHz', width='0.1221004MHz', perchanweightdensity=True, \
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
            nchan=1000, start='354.4452721710GHz', width='0.1221004MHz', perchanweightdensity=True, \
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
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats'] 

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_eph_pcwdT.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_cube_eph_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

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
            
            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')

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
        self.getExpdicts(test_name)

        print("\nSTARTING: iter0 routine")

        # iter0 routine
        tclean(vis=self.msfile, field='21PGiacobini-Zinner', spw=['0', '0'], \
            antenna=['0,1,2,3,4,5,6,7,8,9,10,11', '0,1,2,3,4,5,6,7,8,9'], \
            scan=['7,11,15,19,23','8,12,16,20,24'], \
            intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='data', \
            imagename=file_name+'0', imsize=[80, 80], cell=['0.66arcsec'], \
            phasecenter='TRACKFIELD', stokes='I', specmode='cubesource', \
            nchan=1000, start='354.4452721710GHz', width='0.1221004MHz', perchanweightdensity=True, \
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
            nchan=1000, start='354.4452721710GHz', width='0.1221004MHz', perchanweightdensity=True, \
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
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_eph_briggsbwtaper.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_cube_eph_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

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
            
            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cube_eph_briggsbwtaper
#-------------------------------------------------#
    # Test 5
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mfs_eph(self):
        ''' Standard (single field) ephemeris mfs imaging - central field of Venus (field 2), spw 25 & 45 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')
        self.getExpdicts(test_name)

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

        exp_im_stats = self.exp_dicts['exp_im_stats']

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [144, 266, 0, 0]), \
                (img+'.image', False, [144, 267, 0, 0]), \
                (img+'.image', True, [22, 145, 0, 0]), \
                (img+'.image', False, [21, 145, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[239.37090838deg, -16.96415647deg], [17.8437arcsec, 17.4772arcsec], 90.00000000deg]')

        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.37094084deg, -16.96415506deg], [1.1279arcsec, 0.7875arcsec], 90.00000000deg]')

        exp_psf_stats = self.exp_dicts['exp_psf_stats']


        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-1.05, 1.05])
        self.mom8_creator(img+'.residual', range_list=[-1.05, 1.05])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

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
            
            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_mfs_eph
#-------------------------------------------------#
    # Test 6
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_mtmfs_eph(self):
        ''' Standard (single field) ephemeris mtmfs imaging - central field of Venus (field 2), spw 25 & 45 '''

        test_name = self._testMethodName
        file_name = self.remove_prefix(test_name,'test_')+'.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')
        self.getExpdicts(test_name)

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

        exp_im_stats = self.exp_dicts['exp_im_stats']

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image.tt0', True, [144, 266, 0, 0]), \
                (img+'.image.tt0', False, [144, 267, 0, 0]), \
                (img+'.image.tt0', True, [22, 145, 0, 0]), \
                (img+'.image.tt0', False, [21, 145, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb.tt0', fit_region = \
            'ellipse[[239.37090838deg, -16.96415647deg], [17.8437arcsec, 17.4772arcsec], 90.00000000deg]')

        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf.tt0', fit_region = \
            'ellipse[[239.37094084deg, -16.96415506deg], [1.1279arcsec, 0.7875arcsec], 90.00000000deg]')

        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual.tt0', \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model.tt0', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = self.exp_dicts['exp_model_stats']
 
        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt.tt0')

        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .image.tt1 report
        im1_stats_dict = self.image_stats(img+'.image.tt1', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_im1_stats = self.exp_dicts['exp_im1_stats']
 
        report9 = th.check_dict_vals(exp_im1_stats, im1_stats_dict, '.image.tt1', epsilon=self.epsilon)

        # .residual.tt1 report
        resid1_stats_dict = self.image_stats(img+'.residual.tt1', \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid1_stats = self.exp_dicts['exp_resid1_stats']

        report10 = th.check_dict_vals(exp_resid1_stats, resid1_stats_dict, \
            '.residual.tt1', epsilon=self.epsilon)

        # .model.tt1 report
        model1_stats_dict = self.image_stats(img+'.model.tt1', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model1_stats = self.exp_dicts['exp_model1_stats']

        report11 = th.check_dict_vals(exp_model1_stats, model1_stats_dict, \
            '.model.tt0', epsilon=self.epsilon)

        # .sumwt.tt1 report
        sumwt1_stats_dict = self.image_stats(img+'.sumwt.tt1')

        exp_sumwt1_stats = self.exp_dicts['exp_sumwt1_stats']

        report12 = th.check_dict_vals(exp_sumwt1_stats, sumwt1_stats_dict, \
            '.sumwt.tt1', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9 + report10 + report11 + report12

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image.tt0', range_list=[-1.05, 1.05])
        self.mom8_creator(img+'.residual.tt0', range_list=[-1.05, 1.05])
        test_dict[test_name]['images'].extend( \
            (img+'.image.tt0.moment8.png',img+'.residual.tt0.moment8.png'))

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
            savedict['im1_stats_dict']=im1_stats_dict
            savedict['resid1_stats_dict']=resid1_stats_dict
            savedict['model1_stats_dict']=model1_stats_dict
            savedict['sumwt1_stats_dict']=sumwt1_stats_dict

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')


        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_mtmfs_eph
#-------------------------------------------------#
    # Test 7
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cal(self):
        ''' Calibrator image - field J2258-2758, spw 22 '''

        test_name = self._testMethodName 
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'E2E6.1.00034.S_tclean.ms')
        self.getExpdicts(test_name)

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

        exp_im_stats = self.exp_dicts['exp_im_stats']

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [45, 85, 0, 0]), \
                (img+'.image', False, [45, 86, 0, 0]), \
                (img+'.image', True, [5, 45, 0, 0]), \
                (img+'.image', False, [4, 45, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[344.52476094deg, -27.97251802deg], [34.7828arcsec, 34.7011arcsec], 90.00000000deg]')

        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[344.52484287deg, -27.97253611deg], [8.0802arcsec, 4.8086arcsec], 90.00000000deg]')

        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[344.52480945deg, -27.97253944deg], [12.1076arcsec, 6.2463arcsec], 90.00000000deg]')

        exp_resid_stats = self.exp_dicts['exp_resid_stats']
 
        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[344.52480945deg, -27.97253944deg], [12.1076arcsec, 6.2463arcsec], 90.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # report combination
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8


        failed = self.filter_report(report)

        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-0.015, 2.5])
        self.mom8_creator(img+'.residual', range_list=[-0.015, 2.5])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

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

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

# End of test_standard_cal
#-------------------------------------------------#
    # Test 8
    @stats_dict(test_dict)
    # @unittest.skip("")
    def test_standard_cal_eph(self):
        ''' Standard (single field) ephemeris calibrator imaging - central field of Venus (field 2), spw 25 & 45 '''

        test_name = self._testMethodName 
        file_name = self.remove_prefix(test_name, 'test_')+'.iter'
        img = os.getcwd()+'/'+file_name+'1'
        self.prepData(data_path+'2018.1.00879.S_tclean.ms')
        self.getExpdicts(test_name)

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

        exp_im_stats = self.exp_dicts['exp_im_stats']

        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [144, 266, 0, 0]), \
                (img+'.image', False, [144, 267, 0, 0]), \
                (img+'.image', True, [22, 145, 0, 0]), \
                (img+'.image', False, [21, 145, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', self.epsilon)

        # .mask report (test_standard_cal_eph)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', self.epsilon)

        # .pb report (test_standard_cal_eph)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[239.37090838deg, -16.96415647deg], [17.8437arcsec, 17.4772arcsec], 90.00000000deg]')

        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', self.epsilon)

        # .psf report (test_standard_cal_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.37094084deg, -16.96415506deg], [1.1279arcsec, 0.7875arcsec], 90.00000000deg]')

        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', self.epsilon)

        # .residual report (test_standard_cal_eph)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]')

        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', self.epsilon)

        # .model report (test_standard_cal_eph)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.37089658deg, -16.96414518deg], [12.9657arcsec, 12.4377arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', self.epsilon)

        # .sumwt report (test_standard_cal_eph)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', self.epsilon)

        # report combination (test_standard_cal_eph)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8

        failed = self.filter_report(report)

        add_to_dict(self, output = test_dict, dataset = \
            "2018.1.00879.S_tclean.ms")

        test_dict[test_name]['self.parallel'] = self.parallel
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = []

        img = shutil._basename(img)
        self.mom8_creator(img+'.image', range_list=[-1.05, 1.05])
        self.mom8_creator(img+'.residual', range_list=[-1.05, 1.05])
        test_dict[test_name]['images'].extend( \
            (img+'.image.moment8.png',img+'.residual.moment8.png'))

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

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')

        self.assertTrue(th.check_final(pstr = report), \
            msg = failed)

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
            restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggs', robust=0.5,\
            npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel, verbose=True)

        # retrieve per-channel beam statistics
        bmin_dict, bmaj_dict, pa_dict = \
            self.cube_beam_stats(img+'.psf')

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

        exp_im_stats = self.exp_dicts['exp_im_stats']


        report1 = th.checkall( \
            # checks for image and pb mask movement
            imgmask = [(img+'.image', True, [51, 99, 0, 0]), \
                      (img+'.image', False, [51, 100, 0, 0]), \
                      (img+'.image', True, [9, 56, 0, 0]), \
                      (img+'.image', False, [8, 56, 0, 0])])

        report2 = th.check_dict_vals(exp_im_stats, im_stats_dict, '.image', epsilon=self.epsilon)

        # .mask report (test_mosaic_cube)
        mask_stats_dict = self.image_stats(img+'.mask')

        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        exp_resid_stats = self.exp_dicts['exp_resid_stats']
 
        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube
        exp_wt_stats = self.exp_dicts['exp_wt_stats']

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            exp_bmin_dict = self.exp_dicts['exp_bmin_dict']
            exp_bmaj_dict = self.exp_dicts['exp_bmaj_dict']
            exp_pa_dict = self.exp_dicts['exp_pa_dict']

            report += self.check_dict_vals_beam(exp_bmin_dict, bmin_dict, '.image bmin', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_bmaj_dict, bmaj_dict, '.image bmaj', epsilon=self.epsilon)
            report += self.check_dict_vals_beam(exp_pa_dict, pa_dict, '.image pa', epsilon=self.epsilon)

        failed = self.filter_report(report)

        img = shutil._basename(img)
        add_to_dict(self, output=test_dict, dataset = \
            "E2E6.1.00034.S_tclean.ms")

        self.modify_dict(test_dict, 'test_mosaic_cube', self.parallel)

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
            savedict['wt_stats_dict']=wt_stats_dict
            savedict['bmin_dict']=bmin_stats_dict
            savedict['bmaj_dict']=bmaj_stats_dict
            savedict['pa_dict']=pa_stats_dict

            self.save_dict_to_file(test_name, savedict, test_name+'_cur_stats')

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
        self.getExpdicts(test_name)

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
            restoration=True, restoringbeam='common', pbcor=True, \
            weighting='briggs', robust=0.5,\
            npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel, verbose=True)

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
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']
 
        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube_pcwdT)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        # test_mosaic_cube_pcwdT.exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_pcwdT)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        # test_mosaic_cube_pcwdT.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_pcwdT)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        # test_mosaic_cube_pcwdT.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_pcwdT)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_pcwdT.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube_pcwdT.exp_wt_stats
        exp_wt_stats = self.exp_dicts['exp_wt_stats']

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            # test_mosaic_cube_pcwdT.exp_bmin_dict
            exp_bmin_dict = self.exp_dicts['exp_bmin_dict']

            # test_mosaic_cube_pcwdT.exp_bmaj_dict
            exp_bmaj_dict = self.exp_dicts['exp_bmaj_dict']

            # test_mosaic_cube_pcwdT.exp_pa_dict
            exp_pa_dict = {'*0': 67.42774200439453,'*1': 67.42774200439453,'*10': 67.42865753173828,'*100': 67.46749114990234,'*101': 67.46749114990234,'*102': 67.46749114990234,'*103': 67.46749114990234,'*104': 67.46749114990234,'*105': 67.46749114990234,'*106': 67.46749114990234,'*107': 67.4686050415039,'*108': 67.4686050415039,'*109': 67.4686050415039,'*11': 67.42865753173828,'*110': 67.4686050415039,'*111': 67.46869659423828,'*112': 67.40909576416016,'*113': 67.40939331054688,'*114': 67.40939331054688,'*115': 67.40939331054688,'*116': 67.40939331054688,'*117': 67.40909576416016,'*118': 67.40943145751953,'*119': 67.40943145751953,'*12': 67.42815399169922,'*120': 67.40927124023438,'*121': 67.40985870361328,'*122': 67.40985870361328,'*123': 67.40985870361328,'*124': 67.40985870361328,'*125': 67.40985870361328,'*126': 67.40985870361328,'*127': 67.40949249267578,'*128': 67.40949249267578,'*129': 67.40902709960938,'*13': 67.42815399169922,'*130': 67.40902709960938,'*131': 67.40931701660156,'*132': 67.40931701660156,'*133': 67.40931701660156,'*134': 67.40931701660156,'*135': 67.40931701660156,'*136': 67.40931701660156,'*137': 67.40974426269531,'*138': 67.40974426269531,'*139': 67.40974426269531,'*14': 67.42815399169922,'*140': 67.40974426269531,'*141': 67.40997314453125,'*142': 67.40997314453125,'*143': 67.40997314453125,'*144': 67.40997314453125,'*145': 67.40997314453125,'*146': 67.40936279296875,'*147': 67.40936279296875,'*148': 67.40936279296875,'*149': 67.4095687866211,'*15': 67.42796325683594,'*150': 67.4095687866211,'*151': 67.4095687866211,'*152': 67.40923309326172,'*153': 67.40923309326172,'*154': 67.40921020507812,'*155': 67.40921020507812,'*156': 67.40921020507812,'*157': 67.40921020507812,'*158': 67.40969848632812,'*159': 67.40969848632812,'*16': 67.42796325683594,'*160': 67.40969848632812,'*161': 67.40969848632812,'*162': 67.40969848632812,'*163': 67.40887451171875,'*164': 67.40887451171875,'*165': 67.40887451171875,'*166': 67.40887451171875,'*167': 67.40887451171875,'*168': 67.40816497802734,'*169': 67.40816497802734,'*17': 67.42796325683594,'*170': 67.40816497802734,'*171': 67.40816497802734,'*172': 67.40882873535156,'*173': 67.40799713134766,'*174': 67.40827941894531,'*175': 67.40790557861328,'*176': 67.4001235961914,'*177': 67.39026641845703,'*178': 67.38996124267578,'*179': 67.38996124267578,'*18': 67.42818450927734,'*180': 67.38996124267578,'*181': 67.3896255493164,'*182': 67.38994598388672,'*183': 67.38994598388672,'*184': 67.38994598388672,'*185': 67.38994598388672,'*186': 67.38994598388672,'*187': 67.38994598388672,'*188': 67.38994598388672,'*189': 67.38983917236328,'*19': 67.42818450927734,'*190': 67.38983917236328,'*191': 67.39041137695312,'*192': 67.39041137695312,'*193': 67.39024353027344,'*194': 67.39043426513672,'*195': 67.39043426513672,'*196': 67.39043426513672,'*197': 67.38996124267578,'*198': 67.38996124267578,'*199': 67.38996124267578,'*2': 67.42774200439453,'*20': 67.42818450927734,'*200': 67.38945770263672,'*201': 67.38945770263672,'*202': 67.38945770263672,'*203': 67.38945770263672,'*204': 67.45548248291016,'*205': 67.45548248291016,'*206': 67.47216033935547,'*207': 67.47216033935547,'*208': 67.47383117675781,'*209': 67.47413635253906,'*21': 67.42818450927734,'*210': 67.47413635253906,'*211': 67.47418212890625,'*212': 67.47332763671875,'*213': 67.4735107421875,'*214': 67.4735107421875,'*215': 67.47306823730469,'*216': 67.47306823730469,'*217': 67.47306823730469,'*218': 67.47306823730469,'*219': 67.47328186035156,'*22': 67.4267807006836,'*220': 67.47045135498047,'*221': 67.43890380859375,'*222': 67.43890380859375,'*223': 67.43890380859375,'*224': 67.43890380859375,'*225': 67.4367446899414,'*226': 67.4367446899414,'*227': 67.4367446899414,'*228': 67.4367446899414,'*229': 67.4367446899414,'*23': 67.4267807006836,'*230': 67.4367446899414,'*231': 67.43634796142578,'*232': 67.43634796142578,'*233': 67.43634796142578,'*234': 67.43634796142578,'*235': 67.4362564086914,'*236': 67.4362564086914,'*237': 67.43624877929688,'*238': 67.43585968017578,'*239': 67.43585968017578,'*24': 67.4267807006836,'*240': 67.43636322021484,'*241': 67.43636322021484,'*242': 67.43636322021484,'*243': 67.43710327148438,'*244': 67.43710327148438,'*245': 67.43710327148438,'*246': 67.43762969970703,'*247': 67.43688201904297,'*248': 67.43688201904297,'*249': 67.43688201904297,'*25': 67.4267807006836,'*250': 67.43717956542969,'*251': 67.43717956542969,'*252': 67.43766021728516,'*253': 67.43732452392578,'*254': 67.43732452392578,'*255': 67.43732452392578,'*256': 67.43777465820312,'*257': 67.43777465820312,'*258': 67.43777465820312,'*259': 67.4376449584961,'*26': 67.4267807006836,'*260': 67.44227600097656,'*261': 67.44729614257812,'*262': 67.44729614257812,'*263': 67.4476547241211,'*264': 67.4476547241211,'*265': 67.4476547241211,'*266': 67.4476547241211,'*267': 67.4476547241211,'*268': 67.4476547241211,'*269': 67.44754791259766,'*27': 67.42683410644531,'*270': 67.44754791259766,'*271': 67.44758605957031,'*272': 67.44739532470703,'*273': 67.44739532470703,'*274': 67.44739532470703,'*275': 67.44739532470703,'*276': 67.44766235351562,'*277': 67.44766235351562,'*278': 67.44766235351562,'*279': 67.44766235351562,'*28': 67.42683410644531,'*280': 67.4478759765625,'*281': 67.4478759765625,'*282': 67.44815826416016,'*283': 67.4476547241211,'*284': 67.4476547241211,'*285': 67.4476547241211,'*286': 67.44598388671875,'*287': 67.44554901123047,'*288': 67.44554901123047,'*289': 67.44438934326172,'*29': 67.42701721191406,'*290': 67.47785186767578,'*291': 67.47785186767578,'*292': 67.4777603149414,'*293': 67.4777603149414,'*294': 67.4777603149414,'*295': 67.47737121582031,'*296': 67.47737121582031,'*297': 67.47737121582031,'*298': 67.47763061523438,'*299': 67.47763061523438,'*3': 67.42806243896484,'*30': 67.42701721191406,'*300': 67.47763061523438,'*301': 67.47763061523438,'*302': 67.47763061523438,'*303': 67.47791290283203,'*304': 67.47791290283203,'*305': 67.47791290283203,'*306': 67.47791290283203,'*307': 67.47791290283203,'*308': 67.47820281982422,'*309': 67.47786712646484,'*31': 67.42622375488281,'*310': 67.47786712646484,'*311': 67.47786712646484,'*312': 67.47786712646484,'*313': 67.47786712646484,'*314': 67.47786712646484,'*315': 67.47734832763672,'*316': 67.47640228271484,'*317': 67.47640228271484,'*318': 67.47581481933594,'*319': 67.47581481933594,'*32': 67.42555236816406,'*320': 67.47581481933594,'*321': 67.47581481933594,'*322': 67.47581481933594,'*323': 67.47581481933594,'*324': 67.47581481933594,'*325': 67.47581481933594,'*326': 67.47581481933594,'*327': 67.47581481933594,'*328': 67.47581481933594,'*329': 67.47581481933594,'*33': 67.42555236816406,'*330': 67.47581481933594,'*331': 67.47685241699219,'*332': 67.47685241699219,'*333': 67.47685241699219,'*334': 67.47685241699219,'*335': 67.47685241699219,'*336': 67.47913360595703,'*337': 67.47913360595703,'*338': 67.4788818359375,'*339': 67.4780502319336,'*34': 67.42555236816406,'*340': 67.4780502319336,'*341': 67.4780502319336,'*342': 67.4780502319336,'*343': 67.4780502319336,'*344': 67.48170471191406,'*345': 67.4820556640625,'*346': 67.48235321044922,'*347': 67.4825439453125,'*348': 67.48287200927734,'*349': 67.48287200927734,'*35': 67.42527770996094,'*350': 67.48243713378906,'*351': 67.48243713378906,'*352': 67.48243713378906,'*353': 67.48243713378906,'*354': 67.48257446289062,'*355': 67.48257446289062,'*356': 67.48257446289062,'*357': 67.48297119140625,'*358': 67.48297119140625,'*359': 67.48279571533203,'*36': 67.42559814453125,'*360': 67.48289489746094,'*361': 67.48348999023438,'*362': 67.48348999023438,'*363': 67.48347473144531,'*364': 67.48347473144531,'*365': 67.48347473144531,'*366': 67.48347473144531,'*367': 67.48367309570312,'*368': 67.48367309570312,'*369': 67.52080535888672,'*37': 67.42813110351562,'*370': 67.52080535888672,'*371': 67.52066802978516,'*372': 67.52066802978516,'*373': 67.53384399414062,'*374': 67.53177642822266,'*375': 67.53177642822266,'*376': 67.53177642822266,'*377': 67.53177642822266,'*378': 67.53146362304688,'*379': 67.53092193603516,'*38': 67.42813110351562,'*380': 67.5309829711914,'*381': 67.52033996582031,'*382': 67.5829086303711,'*383': 67.5829086303711,'*384': 67.58301544189453,'*385': 67.58301544189453,'*386': 67.5822525024414,'*387': 67.58255004882812,'*388': 67.58284759521484,'*389': 67.58366394042969,'*39': 67.42864990234375,'*390': 67.58366394042969,'*391': 67.58366394042969,'*392': 67.58366394042969,'*393': 67.58366394042969,'*394': 67.58366394042969,'*395': 67.58395385742188,'*396': 67.58395385742188,'*397': 67.58430480957031,'*398': 67.58397674560547,'*399': 67.58397674560547,'*4': 67.42806243896484,'*40': 67.42864990234375,'*400': 67.58397674560547,'*401': 67.58397674560547,'*402': 67.58397674560547,'*403': 67.58432006835938,'*404': 67.58432006835938,'*405': 67.58432006835938,'*406': 67.5838623046875,'*407': 67.5838623046875,'*408': 67.5838623046875,'*409': 67.5838623046875,'*41': 67.42864990234375,'*410': 67.5838623046875,'*411': 67.5838623046875,'*412': 67.5838623046875,'*413': 67.5837173461914,'*414': 67.5837173461914,'*415': 67.58311462402344,'*416': 67.58311462402344,'*417': 67.58311462402344,'*418': 67.58311462402344,'*419': 67.58311462402344,'*42': 67.42789459228516,'*420': 67.58311462402344,'*421': 67.58311462402344,'*422': 67.58311462402344,'*423': 67.58311462402344,'*424': 67.58311462402344,'*425': 67.58311462402344,'*426': 67.58353424072266,'*427': 67.5837173461914,'*428': 67.58344268798828,'*429': 67.58344268798828,'*43': 67.42789459228516,'*430': 67.58344268798828,'*431': 67.58275604248047,'*432': 67.58275604248047,'*433': 67.58304595947266,'*434': 67.58304595947266,'*435': 67.58304595947266,'*436': 67.58304595947266,'*437': 67.58304595947266,'*438': 67.58304595947266,'*439': 67.58304595947266,'*44': 67.42789459228516,'*440': 67.5838851928711,'*441': 67.5838851928711,'*442': 67.5838851928711,'*443': 67.5838851928711,'*444': 67.5840835571289,'*445': 67.5840835571289,'*446': 67.5840835571289,'*447': 67.5840835571289,'*448': 67.58444213867188,'*449': 67.58444213867188,'*45': 67.42789459228516,'*450': 67.58444213867188,'*451': 67.58393096923828,'*452': 67.58393096923828,'*453': 67.58393096923828,'*454': 67.58477783203125,'*455': 67.5848159790039,'*456': 67.58538055419922,'*457': 67.58538055419922,'*458': 67.58538055419922,'*459': 67.58595275878906,'*46': 67.42789459228516,'*460': 67.58595275878906,'*461': 67.5853042602539,'*462': 67.5853042602539,'*463': 67.5853042602539,'*464': 67.5853042602539,'*465': 67.5853042602539,'*466': 67.5853042602539,'*467': 67.5853042602539,'*468': 67.5853042602539,'*469': 67.5853042602539,'*47': 67.42756652832031,'*470': 67.5853042602539,'*471': 67.5853042602539,'*472': 67.5853042602539,'*473': 67.5853042602539,'*474': 67.5853042602539,'*475': 67.5853042602539,'*476': 67.5853042602539,'*477': 67.58475494384766,'*478': 67.58609771728516,'*479': 67.58577728271484,'*48': 67.4280776977539,'*480': 67.58577728271484,'*481': 67.58577728271484,'*482': 67.58577728271484,'*483': 67.5859603881836,'*484': 67.5859603881836,'*485': 67.58699798583984,'*486': 67.58699798583984,'*487': 67.58699798583984,'*488': 67.58699798583984,'*489': 67.5868911743164,'*49': 67.42770385742188,'*490': 67.5868911743164,'*491': 67.5868911743164,'*492': 67.5868911743164,'*493': 67.5868911743164,'*494': 67.58709716796875,'*495': 67.58709716796875,'*496': 67.58740997314453,'*497': 67.58740997314453,'*498': 67.58740997314453,'*499': 67.58714294433594,'*5': 67.42806243896484,'*50': 67.42770385742188,'*500': 67.58714294433594,'*501': 67.58714294433594,'*502': 67.58714294433594,'*503': 67.58714294433594,'*504': 67.58756256103516,'*505': 67.58756256103516,'*506': 67.58756256103516,'*507': 67.58756256103516,'*51': 67.42770385742188,'*52': 67.42770385742188,'*53': 67.42770385742188,'*54': 67.42770385742188,'*55': 67.42770385742188,'*56': 67.42770385742188,'*57': 67.42770385742188,'*58': 67.42686462402344,'*59': 67.42686462402344,'*6': 67.42777252197266,'*60': 67.42686462402344,'*61': 67.42686462402344,'*62': 67.42632293701172,'*63': 67.42632293701172,'*64': 67.42632293701172,'*65': 67.42557525634766,'*66': 67.42557525634766,'*67': 67.42557525634766,'*68': 67.42557525634766,'*69': 67.42557525634766,'*7': 67.42777252197266,'*70': 67.42498016357422,'*71': 67.42517852783203,'*72': 67.42517852783203,'*73': 67.42532348632812,'*74': 67.42532348632812,'*75': 67.42532348632812,'*76': 67.42532348632812,'*77': 67.42509460449219,'*78': 67.42562866210938,'*79': 67.42562866210938,'*8': 67.42777252197266,'*80': 67.42562866210938,'*81': 67.4272689819336,'*82': 67.4272689819336,'*83': 67.4272689819336,'*84': 67.4272689819336,'*85': 67.4272689819336,'*86': 67.4271011352539,'*87': 67.46137237548828,'*88': 67.46137237548828,'*89': 67.46137237548828,'*9': 67.42777252197266,'*90': 67.46118927001953,'*91': 67.46118927001953,'*92': 67.46118927001953,'*93': 67.46118927001953,'*94': 67.46138000488281,'*95': 67.46138000488281,'*96': 67.46138000488281,'*97': 67.46085357666016,'*98': 67.46085357666016,'*99': 67.46749114990234}
            exp_pa_dict = self.exp_dicts['exp_pa_dicts']

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

            self.save_dict_to_file(test_name,savedict, test_name+'_cur_stats')

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
        self.getExpdicts(test_name)

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
            restoration=True, restoringbeam='common', \
            pbcor=True, weighting='briggsbwtaper', robust=0.5,\
            npixels=0, niter=20000, threshold='0.354Jy', nsigma=0.0, \
            interactive=0, usemask='auto-multithresh', \
            sidelobethreshold=1.25, noisethreshold=5.0, \
            lownoisethreshold=2.0, negativethreshold=0.0, \
            minbeamfrac=0.1, growiterations=75, dogrowprune=True, \
            minpercentchange=1.0, fastnoise=False, restart=True, \
            savemodel='none', calcres=False, calcpsf=False, \
            parallel=self.parallel, verbose=True)

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
        exp_im_stats = self.exp_dicts['exp_im_stats']

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
        exp_mask_stats = self.exp_dicts['exp_mask_stats']

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report (test_mosaic_cube_briggsbwtaper)
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[11.47666677deg, -73.25825652deg], [52.6715arcsec, 52.2589arcsec], 0.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_pb_stats
        exp_pb_stats = self.exp_dicts['exp_pb_stats']

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_briggsbwtaper)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[11.47632032deg, -73.25823681deg], [8.7257arcsec, 8.0720arcsec], 90.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_psf_stats
        exp_psf_stats = self.exp_dicts['exp_psf_stats']

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse [[11.48661818deg, -73.26292371deg], [8.2211arcsec, 7.4698arcsec], 90.00000000deg]')

        # test_mosaic_cube_briggsbwtaper.exp_resid_stats
        exp_resid_stats = self.exp_dicts['exp_resid_stats']

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_briggsbwtaper)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[11.48109199deg, -73.25974151deg], [18.9246arcsec, 17.1916arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_briggsbwtaper.exp_model_stats
        exp_model_stats = self.exp_dicts['exp_model_stats']

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = self.exp_dicts['exp_sumwt_stats']

        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])
        
        #test_mosaic_cube_briggsbwtaper.exp_wt_stats
        exp_wt_stats = self.exp_dicts['exp_wt_stats']

        report9 = th.check_dict_vals(exp_wt_stats, wt_stats_dict, '.weight', epsilon=self.epsilon)

        # report combination (test_mosaic_mfs)
        report = report0 + report1 + report2 + report3 + report4 + report5 + \
            report6 + report7 + report8 + report9

        if self.parallel:
            # test_mosaic_cube_briggsbwtaper.exp_bmin_dict
            exp_bmin_dict = self.exp_dicts['exp_bmin_dict']

            # test_mosaic_cube_briggsbwtaper.exp_bmaj_dict
            exp_bmaj_dict = self.exp_dicts['exp_bmaj_dict']

            # test_mosaic_cube_briggsbwtaper.exp_pa_dict
            exp_pa_dict = self.exp_dicts['exp_pa_dict']

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

            self.save_dict_to_file(test_name,savedict, test_name+'_cur_stats')

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
            #stokes='I', specmode='cubesource', nchan=948, start=5, width=1, \
            stokes='I', specmode='cubesource', nchan=948, start='261.7643758544GHz', width='0.2441755MHz', \
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
            #stokes='I', specmode='cubesource', nchan=948, start=5, width=1, \
            stokes='I', specmode='cubesource', nchan=948, start='261.7643758544GHz', width='0.2441755MHz', \
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
        exp_im_stats = { 'com_bmaj': [False, 0.9347867427516693],
            'com_bmin': [False, 0.7185945667290505],
            'com_pa': [False, -88.16842658816135],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 105408234.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.08379194885492325],
            'max_val_pos': [True, [318, 243, 0, 511]],
            'min_val': [False, -0.07177016884088516],
            'min_val_pos': [True, [222, 212, 0, 674]],
            'im_rms': [False, 0.012333022099039984],
            'rms_per_chan': [False, [0.015614919923233097, 0.012461091353256435, 0.011097045928810322, 0.011531423920752655, 0.012262159210938556, 0.013466869238351801, 0.013850177543492121, 0.01531890407408192, 0.012541673953939087, 0.01208294453325174, 0.010653159036881963, 0.011420173447252088, 0.013334206656096588, 0.011763488390878369, 0.011796359019870296, 0.01104181931344943, 0.011184563712249615, 0.01264541551902955, 0.009626482773124367, 0.012292143836384562, 0.011794937582558012, 0.012491389004109983, 0.012069559609954495, 0.014123734188953532, 0.013010301965684817, 0.010912197489772265, 0.011340332933572663, 0.011450632210528269, 0.013869914134722798, 0.012695694937301936, 0.01235826397254351, 0.01111654248480426, 0.011724524799493456, 0.010250486304411069, 0.012017329992171404, 0.011821215246562619, 0.011990239543801743, 0.01127475496224341, 0.012624475292568882, 0.010907501697973325, 0.012299516262273739, 0.013806443683288433, 0.012011945137801348, 0.01469473390454824, 0.012670678097297744, 0.012239367529021793, 0.010336200749660414, 0.01058418240785468, 0.01065198032813839, 0.01194425004156691, 0.011960787366460288, 0.011746655358855782, 0.010858679957063438, 0.011356591626856372, 0.011500933190081115, 0.011810859539709775, 0.0109385650995274, 0.014603826115089206, 0.01296260008123089, 0.010222393584211103, 0.011094616662570878, 0.011704298759106366, 0.011607587674127054, 0.011684170862616569, 0.010725084592803056, 0.011497519476254164, 0.012040225983484837, 0.012917513730798327, 0.010426226261604447, 0.010492110796809138, 0.013588002797409654, 0.01171482183554535, 0.01090092441045949, 0.009759390302606792, 0.011660577367298764, 0.012708325390558118, 0.014300845004057672, 0.01427551347834974, 0.012620945713508254, 0.011239593988955304, 0.010568868723123064, 0.010629693050992395, 0.011342196981612934, 0.012515710113262043, 0.011315464659780127, 0.012183438662463061, 0.010636150859356899, 0.01113673257810437, 0.010150337895945319, 0.011331541003543102, 0.012032880759962835, 0.013215908777192158, 0.012332406339551951, 0.010789631774502683, 0.011564236200178515, 0.012014659468971945, 0.01112078067318156, 0.011863511040452567, 0.011928775750417956, 0.009863085252111483, 0.012838338216501731, 0.01204164825135906, 0.01302002610270098, 0.01296348344955956, 0.0117135933494222, 0.013699957713875091, 0.015258490235526237, 0.012065440245395739, 0.014396824347923802, 0.011799703652333633, 0.011046110328958482, 0.011110866207965624, 0.012838010542686492, 0.011693991871794889, 0.010796507483093768, 0.012438330715919245, 0.011270437316619008, 0.01004341405864197, 0.01209922488754825, 0.011297404754105185, 0.01191565653258805, 0.01227345556171278, 0.010008409071422076, 0.012291881125164583, 0.010881944096927384, 0.012073427204487796, 0.01048223272432688, 0.01277496205680488, 0.00989571739070131, 0.009664885590356637, 0.011310981311963913, 0.011044841759219924, 0.009971100950937789, 0.011871083805594347, 0.013810310690818138, 0.011921179820603601, 0.011347123656747888, 0.010785356379095263, 0.010868634254483712, 0.010438904036719382, 0.010599769515718388, 0.011778544812755158, 0.011557497774325467, 0.011436781128773981, 0.009949434570622728, 0.012755927531723318, 0.012396754796436255, 0.010195206143916375, 0.012170288210271362, 0.01049476507018805, 0.011531704700641686, 0.011397260559417697, 0.010827813821022391, 0.012474202957147512, 0.012916209760128076, 0.012592720635052698, 0.012538399263307123, 0.012042794123438378, 0.011499164421835947, 0.014581051369908742, 0.012795135272694077, 0.011021949511371662, 0.011275262794722726, 0.010440140139784112, 0.01191172493126342, 0.010522216836243572, 0.010570917001320774, 0.009858926706666618, 0.009868716265793715, 0.011240399216459609, 0.009973215763376798, 0.010712606272862357, 0.010809283492762728, 0.011934895591810128, 0.012738271744860991, 0.012342011017641276, 0.01015607438527581, 0.011954144785020607, 0.011919307668191532, 0.0129253527910233, 0.010485761260138974, 0.011791794054793235, 0.012426006357161115, 0.013754805273426566, 0.010546302059426738, 0.009769139500427614, 0.010824957173257795, 0.010563893431554833, 0.01237529355244673, 0.01000336732175008, 0.010717756216537026, 0.011511965995201196, 0.01373978124866005, 0.00976484090919126, 0.011766359196191671, 0.012086831973818666, 0.011204493289507516, 0.012403792502018491, 0.010839880099183944, 0.012522330102745286, 0.0127581024518738, 0.015063296219699979, 0.013349461943484098, 0.012282151689414506, 0.011892005853184083, 0.0107053855225632, 0.010955042211849907, 0.010811676264760652, 0.011364587120660621, 0.012180348512584014, 0.012565642226737285, 0.011940365703709781, 0.011437165707041659, 0.010689253649605076, 0.01103672624854216, 0.013560307081632748, 0.012392614332869048, 0.011472527677873028, 0.01160956075394297, 0.013511894386541592, 0.013619540768196627, 0.01184638237260289, 0.012299006159816514, 0.014052689891623333, 0.012740358274391734, 0.011503196535383342, 0.012207727807786888, 0.011369503379574734, 0.012331821960951147, 0.011216903654002406, 0.010327434100731244, 0.01349546402293757, 0.012456403872344948, 0.011630844499952893, 0.00933150729634402, 0.011837551660940474, 0.012622488388054375, 0.012408613386189907, 0.011238224166246244, 0.01453295594954949, 0.013227031354695117, 0.009623191707707927, 0.009966921924490758, 0.011379943771563567, 0.011291760958212651, 0.011196475489921979, 0.01148432425242337, 0.01064926468774721, 0.011289629525315011, 0.011875124331931606, 0.014826976127424363, 0.011360857649989892, 0.00984327418440123, 0.013015016844659149, 0.012519426428470828, 0.013462655435703075, 0.01209650838955224, 0.011490279316160798, 0.011573864625482044, 0.011380086626536858, 0.011717790426550249, 0.011134247406686148, 0.013163300585094856, 0.011540481536936997, 0.011684971828873195, 0.012106410337675125, 0.01266509023037908, 0.01097333509820809, 0.012283204054042427, 0.012490144327580036, 0.011410354650791004, 0.011253657707597295, 0.011583323188001718, 0.012576498013555508, 0.010706060134696882, 0.011213523482279967, 0.012686922322681293, 0.01341324314603641, 0.013185860604090317, 0.01192907660352778, 0.010844420105379724, 0.010868162341598157, 0.01113233782846976, 0.009477291430117218, 0.012141607770029811, 0.01139063055273328, 0.012178709658425024, 0.011593212208764312, 0.011736389223319927, 0.013773395086848167, 0.01046112777858449, 0.010780128697386352, 0.010783490945886749, 0.012037734305972995, 0.013077690516182131, 0.010208306820663988, 0.01202608068351915, 0.011517234853496029, 0.012571496255581894, 0.011131019532703674, 0.012941014808113266, 0.011958134415199871, 0.012997941984269392, 0.016396775691009374, 0.016509288411272135, 0.016202371255910016, 0.010597268229777938, 0.011809386296003185, 0.011092750713703557, 0.013147208863983247, 0.011334473512512608, 0.011416053885149314, 0.012930057763437846, 0.0125111195114077, 0.011424229396593207, 0.011743199749722474, 0.011628443493731369, 0.014245859704597325, 0.012166468961620677, 0.012132619531662735, 0.012501279263324474, 0.017489198820383512, 0.011693976586092558, 0.012284243162572758, 0.01179573754482819, 0.013712485841812887, 0.014366924637398424, 0.013397570085984363, 0.01122626935699209, 0.011758874008335958, 0.011562366102793254, 0.013536050834995164, 0.013173922096223895, 0.012383655334239942, 0.012358654755584024, 0.012098794663238537, 0.014693318850918367, 0.014060186282015453, 0.011664145748113608, 0.013341922581101813, 0.01194875803618567, 0.013875308800851275, 0.010538705244952524, 0.011739711602903042, 0.013346871804508319, 0.013768525643763836, 0.012108635017953106, 0.01255326867638848, 0.011017327315357641, 0.012013036170652582, 0.012563026450647668, 0.01451295089961516, 0.01319293499043089, 0.012455816864734294, 0.012107336580030355, 0.011266426491665142, 0.012921921753130632, 0.012485375443257908, 0.013307397680133596, 0.011537620435221779, 0.011810892743791785, 0.011472481031030168, 0.013843883971662791, 0.014100568893439898, 0.011931415556188072, 0.011615745473989206, 0.011414840824634702, 0.012048220475558027, 0.01301754781926097, 0.01159359427325543, 0.012667460964016306, 0.016284994192466942, 0.01432446174929441, 0.01219831434244015, 0.01177188041236153, 0.012688610603965162, 0.011848629133642631, 0.011931951190853246, 0.011336759233115409, 0.010878160267211275, 0.01141514121321475, 0.012156611526245305, 0.012858199676054236, 0.012194170575857947, 0.012591508702076419, 0.011736459933712926, 0.01285966432763003, 0.014248764996529736, 0.011803833912279169, 0.013454592743548455, 0.012348985491283385, 0.01142224994976627, 0.011551462332865834, 0.01194426012146018, 0.012118554845945002, 0.010557828249552698, 0.011746018845110007, 0.01150604915026696, 0.012609955623082439, 0.015380773105875625, 0.012934303675804438, 0.011563441082166438, 0.01252590711936631, 0.013173066740296403, 0.010888916235495698, 0.010921542827725759, 0.012025081856722998, 0.014479920335818878, 0.016195655073482838, 0.010854520873887713, 0.012494244152510567, 0.013502423258460437, 0.013546657797299204, 0.012704779945771238, 0.011847122569214155, 0.013950927503481278, 0.013663446448441944, 0.012260402176464772, 0.011984967247035348, 0.012049197075038286, 0.011514387814418063, 0.010082475305966397, 0.012014084083600067, 0.011490689399150975, 0.012931195873705965, 0.011888144849500024, 0.012328644877063984, 0.013965088359483381, 0.01140236435875301, 0.011953415582090381, 0.01249789318291057, 0.012202798820428355, 0.011511863840289998, 0.013146669887898333, 0.011407720025077764, 0.014273273582137456, 0.012343132826468014, 0.015190977135542192, 0.012371614402965227, 0.0130855996483357, 0.01252536081493026, 0.015060071650653388, 0.013974115829119853, 0.014000424556746247, 0.012760822131736566, 0.012413763331353042, 0.011762125569152591, 0.011587278867276517, 0.01192320093066336, 0.011965155907223706, 0.01519415405001481, 0.014699702014267957, 0.013054479286260016, 0.013915311873957567, 0.013796414554581235, 0.01160080761542055, 0.011609207454324061, 0.011602750098090425, 0.01345509077720275, 0.012969157722007784, 0.014011559652068866, 0.012911988660230928, 0.01234002012192669, 0.012449873779841013, 0.011165624414480699, 0.012654961348232666, 0.011314413619156066, 0.011982121494296987, 0.012337907706602386, 0.012569496575516126, 0.013540575893255948, 0.012827720470942489, 0.012828668156175636, 0.012225694812542712, 0.011926088886978181, 0.012514965514482192, 0.012885424310286216, 0.013343440007252915, 0.011790619141895895, 0.013293158572502045, 0.012504537681088583, 0.011999264759692487, 0.01242516312290984, 0.016680875173125064, 0.012837921023644361, 0.013363687131769457, 0.013979828254468127, 0.013325563148673879, 0.012932798896949486, 0.01234580007350237, 0.012409072309951356, 0.011671781593496093, 0.013419297510797526, 0.013049969215176606, 0.01118798860570238, 0.012462138826771867, 0.013053966485600674, 0.012731729062389415, 0.01579284875782042, 0.014000944967731804, 0.013260665867110093, 0.01376191318750849, 0.011482641854531357, 0.013986208023492985, 0.013343757717123712, 0.013073454019872791, 0.01203325506358, 0.013508184529822777, 0.012501337678861094, 0.012898743862438477, 0.013569797073577639, 0.015127748313051964, 0.013116500931277619, 0.011801720471515658, 0.014334588327358625, 0.014097148850424004, 0.013181284636444404, 0.013970130315954395, 0.013689278350814208, 0.01277205613613693, 0.011845527712848816, 0.01350118052709737, 0.013302490020084047, 0.013840483035346308, 0.014093284283317098, 0.012105389348034519, 0.012792304872702004, 0.01435728765353915, 0.012732787256571974, 0.013012760123425152, 0.011050518445412539, 0.012230697554863319, 0.012034128538795446, 0.011930147471993198, 0.012910748159268577, 0.013372768435942897, 0.013289201752274047, 0.012213442900684142, 0.011811366679249617, 0.013941104399910838, 0.01216823736631352, 0.012017205945075514, 0.012211985933502785, 0.012639996803006853, 0.012586001489890405, 0.013454487414786965, 0.013314632106161848, 0.0149065153080431, 0.013657671650095678, 0.01269419877384009, 0.012605432858571707, 0.013346855803716959, 0.012807869486942432, 0.015028469841526178, 0.014201772295148313, 0.01353312809545639, 0.01256640714555441, 0.01255182549841954, 0.014452440241624415, 0.011886971890394529, 0.015321557935273154, 0.014433761791561271, 0.013273987291610298, 0.013056795394492733, 0.01392481932628729, 0.01217839041697604, 0.014405003885289966, 0.014134180200381593, 0.013821525145573705, 0.01207939809364017, 0.012153125341336541, 0.012818792250357487, 0.0121266531868037, 0.012995216873470761, 0.013042912965242577, 0.011801714078314126, 0.013673486483733165, 0.012938551039197914, 0.011932268713817902, 0.0129492796439923, 0.013126596206617928, 0.0139699019228696, 0.013112440162784214, 0.014096710143933029, 0.013714920443046268, 0.012729377367478803, 0.012217757467013866, 0.01372642318449577, 0.01342268984349765, 0.0161596214859486, 0.01554502816791112, 0.01385303673656151, 0.012664029146470402, 0.01496394766353878, 0.013608648787326606, 0.013743305059356913, 0.012374529165434418, 0.013257131999116076, 0.012245019340600004, 0.012841366681682835, 0.013933656566067043, 0.014671172289196746, 0.011960107721057479, 0.012585475338646731, 0.012579049960133218, 0.012327938359368412, 0.013522115546971263, 0.011278073849888624, 0.012308769668512581, 0.012774716701909617, 0.014245521102087145, 0.013872110283528211, 0.015530797619558747, 0.012629591893160647, 0.012854490702090557, 0.013819957083595982, 0.011293915012501399, 0.013968849122342258, 0.014022336810997172, 0.014821284115710355, 0.012375904758030559, 0.014577526061499108, 0.01335507547825915, 0.013033042777863154, 0.012555524228193098, 0.012166117306323557, 0.01475247999386739, 0.013406250552779251, 0.011014654757624495, 0.01354256975102242, 0.014192444622307081, 0.011272251601273035, 0.013163246916963326, 0.011881432742270076, 0.011291550027323453, 0.013490406563131923, 0.01402643358055753, 0.011807166429960363, 0.012979620685847966, 0.012649809319629219, 0.012041400129844127, 0.014060516503510723, 0.013163368522769956, 0.015377091112262868, 0.01333483490070852, 0.01480072436742792, 0.013006090627617987, 0.01260604214324173, 0.013645567904830245, 0.01314954833409505, 0.011979851376227568, 0.014456320325964367, 0.011627633315446015, 0.013815003645075452, 0.011950662549873965, 0.011626100746086493, 0.013159809892395399, 0.012641421216083285, 0.013783195591672227, 0.011858989790475175, 0.01430188163110481, 0.013331979735069421, 0.015866224086926967, 0.013924755291249133, 0.013821567751028292, 0.01197789675212135, 0.01348702803554566, 0.01240184663677815, 0.013650731554675858, 0.012214549884137386, 0.010450766655865076, 0.010800583902999431, 0.011904444730972377, 0.01255527725016927, 0.01633109854727818, 0.011925740857073706, 0.012934469250114286, 0.011713586308362051, 0.011829318277796697, 0.015953135742004583, 0.012959953115855192, 0.01268726496516468, 0.01325976764678124, 0.010998542893165425, 0.013881986002573577, 0.012493130024531732, 0.013225457712681276, 0.011323354771460145, 0.013512163877253245, 0.011981986979333609, 0.010705414307623645, 0.015411557723988936, 0.012656907165175316, 0.01154548179423097, 0.012436027474269586, 0.012423384490044572, 0.012491492035205595, 0.013350393130031403, 0.013008007336320616, 0.013156253771299203, 0.01239698042898218, 0.012861415096467196, 0.01351158123629954, 0.013558819870434148, 0.013082256416006651, 0.012247610162583132, 0.01330383528652246, 0.012164481658242575, 0.01639354672262505, 0.013784830067551863, 0.013875118049603723, 0.013811469571797933, 0.012266945837148657, 0.01174496096315862, 0.01323540891165691, 0.011341557450151548, 0.011725546717910417, 0.010697575187416571, 0.01241672169357198, 0.013114987292226724, 0.011748110326856147, 0.01176740678294043, 0.012571080348524223, 0.011795458136474401, 0.01250439763927647, 0.01406016257805102, 0.011798710665738074, 0.014669670796195307, 0.01172508680074133, 0.011362078528283887, 0.013347465117762226, 0.013110247988178352, 0.011092010525840032, 0.010641457159832871, 0.011138978740559058, 0.011399785715907563, 0.011944862262296846, 0.012540221985187295, 0.011143136106069776, 0.010580328383115425, 0.015243185738893439, 0.012031375271904045, 0.011378008631125433, 0.010357345821665509, 0.011762147131178579, 0.01302705249292611, 0.013237503346415793, 0.012331057482587782, 0.013489034463890975, 0.015221553129395058, 0.010169070600055983, 0.011274217848748635, 0.01013087180376239, 0.011553815107164992, 0.012713975546552939, 0.012469122236709677, 0.010942563384596843, 0.01295427543429258, 0.012409741405259175, 0.012781418607651994, 0.011255937024567675, 0.01386063768329683, 0.010826526493860181, 0.013589017831079735, 0.017198748686082375, 0.011808705444454717, 0.01069955832202396, 0.012261454983523346, 0.012671935689831215, 0.012610045327698244, 0.013276501031660705, 0.01252924789427766, 0.010322350116349147, 0.012864793630010335, 0.012399919506193529, 0.010587845572973195, 0.013996707225741978, 0.011695893128146603, 0.010578926041378035, 0.010704599931148978, 0.013873710856904942, 0.012901880875233909, 0.013244781218978133, 0.011268582537805884, 0.011414600483590779, 0.012066420191098376, 0.010828856090201424, 0.011511606513911482, 0.011770231627348533, 0.011275941623222357, 0.01312325457629196, 0.010673868980572241, 0.01258879706983797, 0.01077616330363672, 0.012231260183162476, 0.01142243998089413, 0.010719861478157401, 0.014779394774239911, 0.012441962452020179, 0.011363588446218915, 0.010915871749164183, 0.011860005476236708, 0.011121688992330916, 0.0124854220594754, 0.011687095985070787, 0.01213727549975877, 0.010013069331500417, 0.01076526881702059, 0.013035709310736183, 0.011023727589520416, 0.01146800044360023, 0.010667142055079076, 0.011052350040039259, 0.010597300594320135, 0.01217219495334906, 0.011906621775783072, 0.011863753216836513, 0.01198235200593461, 0.012785110569347773, 0.011599442021411518, 0.010104137094898474, 0.010338318088836227, 0.012350156170157703, 0.011398805239719187, 0.010972510820754488, 0.010632512625657713, 0.010999042108302081, 0.012807225391046072, 0.014952117035631016, 0.012046694930262915, 0.012934051011173111, 0.013682025093869574, 0.012800065700554301, 0.010541373640307975, 0.012491246257170997, 0.011935479957140614, 0.01190553209517063, 0.010970700652845241, 0.012366032052383787, 0.011836912961207006, 0.0110789160850661, 0.012515065542790986, 0.01188937744193641, 0.013207219484098455, 0.011234207644791163, 0.011542337653849237, 0.01352324234775997, 0.010766246911482907, 0.009869694668470826, 0.010365458983649418, 0.010466495416702786, 0.012997338835245896, 0.011110917902887762, 0.01111581218294549, 0.011663967696989153, 0.010999428354828157, 0.010164479496510465, 0.011340205385996047, 0.011599590579893996, 0.010210132656324318, 0.012093974962136999, 0.012730299266437097, 0.011206337211825, 0.010797962314059065, 0.009541862423992187, 0.010571002992245213, 0.013799855965355563, 0.01245335365264682, 0.01401486067134829, 0.014105319285909638, 0.012099206460099009, 0.011832274208987553, 0.012626718936209558, 0.011570195932803106, 0.011029280695601135, 0.011649640028011504, 0.010227009574203931, 0.011037193712019255, 0.01108020404325776, 0.012051601084481427, 0.01211705380329035, 0.01214472769584763, 0.012285856955386174, 0.01136090566443648, 0.01362095861863139, 0.012907131010761672, 0.01150248176363178, 0.00995185575082863, 0.010939097781665541, 0.010521029822365217, 0.011352670714686278, 0.012144076291577018, 0.01090090350030123, 0.012305110842593099, 0.009703459977344443, 0.012310931417245924, 0.013132299523332768, 0.011212204857898578, 0.010327658297889251, 0.011595701467957142, 0.009924770034839325, 0.009990488435527281, 0.01148594999967197, 0.010726318730273729, 0.010749969816792, 0.010299898683416563, 0.010247197292435622, 0.011155271661664532, 0.01130465569780169, 0.011199358251136713, 0.01098332921748894, 0.013979778397325475, 0.010205997985730261, 0.010328856176530255, 0.010777294393515293, 0.01106168776347132, 0.011501067992044367, 0.010681141503166292, 0.015197029688680481, 0.011769035856105349, 0.00963268439107862, 0.012011342170742428, 0.009404245490578082, 0.011474873168160446, 0.01041000306802867, 0.01302494720396054, 0.01063578800339353, 0.00995731063562398, 0.012082535821067433, 0.010146329850151303, 0.010081100969426202, 0.013369097216097813, 0.01471394021989113, 0.01105192174046197, 0.01132124872974185, 0.010765414348639797, 0.011532201608984409, 0.009972607317409473, 0.010220128550545394, 0.010692531666863068, 0.011592704080436602, 0.012023593697220197, 0.011559162350337326, 0.01236048774160667, 0.012315253035990155, 0.010575619474081147, 0.01049836160893547, 0.012785385710244344]],
            'im_sum': [False, 4084.4506735376804],
            'regn_sum': [False, 33.93410069151264],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.013772266727603429, 0.012667599322017917, 0.012671293338768432, 0.013400379726227618, 0.012451951398368813, 0.013208167767346695, 0.01276561946530243]],
            'profile': [False, 0.025561013752941966]}


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
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'mask_pix': [False, 7770],
            'mask_regns': [True, 29],
            'npts_real': [True, 191116800]}


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105408234.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000003278255463],
            'im_rms': [False, 0.6370469658095314],
            'npts_0.2': [False, [111186, 111186, 111186, 111186, 111186, 111186, 111186, 111185, 111185, 111184, 111187, 111189, 111190, 111190, 111188, 111189, 111187, 111187, 111188, 111189, 111189, 111187, 111187, 111186, 111185, 111184, 111185, 111185, 111185, 111185, 111184, 111183, 111183, 111184, 111185, 111186, 111187, 111187, 111188, 111189, 111189, 111189, 111189, 111191, 111192, 111193, 111193, 111192, 111192, 111192, 111192, 111192, 111191, 111191, 111190, 111189, 111188, 111189, 111189, 111190, 111191, 111191, 111191, 111192, 111192, 111193, 111193, 111191, 111192, 111191, 111191, 111192, 111192, 111193, 111193, 111194, 111194, 111194, 111194, 111191, 111192, 111192, 111192, 111192, 111192, 111192, 111193, 111193, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111194, 111194, 111193, 111192, 111192, 111191, 111191, 111191, 111191, 111190, 111188, 111188, 111188, 111190, 111189, 111189, 111190, 111188, 111190, 111187, 111187, 111188, 111187, 111188, 111188, 111188, 111188, 111189, 111188, 111188, 111188, 111188, 111189, 111189, 111188, 111189, 111189, 111189, 111188, 111188, 111188, 111191, 111192, 111191, 111191, 111191, 111191, 111189, 111189, 111187, 111187, 111188, 111188, 111187, 111188, 111187, 111187, 111188, 111188, 111188, 111188, 111188, 111188, 111189, 111188, 111188, 111188, 111188, 111187, 111186, 111186, 111186, 111185, 111183, 111182, 111181, 111183, 111182, 111183, 111181, 111182, 111181, 111184, 111184, 111184, 111184, 111185, 111185, 111185, 111184, 111184, 111186, 111187, 111185, 111184, 111184, 111184, 111184, 111185, 111183, 111187, 111186, 111185, 111187, 111186, 111188, 111187, 111187, 111187, 111186, 111187, 111186, 111185, 111185, 111187, 111186, 111187, 111188, 111188, 111187, 111188, 111188, 111187, 111187, 111187, 111186, 111186, 111186, 111186, 111188, 111188, 111187, 111187, 111189, 111189, 111189, 111187, 111188, 111188, 111188, 111189, 111189, 111189, 111189, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111191, 111190, 111189, 111189, 111189, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111188, 111188, 111189, 111190, 111190, 111191, 111191, 111189, 111190, 111190, 111191, 111190, 111190, 111190, 111189, 111190, 111190, 111190, 111190, 111191, 111191, 111191, 111191, 111191, 111191, 111192, 111192, 111192, 111194, 111194, 111196, 111196, 111195, 111196, 111197, 111196, 111196, 111196, 111196, 111196, 111195, 111195, 111194, 111194, 111194, 111194, 111195, 111195, 111195, 111195, 111195, 111195, 111196, 111196, 111196, 111196, 111196, 111196, 111196, 111196, 111197, 111195, 111195, 111195, 111196, 111195, 111195, 111194, 111195, 111195, 111195, 111194, 111194, 111196, 111196, 111195, 111196, 111195, 111195, 111194, 111194, 111195, 111195, 111195, 111196, 111197, 111196, 111198, 111197, 111197, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111195, 111195, 111195, 111195, 111196, 111197, 111197, 111196, 111196, 111198, 111198, 111198, 111197, 111198, 111198, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111199, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111199, 111199, 111199, 111199, 111200, 111200, 111199, 111200, 111200, 111200, 111199, 111199, 111200, 111199, 111199, 111199, 111200, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111200, 111199, 111199, 111199, 111199, 111200, 111200, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111199, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111197, 111198, 111197, 111198, 111199, 111200, 111200, 111200, 111200, 111200, 111201, 111200, 111201, 111200, 111199, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111198, 111199, 111199, 111197, 111197, 111197, 111196, 111196, 111196, 111197, 111197, 111197, 111199, 111199, 111199, 111201, 111200, 111200, 111200, 111199, 111200, 111199, 111199, 111199, 111199, 111200, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111197, 111196, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111195, 111195, 111194, 111194, 111193, 111193, 111193, 111195, 111194, 111194, 111194, 111194, 111194, 111193, 111193, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111191, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111191, 111191, 111192, 111192, 111192, 111192, 111191, 111191, 111191, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111189, 111189, 111189, 111189, 111189, 111188, 111188, 111189, 111189, 111189, 111188, 111188, 111188, 111189, 111189, 111189, 111190, 111190, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111189, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111189, 111189, 111189, 111188, 111188, 111188, 111188, 111188, 111188, 111189, 111189, 111188, 111188, 111187, 111187, 111187, 111187, 111187, 111187, 111187, 111187, 111189, 111189, 111188, 111188, 111187, 111186, 111187, 111187, 111187, 111188, 111186, 111187, 111186, 111185, 111185, 111185, 111184, 111182, 111184, 111184, 111184, 111184, 111184, 111184, 111184, 111185, 111185, 111185, 111185, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111184, 111184, 111184, 111185, 111186, 111184, 111183, 111182, 111184, 111185, 111186, 111187, 111183, 111183, 111182, 111182, 111184, 111184, 111184, 111185, 111185, 111184, 111184, 111184, 111185, 111186, 111185, 111187, 111186, 111184, 111183, 111182, 111184, 111183, 111182, 111182, 111182, 111182, 111183, 111183, 111182, 111182, 111182, 111182, 111182, 111182, 111182, 111182, 111184, 111185, 111185, 111185, 111185, 111185, 111185, 111184, 111183, 111183, 111185, 111185, 111185, 111186, 111186, 111186, 111185, 111186, 111186, 111186, 111186, 111187, 111187, 111188, 111187, 111186, 111187, 111187, 111188, 111187, 111187, 111187, 111187, 111187, 111187, 111186, 111185, 111188, 111187, 111188, 111189, 111186, 111187, 111188, 111188, 111190, 111188, 111188, 111189, 111189, 111185, 111183, 111186, 111185, 111185, 111183, 111183, 111183, 111185, 111185, 111185, 111185, 111185, 111185, 111186, 111186, 111188, 111188, 111185, 111186, 111186, 111186, 111187, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111187, 111187, 111187, 111188, 111188, 111187, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111186, 111187, 111188, 111188, 111188, 111187, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111182, 111182, 111182, 111182, 111183, 111183, 111182, 111182, 111183, 111182, 111182, 111182, 111183, 111183, 111183, 111183, 111183, 111182, 111182, 111181, 111182, 111179, 111182, 111182, 111181, 111182, 111181, 111181, 111181, 111181, 111181, 111181, 111181, 111181, 111180, 111180, 111180, 111180, 111180, 111180, 111180, 111183, 111182, 111182, 111179, 111181, 111179, 111181, 111181, 111181, 111181, 111182, 111180, 111181, 111179, 111181, 111182, 111182]],
            'npts_0.5': [False, [64051, 64051, 64050, 64049, 64050, 64048, 64048, 64048, 64048, 64048, 64049, 64051, 64052, 64052, 64051, 64052, 64052, 64052, 64052, 64052, 64052, 64051, 64051, 64049, 64049, 64049, 64049, 64049, 64049, 64049, 64049, 64049, 64048, 64048, 64049, 64049, 64049, 64049, 64051, 64051, 64051, 64051, 64051, 64053, 64052, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64054, 64054, 64053, 64053, 64053, 64054, 64054, 64052, 64053, 64055, 64054, 64053, 64053, 64052, 64052, 64052, 64052, 64052, 64052, 64053, 64052, 64054, 64055, 64055, 64054, 64054, 64054, 64052, 64055, 64054, 64053, 64052, 64053, 64052, 64053, 64054, 64052, 64052, 64052, 64053, 64053, 64053, 64054, 64055, 64056, 64056, 64056, 64054, 64055, 64054, 64054, 64053, 64053, 64053, 64054, 64053, 64053, 64053, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64052, 64053, 64053, 64053, 64052, 64051, 64048, 64048, 64053, 64053, 64053, 64052, 64054, 64054, 64051, 64051, 64051, 64052, 64054, 64053, 64053, 64053, 64053, 64053, 64051, 64052, 64051, 64051, 64050, 64050, 64051, 64052, 64052, 64050, 64052, 64052, 64051, 64053, 64053, 64054, 64054, 64054, 64053, 64053, 64052, 64052, 64049, 64049, 64048, 64048, 64048, 64046, 64046, 64047, 64046, 64045, 64045, 64045, 64045, 64047, 64047, 64046, 64047, 64049, 64049, 64049, 64049, 64049, 64050, 64051, 64051, 64050, 64050, 64050, 64050, 64050, 64048, 64050, 64050, 64052, 64052, 64052, 64053, 64053, 64052, 64050, 64050, 64050, 64050, 64049, 64049, 64050, 64050, 64050, 64051, 64051, 64051, 64051, 64051, 64052, 64052, 64051, 64051, 64052, 64051, 64051, 64051, 64051, 64051, 64051, 64052, 64053, 64053, 64051, 64052, 64052, 64052, 64053, 64054, 64054, 64054, 64054, 64053, 64053, 64053, 64053, 64053, 64054, 64056, 64055, 64054, 64054, 64054, 64055, 64054, 64055, 64055, 64055, 64055, 64054, 64055, 64054, 64053, 64053, 64053, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64056, 64056, 64055, 64055, 64054, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64055, 64053, 64053, 64055, 64056, 64056, 64057, 64058, 64058, 64059, 64059, 64059, 64060, 64062, 64062, 64062, 64061, 64062, 64061, 64061, 64061, 64058, 64059, 64057, 64057, 64061, 64061, 64060, 64060, 64059, 64061, 64061, 64062, 64060, 64060, 64061, 64061, 64061, 64061, 64060, 64061, 64061, 64060, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64062, 64061, 64063, 64063, 64063, 64063, 64061, 64063, 64062, 64061, 64062, 64062, 64062, 64063, 64062, 64062, 64062, 64062, 64062, 64064, 64064, 64063, 64063, 64062, 64062, 64062, 64065, 64064, 64065, 64065, 64064, 64064, 64064, 64064, 64063, 64063, 64064, 64065, 64065, 64065, 64065, 64068, 64066, 64067, 64067, 64067, 64067, 64067, 64066, 64066, 64066, 64067, 64067, 64068, 64067, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64069, 64068, 64068, 64069, 64069, 64070, 64070, 64070, 64069, 64068, 64069, 64069, 64068, 64068, 64068, 64068, 64069, 64069, 64069, 64068, 64068, 64068, 64069, 64069, 64069, 64069, 64069, 64069, 64068, 64067, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64067, 64067, 64067, 64065, 64066, 64064, 64064, 64064, 64065, 64066, 64066, 64066, 64066, 64067, 64067, 64067, 64067, 64066, 64065, 64064, 64065, 64063, 64063, 64064, 64066, 64065, 64064, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64062, 64062, 64063, 64063, 64063, 64063, 64064, 64064, 64064, 64064, 64065, 64065, 64065, 64065, 64064, 64063, 64063, 64063, 64063, 64062, 64062, 64064, 64064, 64064, 64063, 64063, 64062, 64063, 64063, 64063, 64063, 64063, 64062, 64063, 64062, 64062, 64063, 64063, 64062, 64063, 64064, 64065, 64064, 64064, 64063, 64064, 64063, 64063, 64065, 64063, 64062, 64063, 64063, 64064, 64062, 64062, 64061, 64060, 64059, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64060, 64061, 64060, 64059, 64059, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64057, 64057, 64058, 64059, 64060, 64060, 64059, 64059, 64058, 64057, 64057, 64057, 64057, 64057, 64057, 64058, 64059, 64060, 64059, 64060, 64060, 64063, 64061, 64062, 64063, 64063, 64063, 64063, 64062, 64057, 64058, 64058, 64058, 64060, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64061, 64060, 64061, 64061, 64061, 64061, 64062, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64060, 64059, 64061, 64060, 64061, 64058, 64058, 64058, 64058, 64058, 64059, 64060, 64060, 64059, 64059, 64059, 64060, 64060, 64060, 64060, 64061, 64060, 64060, 64061, 64061, 64061, 64061, 64062, 64062, 64061, 64059, 64061, 64061, 64057, 64058, 64058, 64058, 64059, 64059, 64060, 64056, 64056, 64055, 64055, 64055, 64054, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64055, 64056, 64056, 64056, 64056, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64054, 64053, 64053, 64052, 64052, 64052, 64052, 64054, 64055, 64053, 64052, 64053, 64053, 64052, 64052, 64052, 64052, 64053, 64053, 64053, 64052, 64052, 64052, 64052, 64053, 64052, 64052, 64049, 64050, 64050, 64050, 64051, 64052, 64050, 64054, 64054, 64054, 64052, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64053, 64052, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64055, 64055, 64056, 64054, 64055, 64057, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64056, 64056, 64056, 64056, 64056, 64055, 64056, 64056, 64055, 64056, 64055, 64055, 64055, 64056, 64057, 64058, 64055, 64055, 64054, 64053, 64052, 64052, 64053, 64053, 64053, 64052, 64052, 64052, 64053, 64053, 64052, 64053, 64053, 64054, 64054, 64054, 64054, 64055, 64055, 64054, 64055, 64052, 64051, 64052, 64052, 64052, 64052, 64051, 64053, 64052, 64053, 64053, 64053, 64054, 64053, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64053, 64052, 64053, 64054, 64055, 64054, 64054, 64053, 64054, 64054, 64054, 64053, 64053, 64052, 64051, 64051, 64050, 64050, 64050, 64049, 64050, 64047, 64046, 64046, 64046, 64046, 64045, 64045, 64043, 64045, 64045, 64045, 64045, 64045, 64045, 64045, 64045, 64046, 64045, 64045, 64045, 64045, 64045, 64044, 64044, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64046, 64045, 64044, 64043, 64043, 64043, 64046, 64046, 64046, 64046, 64047, 64046, 64046, 64046, 64046, 64046, 64046]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1043170147097618, 37.05283850477761, 36.76968234028433]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.88011504139996],
            'fit_pix': [False, [240.81620655873363, 209.98247522658707]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.04449443891644478],
            'min_val_pos': [True, [250, 204, 0, 13]],
            'im_rms': [False, 0.012122599757518095],
            'im_sum': [False, 193.4848880386969],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.8853854356294094, 1.0831071465910553, 0.8486133686275924]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462002989996],
            'fit_pix_0': [False, [239.99945901481598, 209.99726416489327]],
            'fit_1': [False, [0.8855283647288508, 1.0833312996382483, 0.851031773306774]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.88011504139996],
            'fit_pix_1': [False, [239.999419222875, 209.9970611505171]],
            'fit_2': [False, [0.8850754412685042, 1.0805755675734443, 0.8486338432784019]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99561005289996],
            'fit_pix_2': [False, [239.99942734104334, 209.9969919997807]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105408234.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.07157289981842041],
            'max_val_pos': [True, [304, 256, 0, 321]],
            'min_val': [False, -0.07177016884088516],
            'min_val_pos': [True, [222, 212, 0, 674]],
            'im_rms': [False, 0.012326444809860093],
            'im_sum': [False, 3748.409866087816],
            'regn_sum': [False, 37.92781735977678],
            'npts_real': [True, 191116800]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_eph)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.03736339882016182],
            'max_val_pos': [True, [315, 229, 0, 511]],
            'min_val': [False, 0.0],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 2.3082033625559958e-05],
            'im_sum': [False, 8.653409476391971],
            'regn_sum': [False, 0.6585759028093889],
            'mask_non0': [True, 0],
            'npts_real': [True, 191116800]}


        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, 
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube_eph)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 45509.83984375],
            'max_val_pos': [True, [0, 0, 0, 0]],
            'min_val': [False, 45195.1875],
            'min_val_pos': [True, [0, 0, 0, 593]],
            'im_rms': [False, 45284.19290422809],
            'npts_real': [True, 948]}


        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_cube_eph)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.3203667402267456],
            'max_val_pos': [True, [240, 209, 0, 944]],
            'min_val': [False, 7.376294524874538e-05],
            'im_rms': [False, 0.12064853677123162],
            'im_sum': [False, 13934226.300038964],
            'npts_0.2': [False, [111186, 111186, 111186, 111186, 111186, 111186, 111186, 111185, 111185, 111184, 111187, 111189, 111190, 111190, 111188, 111189, 111187, 111187, 111188, 111189, 111189, 111187, 111187, 111186, 111185, 111184, 111185, 111185, 111185, 111185, 111184, 111183, 111183, 111184, 111185, 111186, 111187, 111187, 111188, 111189, 111189, 111189, 111189, 111191, 111192, 111193, 111193, 111192, 111192, 111192, 111192, 111192, 111191, 111191, 111190, 111189, 111188, 111189, 111189, 111190, 111191, 111191, 111191, 111192, 111192, 111193, 111193, 111191, 111192, 111191, 111191, 111192, 111192, 111193, 111193, 111194, 111194, 111194, 111194, 111191, 111192, 111192, 111192, 111192, 111192, 111192, 111193, 111193, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111194, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111195, 111194, 111194, 111193, 111192, 111192, 111191, 111191, 111191, 111191, 111190, 111188, 111188, 111188, 111190, 111189, 111189, 111190, 111188, 111190, 111187, 111187, 111188, 111187, 111188, 111188, 111188, 111188, 111189, 111188, 111188, 111188, 111188, 111189, 111189, 111188, 111189, 111189, 111189, 111188, 111188, 111188, 111191, 111192, 111191, 111191, 111191, 111191, 111189, 111189, 111187, 111187, 111188, 111188, 111187, 111188, 111187, 111187, 111188, 111188, 111188, 111188, 111188, 111188, 111189, 111188, 111188, 111188, 111188, 111187, 111186, 111186, 111186, 111185, 111183, 111182, 111181, 111183, 111182, 111183, 111181, 111182, 111181, 111184, 111184, 111184, 111184, 111185, 111185, 111185, 111184, 111184, 111186, 111187, 111185, 111184, 111184, 111184, 111184, 111185, 111183, 111187, 111186, 111185, 111187, 111186, 111188, 111187, 111187, 111187, 111186, 111187, 111186, 111185, 111185, 111187, 111186, 111187, 111188, 111188, 111187, 111188, 111188, 111187, 111187, 111187, 111186, 111186, 111186, 111186, 111188, 111188, 111187, 111187, 111189, 111189, 111189, 111187, 111188, 111188, 111188, 111189, 111189, 111189, 111189, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111191, 111190, 111189, 111189, 111189, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111188, 111188, 111189, 111190, 111190, 111191, 111191, 111189, 111190, 111190, 111191, 111190, 111190, 111190, 111189, 111190, 111190, 111190, 111190, 111191, 111191, 111191, 111191, 111191, 111191, 111192, 111192, 111192, 111194, 111194, 111196, 111196, 111195, 111196, 111197, 111196, 111196, 111196, 111196, 111196, 111195, 111195, 111194, 111194, 111194, 111194, 111195, 111195, 111195, 111195, 111195, 111195, 111196, 111196, 111196, 111196, 111196, 111196, 111196, 111196, 111197, 111195, 111195, 111195, 111196, 111195, 111195, 111194, 111195, 111195, 111195, 111194, 111194, 111196, 111196, 111195, 111196, 111195, 111195, 111194, 111194, 111195, 111195, 111195, 111196, 111197, 111196, 111198, 111197, 111197, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111195, 111195, 111195, 111195, 111196, 111197, 111197, 111196, 111196, 111198, 111198, 111198, 111197, 111198, 111198, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111199, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111199, 111199, 111199, 111199, 111200, 111200, 111199, 111200, 111200, 111200, 111199, 111199, 111200, 111199, 111199, 111199, 111200, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111200, 111199, 111199, 111199, 111199, 111200, 111200, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111199, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111197, 111198, 111197, 111198, 111199, 111200, 111200, 111200, 111200, 111200, 111201, 111200, 111201, 111200, 111199, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111198, 111199, 111199, 111197, 111197, 111197, 111196, 111196, 111196, 111197, 111197, 111197, 111199, 111199, 111199, 111201, 111200, 111200, 111200, 111199, 111200, 111199, 111199, 111199, 111199, 111200, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111197, 111196, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111195, 111195, 111194, 111194, 111193, 111193, 111193, 111195, 111194, 111194, 111194, 111194, 111194, 111193, 111193, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111191, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111192, 111191, 111191, 111192, 111192, 111192, 111192, 111191, 111191, 111191, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111189, 111189, 111189, 111189, 111189, 111188, 111188, 111189, 111189, 111189, 111188, 111188, 111188, 111189, 111189, 111189, 111190, 111190, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111189, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111190, 111189, 111189, 111189, 111188, 111188, 111188, 111188, 111188, 111188, 111189, 111189, 111188, 111188, 111187, 111187, 111187, 111187, 111187, 111187, 111187, 111187, 111189, 111189, 111188, 111188, 111187, 111186, 111187, 111187, 111187, 111188, 111186, 111187, 111186, 111185, 111185, 111185, 111184, 111182, 111184, 111184, 111184, 111184, 111184, 111184, 111184, 111185, 111185, 111185, 111185, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111184, 111184, 111184, 111185, 111186, 111184, 111183, 111182, 111184, 111185, 111186, 111187, 111183, 111183, 111182, 111182, 111184, 111184, 111184, 111185, 111185, 111184, 111184, 111184, 111185, 111186, 111185, 111187, 111186, 111184, 111183, 111182, 111184, 111183, 111182, 111182, 111182, 111182, 111183, 111183, 111182, 111182, 111182, 111182, 111182, 111182, 111182, 111182, 111184, 111185, 111185, 111185, 111185, 111185, 111185, 111184, 111183, 111183, 111185, 111185, 111185, 111186, 111186, 111186, 111185, 111186, 111186, 111186, 111186, 111187, 111187, 111188, 111187, 111186, 111187, 111187, 111188, 111187, 111187, 111187, 111187, 111187, 111187, 111186, 111185, 111188, 111187, 111188, 111189, 111186, 111187, 111188, 111188, 111190, 111188, 111188, 111189, 111189, 111185, 111183, 111186, 111185, 111185, 111183, 111183, 111183, 111185, 111185, 111185, 111185, 111185, 111185, 111186, 111186, 111188, 111188, 111185, 111186, 111186, 111186, 111187, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111186, 111187, 111187, 111187, 111188, 111188, 111187, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111186, 111187, 111188, 111188, 111188, 111187, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111185, 111182, 111182, 111182, 111182, 111183, 111183, 111182, 111182, 111183, 111182, 111182, 111182, 111183, 111183, 111183, 111183, 111183, 111182, 111182, 111181, 111182, 111179, 111182, 111182, 111181, 111182, 111181, 111181, 111181, 111181, 111181, 111181, 111181, 111181, 111180, 111180, 111180, 111180, 111180, 111180, 111180, 111183, 111182, 111182, 111179, 111181, 111179, 111181, 111181, 111181, 111181, 111182, 111180, 111181, 111179, 111181, 111182, 111182]],
            'npts_0.5': [False, [64051, 64051, 64050, 64049, 64050, 64048, 64048, 64048, 64048, 64048, 64049, 64051, 64052, 64052, 64051, 64052, 64052, 64052, 64052, 64052, 64052, 64051, 64051, 64049, 64049, 64049, 64049, 64049, 64049, 64049, 64049, 64049, 64048, 64048, 64049, 64049, 64049, 64049, 64051, 64051, 64051, 64051, 64051, 64053, 64052, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64054, 64054, 64053, 64053, 64053, 64054, 64054, 64052, 64053, 64055, 64054, 64053, 64053, 64052, 64052, 64052, 64052, 64052, 64052, 64053, 64052, 64054, 64055, 64055, 64054, 64054, 64054, 64052, 64055, 64054, 64053, 64052, 64053, 64052, 64053, 64054, 64052, 64052, 64052, 64053, 64053, 64053, 64054, 64055, 64056, 64056, 64056, 64054, 64055, 64054, 64054, 64053, 64053, 64053, 64054, 64053, 64053, 64053, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64052, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64051, 64052, 64053, 64053, 64053, 64052, 64051, 64048, 64048, 64053, 64053, 64053, 64052, 64054, 64054, 64051, 64051, 64051, 64052, 64054, 64053, 64053, 64053, 64053, 64053, 64051, 64052, 64051, 64051, 64050, 64050, 64051, 64052, 64052, 64050, 64052, 64052, 64051, 64053, 64053, 64054, 64054, 64054, 64053, 64053, 64052, 64052, 64049, 64049, 64048, 64048, 64048, 64046, 64046, 64047, 64046, 64045, 64045, 64045, 64045, 64047, 64047, 64046, 64047, 64049, 64049, 64049, 64049, 64049, 64050, 64051, 64051, 64050, 64050, 64050, 64050, 64050, 64048, 64050, 64050, 64052, 64052, 64052, 64053, 64053, 64052, 64050, 64050, 64050, 64050, 64049, 64049, 64050, 64050, 64050, 64051, 64051, 64051, 64051, 64051, 64052, 64052, 64051, 64051, 64052, 64051, 64051, 64051, 64051, 64051, 64051, 64052, 64053, 64053, 64051, 64052, 64052, 64052, 64053, 64054, 64054, 64054, 64054, 64053, 64053, 64053, 64053, 64053, 64054, 64056, 64055, 64054, 64054, 64054, 64055, 64054, 64055, 64055, 64055, 64055, 64054, 64055, 64054, 64053, 64053, 64053, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64056, 64056, 64055, 64055, 64054, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64055, 64053, 64053, 64055, 64056, 64056, 64057, 64058, 64058, 64059, 64059, 64059, 64060, 64062, 64062, 64062, 64061, 64062, 64061, 64061, 64061, 64058, 64059, 64057, 64057, 64061, 64061, 64060, 64060, 64059, 64061, 64061, 64062, 64060, 64060, 64061, 64061, 64061, 64061, 64060, 64061, 64061, 64060, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64062, 64061, 64063, 64063, 64063, 64063, 64061, 64063, 64062, 64061, 64062, 64062, 64062, 64063, 64062, 64062, 64062, 64062, 64062, 64064, 64064, 64063, 64063, 64062, 64062, 64062, 64065, 64064, 64065, 64065, 64064, 64064, 64064, 64064, 64063, 64063, 64064, 64065, 64065, 64065, 64065, 64068, 64066, 64067, 64067, 64067, 64067, 64067, 64066, 64066, 64066, 64067, 64067, 64068, 64067, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64069, 64068, 64068, 64069, 64069, 64070, 64070, 64070, 64069, 64068, 64069, 64069, 64068, 64068, 64068, 64068, 64069, 64069, 64069, 64068, 64068, 64068, 64069, 64069, 64069, 64069, 64069, 64069, 64068, 64067, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64068, 64067, 64067, 64067, 64065, 64066, 64064, 64064, 64064, 64065, 64066, 64066, 64066, 64066, 64067, 64067, 64067, 64067, 64066, 64065, 64064, 64065, 64063, 64063, 64064, 64066, 64065, 64064, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64062, 64062, 64063, 64063, 64063, 64063, 64064, 64064, 64064, 64064, 64065, 64065, 64065, 64065, 64064, 64063, 64063, 64063, 64063, 64062, 64062, 64064, 64064, 64064, 64063, 64063, 64062, 64063, 64063, 64063, 64063, 64063, 64062, 64063, 64062, 64062, 64063, 64063, 64062, 64063, 64064, 64065, 64064, 64064, 64063, 64064, 64063, 64063, 64065, 64063, 64062, 64063, 64063, 64064, 64062, 64062, 64061, 64060, 64059, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64060, 64061, 64060, 64059, 64059, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64057, 64057, 64058, 64059, 64060, 64060, 64059, 64059, 64058, 64057, 64057, 64057, 64057, 64057, 64057, 64058, 64059, 64060, 64059, 64060, 64060, 64063, 64061, 64062, 64063, 64063, 64063, 64063, 64062, 64057, 64058, 64058, 64058, 64060, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64061, 64060, 64061, 64061, 64061, 64061, 64062, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64060, 64059, 64061, 64060, 64061, 64058, 64058, 64058, 64058, 64058, 64059, 64060, 64060, 64059, 64059, 64059, 64060, 64060, 64060, 64060, 64061, 64060, 64060, 64061, 64061, 64061, 64061, 64062, 64062, 64061, 64059, 64061, 64061, 64057, 64058, 64058, 64058, 64059, 64059, 64060, 64056, 64056, 64055, 64055, 64055, 64054, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64055, 64056, 64056, 64056, 64056, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64054, 64053, 64053, 64052, 64052, 64052, 64052, 64054, 64055, 64053, 64052, 64053, 64053, 64052, 64052, 64052, 64052, 64053, 64053, 64053, 64052, 64052, 64052, 64052, 64053, 64052, 64052, 64049, 64050, 64050, 64050, 64051, 64052, 64050, 64054, 64054, 64054, 64052, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64053, 64052, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64055, 64055, 64056, 64054, 64055, 64057, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64056, 64056, 64056, 64056, 64056, 64055, 64056, 64056, 64055, 64056, 64055, 64055, 64055, 64056, 64057, 64058, 64055, 64055, 64054, 64053, 64052, 64052, 64053, 64053, 64053, 64052, 64052, 64052, 64053, 64053, 64052, 64053, 64053, 64054, 64054, 64054, 64054, 64055, 64055, 64054, 64055, 64052, 64051, 64052, 64052, 64052, 64052, 64051, 64053, 64052, 64053, 64053, 64053, 64054, 64053, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64054, 64053, 64052, 64053, 64054, 64055, 64054, 64054, 64053, 64054, 64054, 64054, 64053, 64053, 64052, 64051, 64051, 64050, 64050, 64050, 64049, 64050, 64047, 64046, 64046, 64046, 64046, 64045, 64045, 64043, 64045, 64045, 64045, 64045, 64045, 64045, 64045, 64045, 64046, 64045, 64045, 64045, 64045, 64045, 64044, 64044, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64043, 64046, 64045, 64044, 64043, 64043, 64043, 64046, 64046, 64046, 64046, 64047, 64046, 64046, 64046, 64046, 64046, 64046]],
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

            self.save_dict_to_file(test_name,savedict, test_name+'_stats')
                              

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
            stokes='I', specmode='cubesource', nchan=948, start='261.7643758544GHz', width='0.2441755MHz', \
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
            stokes='I', specmode='cubesource', nchan=948, start='261.7643758544GHz', width='0.2441755MHz', \
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
        exp_im_stats = {'com_bmaj': [False, 0.9379539996102209],
            'com_bmin': [False, 0.7149098922174161],
            'com_pa': [False, -87.11030298442445],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 105415430.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.07559825479984283],
            'max_val_pos': [True, [319, 243, 0, 511]],
            'min_val': [False, -0.06556452810764313],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.011414422783139517],
            'rms_per_chan': [False, [0.013952824376742565, 0.01122374880188617, 0.01036553044594618, 0.010834329676841383, 0.011004250469721491, 0.012220673723321428, 0.012688584549121765, 0.014117962775334711, 0.011720825390122841, 0.01150046125293291, 0.0103864298965176, 0.010206043663966428, 0.012697296801477345, 0.010723566217492262, 0.010225173519223276, 0.009987127133368051, 0.010988351777900801, 0.012634244424198785, 0.008807171800766544, 0.011343850092221799, 0.011825465580970744, 0.01197238247131425, 0.011548195845904923, 0.012509220973762565, 0.011430761938164542, 0.010516294928906189, 0.010556053195249375, 0.010772505499515948, 0.013047874506385932, 0.012144276984028798, 0.010973526429491988, 0.010425155535221001, 0.010592518065153513, 0.009873344777935559, 0.011324096846053648, 0.011235527104471578, 0.010657006518532838, 0.01118696673081834, 0.010441376510623367, 0.010718751435464946, 0.011016019138439663, 0.011685873396758147, 0.011051296889981426, 0.012624531504747455, 0.011495122666352602, 0.010940682875714466, 0.010264885611006013, 0.009859712760843112, 0.010092875748068962, 0.010463782771704004, 0.011300096152831484, 0.011764161175840403, 0.010868168208368192, 0.010959181863245591, 0.010741709188773355, 0.011042677808611756, 0.010660096075986097, 0.014365283533556596, 0.011581756514527888, 0.009245547699758295, 0.009578112722034621, 0.01177349099308397, 0.010156872298378177, 0.010586837316931719, 0.010267663688073065, 0.010587398789755441, 0.011002086760542587, 0.011322845802284855, 0.010571802888845194, 0.009463492687950848, 0.012485991884674745, 0.011654286251208392, 0.010680510333091869, 0.009165114959199092, 0.01055379379540136, 0.012211516574827405, 0.012899995382933438, 0.01268651898262829, 0.012104628004125746, 0.011480697040340734, 0.010351450276380667, 0.010231297135273263, 0.010734718225700834, 0.011480976783324262, 0.010895105705188176, 0.010743715697264495, 0.01046163065539561, 0.010689124300977711, 0.009605882316073025, 0.010850892214708589, 0.011044172446603813, 0.011845893888973565, 0.010956729200407207, 0.009834190529333754, 0.010215185228844547, 0.011530699533545206, 0.011497232472462158, 0.011260029930689159, 0.011017335334340248, 0.009052852724625913, 0.01119258022375134, 0.01200730006025853, 0.012821244509307748, 0.011750440952480788, 0.012067949554595897, 0.01226865939139761, 0.013390593977009245, 0.011391954360365142, 0.01247553014871181, 0.011618457146484221, 0.01057103943451045, 0.010503603547144014, 0.011489545210526548, 0.010639696716389365, 0.01023774144720609, 0.011612913031872072, 0.010602298218108337, 0.010061536507327351, 0.011769640230512876, 0.01114513815176509, 0.009786145266547017, 0.011434759316376563, 0.009614645658584291, 0.0114404042347772, 0.010927713616333915, 0.011406605456654988, 0.010158697258747752, 0.01255648701524223, 0.010390953474902263, 0.009205247483587689, 0.010740889583440584, 0.01040714290630887, 0.00968729070063426, 0.011313585674567998, 0.012295041764605013, 0.010833269101636784, 0.010924049173667369, 0.010763327766352313, 0.010075573878892576, 0.010604514669851524, 0.011153329945884919, 0.010376795145452538, 0.009854785441170289, 0.0101971656766406, 0.008997858698050323, 0.011436041707576339, 0.011157885032837557, 0.009536333234454228, 0.011448591877354017, 0.009761545659605034, 0.010563502414979215, 0.011095646056926442, 0.010593863499552096, 0.010441045391313615, 0.011560380554809626, 0.01082055052946384, 0.011234294149503197, 0.01099869419649269, 0.01036397177069392, 0.013644944544639317, 0.013382997596917573, 0.010835862672442726, 0.009689561225532016, 0.00976430897619941, 0.010417358945731101, 0.009921349092576838, 0.01009193186517613, 0.009139827201116807, 0.009689656205802224, 0.01058547758238007, 0.009391042137992274, 0.010502715087721305, 0.010210266318824265, 0.010724962830436226, 0.012003576561693175, 0.011536067127347383, 0.009624201594632556, 0.011907981244120014, 0.01208076067864276, 0.01183959356655806, 0.010067652438190945, 0.011519835037815833, 0.010954243238264383, 0.012263139383990646, 0.010550259999188079, 0.010250802276482318, 0.011425051094736896, 0.00994842754266715, 0.012159159440724392, 0.010318720323439945, 0.01075689522332167, 0.01031989168014135, 0.011904558939694629, 0.0100317384952476, 0.010628759603048077, 0.01182366591474211, 0.010186582466130276, 0.011700025558958868, 0.01021404960812202, 0.011417438717531925, 0.011999024044278746, 0.01392762086861483, 0.012489437877448223, 0.011034023085982822, 0.01022441315910995, 0.009396806486993136, 0.010049603518199368, 0.010320136400271077, 0.011342298478993896, 0.010956877260669809, 0.011573804456569066, 0.011928855322364727, 0.011162369868918283, 0.010935953225793915, 0.009907023038469125, 0.01213190324208271, 0.011718594894909415, 0.010252741062221396, 0.010697628942382204, 0.012540385481585979, 0.012640707881809363, 0.011064628716590594, 0.011217137005439349, 0.012858901784057773, 0.011865604033374333, 0.010016402134241093, 0.011572786799695743, 0.011180940347114114, 0.01158691506680587, 0.010452168578873398, 0.010015013828739951, 0.012753134186786953, 0.01159742264816528, 0.011140479953007156, 0.009273210113972667, 0.012009152372669881, 0.011597297245766264, 0.011183201585528085, 0.010794712941214729, 0.013432319035404174, 0.011756044963257346, 0.009235916365296604, 0.009550585493158441, 0.010935621577171694, 0.011600852270012972, 0.010171726146430039, 0.010538613417399444, 0.010585259851346408, 0.011328465214403987, 0.011390138403294425, 0.012726058409170848, 0.010431229107882528, 0.009385907035474865, 0.011679787689607923, 0.011588037766780626, 0.012407930784466044, 0.011927334968692697, 0.010270326218799697, 0.010762009157886148, 0.010443602919455827, 0.011123927960244459, 0.01111173018605914, 0.012526600435553984, 0.010808133575665423, 0.011028102659751184, 0.01165680567107892, 0.012096511652251894, 0.010116692211530902, 0.011620223390640536, 0.012331484927412423, 0.010904468913089983, 0.01042956414814794, 0.010706197874821998, 0.011418015034032086, 0.010482062600990922, 0.010702579987894771, 0.012040542417747295, 0.012623704490704539, 0.012404927931175506, 0.01112225502450896, 0.01011890268647475, 0.010821304490882117, 0.011083657944610917, 0.008712699835865027, 0.011142041586668006, 0.011548177387988239, 0.011718402804628053, 0.01146046736782539, 0.011671415860122861, 0.012647959678755254, 0.010111503230536086, 0.010180311952794301, 0.010989613650405006, 0.010980958956414488, 0.01168182312104016, 0.009667052466366928, 0.010663868535818815, 0.01039096710745729, 0.011330755861777427, 0.010146449288215585, 0.012720700984412822, 0.011257444018065938, 0.011833725219232955, 0.015105511275514552, 0.015118315644340102, 0.015026459483228963, 0.010114301662487495, 0.011529798000829197, 0.010603865540796273, 0.013059156339313299, 0.010590663785158584, 0.010836187379582266, 0.01194389040559194, 0.012333781255336603, 0.011851660073013687, 0.011114826717482663, 0.011284888399887987, 0.013798730167180386, 0.01151796566135432, 0.011390679135671104, 0.011291026609377926, 0.0151016486078002, 0.011218273687550894, 0.011125378057198662, 0.01078828584128242, 0.01275516599677738, 0.013629574909518171, 0.011652715077449432, 0.010227344672666882, 0.011246324727313413, 0.010722009871973375, 0.012610173685585, 0.0111171322032092, 0.011752633046504197, 0.011583504981331285, 0.01033974781572764, 0.013460023637319948, 0.012198865352598143, 0.010906560959312252, 0.013254593231151739, 0.01044612581889771, 0.012918388216006268, 0.00949104883405205, 0.01068390726986989, 0.011543287148402682, 0.013097999775280362, 0.010319709504762952, 0.011721139777044602, 0.01002783253581918, 0.01059992482364449, 0.011687465298611938, 0.01244383424699979, 0.012111198028446965, 0.011274718625194424, 0.011229473429933854, 0.01047269603415, 0.01233985006636969, 0.011713563941315284, 0.012235472468679865, 0.011718276624643162, 0.010930794717324438, 0.011378132711342118, 0.013138607613098416, 0.01348495717493486, 0.010435615926960233, 0.011080552522487824, 0.010471769648441117, 0.011589204942684912, 0.01263645875447981, 0.010805760729620164, 0.011668012860705075, 0.012986674140417608, 0.011996709205242533, 0.010530305013564405, 0.012044058292505579, 0.011309029468878441, 0.010812952408266463, 0.011519752425435414, 0.011279202371219486, 0.009632644936531945, 0.010495411553353356, 0.010455028321820263, 0.012238790713492992, 0.010516879654794713, 0.011788515579288837, 0.010610564631407079, 0.012077602742788872, 0.012967894960484744, 0.01051541051095021, 0.012992835871718348, 0.01203428202492438, 0.011097653840889152, 0.01058271081247107, 0.01134733942085389, 0.01118864367682186, 0.009901984490108116, 0.009961053218184702, 0.011142325350306874, 0.011530561488577989, 0.014816798702136317, 0.011400772826957135, 0.011297171049487882, 0.011260085280287836, 0.012349974811933005, 0.011064385685016517, 0.009779353812497203, 0.010976599549148576, 0.013412113014904907, 0.014811874210124549, 0.010844948073749978, 0.011562309636126056, 0.01144612214350528, 0.012561865898821146, 0.01170057574854564, 0.011207345205593425, 0.013079226785468284, 0.012528012528695473, 0.011131949337491873, 0.011544103067186182, 0.010936885538351553, 0.01062273403740401, 0.00901337583106144, 0.010403658781936062, 0.010294436184900314, 0.011338694153499038, 0.011349069992923573, 0.01141893104923401, 0.012404657316721204, 0.010871086097061848, 0.011033323261861766, 0.011390063009330985, 0.011161397895961779, 0.011484704832147597, 0.011403322499188722, 0.010527221499102419, 0.013264762257912547, 0.011434864198259903, 0.013350197401572192, 0.01191736948705002, 0.011744646748731543, 0.011563651219045252, 0.01384967807920074, 0.012943141634616684, 0.012081354045684976, 0.010685649498294232, 0.011261047939547702, 0.011443364497984729, 0.011119637029170248, 0.010863576229039652, 0.011518543007453544, 0.013327302500154575, 0.013326317274475755, 0.011599823241672743, 0.012609385571059924, 0.011720062202005848, 0.010518348239698307, 0.01076332764412709, 0.010031070545281453, 0.012444781300242718, 0.011202314354891581, 0.012748990606449675, 0.011265102148708748, 0.011322592126954056, 0.011477081689336463, 0.010030643468384396, 0.011712640203207936, 0.010690866527092967, 0.01075749040585441, 0.011549709689727983, 0.011187435850903019, 0.012484362273476166, 0.010802199738081946, 0.011365566185466427, 0.010825184897015878, 0.011106979259778382, 0.011543014646103699, 0.011384475352378758, 0.012823588447738373, 0.01096561583839317, 0.011016618210579128, 0.011529829973903068, 0.011236063612136481, 0.011392082368896503, 0.015308664215777578, 0.011327895879957984, 0.012049903696589988, 0.012078936136690943, 0.012003917218173767, 0.01162545164947221, 0.011800910112129283, 0.011989563129743613, 0.010885169003754957, 0.011655885494970849, 0.012368874723909693, 0.010278287498393112, 0.01172492052787343, 0.011678136813822765, 0.011490550950568857, 0.013297366186143666, 0.01216861433668289, 0.01197353346932699, 0.013456311087519322, 0.010473834757480977, 0.012495739982716173, 0.011819050263925646, 0.011971935409738265, 0.010799062880865308, 0.01266347685085823, 0.011800454378469702, 0.011429091125423651, 0.012570084094665483, 0.012998338259497063, 0.012373303428184898, 0.010634443715323157, 0.012613394393797954, 0.012976320130042658, 0.012428751712765147, 0.01374874456056666, 0.013067855423954215, 0.012032385074084728, 0.010182120568787366, 0.013087171483299914, 0.012230477044724044, 0.012156401113828052, 0.012436269025058011, 0.010005210513068687, 0.011325856937439762, 0.01314433897760276, 0.011214165778840093, 0.011853402092926195, 0.010864412862343397, 0.011218941486350994, 0.010899820788606058, 0.011244429446140786, 0.01204535972610005, 0.012073248954623734, 0.012997594140235477, 0.010811100568815062, 0.011034492782254131, 0.011961245724580969, 0.011906288770382517, 0.011134464361566999, 0.011412121510880805, 0.011860156359756823, 0.011518458213456161, 0.011628286110020666, 0.010811581710720258, 0.012430897840000083, 0.011748620204369495, 0.011535229089948314, 0.011629942642785301, 0.011430926606855367, 0.011766634611400033, 0.013550712192092787, 0.012402173122848893, 0.012151904012272523, 0.011126919729675435, 0.011517127284985114, 0.012271514191539275, 0.01056083272848389, 0.012274683541053827, 0.013587642343757938, 0.012210558363567154, 0.01148488347817356, 0.012400496071788409, 0.010838966631267082, 0.012388877439061637, 0.013878707862703451, 0.012275979487678475, 0.011323360170231035, 0.01121870743205782, 0.011178528607032612, 0.011207664238630385, 0.012054333983076697, 0.012740801930889315, 0.010516502337647352, 0.011626373597532264, 0.01158806222576651, 0.011056489775324066, 0.012377951379575826, 0.012368655687074655, 0.013495557379710305, 0.011590472024485267, 0.01245797259670051, 0.013297881181228492, 0.01118529922651134, 0.01057965878434858, 0.012276300760152016, 0.012150321205236144, 0.013953105219985931, 0.014114436179135981, 0.01274382436345107, 0.012281081596996288, 0.013800972680021822, 0.012582517975228314, 0.012181271034299975, 0.011054475605472795, 0.012467326636472894, 0.011057959684919363, 0.011190697634104614, 0.0137403938765662, 0.013173168513310267, 0.010640692934939558, 0.01140539760260584, 0.011389456985471544, 0.01132726374258897, 0.012243194919569279, 0.010209744383518899, 0.010958135798035274, 0.011326545992192997, 0.013465921934065891, 0.012661114030576194, 0.01463229500134103, 0.011356815263979136, 0.0107618210833855, 0.013140946362914829, 0.01083690824830264, 0.013172525764170316, 0.01292701173624945, 0.012958032437609193, 0.010927424833765095, 0.01238002391455297, 0.012142374844486582, 0.012637517586724581, 0.012711500578420475, 0.01117716456070126, 0.013136558203418355, 0.012166307718655997, 0.010563785045088299, 0.012986542492546714, 0.013534043141111252, 0.010712483039058413, 0.012677800115234191, 0.010963257026321243, 0.010218493915542706, 0.013285863198010981, 0.011529610352286326, 0.010752763945362082, 0.011788338556299339, 0.012031806526078178, 0.011693360069662784, 0.012676543716582448, 0.011506670316679127, 0.013149073723102403, 0.012085333675123343, 0.013376491635193005, 0.01111320533604476, 0.010394468117641449, 0.012714584790161415, 0.011906921840545957, 0.01071106827367555, 0.012895878906904017, 0.010881479211369146, 0.012486328813125605, 0.011146752706472468, 0.011403011641073426, 0.012251630810213852, 0.011575987321486422, 0.012300520814192719, 0.011415929167289156, 0.011997599859092881, 0.012336155809177018, 0.0139986705069243, 0.012833653299561033, 0.013243530521819471, 0.011215134946904464, 0.01176907906840581, 0.011247300114259085, 0.012333157771804045, 0.011234534581049166, 0.009636187340416819, 0.010337555453741931, 0.01121452718373124, 0.011663824818025937, 0.014598761444361974, 0.010997415372723404, 0.011772664997746884, 0.011143206745554295, 0.010943425375547907, 0.014809790324021516, 0.0113131848874912, 0.010948777676244364, 0.011848497737600194, 0.010206965249226487, 0.012438666962962424, 0.011689366845999435, 0.011367540627353592, 0.010550650207427622, 0.012253800965697276, 0.011065662593875899, 0.009844975474648092, 0.012576256336498502, 0.01050060341115817, 0.010518522196089118, 0.011282224062765052, 0.012342816323913814, 0.012105562158827054, 0.012263519497937879, 0.011713503695592615, 0.011603600888740944, 0.011992970886539422, 0.011719924591044012, 0.01201361392607509, 0.013025470662629947, 0.0117075264274601, 0.011505856431802979, 0.012139517405937682, 0.010639835381092123, 0.013828392273734275, 0.012473986239082817, 0.012327392574340288, 0.013050214513014611, 0.010582642818901497, 0.010776919707307254, 0.0116128015163686, 0.010435037129990781, 0.011244445807191592, 0.0099457334595068, 0.011245627999041989, 0.011834750409082755, 0.010308843560090767, 0.011923451331761755, 0.012493549377989778, 0.010682440258649381, 0.010475041876924876, 0.01317994741363128, 0.010997139474403644, 0.01292285205044818, 0.010765195874034492, 0.010913278149854792, 0.012149965933537166, 0.011995166203440153, 0.009968868966760918, 0.009859288965251498, 0.009992233063572125, 0.010956240696860069, 0.010628052967280605, 0.012317368384022316, 0.010460242089426683, 0.010377812695525412, 0.0131823239907416, 0.011530812230098518, 0.010455138348376254, 0.009892527849427662, 0.010644347564212557, 0.011984107131278727, 0.01161714570647872, 0.010295509328218053, 0.01171774595258621, 0.014015519055174171, 0.01008130349691525, 0.010210209181435262, 0.009832343555772132, 0.010925458498669043, 0.012718507446211585, 0.011401807428126319, 0.010713287780967865, 0.011393710166341409, 0.011558211436860383, 0.011389539073754604, 0.01041868684277482, 0.014140997935580683, 0.011294804811139486, 0.012644567896654989, 0.014813716037069739, 0.010789204528846867, 0.010523220336355128, 0.012478677433699552, 0.012244910961322112, 0.012360939999442735, 0.011962567085647231, 0.010895040295538411, 0.010352099603701948, 0.01209520015504535, 0.011284333249016745, 0.010504505528710007, 0.011668335946883414, 0.010719416648893577, 0.010662526466092717, 0.010292810347775934, 0.01301544988344848, 0.011229669367442801, 0.012197035018980661, 0.010967359979419512, 0.012088814136393814, 0.012150552411312817, 0.011163207511430446, 0.011061454822150218, 0.012201979837952117, 0.012222140188659438, 0.012669747925406425, 0.009683205538679057, 0.010849905424313843, 0.010008190129822389, 0.010788218346881021, 0.010544110854951673, 0.009497815751181663, 0.013007702174811812, 0.0104723161185551, 0.010798189994942635, 0.01008724934180505, 0.01112061787476007, 0.010494130319224035, 0.010976832275534653, 0.010425632772984762, 0.01193494215628546, 0.010064404866493147, 0.010697646614147852, 0.012418691210891454, 0.011050299565400795, 0.011156427053422657, 0.009263749738842144, 0.011387010385967533, 0.009739884933126648, 0.011092096800857884, 0.011107119395718125, 0.011371437692783724, 0.011355922727726625, 0.012768299897188085, 0.010899173231511367, 0.00960979945669278, 0.010215326580891389, 0.011994841053822135, 0.010082402116272409, 0.010247472871244966, 0.009642791724968364, 0.010367701251711943, 0.012379685699924233, 0.01375774561833414, 0.011070966320980764, 0.011349646202896645, 0.012564226263365142, 0.011588076593878954, 0.010865209525281283, 0.011305972946200745, 0.011448874782984864, 0.011602458570342082, 0.011049277766878713, 0.011270314577347893, 0.010660844000860219, 0.010383828212795732, 0.01212735617945348, 0.010892251540137132, 0.012363704904018879, 0.010459711231242265, 0.01054082892028388, 0.011291045510285912, 0.009463284756444574, 0.009329839249935207, 0.009969865369152661, 0.010374806582412912, 0.01123501605855846, 0.010016857513744605, 0.010995925428895407, 0.01120359185953082, 0.011121279177505237, 0.010041217357694513, 0.010707556498027581, 0.010347100894240534, 0.009660242858766533, 0.01094970071363911, 0.012393773657568496, 0.010516582177694446, 0.009837924440799128, 0.009327982355594285, 0.009762716585448519, 0.012094640829319572, 0.011448380114117775, 0.012212649868948936, 0.01229935958931043, 0.010949366661871488, 0.011090996550054, 0.0120483727064869, 0.011101058406991167, 0.009727322779502374, 0.00999661182447058, 0.009606654563035886, 0.010451457874179315, 0.010167989902002492, 0.01159513049832106, 0.011704691982710753, 0.011067250236874238, 0.011456548903758908, 0.01090104153385783, 0.013002933845534278, 0.012005079216744558, 0.011127458366756866, 0.009188546227543112, 0.010283444687703832, 0.009793711779041458, 0.009968285415044438, 0.010884327221573861, 0.010755688980940557, 0.011036798871365118, 0.009093402264884457, 0.01168164686788081, 0.01223254870607638, 0.01016291888778539, 0.01046307256572417, 0.010359609849231118, 0.009858885682224818, 0.009324334297760346, 0.01048927915844469, 0.010169812427970357, 0.009741990950295076, 0.010053313483894844, 0.009469599185712257, 0.010928813954518774, 0.010561699900345375, 0.01060824740290175, 0.01091379261080372, 0.013100180646622998, 0.010264553846013089, 0.009922009564640219, 0.01067419433475242, 0.010664265783322505, 0.011049101042344456, 0.010870116384194125, 0.013836066487613348, 0.01119181948826399, 0.00970501584862202, 0.010633150282424128, 0.00949470051469533, 0.011954062229303616, 0.010388103528797562, 0.012614311670728257, 0.01061951345866921, 0.009254343841052556, 0.01120456170140556, 0.009956015367818307, 0.009268307974762717, 0.01268117123408987, 0.012851235046575264, 0.010429532569432114, 0.01004488221819719, 0.009782709376884047, 0.01047266077040283, 0.01022848299796577, 0.009847378607697905, 0.010181037936314854, 0.010692860238232162, 0.011229493312116914, 0.010593836713620335, 0.01198429026555785, 0.012437950879111053, 0.010172811077114883, 0.00998401554452982, 0.011244164533793096]],
            'im_sum': [False, 9593.732334884902],
            'regn_sum': [False, 26.102949512250234],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.012651058191565394, 0.011590887207982466, 0.011895721346497457, 0.011957665308940704, 0.011554507601404538, 0.012469100555107608, 0.01181339872829841]],
            'profile': [False, 0.022887460846090554]}

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
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'mask_pix': [False, 8557],
            'mask_regns': [True, 31],
            'npts_real': [True, 191116800]}


        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_pcwdT.exp_pb_stats
        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415430.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000001788139343],
            'im_rms': [False, 0.6370611269211811],
            'npts_0.2': [False, [111207, 111211, 111208, 111208, 111206, 111204, 111203, 111206, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111207, 111206, 111206, 111207, 111206, 111205, 111205, 111206, 111205, 111205, 111207, 111206, 111206, 111206, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111204, 111206, 111204, 111206, 111206, 111206, 111206, 111208, 111209, 111207, 111207, 111208, 111207, 111207, 111207, 111207, 111208, 111207, 111207, 111208, 111209, 111209, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111206, 111206, 111205, 111207, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111209, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111208, 111207, 111205, 111205, 111206, 111206, 111205, 111204, 111206, 111206, 111205, 111205, 111206, 111206, 111205, 111206, 111204, 111206, 111206, 111204, 111203, 111204, 111204, 111203, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111203, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111197, 111198, 111199, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111206, 111205, 111206, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111206, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111203, 111203, 111204, 111202, 111202, 111203, 111202, 111202, 111200, 111203, 111201, 111200, 111200, 111201, 111201, 111201, 111196, 111198, 111196, 111196, 111196, 111197, 111196, 111196, 111196, 111196, 111197, 111197, 111198, 111197, 111199, 111199, 111199, 111198, 111200, 111202, 111200, 111200, 111200, 111201, 111201, 111200, 111199, 111200, 111200, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111200, 111198, 111198, 111198, 111201, 111199, 111199, 111199, 111199, 111199, 111201, 111201, 111200, 111201, 111199, 111199, 111201, 111200, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111202, 111201, 111204, 111204, 111202, 111199, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111204, 111203, 111201, 111199, 111201, 111202, 111203, 111202, 111203, 111205, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111203, 111202, 111201, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111201, 111202, 111198, 111198, 111201, 111199, 111200, 111199, 111200, 111197, 111193, 111194, 111194, 111194, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111199, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111203, 111202, 111203, 111203, 111203, 111204, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111195, 111199, 111199, 111197, 111197, 111196, 111196, 111196, 111199, 111200, 111200, 111198, 111197, 111195, 111196, 111193, 111194, 111196, 111197, 111197, 111196, 111196, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111193, 111194, 111192, 111192, 111192, 111192, 111193, 111190, 111190, 111190, 111191, 111188, 111186, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111184, 111186, 111185, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111197, 111198, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111197, 111196, 111195, 111194, 111194, 111194, 111194, 111193, 111191, 111190, 111192, 111191, 111190, 111188, 111187, 111188, 111190, 111189, 111188, 111186, 111186, 111185, 111185, 111184, 111184, 111185, 111183, 111180, 111179, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111179, 111179, 111180, 111181, 111182, 111179, 111181, 111181, 111185, 111184, 111184, 111185, 111185, 111186, 111186, 111185, 111187, 111187, 111188, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111187, 111186, 111188, 111186, 111185, 111191, 111188, 111187, 111188, 111188, 111188, 111190, 111190, 111192, 111191, 111192, 111190, 111189, 111188, 111189, 111191, 111191, 111192, 111192, 111191, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111194, 111193, 111191, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111197, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111191, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111194, 111195, 111196, 111196, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111196, 111198, 111198, 111197, 111198, 111200, 111198, 111197, 111201, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111196, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111191, 111191, 111194, 111194, 111195, 111193, 111192, 111192, 111190, 111191, 111192, 111190, 111190, 111190, 111192, 111192, 111192, 111192, 111193, 111193, 111193, 111193, 111192, 111193, 111192, 111190, 111197, 111197, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111201, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111198, 111198, 111197, 111197, 111197, 111199, 111199, 111197, 111197, 111198, 111195, 111195, 111193, 111194, 111194, 111193, 111192, 111191, 111191, 111189, 111192, 111192, 111192, 111192, 111192, 111193, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111193, 111193, 111192]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64059, 64060, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64060, 64059, 64059, 64061, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64061, 64061, 64065, 64065, 64063, 64064, 64065, 64063, 64063, 64064, 64064, 64063, 64065, 64066, 64065, 64065, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64065, 64066, 64065, 64065, 64065, 64067, 64067, 64068, 64069, 64069, 64069, 64069, 64070, 64069, 64069, 64069, 64069, 64070, 64070, 64069, 64069, 64068, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64069, 64069, 64069, 64069, 64069, 64068, 64066, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64068, 64066, 64067, 64067, 64067, 64065, 64066, 64065, 64066, 64066, 64067, 64067, 64067, 64066, 64066, 64066, 64064, 64063, 64064, 64065, 64064, 64063, 64059, 64060, 64061, 64064, 64063, 64061, 64058, 64060, 64058, 64059, 64058, 64058, 64059, 64059, 64060, 64059, 64059, 64058, 64058, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64057, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64067, 64065, 64067, 64065, 64064, 64064, 64064, 64064, 64064, 64063, 64065, 64066, 64065, 64068, 64068, 64067, 64067, 64066, 64067, 64067, 64068, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64068, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64060, 64060, 64060, 64060, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64058, 64055, 64055, 64056, 64055, 64054, 64054, 64053, 64055, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64060, 64058, 64059, 64058, 64059, 64059, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64056, 64057, 64056, 64056, 64057, 64059, 64058, 64057, 64058, 64057, 64058, 64058, 64057, 64058, 64059, 64056, 64058, 64058, 64058, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64063, 64063, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64060, 64060, 64060, 64062, 64061, 64060, 64061, 64060, 64057, 64058, 64057, 64057, 64057, 64056, 64055, 64064, 64063, 64060, 64061, 64062, 64063, 64062, 64060, 64059, 64059, 64060, 64061, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64063, 64065, 64064, 64064, 64063, 64062, 64063, 64062, 64063, 64059, 64058, 64060, 64058, 64058, 64059, 64059, 64057, 64055, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64061, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64062, 64063, 64064, 64062, 64061, 64063, 64062, 64061, 64062, 64061, 64059, 64060, 64061, 64062, 64059, 64058, 64059, 64058, 64059, 64057, 64058, 64058, 64060, 64058, 64057, 64056, 64055, 64055, 64054, 64055, 64054, 64055, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64052, 64053, 64053, 64053, 64054, 64053, 64053, 64054, 64054, 64054, 64052, 64050, 64050, 64049, 64048, 64048, 64048, 64049, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64058, 64057, 64057, 64058, 64059, 64057, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64051, 64051, 64051, 64050, 64048, 64049, 64049, 64049, 64049, 64047, 64046, 64047, 64047, 64043, 64043, 64043, 64044, 64044, 64044, 64044, 64043, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64041, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64042, 64041, 64042, 64042, 64044, 64043, 64043, 64045, 64045, 64046, 64047, 64046, 64048, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64047, 64047, 64051, 64048, 64050, 64049, 64050, 64050, 64051, 64051, 64052, 64051, 64051, 64047, 64048, 64048, 64050, 64050, 64049, 64052, 64052, 64049, 64051, 64051, 64052, 64052, 64049, 64049, 64049, 64049, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64053, 64054, 64054, 64054, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64057, 64056, 64056, 64056, 64056, 64055, 64054, 64054, 64055, 64057, 64059, 64059, 64059, 64058, 64057, 64057, 64056, 64056, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64060, 64061, 64061, 64061, 64058, 64057, 64058, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64064, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64062, 64063, 64061, 64060, 64060, 64060, 64059, 64059, 64057, 64055, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64061, 64060, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64064, 64063, 64064, 64062, 64062, 64062, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64061, 64060, 64061, 64061, 64060, 64059, 64058, 64057, 64056, 64055, 64055, 64053, 64053, 64052, 64054, 64053, 64054, 64053, 64052, 64050, 64049, 64049, 64050, 64048, 64049, 64048, 64049, 64050, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64049, 64048, 64048, 64059, 64060, 64059, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64064, 64063, 64065, 64063, 64063, 64063, 64063, 64062, 64061, 64061, 64062, 64060, 64060, 64063, 64065, 64064, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64058, 64059, 64060, 64059, 64059, 64058, 64059, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64062, 64062, 64061, 64061, 64061, 64060, 64060, 64059, 64060, 64061, 64059, 64058, 64058, 64058, 64057, 64056, 64055, 64053, 64053, 64055, 64053, 64050, 64050, 64050, 64052, 64053, 64051, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64054, 64056, 64055, 64055]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1042808716193668, 37.06017594715826, 36.758467524557425]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.88011504139996],
            'fit_pix': [False, [240.7427698068262, 210.14394033791376]]}


        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_pcwdT.exp_psf_stats
        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.0549619235098362],
            'min_val_pos': [True, [250, 204, 0, 946]],
            'im_rms': [False, 0.012125671956334268],
            'im_sum': [False, 75.75408568766147],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.8864100738559655, 1.0840855127823963, 0.8443649865524704]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462002989996],
            'fit_pix_0': [False, [239.99836406712743, 209.9973158998747]],
            'fit_1': [False, [0.8870507306043662, 1.0832634520434032, 0.844194602156315]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.88011504139996],
            'fit_pix_1': [False, [239.99851823692387, 209.99724648184971]],
            'fit_2': [False, [0.8876820160097865, 1.083932340681245, 0.8429234090714239]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99561005289996],
            'fit_pix_2': [False, [239.99845067901538, 209.9971661763072]]}


        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        # test_mosaic_cube_eph_pcwdT.exp_resid_stats
        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415430.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.06294216960668564],
            'max_val_pos': [True, [320, 214, 0, 617]],
            'min_val': [False, -0.06556452810764313],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.011408315963598554],
            'im_sum': [False, 9275.945943821356],
            'regn_sum': [False, -16.121912690320187],
            'npts_real': [True, 191116800]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_eph_pcwdT)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_eph_pcwdT.exp_model_stats
        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.039518825709819794],
            'max_val_pos': [True, [149, 256, 0, 717]],
            'min_val': [False, -0.0011229703668504953],
            'min_val_pos': [True, [322, 245, 0, 589]],
            'im_rms': [False, 2.176988118839087e-05],
            'im_sum': [False, 8.197739126393571],
            'regn_sum': [False, 1.1301045540021732],
            'mask_non0': [True, 0],
            'npts_real': [True, 191116800]}


        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, 
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube_eph_briggsbwtaper)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_eph_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 45346.8984375],
            'max_val_pos': [True, [0, 0, 0, 283]],
            'min_val': [False, 45300.84765625],
            'min_val_pos': [True, [0, 0, 0, 856]],
            'im_rms': [False, 45327.483918476515],
            'npts_real': [True, 948]}


        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_cube_eph_pcwdT)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        # test_mosaic_cube_eph_pcwdT.exp_wt_stats
        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.3204350769519806],
            'max_val_pos': [True, [240, 210, 0, 595]],
            'min_val': [False, 7.413938146783039e-05],
            'im_rms': [False, 0.1206465686186194],
            'im_sum': [False, 13934226.301663086],
            'npts_0.2': [False, [111207, 111211, 111208, 111208, 111206, 111204, 111203, 111206, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111207, 111206, 111206, 111207, 111206, 111205, 111205, 111206, 111205, 111205, 111207, 111206, 111206, 111206, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111204, 111206, 111204, 111206, 111206, 111206, 111206, 111208, 111209, 111207, 111207, 111208, 111207, 111207, 111207, 111207, 111208, 111207, 111207, 111208, 111209, 111209, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111206, 111206, 111205, 111207, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111209, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111208, 111207, 111205, 111205, 111206, 111206, 111205, 111204, 111206, 111206, 111205, 111205, 111206, 111206, 111205, 111206, 111204, 111206, 111206, 111204, 111203, 111204, 111204, 111203, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111203, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111197, 111198, 111199, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111206, 111205, 111206, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111206, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111206, 111206, 111203, 111203, 111204, 111202, 111202, 111203, 111202, 111202, 111200, 111203, 111201, 111200, 111200, 111201, 111201, 111201, 111196, 111198, 111196, 111196, 111196, 111197, 111196, 111196, 111196, 111196, 111197, 111197, 111198, 111197, 111199, 111199, 111199, 111198, 111200, 111202, 111200, 111200, 111200, 111201, 111201, 111200, 111199, 111200, 111200, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111200, 111198, 111198, 111198, 111201, 111199, 111199, 111199, 111199, 111199, 111201, 111201, 111200, 111201, 111199, 111199, 111201, 111200, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111202, 111201, 111204, 111204, 111202, 111199, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111204, 111203, 111201, 111199, 111201, 111202, 111203, 111202, 111203, 111205, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111203, 111202, 111201, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111202, 111201, 111202, 111198, 111198, 111201, 111199, 111200, 111199, 111200, 111197, 111193, 111194, 111194, 111194, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111199, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111203, 111202, 111203, 111203, 111203, 111204, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111195, 111199, 111199, 111197, 111197, 111196, 111196, 111196, 111199, 111200, 111200, 111198, 111197, 111195, 111196, 111193, 111194, 111196, 111197, 111197, 111196, 111196, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111193, 111194, 111192, 111192, 111192, 111192, 111193, 111190, 111190, 111190, 111191, 111188, 111186, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111184, 111186, 111185, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111197, 111198, 111198, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111197, 111196, 111195, 111194, 111194, 111194, 111194, 111193, 111191, 111190, 111192, 111191, 111190, 111188, 111187, 111188, 111190, 111189, 111188, 111186, 111186, 111185, 111185, 111184, 111184, 111185, 111183, 111180, 111179, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111179, 111179, 111180, 111181, 111182, 111179, 111181, 111181, 111185, 111184, 111184, 111185, 111185, 111186, 111186, 111185, 111187, 111187, 111188, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111187, 111186, 111188, 111186, 111185, 111191, 111188, 111187, 111188, 111188, 111188, 111190, 111190, 111192, 111191, 111192, 111190, 111189, 111188, 111189, 111191, 111191, 111192, 111192, 111191, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111194, 111193, 111191, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111197, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111191, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111194, 111195, 111196, 111196, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111196, 111198, 111198, 111197, 111198, 111200, 111198, 111197, 111201, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111196, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111198, 111198, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111191, 111191, 111194, 111194, 111195, 111193, 111192, 111192, 111190, 111191, 111192, 111190, 111190, 111190, 111192, 111192, 111192, 111192, 111193, 111193, 111193, 111193, 111192, 111193, 111192, 111190, 111197, 111197, 111198, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111201, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111198, 111198, 111197, 111197, 111197, 111199, 111199, 111197, 111197, 111198, 111195, 111195, 111193, 111194, 111194, 111193, 111192, 111191, 111191, 111189, 111192, 111192, 111192, 111192, 111192, 111193, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111193, 111193, 111192]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64059, 64060, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64060, 64059, 64059, 64061, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64061, 64061, 64065, 64065, 64063, 64064, 64065, 64063, 64063, 64064, 64064, 64063, 64065, 64066, 64065, 64065, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64065, 64066, 64065, 64065, 64065, 64067, 64067, 64068, 64069, 64069, 64069, 64069, 64070, 64069, 64069, 64069, 64069, 64070, 64070, 64069, 64069, 64068, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64069, 64069, 64069, 64069, 64069, 64068, 64066, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64068, 64066, 64067, 64067, 64067, 64065, 64066, 64065, 64066, 64066, 64067, 64067, 64067, 64066, 64066, 64066, 64064, 64063, 64064, 64065, 64064, 64063, 64059, 64060, 64061, 64064, 64063, 64061, 64058, 64060, 64058, 64059, 64058, 64058, 64059, 64059, 64060, 64059, 64059, 64058, 64058, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64057, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64067, 64065, 64067, 64065, 64064, 64064, 64064, 64064, 64064, 64063, 64065, 64066, 64065, 64068, 64068, 64067, 64067, 64066, 64067, 64067, 64068, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64068, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64060, 64060, 64060, 64060, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64058, 64055, 64055, 64056, 64055, 64054, 64054, 64053, 64055, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64060, 64058, 64059, 64058, 64059, 64059, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64056, 64057, 64056, 64056, 64057, 64059, 64058, 64057, 64058, 64057, 64058, 64058, 64057, 64058, 64059, 64056, 64058, 64058, 64058, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64063, 64063, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64060, 64060, 64060, 64062, 64061, 64060, 64061, 64060, 64057, 64058, 64057, 64057, 64057, 64056, 64055, 64064, 64063, 64060, 64061, 64062, 64063, 64062, 64060, 64059, 64059, 64060, 64061, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64063, 64065, 64064, 64064, 64063, 64062, 64063, 64062, 64063, 64059, 64058, 64060, 64058, 64058, 64059, 64059, 64057, 64055, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64061, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64062, 64063, 64064, 64062, 64061, 64063, 64062, 64061, 64062, 64061, 64059, 64060, 64061, 64062, 64059, 64058, 64059, 64058, 64059, 64057, 64058, 64058, 64060, 64058, 64057, 64056, 64055, 64055, 64054, 64055, 64054, 64055, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64052, 64053, 64053, 64053, 64054, 64053, 64053, 64054, 64054, 64054, 64052, 64050, 64050, 64049, 64048, 64048, 64048, 64049, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64058, 64057, 64057, 64058, 64059, 64057, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64051, 64051, 64051, 64050, 64048, 64049, 64049, 64049, 64049, 64047, 64046, 64047, 64047, 64043, 64043, 64043, 64044, 64044, 64044, 64044, 64043, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64041, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64042, 64041, 64042, 64042, 64044, 64043, 64043, 64045, 64045, 64046, 64047, 64046, 64048, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64047, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64047, 64047, 64051, 64048, 64050, 64049, 64050, 64050, 64051, 64051, 64052, 64051, 64051, 64047, 64048, 64048, 64050, 64050, 64049, 64052, 64052, 64049, 64051, 64051, 64052, 64052, 64049, 64049, 64049, 64049, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64053, 64054, 64054, 64054, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64057, 64056, 64056, 64056, 64056, 64055, 64054, 64054, 64055, 64057, 64059, 64059, 64059, 64058, 64057, 64057, 64056, 64056, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64060, 64061, 64061, 64061, 64058, 64057, 64058, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64064, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64062, 64063, 64061, 64060, 64060, 64060, 64059, 64059, 64057, 64055, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64061, 64060, 64063, 64063, 64063, 64063, 64064, 64063, 64063, 64063, 64064, 64063, 64064, 64062, 64062, 64062, 64063, 64062, 64062, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64061, 64060, 64061, 64061, 64060, 64059, 64058, 64057, 64056, 64055, 64055, 64053, 64053, 64052, 64054, 64053, 64054, 64053, 64052, 64050, 64049, 64049, 64050, 64048, 64049, 64048, 64049, 64050, 64051, 64050, 64049, 64049, 64049, 64050, 64049, 64049, 64048, 64048, 64059, 64060, 64059, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64064, 64063, 64065, 64063, 64063, 64063, 64063, 64062, 64061, 64061, 64062, 64060, 64060, 64063, 64065, 64064, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64058, 64059, 64060, 64059, 64059, 64058, 64059, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64062, 64062, 64061, 64061, 64061, 64060, 64060, 64059, 64060, 64061, 64059, 64058, 64058, 64058, 64057, 64056, 64055, 64053, 64053, 64055, 64053, 64050, 64050, 64050, 64052, 64053, 64051, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64054, 64056, 64055, 64055]],
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
            stokes='I', specmode='cubesource', nchan=948, start='261.7643758544GHz', width='0.2441755MHz', \
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
            stokes='I', specmode='cubesource', nchan=948, start='261.7643758544GHz', width='0.2441755MHz', \
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
        exp_im_stats = {'com_bmaj': [False, 0.9369754764795376],
            'com_bmin': [False, 0.714328922231186],
            'com_pa': [False, -87.11352315281212],
            'npts': [True, 191116800],
            'npts_unmasked': [False, 105415354.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.07545597106218338],
            'max_val_pos': [True, [319, 243, 0, 511]],
            'min_val': [False, -0.0654265433549881],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.01138825435221424],
            'rms_per_chan': [False, [0.01391551567345628, 0.011198813123537344, 0.010345340419894959, 0.010811274412962498, 0.010980385487975988, 0.012191143283614564, 0.012655534175377578, 0.014079166218826298, 0.011693246507781589, 0.011473161953430535, 0.01036509324966158, 0.010185600376660019, 0.012665564711609607, 0.010700887971683245, 0.010204218275929773, 0.009966093553230537, 0.010963421132572218, 0.012603110329282182, 0.008792978601938796, 0.011317940881447585, 0.011797013644227599, 0.011942828512409687, 0.011521270528107433, 0.01247837136416052, 0.011404686740893335, 0.010493928899296148, 0.010533345727638473, 0.010748676360497075, 0.013014328613251656, 0.012113374360662368, 0.010949392947331485, 0.0104040724325028, 0.010570279485877935, 0.009853135727721589, 0.011297806905893708, 0.011209623201283268, 0.0106340973552856, 0.011160040827923315, 0.010418842531732761, 0.010695826361701173, 0.010991351793005264, 0.011658097363124797, 0.011025582495098267, 0.012592575648849613, 0.011468895561306314, 0.010916062854478307, 0.010244011995244775, 0.009840442611445854, 0.010074028270638033, 0.010442778110347897, 0.011273588111584297, 0.011736170779119324, 0.010844709387112735, 0.010936670112463518, 0.010720025449490744, 0.011017391510606475, 0.01063744025823964, 0.014325808895565496, 0.0115557081884916, 0.00923052251363914, 0.009561297096677555, 0.01174533333395733, 0.010136810058581108, 0.010565585125265413, 0.010246552076802926, 0.01056456720492672, 0.010977619843268623, 0.011296872184220713, 0.010549875904601417, 0.009447489511758254, 0.012454700907368103, 0.011625889249378475, 0.010657897093555515, 0.009149386889594392, 0.010530641723948242, 0.012180883491652838, 0.012866352369947681, 0.012653690321379574, 0.012075212105291451, 0.01145391073954902, 0.010329345287953384, 0.010210823768431073, 0.010712144642626212, 0.011453256770831805, 0.010870741976014217, 0.01072160369091619, 0.01044009423911047, 0.010666232894270758, 0.009588296187507343, 0.010827533142129158, 0.011019715788050544, 0.011816983805996975, 0.010932443307819282, 0.009815306845600016, 0.010194412850479275, 0.011503936803220299, 0.011470414777124219, 0.011235259536549037, 0.010993846172661236, 0.00903766289885872, 0.011166712018514637, 0.011978779079066462, 0.012787919649476436, 0.011722765399185715, 0.012038142657925457, 0.012238668733087056, 0.013354472436258896, 0.01136602311610165, 0.012444105036993268, 0.011591124374724322, 0.010549276880037924, 0.010481459506565658, 0.011463814602880623, 0.01061697870999084, 0.010218277881268894, 0.011585800843130207, 0.010579855240934197, 0.010043023725294404, 0.011742356838586324, 0.011120549918916912, 0.009768703574725066, 0.011408413544833518, 0.009595849589896756, 0.011414755024059987, 0.010904002508659708, 0.011381100891792941, 0.010138297997786417, 0.012524523642418283, 0.010368857237119633, 0.009189210260360159, 0.010717792727978586, 0.010385047802135778, 0.009668636919135873, 0.011288628646689006, 0.01226453677286331, 0.010810164828898535, 0.010900410347242309, 0.01073834085489203, 0.010054557996776342, 0.01058190139223016, 0.01112769307494144, 0.010355479854678822, 0.009836085362232252, 0.010176917852749, 0.0089835084142114, 0.01140878518608462, 0.011132856825588405, 0.009518558455639786, 0.011420761430388523, 0.009742653641841278, 0.010541266474157504, 0.011071649044048996, 0.010572080379685642, 0.010419287576888581, 0.011533118229040589, 0.01079756738002291, 0.011208088280964005, 0.010973999363405545, 0.010342535577589197, 0.01360858079298801, 0.013346554388556493, 0.01081165486577898, 0.009671134298033126, 0.00974645329159148, 0.010395247664496806, 0.009902441097076624, 0.010071512650368323, 0.009123944267933446, 0.009670367381145469, 0.010561708075660408, 0.009374145734522306, 0.010480613589791602, 0.010188881526182294, 0.010701719237319189, 0.011973884113993664, 0.011508903889746258, 0.009606412277797797, 0.011879767711977917, 0.012051196773676552, 0.011811551483983138, 0.010048625839022163, 0.011492812761157266, 0.010929950311566461, 0.01223238383113695, 0.01052752405439503, 0.010230198195914286, 0.01139858392141661, 0.009929041715717861, 0.012128463571431863, 0.0102966296601799, 0.010733421630739069, 0.010298427016025418, 0.01187617441191375, 0.010012077602361613, 0.010606446457556056, 0.011796736959529392, 0.010165267691489592, 0.011672502651117342, 0.010193309032991338, 0.011390952452404632, 0.011970661120405518, 0.013890279243384427, 0.01245774213057081, 0.011009134249838832, 0.010203358081979731, 0.00937986973959424, 0.010029155665712653, 0.01029865205962035, 0.011316283888592031, 0.01093328067262625, 0.011546369159699154, 0.011900073114081743, 0.011136586569292822, 0.010910845781053766, 0.009887444991602188, 0.01210332730597261, 0.01169126529722959, 0.010231104150216084, 0.010674369291115909, 0.012508523181995508, 0.012608837843702496, 0.011041149758156338, 0.011192281551465883, 0.012825393640233185, 0.011836889418171414, 0.009997617957231997, 0.011545436933067615, 0.011156641657883177, 0.011559415957711613, 0.010430718794276825, 0.009994723339651105, 0.012719483986911648, 0.011570039242933336, 0.011115070523757615, 0.009256703143477844, 0.011979715844125823, 0.01157055108823778, 0.01115699215897198, 0.01077075613713458, 0.013396495730722524, 0.011727621917112801, 0.009219980862272792, 0.00953358502102843, 0.010912021877420653, 0.011574615990497597, 0.010152294064037795, 0.01051687628723948, 0.01056251557771618, 0.011302592044929533, 0.01136320638487575, 0.01269415147830265, 0.01040936046532526, 0.00936827588400146, 0.01165178441401517, 0.011561000909089937, 0.012377491714670702, 0.011898334521497491, 0.0102495022172728, 0.01073966773082043, 0.010422736316839277, 0.011098521560582735, 0.011086588672418058, 0.012496652625511653, 0.010784776326271088, 0.011004579018834703, 0.011628658667145072, 0.0120665759226352, 0.010096402548412669, 0.011592984085871566, 0.012300668050443006, 0.010880753232816672, 0.010407554889909798, 0.010682966499209853, 0.011391599367702702, 0.010459876135538065, 0.010679522774863785, 0.012010982549606468, 0.012591121909487682, 0.012373586030356034, 0.011096697003823276, 0.010098737365759863, 0.01079718160149053, 0.011059175772192799, 0.008699436224291472, 0.011117258550168181, 0.011521298995493125, 0.0116907022709856, 0.01143472344879642, 0.011644249457362595, 0.012615807414612969, 0.010091562674338979, 0.01016033781618809, 0.010965081423838976, 0.010957455775951107, 0.01165532279623697, 0.00964942467902527, 0.010641495187787548, 0.010368737724126682, 0.011305078692750416, 0.01012587036515205, 0.012687409477921205, 0.011231199903396915, 0.011806009103003515, 0.01506330514919511, 0.015075603541390946, 0.014984096748374406, 0.010093661661164186, 0.011503207407057592, 0.010580984678778665, 0.01302805928019781, 0.010568295375048668, 0.010813162725367507, 0.011914499539068141, 0.012303520180347303, 0.01182346881689707, 0.01108994484072009, 0.011258519264891173, 0.01376157080278537, 0.011491594691571207, 0.01136533616605801, 0.011265139203673084, 0.015059264139043956, 0.011191801261909358, 0.011101309583405152, 0.010766794482062077, 0.01272392361288182, 0.013594050568104631, 0.0116258354793538, 0.010207657821216775, 0.011221513210258821, 0.010699688518146865, 0.01257831136093649, 0.01109329974751451, 0.011725835674511464, 0.011556158403140908, 0.01031770183476827, 0.013423810023143468, 0.012169170240339172, 0.010884181056196568, 0.013219747679734682, 0.010424564368164566, 0.01288569704155284, 0.00947356549182442, 0.010661602021293416, 0.011516515575887253, 0.01306304685508459, 0.010298022620306848, 0.011692784249875054, 0.010008629208782203, 0.010577465457789566, 0.011660345993726041, 0.012412983436878435, 0.012081989543155447, 0.011248495542875618, 0.011203927610165042, 0.010451326863818826, 0.012309473277926504, 0.011686218679865537, 0.012205192543136487, 0.011690905887865975, 0.01090712213526366, 0.011352194732435289, 0.013105331937951948, 0.013449661754540916, 0.0104145758665503, 0.01105667310343271, 0.010449982772855117, 0.011561698977045239, 0.012604025769663142, 0.010782981652464312, 0.011641618848162916, 0.012954189541278683, 0.011967775316031427, 0.010507995510634813, 0.012014661812115036, 0.01128242812999333, 0.010790432467561338, 0.011494794254532818, 0.011254069964641665, 0.009614709928895402, 0.01047387885692144, 0.010434069944095426, 0.012208279425905553, 0.01049515603390529, 0.011760783954543573, 0.01058840194134269, 0.012048299295972292, 0.012935624163251395, 0.010493419781751674, 0.012958610519237282, 0.01200499263318809, 0.011072342770609367, 0.010560716266239585, 0.011321540764256108, 0.011163118754599476, 0.009882703631943743, 0.009940915626621335, 0.011117996900318469, 0.011504650877207725, 0.014774711597258015, 0.011374638966577392, 0.01127153301735087, 0.01123578325598108, 0.012320248016013547, 0.011040033388996088, 0.00976207700781444, 0.010952589142997234, 0.013376855026179262, 0.014771197069397246, 0.010821479927248953, 0.011536002156779161, 0.011420030412244745, 0.012530904455104111, 0.011673877887152935, 0.011182299948114325, 0.013046223919109444, 0.012497582674964534, 0.011107827833385065, 0.011518419907441606, 0.010913673318033898, 0.01060017708795507, 0.00899894907135792, 0.01038245697515491, 0.010273633958549952, 0.011313525903307774, 0.011323165045246672, 0.01139332979415624, 0.012374088003489615, 0.010846875329190117, 0.011009469191714416, 0.011364296685645613, 0.011137214994526303, 0.011459177716679913, 0.011378136379264674, 0.010506628387558576, 0.01323121123896279, 0.011408247969116499, 0.013315475615833672, 0.01188893441512642, 0.011717811820115756, 0.011536696713093687, 0.01381308497552019, 0.012909669223975002, 0.012052432776570917, 0.010663297919072178, 0.011236748557025281, 0.011417204578490225, 0.011095209740528637, 0.010839760972536343, 0.01149151405819477, 0.013292450666728076, 0.013290942443756752, 0.01157408796605491, 0.012578005105955171, 0.011693245886345044, 0.010496704688724432, 0.010740626671809101, 0.010010859488681986, 0.012413441211024199, 0.0111771015375279, 0.012717297337893138, 0.011239869815280454, 0.011295972165029521, 0.011451193895247251, 0.010010935624683026, 0.011685132675417445, 0.01066902459651249, 0.010734910511270212, 0.011522490602286457, 0.0111625617236132, 0.012453652154402507, 0.01077834604770361, 0.01133891022845769, 0.01080167632551295, 0.011081788903629623, 0.01151597543815762, 0.011359160652902326, 0.012790343301677243, 0.010942661113587283, 0.010993110777618225, 0.011503303126971972, 0.011211381123856666, 0.01136618108236803, 0.01526613213144868, 0.011301931623827068, 0.012021306396428837, 0.012050718223203413, 0.011975360347810615, 0.011599033183117623, 0.011773678079377665, 0.011962216441033048, 0.010862671183139973, 0.0116288915814634, 0.012338926394263355, 0.010258065701063506, 0.01169652677840939, 0.011650658142161052, 0.011464672233430397, 0.013263557550491715, 0.012140503025387309, 0.01194531784247912, 0.013422603089499344, 0.010452626104052145, 0.012466031418590038, 0.011792624772612847, 0.011942971840778874, 0.010774727251425087, 0.012631943286551622, 0.011771710373427821, 0.0114037432592323, 0.012538754797458255, 0.012966952553430005, 0.012342390461358056, 0.010612501763890395, 0.012582420042745459, 0.01294395834067153, 0.012398621986183403, 0.013712645273892156, 0.013034791169337609, 0.01200322609060809, 0.010161505762754692, 0.013053053524129724, 0.012201202735098574, 0.012127196336283603, 0.012404880207630031, 0.009985547344910972, 0.011299663363624748, 0.013111067507206713, 0.011188559302689692, 0.011824333574171177, 0.010839952738754715, 0.011193187589903264, 0.010876358797037692, 0.011219974059371918, 0.01201633585948114, 0.012044551730118282, 0.01296471624433944, 0.010788423266819827, 0.011009711238469836, 0.011933472009159209, 0.011877563687803601, 0.011109406021886366, 0.011386069046057871, 0.011831740573136865, 0.01149286871374616, 0.011603462340601943, 0.010788340834009983, 0.012400325913548574, 0.01172042951924379, 0.011508336690804491, 0.011603863620446065, 0.01140454758188262, 0.011738446323379096, 0.013515722571110306, 0.012372465991333601, 0.0121227469614165, 0.011102598343154944, 0.011491254228204672, 0.012241753399869595, 0.010539041691398094, 0.01224561210231814, 0.013550939413897964, 0.012181611935657527, 0.011459087347053958, 0.012369887184443956, 0.010815942800915879, 0.012359833407768046, 0.013842524915639268, 0.012245748804495747, 0.011298177134238358, 0.011194627818937387, 0.011153466016506314, 0.011182445216799963, 0.012024873331244057, 0.01270776619936928, 0.010494129612726808, 0.011599110075928324, 0.011561517779585946, 0.011031652460864823, 0.012348419825214113, 0.012338799395344395, 0.013460152696158127, 0.011563583231934085, 0.012428034338536422, 0.013263765175452425, 0.011160583396022295, 0.010557806307679033, 0.01224672302896639, 0.012121309130817264, 0.013916453902819157, 0.014077298279752048, 0.012711605066929196, 0.012250003914288488, 0.013764666026541743, 0.012551498639013161, 0.012151557472141155, 0.011029666905482212, 0.012436344694198326, 0.01103359090008633, 0.011165879337660611, 0.013703783935326216, 0.013138992153397698, 0.01061770802488673, 0.011379580871463194, 0.011363558769280465, 0.011301655797430445, 0.012213251940275152, 0.010188687808435192, 0.01093307948833843, 0.011300172604409197, 0.013429973366686596, 0.012629334827924249, 0.014592125389570249, 0.011330657609183274, 0.010738440896738487, 0.013106740883559953, 0.010813565199318004, 0.013139196624239438, 0.01289408201147847, 0.012925366980481644, 0.010903130670237137, 0.012349194729709405, 0.012112273832888274, 0.012606279152308789, 0.01268010835742184, 0.011152480237942143, 0.0131034147826204, 0.012135949055839654, 0.01054169972935904, 0.012954294049810185, 0.013498063316998349, 0.010689926036177697, 0.012645585515160799, 0.010939852305324379, 0.010198068807222845, 0.013250607762112484, 0.011503714882512108, 0.010729180287361308, 0.0117603175714084, 0.012003440050524547, 0.011666573293905795, 0.012644084576131592, 0.011481595584960601, 0.013115499459816035, 0.012056553370898936, 0.013342524059805629, 0.011089557500685748, 0.010373673319258344, 0.012682253302271254, 0.011878806051804332, 0.010689626984112435, 0.012863815083887164, 0.01085789514807431, 0.012455020887258656, 0.011122446017804741, 0.01137685400909025, 0.012221378234692873, 0.011550078616777907, 0.01227001769463364, 0.011390940813890495, 0.0119687132200359, 0.012306802523615809, 0.013961020939233947, 0.012800009293421878, 0.013208848259047504, 0.011190368947139139, 0.01174188589309877, 0.011222246238821105, 0.012303348392178856, 0.011209902267838696, 0.009619180150990737, 0.010316924946913893, 0.011189554187211436, 0.011635960931678534, 0.014558189535781018, 0.010972710109294065, 0.011745480701183461, 0.011119071848271713, 0.010919774271012239, 0.014769272715155945, 0.011287044162453371, 0.0109246067254714, 0.011820284065314601, 0.010186394935158796, 0.012409454592021741, 0.011663031392865755, 0.011342310110051003, 0.010527865700960173, 0.012224475693825062, 0.01104268732124472, 0.009826944684556263, 0.012545432034855837, 0.010480227461350721, 0.010496787531049849, 0.01125771922624902, 0.012312540487719218, 0.01207704300197997, 0.012234200225847987, 0.011685728360191504, 0.01157707846019137, 0.01196572637307069, 0.011693324337565878, 0.011986360330347943, 0.012992981866439324, 0.011680706871814833, 0.0114791720673399, 0.012110186278473838, 0.010618820888381186, 0.013792118712727987, 0.012443766995919946, 0.012297074409095122, 0.01301595775530087, 0.010561100212335944, 0.010754847935992918, 0.011585717748904531, 0.01041456820730729, 0.0112199249155143, 0.009926342690918678, 0.011219993363492563, 0.011806962512050044, 0.010287643219072226, 0.011895677416971553, 0.012462596288443912, 0.010659439525448759, 0.010453644357873769, 0.013145997469789727, 0.010972772727225468, 0.01288923325566898, 0.010742121125745957, 0.010889638693321782, 0.012120258362527926, 0.011966472160895915, 0.009949471030270716, 0.009839749262065328, 0.009973078454802603, 0.010933075471764827, 0.01060638993098246, 0.012286990439496215, 0.010438470622262108, 0.010356165575680594, 0.013147829987749163, 0.011503997929352479, 0.010433809094295407, 0.009873918594700304, 0.010621956831719652, 0.0119549965604376, 0.011589897569645406, 0.010274629111251104, 0.011691095946987104, 0.013977870188447398, 0.010060873324879281, 0.010190241744687996, 0.009814124490438883, 0.010901623127690099, 0.012685857237233703, 0.011375842848159701, 0.01069010023721738, 0.011367480110733238, 0.011531773290414814, 0.011363141945805766, 0.010396688594581247, 0.014102047656936085, 0.011268816423457692, 0.012612800193041639, 0.014773048045272718, 0.01076767128510351, 0.010501749810623658, 0.012447988623926176, 0.012214461134542014, 0.012330200793774638, 0.011934973211024087, 0.010871025246303349, 0.010330734021745763, 0.012066539477031957, 0.011258502590474697, 0.010481886035595024, 0.011641613815714226, 0.010697465976347997, 0.01063974836908695, 0.010272304039852421, 0.012982026257443772, 0.011204399510922287, 0.012167644387189354, 0.010943922826048561, 0.012059322542934972, 0.01212136785938398, 0.01113752353832333, 0.01103761873105044, 0.012172340521742936, 0.012192280403185186, 0.012638606311108072, 0.009666082928467531, 0.010827198995898672, 0.009989254933653461, 0.010766343151204736, 0.010523079936086662, 0.009480797063413782, 0.012974173906662887, 0.010451118264427971, 0.010774942144578742, 0.010068511119076301, 0.011096479735761862, 0.010473014479068682, 0.010952564841578671, 0.010404561954280443, 0.011906488710142158, 0.010044360213585223, 0.010675026057278372, 0.01238813057303107, 0.01102580504548289, 0.011130866137614319, 0.009247742813543161, 0.011360260958360121, 0.009722381180439672, 0.011067493252939969, 0.011083008453611695, 0.011347008857267884, 0.011330175731110484, 0.012735751826678406, 0.010875781053330653, 0.009592931784560565, 0.010195486601363851, 0.011966287452965083, 0.010062644233046565, 0.010226784832337459, 0.009625918296838681, 0.010347269612854084, 0.012349513190613302, 0.013722032275909152, 0.011046561418105911, 0.011323547931902904, 0.012533637327172071, 0.011562137265752967, 0.010841189478950087, 0.01128045928398688, 0.011423024199729596, 0.01157569505918874, 0.011025249697811813, 0.011245214145903617, 0.010637604657311882, 0.010362950547746, 0.012097746021352416, 0.01086892499251379, 0.012332641800251752, 0.010437202052044289, 0.010518430848280405, 0.011265188002607257, 0.009446649170151047, 0.009313443411281427, 0.009950020938308227, 0.010352970403371981, 0.01120856919648851, 0.00999711443189775, 0.010970975051248241, 0.011178070257493464, 0.01109598885398961, 0.010020899626209173, 0.01068431468070205, 0.01032706506180335, 0.009642814392534374, 0.010926095321583337, 0.012362877682361708, 0.010494021286914033, 0.009819067073263311, 0.00931149262950656, 0.009743758926703526, 0.012064019100826428, 0.01142151547989276, 0.01218306069157812, 0.012267823706118366, 0.010923883485025659, 0.011066397991132858, 0.01201974828647091, 0.011075253829816417, 0.009709387561967124, 0.009977071409927643, 0.00958955819211091, 0.0104295434121177, 0.010148550862387666, 0.01156742075660445, 0.01167526642424578, 0.011041524105200186, 0.011429319422665804, 0.010875854302731846, 0.012968677540407597, 0.011974985417507491, 0.011102470907026483, 0.009173031798793892, 0.010262012294932885, 0.009774694126784163, 0.009948520531037628, 0.010860247819736294, 0.010732870358521532, 0.01101143231566375, 0.009078736041032105, 0.011654313259812988, 0.012201486714527494, 0.010142054258460756, 0.010440632541945435, 0.01033896005505243, 0.009840414820710316, 0.009307337990383946, 0.010468262860276429, 0.01015023193353219, 0.00972346195351201, 0.010032669260383301, 0.009452312844075035, 0.010903728244818642, 0.01053900088745352, 0.01058487381303189, 0.01088981282562214, 0.013068034166664003, 0.01024335498989072, 0.0099025550259023, 0.0106510439612941, 0.010640736068391972, 0.011024260755162064, 0.010845663505606009, 0.013798681901385973, 0.01116586973647823, 0.00968625348566756, 0.010610949217243839, 0.009476606643439772, 0.011924206044332279, 0.010365470663524905, 0.012581502692677621, 0.010596722985717436, 0.009237430140717045, 0.01117810832618649, 0.009935962754534854, 0.00925250109079027, 0.012648952116913947, 0.012817826755854337, 0.010407021096079846, 0.01002431077596234, 0.009763642939969097, 0.010450134341096633, 0.010207801770397316, 0.009826923540041877, 0.010159473748469525, 0.010669894420236555, 0.011203266954386028, 0.01057176040813721, 0.011955848725536922, 0.012407740563333845, 0.010151934263707767, 0.009964868788700714, 0.011218820859613681]],
            'im_sum': [False, 9556.276200833598],
            'regn_sum': [False, 25.980929710112832],
            'npts_real': [True, 191116800],
            'rms_per_field': [False, [0.012620097693620318, 0.01156422398409204, 0.011867619720565525, 0.01192928179151625, 0.011527934240669084, 0.012438354923861005, 0.011785502633137122]],
            'profile': [False, 0.022832596793556984]}

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
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'mask_pix': [False, 8534],
            'mask_regns': [True, 31],
            'npts_real': [True, 191116800]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse [[239.37091637deg, -16.96407526deg], [28.1142arcsec, 27.0960arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415354.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, 0.20000001788139343],
            'im_rms': [False, 0.6370609562630469],
            'npts_0.2': [False, [111206, 111210, 111208, 111207, 111205, 111204, 111201, 111204, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111208, 111206, 111206, 111207, 111205, 111205, 111205, 111206, 111204, 111205, 111206, 111205, 111205, 111205, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111203, 111206, 111204, 111206, 111206, 111206, 111206, 111207, 111208, 111207, 111207, 111207, 111206, 111207, 111207, 111206, 111207, 111207, 111207, 111207, 111209, 111208, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111207, 111206, 111205, 111207, 111205, 111206, 111206, 111206, 111205, 111205, 111207, 111208, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111205, 111206, 111206, 111206, 111206, 111204, 111205, 111205, 111205, 111205, 111206, 111205, 111205, 111206, 111204, 111206, 111206, 111203, 111203, 111204, 111204, 111204, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111203, 111202, 111202, 111202, 111202, 111203, 111203, 111201, 111202, 111202, 111203, 111203, 111202, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111196, 111198, 111198, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111207, 111205, 111205, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111206, 111204, 111205, 111206, 111206, 111205, 111206, 111206, 111205, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111206, 111203, 111203, 111203, 111201, 111202, 111203, 111201, 111203, 111200, 111202, 111201, 111200, 111200, 111201, 111201, 111201, 111197, 111198, 111197, 111196, 111197, 111196, 111196, 111196, 111196, 111196, 111197, 111197, 111197, 111197, 111198, 111199, 111198, 111197, 111200, 111202, 111200, 111200, 111200, 111200, 111200, 111200, 111199, 111200, 111199, 111197, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111200, 111198, 111198, 111198, 111201, 111199, 111199, 111199, 111198, 111199, 111201, 111200, 111200, 111201, 111198, 111199, 111201, 111199, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111201, 111201, 111203, 111204, 111202, 111198, 111199, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111204, 111203, 111201, 111199, 111201, 111201, 111203, 111202, 111202, 111204, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111202, 111201, 111201, 111202, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111200, 111202, 111199, 111199, 111201, 111199, 111199, 111200, 111200, 111197, 111193, 111194, 111194, 111195, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111198, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111196, 111198, 111199, 111197, 111196, 111197, 111196, 111196, 111199, 111200, 111199, 111197, 111197, 111195, 111195, 111193, 111194, 111196, 111198, 111196, 111196, 111195, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111192, 111193, 111192, 111192, 111192, 111191, 111193, 111190, 111190, 111190, 111191, 111188, 111185, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111185, 111186, 111184, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111197, 111196, 111195, 111194, 111194, 111194, 111193, 111193, 111191, 111190, 111191, 111191, 111191, 111188, 111186, 111187, 111190, 111188, 111188, 111187, 111186, 111185, 111185, 111184, 111184, 111184, 111184, 111181, 111181, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111179, 111179, 111180, 111181, 111181, 111178, 111181, 111181, 111185, 111184, 111184, 111185, 111185, 111186, 111185, 111185, 111188, 111187, 111189, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111187, 111186, 111187, 111186, 111186, 111191, 111186, 111187, 111189, 111188, 111188, 111190, 111190, 111192, 111191, 111192, 111190, 111190, 111188, 111189, 111191, 111191, 111192, 111192, 111190, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111193, 111191, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111196, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111191, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111193, 111195, 111195, 111195, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111197, 111198, 111196, 111198, 111198, 111198, 111198, 111200, 111198, 111197, 111201, 111200, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111197, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111197, 111199, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111192, 111191, 111194, 111194, 111194, 111193, 111191, 111192, 111190, 111190, 111191, 111190, 111190, 111190, 111192, 111192, 111191, 111192, 111192, 111193, 111193, 111192, 111192, 111192, 111191, 111190, 111197, 111198, 111197, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111202, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111196, 111197, 111197, 111199, 111199, 111197, 111197, 111198, 111195, 111195, 111193, 111193, 111194, 111193, 111191, 111191, 111191, 111190, 111192, 111192, 111192, 111192, 111192, 111193, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111194, 111193, 111192]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64058, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64060, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64062, 64061, 64065, 64065, 64063, 64062, 64064, 64063, 64061, 64063, 64064, 64063, 64065, 64065, 64065, 64064, 64064, 64064, 64064, 64065, 64064, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64064, 64066, 64064, 64065, 64066, 64067, 64067, 64068, 64068, 64069, 64070, 64069, 64070, 64069, 64070, 64069, 64069, 64070, 64070, 64069, 64069, 64067, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64067, 64069, 64069, 64069, 64068, 64068, 64065, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64066, 64066, 64068, 64066, 64065, 64065, 64065, 64065, 64066, 64066, 64067, 64067, 64066, 64066, 64066, 64066, 64065, 64063, 64064, 64065, 64064, 64062, 64059, 64060, 64061, 64064, 64063, 64060, 64058, 64060, 64058, 64059, 64058, 64058, 64058, 64059, 64059, 64059, 64059, 64058, 64057, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64066, 64065, 64067, 64065, 64064, 64064, 64064, 64063, 64064, 64063, 64064, 64066, 64064, 64068, 64068, 64066, 64066, 64066, 64067, 64067, 64068, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64067, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64059, 64060, 64060, 64059, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64057, 64055, 64056, 64056, 64055, 64054, 64054, 64053, 64054, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64055, 64056, 64055, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64058, 64058, 64057, 64058, 64059, 64056, 64057, 64058, 64057, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64059, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64063, 64062, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64061, 64060, 64060, 64062, 64061, 64060, 64060, 64060, 64057, 64058, 64057, 64057, 64056, 64055, 64056, 64064, 64063, 64059, 64061, 64062, 64062, 64061, 64060, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64061, 64062, 64062, 64062, 64063, 64064, 64064, 64064, 64063, 64061, 64062, 64062, 64063, 64059, 64058, 64059, 64058, 64058, 64059, 64059, 64057, 64055, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64054, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64060, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64063, 64063, 64064, 64062, 64062, 64062, 64062, 64061, 64061, 64061, 64059, 64060, 64061, 64060, 64059, 64057, 64058, 64058, 64059, 64058, 64058, 64058, 64060, 64058, 64057, 64054, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64054, 64053, 64051, 64050, 64050, 64050, 64048, 64048, 64048, 64050, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64059, 64057, 64057, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64051, 64051, 64051, 64050, 64048, 64049, 64049, 64050, 64049, 64047, 64046, 64046, 64046, 64044, 64044, 64043, 64044, 64044, 64044, 64044, 64042, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64042, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64042, 64040, 64042, 64042, 64044, 64044, 64043, 64044, 64045, 64045, 64046, 64044, 64048, 64048, 64047, 64047, 64048, 64048, 64046, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64047, 64047, 64051, 64049, 64049, 64049, 64051, 64050, 64051, 64050, 64052, 64051, 64051, 64047, 64048, 64047, 64049, 64049, 64049, 64051, 64051, 64050, 64051, 64051, 64051, 64051, 64049, 64049, 64048, 64048, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64055, 64055, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64056, 64056, 64055, 64054, 64055, 64055, 64058, 64058, 64058, 64058, 64058, 64057, 64057, 64056, 64057, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64059, 64060, 64061, 64061, 64058, 64058, 64058, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64063, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64064, 64064, 64062, 64063, 64061, 64060, 64060, 64060, 64059, 64058, 64057, 64054, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64060, 64060, 64063, 64064, 64062, 64063, 64064, 64063, 64063, 64062, 64064, 64063, 64063, 64062, 64062, 64062, 64063, 64063, 64061, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64061, 64061, 64061, 64061, 64060, 64059, 64058, 64056, 64056, 64055, 64055, 64053, 64052, 64052, 64054, 64051, 64052, 64053, 64052, 64050, 64049, 64048, 64048, 64048, 64049, 64048, 64049, 64049, 64050, 64049, 64049, 64049, 64049, 64050, 64049, 64048, 64048, 64048, 64059, 64059, 64058, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64065, 64063, 64062, 64063, 64062, 64062, 64061, 64061, 64062, 64061, 64060, 64064, 64064, 64065, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64059, 64058, 64060, 64059, 64058, 64059, 64059, 64061, 64059, 64060, 64060, 64060, 64061, 64061, 64062, 64062, 64060, 64061, 64061, 64060, 64060, 64059, 64060, 64061, 64058, 64058, 64058, 64058, 64057, 64056, 64055, 64054, 64053, 64055, 64053, 64050, 64050, 64051, 64052, 64053, 64052, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64054, 64056, 64056, 64055]],
            'npts_real': [True, 191116800],
            'fit': [False, [1.1042804178790713, 37.06023555493437, 36.758390484594315]],
            'fit_loc_chan': [True, 474],
            'fit_loc_freq': [1e-10, 261.88011504139996],
            'fit_pix': [False, [240.74284782362395, 210.1445947261468]]}


        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report (test_mosaic_cube_eph)
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[239.36977532deg, -16.96391179deg], [1.0415arcsec, 0.9313arcsec], 90.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [240, 210, 0, 0]],
            'min_val': [False, -0.05504078045487404],
            'min_val_pos': [True, [250, 204, 0, 946]],
            'im_rms': [False, 0.012116211613735845],
            'im_sum': [False, 75.47794335577694],
            'npts_real': [True, 191116800],
            'fit_0': [False, [0.886326655820083, 1.0827417193282847, 0.8432044545236995]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 261.76462002989996],
            'fit_pix_0': [False, [239.99836883975055, 209.99733278296767]],
            'fit_1': [False, [0.886964256657234, 1.0819245159296407, 0.8430338914805372]],
            'fit_loc_chan_1': [True, 474],
            'fit_loc_freq_1': [1e-10, 261.88011504139996],
            'fit_pix_1': [False, [239.9985221259162, 209.99726354985233]],
            'fit_2': [False, [0.8875922689198896, 1.08258955698927, 0.8417657253630042]],
            'fit_loc_chan_2': [True, 947],
            'fit_loc_freq_2': [1e-10, 261.99561005289996],
            'fit_pix_2': [False, [239.99845501846443, 209.99718365233807]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report (test_mosaic_cube_eph_briggsbwtaper)
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]')

        # test_mosaic_cube_eph_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 191116800],
            'npts_unmasked': [False, 105415354.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.06285478919744492],
            'max_val_pos': [True, [320, 214, 0, 617]],
            'min_val': [False, -0.0654265433549881],
            'min_val_pos': [True, [268, 151, 0, 57]],
            'im_rms': [False, 0.011382174483049073],
            'im_sum': [False, 9240.214643091973],
            'regn_sum': [False, -16.069532327876914],
            'npts_real': [True, 191116800]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report (test_mosaic_cube_eph_briggsbwtaper)
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[239.36796846deg, -16.96307757deg], [5.3756arcsec, 3.3987arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_mosaic_cube_eph_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.03946136683225632],
            'max_val_pos': [True, [149, 256, 0, 717]],
            'min_val': [False, -0.001125899376347661],
            'min_val_pos': [True, [322, 245, 0, 589]],
            'im_rms': [False, 2.1682798001723417e-05],
            'im_sum': [False, 8.168397554662079],
            'regn_sum': [False, 1.1245928519638255],
            'mask_non0': [True, 0],
            'npts_real': [True, 191116800]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, 
            '.model', epsilon=self.epsilon)

        # .sumwt report (test_mosaic_cube_eph_briggsbwtaper)
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_mosaic_cube_eph_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 948],
            'npts_unmasked': [True, 948.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 45549.265625],
            'max_val_pos': [True, [0, 0, 0, 283]],
            'min_val': [False, 45503.4375],
            'min_val_pos': [True, [0, 0, 0, 856]],
            'im_rms': [False, 45529.93611237143],
            'npts_real': [True, 948]}


        report8 = th.check_dict_vals(exp_sumwt_stats, sumwt_stats_dict, \
            '.sumwt', epsilon=self.epsilon)

        # .weight report (test_mosaic_cube_eph_briggsbwtaper)
        wt_stats_dict = self.image_stats(img+'.weight', masks=[ \
            pb_stats_dict['pb_mask_0.2'], pb_stats_dict['pb_mask_0.5']])

        # test_mosaic_cube_eph_briggsbwtaper.exp_wt_stats
        exp_wt_stats = {'npts': [True, 191116800],
            'npts_unmasked': [True, 191116800.0],
            'freq_bin': [1e-10, 244175.5],
            'start': [True, 261764400000.0],
            'end': [True, 261995600000.0],
            'start_delta': [False, 261764400000.0],
            'end_delta': [False, 261995600000.0],
            'nchan': [True, 948],
            'max_val': [False, 0.32043445110321045],
            'max_val_pos': [True, [240, 210, 0, 595]],
            'min_val': [False, 7.414929132210091e-05],
            'im_rms': [False, 0.12064659331008178],
            'im_sum': [False, 13934226.308718286],
            'npts_0.2': [False, [111206, 111210, 111208, 111207, 111205, 111204, 111201, 111204, 111203, 111205, 111202, 111202, 111204, 111204, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111205, 111204, 111205, 111205, 111205, 111205, 111205, 111205, 111204, 111204, 111205, 111205, 111205, 111205, 111206, 111206, 111207, 111208, 111206, 111206, 111207, 111205, 111205, 111205, 111206, 111204, 111205, 111206, 111205, 111205, 111205, 111206, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111203, 111206, 111204, 111206, 111206, 111206, 111206, 111207, 111208, 111207, 111207, 111207, 111206, 111207, 111207, 111206, 111207, 111207, 111207, 111207, 111209, 111208, 111209, 111209, 111208, 111208, 111208, 111207, 111206, 111207, 111207, 111207, 111206, 111205, 111207, 111205, 111206, 111206, 111206, 111205, 111205, 111207, 111208, 111208, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111207, 111206, 111205, 111206, 111206, 111206, 111206, 111204, 111205, 111205, 111205, 111205, 111206, 111205, 111205, 111206, 111204, 111206, 111206, 111203, 111203, 111204, 111204, 111204, 111204, 111204, 111202, 111204, 111203, 111203, 111203, 111202, 111202, 111202, 111202, 111202, 111203, 111202, 111202, 111202, 111202, 111203, 111203, 111201, 111202, 111202, 111203, 111203, 111202, 111198, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111197, 111197, 111195, 111195, 111196, 111196, 111198, 111198, 111200, 111202, 111202, 111202, 111200, 111200, 111200, 111206, 111207, 111205, 111205, 111206, 111208, 111208, 111208, 111208, 111208, 111207, 111208, 111207, 111207, 111206, 111205, 111206, 111204, 111205, 111206, 111206, 111205, 111206, 111206, 111205, 111206, 111205, 111205, 111205, 111205, 111205, 111205, 111206, 111206, 111206, 111206, 111206, 111206, 111206, 111205, 111205, 111206, 111206, 111203, 111203, 111203, 111201, 111202, 111203, 111201, 111203, 111200, 111202, 111201, 111200, 111200, 111201, 111201, 111201, 111197, 111198, 111197, 111196, 111197, 111196, 111196, 111196, 111196, 111196, 111197, 111197, 111197, 111197, 111198, 111199, 111198, 111197, 111200, 111202, 111200, 111200, 111200, 111200, 111200, 111200, 111199, 111200, 111199, 111197, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111200, 111198, 111198, 111198, 111201, 111199, 111199, 111199, 111198, 111199, 111201, 111200, 111200, 111201, 111198, 111199, 111201, 111199, 111200, 111201, 111200, 111200, 111200, 111199, 111199, 111199, 111200, 111198, 111199, 111199, 111199, 111202, 111203, 111203, 111201, 111201, 111203, 111204, 111202, 111198, 111199, 111198, 111198, 111198, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111195, 111195, 111196, 111195, 111195, 111194, 111194, 111204, 111203, 111200, 111201, 111205, 111204, 111203, 111201, 111199, 111201, 111201, 111203, 111202, 111202, 111204, 111202, 111205, 111205, 111205, 111206, 111204, 111202, 111202, 111201, 111201, 111202, 111201, 111202, 111202, 111202, 111202, 111202, 111202, 111200, 111202, 111199, 111199, 111201, 111199, 111199, 111200, 111200, 111197, 111193, 111194, 111194, 111195, 111195, 111191, 111193, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111197, 111199, 111198, 111200, 111200, 111200, 111200, 111200, 111198, 111199, 111199, 111199, 111198, 111199, 111197, 111199, 111198, 111198, 111198, 111198, 111198, 111200, 111200, 111200, 111200, 111200, 111200, 111200, 111201, 111201, 111201, 111202, 111203, 111202, 111202, 111203, 111203, 111203, 111203, 111203, 111203, 111201, 111200, 111200, 111198, 111196, 111196, 111198, 111199, 111197, 111196, 111197, 111196, 111196, 111199, 111200, 111199, 111197, 111197, 111195, 111195, 111193, 111194, 111196, 111198, 111196, 111196, 111195, 111195, 111194, 111194, 111195, 111195, 111194, 111194, 111192, 111193, 111192, 111192, 111192, 111191, 111193, 111190, 111190, 111190, 111191, 111188, 111185, 111189, 111191, 111191, 111191, 111191, 111189, 111189, 111189, 111188, 111188, 111189, 111187, 111187, 111190, 111188, 111187, 111185, 111186, 111184, 111184, 111184, 111184, 111184, 111184, 111186, 111188, 111194, 111194, 111194, 111194, 111196, 111196, 111196, 111196, 111196, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111197, 111197, 111198, 111197, 111196, 111195, 111194, 111194, 111194, 111193, 111193, 111191, 111190, 111191, 111191, 111191, 111188, 111186, 111187, 111190, 111188, 111188, 111187, 111186, 111185, 111185, 111184, 111184, 111184, 111184, 111181, 111181, 111179, 111180, 111180, 111180, 111180, 111180, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111179, 111180, 111179, 111179, 111180, 111181, 111181, 111178, 111181, 111181, 111185, 111184, 111184, 111185, 111185, 111186, 111185, 111185, 111188, 111187, 111189, 111187, 111188, 111188, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111189, 111188, 111187, 111186, 111187, 111186, 111186, 111191, 111186, 111187, 111189, 111188, 111188, 111190, 111190, 111192, 111191, 111192, 111190, 111190, 111188, 111189, 111191, 111191, 111192, 111192, 111190, 111191, 111191, 111192, 111191, 111191, 111191, 111191, 111191, 111190, 111193, 111193, 111193, 111194, 111195, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111193, 111191, 111191, 111191, 111191, 111192, 111192, 111196, 111196, 111196, 111195, 111195, 111195, 111196, 111196, 111194, 111194, 111194, 111193, 111193, 111192, 111191, 111191, 111195, 111194, 111194, 111196, 111195, 111195, 111196, 111196, 111195, 111195, 111195, 111193, 111193, 111194, 111193, 111195, 111195, 111195, 111195, 111195, 111197, 111196, 111196, 111198, 111198, 111198, 111200, 111199, 111198, 111199, 111199, 111198, 111197, 111198, 111196, 111198, 111198, 111198, 111198, 111200, 111198, 111197, 111201, 111200, 111200, 111200, 111200, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111198, 111198, 111198, 111198, 111197, 111195, 111196, 111199, 111199, 111199, 111199, 111199, 111198, 111197, 111197, 111197, 111197, 111198, 111198, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111197, 111199, 111199, 111198, 111199, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111197, 111196, 111196, 111196, 111196, 111195, 111195, 111195, 111194, 111194, 111194, 111194, 111192, 111191, 111194, 111194, 111194, 111193, 111191, 111192, 111190, 111190, 111191, 111190, 111190, 111190, 111192, 111192, 111191, 111192, 111192, 111193, 111193, 111192, 111192, 111192, 111191, 111190, 111197, 111198, 111197, 111197, 111197, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111199, 111202, 111202, 111202, 111202, 111202, 111202, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111198, 111197, 111199, 111199, 111199, 111199, 111198, 111199, 111199, 111198, 111198, 111198, 111198, 111198, 111197, 111198, 111197, 111198, 111198, 111198, 111198, 111198, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111199, 111199, 111199, 111198, 111198, 111199, 111198, 111198, 111198, 111196, 111197, 111197, 111199, 111199, 111197, 111197, 111198, 111195, 111195, 111193, 111193, 111194, 111193, 111191, 111191, 111191, 111190, 111192, 111192, 111192, 111192, 111192, 111193, 111192, 111193, 111193, 111193, 111193, 111194, 111193, 111193, 111194, 111193, 111192]],
            'npts_0.5': [False, [64063, 64063, 64064, 64063, 64063, 64062, 64058, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64059, 64058, 64058, 64059, 64059, 64059, 64060, 64059, 64059, 64060, 64059, 64059, 64059, 64059, 64059, 64062, 64061, 64065, 64065, 64063, 64062, 64064, 64063, 64061, 64063, 64064, 64063, 64065, 64065, 64065, 64064, 64064, 64064, 64064, 64065, 64064, 64064, 64064, 64064, 64065, 64064, 64064, 64065, 64065, 64064, 64066, 64064, 64065, 64066, 64067, 64067, 64068, 64068, 64069, 64070, 64069, 64070, 64069, 64070, 64069, 64069, 64070, 64070, 64069, 64069, 64067, 64068, 64068, 64069, 64069, 64069, 64067, 64067, 64067, 64069, 64069, 64069, 64068, 64068, 64065, 64066, 64066, 64066, 64067, 64067, 64067, 64069, 64067, 64065, 64068, 64066, 64067, 64066, 64066, 64066, 64068, 64066, 64065, 64065, 64065, 64065, 64066, 64066, 64067, 64067, 64066, 64066, 64066, 64066, 64065, 64063, 64064, 64065, 64064, 64062, 64059, 64060, 64061, 64064, 64063, 64060, 64058, 64060, 64058, 64059, 64058, 64058, 64058, 64059, 64059, 64059, 64059, 64058, 64057, 64059, 64058, 64057, 64058, 64058, 64057, 64057, 64057, 64057, 64057, 64056, 64055, 64055, 64055, 64054, 64054, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64056, 64056, 64058, 64058, 64058, 64061, 64061, 64061, 64065, 64065, 64064, 64064, 64064, 64069, 64067, 64068, 64067, 64066, 64065, 64067, 64065, 64064, 64064, 64064, 64063, 64064, 64063, 64064, 64066, 64064, 64068, 64068, 64066, 64066, 64066, 64067, 64067, 64068, 64065, 64065, 64064, 64064, 64065, 64066, 64067, 64067, 64063, 64061, 64063, 64065, 64062, 64062, 64061, 64060, 64059, 64060, 64060, 64059, 64059, 64059, 64060, 64059, 64058, 64058, 64059, 64058, 64057, 64055, 64056, 64056, 64055, 64054, 64054, 64053, 64054, 64055, 64055, 64056, 64057, 64057, 64058, 64058, 64057, 64058, 64058, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64058, 64057, 64055, 64055, 64056, 64056, 64055, 64056, 64055, 64056, 64057, 64059, 64057, 64057, 64058, 64056, 64058, 64058, 64057, 64058, 64059, 64056, 64057, 64058, 64057, 64057, 64059, 64057, 64056, 64057, 64059, 64057, 64057, 64059, 64056, 64057, 64057, 64057, 64062, 64064, 64063, 64063, 64062, 64062, 64062, 64062, 64062, 64061, 64061, 64062, 64061, 64060, 64060, 64062, 64061, 64060, 64060, 64060, 64057, 64058, 64057, 64057, 64056, 64055, 64056, 64064, 64063, 64059, 64061, 64062, 64062, 64061, 64060, 64059, 64059, 64060, 64060, 64060, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64062, 64062, 64061, 64062, 64062, 64062, 64063, 64064, 64064, 64064, 64063, 64061, 64062, 64062, 64063, 64059, 64058, 64059, 64058, 64058, 64059, 64059, 64057, 64055, 64055, 64053, 64053, 64052, 64053, 64053, 64053, 64053, 64054, 64055, 64055, 64055, 64055, 64056, 64056, 64055, 64055, 64056, 64056, 64056, 64058, 64057, 64058, 64060, 64061, 64061, 64061, 64060, 64059, 64060, 64061, 64061, 64062, 64062, 64062, 64063, 64063, 64063, 64064, 64065, 64064, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64065, 64066, 64067, 64066, 64064, 64062, 64063, 64063, 64064, 64062, 64062, 64062, 64062, 64061, 64061, 64061, 64059, 64060, 64061, 64060, 64059, 64057, 64058, 64058, 64059, 64058, 64058, 64058, 64060, 64058, 64057, 64054, 64054, 64054, 64054, 64054, 64054, 64055, 64055, 64055, 64054, 64054, 64054, 64054, 64054, 64053, 64052, 64052, 64052, 64052, 64051, 64051, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64053, 64054, 64053, 64051, 64050, 64050, 64050, 64048, 64048, 64048, 64050, 64049, 64051, 64054, 64055, 64055, 64056, 64056, 64057, 64059, 64057, 64059, 64059, 64057, 64057, 64058, 64059, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64058, 64059, 64058, 64058, 64057, 64057, 64057, 64057, 64056, 64056, 64055, 64055, 64055, 64055, 64055, 64055, 64055, 64054, 64053, 64052, 64051, 64051, 64050, 64051, 64051, 64051, 64050, 64048, 64049, 64049, 64050, 64049, 64047, 64046, 64046, 64046, 64044, 64044, 64043, 64044, 64044, 64044, 64044, 64042, 64043, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64041, 64040, 64042, 64042, 64042, 64042, 64044, 64044, 64044, 64044, 64042, 64040, 64042, 64042, 64044, 64044, 64043, 64044, 64045, 64045, 64046, 64044, 64048, 64048, 64047, 64047, 64048, 64048, 64046, 64048, 64047, 64047, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64048, 64047, 64047, 64051, 64049, 64049, 64049, 64051, 64050, 64051, 64050, 64052, 64051, 64051, 64047, 64048, 64047, 64049, 64049, 64049, 64051, 64051, 64050, 64051, 64051, 64051, 64051, 64049, 64049, 64048, 64048, 64049, 64052, 64052, 64052, 64051, 64054, 64054, 64054, 64053, 64054, 64054, 64054, 64054, 64053, 64054, 64054, 64055, 64055, 64057, 64056, 64056, 64056, 64055, 64055, 64055, 64055, 64056, 64056, 64056, 64056, 64056, 64055, 64054, 64055, 64055, 64058, 64058, 64058, 64058, 64058, 64057, 64057, 64056, 64057, 64057, 64057, 64057, 64057, 64058, 64058, 64058, 64057, 64058, 64058, 64058, 64058, 64058, 64058, 64058, 64059, 64059, 64059, 64061, 64059, 64060, 64061, 64061, 64058, 64058, 64058, 64059, 64062, 64063, 64061, 64063, 64062, 64062, 64062, 64065, 64065, 64065, 64063, 64062, 64061, 64061, 64061, 64061, 64062, 64062, 64063, 64063, 64063, 64063, 64063, 64063, 64064, 64063, 64064, 64064, 64062, 64063, 64061, 64060, 64060, 64060, 64059, 64058, 64057, 64054, 64059, 64060, 64061, 64060, 64060, 64060, 64060, 64060, 64060, 64060, 64063, 64064, 64062, 64063, 64064, 64063, 64063, 64062, 64064, 64063, 64063, 64062, 64062, 64062, 64063, 64063, 64061, 64062, 64062, 64062, 64062, 64064, 64064, 64062, 64061, 64061, 64061, 64061, 64060, 64059, 64058, 64056, 64056, 64055, 64055, 64053, 64052, 64052, 64054, 64051, 64052, 64053, 64052, 64050, 64049, 64048, 64048, 64048, 64049, 64048, 64049, 64049, 64050, 64049, 64049, 64049, 64049, 64050, 64049, 64048, 64048, 64048, 64059, 64059, 64058, 64058, 64059, 64062, 64061, 64062, 64061, 64061, 64061, 64061, 64061, 64063, 64063, 64063, 64063, 64065, 64063, 64062, 64063, 64062, 64062, 64061, 64061, 64062, 64061, 64060, 64064, 64064, 64065, 64063, 64064, 64062, 64061, 64061, 64064, 64062, 64062, 64061, 64061, 64061, 64061, 64061, 64061, 64061, 64060, 64060, 64059, 64058, 64060, 64059, 64058, 64059, 64059, 64061, 64059, 64060, 64060, 64060, 64061, 64061, 64062, 64062, 64060, 64061, 64061, 64060, 64060, 64059, 64060, 64061, 64058, 64058, 64058, 64058, 64057, 64056, 64055, 64054, 64053, 64055, 64053, 64050, 64050, 64051, 64052, 64053, 64052, 64052, 64053, 64054, 64054, 64054, 64054, 64054, 64055, 64056, 64054, 64056, 64056, 64055]],
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

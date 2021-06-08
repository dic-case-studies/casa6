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

            self.save_dict_to_file(test_name,savedict, test_name+'_current_metrics')

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

       
        test_name = self._testMethodName
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

        exp_im_stats = {'com_bmaj': [False, 4.496111869844004],
            'com_bmin': [False, 3.3033772209261785],
            'com_pa': [False, 86.47054060326587],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 3.1939454078674316],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.3596886694431305],
            'min_val_pos': [True, [58, 39, 0, 490]],
            'im_rms': [False, 0.055886518611949515],
            'rms_per_chan': [False, [0.046436697869152804, 0.04934544120928089, 0.04896199177473545, 0.045242469255244916, 0.04076142987134027, 0.04066943677263002, 0.04376487990076443, 0.043694037974636946, 0.04459647157043778, 0.048017215146382866, 0.04627143010988964, 0.04420360428221782, 0.04156732965986629, 0.04050710499664779, 0.039057803588974566, 0.042761211802002805, 0.04252217160173592, 0.04359491616287709, 0.042110033059135944, 0.04292671323176175, 0.04664013025271547, 0.044981539437139055, 0.04665173264254359, 0.04681796405643853, 0.046889986578062995, 0.04456627223782911, 0.043437573543890357, 0.04696837963609141, 0.04119206881920388, 0.04253605374377071, 0.04488390666581614, 0.046219678502766974, 0.05084848031791911, 0.05094348922893772, 0.051587061193443084, 0.053067146840733756, 0.049471381348813376, 0.04595651140680717, 0.04495376242396089, 0.04038566412630896, 0.03864010620979754, 0.04199809692448168, 0.04636255623563074, 0.04548627909545541, 0.0416926662524104, 0.04004642564522245, 0.04359098502992989, 0.048897197549261646, 0.05023278382163298, 0.04587088821848791, 0.04197161169982609, 0.04629835772440621, 0.04930300708497102, 0.051480523822865304, 0.0444920176542373, 0.04123743560674529, 0.04782542643935964, 0.04746206516205945, 0.04416113994824711, 0.042678569765962016, 0.04258109174037231, 0.0415876900426464, 0.0460845913679538, 0.05012091495742385, 0.04691654308190901, 0.0476245504240818, 0.04773492099109901, 0.04495658663987826, 0.045487110777595444, 0.04562602278476289, 0.04696766248115374, 0.047969470203093195, 0.05112473591159104, 0.05068971697370774, 0.049147417386991965, 0.04507639438223545, 0.041421888684340255, 0.0425226081122096, 0.04666948693254447, 0.05137830022927088, 0.0487975705620606, 0.046813770813629384, 0.04353027014416515, 0.04118253262925014, 0.0434871216800935, 0.04636335671805942, 0.05095272918530583, 0.04726479357778318, 0.047112086285304035, 0.04810445474853609, 0.044180807950651574, 0.042093973726747065, 0.04295999527148268, 0.04124848250926838, 0.04070068168513572, 0.043414986596181405, 0.04661471496530699, 0.050460952779931165, 0.05008705270284653, 0.041392912138652954, 0.039904238787250336, 0.04679096971894066, 0.05151217574354204, 0.049486972953524395, 0.04519692999924898, 0.04440251398328859, 0.04924052476880247, 0.05214177625986392, 0.04692878097318393, 0.046059627017696655, 0.04612626064171623, 0.043618412970095506, 0.03976377438184463, 0.042115253722497, 0.04918542098706065, 0.04545362803301337, 0.044580289917133406, 0.043606589149843666, 0.04589475013866274, 0.046721917469477514, 0.045886670071669494, 0.04843013834859285, 0.0468518391646326, 0.04374076376976529, 0.04337661028501193, 0.04619439910530295, 0.05030786730965178, 0.04708647823192546, 0.04352457917341813, 0.04329424090647755, 0.04322675306877749, 0.044531887145507425, 0.04749786740181839, 0.05084068642775522, 0.04816327143578811, 0.04113284396067957, 0.0458824702579652, 0.04457602600259906, 0.042673950698155896, 0.04135518185294876, 0.04363527480521624, 0.045915816324020906, 0.042947043718610724, 0.042114974596657626, 0.040225926243584544, 0.03782476669855668, 0.040947824695726054, 0.04104108866239226, 0.04075579015241123, 0.03995035368137553, 0.04193807816156728, 0.04180721762719505, 0.047504768816046046, 0.04736784374922038, 0.044766797754126904, 0.04388574210339384, 0.04105473583139812, 0.046076845774234994, 0.0467644249374563, 0.048440853239654656, 0.05244702471581019, 0.053676932789258246, 0.045041508474699773, 0.044465034511806445, 0.050808335496450474, 0.04613109382684326, 0.04459087256845239, 0.04784074724761662, 0.04518331531207736, 0.045309429358903475, 0.04667744939767361, 0.039019394101063434, 0.03614463130987937, 0.041436322331324854, 0.04292124772257666, 0.04052206582901182, 0.04111202794825397, 0.04516745212774294, 0.047704624181569434, 0.047580045401103596, 0.04717948381553734, 0.04412987455554941, 0.04256697937153531, 0.04181606458168034, 0.04907721236107465, 0.049563858911260326, 0.0421301373909515, 0.04059486020625589, 0.04019192818638815, 0.0447132231289818, 0.047353768692884754, 0.04758399846175147, 0.046136697711113436, 0.040747266079874625, 0.04076262287036722, 0.04067101148976313, 0.040799759360021115, 0.04299405385428427, 0.04377280211034057, 0.04345150159208221, 0.044980372184520695, 0.04392953325826627, 0.04367090267082839, 0.047776860736536246, 0.048070363879887126, 0.04672031680965515, 0.04724325319701497, 0.04407948891411036, 0.04546998614564801, 0.04455527757983371, 0.046142136341844255, 0.04320616054281921, 0.045534956659338044, 0.04859329360282889, 0.048900190734707506, 0.04929999273771512, 0.04632358475481116, 0.04428851750274284, 0.04784975169537624, 0.04765093319911423, 0.043881745818891524, 0.043324369056045746, 0.04514575976732948, 0.045891383907318166, 0.047016150372676854, 0.045037014585648856, 0.04010251969147668, 0.04159175015205475, 0.04137216214036925, 0.04275113249661427, 0.04186177273686212, 0.042649167994782675, 0.05214059329636489, 0.04997177400191877, 0.04413320467298178, 0.042953912737445386, 0.041919859361076166, 0.04215158603777731, 0.04172864652945504, 0.045968321694947474, 0.04875739844937148, 0.045170973253897566, 0.04516555305609894, 0.045270464081727156, 0.04802056986652718, 0.04808519866194299, 0.049062815490268495, 0.049016631820180526, 0.045425672802811776, 0.04732161236976712, 0.04346775339199203, 0.04153855085303185, 0.04206187024253007, 0.04332446409616924, 0.049482501906829714, 0.0507054466430013, 0.044307426639965806, 0.04543831583909989, 0.047296688638528, 0.04643888426830375, 0.044806402864451686, 0.040267103806810325, 0.04115060022819575, 0.04706710046145555, 0.046288412762023334, 0.04357966717296147, 0.04520833664590388, 0.0445773127076585, 0.042527732196947624, 0.04125083538851448, 0.0464798273652813, 0.0477061093380914, 0.04734605419206378, 0.04816610067157907, 0.04420919983233172, 0.04293696948203265, 0.04200859180941565, 0.04215638496043083, 0.04388855282240666, 0.04582580951984011, 0.046035236630698695, 0.041593507191228, 0.04355778141475355, 0.05088599406941826, 0.04940169882218437, 0.04859676205214566, 0.04856506220858988, 0.04463197009403143, 0.0469928195631672, 0.04892660260584388, 0.04188083544330787, 0.04133702742116763, 0.04165002972255051, 0.04465837542426475, 0.04524045036081385, 0.04333489557469641, 0.036989245194290137, 0.04199435378560823, 0.0481421046288897, 0.04716445679719833, 0.04154312168845833, 0.04176322011685239, 0.04046233927208671, 0.039883571656071926, 0.03765969072148476, 0.041965142187759395, 0.04559238206427811, 0.0453795609369959, 0.045403985370455834, 0.043728555125363625, 0.04582630825421479, 0.04574769870819116, 0.041828903815262365, 0.044431679448605176, 0.04407619669845318, 0.044814214081711067, 0.04848916851754347, 0.04465432486321319, 0.04394693625127326, 0.04702186055677961, 0.04973534695469247, 0.049146429948471695, 0.049174047326561644, 0.04641722732162844, 0.0413493581116718, 0.040667117853418484, 0.038859433299798685, 0.04364937330851356, 0.04420011801442241, 0.044643107270842856, 0.04811818762813673, 0.04687352375691469, 0.04256380947059253, 0.04359082466003849, 0.04230980866473369, 0.04281115946568443, 0.04422596088989216, 0.042246858498182424, 0.040077278402510984, 0.042549874146797714, 0.04405511204341645, 0.045906337049994056, 0.0465523248225122, 0.046191054551640416, 0.044661499298400485, 0.042991036849962264, 0.043753909902568, 0.04161174584570825, 0.04519950697840892, 0.046265539276719315, 0.04268460646762688, 0.04193615173426503, 0.039888183799370425, 0.04520227446386294, 0.044925270927308386, 0.04304844469474453, 0.043324892275943376, 0.042493393018593686, 0.044060400863577764, 0.04833598244473396, 0.04667128831700136, 0.04400275471906582, 0.044034772926916044, 0.0449817248307724, 0.044117220953356004, 0.04653280231209841, 0.045651252724683906, 0.047787009258557236, 0.04695171563115938, 0.04586997367936294, 0.04653594966524015, 0.046746865286752294, 0.04596274133397725, 0.043423193820342536, 0.0429095051365923, 0.04261310963062516, 0.04079337018589888, 0.044095264125143696, 0.046447871571280504, 0.04604306292098572, 0.04133170707137652, 0.0420198631265, 0.045862996429123136, 0.042615400719355566, 0.04152397719066135, 0.041999641834113935, 0.044506562146456125, 0.0489076227548345, 0.050881675342586384, 0.04865983581349746, 0.0466263016917132, 0.039185250702358655, 0.039570364674090364, 0.040831878878307054, 0.04122134099737277, 0.042839947818260435, 0.04505271727089972, 0.04371049447279272, 0.04496821842378058, 0.043527303075767916, 0.04146834265759304, 0.0413392604938528, 0.04095808278187557, 0.04185026895871849, 0.042819533882487534, 0.044015774767846225, 0.04472485535881441, 0.04597126224550101, 0.04755869539738559, 0.04507188732391839, 0.04305694521144407, 0.0407981952355781, 0.04368203628914822, 0.046542715674692335, 0.05031311018545435, 0.049150553827940056, 0.050420119426495796, 0.05042319018370619, 0.048208632662477956, 0.041953421435925224, 0.04338425755527854, 0.04425432187285405, 0.0469656376516665, 0.04772128602161581, 0.04422046207371832, 0.0431820595706624, 0.041876735808617165, 0.048004314058412025, 0.04900912687093005, 0.04344985549827121, 0.0472484655849604, 0.04287398199229198, 0.044514729537670016, 0.043700386055600185, 0.04342698342556917, 0.04791522096395403, 0.04903780430384516, 0.045944967610429314, 0.0431244319437859, 0.042030569388848776, 0.046092789879613134, 0.04236907961853478, 0.039984657655645614, 0.04224856764482613, 0.04382141692853359, 0.04496437538277435, 0.051192810669175526, 0.05216376798749919, 0.04532434451335352, 0.0431429415470322, 0.035891139856947835, 0.03777495237969397, 0.04456327531280985, 0.04843266965180281, 0.043500726406838025, 0.0395594486007869, 0.04261777241666637, 0.04579923208972242, 0.045738668555702314, 0.04372132290200405, 0.044647710224980056, 0.04249636238690276, 0.04015298023581446, 0.04619275637175317, 0.04612509644727192, 0.04810497095694837, 0.048132106400979105, 0.04395205678082682, 0.04458138815599986, 0.0431716555054271, 0.04391610135390437, 0.04570699536249715, 0.04330642753889957, 0.040729093816539746, 0.04974366860566736, 0.045342501262102326, 0.04167534163587883, 0.042307335955378014, 0.04232635407330191, 0.04414200953135981, 0.04588069278483478, 0.044283750950693354, 0.04369716817887697, 0.04179124447894311, 0.0466318762331663, 0.07806544022832206, 0.15776378798837767, 0.27120535650240335, 0.35558286881232487, 0.36936842458936664, 0.33274253987994623, 0.29442340895699953, 0.27098045350482386, 0.25221602783952335, 0.24069468654698972, 0.2484813999384745, 0.2589627963754489, 0.2717264009097436, 0.2841120694868153, 0.2632075644478858, 0.19501431775504555, 0.1153205374797238, 0.05942547189735459, 0.044180540316820296, 0.04319301077871217, 0.0431325767114206, 0.04085491190301253, 0.0440834670662473, 0.04857361140230655, 0.045041767961530546, 0.04179778558721278, 0.040806152865383055, 0.04204786435876296, 0.045247925684986014, 0.04734484038570823, 0.04844246502217243, 0.04744621613609996, 0.047019488078712035, 0.045119054192661455, 0.0428316286870027, 0.04277679574653031, 0.04587281553677757, 0.04624714875887038, 0.04754205960961853, 0.0506659414706044, 0.046054933440241363, 0.04211559712901001, 0.04627071885689092, 0.049834614275276155, 0.05046074730273139, 0.047724282773597866, 0.041494289753345064, 0.04431996966311481, 0.04843417805839743, 0.04814998715103537, 0.04907738913499242, 0.04944282864619265, 0.04521747648573062, 0.04219988564277459, 0.04285237591833291, 0.038456315088847524, 0.04011490406736227, 0.043681775665527174, 0.045537973609178314, 0.045569222781213066, 0.04354643976503561, 0.04368748023632953, 0.04419588815529541, 0.046819346881899856, 0.04462637334638743, 0.04521544542145996, 0.04557127316946643, 0.04505243150891467, 0.04323891857602188, 0.04467411762339298, 0.04361373829523538, 0.039607602510523594, 0.04052242922576435, 0.046030992267336805, 0.04662190935211449, 0.043745211699971386, 0.04350194213269832, 0.04350759565774762, 0.04275327682417688, 0.04543882275179457, 0.043127001714826606, 0.044817849789228176, 0.04545663518730103, 0.04647914055836658, 0.04096308925012331, 0.04179146339703641, 0.045368231933929415, 0.047729927579010734, 0.05029062044525315, 0.046690435074565884, 0.0501755327667634, 0.047333358370090374, 0.04001043698156056, 0.04595460077786111, 0.050040341006851105, 0.046937775103807214, 0.04531478433611846, 0.04566360634952669, 0.044947005667432906, 0.044808455821824904, 0.041235928977575415, 0.03968341783187546, 0.040143871244530484, 0.04375295607460526, 0.04904403283029426, 0.049564362573488084, 0.04443334861343526, 0.04274910708086753, 0.040725183572862526, 0.04583017072214125, 0.04878966115305859, 0.04128730507404871, 0.04412335208453002, 0.046495354201513926, 0.04138182706024895, 0.04217817788506256, 0.04624709060618901, 0.04431090481619349, 0.044972112128921334, 0.046304769603916586, 0.0485513869194507, 0.04424903407642042, 0.04428396754117315, 0.04404354841945006, 0.045999757978571455, 0.04579974483101019, 0.04493050206105136, 0.04782962784941458, 0.04727426490033995, 0.04380589997033528, 0.04485366326335109, 0.04807820755332521, 0.04920756428479516, 0.04937645285679352, 0.049829580867415915, 0.048764293156774385, 0.04829510191146219, 0.04693900499211445, 0.04633655498260907, 0.04175610883019444, 0.03708686949259042, 0.0445156121518753, 0.049115336033239775, 0.04712301573002867, 0.047231432170403646, 0.04348383422002165, 0.04321561745224448, 0.04382385998504035, 0.04040596100029039, 0.04224179133545749, 0.04468346149294344, 0.03988319187104351, 0.03828601116171196, 0.04419977254195108, 0.042373853288575375, 0.03831681970812281, 0.03655782091817188, 0.042670253436652504, 0.046708351984574965, 0.04473726624213793, 0.04612807135571747, 0.051406361767195756, 0.05297824157494852, 0.04513533169978514, 0.03839869776421337, 0.04053378357153958, 0.04163322724858157, 0.04182891691882434, 0.04842081242381634, 0.04879504210169032, 0.046836217179586116, 0.04290682817489914, 0.044582341784591105, 0.04490003205842483, 0.04602142844487166, 0.045481230956658146, 0.04538314395623165, 0.04223445460877744, 0.04402352206900099, 0.045266977321175725, 0.0436593613008295, 0.041516160487384525, 0.039738897911894416, 0.04569307821046952, 0.04733008502119898, 0.041957993557838885, 0.042021219098555405, 0.04792601786683623, 0.05004385861635714, 0.04833875608957245, 0.04665437443805753, 0.044992448559262266, 0.03975975299026254, 0.043846557978077016, 0.04828840176505269, 0.04652716565752805, 0.04730529551337945, 0.04498167597122633, 0.04216975907784085, 0.03886523559977788, 0.04194369924654704, 0.045727387422739754, 0.044506117005971034, 0.04410919693599897, 0.04062528253553941, 0.039909011598196674, 0.03746370169882932, 0.041062210165122744, 0.04562615759438729, 0.045414085338632955, 0.041010795636350655, 0.04286251865671152, 0.03797225257453629, 0.04330598203899954, 0.04610718044218543, 0.04625332381908609, 0.04588018597990095, 0.04098233795751277, 0.04075419076561464, 0.04131496593419855, 0.04202835123632931, 0.045677278987934396, 0.04427965467407332, 0.04321194076704215, 0.04039862653270404, 0.044175641167403176, 0.046434921571469986, 0.04612819801366069, 0.04694776780472301, 0.04099383222212228, 0.03917241219341224, 0.042068449460125036, 0.04370712675052841, 0.04566609615291663, 0.04554385827900137, 0.04281206035046255, 0.04347090687004217, 0.04112938716216996, 0.04202862504150548, 0.04295046737687642, 0.04545657385475635, 0.046957344320490174, 0.044930510000153544, 0.04168760492953715, 0.04090601771168612, 0.0455630776972825, 0.04597795049825335, 0.04133169705894654, 0.04120487628598337, 0.04758132208426841, 0.04750359649922572, 0.043746754411552416, 0.039969611653241026, 0.0370851372143889, 0.04089311315935916, 0.04083945567726639, 0.0428484838520684, 0.039332551408404214, 0.03858436205685368, 0.040934957347207916, 0.050654546446170026, 0.05097260586321339, 0.041853153510787375, 0.04338997701915169, 0.04678955322620283, 0.045415538614583896, 0.0442047946569195, 0.04389453231034864, 0.04625370428755419, 0.04914433760878362, 0.04640838863763117, 0.04675583353174989, 0.048344256903468295, 0.04769398098282866, 0.04409896357510309, 0.039703360252422194, 0.04255562166678417, 0.04217126611200596, 0.04003525253301322, 0.042949092700844, 0.04162791047900747, 0.0413162126868823, 0.04016158990534426, 0.04080225782419112, 0.04504404215012396, 0.047715803434991544, 0.04459195804658081, 0.04758174926498982, 0.04747017594826984, 0.0475656917361288, 0.046386447740224386, 0.04725118171029761, 0.04808790481799425, 0.04858139145256644, 0.046394274503655984, 0.04546742468429808, 0.04473916398152629, 0.04091007186311367, 0.040473323570406845, 0.04149693391689712, 0.04489012854384838, 0.04568975487582137, 0.04308233314534105, 0.04428562845682951, 0.04811203918428873, 0.047285389364052244, 0.05154243012970362, 0.04874546793919434, 0.046481953117862844, 0.0467410772482252, 0.04540861377997277, 0.04768072021944487, 0.04668965399420149, 0.04534680854356319, 0.045717731495190275, 0.04352592050803243, 0.042246781377357774, 0.04297653390841251, 0.04482593340392208, 0.04213803084079387, 0.040341402839927605, 0.041202696841319074, 0.04455581844575452, 0.04583316901662954, 0.04476653987017234, 0.04404044993113527, 0.04323584216597428, 0.0407148005579045, 0.039810485948214555, 0.04344594505571544, 0.04735161553566336, 0.04856055968640867, 0.04691866083987017, 0.047895854140470075, 0.04881197291412326, 0.045525706148063026, 0.046417520159948655, 0.04474300103353705, 0.04155988246327171, 0.04377189815845601, 0.04795878974600023, 0.04525235421707231, 0.0387289509478405, 0.03922893116237577, 0.03803757186243102, 0.03949120817560298, 0.040143851753887225, 0.043114834889568, 0.04602567054512971, 0.0437991817004047, 0.04106584766908839, 0.042641255317866854, 0.04550032987759433, 0.04234298059824302, 0.041622781448595915, 0.04692339278596035, 0.044528609287418454, 0.041743387912223796, 0.042208493978792685, 0.041125083988986934, 0.0405759128144399, 0.04467213327019776, 0.043214303511136115, 0.041508134050783296, 0.041428173047111556, 0.040282008128922166, 0.040203688250344996, 0.039737248806300016, 0.04051657728302586, 0.041622680469507414, 0.04221313659186694, 0.040082710053825375, 0.039719538697112, 0.04411726176802357, 0.04282786046798159, 0.04242926571118004, 0.041994588868336645, 0.046654630912404196, 0.049532951890798115, 0.04732551195144287, 0.04615435180286521, 0.04630096966294292, 0.04405921913621422, 0.04057920168856714, 0.043127866489089076, 0.04159706690252945, 0.0409040743352284, 0.044183587281843774, 0.04617277543544684, 0.04745291146776343, 0.05001258828501942, 0.044527038550667745, 0.04147790210671815, 0.04015206913630134, 0.0424791772606839, 0.045002347470638235, 0.04536894399192883, 0.043285293689377224, 0.04377719366412504, 0.04372734733000784, 0.042147854310142716, 0.038394077278512836, 0.03841638903877585, 0.040738262575049415, 0.04406275824470088, 0.04561539364919089, 0.04815410288701247, 0.04380385040469506, 0.03927703898847133, 0.041332503806580274, 0.04346874721018446, 0.04111359760443814, 0.04189287559788935, 0.03992992170852107, 0.041589487789998575, 0.04253549033084363, 0.045843748226115796, 0.04003225597131588, 0.03787007823901111, 0.03799338617254116, 0.04285197349285763, 0.045624047198232466, 0.04472829925619543, 0.047179414501805184, 0.04686834076548325, 0.042443768278091494, 0.047928754185026366, 0.049687264773723455, 0.04154905997274937, 0.04218387341604758, 0.04520406326908012, 0.04139958233783231, 0.04019798332273786, 0.04190762617444021, 0.04335420732360746, 0.04237465977091424, 0.04628070745735236, 0.04495196078901107, 0.046447062304978685, 0.050760228472102153, 0.04652902848369043, 0.04255680214098841, 0.04008157254800408, 0.03890717013039061, 0.04252967926532384, 0.04688988230610397, 0.04741546640409401, 0.046850148636497245, 0.041072400781220014, 0.03556251474091502, 0.03901253247042472, 0.03899358613011951, 0.045269905619281584, 0.04224687380267042, 0.03805147308971664, 0.041744101599233294, 0.04474441172714239, 0.04724508885906723, 0.049218613302406414, 0.04768860292698267, 0.04518556813818921, 0.04405214586861989, 0.04411513272767838, 0.0481235909145108, 0.04739636465817227, 0.04211249577864894, 0.04266132751769585, 0.04161774305336725, 0.043100364700299126, 0.04552503274158557, 0.04993780254601025, 0.04782835142613495, 0.04378025606649709, 0.042449558484278845, 0.03810731472722694, 0.039204816964426396, 0.04219617471250491, 0.045762028874738465, 0.04494528717024618, 0.04589857251270448, 0.0472983276430976, 0.04757687497726176, 0.04903733760320237, 0.04938431289090874, 0.04204901564635354, 0.04359603587681168, 0.048005009794489566, 0.04582925273166158, 0.0425409211761743, 0.04749963733807478, 0.047810341596607894, 0.046306134109708144, 0.04328346075754866, 0.045757773687285064, 0.045906416400406, 0.048505022928582706, 0.04447474104571062, 0.047842325597475706, 0.050649418080458926, 0.043857959451520837, 0.03573407062144151, 0.04277427524098398, 0.04295178228999916, 0.047484579340476736, 0.04679772086423353, 0.043423465631174406, 0.0458444488504807, 0.04353395605574517, 0.041243884477055695, 0.04250780220611674, 0.045617517133681645, 0.04877347310740958, 0.04461816696526135, 0.04073852230707331, 0.04631327857005907]],
            'im_sum': [False, 2271.747356080856],
            'regn_sum': [False, 252.51129265129566],
            'npts_real': [True, 6400000],
            'profile': [False, 2.7224358130019817],
            'fit': [False, [3.0141833841972825, 6.00061793161013, 5.668161770679299]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.50497926661194],
            'fit_pix': [False, [45.783391223533, 41.07442247498776]]}

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
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'mask_pix': [False, 9642],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        exp_pb_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20036059617996216],
            'im_rms': [False, 0.576840993910233],
            'npts_0.2': [False, [3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233]],
            'npts_0.5': [False, [1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549]],
            'npts_real': [True, 6400000],
            'fit': [False, [1.046847676114786, 28.075049566294457, 28.075049566280292]],
            'fit_loc_chan': [True, 500],
            'fit_loc_freq': [1e-10, 354.50632237101223],
            'fit_pix': [False, [40.0, 40.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        exp_psf_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.16418913006782532],
            'min_val_pos': [True, [36, 35, 0, 993]],
            'im_rms': [False, 0.08711620311345131],
            'im_sum': [False, 2742.739240222416],
            'npts_real': [True, 6400000],
            'fit_0': [False, [1.1041815481079804, 3.942786862512436, 2.861422654767207]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.44539427140006],
            'fit_pix_0': [False, [39.99991697297065, 39.997612723971876]],
            'fit_1': [False, [1.1041897949047454, 3.943995482479015, 2.858120721885849]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.50632237101223],
            'fit_pix_1': [False, [39.99991481982509, 39.99762628119395]],
            'fit_2': [False, [1.1041649735486678, 3.946028634858572, 2.8548012143390697]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672504706244],
            'fit_pix_2': [False, [39.999910983933844, 39.997642687543916]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 0.3189396262168884],
            'max_val_pos': [True, [59, 65, 0, 492]],
            'min_val': [False, -0.32551121711730957],
            'min_val_pos': [True, [32, 47, 0, 491]],
            'im_rms': [False, 0.046506649333909536],
            'im_sum': [False, 236.25628862426532],
            'regn_sum': [False, 51.38772555207834],
            'npts_real': [True, 6400000]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.3655270338058472],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.07262927293777466],
            'min_val_pos': [True, [58, 39, 0, 490]],
            'im_rms': [False, 0.001713471051900203],
            'im_sum': [False, 52.686264319345355],
            'regn_sum': [False, 53.27634758129716],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400000]}

        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1009.5013427734375],
            'max_val_pos': [True, [0, 0, 0, 979]],
            'min_val': [False, 1007.9752807617188],
            'min_val_pos': [True, [0, 0, 0, 0]],
            'im_rms': [False, 1008.7183259897186],
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

            self.save_dict_to_file(test_name, savedict, test_name+'_stats')



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
        exp_im_stats = {'com_bmaj': [False, 4.555337133782812],
            'com_bmin': [False, 3.366170346656177],
            'com_pa': [False, 84.29991156839789],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 3.293868064880371],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.3665587902069092],
            'min_val_pos': [True, [32, 47, 0, 490]],
            'im_rms': [False, 0.05691195242659064],
            'rms_per_chan': [False, [0.0511710509742412, 0.0491058318709313, 0.047427317804676385, 0.04594435148674381, 0.0415009859700181, 0.04235219163430191, 0.044665370448235396, 0.04396920906263267, 0.044964299741928296, 0.046866284052463285, 0.04511061287632169, 0.04513810215355222, 0.04341798783967012, 0.043050257310306664, 0.042129786010037105, 0.04636869210044085, 0.044059859711359506, 0.0431005878298246, 0.040842230625417586, 0.044253567123580514, 0.04788789162384652, 0.046765812215769344, 0.04847899075342864, 0.04515956443260174, 0.04538329846919214, 0.04399170567475163, 0.04102228636632494, 0.04847222769658979, 0.043927446546726236, 0.04390410197613647, 0.04780115032582894, 0.0466626671608103, 0.04777581073270312, 0.046921663192160436, 0.04744313453465746, 0.04824639384380612, 0.04788860970633028, 0.04880504778837709, 0.046702319349569796, 0.04136799544094345, 0.040308435767093816, 0.04260226577052781, 0.04547176054951916, 0.04286001122690497, 0.040840505618497565, 0.04128661785358394, 0.0404748729255245, 0.042631249975814, 0.04540595536141528, 0.04630780990356906, 0.04297565065452147, 0.042575767914947674, 0.04552191770673945, 0.050248180795149845, 0.046241521425769844, 0.04330203798503876, 0.04768733950182284, 0.04847564341337296, 0.044981581100583645, 0.04189074558356269, 0.04188414937705219, 0.04079794880186283, 0.0470115336511376, 0.04986855588643256, 0.047600572328423164, 0.04718238766642814, 0.047370448573660746, 0.04499853260000665, 0.04583102311333413, 0.04539363938016473, 0.04264057951288134, 0.04583432526287408, 0.05071215052060615, 0.051313592516402856, 0.04942251962800088, 0.0452659277052202, 0.0458722320785561, 0.04613611570301598, 0.04465261606782503, 0.04829763157869588, 0.04782813524298008, 0.04494709830097677, 0.042381751521326096, 0.04141890895233068, 0.04467059934599844, 0.04559961356719576, 0.04547567220115387, 0.04435778612501103, 0.046562842912044225, 0.04578554313014109, 0.04331591067646111, 0.04242197635982062, 0.04343724664810054, 0.04163437621955441, 0.04287598365908082, 0.04555346099364086, 0.04597304625448051, 0.05200878188299578, 0.05414682485856498, 0.045479063677547865, 0.04216966989655045, 0.041714049861264584, 0.045120142864486165, 0.04857578684261745, 0.04716268927025188, 0.046611049794305984, 0.048399908635185306, 0.04744089384493704, 0.045090408799761825, 0.049220554347253656, 0.047350986737763306, 0.039940624100381614, 0.03680860477408705, 0.04059962123891759, 0.04907408990430075, 0.04664142745221579, 0.046223478522023874, 0.04467307845563109, 0.048889824810990946, 0.04749560034065197, 0.04585869037339095, 0.046341359374244856, 0.04496882756089031, 0.04347350826006669, 0.04366455133429647, 0.04510413768546546, 0.04767605563780636, 0.04628777540015753, 0.04383844882584038, 0.04557909476816198, 0.04448225097973088, 0.04461731749963628, 0.04871954222966844, 0.05128764741906331, 0.048084815499627044, 0.042300818407911896, 0.04625528537619014, 0.04336890388492618, 0.04038118961870588, 0.04047602846973458, 0.042567214889519256, 0.04515217595189108, 0.044129137309465494, 0.04307067417311321, 0.04043559745933244, 0.03818958694565964, 0.039706852487315224, 0.042979483899121214, 0.043916140730896176, 0.040544206455762705, 0.043326762913106526, 0.04206060328268939, 0.046981420807857076, 0.04845391512559523, 0.04672629320983973, 0.04404219681984842, 0.04172994421619655, 0.04459099831496948, 0.043866258959796615, 0.047358150403602146, 0.049035864489749684, 0.05094280776928182, 0.04492396944906997, 0.04400046101968693, 0.04811706718667277, 0.0426986402937336, 0.04467005664630853, 0.04954293591254802, 0.04443597582811842, 0.04545636238646333, 0.048839934131079746, 0.04315819427142755, 0.0386181639950708, 0.04402914248162857, 0.046187749772956525, 0.041035310855728035, 0.04099317590057528, 0.04441207051714462, 0.046235847532745226, 0.0482010704543739, 0.04392063938988062, 0.04064756145475974, 0.04138550275559402, 0.040243183655705554, 0.04518851623269384, 0.045694083704164334, 0.040105531014094836, 0.03987778125691593, 0.04093442367509905, 0.047282985328581696, 0.04806981059082117, 0.04716012661373992, 0.04561097088491241, 0.04064006261976807, 0.04116720745451455, 0.0406780229779118, 0.04265757275994845, 0.04376732363903251, 0.04454225656204078, 0.043386922587584184, 0.04217701011036419, 0.04332828002160028, 0.04266178634176287, 0.04551570766907826, 0.0480457262655882, 0.04561792734222954, 0.04409469369039748, 0.04057848039162825, 0.04236357519774322, 0.046349524387768265, 0.05030629016867036, 0.04485247188838552, 0.04449155134718874, 0.04992276056798944, 0.04940551822489045, 0.04868295725037189, 0.04527712738694331, 0.04186862261555057, 0.0429630524634968, 0.04499982893318418, 0.04304105923032144, 0.044339172561306595, 0.04739996899328926, 0.04501416069400519, 0.044652280084299306, 0.04412519586902521, 0.04135166306993902, 0.044132282709036866, 0.04231884284017298, 0.04469923871251993, 0.04473574988832509, 0.042215441874413284, 0.04990394042300807, 0.0494353630705107, 0.044832618004111466, 0.045785099166762495, 0.045199260880778115, 0.04508466465276481, 0.043361469974022196, 0.04658883658994638, 0.050564414356763086, 0.046578192297974046, 0.04585650570372835, 0.04586022834671364, 0.04464739885130335, 0.046741192518259464, 0.04996942633464429, 0.047683534337113816, 0.042319162274250886, 0.04671506777733548, 0.04409193227846734, 0.0422249069485366, 0.04467846094556274, 0.04509045753574932, 0.04804217831929826, 0.047362437062488175, 0.041702495866345, 0.04208843910529905, 0.04386681352310239, 0.04527932894324497, 0.04508196034274706, 0.038003708577452125, 0.04132926646915876, 0.0465432753615036, 0.04560711306130088, 0.04576331308915141, 0.04656752899110721, 0.04360772210798942, 0.0399468128147913, 0.040997556628144825, 0.04935035279770661, 0.050949630321551524, 0.04495109163647917, 0.04127036696343021, 0.041545778001206186, 0.04298065481299039, 0.04180053472952686, 0.04053515744138335, 0.044116665636433926, 0.047682248394217, 0.04578460093590072, 0.042317505080424145, 0.04112408141815903, 0.04950290975819657, 0.04906809053165872, 0.04712610814404197, 0.04483025733112527, 0.04220710327441495, 0.04689904031040545, 0.04770977093977536, 0.04061446274621193, 0.03855395424343832, 0.04103417297823304, 0.04270719013150538, 0.04311744934020913, 0.041787918368279865, 0.03950595765143592, 0.044205208171978934, 0.04804693977027296, 0.047143593610045664, 0.044883052091199194, 0.04671632857686066, 0.04610526965682995, 0.04074178180069125, 0.0381349688144203, 0.042360964011485404, 0.04382490315637503, 0.044781962713744, 0.045960278135158655, 0.04416267514227893, 0.045058342659130744, 0.04408765168781989, 0.04308487473515519, 0.04610518692950167, 0.04302606342567099, 0.04411797662376259, 0.047592576874859766, 0.04479946723379238, 0.04314016463758305, 0.045018361314295346, 0.048779415499810325, 0.04719537393597005, 0.04686169568890435, 0.04573342413095565, 0.04136888072681148, 0.04239595421044405, 0.040295033961612645, 0.043324138427140765, 0.042886299345549495, 0.04370150121377847, 0.04726720119047192, 0.04688152497663171, 0.04216566673700439, 0.040303789036477126, 0.040888455245350994, 0.04056560828708437, 0.041015157238264835, 0.03831574147172779, 0.03776969003744252, 0.04341599590916081, 0.044613221687851944, 0.04378484155531375, 0.044154545418134135, 0.04405293992229615, 0.04555903572340968, 0.0450423604604058, 0.04560720061928255, 0.045134063340471854, 0.048497031964785664, 0.048328710055805815, 0.04369752093987805, 0.04387499374878124, 0.045679250975076385, 0.04828545607378901, 0.04503792573426137, 0.04394604355265449, 0.042199847559304636, 0.04314833855323105, 0.04483643841068481, 0.048007724008618906, 0.04581600864856695, 0.04258056593979621, 0.0420584098318127, 0.04355085221133504, 0.042962592132744094, 0.046111055050587504, 0.04395769623502132, 0.04468357289629069, 0.04697969938166983, 0.04595402000023419, 0.04697419409395668, 0.049219344252601376, 0.048005086759996325, 0.05029377733754566, 0.049082308875115686, 0.04267820296145758, 0.04094884844716048, 0.04297025758607257, 0.046713378508162706, 0.048087964531066105, 0.04444586038351215, 0.04368701270899749, 0.04742946140995044, 0.04514395002893858, 0.04391551654271971, 0.04123374048433164, 0.0457787324442786, 0.049064605477177343, 0.04490913205635445, 0.04490228337451022, 0.047401908608303144, 0.041911945003784895, 0.040215578543547716, 0.041307267340881325, 0.045368858694375895, 0.04628760318753789, 0.0452571349151798, 0.04299040266972557, 0.04566401299145563, 0.042186218440877665, 0.039590605280172816, 0.040980868671193174, 0.04163775576388983, 0.04130413769013361, 0.043945668906054836, 0.04561215824721322, 0.04445602404327032, 0.045783098980149074, 0.04778146396244574, 0.04622256875039857, 0.04429939325010576, 0.04222737410728622, 0.04480220999980045, 0.04573713589666081, 0.04807642278115695, 0.05012661487110346, 0.051057847293619764, 0.049202903288461886, 0.044066078590053515, 0.03869330697784115, 0.04005750370613939, 0.042583087272016096, 0.04646433953331963, 0.04843656672825028, 0.04337159868478386, 0.04129161674523401, 0.04159979426678676, 0.04783796708212357, 0.04803509551594159, 0.04271292352623876, 0.0455489830132127, 0.04155378961917943, 0.04315514135818257, 0.04390220146264162, 0.04410930681949159, 0.0472019930773177, 0.04874645246192576, 0.044078769591148086, 0.04046239100958424, 0.04323542240548025, 0.046749951792652615, 0.04413279360486385, 0.03991442364566927, 0.042099668239468174, 0.04211068440864971, 0.043505387166735715, 0.047323848003642095, 0.04745554969935124, 0.0425372910672586, 0.043669644783594995, 0.03865510143522782, 0.04013774912311121, 0.042941638685498495, 0.044228962278241336, 0.042187197445406895, 0.040637885246025025, 0.04176379743472241, 0.046266099841557926, 0.04835892509809472, 0.04440581513278651, 0.044848534784280254, 0.04332244038665112, 0.04095861338972233, 0.04477945828350456, 0.04474446370841756, 0.04672014871005834, 0.04890087092035712, 0.044920945647659714, 0.04353612087720862, 0.04285252796775038, 0.04154177855921692, 0.04569167250444369, 0.043389563137158964, 0.04108722155650539, 0.04890667415484547, 0.04430280757082574, 0.042555447692747425, 0.043954675540407925, 0.04473920869437501, 0.045953716831692866, 0.04611529686424024, 0.045527063359704516, 0.04437703908042696, 0.04152190513121626, 0.04588997114091835, 0.08070828360453258, 0.17011597131696088, 0.2873350500968064, 0.37507381962364805, 0.38993932248535357, 0.34765529897938036, 0.30792990574742735, 0.28128455623042625, 0.26295627858727993, 0.2544366049193439, 0.2596753362283072, 0.27214812217106316, 0.2825564132466614, 0.3015568378388069, 0.2758339466259943, 0.2056344972324021, 0.12362288409410535, 0.06400462469485993, 0.04362897942331456, 0.044045856282146134, 0.043864817316223895, 0.04108521207383964, 0.04507105433252286, 0.04895275229699207, 0.0473583518597183, 0.04168874503456669, 0.04202879335770104, 0.04473483090647776, 0.045188066866670505, 0.04676576183010721, 0.05038043826705083, 0.04588896571630465, 0.04622699335322239, 0.04592513568585445, 0.04195009884133858, 0.04137032401295382, 0.04521559106259956, 0.04966261459892102, 0.04758965817729824, 0.045636840033258866, 0.04417985445759443, 0.043269319580502646, 0.04935094731153835, 0.05370886483877212, 0.052221461037257215, 0.04831597277551329, 0.04271439613298134, 0.04209441992597077, 0.04650469698846066, 0.04962336274206145, 0.051016831892572785, 0.04743221239707722, 0.04568975014812525, 0.04510746645325057, 0.04400632333172186, 0.039998180475841076, 0.044621910311143, 0.04609853642164451, 0.04522393365538562, 0.04490998250978899, 0.043346339569221785, 0.043457016809540165, 0.0425993943925095, 0.046717444841513855, 0.04541188524705667, 0.045851992759998454, 0.04519253299037234, 0.04346783288282716, 0.043424134416619686, 0.04374210351649922, 0.04211478761613377, 0.039450140949390845, 0.0405040633318485, 0.045318930207222954, 0.04435671796233326, 0.043546793546948286, 0.043575865027333495, 0.045246238459986016, 0.041681659896991234, 0.043331177539123356, 0.042239183489564694, 0.044955958036043096, 0.045374649862824266, 0.047828787509047156, 0.04232759040225429, 0.04110862159028894, 0.044151143962815505, 0.047682404696488066, 0.04996482173530259, 0.046759940781707345, 0.05150270065420599, 0.048914591535924874, 0.040131921828050365, 0.04723972314983509, 0.0506400595414927, 0.0478873987229166, 0.047566739808015825, 0.04723356186385765, 0.045710170646954096, 0.046457969765637104, 0.0442617450705403, 0.040953448770182484, 0.040839513639018374, 0.042784628654277394, 0.04464862031842329, 0.04597979973084751, 0.04441919498390057, 0.042193966876797266, 0.04198679442123217, 0.04320787012340684, 0.0484379621838361, 0.04183535948295165, 0.04494146633002023, 0.048966414872948505, 0.042046722384187205, 0.03891622722221963, 0.042053695569046486, 0.039838669658182654, 0.04141147790843661, 0.04458325311829123, 0.04735379934500237, 0.044437525423351526, 0.04510354102167274, 0.04412668192463012, 0.04310266863107471, 0.04275755168127876, 0.04246356363697789, 0.046138842164377226, 0.04367898256061182, 0.04065078652072923, 0.041095669529102166, 0.04266922755691068, 0.043871990393221365, 0.04538459538693181, 0.04657986803965163, 0.04591004079602359, 0.046475106624329146, 0.04524991461308665, 0.04294724823755841, 0.040834275730254596, 0.037567738401474766, 0.04488960175150954, 0.04999511866122657, 0.04626706093539809, 0.045619597740860114, 0.044068027471420856, 0.04624198737247629, 0.04430221593730956, 0.03992111296904295, 0.04041328132554766, 0.04291991423048277, 0.04070631516567948, 0.03967701840711301, 0.042280483357614966, 0.03807842858645987, 0.038253135332009185, 0.04031507349837287, 0.04253909225308235, 0.04502400749590674, 0.04657554806832086, 0.04695604364284647, 0.048191751162209855, 0.051234660662775516, 0.0451858743732422, 0.040620794530063314, 0.04250481945172377, 0.04503577071485718, 0.04434892065177995, 0.04986518589887667, 0.04805347503299002, 0.04663912561261378, 0.045511813628137686, 0.04539590766902981, 0.04619076921941353, 0.048657856053205266, 0.04835269980457421, 0.045699917213065745, 0.04176513491526662, 0.043166407129944774, 0.04274112383516181, 0.03974718655782916, 0.04306056539817251, 0.04452636935738953, 0.04850082064279787, 0.04708455557123165, 0.041004121067491205, 0.04212389205777563, 0.04813538396706355, 0.050250948335518686, 0.04676145050138597, 0.043736102971985595, 0.040845160557585684, 0.03748637318500156, 0.04337656627557361, 0.04749475263440817, 0.046295219404072135, 0.048143619541977684, 0.04488974005750239, 0.040573526719343694, 0.037986684846188605, 0.0411242047327654, 0.045419802025673355, 0.04399195151277112, 0.04481399144552655, 0.042982118616007044, 0.045008250330539894, 0.04529351269873727, 0.04784826443602788, 0.04692498277483518, 0.04421288548266149, 0.04110530878634437, 0.04267277479799294, 0.038170622192963966, 0.044482197900083356, 0.04675206064562118, 0.04375449822428938, 0.04333647628267722, 0.04108680149011656, 0.04003212690634986, 0.040520752932678726, 0.04126731341982271, 0.04332949788908975, 0.041824390226529035, 0.04222385525269832, 0.03989487789220953, 0.04138708970795852, 0.04245894216308035, 0.044725873951922124, 0.04721825696276477, 0.04202813192653578, 0.0399460528657829, 0.044571106380753886, 0.0447323229086075, 0.044050067662238006, 0.04690696307911725, 0.04312504591995337, 0.04535449794212105, 0.04175185169290581, 0.04157695497192835, 0.042609494573312895, 0.04625807567746221, 0.050335731719151525, 0.04726435688804032, 0.04234447811560612, 0.03982635414749611, 0.042865048792595194, 0.04545209156888121, 0.042407326449704304, 0.043206572383087435, 0.04800092199739377, 0.04829777325509949, 0.04482363982150042, 0.04206650083122971, 0.04255687985599084, 0.04429451586584098, 0.041962086054368866, 0.043919826900260674, 0.03934343424164113, 0.037813965833677535, 0.04142744265610587, 0.051316829486620696, 0.049805988318551306, 0.041625950162602135, 0.04338778236760687, 0.048651992345158844, 0.0458107323638186, 0.04177986842250171, 0.041755522528424016, 0.04459573786058077, 0.04726495154821514, 0.04456404711820412, 0.04605781488693276, 0.04820277873090167, 0.046173531413243735, 0.04330581712817806, 0.040119221546603766, 0.0418880880070373, 0.041528662619733044, 0.043534078906845075, 0.04659356346381995, 0.044376291876952176, 0.043186166807959975, 0.040499731261696446, 0.04155453510073291, 0.04476806930534024, 0.04608939080278736, 0.04246457199963066, 0.044811740665129836, 0.04636203862946733, 0.04809773817125487, 0.047074582717482175, 0.04706955716115251, 0.04881264031749564, 0.048407027251639534, 0.04709565523751361, 0.04547126227968264, 0.042303449988877063, 0.04136077275872563, 0.04305732524780119, 0.04127375854913234, 0.04630413725126096, 0.04759718637329068, 0.04467009557051697, 0.044047586246901194, 0.046733636823707424, 0.04472946107601905, 0.050777975000485906, 0.05125934174063112, 0.045936831440501924, 0.043707393859863754, 0.04311598030391774, 0.04599191804260334, 0.04798843524134782, 0.04490284314987071, 0.04777219395782297, 0.04969023200395486, 0.04740904203539907, 0.0440009652574982, 0.0440677574873923, 0.040270676991024414, 0.041189248650477454, 0.03996385955441793, 0.043444626747977264, 0.04732738458665101, 0.04474067829417449, 0.043869244275380526, 0.041537677380427324, 0.03997971990014071, 0.03912187617621613, 0.041581602406675335, 0.04628875507974483, 0.04787231535426587, 0.04662323393957114, 0.05003078944788062, 0.049874018140278344, 0.04353304807270341, 0.04369661636198581, 0.04392774186515209, 0.04166153900515943, 0.041484546214130044, 0.04502193469191308, 0.04427284993262843, 0.03912138117134029, 0.03919150643985177, 0.04087309848399349, 0.04255661402141116, 0.04411047509452905, 0.04365853958835935, 0.04434415264072045, 0.04376116550989087, 0.04057536714967884, 0.04354290123808518, 0.04571809685548106, 0.042792338786651646, 0.043098490502763244, 0.049854409965534925, 0.04618636534124665, 0.044123660069196165, 0.042086391478348084, 0.04048235134370972, 0.04251342608958859, 0.046362598362469004, 0.045539527217078185, 0.042719248349796464, 0.042142071143352035, 0.04027273308296864, 0.03955791363749769, 0.03944671741997897, 0.039158461234095734, 0.04266214242670252, 0.042423747609326996, 0.04125398444071295, 0.040742964126010424, 0.04189489106697355, 0.04037726362998526, 0.04231378678619918, 0.04290878280725253, 0.0465167727837056, 0.04790521430006423, 0.04919630398490165, 0.05001619884001604, 0.04542644389291285, 0.040298221493403064, 0.04127718613267144, 0.04354987117807297, 0.04080881213017432, 0.04328542626191592, 0.04643675404306652, 0.048183512607549876, 0.049184734318566646, 0.04840071684592878, 0.044204873534443294, 0.0420089699437857, 0.039827562677114736, 0.04178428137270748, 0.044483539959811885, 0.045202770846765855, 0.04282458887249053, 0.04311984281040187, 0.043317366061525156, 0.040266420533826563, 0.03746444088907071, 0.03876012103208295, 0.041242722451112476, 0.0453531604402292, 0.04651219150499251, 0.049103871858467436, 0.046289432392027416, 0.0418003626070388, 0.042559230557822264, 0.041694666797074585, 0.03977286339079858, 0.042047934721173366, 0.04325667764783359, 0.04478087241686524, 0.043467564075746284, 0.046467472384519405, 0.040069313414233365, 0.03932669033135957, 0.036674081932373, 0.04035573321611562, 0.043825746315377785, 0.04226006500936928, 0.04426403927091294, 0.04428239796857482, 0.04007931179126735, 0.046470043730921604, 0.04861436337918293, 0.041802373836734054, 0.041731738317003275, 0.04327330997242376, 0.03962533653138047, 0.039946629679497796, 0.04350867084763192, 0.043463668193889166, 0.04007075045958809, 0.04561667281300996, 0.04574049416764589, 0.044864232866475985, 0.048466440806165946, 0.043636359901250456, 0.03984589129147291, 0.03908402303474198, 0.03785282550309543, 0.04522641562273479, 0.05128945098693023, 0.05052140630211563, 0.04569754159001448, 0.04224544042769382, 0.037523115503082, 0.04088410465817124, 0.04123471656317448, 0.04274657538287841, 0.04276422107388948, 0.04091079128023806, 0.04231664731089127, 0.04755323530150363, 0.04936026674954882, 0.05017393358130251, 0.0480175491831683, 0.046315161255912146, 0.044508746245854595, 0.04409176859115316, 0.04652053544581553, 0.04379336598194394, 0.04157993210004191, 0.04474370008153558, 0.046010065609617984, 0.04716694473114898, 0.04641569830102659, 0.04849923245678476, 0.04784957531124401, 0.04220093989806986, 0.04083499581574191, 0.036692163257465785, 0.039565171768156694, 0.04390553432840742, 0.04396788774254892, 0.04207350687266713, 0.0424163824454376, 0.04468932647219615, 0.04857659362984081, 0.05032269061240432, 0.049439548889266965, 0.0418949978731385, 0.04164303251149575, 0.04715653986013338, 0.04760769719247853, 0.04227993046284273, 0.046795607676726574, 0.048698666879539707, 0.049068165913337534, 0.04410003502982793, 0.04711619369250108, 0.04587679228814314, 0.0475623617876886, 0.04619988912165085, 0.04922100826271814, 0.052913329773182534, 0.0467717887991848, 0.039225326729021785, 0.04433580305969694, 0.04280667653011399, 0.046413255513419543, 0.04715064792703495, 0.044548670643813876, 0.04796392427546589, 0.04527100121963938, 0.04321478812841821, 0.042187416713115364, 0.04305771697722239, 0.04569364972904935, 0.04310032933326111, 0.0381080331502587, 0.04519094524636647]],
            'im_sum': [False, 2541.550916265518],
            'regn_sum': [False, 269.4618835449219],
            'npts_real': [True, 6400000],
            'profile': [False, 2.806910952210095],
            'fit': [False, [3.110410414279016, 6.176995156924437, 5.7534551319905995]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.50497926661194],
            'fit_pix': [False, [45.65467797815349, 41.05254434252133]]}

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
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'mask_pix': [False, 9852],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_pb_stats
        exp_pb_stats = {
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20036059617996216],
            'im_rms': [False, 0.576840993910233],
            'npts_0.2': [False, [3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233]],
            'npts_0.5': [False, [1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549]],
            'npts_real': [True, 6400000],
            'fit': [False, [1.046847676114786, 28.075049566294457, 28.075049566280292]],
            'fit_loc_chan': [True, 500],
            'fit_loc_freq': [1e-10, 354.50632237101223],
            'fit_pix': [False, [40.0, 40.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_psf_stats
        exp_psf_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.16455988585948944],
            'min_val_pos': [True, [34, 36, 0, 949]],
            'im_rms': [False, 0.08653503618208624],
            'im_sum': [False, 2857.8199424562085],
            'npts_real': [True, 6400000],
            'fit_0': [False, [1.095566785329151, 4.020892595865358, 2.9706317101525186]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.44539427140006],
            'fit_pix_0': [False, [40.00007105295384, 39.99750303438477]],
            'fit_1': [False, [1.0955778815758896, 4.020331875996908, 2.9703856624539005]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.50632237101223],
            'fit_pix_1': [False, [40.00007131349202, 39.99750438784483]],
            'fit_2': [False, [1.0956200407624843, 4.019976254809749, 2.9692129137171084]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672504706244],
            'fit_pix_2': [False, [40.00007226924553, 39.99750678468819]]}

        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_pcwdT.exp_resid_stats
        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 0.3233320415019989],
            'max_val_pos': [True, [20, 53, 0, 497]],
            'min_val': [False, -0.3440611958503723],
            'min_val_pos': [True, [32, 47, 0, 490]],
            'im_rms': [False, 0.04664051079977191],
            'im_sum': [False, 331.069520983562],
            'regn_sum': [False, 51.19770630635321],
            'npts_real': [True, 6400000]}


        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_eph_pcwdT.exp_model_stats
        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.4161585569381714],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.07062584906816483],
            'min_val_pos': [True, [58, 39, 0, 491]],
            'im_rms': [False, 0.001756862244366767],
            'im_sum': [False, 55.41836966946721],
            'regn_sum': [False, 56.15662528015673],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400000]}


        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_cube_eph_pcwdT.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1011.9310302734375],
            'max_val_pos': [True, [0, 0, 0, 70]],
            'min_val': [False, 1011.7631225585938],
            'min_val_pos': [True, [0, 0, 0, 475]],
            'im_rms': [False, 1011.8442614231698],
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
        exp_im_stats = {'com_bmaj': [False, 4.555308268539054],
            'com_bmin': [False, 3.365998340500264],
            'com_pa': [False, 84.3021528064608],
            'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 3.2901856899261475],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.3735768795013428],
            'min_val_pos': [True, [32, 47, 0, 490]],
            'im_rms': [False, 0.05690275351685566],
            'rms_per_chan': [False, [0.05116729711532761, 0.04910328853904836, 0.047426418701247246, 0.04594330791478071, 0.0414993552853678, 0.0423505809213826, 0.04466391431221988, 0.043967366538047975, 0.044963184071386615, 0.046865525000930226, 0.04510974334062459, 0.045136820839991335, 0.04341679043913987, 0.04304826549628102, 0.042126962919566205, 0.04636587345482038, 0.04405803763461302, 0.04309976376065565, 0.040840793499483045, 0.044252013742321836, 0.047887870296788884, 0.046764589337172605, 0.04847735064746266, 0.04515817689816514, 0.04538150763387064, 0.04398946623855847, 0.04102177779879191, 0.04847160231514866, 0.04392692511659066, 0.04390215584266637, 0.04779868528096379, 0.04666063084450012, 0.04777377982019708, 0.04691930613995507, 0.047441610926343136, 0.048245599052981564, 0.04788910491952359, 0.048803735800761784, 0.04670042526063481, 0.04136639616777576, 0.04030704600712542, 0.042601238108909305, 0.045470583305565054, 0.04285820957438361, 0.04083882505699726, 0.041285577562955936, 0.040474623940935296, 0.04263010153067033, 0.04540432070607112, 0.04630602417284128, 0.0429746424844384, 0.04257500816494042, 0.04552042269976565, 0.05024578770460155, 0.046238668317632826, 0.043299559215245594, 0.04768508511368846, 0.048473069808366215, 0.044979066244736715, 0.041888204543854675, 0.04188138044589389, 0.04079656242788451, 0.04701025880153929, 0.04986681145895155, 0.04759934598593846, 0.047182637071428786, 0.04737029057006218, 0.04499798886886377, 0.045830492437916384, 0.045393173538764306, 0.04263933129597192, 0.04583134301637154, 0.050710037763646236, 0.05131220638333853, 0.049421180618205156, 0.04526535786873854, 0.04587137359358024, 0.04613435085464911, 0.04465050414869148, 0.048295791554697216, 0.04782637095324365, 0.044945369045715944, 0.04238044018406295, 0.04141802997818246, 0.044669561101233254, 0.04559871326629393, 0.045475820004178415, 0.04435791148011699, 0.04656188251072096, 0.04578388853573801, 0.043313908016516044, 0.042420522095935005, 0.04343673451507639, 0.04163288329698661, 0.04287399631600737, 0.04555165210911033, 0.045972303866938506, 0.0520067139548869, 0.05414397034129014, 0.045477161769417, 0.04216805692639514, 0.04171368363497993, 0.04511947281370644, 0.04857451421925222, 0.047162626066362706, 0.046610426623432503, 0.04839892239719139, 0.04744000855074887, 0.04508954843254543, 0.049218898729637375, 0.04735052852414996, 0.03994086525734225, 0.03680807594182824, 0.04059802557102472, 0.04907177102581404, 0.04663967226487628, 0.046221688188371045, 0.04467138420835264, 0.0488870498699097, 0.04749303063784769, 0.04585638384646808, 0.046340088557040665, 0.0449692535448143, 0.043473232734277015, 0.04366214249309008, 0.04510211367015539, 0.0476739565677098, 0.046285425210523194, 0.0438362110941993, 0.0455771382039262, 0.04448063725848597, 0.04461522245603117, 0.048717303055004504, 0.05128674815534695, 0.04808358290360343, 0.04229888417203775, 0.046253408410955064, 0.04336772593087985, 0.04038111190719556, 0.04047644797468002, 0.042566326017542774, 0.045150770987751505, 0.04412659675607936, 0.043067842482071066, 0.040433927490619075, 0.038188013744777496, 0.03970523796212651, 0.04297749439458233, 0.04391372479960849, 0.040542730072032615, 0.04332445414084403, 0.04205791715402477, 0.04697902315753677, 0.04845145900248261, 0.04672439733294359, 0.04404054749582631, 0.0417283901460653, 0.044590691437377074, 0.04386615482754906, 0.04735643378516962, 0.0490346467112115, 0.05094191025864634, 0.04492263012597952, 0.0439987963665548, 0.04811597080599062, 0.042698213141347904, 0.04466812382968773, 0.049541451489508656, 0.044434733284454575, 0.04545464443745161, 0.04883760999443205, 0.04315560270187941, 0.03861569202243875, 0.04402629591850435, 0.046185667414745536, 0.04103342752145166, 0.04099088624288794, 0.04440990529811676, 0.04623431830922685, 0.048199762751454454, 0.04391946116703184, 0.040646146268483244, 0.041383999065549446, 0.040241623506358605, 0.0451877006739886, 0.045693713559264354, 0.04010548773507098, 0.03987747934721885, 0.040933220036877924, 0.04728224506368913, 0.04806952295864596, 0.04715820775205383, 0.045609080346869224, 0.04063899648712911, 0.04116553360560871, 0.040677424962669656, 0.042656528116295035, 0.04376496689769152, 0.04454045209136651, 0.043385854797762864, 0.042176955213025315, 0.043327529527520246, 0.04266039003469953, 0.04551353196986531, 0.048042476018018455, 0.04561557483186576, 0.044093983908934556, 0.040577076286147255, 0.04236165958012301, 0.04634772243116468, 0.05030507789509709, 0.04485043075625129, 0.04448957351936966, 0.049921236379178904, 0.04940394416253078, 0.04868163632266521, 0.045276159371444934, 0.04186768961609199, 0.042962619401706266, 0.04499874246747334, 0.04303980078643947, 0.04433739144122113, 0.04739878279484431, 0.045013070236028435, 0.04465090677750713, 0.04412373599412178, 0.041350761968379425, 0.04413152043979249, 0.04231813977434867, 0.044697973066269284, 0.044733079449081656, 0.04221292312171885, 0.049902738947518704, 0.04943451574019778, 0.0448310339601354, 0.04578332817064797, 0.045197004135803624, 0.045081348335432175, 0.04335874376181495, 0.0465869499721449, 0.05056216409040827, 0.046575568509305186, 0.04585429096259659, 0.04585842605312868, 0.04464686522739467, 0.0467400768167714, 0.04996796041008005, 0.04768328026539277, 0.042319349737677474, 0.04671450416204197, 0.04409012943641373, 0.04222298996228015, 0.04467634333740302, 0.04508817276600825, 0.04804137837347428, 0.047362208399825044, 0.04170161594378143, 0.04208715741910808, 0.043864629804226715, 0.045276986319875305, 0.04508036760975704, 0.03800269893013875, 0.0413272210638029, 0.046540409201921186, 0.045605298741620885, 0.04576253485466565, 0.04656704271109055, 0.0436069696252561, 0.039945647173359786, 0.040996280688094786, 0.049348450214585664, 0.05094742705051873, 0.04494995990133932, 0.04126985904652916, 0.04154431030703643, 0.04297880294700325, 0.04179852472526678, 0.04053303795393168, 0.04411398235894323, 0.047679531818191775, 0.04578200966604946, 0.042314939144074366, 0.04112240353315196, 0.049501223897322325, 0.04906628288924678, 0.04712503578312127, 0.04482958555846782, 0.04220731130879237, 0.0468995045390844, 0.04770981318513866, 0.04061386447342675, 0.038552792876544714, 0.04103227158168726, 0.042704742999142724, 0.043115723575801096, 0.04178645552433733, 0.039503696868969373, 0.04420390820893296, 0.04804491853337866, 0.04714149981886613, 0.044881146395607974, 0.04671514999059933, 0.04610424955968707, 0.04074110041067899, 0.03813303121310113, 0.042359326620172735, 0.04382406797901019, 0.04478028970717848, 0.04595887326111854, 0.04416052365507004, 0.04505544088060466, 0.04408602625168516, 0.04308312548287552, 0.046102623377645584, 0.04302475305749556, 0.04411744969475481, 0.047592987643520976, 0.04479884234058751, 0.04313888360782502, 0.04501734215777683, 0.048777900659180716, 0.04719463324073782, 0.04686053113486287, 0.0457318123274814, 0.041366972887140384, 0.042393372513732125, 0.04029250979526522, 0.04332321589885846, 0.0428860882004734, 0.043700173470180105, 0.04726499818834088, 0.04687906462198162, 0.042164448744996974, 0.04030279837757205, 0.04088750010051021, 0.04056480312057034, 0.04101396020350358, 0.03831494269839726, 0.03776852241273864, 0.04341461869097921, 0.04461174731501267, 0.04378303295871123, 0.04415256586733245, 0.044051095937795125, 0.0455573185689383, 0.04504084407299739, 0.045605725537191104, 0.045131990916430144, 0.048495478200801766, 0.04832721933337852, 0.04369600973526794, 0.04387415810530338, 0.0456787556424287, 0.0482841134168499, 0.04503638445241318, 0.04394433607791088, 0.042197338669291866, 0.043146801869541264, 0.04483405695685402, 0.048004073064772786, 0.04581415761436921, 0.042578986065324995, 0.04205705354025357, 0.04354942578034775, 0.042961454434658854, 0.046110373876017546, 0.04395638231047666, 0.044681857041917294, 0.04697752511454902, 0.04595232785227561, 0.04697195063233949, 0.04921719863300786, 0.04800311683691384, 0.0502914658012708, 0.049080071305503306, 0.04267724319831937, 0.040947641550522454, 0.04296850650970498, 0.04671051762711886, 0.048084258779988225, 0.04444361787740532, 0.04368598990051458, 0.047427486020162994, 0.04514122918227242, 0.04391312991625713, 0.0412325979086211, 0.04577780367515201, 0.04906199418776229, 0.04490787040577596, 0.044902305631628975, 0.04740113122963352, 0.0419103137563522, 0.04021369447180656, 0.041305998775877324, 0.04536651777353465, 0.046285215913978696, 0.04525552054177476, 0.042989346441429474, 0.045661716737883815, 0.04218495026990458, 0.03959053147507961, 0.040980487248962166, 0.04163741286266258, 0.04130334836380524, 0.043943763617896044, 0.04561039685446372, 0.04445528488105397, 0.04578242647887811, 0.047778807047666176, 0.04622034706354849, 0.04429850830245711, 0.04222642385165843, 0.04480015058060434, 0.045734559252952924, 0.04807404480263517, 0.05012428367681061, 0.05105605554538161, 0.04920209840684151, 0.044065329806392196, 0.038692138695392, 0.04005610025336422, 0.04258073437893235, 0.04646186055215285, 0.048433928146113206, 0.043370226538408975, 0.04129038654167129, 0.041598121526794776, 0.047835863101601664, 0.048033875461006764, 0.04271221281159591, 0.04554792952420345, 0.0415522529932429, 0.043152884515414526, 0.04389968952379298, 0.044107684979145405, 0.047201018090718146, 0.0487453858125432, 0.044077711881765236, 0.0404606880973439, 0.04323348628569985, 0.04674907384540912, 0.04413263166289932, 0.03991372966969408, 0.0420977890669076, 0.04210899051239342, 0.04350308561141041, 0.04732264024476415, 0.04745564274457208, 0.042536306726198146, 0.043667732844827244, 0.03865358412673667, 0.04013637952908064, 0.04294116972023388, 0.044228597443581055, 0.042185109956491204, 0.04063573130945376, 0.04176346385714941, 0.04626539531938366, 0.04835732946584149, 0.044405109577706506, 0.044847545886691334, 0.04332039301742537, 0.04095626744624328, 0.04477842957323933, 0.04474270966369732, 0.04671807460864309, 0.04889853448749098, 0.04491759805337659, 0.04353300829061404, 0.04284958351691704, 0.041539539721888986, 0.0456899551854273, 0.04338830634916568, 0.041085728750192174, 0.04890482347719496, 0.04430101176815802, 0.04255564991587231, 0.04395392374396662, 0.04473671961834208, 0.04595060411498081, 0.04611206145827532, 0.04552378431768198, 0.04437480149200546, 0.041519953084110245, 0.045887086449703775, 0.08070184930184159, 0.17010130888346403, 0.28731363235925883, 0.37518382871348555, 0.3892949430156687, 0.34719364976424666, 0.3078523995738813, 0.2812591806886139, 0.26293234156257417, 0.25441423474359065, 0.2596529746362085, 0.27212458938195516, 0.2825308396246516, 0.3015067161449368, 0.27588177170354, 0.20557769337275264, 0.12361313059557615, 0.06399998855090439, 0.04362733856685197, 0.04404372601386455, 0.04386212073928898, 0.041083416718789054, 0.0450698028920436, 0.0489506651610377, 0.04735644686143907, 0.04168738281845572, 0.04202630918442599, 0.04473263738953591, 0.04518664952293885, 0.04676527885711022, 0.05037878502988409, 0.045887840411189804, 0.046225890916499965, 0.045923548272185116, 0.04194891373012146, 0.04136925844204895, 0.04521463179560795, 0.049660081211829445, 0.04758793776748543, 0.04563749265618054, 0.04417805327634696, 0.04326665821138246, 0.04934834142218267, 0.05370581757698511, 0.05221953158343308, 0.04831487631584759, 0.04271272447294109, 0.04209199276836319, 0.04650228903648591, 0.04962175749457209, 0.051015647237626215, 0.04743232128091371, 0.045689953379644016, 0.04510627113800185, 0.04400429197580577, 0.03999592219790395, 0.04462047027630851, 0.04609776347140558, 0.04522302504042985, 0.044908773531102134, 0.04334460653337687, 0.04345562043831261, 0.04259918849825984, 0.046718062120882406, 0.045411412550706076, 0.04585131104052212, 0.04519123184158126, 0.04346589994978304, 0.04342194896693256, 0.04374109326333788, 0.04211361929092754, 0.03944878888166286, 0.04050334438619292, 0.04531807124028844, 0.04435500695751117, 0.04354473464530376, 0.04357456372210875, 0.04524479441767891, 0.04168088158063653, 0.04333036517986282, 0.042238130696952515, 0.04495482588819793, 0.045373235462835, 0.04782733021721891, 0.042325836165534265, 0.04110668733004752, 0.04414955618172679, 0.047681669277635834, 0.04996316231863413, 0.04675767429764171, 0.05150060758469386, 0.048912136361160846, 0.040130035480996126, 0.04723787618862338, 0.050638729217014844, 0.04788674757532179, 0.04756502310460671, 0.047231228303738104, 0.04570845894664183, 0.04645639459337434, 0.04426033263371052, 0.04095195034899301, 0.04083790216639175, 0.04278262046018793, 0.044646211130660994, 0.04597753292560869, 0.044417165169960214, 0.04219176006863156, 0.04198454298338819, 0.04320673978883932, 0.04843640780502863, 0.04183375986968685, 0.04494023711710365, 0.04896550858783924, 0.04204504784589675, 0.03891437574954094, 0.04205280481680479, 0.03983740222405669, 0.04140986031086472, 0.04458156617725206, 0.04735151997510863, 0.04443560282422139, 0.04510204799965698, 0.0441253377336145, 0.043101549595038965, 0.042756441514509747, 0.04246240882663583, 0.046137880474099266, 0.043678039065423765, 0.040650834802730185, 0.0410956660194223, 0.04266933916516755, 0.043871678059890624, 0.04538342487996103, 0.046579600467178, 0.04590930188815202, 0.04647332343303196, 0.04524825618563781, 0.04294660021725802, 0.040833579003721955, 0.037566145809399806, 0.044888893446522936, 0.04999472047604279, 0.046266559267686524, 0.045619635720472854, 0.044067332150996355, 0.04623960384585643, 0.044300532421850114, 0.039920421649818905, 0.04041223789442004, 0.042918749466879444, 0.040705458404786815, 0.03967494081583607, 0.04227877507665494, 0.038077198119929184, 0.03825145402108632, 0.040312613834383594, 0.04253746903694269, 0.045022399617389554, 0.046572992277021406, 0.04695284513078954, 0.04818922665529497, 0.05123339319546236, 0.04518482124935362, 0.04061884207713111, 0.04250283264840193, 0.04503427398657955, 0.04434723268132182, 0.04986428424921361, 0.04805228945178296, 0.04663744846311706, 0.04550970679929021, 0.04539356025786973, 0.04618918417072768, 0.0486586927397824, 0.04835199058089579, 0.04569785627605178, 0.04176343558839675, 0.043164788725333796, 0.04273979825665297, 0.03974652444574154, 0.04305890119936653, 0.044524312426168694, 0.04849818613205884, 0.04708282979857598, 0.04100243099883566, 0.04212131087248006, 0.048132645471550425, 0.05024909613108227, 0.04676084814207293, 0.043735361130646015, 0.040843740872219425, 0.03748511419802302, 0.0433744333017367, 0.047492457425993104, 0.046294000690351354, 0.04814330748940135, 0.04488926512527577, 0.0405720825559782, 0.03798507429819145, 0.04112199698827656, 0.04541800693427667, 0.04399103227789008, 0.044811924939795536, 0.04297995548613756, 0.04500612419237003, 0.04528999135299485, 0.04784543262526196, 0.04692332440061222, 0.044211520860239995, 0.041103166249598255, 0.042670726426144556, 0.038169338024399216, 0.04448112369429627, 0.04675022431271521, 0.04375299422614655, 0.0433349760053098, 0.04108428453875046, 0.04003054323237877, 0.04051890823955529, 0.04126487998966233, 0.04332784406450002, 0.04182300821482669, 0.042223373931770954, 0.03989379513887651, 0.04138599158886148, 0.042457626246167175, 0.04472529482201079, 0.04721720039101929, 0.04202689797253457, 0.039944641948859584, 0.04456923698524046, 0.04473106791374529, 0.04404942986606497, 0.046905347594676565, 0.04312306523518605, 0.04535138950795937, 0.04174908896764267, 0.04157587339533971, 0.042608430372081355, 0.046255857137618026, 0.050334019816933494, 0.047262673974153624, 0.04234290159930733, 0.039825039860077925, 0.04286342904753436, 0.04544938622257436, 0.042405256687195414, 0.04320519061617756, 0.047999415130281546, 0.04829639741606494, 0.04482090085083877, 0.04206381528350167, 0.04255443908864056, 0.0442921434882811, 0.041960066400263896, 0.04391831361706683, 0.03934235574357688, 0.037812808264869054, 0.041425603196271324, 0.05131482840843496, 0.049804711443060444, 0.041624710266980496, 0.04338621809312778, 0.048650060299028484, 0.04580910452338551, 0.041778427124963075, 0.04175323739724728, 0.0445940895206801, 0.047263607157656456, 0.044562619096912295, 0.04605684752999258, 0.048202720997610034, 0.04617320157810235, 0.0433045049680528, 0.04011819867415498, 0.04188806053778063, 0.041528772970416516, 0.04353236779918611, 0.046591376195629995, 0.04437447229596645, 0.043184541900372894, 0.040498547619413926, 0.04155340684034046, 0.044766101868641356, 0.046087324852919956, 0.04246403404491755, 0.044810686166554974, 0.04636038757292936, 0.04809695408872083, 0.04707321576655934, 0.04706761650373703, 0.048810773927350976, 0.04840533161464197, 0.047094149804909835, 0.045470281325987655, 0.04230194461342957, 0.041359611137599804, 0.043056218062010744, 0.041271484867300444, 0.046301688142634724, 0.04759489106052545, 0.04466810209250278, 0.04404609196180724, 0.04673197095881622, 0.044727800835898175, 0.05077610106156217, 0.05125762797507319, 0.045935416907246375, 0.04370548728589216, 0.043114325710940156, 0.04598968711062229, 0.047986350587714124, 0.04490135736985996, 0.047769969596535934, 0.04968747192234251, 0.04740765994204652, 0.04399974079604498, 0.044065997308604554, 0.040269053288832964, 0.04118791833682217, 0.039964015212371945, 0.043444497280207044, 0.047326487416182644, 0.04473860791411215, 0.04386694154450048, 0.041537232314026916, 0.03997830429162584, 0.039120541700435656, 0.041581151936315563, 0.04628740285788258, 0.047870418269909634, 0.04662154968351817, 0.05002896415585355, 0.04987150472890811, 0.04353202047827673, 0.043696275087509985, 0.04392725797938837, 0.04166027880513969, 0.041483596834816776, 0.045021385292791635, 0.04427139201545601, 0.03911876423397845, 0.03918975284429804, 0.040871536108388624, 0.04255490957923923, 0.04410847737839581, 0.04365755111507619, 0.04434324899898308, 0.043760502574852596, 0.04057433497304218, 0.043540525620750274, 0.045716653046709535, 0.04279001072160254, 0.04309528912055374, 0.0498528050112827, 0.046185296176492675, 0.044122197472394216, 0.04208540126775135, 0.04048122509041643, 0.04251179186239604, 0.046360950532881896, 0.04553776498283139, 0.04271723000127095, 0.04214051069373384, 0.04027143534511742, 0.039556138489226705, 0.0394454767577124, 0.03915745308048971, 0.0426604031784745, 0.042422507532394564, 0.04125278523542767, 0.04074152532316104, 0.041894423444072335, 0.04037606399142927, 0.042311893718620465, 0.04290735486069182, 0.046516723700828634, 0.04790502752344351, 0.04919419602601001, 0.05001314822304484, 0.045424186087659205, 0.04029647269131899, 0.04127462412561173, 0.043548064056534753, 0.040808158667153566, 0.043283981292258696, 0.0464351497056676, 0.04818190251415473, 0.049183400932093654, 0.048399220361180106, 0.04420213851078221, 0.04200595490637232, 0.03982505339585597, 0.04178307840711175, 0.044481540408521744, 0.04520050300192068, 0.04282290867811043, 0.04311880424078224, 0.043315776710366856, 0.040264894608854636, 0.037462886870982474, 0.038757597550494155, 0.0412408704365306, 0.04535179324365991, 0.04651100197075401, 0.049101883946623036, 0.046286993356481856, 0.041798108753039694, 0.0425566794873976, 0.04169336839195659, 0.03977126848163422, 0.04204580601534974, 0.043254185547443355, 0.044778931491746145, 0.04346620053196143, 0.046465892522302564, 0.04006725050440246, 0.039324297098618045, 0.03667246919554188, 0.0403557210787201, 0.043825950075946356, 0.0422603088658769, 0.044263827977230434, 0.04428082489573264, 0.04007741767110468, 0.04646821240560842, 0.04861327790634655, 0.04180092676958734, 0.04173035337541241, 0.043271644332236445, 0.03962366788978274, 0.03994462052197399, 0.043506825274811055, 0.04346198901082906, 0.040069535811447766, 0.04561720126686115, 0.04574076939583449, 0.04486419532169502, 0.04846571817577178, 0.04363466501242077, 0.03984455704792209, 0.039082349301494984, 0.037851259206999756, 0.04522511481330434, 0.05128732194609073, 0.05051850833022727, 0.04569586487957692, 0.04224457362090452, 0.03752153947177082, 0.040881507368870584, 0.04123285085935292, 0.04274551603208298, 0.04276409424482033, 0.04091018414739487, 0.04231569406471316, 0.04755138819432719, 0.04935750437017195, 0.05017080509792138, 0.04801575844214858, 0.04631413253043159, 0.04450693552005502, 0.04409014954821563, 0.04652019721821898, 0.04379161410269631, 0.04157737100102083, 0.044741352367445815, 0.04600775206181546, 0.04716394296707482, 0.04641384861812514, 0.04849827792669965, 0.0478487287138147, 0.04220017321454819, 0.040834003198488567, 0.03669124831424792, 0.0395635178514907, 0.04390370738594051, 0.043966603327322884, 0.042073404323273014, 0.04241675826602162, 0.04468830730576528, 0.04857588956877838, 0.05032217001306757, 0.049438664972429734, 0.041893852088812134, 0.04164206124160698, 0.04715555912360174, 0.04760571781241883, 0.04227797806939574, 0.046793582060116454, 0.04869569838080588, 0.049064610435668725, 0.044097485690700176, 0.04711472626671868, 0.0458748920588632, 0.047560856432347066, 0.04619807952679293, 0.04921976937205322, 0.0529110537758894, 0.04676898017160896, 0.03922334691038973, 0.04433506968622178, 0.04280558422000816, 0.04641031000443531, 0.04714838523779313, 0.0445472715801689, 0.04796292047299996, 0.04526961202175698, 0.043212448674053154, 0.04218481352743408, 0.04305622655668058, 0.04569258078022905, 0.043098888977329636, 0.03810730185094095, 0.04518896225414396]],
            'im_sum': [False, 2537.922123380079],
            'regn_sum': [False, 268.8044554442167],
            'npts_real': [True, 6400000],
            'profile': [False, 2.806226167884636],
            'fit': [False, [3.1062863101718547, 6.1773240682602095, 5.744644158900882]],
            'fit_loc_chan': [True, 489],
            'fit_loc_freq': [1e-10, 354.50497926661194],
            'fit_pix': [False, [45.65720045662228, 41.04610803971941]]}

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
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'mask_pix': [False, 9849],
            'mask_regns': [True, 1],
            'npts_real': [True, 6400000]}

        report3 = th.check_dict_vals(exp_mask_stats, mask_stats_dict, '.mask', epsilon=self.epsilon)

        # .pb report
        pb_stats_dict = self.image_stats(img+'.pb', fit_region = \
            'ellipse[[93.24066271deg, 22.57227012deg], [22.4071arcsec, 21.7699arcsec], 90.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_pb_stats
        exp_pb_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, 0.20036059617996216],
            'im_rms': [False, 0.576840993910233],
            'npts_0.2': [False, [3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233, 3233]],
            'npts_0.5': [False, [1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549, 1549]],
            'npts_real': [True, 6400000],
            'fit': [False, [1.046847676114786, 28.075049566294457, 28.075049566280292]],
            'fit_loc_chan': [True, 500],
            'fit_loc_freq': [1e-10, 354.50632237101223],
            'fit_pix': [False, [40.0, 40.0]]}

        report4 = th.check_dict_vals(exp_pb_stats, pb_stats_dict, '.pb', epsilon=self.epsilon)

        # .psf report
        psf_stats_dict = self.image_stats(img+'.psf', fit_region = \
            'ellipse[[93.24070675deg, 22.57229969deg], [4.6173arcsec, 3.4331arcsec], 90.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_psf_stats
        exp_psf_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.0],
            'max_val_pos': [True, [40, 40, 0, 0]],
            'min_val': [False, -0.16458643972873688],
            'min_val_pos': [True, [34, 36, 0, 949]],
            'im_rms': [False, 0.086535984124565],
            'im_sum': [False, 2855.491490987678],
            'npts_real': [True, 6400000],
            'fit_0': [False, [1.0955776643770092, 4.02083883935996, 2.9704405181863174]],
            'fit_loc_chan_0': [True, 1],
            'fit_loc_freq_0': [1e-10, 354.44539427140006],
            'fit_pix_0': [False, [40.000070982740574, 39.997503447253926]],
            'fit_1': [False, [1.0955889220173731, 4.020278276279631, 2.970194228997147]],
            'fit_loc_chan_1': [True, 500],
            'fit_loc_freq_1': [1e-10, 354.50632237101223],
            'fit_pix_1': [False, [40.00007124298882, 39.997504801396616]],
            'fit_2': [False, [1.0956310542682013, 4.019922686595429, 2.96902152077006]],
            'fit_loc_chan_2': [True, 999],
            'fit_loc_freq_2': [1e-10, 354.5672504706244],
            'fit_pix_2': [False, [40.000072198599625, 39.99750719898725]]}


        report5 = th.check_dict_vals(exp_psf_stats, psf_stats_dict, '.psf', epsilon=self.epsilon)

        # .residual report
        resid_stats_dict = self.image_stats(img+'.residual', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]')

        # test_standard_cube_eph_briggsbwtaper.exp_resid_stats
        exp_resid_stats = {'npts': [True, 6400000],
            'npts_unmasked': [False, 3233000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 0.32331418991088867],
            'max_val_pos': [True, [20, 53, 0, 497]],
            'min_val': [False, -0.3310726284980774],
            'min_val_pos': [True, [32, 47, 0, 490]],
            'im_rms': [False, 0.046640536710825437],
            'im_sum': [False, 331.0392044886937],
            'regn_sum': [False, 51.19282157626003],
            'npts_real': [True, 6400000]}

        report6 = th.check_dict_vals(exp_resid_stats, resid_stats_dict, \
            '.residual', epsilon=self.epsilon)

        # .model report
        model_stats_dict = self.image_stats(img+'.model', fit_region = \
            'ellipse[[93.23955772deg, 22.57243258deg], [6.2806arcsec, 5.4290arcsec], 0.00000000deg]', masks=mask_stats_dict['mask'])

        # test_standard_cube_eph_briggsbwtaper.exp_model_stats
        exp_model_stats = {'npts': [True, 6400000],
            'npts_unmasked': [True, 6400000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1.4160760641098022],
            'max_val_pos': [True, [46, 41, 0, 489]],
            'min_val': [False, -0.070622518658638],
            'min_val_pos': [True, [58, 39, 0, 491]],
            'im_rms': [False, 0.0017559473893058901],
            'im_sum': [False, 55.33133029192686],
            'regn_sum': [False, 56.097182389348745],
            'mask_non0': [True, 0],
            'npts_real': [True, 6400000]}


        report7 = th.check_dict_vals(exp_model_stats, model_stats_dict, \
            '.model', epsilon=self.epsilon)

        # .sumwt report
        sumwt_stats_dict = self.image_stats(img+'.sumwt')

        # test_standard_cube_eph_briggsbwtaper.exp_sumwt_stats
        exp_sumwt_stats = {'npts': [True, 1000],
            'npts_unmasked': [True, 1000.0],
            'freq_bin': [1e-10, 122100.40002441406],
            'start': [True, 354445000000.0],
            'end': [True, 354567300000.0],
            'start_delta': [False, 354445000000.0],
            'end_delta': [False, 354567300000.0],
            'nchan': [True, 1000],
            'max_val': [False, 1012.2438354492188],
            'max_val_pos': [True, [0, 0, 0, 70]],
            'min_val': [False, 1012.075927734375],
            'min_val_pos': [True, [0, 0, 0, 475]],
            'im_rms': [False, 1012.1570928591075],
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

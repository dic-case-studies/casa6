# test_stk_vlass_pipeline_imaging.py
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
# [https://open-jira.nrao.edu/browse/CAS-12427]
#
# Based on the requirements listed here:
# 
#
##########################################################################


'''
Datasets :
---------------
VLASS 1.1 : 
    - New/Old ACUs : Pointing offset between two sets of antennas is about 45 arcsec (upto 1 arcmin)
    - Pointing nods between integrations : 1.5 arcmin between integrations.
VLASS 1.2 : 
     - Pointing nods between integrations : 1.5 arcmin between integrations

For both observation epochs : 3C286  :    
               - 1  pointing, 2 integrations
               - 5 pointings, 5x2=10 integrations

TODO : Prune the pointing table.
             Get VLASS 1.2 dataset of 3c286


Algorithm Options : 
-------------------------

New parameter :  pointingsigdev : [ bin width,  mean shift width ]

For VLASS data, [300, 300] will treat all antennas together, and no time-dependent pointing.
        - This will be similar to do-pointing=False except that the time-step will be every 300 sec. 

For VLASS 1.1, a default of [ 30 arcsec,  30 arcsec ] is good.  Factor of 2 slower than [300,300]
For VLASS 1.2,  we can do [ 300arcsec, 30 arcsec ] to treat all antennas together, but to keep the time-dependence corrections.

If you sub-select antennas for new or old ACUs, then pointingsigdev applies to the selected antennas only. So, if selecting only new ACUs from VLASS 1.1 data, can set [ 300.0, 30.0 ]. 


'''

fieldnames = {'18952':'1447214-161911',
              '18953':'1447338-161911',
              '18954':'1447463-161911', 
              '18955':'1447587-161911', 
              '18956':'1448111-161911', 
              '18957':'1448236-161911', 
              '18958':'1448360-161911', 
              '18959':'1448484-161911', 
              '18960':'1449009-161911'}

fieldcenters = {'18952':'J2000 14:47:21.451021 -16.19.12.00000',
                '18953':'J2000 14:47:33.886250 -16.19.12.00000',
                '18954':'J2000 14:47:46.321479 -16.19.12.00000', 
                '18955':'J2000 14:47:58.756708 -16.19.12.00000', 
                '18956':'J2000 14:48:11.191937 -16.19.12.00000', 
                '18957':'J2000 14:48:23.627167 -16.19.12.00000', 
                '18958':'J2000 14:48:36.062396 -16.19.12.00000', 
                '18959':'J2000 14:48:48.497625 -16.19.12.00000', 
                '18960':'J2000 14:49:00.932854 -16.19.12.00000', 
                '18961':'J2000 14:49:13.368083 -16.19.12.00000'}

alltimes = {'18952': {'one':'<13:34:14.7',  'two':'>13:34:14.7'} , 
            '18953': {'one':'<13:34:15.6',  'two':'>13:34:15.6'} , 
            '18954': {'one':'<13:34:16.5',  'two':'>13:34:16.5'} , 
            '18955': {'one':'<13:34:17.4',  'two':'>13:34:17.4'} , 
            '18956': {'one':'<13:34:18.3',  'two':'>13:34:18.3'} , 
            '18957': {'one':'<13:34:19.2',  'two':'>13:34:19.2'} , 
            '18958': {'one':'<13:34:20.1',  'two':'>13:34:20.1'} ,
            '18959': {'one':'<13:34:21.0',  'two':'>13:34:21.0'} ,
            '18960': {'one':'<13:34:21.9',  'two':'>13:34:21.9'} }

pixlocs = {'18953':[3025,4757,0,0],
               '18954':[3621,4758,0,0],
               '18955':[4218,4758,0,0],
               '18956':[4815,4758,0,0],
               '18957':[5411,4758,0,0],
               '18958':[6008,4758,0,0],
               '18959':[6604,4758,0,0],
               '18960':[7201,4757,0,0] }

##########################################################################
### NOTES

## pointingsigdev parameter needs to now show a warning for the default of usepointing=False !!! 

##########################################################################

# Imports #
import os
import subprocess
import glob
import time
import sys
import shutil
import unittest
import inspect
import pprint
import numpy
import pylab as pl

CASA6 = False
 
try:
    from casatools import ctsys
    from casatools import image as iatool
    from casatasks import tclean, casalog
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    CASA6 = True

    # CASA6 doesn't need default
    def default(atask):
        pass

    ctsys_resolve = ctsys.resolve
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from parallel.parallel_task_helper import ParallelTaskHelper

    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'data')
        return os.path.join(dataPath,apath)

import casatestutils as ctu
from casatestutils.imagehelpers import TestHelpers
imh = TestHelpers()


# Path to data
#data_path = '/lustre/rurvashi/StakeholderTests/VLASS/DATA_FOR_CAS12427/'
data_path = ctsys_resolve('stakeholders/vlass/')

##############################################################
############# Test settings #######################
## CFCAche settings. Default oversampling=20. Set to 2 to make the CFCache smaller, at the cost of accuracy
os.environ['ATerm_CONVSIZE'] = '512'
os.environ['ATerm_OVERSAMPLING'] = '2'
cfcache_path  = data_path+'CFCACHE_OS2'

## Number of iterations. This will test major and minor cycles and get better numerical results. (About 3 major cycles each)
#niter=100
## ( Set to niter=0 for a faster run and approximate results. It will still test for code changes. 
niter=0

## A flag to run the tclean commands (True) as well as make the reports.
## False : Just remake the reports from existing outputs (False, useful for debugging the test scripts).
do_row=True


#(1)# Image all the fields in the MS to generate the full 'rowplot' plots and show full trends. 
#fields = ['18954','18955','18956','18957','18958','18959']
#valcheck = 3  # index of field id 18957. Test for pass/fail on numbers from one field only.

#(2)# For a faster test run, do only one field.  No plots will be made, but the same pass/fail tests/checks apply.
fields = ['18957']
valcheck = 0  # index of field id 18957

#################################################################

## Base Test class with Utility functions
class TestTcleanBase(unittest.TestCase):

    def setUp(self):
        self._myia = iatool()
        self.epsilon = 0.1 # sets epsilon as a percentage (10%)
#        self.msfile = ""
        self.img = "tst"

        self.newACU = 'ea01, ea02, ea08, ea14, ea17, ea21' 
        ##antenna = '0,1,6,12,15,19'  ## new ACU ants for this test dataset.
        ##antenna = '1,2,7,8,14,16,17,21,27,28' # Full list of ants with new ACUs (Summer 2018)
        self.oldACU = '!'+self.newACU

        # To use subdir in the output image names in some tests (CAS-10937)
        self.img_subdir = 'refimager_tst_subdir'
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

    def tearDown(self):
        ctu.generate_weblog("tclean_VLASS_pipeline",test_dict)
        self._myia.done()
        """ don't delete it all """
#        self.delData()

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self,msname=""):
#        if msname != "":
#             self.msfile=msname
        self.dataset1 = 'vlass1.1_rowtest_6fields.ms' #'PcorrVis_field18954.ms'
        self.fieldname1 = fieldnames['18954']
        self.timerange1two = alltimes['18954']['two']

        self.dataset2 =  'vlass1.1_rowtest_6fields.ms' #'PcorrVis_field18957.ms'
        self.fieldname2 = fieldnames['18957']


    def delData(self,msname=""):
        return
#        if msname != "":
#             self.msfile=msname
#        if (os.path.exists(self.msfile)):
#             os.popen('rm -rf ' + self.msfile)
#        os.popen('rm -rf ' + self.img_subdir)
#        os.popen('rm -rf ' + self.img+'*')

    def getNameDoc(self):
#        print inspect.stack()[1]
        testname=inspect.stack()[1][3]
        print("Test name  : " + testname)
#        tname = inspect.getframeinfo(inspect.currentframe()).function
        doc = eval('self.'+testname + '.__doc__')
        print("Doc : " +  doc)
        return testname, doc


    def image_stats(self,img,suf):
        self._myia.open(img+suf)
        rms = self._myia.statistics()['rms'][0]  # reads RMS value for image
        int_flux = self._myia.statistics()['sum'][0]  # reads integrated 
                                                      # flux value for image
        major = self._myia.commonbeam()['major']['value']
        minor = self._myia.commonbeam()['minor']['value']
        pa = self._myia.commonbeam()['pa']['value']
        self._myia.close()

        return rms, int_flux, major, minor, pa

    def image_metrics(self, img, loc):
        # Read intensity, alpha and PB at source location
        intval = imh.get_pix(img+'.image.tt0', loc)
        spxval = imh.get_pix(img+'.alpha', loc)
        pbval = imh.get_pix(img+'.pb.tt0', loc)

        # PB-correction
        if pbval>0.0:
            pbcor_intval = intval/pbval
        else:
            pbcor_intval=0.0 # null value

        return pbcor_intval, spxval, pbval

    def image_list(self,img,mode):

        imlist0 = [img+'.alpha', img+'.alpha.error',\
                   img+'.psf.tt0', img+'.psf.tt1', img+'.psf.tt2', \
                   img+'.residual.tt0', img+'.residual.tt1', \
                   img+'.pb.tt0', \
                   img+'.sumwt.tt0', img+'.sumwt.tt1', img+'.sumwt.tt2']

        imlist1 = [img+'.model.tt0', img+'.model.tt1', \
                   img+'.image.tt0', img+'.image.tt1']

        if mode == 'niterg0':
            return imlist0 + imlist1
        else:
            return imlist0

    def make_png(self, imname, out):
        self._myia.open(imname)
        pix = self._myia.getchunk()[:,:,0,0]
        self._myia.close()
        
        pl.imsave(fname=out, arr=numpy.transpose(pix), cmap='jet',origin='lower')

    def make_rowplot(self, fields, intvals, spxvals, pbvals, fname):
        pl.ioff()
        pl.clf()
        pl.figure(figsize=(8,6))
        fids=[]
        for ff in fields:
            fids.append(eval(ff))
        pl.plot( fids, intvals , 'o-', label='Sky Intensity' )
        pl.plot( fids, spxvals,'o-', label='Sky Alpha' )
        pl.plot( fids, pbvals,'o-', label='PB gain' )
        pl.hlines( y= 1.38 ,xmin=fids[0],xmax=fids[len(fids)-1],color='m')
        pl.hlines( y=-0.46 ,xmin=fids[0],xmax=fids[len(fids)-1],color='c')
        pl.axis([fids[0],fids[len(fids)-1],-2.5,2.5])
        pl.vlines(18956.7,-2.5,2.5,color='k',linestyle='--')
        pl.legend(loc='lower left')

        pl.savefig(fname)

##############################################
test_dict = {}
####    Tests     ####
class Test_vlass_1p1_row(TestTcleanBase):

    @ctu.stats_dict(test_dict)
    def test_vlass_1p1_row_pcorr0(self):
        '''
        ---------------------------------------------------------------------------------------------------------------------------------------------
        VLASS 1.1 : No Pointing Corrections. MS field phasecenter is treated as the pointing center for all antennas. 
                           Image each OTF field separately (to make a scan plot)
                           Use mtmfs with conjbeams True and test pb-corrected intensity and spectral index.

        Test that the intensity and spectral index are wrong (pick known wrong values). 
        The plot should show a slope in intensity and spectral index arising from the average pointing offset
        in the majority of the antennas (ie. 20 old ACU ants) between the true pointing and the recorded phase center.

        This is the state of CASA 5.6 for usepointing=False
        ---------------------------------------------------------------------------------------------------------------------------------------------
        '''
        testname, testdoc = self.getNameDoc()

        intvals=[]
        spxvals=[]
        pbvals=[]

        #cfcache_oversampling1='CFCACHE_OS1'

        for field in fields:
            file_name = 'im_'+testname+'_'+field
            img = os.getcwd()+'/'+file_name
            #msname = data_path + 'PcorrVis_field'+field+'.ms'
            msname = data_path+'vlass1.1_rowtest_6fields.ms'
            
            if do_row==True:
                os.system('rm -rf '+img+'.*')
                tclean(vis=msname, \
                       imagename=img, spw='2~17',antenna='', \
                       timerange='', \
                       field=fieldnames[field] ,uvrange='<12.0km',\
                       cell='0.3arcsec', imsize=10000, phasecenter="J2000 14:48:31.561 -16.18.56.120", \
                       specmode='mfs',nterms=2, reffreq='3.0GHz',\
                       deconvolver='mtmfs',weighting='uniform',\
                       gridder='awproject', cfcache=cfcache_path,wbawp=True,\
                       pblimit=0.02,psterm=False, conjbeams=True,\
                       wprojplanes=64,niter=niter,datacolumn='data',\
                       mosweight=False, usepointing=False, 
                       pointingoffsetsigdev=300.0,\
                       parallel=self.parallel)

            intval, spxval, pbval = self.image_metrics(img=img, loc=[5792,4705,0,0])

            intvals.append(intval)
            spxvals.append(spxval)
            pbvals.append(pbval)

        
        if niter==0:   ## Approximate - initial residual image
            truth_int = 1.19
            truth_spx = -0.78
            truth_pb = 0.92
        else:   ## A bit more accurate after a few major cycles
            truth_int = 1.19
            truth_spx = -0.78
            truth_pb = 0.92
        
        out, report1 = imh.check_val(  intvals[valcheck], truth_int, valname = 'PBCor Intensity', epsilon = truth_int*self.epsilon)
        out, report2 = imh.check_val(  spxvals[valcheck], truth_spx, valname = 'Spectral Index', epsilon = self.epsilon)
        out, report3 = imh.check_val(  pbvals[valcheck], truth_pb, valname = 'PB gain', epsilon = truth_pb*self.epsilon)
        report = report1 + report2 + report3
        report = report + 'Test values for field '+ fields[valcheck] + ' with name ' + fieldnames[fields[valcheck]] + '\n'

        test_dict[testname]['Test Description'] = testdoc
        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['Report'] = report

        if len(fields)>1:
            self.make_rowplot(fields, intvals, spxvals, pbvals, testname+'_rowplot.png')
            test_dict[testname]['Figures'] = [testname+'_rowplot.png']

        ctu.add_to_dict(self, output = test_dict, dataset = msname)
        
        self.assertTrue(imh.check_final(pstr = report), msg = report)


    @ctu.stats_dict(test_dict)
    def test_vlass_1p1_row_pcorr1_oneint(self):
        '''
        ---------------------------------------------------------------------------------------------------------------------------------------------
        VLASS 1.1 :  Corrections for antenna-groups (new/old ACUs). psd = [30,300]
                           Image each OTF field separately (to make a scan plot)
                           Use mtmfs with conjbeams True and test pb-corrected intensity and spectral index.
                           Image only one integration (first one per scan).

        Test that the intensity and spectral index are correct.  
        This tests the antenna clustering using just one timestep of data (i.e. no time-dependent errors)

        [ If the second timestep is included (or chosen) the slope re-appears
           because with psf=[x,300], the pointing offset is read only from the first timestep in a VLASS scan. ]
        ---------------------------------------------------------------------------------------------------------------------------------------------
        '''
        testname, testdoc = self.getNameDoc()

#        fields = ['18954','18955','18956','18957','18958','18959']

        intvals=[]
        spxvals=[]
        pbvals=[]

        #cfcache_oversampling1='CFCACHE_OS1'

        for field in fields:
            file_name = 'im_'+testname+'_'+field
            img = os.getcwd()+'/'+file_name
            #msname = data_path + 'PcorrVis_field'+field+'.ms'
            msname = data_path+'vlass1.1_rowtest_6fields.ms'
            
            if do_row==True:
                os.system('rm -rf '+img+'.*')
                tclean(vis=msname, \
                       imagename=img, spw='2~17',antenna='', \
                       timerange=alltimes[field]['one'], \
                       field=fieldnames[field] ,uvrange='<12.0km',\
                       cell='0.3arcsec', imsize=10000, phasecenter="J2000 14:48:31.561 -16.18.56.120", \
                       specmode='mfs',nterms=2, reffreq='3.0GHz',\
                       deconvolver='mtmfs',weighting='uniform',\
                       gridder='awproject', cfcache=cfcache_path,wbawp=True,\
                       pblimit=0.02,psterm=False, conjbeams=True,\
                       wprojplanes=64,niter=niter,datacolumn='data',\
                       mosweight=False, usepointing=True, 
                       pointingoffsetsigdev=[30.0, 300.0],\
                       parallel=self.parallel)

            intval, spxval, pbval = self.image_metrics(img=img, loc=[5792,4705,0,0])

            intvals.append(intval)
            spxvals.append(spxval)
            pbvals.append(pbval)

        if niter==0:   ## Approximate - initial residual image
            truth_int = 1.38
            truth_spx = -0.6
            truth_pb = 0.93
        else:   ## A bit more accurate after a few major cycles
            truth_int = 1.38
            truth_spx = -0.6
            truth_pb = 0.93
        
        out, report1 = imh.check_val(  intvals[valcheck], truth_int, valname = 'PBCor Intensity', epsilon = truth_int*self.epsilon)
        out, report2 = imh.check_val(  spxvals[valcheck], truth_spx, valname = 'Spectral Index', epsilon = self.epsilon)
        out, report3 = imh.check_val(  pbvals[valcheck], truth_pb, valname = 'PB gain', epsilon = truth_pb*self.epsilon)

        report = report1 + report2 + report3
        report = report + 'Test values for field '+ fields[valcheck] + ' with name ' + fieldnames[fields[valcheck]] + '\n'

        test_dict[testname]['Test Description'] = testdoc
        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['Report'] = report

        if len(fields)>1:
            self.make_rowplot(fields, intvals, spxvals, pbvals, testname+'_rowplot.png')
            test_dict[testname]['Figures'] = [testname+'_rowplot.png']

        ctu.add_to_dict(self, output = test_dict, dataset = msname)
        
        self.assertTrue(imh.check_final(pstr = report), msg = report)



    @ctu.stats_dict(test_dict)
    def test_vlass_1p1_row_pcorr1_twoint(self):
        '''
        ---------------------------------------------------------------------------------------------------------------------------------------------
       VLASS 1.1 :  Corrections for antenna-groups (new/old ACUs). psd = [30,300]
                            Image each OTF field separately (to make a scan plot)
                           Use mtmfs with conjbeams True and test pb-corrected intensity and spectral index.
                           Image both integrations in each scan.

        Test that the intensity and spectral index are wrong. Antenna clustering is applied, but time dependence it not.
        
        If the second timestep is included (or chosen) the slope re-appears
        because with psf=[x,300], the pointing offset is read only from the first timestep in a VLASS scan.

        This is a redundant test.  It need not be included in stakeholder verification tests, but is useful to illustrate
        an effect. 
        ---------------------------------------------------------------------------------------------------------------------------------------------
        '''
        testname, testdoc = self.getNameDoc()

#        fields = ['18954','18955','18956','18957','18958','18959']

        intvals=[]
        spxvals=[]
        pbvals=[]

        #cfcache_oversampling1='CFCACHE_OS1'

        for field in fields:
            file_name = 'im_'+testname+'_'+field
            img = os.getcwd()+'/'+file_name
            #msname = data_path + 'PcorrVis_field'+field+'.ms'
            msname = data_path+'vlass1.1_rowtest_6fields.ms'
            
            if do_row==True:
                os.system('rm -rf '+img+'.*')
                tclean(vis=msname, \
                       imagename=img, spw='2~17',antenna='', \
                       timerange='', \
                       field=fieldnames[field] ,uvrange='<12.0km',\
                       cell='0.3arcsec', imsize=10000, phasecenter="J2000 14:48:31.561 -16.18.56.120", \
                       specmode='mfs',nterms=2, reffreq='3.0GHz',\
                       deconvolver='mtmfs',weighting='uniform',\
                       gridder='awproject', cfcache=cfcache_path,wbawp=True,\
                       pblimit=0.02,psterm=False, conjbeams=True,\
                       wprojplanes=64,niter=niter,datacolumn='data',\
                       mosweight=False, usepointing=True, 
                       pointingoffsetsigdev=[30.0, 300.0],\
                       parallel=self.parallel)

            intval, spxval, pbval = self.image_metrics(img=img, loc=[5792,4705,0,0])

            intvals.append(intval)
            spxvals.append(spxval)
            pbvals.append(pbval)

        if niter==0:   ## Approximate - initial residual image
            truth_int = 1.17
            truth_spx = -0.77
            truth_pb = 0.93
        else:   ## A bit more accurate after a few major cycles
            truth_int = 1.17
            truth_spx = -0.77
            truth_pb = 0.93

        
        out, report1 = imh.check_val(  intvals[valcheck], truth_int, valname = 'PBCor Intensity', epsilon = truth_int*self.epsilon)
        out, report2 = imh.check_val(  spxvals[valcheck], truth_spx, valname = 'Spectral Index', epsilon = self.epsilon)
        out, report3 = imh.check_val(  pbvals[valcheck], truth_pb, valname = 'PB gain', epsilon = truth_pb*self.epsilon)
        report = report1 + report2 + report3

        report = report + 'Test values for field '+ fields[valcheck] + ' with name ' + fieldnames[fields[valcheck]] + '\n'

        test_dict[testname]['Test Description'] = testdoc
        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['Report'] = report

        if len(fields)>1:
            self.make_rowplot(fields, intvals, spxvals, pbvals, testname+'_rowplot.png')
            test_dict[testname]['Figures'] = [testname+'_rowplot.png']

        ctu.add_to_dict(self, output = test_dict, dataset = msname)
        
        self.assertTrue(imh.check_final(pstr = report), msg = report)


    @ctu.stats_dict(test_dict)
    def test_vlass_1p1_row_pcorr1and2_twoint(self):
        '''
        ---------------------------------------------------------------------------------------------------------------------------------------------
        VLASS 1.1 :  Corrections for antenna-groups (new/old ACUs). and time-dependent nodding :psd = [30,30]
                           Image each OTF field separately (to make a scan plot)
                           Use mtmfs with conjbeams True and test pb-corrected intensity and spectral index.
                           Use both integrations per scan        

        Test that the intensity and spectral index are correct.  

        This is the setting most optimal for VLASS.  Perhaps the stakeholder test needs just this one run. 
        ---------------------------------------------------------------------------------------------------------------------------------------------
        '''
        testname, testdoc = self.getNameDoc()

#        fields = ['18954','18955','18956','18957','18958','18959']

        intvals=[]
        spxvals=[]
        pbvals=[]

        #cfcache_oversampling1='CFCACHE_OS1'

        for field in fields:
            file_name = 'im_'+testname+'_'+field
            img = os.getcwd()+'/'+file_name
            #msname = data_path + 'PcorrVis_field'+field+'.ms'
            msname = data_path+'vlass1.1_rowtest_6fields.ms'
            
            if do_row==True:
                os.system('rm -rf '+img+'.*')
                tclean(vis=msname, \
                       imagename=img, spw='2~17',antenna='', \
                       timerange='', \
                       field=fieldnames[field] ,uvrange='<12.0km',\
                       cell='0.3arcsec', imsize=10000, phasecenter="J2000 14:48:31.561 -16.18.56.120", \
                       specmode='mfs',nterms=2, reffreq='3.0GHz',\
                       deconvolver='mtmfs',weighting='uniform',\
                       gridder='awproject', cfcache=cfcache_path,wbawp=True,\
                       pblimit=0.02,psterm=False, conjbeams=True,\
                       wprojplanes=64,niter=niter,datacolumn='data',\
                       mosweight=False, usepointing=True, 
                       pointingoffsetsigdev=[30.0, 30.0],\
                       parallel=self.parallel)

            intval, spxval, pbval = self.image_metrics(img=img, loc=[5792,4705,0,0])

            intvals.append(intval)
            spxvals.append(spxval)
            pbvals.append(pbval)

        if niter==0:   ## Approximate - initial residual image
            truth_int = 1.32
            truth_spx = -0.66
            truth_pb = 0.84
        else:   ## A bit more accurate after a few major cycles
            truth_int = 1.32
            truth_spx = -0.58
            truth_pb = 0.84

        
        out, report1 = imh.check_val(  intvals[valcheck], truth_int, valname = 'PBCor Intensity', epsilon = truth_int*self.epsilon)
        out, report2 = imh.check_val(  spxvals[valcheck], truth_spx, valname = 'Spectral Index', epsilon = self.epsilon)
        out, report3 = imh.check_val(  pbvals[valcheck], truth_pb, valname = 'PB gain', epsilon = truth_pb*self.epsilon)

        report = report1 + report2 + report3
        report = report + 'Test values for field '+ fields[valcheck] + ' with name ' + fieldnames[fields[valcheck]] + '\n'

        test_dict[testname]['Test Description'] = testdoc
        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['Report'] = report

        if len(fields)>1:
            self.make_rowplot(fields, intvals, spxvals, pbvals, testname+'_rowplot.png')
            test_dict[testname]['Figures'] = [testname+'_rowplot.png']

        ctu.add_to_dict(self, output = test_dict, dataset = msname)
        
        self.assertTrue(imh.check_final(pstr = report), msg = report)

    @ctu.stats_dict(test_dict)
    def test_vlass_1p1_jointmosaic_pcorr1and2(self):
        '''
        ---------------------------------------------------------------------------------------------------------------------------------------------
        VLASS 1.1 : With both types of pointing correction ( antenna clusters and time-dependent nodding )
                                        : posigdev=[30,30] 
                           Use two scans/fields and all antennas together to test the field mosaic pointing offset.
                           Use mtmfs with conjbeams True and test pb-corrected intensity and spectral index.

        Test that the intensity and spectral index are correct. 

        This is the setting most optimal for VLASS.
        ---------------------------------------------------------------------------------------------------------------------------------------------
        '''
        testname, testdoc = self.getNameDoc()

        file_name = 'im_'+testname
        self.prepData()

        img = os.getcwd()+'/'+file_name
        msname = data_path+'vlass1.1_rowtest_6fields.ms'

        if do_row==True:
            os.system('rm -rf '+img+'.*')
            tclean(vis=msname, \
                   imagename=img, spw='2~17',antenna='', \
                   timerange='', \
                   field='1~4' ,uvrange='<12.0km',\
                   cell='0.3arcsec', imsize=10000, phasecenter="J2000 14:48:31.561 -16.18.56.120", \
                   specmode='mfs',nterms=2, reffreq='3.0GHz',\
                   deconvolver='mtmfs',weighting='uniform',\
                   gridder='awproject', cfcache=cfcache_path,wbawp=True,\
                   pblimit=0.02,psterm=False, conjbeams=True,\
                   wprojplanes=64,niter=niter,datacolumn='data',\
                   mosweight=False, usepointing=True, 
                   pointingoffsetsigdev=[30.0, 30.0],\
                   parallel=self.parallel)
            

        intval, spxval, pbval = self.image_metrics(img=img, loc=[5792,4705,0,0])

        if niter==0:   ## Approximate - initial residual image
            truth_int = 1.3
            truth_spx = -0.67
            truth_pb = 0.94
        else:   ## A bit more accurate after a few major cycles
            truth_int = 1.3
            truth_spx = -0.5
            truth_pb = 0.94

       
        #report0 = imh.check_imexist( self.image_list( img, 'niter0') )
        out, report1 = imh.check_val(  intval, truth_int, valname = 'allants : PBCor Intensity', epsilon = truth_int*self.epsilon)
        out, report2 = imh.check_val(  spxval, truth_spx, valname = 'allants : Spectral Index', epsilon = self.epsilon)
        out, report3 = imh.check_val(  pbval, truth_pb, valname = 'allants : PB gain', epsilon = truth_pb*self.epsilon)


        report = report1 + report2 +report3

        test_dict[testname]['Test Description'] = testdoc

        test_dict[testname]['self.parallel'] = self.parallel
        test_dict[testname]['Report'] = report

        self.make_png(img+'.pb.tt0', out = img+'.pb.tt0.png')
        test_dict[testname]['Figures'] = [img+'.pb.tt0.png']

        ctu.add_to_dict(self, output = test_dict, dataset = self.dataset1)

        self.assertTrue(imh.check_final(pstr = report), msg = report)

def suite():
     return [Test_vlass_1p1_row, Test_standard] #[Test_tclean_ALMA]

# Main #
if __name__ == '__main__':
    unittest.main()



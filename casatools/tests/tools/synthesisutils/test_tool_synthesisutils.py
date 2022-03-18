##########################################################################
# test_tool_synthesisutils.py
#
# Copyright (C) 2021
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
# [https://open-jira.nrao.edu/browse/CAS-13590]
#
# Based on the requirements listed in CASADocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatools.synthesisutils.html
#
# Test case: requirement
# methods tested: getoptimumsize, fitpsfbeam, advisechansel 
#
# test_fitpsfbeam
# test single term psf fitting 
# test a wrong nterms (>1) for single term psf fitting
# test a wrong psfcutoff (<0, 1) 
# test mulit-term psf fitting (nterms=2)
# 
##########################################################################
 
 
####    Imports     ####
import os
import sys
import unittest
import copy
import shutil
import numpy
 
from casatestutils import testhelper as th
from casatestutils import compare  
 
from casatools import ctsys, synthesisutils, image, msmetadata, measures, table, quanta
su = synthesisutils() 
_qa = quanta()

class sutest_base(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Location of input data
        cls.datapath = ctsys.resolve('unittest/synthesisutils/')
        cls.ia = image()
        cls.tb = table()
        cls.msmd = msmetadata()
    @classmethod
    def tearDownClass(cls):
        cls.ia.done()
        cls.tb.done()
        cls.msmd.done()

    def checklistval(kind, inlist, reflist):
        state = False
        pstr = ''
        if inlist==reflist:
            state = True
        else:
            state = False
            pstr = 'Fail : {} {} expected {}'.format(kind,inlist, reflist)
        return (state, pstr)

    def checkdict(self,indict,refdict, testdesc=''):
        """
        Compare dictionaries between indict and refdict
        testesc: optional test description printed with the message for failed case
        
        Returns True/False and  messages in cases of failures 
        also print the failure messages if at least one of the tests fails
        """
        tests =[]
        msg = ''
        if testdesc!='':
           msg += '[ {} ]'.format(testdesc)
        for key in indict:
            if key in refdict:
                if isinstance(indict[key], numpy.ndarray):
                    if isinstance(refdict[key], numpy.ndarray):
                        if numpy.array_equal(indict[key],refdict[key]): 
                            tests.append(True)
                        else: 
                            #print("{0}:{1} (current) != {0}:{2} (reference) ".format(key,indict[key],refdict[key]))
                            msg += "Fail: {0}:{1} (current) != {0}:{2} (reference)\n ".format(key,indict[key],refdict[key])
                            tests.append(False)
                         
                    else:
                         msg += "Fail: Type mismatch for {}: current {}, reference {} \n".format(key,type(indict[key]), type(refdict[key]))
                elif isinstance(indict[key], list):
                     if isinstance(refdict[key], list):
                         if indict[key]==refdict[key]:
                             tests.append(True)
                         else:
                             tests.append(False)
                             msg += "Fail: {0}:{1} (current) != {0}:{2}(reference) \n".format(key, indict[key], refdict[key])
                else:
                    if indict[key]==refdict[key]:
                        tests.append(True)
                    else:
                        tests.append(False)
                        msg += "Fail {0}:{1} (current) != {0}:{2} (reference) \n".format(key, indict[key], refdict[key])
            else:
                pass
 
        allpass = all(tests)
        if not allpass:
            print(msg)

        return (allpass, msg)

    def topof_to_outframef(self, vis, infreq, outframe, fieldid):
        """
        Calculate corresponding outframe frequencies for the input TOPO
        frequency for a specific source for the given observing time stamps
        as the outframe frequency shift.
        (returns a list of measure frequency)
        """
        _me = measures() 
        # source direction
        self.tb.open(vis+'/FIELD')
        dir = self.tb.getcell('PHASE_DIR',fieldid)
        self.tb.done()
        mdir = _me.direction('J2000', str(dir[0][0])+'rad',str(dir[1][0])+'rad')
    
        # location
        self.tb.open(vis+'/OBSERVATION')
        tel = self.tb.getcol('TELESCOPE_NAME')
        self.tb.done() 
   
        # get time range for the source
        self.msmd.open(vis)
        tms = self.msmd.timesforfield(fieldid)
        self.msmd.done()
   
        # set location
        _me.doframe(_me.observatory(tel[0])) 
        # set direction
        _me.doframe(mdir)

        # input freq
        inmf = _me.frequency('TOPO',infreq)
    
        outflist=[]
        for tm in tms:
            ep = _me.epoch('UTC', str(tm)+'s')
            # set frame for a specific epoch
            _me.doframe(ep)
            #me.showframe()
    
            outf = _me.measure(inmf,outframe)

            outflist.append(outf)

        return outflist
 
####    Tests     ####
class fitpsfbeam_test(sutest_base):
    ### su.fitpsBeam has three parameters, imagename, nterms, and psfcutoff.
    ### nterms and psfcutoff have the default values but other values need to be tested.
    ### Appropriate multiterm psfs need to be present for nterms>1
   
     
    ### Set Up
    def setUp(self):
        # input images
        self.inputims=[]
        psfim1 = 'su_fitpsfbeam_test_mfs.psf'
        self.inputims.append(psfim1)
        shutil.copytree(os.path.join(self.datapath, psfim1),psfim1)
        # base name for multiterm 
        psfim2base = 'su_fitpsfbeam_test_mtmfs.psf.tt'
        for ext in ['0','1','2']:
            psfname = psfim2base+ext
            psf2ttfullpath = os.path.join(self.datapath, psfname) 
            shutil.copytree(os.path.join(self.datapath, psfname),psfname)
            self.inputims.append(psfname)
            
    ### Teardown
    def tearDown(self):
        for img in self.inputims:
            shutil.rmtree(img)
 
    ### Test Cases
    def test_mfs(self):
        '''Test that fitting of mfs psf works '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]
        self.ia.open(psfim)
        origbeam = self.ia.restoringbeam()
        # Set different values for the beam
        resetbeam = copy.deepcopy(origbeam)
        resetbeam['major']['value']=1.0
        resetbeam['minor']['value']=1.0
        resetbeam['positionangle']['value']=0.0
        self.ia.setrestoringbeam(major=resetbeam['major'], minor=resetbeam['minor'], pa=resetbeam['positionangle'])
        self.ia.done()
 
        # expected values
        bmref= {'major': {'unit': 'arcsec', 'value': 52.93011474609375}, 'minor': {'unit': 'arcsec', 'value': 49.147438049316406}, 'positionangle': {'unit': 'deg', 'value': -87.39420318603516}}
 
        su = synthesisutils()
        ret = su.fitPsfBeam(imagename=psfim[:-4])   
        del su
        self.ia.open(psfim)
        newbeam = self.ia.restoringbeam()
        self.ia.done()

        self.assertTrue(ret)
        # test against expected values
        # - original beam in psf image is different from fitted values
        self.assertDictContainsSubset(newbeam, bmref)

    def test_mfs_wrong_nterms(self):
        '''Test that it catches if nterms is inconsistent with input psf (nterms=2, for a single term  psf) '''
        su = synthesisutils()
        psfim = self.inputims[0]

        self.assertRaises(RuntimeError, su.fitPsfBeam, imagename=psfim[:-4],nterms=2)   
        del su

    def test_mfs_wrong_psfcutoff(self):
        '''Test that psfcutoff is given outside the allowed range (1.0) '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]

        su = synthesisutils()
        ret1 = su.fitPsfBeam(imagename=psfim[:-4],psfcutoff=1.0)   
        ret2 = su.fitPsfBeam(imagename=psfim[:-4],psfcutoff=-1.0)   
        del su
        self.assertFalse(ret1)
        self.assertFalse(ret2)

    def test_mfs_largerpsfcutoff(self):
        '''Test that psfcutoff with a valid (larger) number  works '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfim = self.inputims[0]
        self.ia.open(psfim)
        origbeam = self.ia.restoringbeam()
        # Set different values for the beam
        resetbeam = copy.deepcopy(origbeam)
        resetbeam['major']['value']=1.0
        resetbeam['minor']['value']=1.0
        resetbeam['positionangle']['value']=0.0
        self.ia.setrestoringbeam(major=resetbeam['major'], minor=resetbeam['minor'], pa=resetbeam['positionangle'])
        self.ia.done()
 
        # expected values
        bmref= {'major': {'unit': 'arcsec', 'value': 49.89558029174805}, 'minor': {'unit': 'arcsec', 'value': 46.80238342285156}, 'positionangle': {'unit': 'deg', 'value': -88.28898620605469}}
 
        su = synthesisutils()
        ret = su.fitPsfBeam(imagename=psfim[:-4],psfcutoff=0.5)   
        del su
        self.ia.open(psfim)
        newbeam = self.ia.restoringbeam()
        self.ia.done()

        self.assertTrue(ret)
        # test against expected values
        # - original beam in psf image is different from fitted values
        self.assertDictContainsSubset(newbeam, bmref)

    def test_mtmfs_nterms2(self):
        '''Test that fitting of multiterm  psf works '''
        # prep psf image
        # Save the origin restoringbeam for reference
        psfimtt0 = self.inputims[1]
        self.ia.open(psfimtt0)
        origbeam = self.ia.restoringbeam()
        # Set different values for the beam
        resetbeam = copy.deepcopy(origbeam)
        resetbeam['major']['value']=1.0
        resetbeam['minor']['value']=1.0
        resetbeam['positionangle']['value']=0.0
        self.ia.setrestoringbeam(major=resetbeam['major'], minor=resetbeam['minor'], pa=resetbeam['positionangle'])
        self.ia.done()
 
        # expected values
        bmref= {'major': {'unit': 'arcsec', 'value': 52.93011474609375}, 'minor': {'unit': 'arcsec', 'value': 49.147438049316406}, 'positionangle': {'unit': 'deg', 'value': -87.39420318603516}}
 
        su = synthesisutils()
        ret = su.fitPsfBeam(imagename=psfimtt0[:-8],nterms=2)   
        del su
        self.ia.open(psfimtt0)
        newbeam = self.ia.restoringbeam()
        self.ia.done()

        self.assertTrue(ret)
        # test against expected values
        # - original beam in psf image is different from fitted values
        self.assertDictContainsSubset(newbeam, bmref)
 
class getoptimumsize_test(sutest_base):
    ### su.getoptimumsize takes a single parameter, size.

    ### Set Up
    def setUp(self):
        pass
    ### Teardown
    def tearDown(self):
        pass
    ### Test Cases
    def test_default(self):
        '''Test default size'''
        su = synthesisutils()
        ret = su.getOptimumSize()
        del su
        self.assertEqual(ret,100)

    def test_oddimsize(self):
        '''Test odd non-optimal number'''
        su = synthesisutils()
        ret = su.getOptimumSize(size=501)
        del su
        self.assertEqual(ret,512)

    def test_evenimsize(self):
        '''Test even non-optimal number '''
        su = synthesisutils()
        ret = su.getOptimumSize(size=510)
        del su
        self.assertEqual(ret,512)

class advisechansel_test(sutest_base):
    ### su.advisechansel returns data selections or frequency range for spectral imaging

    ### Setup
    def setUp(self):
        # input MSes
        self.inputmses=[]
        ms1 = 'twhya.short.ms'
        ms2 = 'alma_ephemobj_icrs.ms'
        self.inputmses.append(ms1)
        self.inputmses.append(ms2)
        for inms in self.inputmses:
            if not os.path.exists(inms):
                shutil.copytree(os.path.join(self.datapath, inms),inms)

    ### Teardown
    def TearDown(self):
        for inpdata in self.inputmses:
            if os.path.exists(inpdata):
                shutil.rmtree(indata)

    ### Test cases
    def test_advisechanelsel_datasel(self):
        '''Test that data selection parameters or given frequency range, etc is returned correctly'''
        # Run multiple tests with slightly different input but should yield identical results.

        # First, try it in the data frame 
        su = synthesisutils()

        inputms ='twhya.short.ms'
        # start: Chan 500 of spw0, end: Chan 3840 of spw1
        fsttopo = '356.55897119140625GHz'
        fentopo = '358.2029418946312GHz'
        cw = '122.070kHz'
        res = su.advisechansel(freqstart=fsttopo,freqend=fentopo, freqstep=cw, freqframe='TOPO', fieldid=2, msname=inputms)   
        #print('res=',res)
        # res['nchan'] = [3340, 3840], res['spw'] = [0,1], res['start']=[500,0]
        refdict = {'nchan': numpy.array([3340, 3840]), 'spw': numpy.array([0, 1]), 'start': numpy.array([500,   0])}

        (pof, errmsg) = self.checkdict(res,refdict, 'TOPO test')

        # calculate input frequencies in LSRK in the time range for a specific source
        listoffreqs = self.topof_to_outframef(inputms, fsttopo, 'LSRK', 2)
        fval = []
        for f in listoffreqs:
            fval.append(f['m0']['value'])
        fstlsrk = str(_qa.convert(str(min(fval))+'Hz', 'GHz')['value'])+'GHz' 

        listoffreqs = self.topof_to_outframef(inputms, fentopo, 'LSRK', 2)
        fval.clear()
        for f in listoffreqs:
            fval.append(f['m0']['value'])
        fenlsrk = str(_qa.convert(str(min(fval))+'Hz', 'GHz')['value'])+'GHz' 
        #fstlsrk= 356.58364232412697GHz
        #fenlsrk= 358.22772677745405GHz
        # extra adding makes slightly larger chan range?
        reflsrkdict = {'nchan': numpy.array([3341, 3840]), 'spw': numpy.array([0, 1]), 'start': numpy.array([499,   0])}
        # run advisechansel with LSRK 
        reslsrk = su.advisechansel(freqstart=fstlsrk,freqend=fenlsrk, freqstep=cw, freqframe='LSRK', fieldid=2, msname=inputms)   
        del su

        (poflsrk, errmsglsrk) = self.checkdict(reslsrk,reflsrkdict,'LSRK test')
        self.assertTrue(all([pof,poflsrk]))

    def test_advisechanelsel_datasel_ephem(self):
        '''Test that data selection parameters for given frequency range,etc for ephemeris object is returned correctly'''

        su = synthesisutils()
        #ret = su.advisechansel(xxxx) 
        inputms = 'alma_ephemobj_icrs.ms'
        # spw2 descending freq order 
        fst = '220.34439739444292GHz'
        fen = '220.29556926944292GHz' 
        cw  = '122.0703125kHz'
        res = su.advisechansel(freqstart=fst,freqend=fen, freqstep=cw, freqframe='SOURCE', fieldid=1, ephemtable='TRACKFIELD', msname=inputms)   
        #print('res=',res)
        # res['nchan'] = [3340, 3840], res['spw'] = [0,1], res['start']=[500,0]
        refdict = {'nchan': numpy.array([402]), 'spw': numpy.array([2]), 'start': numpy.array([251])}

        (pof, errmsg) = self.checkdict(res,refdict, 'internal ephem table test')

        # alternatively use an external table
        shutil.copytree(inputms+'/FIELD/EPHEM0_Uranus_57362.91000000.tab', 'external_EPHEM0_Uranus_57362.91000000.tab')
        resext = su.advisechansel(freqstart=fst,freqend=fen, freqstep=cw, freqframe='SOURCE', fieldid=1, ephemtable='external_EPHEM0_Uranus_57362.91000000.tab', msname=inputms)   
        (pofext, errmsgext) = self.checkdict(resext,refdict,'external ephem table test')
        shutil.rmtree('external_EPHEM0_Uranus_57362.91000000.tab')

        # or alternatively use the default table of the object known by CASA
        resdef = su.advisechansel(freqstart=fst,freqend=fen, freqstep=cw, freqframe='SOURCE', fieldid=1, ephemtable='Uranus', msname=inputms)   
        (pofdef, errmsgdef) = self.checkdict(resdef,refdict,'default ephem table test')

        del su

        self.assertTrue(all([pof,pofext,pofdef]))
        

    def test_su_adivsechansel_getfreqrange(self):
        '''Test that frequency range for given data selections is returned correctly'''
        refdict = {'freqend': {'unit': 'Hz', 'value': 358203002929.7875}, 'freqstart': {'unit': 'Hz', 'value': 356558910156.25}}
        su = synthesisutils()
        inputms ='twhya.short.ms'
        res = su.advisechansel(freqframe='TOPO', getfreqrange=True,spwselection='0:500~3839,1',fieldid=2, msname=inputms)
        
        (pof, errmsg) = self.checkdict(res, refdict, 'TOPO test') 
        #fstlsrk= '356.58364232412697GHz' fenlsrk= '358.22772677745405GHz'
        # corresponding LSRK frequencies
        reslsrk = su.advisechansel(freqframe='LSRK', getfreqrange=True,spwselection='0:500~3839,1', fieldid=2, msname=inputms)

        del su

        reflsrkdict = {'freqend': {'unit': 'Hz', 'value': 358227726777.45405}, 'freqstart': {'unit': 'Hz', 'value': 356583642324.12697}}
        diff_freqstart = reslsrk['freqstart']['value']-reflsrkdict['freqstart']['value']
        diff_freqend = reslsrk['freqend']['value']-reflsrkdict['freqend']['value']
        # there may be some padding added so allow difference of ~the channel width (122.07 kHz) 
        tol = 122070.0 
        self.assertTrue(abs(diff_freqstart) < tol and diff_freqstart < 0)
        self.assertTrue(abs(diff_freqend) < tol and diff_freqend >0)
        self.assertTrue(pof)

    def test_su_advisechansel_getfreqrange_ephem(self):
        '''Test that frequency range for given data selections for an ephemeris object is returned correctly'''
        su = synthesisutils()
        inputms ='alma_ephemobj_icrs.ms'
        tol  = 122070.3125 # chan width in Hz
        refdict = {'freqend': {'unit': 'Hz', 'value': 220344661518.86124}, 'freqstart': {'unit': 'Hz', 'value': 220295510184.3427}}
        res = su.advisechansel(freqframe='SOURCE', fieldid=1, getfreqrange=True, spwselection='2:251~652', ephemtable='TRACKFIELD', msname=inputms)   
        shutil.copytree(inputms+'/FIELD/EPHEM0_Uranus_57362.91000000.tab', 'external_EPHEM0_Uranus_57362.91000000.tab')
        resext = su.advisechansel(freqframe='SOURCE', fieldid=1, getfreqrange=True, spwselection='2:251~652',  ephemtable='external_EPHEM0_Uranus_57362.91000000.tab', msname=inputms)   
     
        shutil.rmtree('external_EPHEM0_Uranus_57362.91000000.tab')

        # res and resext should be identical as it effectively uses the same ephem table
        self.assertEqual(res['freqstart']['value'], refdict['freqstart']['value'])
        self.assertEqual(res['freqend']['value'], refdict['freqend']['value'])
        self.assertEqual(resext['freqstart']['value'], refdict['freqstart']['value'])
        self.assertEqual(resext['freqend']['value'], refdict['freqend']['value'])

        # or alternatively use the default table of the object known by CASA
        resdef = su.advisechansel(freqframe='SOURCE', fieldid=1, getfreqrange=True, spwselection='2:251~652', ephemtable='Uranus', msname=inputms)   
        del su

        diff_freqstart_def = resdef['freqstart']['value'] - refdict['freqstart']['value']
        diff_freqend_def = resdef['freqend']['value'] - refdict['freqend']['value']
        self.assertTrue(abs(diff_freqstart_def) < tol )
        self.assertTrue(abs(diff_freqend_def) < tol )

    def test_su_adivsechanel_defaults(self):
        '''Test non specified parameter case for proper error/warning message '''
        su = synthesisutils()
        self.assertRaises(Exception, su.advisechansel())


if __name__ == '__main__':
    unittest.main()

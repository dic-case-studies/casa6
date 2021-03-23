##########################################################################
# test_req_task_visstat.py
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_visstat/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import visstat
    CASA6 = True
    
    tb = casatools.table()
    qa = casatools.quanta()
    me = casatools.measures()
    
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
import sys
import os
import unittest
import copy
import shutil
import numpy as np

if CASA6:
    loc_zip = zip
else:
    from itertools import izip
    loc_zip = izip
    
    
epsilon = 0.0001

### Data ###
if CASA6:
    datapath = casatools.ctsys.resolve('unittest/visstat/outlier_ut.ms/')
    mms_data = casatools.ctsys.resolve('unittest/visstat/outlier_mms.mms/')
    selectiondata = casatools.ctsys.resolve('unittest/visstat/uid___X02_X3d737_X1_01_small.ms/')
    mms_select = casatools.ctsys.resolve('unittest/visstat/uid___X02_X3d737_X1_01_small.mms')
    singledish = casatools.ctsys.resolve('unittest/visstat/analytic_spectra_tsys.ms')
    # Data for merged test
    merged_data_path = casatools.ctsys.resolve('unittest/visstat/')
    
else:
    datapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/visstat/outlier_ut.ms/'
    mms_data = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/visstat/outlier_mms.mms/'
    selectiondata = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/visstat/uid___X02_X3d737_X1_01_small.ms/'
    mms_select = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/visstat/uid___X02_X3d737_X1_01_small.mms'
    singledish = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/visstat/analytic_spectra_tsys.ms'
    # Data from merged test
    merged_data_path = os.path.join(os.environ.get('CASAPATH').split()[0], 'casatestdata/unittest/visstat/')
    
axislist = ['flag', 'antenna1', 'antenna2', 'feed1', 'feed2', 'field_id', 'array_id', 'data_desc_id', 'flag_row', 'interval', 'scan', 'scan_number', 'time', 'weight_spectrum', 'amp', 'amplitude', 'phase', 'real', 'imag', 'imaginary', 'uvrange']
 
keylist = ['firstquartile', 'isMasked', 'isWeighted', 'max', 'maxDatasetIndex', 'maxIndex', 'mean', 'medabsdevmed', 'median', 'min', 'minDatasetIndex', 'minIndex', 'npts', 'rms', 'stddev', 'sum', 'sumOfWeights', 'sumsq', 'thirdquartile', 'variance']

# NOTE: I also need to get and/or make a multimesset to test this data from. I could usr the data from listobs to do this.

if CASA6:
    ###
    ### this test uses a sort(...) "cmp" function, but python 3 uses a key function...
    ### this function from the python 3 documentation converts...
    ###
    def cmp_to_key(mycmp):
        'Convert a cmp= function into a key= function'
        class K:
            def __init__(self, obj, *args):
                self.obj = obj
            def __lt__(self, other):
                return mycmp(self.obj, other.obj) < 0
            def __gt__(self, other):
                return mycmp(self.obj, other.obj) > 0
            def __eq__(self, other):
                return mycmp(self.obj, other.obj) == 0
            def __le__(self, other):
                return mycmp(self.obj, other.obj) <= 0
            def __ge__(self, other):
                return mycmp(self.obj, other.obj) >= 0
            def __ne__(self, other):
                return mycmp(self.obj, other.obj) != 0
        return K


nostat = visstat(selectiondata)
nostatmms = visstat(mms_select)

class visstat_test(unittest.TestCase):
     
    def setUp(self):
        if not CASA6:
            default(visstat)
            
        # Data set up from merged test
        self.msfile = "ngc5921_add_corect_model.ms"
        self.msfile2 ="OrionS_rawACSmod_calave.ms"
        self.msfile2_asap="OrionS_rawACSmod_calave.asap"
        self.msfile3="OrionS_rawACSmod_calave_intent_on1_off3.ms"
        # Never used:
        # self.msfile4="OrionS_rawACSmod_calave_intent_on1.ms"
        self.msfile5="OrionS_rawACSmod_calave_intent_on3_off1.ms"
        # Never used:
        # self.msfile6="OrionS_rawACSmod_calave_intent_off1.ms"
        self.msfile_flagged_spw1='OrionS_rawACSmod_calave_flagged_spw1.ms'
        self.msfile7="visstat2_test6_scan.txt"
        self.msfile8="visstat2_test6_amp.txt"
        self.msfile9="visstat2_test5.txt"
        self.msfile10="visstat2_test7.txt"
        self.msfile11="visstat2_test8.txt"
        self.msfile13='visstat2_test10_check_on.txt'
        self.msfile14='visstat2_test10_check_off.txt'

        shutil.copytree(os.path.join(merged_data_path,self.msfile), self.msfile, symlinks=True)
        shutil.copytree(os.path.join(merged_data_path,self.msfile2), self.msfile2, symlinks=True)
        shutil.copytree(os.path.join(merged_data_path,self.msfile2_asap), self.msfile2_asap,
                        symlinks=True)
        shutil.copytree(os.path.join(merged_data_path,self.msfile3), self.msfile3, symlinks=True)
        # shutil.copytree(os.path.join(datapath,self.msfile4), self.msfile4, symlinks=True)
        shutil.copytree(os.path.join(merged_data_path,self.msfile5), self.msfile5, symlinks=True)
        # shutil.copytree(os.path.join(datapath,self.msfile6), self.msfile6, symlinks=True)
        shutil.copytree(os.path.join(merged_data_path,self.msfile_flagged_spw1),
                        self.msfile_flagged_spw1, symlinks=True)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile7), self.msfile7)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile8), self.msfile8)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile9), self.msfile9)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile10), self.msfile10)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile11), self.msfile11)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile13), self.msfile13)
        shutil.copyfile(os.path.join(merged_data_path,self.msfile14), self.msfile14)
        
        if not CASA6:
            default('visstat')
            
            
        
    def tearDown(self):
        # Teardown from merged test
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.msfile2)
        shutil.rmtree(self.msfile2_asap)
        shutil.rmtree(self.msfile3)
        # shutil.rmtree(self.msfile4)
        shutil.rmtree(self.msfile5)
        # shutil.rmtree(self.msfile6)
        shutil.rmtree(self.msfile_flagged_spw1)
        os.remove(self.msfile7)
        os.remove(self.msfile8)
        os.remove(self.msfile9)
        os.remove(self.msfile10)
        os.remove(self.msfile11)
        os.remove(self.msfile13)
        os.remove(self.msfile14)
   
   
    def test_axis(self):
        '''
            test_axis
            ------------------------------
            
            Test the axis parameter values.
            visstat should return a dict and the keys should match the provided key list.
            
            This test iterates over all the possible axis values
        '''
        for axis in axislist:

            axisstat = visstat(datapath, axis=axis)
            axismms = visstat(mms_data, axis=axis)
            
            self.assertTrue( type(axisstat) == type(dict()), msg='output is not a dict for axis {}'.format(axis) )
            self.assertTrue( type(axismms) == type(dict()), msg='output is not a dict for axis {} on mms'.format(axis) )
        
            self.assertTrue(sorted(list(axisstat[list(axisstat.keys())[0]].keys())) == keylist, msg='keys do not match the key list for axis {}'.format(axis))
            self.assertTrue(sorted(list(axismms[list(axismms.keys())[0]].keys())) == keylist, msg='keys do not match the key list for axis {} for mms'.format(axis))
    
            
    def test_reportingaxes(self):
        '''
            test_reportingaxes
            -----------------------------
            
            Test the reportingaxes parameter.
            The output should be a dict and contain all the expected keys.
            
            Iterate over all the possible values.
        '''
        
        for axes in ['ddid', 'field', 'integration']:
            
            reportstat = visstat(datapath, reportingaxes=axes)
            reportmms = visstat(mms_data, reportingaxes=axes)
            
            self.assertTrue( type(reportstat) == type(dict()) )
            self.assertTrue( type(reportmms) == type(dict()) )
            
            self.assertTrue(sorted(list(reportstat[list(reportstat.keys())[0]].keys())) == keylist)
            self.assertTrue(sorted(list(reportmms[list(reportmms.keys())[0]].keys())) == keylist)
            
    def test_useflags(self):
        '''
            test_useflags
            ----------------------
            
            Check that the useflags parameter produces different results then when useflags = False.
        '''
        
        withflags = visstat(datapath, useflags=True)
        withoutflags = visstat(datapath, useflags = False)
        
        flagsmms = visstat(mms_select, useflags=False)
        noflagmms = visstat(mms_select, useflags=True)
        
        self.assertTrue(withflags != withoutflags)
        self.assertTrue(flagsmms != noflagmms)
        
    def test_datacolumn(self):
        '''
            test_datacolumn
            ----------------------------
            
            Check the data column parameter.
            
            Iterate over possible data column inputs and check that a dictionary is created by visstat
            also check that all the keys that should be present are there.
            
            (This last step may not be nessisary for this test)
        '''
        
        for col in [ 'data', 'corrected', 'model' ]:
            colstat = visstat(datapath, datacolumn=col)
            colmms = visstat(mms_data, datacolumn=col)
            
            self.assertTrue( type(colstat) == type(dict()), msg = 'Fails for column: {}'.format(col) )
            self.assertTrue( type(colmms) == type(dict()), msg = 'Fails for column: {} on mms'.format(col) )
            
            self.assertTrue(sorted(list(colstat[list(colstat.keys())[0]].keys())) == keylist, msg = 'Fails for column: {}'.format(col))
            self.assertTrue(sorted(list(colmms[list(colmms.keys())[0]].keys())) == keylist, msg = 'Fails for column: {} on mms'.format(col))
            
        floatstat = visstat(singledish, datacolumn='float_data')
        self.assertTrue( type(floatstat) == type(dict()) )
            
    def test_spw(self):
        '''
            test_spw
            ---------------
            
            Test the spectral window selection parameter.
            
            Assert that a selection using the spw parameter returns a different result than no selection.
        '''
        
        spwstat = visstat(selectiondata, spw='0')
        
        spwmms = visstat(mms_select, spw='0')
        
        self.assertTrue(spwstat != nostat)
        self.assertTrue(spwmms != nostatmms)
        
    def test_field(self):
        '''
            test_field
            --------------
            
            Test the field selection parameter.
            
            Assert that a selection using the field parameter returns a different result than no selection.
        '''
        
        fieldstat = visstat(selectiondata, field='0')
        
        fieldmms = visstat(mms_select, field='0')
        
        self.assertTrue( fieldstat != nostat )
        self.assertTrue( fieldmms != nostatmms )
        
    def test_selectdata(self):
        '''
            test_selectdata
            -----------------------
            
            Test the selectdata parameter
            
            Assert that the select data parameter prevents other selection fields from having an affect
            
            Assert that with selectdata=False and an active selection produces the same results as the task with no selection
        '''
        
        selectstat = visstat(selectiondata, selectdata=False, antenna='DV01')
        
        selectmms = visstat(mms_select, selectdata=False, antenna='DV01')
        
        self.assertTrue( selectstat == nostat )
        self.assertTrue( selectmms == nostatmms )
        
    def test_antenna(self):
        '''
            test_antenna
            ---------------------------
            
            Test the antenna selection parameter
            
            Assert that selection with this parameter will return a different result than no selection.
        '''
        
        antennastat = visstat(selectiondata, antenna='DV01')
        
        antennamms = visstat(mms_select, antenna='DV01')
        
        self.assertTrue( antennastat != nostat )
        self.assertTrue( antennamms != nostatmms )
        
    def test_uvange(self):
        '''
            test_uvrange
            --------------------
            
            Test the uvrange selection parameter
            
            Assert that the selection with this parameter will return a different value than no selection.
        '''
        
        uvRangeStat = visstat(selectiondata, uvrange='0~10')
        
        uvRangemms = visstat(mms_select, uvrange='0~10')
        
        self.assertTrue( uvRangeStat != nostat )
        self.assertTrue( uvRangemms != nostatmms)
        
    def test_timerange(self):
        '''
            test_timerange
            ------------------
            
            Test the timerange selection parameter
            
            Assert that the selection with this parameter will return a different value than no selection.
        '''
        
        timerangeStat = visstat(selectiondata, timerange='03:01:30~03:05:00')
        
        timerangemms = visstat(mms_select, timerange='03:01:30~03:05:00')
        
        self.assertTrue( timerangeStat != nostat )
        self.assertTrue( timerangemms != nostatmms )
        
    def test_correlation(self):
        '''
            test_correlation
            -------------------
            
            Test the correlation parameter
            
            Assert that the selection with this parameter will return a different value than no selection.
        '''
        
        corrStat = visstat(selectiondata, correlation='XX')
        
        corrmms = visstat(mms_select, correlation='XX')
        
        self.assertTrue( corrStat != nostat )
        self.assertTrue( corrmms != nostatmms )
        
    def test_scan(self):
        '''
            test_scan
            ------------
        
            Test the scan selection parameter
        
            Assert that the selction with this parameter will return a different result than no selection.
        '''
    
        scanStat = visstat(selectiondata, scan='1')
        
        scanmms = visstat(mms_select, scan='1')
        
        self.assertTrue( scanStat != nostat )
        self.assertTrue( scanmms != nostatmms )
        
    def test_array(self):
        '''
            test_array
            -------------
            
            Test the array selection parameter.
            
            Assert that checking an out of range array returns a NoneType, and valid selections retrun a dictionary
        '''
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(RuntimeError):
                arrayFail = visstat(selectiondata, array='1')
            with self.assertRaises(RuntimeError):
                arrayFailmms = visstat(mms_select, array='1')
        else:
            arrayFail = visstat(selectiondata, array='1')
            arrayFailmms = visstat(mms_select, array='1')
            
            self.assertTrue( type(arrayFail) == type(None) )
            self.assertTrue( type(arrayFailmms) == type(None) )
        
        arrayPass = visstat(selectiondata, array='0')
        arrayPassmms = visstat(mms_select, array='0')
        
        self.assertTrue( type(arrayPass) == type(dict()) )
        self.assertTrue( type(arrayPassmms) == type(dict()) )
        
    def test_observation(self):
        '''
            test_observation
            -----------------------
            
            Test the observation selection parameter
            
            Assert that checking an out of range observation ID returns a NoneType, and valid selections return a dictionary
        '''
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(RuntimeError):
                observationFail = visstat(selectiondata, observation=1)
            with self.assertRaises(RuntimeError):
                observationFailmms = visstat(mms_select, observation=1)
        else:
            observationFail = visstat(selectiondata, observation=1)
            observationFailmms = visstat(mms_select, observation=1)
            self.assertTrue( type(observationFail) == type(None) )
            self.assertTrue( type(observationFailmms) == type(None) )
        
        observationPass = visstat(selectiondata, observation=0)
        observationPassmms = visstat(mms_select, observation=0)
    
        self.assertTrue( type(observationPass) == type(dict()) )
        self.assertTrue( type(observationPassmms) == type(dict()) )
        
    def test_timeavg(self):
        '''
            test_timeaverage
            ---------------------
            
            Test the timeaverage parameter
            
            Assert that the dict produced when timeaverage = True is different from the one produced when timeaverage=False
        '''

        timeavgTrue = visstat(selectiondata, timebin='1s', timeaverage=True)
        timeavgFalse = visstat(selectiondata, timeaverage=False)
        
        timeavgTruemms = visstat(mms_select, timebin='1s', timeaverage=True)
        timeavgFalsemms = visstat(mms_select, timeaverage=False)
        
        self.assertTrue( timeavgFalse != timeavgTrue )
        self.assertTrue( timeavgFalsemms != timeavgTruemms)
        
    def test_timebin(self):
        '''
            test_timebin
            -------------------
            
            Test the timebin parameter
            
            Assert that the result when given a bin width for averaging is different than when none is given.
        '''
        
        timebinSelect = visstat(selectiondata, timeaverage=True, timebin='10s')
        
        timebinmms = visstat(selectiondata, timeaverage=True, timebin='10s')
        
        self.assertTrue( timebinSelect != nostat )
        self.assertTrue( timebinmms != nostatmms )
        
    def test_timespan(self):
        '''
            test_timespan
            ------------------
            
            Test the timespan parameter
            
            Assert that all parameter settings give different results than the default output
        '''

        scanSelect = visstat(selectiondata, timeaverage=True, timebin='1s', timespan='scan')
        stateSelect = visstat(selectiondata, timeaverage=True, timebin='1s',
                              timespan='state')
        bothSelect = visstat(selectiondata, timeaverage=True,  timebin='1s',
                             timespan='scan, state')
        
        scanSelectmms = visstat(mms_select, timeaverage=True,  timebin='1s', timespan='scan')
        stateSelectmms = visstat(mms_select, timeaverage=True,  timebin='1s',
                                 timespan='state')
        bothSelectmms = visstat(mms_select, timeaverage=True,  timebin='1s',
                                timespan='scan, state')
        
        for item in [scanSelect, stateSelect, bothSelect]:
            self.assertTrue( item != nostat )
            
        for item in [scanSelectmms, stateSelectmms, bothSelectmms]:
            self.assertTrue( item != nostatmms )
        
    def test_maxuvwdistance(self):
        '''
            test_maxuvwdistance
            -----------------------
            
            Test the maxuvwdistance parameter
            
            Assert that the output is a python dict. Once again this selection seems to not change the values that are returned
        '''

        uvwSelect = visstat(selectiondata, timeaverage=True, timebin='1s',
                            maxuvwdistance=10.0)
        
        uvwSelectmms = visstat(mms_select, timeaverage=True, timebin='1s',
                               maxuvwdistance=10.0)
        
        self.assertTrue( uvwSelect != nostat )
        self.assertTrue( uvwSelectmms != nostatmms )
        
    def test_intent(self):
        '''
            test_intent
            -----------------
            
            Test the intent parameter
            
            Assert that the specified selection creates a different dict than the default values.
        '''
        
        intentSelect = visstat(selectiondata, intent='CALIBRATE_AMPLI.ON_SOURCE')
        
        intentSelectmms = visstat(mms_select, intent='CALIBRATE_AMPLI.ON_SOURCE')
        
        self.assertTrue( intentSelect != nostat )
        self.assertTrue( intentSelectmms != nostatmms )
        
    # The merged test cases from test_visstat begin here
    
    def compare(self, a, b):
        for d1, d2 in loc_zip(a,b):
            if(d1.split(':')[0]==d2.split(':')[0]):
                if(not np.allclose(np.array([float(d1.split(':')[1])]), np.array([float(d2.split(':')[1])]))):
                    raise Exception(d1.split(':')[0] + ' ' + 'values are not consistent!')
                    
    def test_defaultValues(self):
        '''Visstat 01: Default values'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        expected = {self.msfile:
                    {'DATA_DESC_ID=0': {'firstquartile': 0.024146264419,
                                        'isMasked': True,
                                        'isWeighted': False,
                                        'max': 73.75,
                                        'maxDatasetIndex': 12,
                                        'maxIndex': 2408,
                                        'mean': 4.837103133618731,
                                        'medabsdevmed': 0.04501341888681054,
                                        'median': 0.05355948396027088,
                                        'min': 2.2130521756480448e-05,
                                        'minDatasetIndex': 54,
                                        'minIndex': 8692,
                                        'npts': 2660994.0,
                                        'thirdquartile': 0.3291134536266327,
                                        'rms': 17.081207832906546,
                                        'stddev': 16.382008276126726,
                                        'sum': 12871502.415939873,
                                        'sumsq': 776391995.3973862,
                                        'variance': 268.3701951590845}}}

        v2 = visstat(vis=self.msfile, axis='amp', datacolumn='data', reportingaxes='ddid')

        if v2.keys() != expected[self.msfile].keys():
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Wrong dictionary keys. Expected %s, got %s" % \
                            (expected[self.msfile], v2)
            raise Exception("Wrong dictionary keys. Expected %s, got %s" % \
                            (expected[self.msfile], v2))




        if 'DATA_DESC_ID=0' not in v2:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Dictionary returned from visstat does not have key DATA_DESC_ID=0"
            raise Exception("Dictionary returned from visstat does not have key DATA_DESC_ID=0")

        for e in expected[self.msfile]['DATA_DESC_ID=0'].keys():
            print('')
            print("Checking %s: %s vs %s" % \
                   (e, expected[self.msfile]['DATA_DESC_ID=0'][e], v2['DATA_DESC_ID=0'][e]))
            failed = False

            if type(expected[self.msfile]['DATA_DESC_ID=0'][e])==bool:
                if expected[self.msfile]['DATA_DESC_ID=0'][e] != v2['DATA_DESC_ID=0'][e]:
                    failed = True
            elif abs((expected[self.msfile]['DATA_DESC_ID=0'][e] - v2['DATA_DESC_ID=0'][e])/expected[self.msfile]['DATA_DESC_ID=0'][e]) > epsilon:
                failed = True

            if failed:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                        +"\nError: Numbers differ, expected %s, got %s" % \
                        ( str(v2['DATA_DESC_ID=0'][e]), str(expected[self.msfile]['DATA_DESC_ID=0'][e]) )

        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
    def test_channelSelectFlags(self):
        '''Visstat 02: Check channel selections, useflags=True, reportingaxes='ddid',correlation=corr, datacolumn=data, axis=amp'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        for ch in [1, 2, 4, 7, 13, 62]:
          for corr in ['ll', 'rr', 'll,rr']:
            print("Call with spw='0:1~"+str(ch)+"', correlation="+corr)
            s2 = visstat(vis=self.msfile, axis='amp', datacolumn='data', spw='0:1~'+str(ch), correlation=corr, reportingaxes='ddid', useflags=True)
            print('')
            print('s2', s2)
            n_expected = 2660994/63 * ch
            if corr in ['ll', 'rr']:
                n_expected /= 2
            n = int(s2['DATA_DESC_ID=0']['npts'])
            print("Checking npts: %s vs %s" % (n, n_expected))
            if n != n_expected:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError:"+str(n_expected) + " points expected, but npts = " + str(n)
                raise Exception(str(n_expected) + " points expected, but npts = " + str(n))
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
    def test_datacolModel(self):
        '''Visstat 03: Default values with datacolum=model, reportingaxis=ddid'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        expected = {self.msfile:
                    {'DATA_DESC_ID=0': {'firstquartile': 1.0,
                                        'isMasked': True,
                                        'isWeighted': False,
                                        'max': 1.0,
                                        'maxDatasetIndex': 0,
                                        'maxIndex': 0,
                                        'mean': 1.0,
                                        'medabsdevmed': 0.0,
                                        'median': 1.0,
                                        'min': 1.0,
                                        'minDatasetIndex': 0,
                                        'minIndex': 0,
                                        'npts': 2660994.0,
                                        'thirdquartile': 1.0,
                                        'rms': 1.0,
                                        'stddev': 0.0,
                                        'sum': 2660994.0,
                                        'sumsq': 2660994.0,
                                        'variance': 0.0}}}

        v2 = visstat(vis=self.msfile, axis='amp', datacolumn='model', reportingaxes='ddid')
        if v2.keys() != expected[self.msfile].keys():
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Wrong dictionary keys. Expected %s, got %s" % \
                            (expected[self.msfile], v2)

        if 'DATA_DESC_ID=0' not in v2:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Dictionary returned from visstat does not have key DATA_DESC_ID=0"
            raise Exception("Dictionary returned from visstat does not have key DATA_DESC_ID=0")

        for e in expected[self.msfile]['DATA_DESC_ID=0'].keys():
            print('')
            print("Checking %s: %s vs %s" % \
                   (e, expected[self.msfile]['DATA_DESC_ID=0'][e], v2['DATA_DESC_ID=0'][e]))
            failed = False
            if expected[self.msfile]['DATA_DESC_ID=0'][e] == 0:
                if v2['DATA_DESC_ID=0'][e] != 0:
                    failed = True
            else:
                if abs((expected[self.msfile]['DATA_DESC_ID=0'][e] - v2['DATA_DESC_ID=0'][e])/expected[self.msfile]['DATA_DESC_ID=0'][e]) > epsilon:
                    failed = True
            if failed:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Numbers differ, expected2 %s, got %s" % \
                      (str(expected[self.msfile]['DATA_DESC_ID=0'][e]), str(v2['DATA_DESC_ID=0'][e]))

        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
    def test_special_cases(self):
        '''Visstat 04: Test of special cases'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        for a in range(1, 5):
            s2 = visstat(vis=self.msfile, axis='antenna1', antenna=str(a)+'&26')
            print('')
            print("antenna =", a, "; mean = ", s2['DATA_DESC_ID=0']['mean'])

            # Note there's a counting from 0 or 1 issue here
            # with the antenna numbering
            if self.msfile == self.msfile:
                offset = 1
            else:
                offset = 0

            if abs((s2['DATA_DESC_ID=0']['mean']+offset) - a) > epsilon:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                        +"\nError: Failed when antenna= "+str(a)+'&26'
                raise Exception("Error!")

        for scan in range(1, 8):
            s2 = visstat(vis=self.msfile, axis='scan_number', scan=str(scan))
            if abs(s2['DATA_DESC_ID=0']['mean'] - scan) > epsilon:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                +"\nError: Failed for scan = "+str(scan)
                raise Exception("Error!")

        self.assertTrue(retValue['success'],retValue['error_msgs'])



    def test_reportAxisInt(self):
        '''Visstat 05: Test using reportingaxes=integration, datacolumn=float_data'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        correlation_type=['RR', 'LL']
        sd_correlation_type=['0', '1']
        datacolumn_list=['float_data']
        spw_list=['0', '1', '2', '3']
        reporting_axes=['integration']
        check_list=['rms', 'max', 'min', 'sum', 'median', 'stddev', 'mean']
        intent_list=['']

        tb.open(self.msfile2)
        ref=tb.getcolkeyword('TIME','MEASINFO')['Ref']
        tt=tb.getcol('TIME')
        tb.close()

        result_list=[]

        for dt in datacolumn_list:
            for col, sd_pol in loc_zip(correlation_type, sd_correlation_type):
                num_tt=0
                for time in tt:
                    for spwin in spw_list:
                        trange = qa.time(me.epoch('ref','%fs' % time)['m0'], prec=8, form='ymd')[0]
                        v2 = visstat(vis=self.msfile2, axis='amp', timerange=str(trange),reportingaxes=reporting_axes[0],
                                      correlation=col, datacolumn=dt, spw=spwin, intent=intent_list[0])
                        v2_keys=list(v2.keys())
                        for check in check_list:
                            print(check, v2[str(v2_keys[0])][check])
                            result_list.append(check+':' + str(v2[str(v2_keys[0])][check]))

                        num_tt +=1
                        if num_tt==3:
                            break

        with open('visstat2_test5.txt','r') as fout:
            f = fout.read()
        self.compare(np.array(result_list), np.array(f[:-1].split(' ')))
        



    def test_reportAxisField(self):
        '''Visstat 06: Test using reportingaxes=field'''

        datacolumn_list=['data', 'corrected', 'model']
        correlation_type=['RR', 'LL']
        field_list=['0','1','2']
        spw_list=['0']
        axis_list=['amp', 'scan_number']
        field_list=['0','1','2']
        useflags_list=[True, False]
        reporting_axes=['field']
        check_list=['rms', 'max', 'min', 'sum', 'median', 'stddev', 'mean']

        ax_scan=[]
        ax_amp=[]

        for ax in axis_list:
            for col in correlation_type:
                for fd in field_list:
                    for spwin in spw_list:
                        for fg in useflags_list:

                            if(ax=='scan_number'):
                                v2 = visstat(vis=self.msfile, axis=ax, useflags=fg, correlation=col,
                                      spw=spwin, reportingaxes=reporting_axes[0])
                                for check in check_list:
                                    ax_scan.append(check+':'+str(v2['FIELD_ID='+ fd][check]))

                            if(ax=='amp'):
                                for dt in datacolumn_list:
                                    v2 = visstat(vis=self.msfile, axis=ax, useflags=fg, datacolumn=dt, correlation=col,
                                      spw=spwin, reportingaxes=reporting_axes[0])
                                    for check in check_list:
                                        ax_amp.append(check+':'+str(v2['FIELD_ID='+ fd][check]))

        with open('visstat2_test6_scan.txt','r') as fout:
            f_scan = fout.read()
        with open('visstat2_test6_amp.txt','r') as fout:
            f_amp = fout.read()

        #check when ax=scan
        self.compare(np.array(ax_scan),np.array(f_scan[:-1].split(' ')))
        #check when ax~amp
        self.compare(np.array(ax_amp),np.array(f_amp[:-1].split(' ')))



    def test_datacolCorrected(self):
        '''Visstat 07: Default values with datacolum=corrected, reportingaxis=ddid'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        expected = {self.msfile:
                    {'DATA_DESC_ID=0': {'firstquartile': 0.024146264419,
                                        'isMasked': True,
                                        'isWeighted': False,
                                        'max': 73.75,
                                        'maxDatasetIndex': 12,
                                        'maxIndex': 2408,
                                        'mean': 4.837103133618731,
                                        'medabsdevmed': 0.04501341888681054,
                                        'median': 0.05355948396027088,
                                        'min': 2.2130521756480448e-05,
                                        'minDatasetIndex': 54,
                                        'minIndex': 8692,
                                        'npts': 2660994.0,
                                        'thirdquartile': 0.3291134536266327,
                                        'rms': 17.081207832906546,
                                        'stddev': 16.382008276126726,
                                        'sum': 12871502.415939873,
                                        'sumsq': 776391995.3973862,
                                        'variance': 268.3701951590845}}}

        v2 = visstat(vis=self.msfile, axis='amp', datacolumn='corrected', reportingaxes='ddid')
        if v2.keys() != expected[self.msfile].keys():
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Wrong dictionary keys. Expected %s, got %s" % \
                            (expected[self.msfile], v2)
            raise Exception("Wrong dictionary keys. Expected %s, got %s" % \
                            (expected[self.msfile], v2))

        if 'DATA_DESC_ID=0' not in v2:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Dictionary returned from visstat does not have key DATA_DESC_ID=0"
            raise Exception("Dictionary returned from visstat does not have key DATA_DESC_ID=0")

        for e in expected[self.msfile]['DATA_DESC_ID=0'].keys():
            print('')
            print("Checking %s: %s vs %s" % \
                   (e, expected[self.msfile]['DATA_DESC_ID=0'][e], v2['DATA_DESC_ID=0'][e]))
            failed = False
            if expected[self.msfile]['DATA_DESC_ID=0'][e] == 0:
                if v2['DATA_DESC_ID=0'][e] != 0:
                    failed = True
            else:
                if abs((expected[self.msfile]['DATA_DESC_ID=0'][e] - v2['DATA_DESC_ID=0'][e])/expected[self.msfile]['DATA_DESC_ID=0'][e]) > epsilon:
                    failed = True
            if failed:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Numbers differ, expected2 %s, got %s" % \
                      (str(expected[self.msfile]['DATA_DESC_ID=0'][e]), str(v2['DATA_DESC_ID=0'][e]))


        self.assertTrue(retValue['success'],retValue['error_msgs'])



    def test_datacolMulti(self):
        '''Visstat 08: Test when using reportingaxes='integration, datacolumn=data,corrected,model'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        #clearcal(self.msfile, addmodel=True)
        correlation_type=['','LL','RR']
        datacolumn_list=['data', 'corrected', 'model']

        tb.open(self.msfile)
        ref=tb.getcolkeyword('TIME','MEASINFO')['Ref']
        tt=tb.getcol('TIME')
        tb.close()

        trange = qa.time(me.epoch('ref','%fs' % tt[0])['m0'], prec=8, form='ymd')[0]
        s2 = visstat(vis=self.msfile, axis='amp', timerange=str(trange),reportingaxes='integration')
        s2_keys=list(s2.keys())

        check_list=['rms', 'medabsdevmed', 'min', 'max', 'sum', 'median', 'sumsq', 'stddev', 'variance', 'npts', 'mean']
        result_list=[]

        for dt in datacolumn_list:
            for col in correlation_type:
                num_tt=0
                for time in tt:
                    trange = qa.time(me.epoch('ref','%fs' % time)['m0'], prec=8, form='ymd')[0]
                    s2 = visstat(vis=self.msfile, axis='amp', timerange=str(trange),reportingaxes='integration', correlation=col, datacolumn=dt)
                    for check in check_list:
                        result_list.append(check+':' + str(s2[str(s2_keys[0])][check]))

                    num_tt +=1
                    if num_tt==10:
                        break

        with open('visstat2_test8.txt','r') as fout:
            f = fout.read()
        self.compare(np.array(result_list), np.array(f[:-1].split(' ')))



    def test_corrLLRR(self):
        '''Visstat 09: Test using reportingaxes=ddid, correlation=[LL,RR], datacolumn=float_data spw=[0,1,2,3]'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        correlation_type=['RR', 'LL']
        sd_correlation_type=['0', '1']
        datacolumn_list=['float_data']
        spw_list=['0', '1', '2', '3']
        reporting_axes=['ddid']
        check_list=['rms', 'max', 'min', 'sum', 'median', 'stddev', 'mean']
        intent_list=['']

        tb.open(self.msfile2)
        ref=tb.getcolkeyword('TIME','MEASINFO')['Ref']
        tt=tb.getcol('TIME')
        tb.close()

        result_list=[]

        for dt in datacolumn_list:
            for col, sd_pol in loc_zip(correlation_type, sd_correlation_type):
                num_tt=0
                for time in tt:
                    for spwin in spw_list:
                        trange = qa.time(me.epoch('ref','%fs' % time)['m0'], prec=8, form='ymd')[0]
                        v2 = visstat(vis=self.msfile2, axis='amp', timerange=str(trange),reportingaxes=reporting_axes[0],
                                      correlation=col, datacolumn=dt, spw=spwin, intent=intent_list[0])
                        v2_keys=list(v2.keys())
                        for check in check_list:
                            result_list.append(check+':' + str(v2[str(v2_keys[0])][check]))
                        num_tt +=1
                        if num_tt==3:
                            break

        with open('visstat2_test7.txt','r') as fout:
            f = fout.read()
        self.compare(np.array(result_list), np.array(f[:-1].split(' ')))



    def test_intentOnOff(self):
        '''Visstat 10: Test using reportingaxes=field, datacolumn=corrected, intent=[on,off]'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        correlation_type=['RR','LL']
        sd_correlation_type=['0', '1']
        datacolumn_list=['corrected']
        spw_list=['0']
        field_list=['1']
        reporting_axes=['field']
        check_list=['rms', 'max', 'min', 'sum', 'median', 'stddev', 'mean']

        intent_on=[]
        intent_off=[]

        for dt in datacolumn_list:
            for col, sd_pol in loc_zip(correlation_type, sd_correlation_type):
                    for fd in field_list:
                        v2_intent_on = visstat(vis=self.msfile3, axis='real',reportingaxes=reporting_axes[0],
                                         correlation=col, datacolumn=dt, intent='OBSERVE_TARGET#ON_SOURCE',field=fd)
                        v2_intent_off = visstat(vis=self.msfile5, axis='real',reportingaxes=reporting_axes[0],
                                         correlation=col, datacolumn=dt, intent='OBSERVE_TARGET#OFF_SOURCE',field=fd)
                        for check in check_list:
                            print('')
                            print('check intent on', check)
                            print(v2_intent_on['FIELD_ID='+ fd][check])

                            print('check intent off', check)
                            print(v2_intent_off['FIELD_ID='+ fd][check])

                            intent_on.append(check+':'+str(v2_intent_on['FIELD_ID='+ fd][check]))
                            intent_off.append(check+':'+str(v2_intent_off['FIELD_ID='+ fd][check]))

        with open('visstat2_test10_check_on.txt','r') as fout:
            f_intent_on = fout.read()
        f_on_split=f_intent_on[:-1].split(' ')
        self.compare(np.array(f_on_split), np.array(intent_on))

        with open('visstat2_test10_check_off.txt','r') as fout:
            f_intent_off = fout.read()
        f_off_split=f_intent_off[:-1].split(' ')
        self.compare(np.array(f_off_split), np.array(intent_off))

    def test_timeAvgWithinScans(self):
        '''Visstat 11: Test of time averaging within scans'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        period = 120

        # get list of scan numbers in test ms
        tb.open(self.msfile)
        all_scans = list(set(tb.getcol('SCAN_NUMBER')))
        all_scans.sort()
        tb.close()

        v2 = visstat(vis=self.msfile, axis='amp', reportingaxes='ddid',
                      datacolumn='data', useflags=False, timeaverage=True,
                      timebin='%ds' % period)
        def scanTime(x):
            '''Get scan number and time from a single statistics result set
            in dictionary returned from visstat'''
            xparts = x.split(',')
            scan = [v for v in xparts if v.startswith('SCAN_NUMBER=')][0]
            time = [v for v in xparts if v.startswith('TIME=')][0]
            return (int(scan[len('SCAN_NUMBER='):]), float(time[len('TIME='):]))

        def compareKeys(x, y):
            '''Comparison of visstat dictionary keys: scan number, then time'''
            if CASA6:
                x = scanTime(x)
                y = scanTime(y)
                return (x>y)-(x<y)
            else:
                return cmp(scanTime(x), scanTime(y))

        # sort the dictionary keys, and create an ordered list of dictionary
        # elements (i.e, statistics for every averaging interval)
        v2_keys = v2.keys()
        if CASA6:
            v2_keys = sorted(v2_keys,key=cmp_to_key(compareKeys))
        else:
            v2_keys.sort(cmp=compareKeys)
        ordered_v2 = [(scanTime(k), v2[k]) for k in v2_keys]

        # get ordered list of scan numbers in visstat results
        scans = list(set(s for ((s, t), v) in ordered_v2))
        scans.sort()

        # since averaging does not cross scan boundaries, the "scans" and
        # "all_scans" list should be the same
        if scans != all_scans:
            retValue['success'] = False
            retValue['error_msgs'] += "\nError: Expected statistics for all " \
                                      "scans %s, got scans %s" % \
                                      (str(all_scans), str(scans))

        # within every scan, the reported statistics should have timestamps that
        # differ by the averaging period
        for scan in scans:
            times = [int(round(t)) for ((s, t), v) in ordered_v2 if s == scan]
            if len(times) > 1:
                for i in range(2, len(times)):
                    diff = times[i] - times[i - 1]
                    if diff != period:
                        retValue['success'] = False
                        retValue['error_msgs'] += "\nError: Timestamp difference for scan %d (%d) " \
                                                  "are not equal to averaging period (%d)" % \
                                                  (scan, diff, period)

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test_timeAvgAcrossScans(self):
        '''Visstat 12: Test of time averaging across scans'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        period = 120

        v2 = visstat(vis=self.msfile, axis='amp', reportingaxes='ddid',
                      datacolumn='data', useflags=False, timeaverage=True,
                      timebin='%ds' % period, timespan='scan')

        def statTime(x):
            '''Get the time from a single statistics result set in dictionary returned from
            visstat'''
            xparts = x.split(',')
            time = [v for v in xparts if v.startswith('TIME=')][0]
            return float(time[len('TIME='):])

        def compareKeys(x, y):
            '''Comparison of visstat dictionary keys by time'''
            if CASA6:
                x = statTime(x)
                y = statTime(y)
                return (x>y)-(x<y)
            else:
                return cmp(statTime(x), statTime(y))

        # sort the dictionary keys, and create an ordered list of dictionary
        # elements (i.e, statistics for every averaging interval)
        v2_keys = list(v2.keys())
        if CASA6:
            v2_keys = sorted(v2_keys,key=cmp_to_key(compareKeys))
        else:
            v2_keys.sort(cmp=compareKeys)
        ordered_v2 = [(statTime(k), v2[k]) for k in v2_keys]

        # reported statistics should have timestamps that differ by, at least,
        # the averaging period (don't check that differences are equal to the
        # averaging period because of time gaps in the data)
        times = [int(round(t)) for (t, v) in ordered_v2]
        for i in range(2, len(times)):
            diff = times[i] - times[i - 1]
            if diff < period:
                retValue['success'] = False
                retValue['error_msgs'] += "\nError: Timestamp difference (%d) " \
                                          "not at least averaging period (%d)" % \
                                          (diff, period)

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test_handle_all_flagged_groups(self):
        '''visstat 13: handle all-flagged sub-selections in reportingaxes, CAS-12857'''
            
        # The MS used in this test is the result of flagging spw=1, as follows:
        # flagdata(vis='OrionS_rawACSmod_calave.ms', mode='manual', spw='1')
        res = visstat(vis=self.msfile_flagged_spw1, axis='amp', datacolumn='corrected',
                      useflags=True, reportingaxes='ddid')
            
        # Check output dict for empty groups/sub-selection across reportingaxes (spw)
        flagged_spw = 'DATA_DESC_ID=1'
        self.assertEqual(res[flagged_spw]['npts'], 0)
        for entry in ['firstquartile', 'mean', 'medabsdevmed', 'median', 'rms', 'stddev',
                        'sum', 'sumOfWeights', 'sumsq', 'thirdquartile', 'variance']:
            self.assertTrue(np.isnan(res[flagged_spw][entry]))
        for entry in ['isMasked', 'isWeighted']:
            self.assertEqual(res[flagged_spw][entry], False)
        # Basic check on other SPWs
        for entry in ['DATA_DESC_ID=0', 'DATA_DESC_ID=2', 'DATA_DESC_ID=3']:
            self.assertEqual(res[entry]['npts'], 16384)

    def test_doquantiles(self):
        """test doquantiles parameter"""
        # doquantiles=True, all stats computed
        res = visstat(
            vis=self.msfile, axis='amp', datacolumn='data',
            useflags=True, reportingaxes='ddid', doquantiles=True
        )
        self.assertTrue(
            res['DATA_DESC_ID=0'] and type(res['DATA_DESC_ID=0']) == dict,
            'res["DATA_DESC_ID=0"] does not exist or is not a dict' 
        )
        res2 = res['DATA_DESC_ID=0']
        expected_keys = [
            'firstquartile', 'isMasked', 'isWeighted', 'max', 'maxDatasetIndex',
            'maxIndex', 'mean', 'medabsdevmed', 'median', 'min', 'minDatasetIndex',
            'minIndex', 'npts', 'rms', 'stddev', 'sum', 'sumOfWeights', 'sumsq',
            'thirdquartile', 'variance'
        ]
        # for python3 must cast dict.keys() to list
        got_keys = list(res2.keys())
        got_keys.sort()
        self.assertTrue(
            got_keys == expected_keys,
            'list of keys for doquantiles=True is incorrect'
        )

        # doquantiles=False, no quantiles computed
        res = visstat(
            vis=self.msfile, axis='amp', datacolumn='data',
            useflags=True, reportingaxes='ddid', doquantiles=False
        )
        self.assertTrue(
            res['DATA_DESC_ID=0'] and type(res['DATA_DESC_ID=0']) == dict,
            'res["DATA_DESC_ID=0"] does not exist or is not a dict'
        )
        res2 = res['DATA_DESC_ID=0']
        expected_keys = [
            'isMasked', 'isWeighted', 'max', 'maxDatasetIndex', 'maxIndex',
            'mean', 'min', 'minDatasetIndex', 'minIndex', 'npts', 'rms',
            'stddev', 'sum', 'sumOfWeights', 'sumsq', 'variance'
        ]
        # for python3 must cast dict.keys() to list
        got_keys = list(res2.keys())
        got_keys.sort()
        self.assertTrue(
            got_keys == expected_keys,
            'list of keys for doquantiles=False is incorrect'
        )
        
def suite():
    return[visstat_test]

if __name__ == '__main__':
    unittest.main()

#########################################################################
# test_task_hanningsmooth.py
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
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.hanningsmooth.html
#
##########################################################################
import os
import shutil
import unittest

from casatasks.private.parallel.parallel_data_helper import ParallelDataHelper
from casatasks import hanningsmooth, mstransform, partition, cvel, split, clearcal
from casatools import ctsys, ms
from casatestutils import testhelper as th

# Path for data
datapath = ctsys.resolve('unittest/hanningsmooth/')

'''
functional tests for task hanningsmooth
'''

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:   
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/hanningsmooth/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR

print('hanningsmooth tests will use data from %s' % datapath)

class test_base(unittest.TestCase):

    def setUp_ngc5921(self):
        # only DATA column
        self.msfile = 'ngc5921_ut.ms'
        if testmms:
            self.msfile = 'ngc5921_ut.mms'
            
        if (not os.path.exists(self.msfile)):
            shutil.copytree(os.path.join(datapath,self.msfile), self.msfile)

    def setUp_almams(self):
        # MS with DATA and CORRECTED_DATA
        self.msfile = 'ALMA-data-mst-science-testing-CAS-5013-one-baseline-one-timestamp.ms'
        if testmms:
            self.msfile = 'ALMA-data-mst-science-testing-CAS-5013-one-baseline-one-timestamp.mms'
            
        if (not os.path.exists(self.msfile)):
            shutil.copytree(os.path.join(datapath,self.msfile), self.msfile)

    def createMMS(self, msfile, column='data', axis='auto',scans='',spws=''):
        '''Create MMSs for tests with input MMS'''
        prefix = msfile.rstrip('.ms')
        if not os.path.exists(msfile):
            os.system('cp -RL '+os.path.join(datapath,msfile)+' '+ msfile)
        
        # Create an MMS for the tests
        self.testmms = prefix + ".test.mms"

        if os.path.exists(self.testmms):
            os.system("rm -rf " + self.testmms)
            
        print("................. Creating test MMS ..................")
        partition(vis=msfile, outputvis=self.testmms, datacolumn=column,
                    createmms=True,separationaxis=axis, scan=scans, spw=spws)


class hanningsmooth_test1(test_base):
    
    def setUp(self):
        self.setUp_ngc5921()

    def tearDown(self):
        if (os.path.exists(self.outputms)):
            shutil.rmtree(self.outputms,ignore_errors=True)        
        
    def test1(self):
        """hanningsmooth - Test 1: Wrong input MS should raise an exception"""
        msfile = 'badmsfile'
        self.outputms = 'none.ms'
        passes = False
        # CASA5 tasks do not throw exceptions, they return a value of False
        try:
            ret = hanningsmooth(vis=msfile)
            if not ret:
                passes = True
        except:
            passes = True
        self.assertTrue(passes)
        
    def test2(self):
        '''hanningsmooth - Test 2: Check that output MS is created'''
        self.outputms = 'hann2.ms'
        hanningsmooth(vis=self.msfile, outputvis=self.outputms, datacolumn='corrected')
        # Smoothed data should be saved in DATA column of outupt MS
        self.assertTrue(os.path.exists(self.outputms))
                
    def test3(self):
        '''hanningsmooth - Test 3: Check theoretical and calculated values on non-existing CORRECTED column'''
        self.outputms = 'hann3.ms'

      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.msfile, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [False])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [False])

        # It should fall-back and use the input DATA column
        hanningsmooth(vis=self.msfile, outputvis=self.outputms, datacolumn='corrected')

      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.outputms, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [True])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [True])

        data_col = th.getVarCol(self.msfile, 'DATA')
        corr_col = th.getVarCol(self.outputms, 'DATA')
        nrows = len(corr_col)
        
      # Loop over every 2nd row,pol and get the data for each channel
        max = 1e-05
        for i in range(1,nrows,2) :
            row = 'r%s'%i            
            # polarization is 0-1
            for pol in range(0,2) :
                # array's channels is 0-63
                for chan in range(1,62) :
                    # channels must start from second and end before the last
                    data = data_col[row][pol][chan]
                    dataB = data_col[row][pol][chan-1]
                    dataA = data_col[row][pol][chan+1]
        
                    Smoothed = th.calculateHanning(dataB,data,dataA)
                    CorData = corr_col[row][pol][chan]
                    
                    # Check the difference
                    self.assertTrue(abs(CorData-Smoothed) < max )

    def test4(self):
        '''hanningsmooth - Test 4: Theoretical and calculated values should be the same for MMS-case'''
    
        # Split the input to decrease the running time
        split(self.msfile, outputvis='splithan.ms',scan='1,2',datacolumn='data')
        self.msfile = 'splithan.ms'
        
        # create a test MMS. It creates self.testmms
        self.createMMS(self.msfile)
        self.outputms = 'hann4.mms'
        
      # check correct flagging (just for one row as a sample)
        mslocal = ms()
        mslocal.open(self.msfile)
        mslocal.sort('sorted.ms',['OBSERVATION_ID','ARRAY_ID','SCAN_NUMBER','FIELD_ID','DATA_DESC_ID','ANTENNA1','ANTENNA2','TIME'])
        mslocal.close()
        self.msfile = 'sorted.ms'
        flag_col = th.getVarCol(self.msfile, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [False])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [False])
        
        data_col = th.getVarCol(self.msfile, 'DATA')        
        hanningsmooth(vis=self.testmms, outputvis=self.outputms, datacolumn='data', keepmms=True)
        self.assertTrue(ParallelDataHelper.isParallelMS(self.outputms), 'Output should be an MMS')

      # Sort the MMS
        mslocal.open(self.outputms)
        mslocal.sort('sorted.mms',['OBSERVATION_ID','ARRAY_ID','SCAN_NUMBER','FIELD_ID','DATA_DESC_ID','ANTENNA1','ANTENNA2','TIME'])
        mslocal.close()
        self.outputms = 'sorted.mms'
        
        corr_col = th.getVarCol(self.outputms, 'DATA')
        nrows = len(corr_col)

      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.outputms, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [True])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [True])
        
      # Loop over every 2nd row,pol and get the data for each channel
        max = 1e-05
        for i in range(1,nrows,2) :
            row = 'r%s'%i            
            # polarization is 0-1
            for pol in range(0,2) :
                # array's channels is 0-63
                for chan in range(1,62) :
                    # channels must start from second and end before the last
                    data = data_col[row][pol][chan]
                    dataB = data_col[row][pol][chan-1]
                    dataA = data_col[row][pol][chan+1]
        
                    Smoothed = th.calculateHanning(dataB,data,dataA)
                    CorData = corr_col[row][pol][chan]
                    
                    # Check the difference
                    self.assertTrue(abs(CorData-Smoothed) < max )

    def test6(self):
        '''hanningsmooth - Test 6: Flagging should be correct with datacolumn==ALL'''
        self.outputms = 'hann6.ms'
        
      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.msfile, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [False])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [False])

        hanningsmooth(vis=self.msfile,outputvis=self.outputms, datacolumn='all')

      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.outputms, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [True])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [True])

    def test7(self):
        '''hanningsmooth - Test 7: Flagging should be correct when hanning smoothing within cvel (no transform)'''
        self.outputms = 'cvelngc.ms'
        clearcal(vis=self.msfile)
        
      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.msfile, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [False])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [False])

        cvel(vis=self.msfile, outputvis=self.outputms, hanning=True)

      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.outputms, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [True])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [True])

    def test8(self):
        '''hanningsmooth - Test 8: Flagging should be correct when hanning smoothing within mstransform (with regrid)'''
        self.outputms = 'cvelngc.ms'
        
      # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.msfile, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [False])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [False])
        self.assertTrue(flag_col['r1'][0][62] == [False])

        # CAS-4114 cvel doesn't support MMS. Compare with mstransform instead.
#        cvel(vis=self.msfile, outputvis=self.outputms, hanning=True, outframe='cmb')
        mstransform(vis=self.msfile, outputvis=self.outputms, datacolumn='data',
                    hanning=True, regridms=True, outframe='cmb')

        # check correct flagging (just for one row as a sample)
        flag_col = th.getVarCol(self.outputms, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [True])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][2] == [False])
        self.assertTrue(flag_col['r1'][0][60] == [False])
        self.assertTrue(flag_col['r1'][0][61] == [True])
        self.assertTrue(flag_col['r1'][0][62] == [True])


class hanningsmooth_test2(test_base):
    
    def setUp(self):
        self.setUp_almams()

    def tearDown(self):
        if (os.path.exists(self.outputms)):
            shutil.rmtree(self.outputms,ignore_errors=True)       
             
    def test_default_cols(self):
        '''hanningsmooth: Default datacolumn=all and MMS output'''
        
        self.createMMS(self.msfile,column='all')
        self.outputms = 'hannall.ms'

        hanningsmooth(vis=self.testmms, outputvis=self.outputms)
        self.assertTrue(ParallelDataHelper.isParallelMS(self.outputms), 'Output should be an MMS')
        
        # Should have all scratch columns in output
        cd = th.getColDesc(self.outputms, 'DATA')
        self.assertGreater(len(cd), 0, 'DATA column does not exist')
        cc = th.getColDesc(self.outputms, 'CORRECTED_DATA')
        self.assertGreater(len(cc), 0, 'CORRECTED_DATA does not exist')
        
        # Now repeat the above steps but create an output MS by setting keepmms=False
        os.system('rm -rf '+self.outputms)
        hanningsmooth(vis=self.testmms, outputvis=self.outputms, keepmms=False)
        self.assertFalse(ParallelDataHelper.isParallelMS(self.outputms), 'Output should be a normal MS')
        
        # Should have all scratch columns in output
        cd = th.getColDesc(self.outputms, 'DATA')
        self.assertGreater(len(cd), 0, 'DATA column does not exist')
        cc = th.getColDesc(self.outputms, 'CORRECTED_DATA')
        self.assertGreater(len(cc), 0, 'CORRECTED_DATA does not exist')

    def test_corrected_col(self):
        '''hanningsmooth: Apply smoothing in CORRECTED column'''
        self.outputms = 'hanncorr.ms'

      # check correct flagging before (just for one row as a sample)
        flag_col = th.getVarCol(self.msfile, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [False])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][3838] == [False])
        self.assertTrue(flag_col['r1'][0][3839] == [False])
        
        # input column
        data_col = th.getVarCol(self.msfile, 'CORRECTED_DATA') 
               
        hanningsmooth(vis=self.msfile, outputvis=self.outputms, datacolumn='corrected')
        
        # output smoothed column
        corr_col = th.getVarCol(self.outputms, 'DATA')
        nrows = len(corr_col)
        
      # check correct flagging after (just for one row as a sample)
        flag_col = th.getVarCol(self.outputms, 'FLAG')
        self.assertTrue(flag_col['r1'][0][0] == [True])
        self.assertTrue(flag_col['r1'][0][1] == [False])
        self.assertTrue(flag_col['r1'][0][3838] == [False])
        self.assertTrue(flag_col['r1'][0][3839] == [True])

      # Loop over every 2nd row,pol and get the data for each channel
        max = 1e-04
        for i in range(1,nrows,2) :
            row = 'r%s'%i            
            # polarization is 0-1
            for pol in range(0,2) :
                # array's channels is 0-3840
                for chan in range(1,3839) :
                    # channels must start from second and end before the last
                    data = data_col[row][pol][chan]
                    dataB = data_col[row][pol][chan-1]
                    dataA = data_col[row][pol][chan+1]
        
                    Smoothed = th.calculateHanning(dataB,data,dataA)
                    CorData = corr_col[row][pol][chan]
                    
                    # Check the difference
                    self.assertTrue(abs(CorData-Smoothed) < max, 
                                    'CorData=%s Smoothed=%s in row=%s pol=%s chan=%s'%(CorData,Smoothed,row,pol,chan))


class Cleanup(test_base):

    def tearDown(self):
        shutil.rmtree('ngc5921_ut.ms', ignore_errors=True)
        shutil.rmtree('ALMA-data-mst-science-testing-CAS-5013-one-baseline-one-timestamp.ms', ignore_errors=True)
        
    def test_runTest(self):
        '''hanningsmooth: Cleanup'''
        pass
            
if __name__ == '__main__':
    unittest.main()

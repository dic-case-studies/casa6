##########################################################################
# test_req_task_gaincal.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_gaincal/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import gaincal, casalog
    CASA6 = True
    tb = casatools.table()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import os
import unittest
import shutil
import numpy as np
import pylab as pl

if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/vla/gaincaltest2.ms')
    compCal = casatools.ctsys.resolve('caltables/gaincaltest2.ms.G0')
    tCal = casatools.ctsys.resolve('caltables/gaincaltest2.ms.T0')
    # Reference Cals
    combinedRef = casatools.ctsys.resolve('caltables/genDataCombine.G0')
    preTRef = casatools.ctsys.resolve('caltables/genDataPreT.G0')
    preGRef = casatools.ctsys.resolve('caltables/genDataPreG.T0')
    calModeP = casatools.ctsys.resolve('caltables/calModeTest.G0')
    calModeA = casatools.ctsys.resolve('caltables/calModeTest.G1')
    typeCalK = casatools.ctsys.resolve('caltables/gaintypek.G0')
    typeCalSpline = casatools.ctsys.resolve('caltables/gaintypeSpline.G0')
    spwMapCal = casatools.ctsys.resolve('caltables/spwMap.G0')
    
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/gaincaltest2.ms'
        compCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/gaincaltest2.ms.G0'
        tCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/gaincaltest2.ms.T0'
        
        combinedRef = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/genDataCombine.G0'
        preTRef = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/genDataPreT.G0'
        preGRef = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/genDataPreG.T0'
        calModeP = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/calModeTest.G0'
        calModeA = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/calModeTest.G1'
        typeCalK = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/gaintypek.G0'
        typeCalSpline = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/gaintypeSpline.G0'
        spwMapCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/spwMap.G0'
        
        
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/gaincaltest2.ms'
        compCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/gaincaltest2.ms.G0'
        tCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/gaincaltest2.ms.T0'
        
        combinedRef = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/genDataCombine.G0'
        preTRef = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/genDataPreT.G0'
        preGRef = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/genDataPreG.T0'
        calModeP = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/calModeTest.G0'
        calModeA = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/calModeTest.G1'
        typeCalK = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/gaintypek.G0'
        typeCalSpline = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/gaintypeSpline.G0'
        spwMapCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/spwMap.G0'
        
fullRangeCal = 'testgaincal.cal'
maxScanCal = 'testScan.cal'
int70Cal = 'int70.cal'
int30Cal = 'int30.cal'

tempCal = 'temp.cal'
tempCal2 = 'temp2.cal'
selectCal = 'select.cal'


def getparam(caltable, colname='CPARAM'):
    ''' Open a caltable and get the provided column '''

    tb.open(caltable)
    outtable = tb.getcol(colname)
    tb.close()

    return outtable


def tableComp(table1, table2, cols=[], rtol=8e-7, atol=1e-8):
    ''' Compare two caltables '''

    tableVal1 = {}
    tableVal2 = {}

    tb.open(table1)
    colname1 = tb.colnames()

    for col in colname1:
        try:
            tableVal1[col] = tb.getcol(col)
        except RuntimeError:
            pass
    tb.close()

    tb.open(table2)
    colname2 = tb.colnames()

    for col in colname2:
        try:
            tableVal2[col] = tb.getcol(col)
        except RuntimeError:
            pass
    tb.close()

    truthDict = {}

    for col in tableVal1.keys():

        try:
            truthDict[col] = np.isclose(tableVal1[col], tableVal2[col], rtol=rtol, atol=atol)
        except TypeError:
            print(col, 'ERROR in finding truth value')
            casalog.post(message=col+': ERROR in determining the truth value')

    if len(cols) == 0:
        
        truths = [[x, np.all(truthDict[x] == True)] for x in truthDict.keys()]

    else:

        truths = [[x, np.all(truthDict[x] == True)] for x in cols]

    return np.array(truths)


class gaincal_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        gaincal(vis=datapath, caltable=fullRangeCal, combine='scan', solint='inf', field='0', refant='0',
                smodel=[1, 0, 0, 0], scan='0~9')

        gaincal(vis=datapath, caltable=maxScanCal, solint='inf', field='0', refant='0',
                smodel=[1, 0, 0, 0], scan='0~9')

        gaincal(vis=datapath, caltable=int70Cal, solint='70s', field='0', refant='0',
                smodel=[1, 0, 0, 0], scan='0~9')

        gaincal(vis=datapath, caltable=int30Cal, solint='30s', field='0', refant='0',
                smodel=[1, 0, 0, 0], scan='0~9')
        
        gaincal(vis=datapath, caltable=selectCal, solint='inf', field='0', refant='0',
                smodel=[1, 0, 0, 0], scan='2', spw='2')

    def tearDown(self):
        if os.path.exists(tempCal):
            shutil.rmtree(tempCal)
            
        if os.path.exists(tempCal2):
            shutil.rmtree(tempCal2)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(fullRangeCal):
            shutil.rmtree(fullRangeCal)

        if os.path.exists(maxScanCal):
            shutil.rmtree(maxScanCal)

        if os.path.exists(int70Cal):
            shutil.rmtree(int70Cal)

        if os.path.exists(int30Cal):
            shutil.rmtree(int30Cal)
            
        if os.path.exists(selectCal):
            shutil.rmtree(selectCal)

    def test_correctGains(self):
        '''
            test_correctGains
            -------------------
            
            Check that the gaincal results match a reference gaincal table
        '''

        self.assertTrue(np.all(tableComp(fullRangeCal, combinedRef)[:,1] == 'True'))


    def test_intervalSNR(self):
        '''
            test_intervalSNR
            ------------------
            
            Check that shorter interval times result in a lower signal to noise
        '''

        snrCombine = np.mean(getparam(fullRangeCal, 'SNR'))
        snrScans = np.mean(getparam(maxScanCal, 'SNR'))
        int70Snr = np.mean(getparam(int70Cal, 'SNR'))
        int30Snr = np.mean(getparam(int30Cal, 'SNR'))

        self.assertTrue(int30Snr < int70Snr < snrScans < snrCombine)

    def test_minSNR(self):
        '''
            test_minSNR
            -------------
            
            Check that values below the provided SNR threshold are flagged
        '''

        gaincal(vis=datapath, caltable=tempCal, solint='30s', field='0', refant='0',
                smodel=[1, 0, 0, 0], minsnr=1000)

        flagged = getparam(tempCal, 'FLAG')

        self.assertTrue(np.all(flagged == 1))

    def test_fieldSelect(self):
        '''
            test_fieldSelect
            ------------------
            
            Check that the field selection parameter functions properly
        '''

        fields = getparam(fullRangeCal, 'FIELD_ID')

        self.assertTrue(np.all(fields == 0))

    def test_refantSelect(self):
        '''
            test_refantSelect
            -------------------
            
            Check that the refant selection functions properly
        '''

        refants = getparam(fullRangeCal, 'ANTENNA2')

        self.assertTrue(np.all(refants == 0))

    def test_scanSelect(self):
        '''
            test_scanSelect
            -----------------
            
            Check that the scan selection functions properly
        '''

        scans = getparam(selectCal, 'SCAN_NUMBER')

        self.assertTrue(np.all(scans == 2))

    def test_spwSelect(self):
        '''
            test_spwSelect
            ----------------
            
            Check that the spw selection parameter functions properly
        '''

        spws = getparam(selectCal, 'SPECTRAL_WINDOW_ID')

        self.assertTrue(np.all(spws == 2))

    def test_refantDiff(self):
        '''
            test_refantDiff
            -----------------
            
            Check that selecting refant will cause that refant to be set at 0
        '''

        gaincal(vis=datapath, caltable=tempCal, solint='inf', field='0', combine='scan', refant='1',
                smodel=[1, 0, 0, 0])

        gAmp = getparam(tempCal)
        refs = [[np.mean(gAmp.imag[j,0,i::10]) for i in range(10)] for j in range(2)]

        self.assertTrue(refs[0][1] == 0 and refs[1][1] == 0)

    def test_preapplyT0(self):
        '''
            test_preapplyT0
            -----------------
            
            Check that pre applying the T table results in the regular G table calibration
        '''

        gaincal(vis=datapath, caltable=tempCal, refant='0', solint='inf', smodel=[1, 0, 0, 0], gaintype='G', field='0', combine='scan',
                  gaintable=[tCal])
        
        self.assertTrue(np.all(tableComp(preTRef, tempCal)[:,1] == 'True'))

    def test_preapplyG0(self):
        '''
            test_preapplyG0
            -----------------
            
            Check that pre applying the G table results in the regular T table calibration
        '''

        gaincal(vis=datapath, caltable=tempCal, refant='0', solint='int', smodel=[1, 0, 0, 0],
                gaintype='T', gaintable=[compCal])
        
        self.assertTrue(np.all(tableComp(preGRef, tempCal)[:,1] == 'True'))
        
    def test_antennaSelect(self):
        '''
            test_antennaSelect
            --------------------
            
            Check that antennas that aren't selected are flagged
        '''
        
        gaincal(vis=datapath, caltable=tempCal, refant='0', field='0', solint='inf', combine='scan', antenna='0~5&', smodel=[1,0,0,0], gaintype='G')
        
        flags = getparam(tempCal, colname='FLAG')
        flagnum = np.sum(flags)
        
        self.assertTrue(flagnum == 32)
        
    def test_minBl(self):
        '''
            test_minBl
            ------------
            
            Check that if the min baseline threshold isn't met those antennas aren't used. If no antennas have enough then a file is not written.
        '''
        
        gaincal(vis=datapath, caltable=tempCal, refant='0', solint='int', smodel=[1,0,0,0], gaintype='G', combine='scan', antenna='0~5&', minblperant=6)
        
        self.assertFalse(os.path.exists(tempCal))
        
    def test_preboth(self):
        '''
            test_preboth
            --------------
            
            Check that when preapplying both then gaintype T increases the SNR
            
        '''
        
        gaincal(vis=datapath, caltable=tempCal, refant='0', solint='inf', smodel=[1, 0, 0, 0], gaintype='G', field='0',
                  gaintable=[tCal, compCal], gainfield=['0','0'], interp=[''])
        
        gaincal(vis=datapath, caltable=tempCal2, refant='0', solint='inf', smodel=[1, 0, 0, 0], gaintype='T', field='0',
                  gaintable=[tCal, compCal])
        
        SNR1 = np.mean(getparam(tempCal, colname='SNR'))
        SNR2 = np.mean(getparam(tempCal2, colname='SNR'))
        
        self.assertTrue(SNR1 < SNR2)
        
    def test_calModeP(self):
        '''
            test_calModeP
            ---------------
            
            Check that the output with calmode 'p' is equal to a reference calibration table
        '''
        
        gaincal(vis=datapath, caltable=tempCal, field='0', smodel=[1,0,0,0], solint='inf', combine='scan', calmode='p')
        
        self.assertTrue(np.all(tableComp(tempCal, calModeP)[:,1] == 'True'))
        
    def test_calModeA(self):
        '''
            test_calModeA
            ---------------
            
            Check that the output with calmode 'a' is equal to a reference calibration table
        '''
        
        gaincal(vis=datapath, caltable=tempCal, field='0', smodel=[1,0,0,0], solint='inf', combine='scan', calmode='a')
        
        self.assertTrue(np.all(tableComp(tempCal, calModeA)[:,1] == 'True'))
        
    def test_gainTypeK(self):
        '''
            test_gainTypeK
            ----------------
            
            Check that the output with gaintype k is equal to a reference calibration table
        '''
        
        gaincal(vis=datapath, caltable=tempCal, field='0', smodel=[1,0,0,0], solint='inf', combine='scan', gaintype='KCROSS', refant='0')
        
        self.assertTrue(np.all(tableComp(tempCal, typeCalK)[:,1] == 'True'))
        
    def test_gainTypeSpline(self):
        '''
            test_gainTypeK
            ----------------
            
            Check that the output with gaintype GSPLINE is equal to a reference calibration table
        '''
        
        gaincal(vis=datapath, caltable=tempCal, field='0', smodel=[1,0,0,0], solint='inf', combine='scan', gaintype='GSPLINE', refant='0')
        
        self.assertTrue(np.all(tableComp(tempCal, typeCalSpline)[:,1] == 'True'))
        
    def test_spwMap(self):
        '''
            test_spwMap
            -------------
            
            Check that the output with spwMap matches to a reference calibration table
        '''
        
        gaincal(vis=datapath, caltable=tempCal, field='0', smodel=[1,0,0,0], solint='inf', combine='scan', refant='0',spwmap=[0,0,1,1])
        
        self.assertTrue(np.all(tableComp(tempCal, spwMapCal)[:,1] == 'True'))
        


def suite():
    return[gaincal_test]


if __name__ == '__main__':
    unittest.main()

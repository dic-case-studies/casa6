##########################################################################
# test_req_task_applycal.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_applycal/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import applycal
    CASA6 = True
    tb = casatools.table()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy as np


if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/vla/gaincaltest2.ms')
    tCal = casatools.ctsys.resolve('caltables/gaincaltest2.ms.T0')
    gCal = casatools.ctsys.resolve('caltables/gaincaltest2.ms.G0')
    
    altdata = casatools.ctsys.resolve('visibilities/vla/ngc5921.ms/')
    altCal = casatools.ctsys.resolve('caltables/ngc5921.ref1a.gcal')
    callibfile = casatools.ctsys.resolve('text/refcallib.txt')
    
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/gaincaltest2.ms/'
        tCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/gaincaltest2.ms.T0'
        gCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/gaincaltest2.ms.G0'
        
        altdata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/ngc5921.ms/'
        altCal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ngc5921.ref1a.gcal/'
        callibfile = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/text/refcallib.txt'
        
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/gaincaltest2.ms/'
        tCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/gaincaltest2.ms.T0'
        gCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/gaincaltest2.ms.G0'
        
        altdata = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/ngc5921.ms/'
        altCal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ngc5921.ref1a.gcal/'
        callibfile = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/text/refcallib.txt'
    
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

datacopy = 'applycalcopy.ms'
dataref = 'referencedata.ms'
altcopy = 'altcopy.ms'
temptcal = 'temptcal.T0'

class applycal_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def setUp(self):
        shutil.copytree(datapath, datacopy)
        shutil.copytree(datapath, dataref)
    
    def tearDown(self):
        shutil.rmtree(datacopy)
        if os.path.exists('applycalcopy.ms.flagversions'):
            shutil.rmtree('applycalcopy.ms.flagversions')
            
        if os.path.exists(dataref):
            shutil.rmtree(dataref)
            
        if os.path.exists('referencedata.ms.flagversions'):
            shutil.rmtree('referencedata.ms.flagversions')
            
        if os.path.exists(altcopy):
            shutil.rmtree(altcopy)
            
        if os.path.exists('altcopy.ms.flagversions'):
            shutil.rmtree('altcopy.ms.flagversions')
            
        if os.path.exists(temptcal):
            shutil.rmtree(temptcal)

        
    @classmethod
    def tearDownClass(cls):
        pass
    
    def test_corrected(self):
        '''
            test_corrected
            ----------------
            
            Check that the CORRECTED_DATA column is added to the MS
        '''
        
        applycal(vis=datacopy, gaintable=[gCal])
        try:
            corrected = getparam(datacopy, 'CORRECTED_DATA')
            colmade = True
        except RuntimeError:
            colmade = False
            
        self.assertTrue(colmade, msg='No CORRECTED_DATA colun was made')
        
        
    def test_overwritten(self):
        '''
            test_overwritten
            ------------------
            
            Check that the CORRECTED_DATA column is overwritten with consecutive runs of applycal
        '''
        
        applycal(vis=datacopy, gaintable=[gCal])
        correctedFirst = getparam(datacopy, 'CORRECTED_DATA')
        applycal(vis=datacopy, gaintable=[tCal])
        correctedSecond = getparam(datacopy, 'CORRECTED_DATA')
        
        self.assertFalse(np.all(correctedFirst == correctedSecond))
        
    def test_fieldSelect(self):
        '''
            test_fieldSelect
            ------------
            
            Check that fields are properly selected
        '''
        # TODO: Need a dataset with more fields
        
        applycal(vis=datacopy, gaintable=[gCal], field='1')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.33317375198725724+0.028998389148771512j)
        
    def test_spwSelect(self):
        '''
            test_spwSelect
            ----------------
            
            Check that spws are properly selected
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], spw='1')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.28994172236678006+0.03374482270467842j)
        
    def test_intentSelect(self):
        '''
            test_intentSelect
            -------------------
            
            Check that intents are properly selected
        '''

        applycal(vis=datacopy, gaintable=[gCal], intent='*AMPLI*')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.40933912826626284+0.010524099430369143j)
        
    def test_timerangeSelect(self):
        '''
            test_timerangeSelect
            ----------------------
            
            Check that the timerange parameter properly selects a subset of the data
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], timerange='04:33:23~04:38:23')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.26859135173414944+0.043249139614812186j)
        
    def test_uvrangeSelect(self):
        '''
            test_uvrangeSelect
            --------------------
            
            Check that the uvrange parameter properly selects a subset of the data
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], uvrange='>100klambda')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.40898735923515789+0.010907960569579565j)
        
    def test_antennaSelect(self):
        '''
            test_antennaSelect
            --------------------
            
            Check that the antenna parameter properly selects a subset of the data
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], antenna='0~5&')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.28547721979290314+0.024605255460554348j)
        
    def test_scanSelect(self):
        '''
            test_scanSelect
            -----------------
            
            Check that the scan parameter properly selects a subset of the data
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], scan='10')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.26857740451349088+0.04325810174871688j)
        
    def test_callib(self):
        '''
            test_callib
            -------------
            
            Check that a callib file is taken and the calibration tables are pre-applied
        '''
        shutil.copytree(tCal, temptcal)
        
        applycal(vis=datacopy, docallib=True, callib=callibfile)
        #applycal(vis=dataref, gaintable=[tCal])
        
        #noSelect = getparam(dataref, 'CORRECTED_DATA')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.31622877260041754+0.049053979716640904j)
        
    def test_gaintable(self):
        '''
            test_gaintable
            ----------------
            
            Check that the gaintable parameter selects caltables to be applied
        '''
        
        applycal(vis=datacopy, gaintable=[gCal,tCal])
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.50000433754678164+3.7223202622521605e-05j)
        
    def test_gainfield(self):
        '''
            test_gainfield
            ----------------
            
            Check that the gainfield parameter selects a subset of the gaintables
        '''
        
        applycal(vis=datacopy, gaintable=[tCal], gainfield='1')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.31430447832422664+0.048264338579282674j)
        
    def test_interp(self):
        '''
            test_interp
            -------------
            
            Check that the interp parameter changes the interpolation type for each gaintable
        '''
        
        applycal(vis=datacopy, gaintable=[tCal], interp='cubicflag')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.31622877260293381+0.049053979715317733j)
        
    def test_spwmap(self):
        '''
            test_spwmap
            -------------
            
            Check that the spwmap parameter applys a subset of the gaintable(s)
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], spwmap=[0,0,1,1])
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.39223749429011501+0.033263095996865069j)
        
    def test_calwt(self):
        '''
            test_calwt
            ------------
            
            Check that calwt = True calibrates the data weights per gaintable
        '''
        
        #NOTE All boolear parameter values must be in an array/list?
        
        applycal(vis=datacopy, gaintable=[gCal], calwt=[False])
        datamean = np.mean(getparam(datacopy, 'WEIGHT'))
        
        self.assertTrue(datamean == 576072000.0)
        
    def test_parang(self):
        '''
            test_parang
            -------------
            
            Check that the parallactic angle correction is applied if parang = True
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], parang=True)
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.40936784517874114+0.010518063253447309j)
        
    def test_calflag(self):
        '''
            test_calflag
            ----------------
            
            Check that calflag calibrated the data and the flags
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], applymode='calflag')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.40933912826626284+0.010524099430369143j)
        
    def test_trial(self):
        '''
            test_trial
            ------------
            
            Check that the trial mode leaves the dataset unchanged
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], applymode='trial')
        data = getparam(datacopy, 'DATA')
        corrected = getparam(datacopy, 'CORRECTED_DATA')
        
        self.assertTrue(np.all(data == corrected))
    
    def test_calflagstrict(self):
        ''' 
            test_calflagstrict
            --------------------
            
            Check that all selected data that have no solutions will be flagged
        '''
        
        #TODO: Come back to this, same as non strict opion right now
        
        applycal(vis=datacopy, gaintable=[gCal], applymode='calflagstrict')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.40933912826626284+0.010524099430369143j)
        
        
    def test_flagonly(self):
        '''
            test_flagonly
            ---------------
            
            Check that flags are applied, but not the calibration itself
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], applymode='flagonly')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.25983810264064156+0.045579181018852881j)
        
        
    def test_flagonlystrict(self):
        '''
            test_flagonlystrict
            ---------------------
            
            Check that all the selected data that have no solutions will be flagged
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], applymode='flagonlystrict')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.25983810264064156+0.045579181018852881j)
        
        
    def test_calonly(self):
        '''
            test_calonly
            --------------
            
            Check that only calibration and weights are applied, not flags
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], applymode='calonly')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        
        self.assertTrue(datamean == 0.40933912826626284+0.010524099430369143j)
        
        
    def test_flagbackup(self):
        '''
            test_flagbackup
            -----------------
            
            Check to see that a backup for the flags was made
        '''
        
        applycal(vis=datacopy, gaintable=[gCal], flagbackup=False)
        flagpath = datacopy + '.flagversions'
        
        self.assertFalse(os.path.exists(flagpath))
        
    
    
def suite():
    return[applycal_test]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--strict', action="store_true", default=False, dest='strict_flag')
    parser.add_argument('unittest_args', nargs='*')

    args = parser.parse_args()
    strict = False
    if args.strict_flag:
        strict = True
    sys.argv[1:] = args.unittest_args
    unittest.main() 

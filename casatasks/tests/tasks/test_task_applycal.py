##########################################################################
# test_task_applycal.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.calibration.applycal.html
#
#
##########################################################################
import sys
import os
import unittest
import shutil
import numpy as np

import casatools
from casatasks import applycal, mstransform, gaincal, casalog, clearcal
from casatasks.private.callibrary import callibrary
tb = casatools.table()
calibrater = casatools.calibrater

datapath = casatools.ctsys.resolve('unittest/applycal/')

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

def change_perms(path):
    os.chmod(path, 0o777)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root,d), 0o777)
        for f in files:
            os.chmod(os.path.join(root,f), 0o777)


def drop_solution_for_antenna(caltable, antenna):
    """Drop calibration solution for given antenna

    Arguments:
        caltable {string} -- caltable name
        antenna {int} -- antenna ID
    """
    tb.open(caltable, nomodify=False)
    try:
        tsel = tb.query('ANTENNA1={}'.format(antenna))
        rows = tsel.rownumbers()
        tsel.close()
        tb.removerows(rows)
    finally:
        tb.close()

# Input data files
msdata = 'gaincaltest2.ms'
altdata = 'ngc5921.ms'
tCal = 'gaincaltest2.ms.T0'
gCal = 'gaincaltest2.ms.G0'
altCal = 'ngc5921.ref1a.gcal'
callibfile = 'refcallib.txt'
mmsdata = 'ngc5921.applycal.mms'
mmsbcal = 'ngc5921.bcal'
mmsgcal = 'ngc5921.gcal'
mmsflux = 'ngc5921.fluxscale'
vlbadata = 'ba123a.ms'
vlbaGCCal = 'ba123a.gc'


# Locally copied data files
datacopy = 'applycalcopy.ms'
dataref = 'referencedata.ms'
altcopy = 'altcopy.ms'
temptcal = 'temptcal.T0'

mmscopy = 'mmsapplycalcopy.mms'
mmsbcalcopy = 'mmsbcalcopy.cal'
mmsgcalcopy = 'mmsgcalcopy.cal'
mmsfluxcopy = 'mmsfluxcopy.cal'

vlbacopy = 'vlbacopy.ms'

class applycal_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        shutil.copytree(os.path.join(datapath,vlbadata), vlbacopy)
        pass

    def setUp(self):
        shutil.copytree(os.path.join(datapath, msdata), datacopy)
        shutil.copytree(os.path.join(datapath, msdata), dataref)
        shutil.copytree(os.path.join(datapath, mmsdata), mmscopy)

        shutil.copytree(os.path.join(datapath, mmsbcal), mmsbcalcopy)
        shutil.copytree(os.path.join(datapath, mmsgcal), mmsgcalcopy)
        shutil.copytree(os.path.join(datapath, mmsflux), mmsfluxcopy)

        change_perms(mmsbcalcopy)
        change_perms(mmsgcalcopy)
        change_perms(mmsfluxcopy)
        
        # Link  cal tables to local directory
        if not os.path.exists(tCal):
            os.system('ln -s '+ datapath + tCal + ' .')
        if not os.path.exists(gCal):
            os.system('ln -s '+ datapath + gCal + ' .')
        if not os.path.exists(altCal):
            os.system('ln -s '+ datapath + altCal + ' .')
        if not os.path.exists(callibfile):
            os.system('ln -s '+ datapath + callibfile + ' .')
        if not os.path.exists(vlbaGCCal):
            os.system('ln -s '+ datapath + vlbaGCCal + ' .')

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

        if os.path.exists(mmscopy):
            shutil.rmtree(mmscopy)

        if os.path.exists(mmsbcalcopy):
            shutil.rmtree(mmsbcalcopy)

        if os.path.exists(mmsgcalcopy):
            shutil.rmtree(mmsgcalcopy)

        if os.path.exists(mmsfluxcopy):
            shutil.rmtree(mmsfluxcopy)

        if os.path.exists('cl_fldmap_test.ms.flagversions'):
            shutil.rmtree('cl_fldmap_test.ms.flagversions')
                
        if os.path.exists('mmsapplycalcopy.mms.flagversions'):
            shutil.rmtree('mmsapplycalcopy.mms.flagversions')
        
        
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(vlbacopy)
        shutil.rmtree('cl_fldmap_test.ms', ignore_errors=True)
        shutil.rmtree('cl_fldmap_test.Gf0', ignore_errors=True)
        shutil.rmtree('cl_fldmap_test.Gf01', ignore_errors=True)
        caldata = [tCal,gCal,altCal,callibfile,vlbaGCCal]
        for ff in caldata:
            if os.path.exists(ff):
                os.unlink(ff)
        outlist = ['callib_f0.txt','callib_f01_m0.txt','callib_f01_s0.txt','callib_f01.txt','callib_f01_s01.txt','callib_f01_m01.txt']
        for f in outlist:
            if os.path.exists(f):
                os.remove(f)
 
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

        self.assertTrue(np.isclose(datamean, 0.33317375198725724+0.028998389148771512j))

    def test_spwSelect(self):
        '''
            test_spwSelect
            ----------------

            Check that spws are properly selected
        '''

        applycal(vis=datacopy, gaintable=[gCal], spw='1')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.28994172236678006+0.03374482270467842j))

    def test_intentSelect(self):
        '''
            test_intentSelect
            -------------------

            Check that intents are properly selected
        '''

        applycal(vis=datacopy, gaintable=[gCal], intent='*AMPLI*')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.40933912826626284+0.010524099430369143j))

    def test_timerangeSelect(self):
        '''
            test_timerangeSelect
            ----------------------

            Check that the timerange parameter properly selects a subset of the data
        '''

        applycal(vis=datacopy, gaintable=[gCal], timerange='04:33:23~04:38:23')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.26859135173414944+0.043249139614812186j))

    def test_uvrangeSelect(self):
        '''
            test_uvrangeSelect
            --------------------

            Check that the uvrange parameter properly selects a subset of the data
        '''

        applycal(vis=datacopy, gaintable=[gCal], uvrange='>100klambda')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.40898735923515789+0.010907960569579565j))

    def test_antennaSelect(self):
        '''
            test_antennaSelect
            --------------------

            Check that the antenna parameter properly selects a subset of the data
        '''

        applycal(vis=datacopy, gaintable=[gCal], antenna='0~5&')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.28547721979290314+0.024605255460554348j))

    def test_scanSelect(self):
        '''
            test_scanSelect
            -----------------

            Check that the scan parameter properly selects a subset of the data
        '''

        applycal(vis=datacopy, gaintable=[gCal], scan='10')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.26857740451349088+0.04325810174871688j))

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

        self.assertTrue(np.isclose(datamean, 0.31622877260041754+0.049053979716640904j))
        
    def test_callib_field(self):
        '''
            test_callib_field
            -----------------
            
            Exercise fldmap-related maneuvers in cal library
        '''

        clfldmaptest='cl_fldmap_test.ms'

        # slice out just spw=0
        mstransform(vis=datacopy,outputvis=clfldmaptest,datacolumn='data',spw='0')

        # Create cal tables to apply in different ways below
        #  solutions from field=0 only
        Gf0='cl_fldmap_test.Gf0'
        gaincal(vis=clfldmaptest,caltable=Gf0,field='0',solint='inf',refant='0',smodel=[1,0,0,0])
        #  solutions from fields 0 & 1
        Gf01='cl_fldmap_test.Gf01'
        gaincal(vis=clfldmaptest,caltable=Gf01,field='0,1',solint='inf',refant='0',smodel=[1,0,0,0])

        # a cal library object to use variously below
        clib=callibrary()

        # traditional apply using field 0 -only table
        clearcal(vis=clfldmaptest)
        applycal(vis=clfldmaptest,field='0,1',gaintable=[Gf0],interp=['linear'],flagbackup=False)
        cdref=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]
        self.assertTrue(np.isclose(np.abs(np.mean(cdref)),1.00150795))

        # cal library apply using field 0 -only table, no cal selection
        clearcal(vis=clfldmaptest)
        clib.clear()
        f0='callib_f0.txt'
        clib.add(caltable=Gf0,tinterp='linear')
        clib.write(f0)
        applycal(vis=clfldmaptest,field='0,1',docallib=True,callib=f0,flagbackup=False)
        cd0=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cd0-cdref)),0.0))

        # cal library apply using field 0,1 table, *mapping* field 0 only from cal
        clearcal(vis=clfldmaptest)
        clib.clear()
        f01m0='callib_f01_m0.txt'
        clib.add(caltable=Gf01,tinterp='linear',fldmap=[0,0])
        clib.write(f01m0)
        applycal(vis=clfldmaptest,field='0,1',docallib=True,callib=f01m0,flagbackup=False)
        cd01m0=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cd01m0-cdref)),0.0))

        # cal library apply using field 0,1 table, *selecting* field 0 only from cal
        clearcal(vis=clfldmaptest)
        clib.clear()
        f01s0='callib_f01_s0.txt'
        clib.add(caltable=Gf01,tinterp='linear',fldmap='0')
        clib.write(f01s0)
        applycal(vis=clfldmaptest,field='0,1',docallib=True,callib=f01s0,flagbackup=False)
        cd01s0=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cd01s0-cdref)),0.0))


        # traditional apply using field 0,1 table
        clearcal(vis=clfldmaptest)
        applycal(vis=clfldmaptest,field='0,1',gaintable=[Gf01],interp=['linear'])
        cdref=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cdref)),1.00463621186))

        # cal library apply using field 0,1 -only table, no cal selection
        clearcal(vis=clfldmaptest)
        clib.clear()
        f01='callib_f01.txt'
        clib.add(caltable=Gf01,tinterp='linear')
        clib.write(f01)
        applycal(vis=clfldmaptest,field='0,1',docallib=True,callib=f01,flagbackup=False)
        cd01=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cd01-cdref)),0.0))

        # cal library apply using field 0,1 table, *selecting* field 0,1 explicitly
        clearcal(vis=clfldmaptest)
        clib.clear()
        f01s01='callib_f01_s01.txt'
        clib.add(caltable=Gf01,tinterp='linear',fldmap='0,1')
        clib.write(f01s01)
        applycal(vis=clfldmaptest,field='0,1',docallib=True,callib=f01s01,flagbackup=False)
        cd01s01=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cd01s01-cdref)),0.0))


        # traditional apply using field 0,1 table, w/ gainfield='nearest'
        clearcal(vis=clfldmaptest)
        applycal(vis=clfldmaptest,field='0,1',gaintable=[Gf01],gainfield=['nearest'],interp=['linear'])
        cdref=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cdref)),1.00319536))

        # cal library apply using field 0,1 table, *mapping* field 0,1 (self!) only from cal
        clearcal(vis=clfldmaptest)
        clib.clear()
        f01m01='callib_f01_m01.txt'
        clib.add(caltable=Gf01,tinterp='linear',fldmap=[0,1])
        clib.write(f01m01)
        applycal(vis=clfldmaptest,field='0,1',docallib=True,callib=f01m01,flagbackup=False)
        cd01m01=getparam(clfldmaptest,'CORRECTED_DATA')[0:4:3,:,:]  # p-hands only
        self.assertTrue(np.isclose(np.abs(np.mean(cd01m01-cdref)),0.0))

        
    def test_callib_missing_antenna0(self):
        '''
            test_callib_missing_antenna0
            -------------

            Check that a callib file is taken and the calibration tables are pre-applied
            Verify bug fix for CAS-12881
        '''
        shutil.copytree(tCal, temptcal)

        # flag solution for ANTENNA 0
        drop_solution_for_antenna(temptcal, 0)

        applycal(vis=datacopy, docallib=True, callib=callibfile)
        #applycal(vis=dataref, gaintable=[tCal])

        #noSelect = getparam(dataref, 'CORRECTED_DATA')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))
        #print(f'datamean = {datamean}')

        self.assertTrue(np.isclose(datamean, 0.31622877260041754+0.049053979716640904j))


    def test_gaintable(self):
        '''
            test_gaintable
            ----------------

            Check that the gaintable parameter selects caltables to be applied
        '''

        applycal(vis=datacopy, gaintable=[gCal,tCal])
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.50000433754678164+3.7223202622521605e-05j))

    def test_gainfield(self):
        '''
            test_gainfield
            ----------------

            Check that the gainfield parameter selects a subset of the gaintables
        '''

        applycal(vis=datacopy, gaintable=[tCal], gainfield='1')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.31430447832422664+0.048264338579282674j))

    def test_interp(self):
        '''
            test_interp
            -------------

            Check that the interp parameter changes the interpolation type for each gaintable
        '''

        applycal(vis=datacopy, gaintable=[tCal], interp='cubicflag')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.31622877260293381+0.049053979715317733j))

    def test_spwmap(self):
        '''
            test_spwmap
            -------------

            Check that the spwmap parameter applys a subset of the gaintable(s)
        '''

        applycal(vis=datacopy, gaintable=[gCal], spwmap=[0,0,1,1])
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.39223749429011501+0.033263095996865069j))

    def test_spwmapMulti(self):
        '''
            test_spwmapMulti
            -------------------

            Check the function of spwmap when provided a list of lists
        '''

        applycal(vis=datacopy, gaintable=[gCal,tCal], spwmap=[[0,0,1,1],[0,0,1,1]])
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.47988006537755368+0.028796507242145018j))

    def test_calwt(self):
        '''
            test_calwt
            ------------

            Check that calwt = True calibrates the data weights per gaintable
        '''

        #NOTE All boolear parameter values must be in an array/list?

        applycal(vis=datacopy, gaintable=[gCal], calwt=[False])
        datamean = np.mean(getparam(datacopy, 'WEIGHT'))

        self.assertTrue(np.isclose(datamean, 576072000.0))

    def test_parang(self):
        '''
            test_parang
            -------------

            Check that the parallactic angle correction is applied if parang = True
        '''

        applycal(vis=datacopy, gaintable=[gCal], parang=True)
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.40936784517874114+0.010518063253447309j))

    def test_calflag(self):
        '''
            test_calflag
            ----------------

            Check that calflag calibrated the data and the flags
        '''

        applycal(vis=datacopy, gaintable=[gCal], applymode='calflag')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.40933912826626284+0.010524099430369143j))

    def test_trial(self):
        '''
            test_trial
            ------------

            Check that the trial mode leaves the dataset unchanged
        '''

        applycal(vis=datacopy, gaintable=[gCal], applymode='trial')
        
        tb.open(datacopy)
        columns = tb.colnames()
        tb.close()
        
        self.assertFalse('CORRECTED_DATA' in columns)

    def test_calflagstrict(self):
        '''
            test_calflagstrict
            --------------------

            Check that all selected data that have no solutions will be flagged
        '''

        #TODO: Come back to this, same as non strict opion right now

        applycal(vis=datacopy, gaintable=[gCal], applymode='calflagstrict')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.40933912826626284+0.010524099430369143j))


    def test_flagonly(self):
        '''
            test_flagonly
            ---------------

            Check that flags are applied, but not the calibration itself
        '''

        # CORRECTED_DATA column doesn't exist if using flagonly
        applycal(vis=datacopy, gaintable=[gCal], applymode='flagonly')
        
        tb.open(datacopy)
        columns = tb.colnames()
        tb.close()
        
        self.assertFalse('CORRECTED_DATA' in columns)
        self.assertTrue(os.path.exists('applycalcopy.ms.flagversions'))


    def test_flagonlystrict(self):
        '''
            test_flagonlystrict
            ---------------------

            Check that all the selected data that have no solutions will be flagged
        '''

        applycal(vis=datacopy, gaintable=[gCal], applymode='flagonlystrict')
        
        tb.open(datacopy)
        columns = tb.colnames()
        tb.close()
        
        self.assertFalse('CORRECTED_DATA' in columns)
        self.assertTrue(os.path.exists('applycalcopy.ms.flagversions'))


    def test_calonly(self):
        '''
            test_calonly
            --------------

            Check that only calibration and weights are applied, not flags
        '''

        applycal(vis=datacopy, gaintable=[gCal], applymode='calonly')
        datamean = np.mean(getparam(datacopy, 'CORRECTED_DATA'))

        self.assertTrue(np.isclose(datamean, 0.40933912826626284+0.010524099430369143j))


    def test_flagbackup(self):
        '''
            test_flagbackup
            -----------------

            Check to see that a backup for the flags was made
        '''

        applycal(vis=datacopy, gaintable=[gCal], flagbackup=True)
        flagpath = datacopy + '.flagversions'
        files = os.listdir(flagpath)
        result = any('flags.applycal' in i for i in files)

        self.assertTrue(result)

    ### Merged Test

    def test1_applycal_fluxscale_gcal_bcal(self):
        """Test: Apply calibration using fluxscal gcal and bcal tables. Create flagbackup for an MMS"""

        aux = [mmsbcalcopy, mmsgcalcopy, mmsfluxcopy]

        # Repository caltables are pre-v4.1, and we
        # must update them _before_ applycal to avoid contention
        casalog.post("Updating pre-v4.1 caltables: %s" % str(aux),"WARN","test1_applycal_fluxscale_gcal_bcal")
        cblocal = calibrater()
        for oldct in aux:
            cblocal.updatecaltable(oldct)
        casalog.post("Pre-v4.1 caltables updated","INFO","test1_applycal_fluxscale_gcal_bcal")

        # Run applycal in MMS mode. Verify that the flagbackup is correctly created for the top-level MMS only
        applycal(vis=mmscopy,field='',spw='',selectdata=False,gaintable=aux,
                 gainfield=['nearest','nearest','0'],
                 interp=['linear', 'linear','nearest'],spwmap=[], flagbackup=True)

        # Verify that flagbackup works
        self.assertTrue(os.path.exists(mmscopy+'.flagversions'), 'Backup of flags was not created')
        files = os.listdir(mmscopy+'/SUBMSS')
        print(files)
        for ff in files:
            self.assertFalse(ff.__contains__('flagversions'))


    def test_gaincurve(self):
        applycal(vis=vlbacopy, gaintable=[vlbaGCCal], applymode='calonly')
        datamean = np.mean(getparam(vlbacopy, 'DATA'))
        correctedmean = np.mean(getparam(vlbacopy, 'CORRECTED_DATA'))

        self.assertTrue(correctedmean != datamean)
        self.assertTrue(np.isclose(correctedmean, 2.68767602411+4.91467419904e-06j))

if __name__ == '__main__':
    unittest.main()

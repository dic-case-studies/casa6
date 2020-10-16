##########################################################################
# test_req_task_polcal.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_polcal/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import polcal, casalog
    tb = casatools.table()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy as np

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/vla/polcalTestSet.ms')
    calpath = casatools.ctsys.resolve('caltables/polcalTestSet.ms.Dtrue')
    # circular data
    datapathCirc = casatools.ctsys.resolve('visibilities/vla/polcalTestSetCirc.ms')
    calpathCirc = casatools.ctsys.resolve('caltables/polcalTestSetCirc.ms.Dtrue')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/polcalTestSet.ms'
        calpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/polcalTestSet.ms.Dtrue'
        # circular data
        datapathCirc = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/polcalTestSetCirc.ms'
        calpathCirc = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/polcalTestSetCirc.ms.Dtrue'

    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/polcalTestSet.ms'
        calpath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/polcalTestSet.ms.Dtrue'
        # circular data
        datapathCirc = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/polcalTestSetCirc.ms'
        calpathCirc = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/polcalTestSetCirc.ms.Dtrue'


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
            casalog.post(message=col + ': ERROR in determining the truth value')

    if len(cols) == 0:

        truths = [[x, np.all(truthDict[x] == True)] for x in truthDict.keys()]

    else:

        truths = [[x, np.all(truthDict[x] == True)] for x in cols]

    return np.array(truths)


def getSelection(table, field='0', spw='0'):

    output = []

    tb.open(table)
    fields = tb.getcol('FIELD_ID')
    spws = tb.getcol('SPECTRAL_WINDOW_ID')
    results = tb.getcol('CPARAM')

    # give the cparam result for selected field and spw
    for i in range(tb.nrows()):
        if fields[i] == int(field) and spws[i] == int(spw):
            output.append(results[:,:,i])

    tb.close()

    return output


outcal = 'polcalTestOutput.ms.Df'
datacopy = 'polcalTestCopy.ms'
datacopyCirc = 'polcalTestCopyCirc.ms'


class polcal_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        # Use a copy of the MS
        shutil.copytree(datapath, datacopy)
        shutil.copytree(datapathCirc, datacopyCirc)

    def tearDown(self):
        # After each test remove the copy and the output cal table
        if os.path.exists(datacopy):
            shutil.rmtree(datacopy)
        if os.path.exists(datacopyCirc):
            shutil.rmtree(datacopyCirc)
        if os.path.exists(outcal):
            shutil.rmtree(outcal)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_noPolDflls(self):
        ''' Test on field 0 and spw 0. Use channelized solutions '''

        polcal(vis=datacopy, caltable=outcal, field='0', spw='0', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Dflls',
               smodel=[1.0, 0.0, 0.0, 0.0], refant='5')

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal

        self.assertTrue(np.all(np.isclose(calresult, 0)))

    def test_noPolDlls(self):
        ''' Test on Field 0 and spw 0 there should be no source or instrumental pol. This is probably not a useful case '''

        polcal(vis=datacopy, caltable=outcal, field='0', spw='0', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Dlls',
               smodel=[1.0, 0.0, 0.0, 0.0], refant='5')

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, :10]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    def test_sourcePolDflls(self):
        ''' Test one Field 1 and spw 0. Channelized solution '''

        polcal(vis=datacopy, caltable=outcal, field='1', spw='0', solint='inf', combine='scan',
               preavg=101.0, minsnr=0.0, poltype='Dflls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # Need an unchannelized version of the cal to compare to?

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal

        self.assertTrue(np.all(np.isclose(calresult, 0)))

    def test_sourcePolDlls(self):
        ''' Test on Field 1 and spw 0. Should return the source pol '''

        polcal(vis=datacopy, caltable=outcal, field='1', spw='0', solint='inf', combine='scan',
               preavg=101.0, minsnr=0.0, poltype='Dlls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # Need an unchannelized version of the cal to compare to?

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, :10]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    def test_onlyInstPolDflls(self):
        ''' Test on Field 0 and spw 1. Use channelized solutions '''

        polcal(vis=datacopy, caltable=outcal, field='0', spw='1', solint='inf', combine='scan',
               preavg=101.0, minsnr=0.0, poltype='Dflls', refant='5',
               smodel=[1.0, 0.0, 0.0, 0.0])

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 10:20]
        tb.close()

        resDiff = calresult-refresult

        self.assertTrue(np.all(np.isclose(resDiff, 0)))

    def test_onlyInstPolDlls(self):
        ''' Test on Field 0 and spw 1. Should clear inst pol and see no source pol. needs refant '''

        polcal(vis=datacopy, caltable=outcal, field='0', spw='1', solint='inf', combine='scan',
               preavg=101.0, minsnr=0.0, poltype='Dlls', refant='5',
               smodel=[1.0, 0.0, 0.0, 0.0])

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 10:20]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    def test_sourceInstPolDflls(self):
        ''' Test on Field 1 and spw 1. Use channelized solutions '''

        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf', combine='scan',
               preavg=101.0, minsnr=0.0, poltype='Dflls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 10:20]
        tb.close()

        np.all(np.isclose(calresult, refresult))

    def test_sourceInstPolDlls(self):
        ''' Test on Field 1 and spw 1. Should see source pol '''

        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf', combine='scan',
               preavg=101.0, minsnr=0.0, poltype='Dlls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 10:20]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    def test_sourceInstPolNoise(self):
        ''' Test on Field 1 and spw 2 or 3. Should see source pol with noise '''
        polcal(vis=datacopy, caltable=outcal, field='1', spw='2', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Dflls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 20:30]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult, atol=9e-4)))

    def test_combineSpw(self):
        ''' Test that combine spws increases the SNR/ decreases the noise '''

        polcal(vis=datacopy, caltable=outcal, field='1', spw='2', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Dflls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # TODO come back to this us df?

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 20:30]
        tb.close()

        allChanMeanDiff = (np.mean(calresult) - np.mean(refresult)) / np.sqrt(8)
        singleChanMeanDiff = np.mean(refresult) - np.mean(calresult[:, 0, :])
        meanDiffs = singleChanMeanDiff - allChanMeanDiff

        self.assertTrue(np.isclose(meanDiffs, 0, atol=2e-4))

    def test_noSourceInstPolNoise(self):
        ''' Test on Field 0 and spw 2 or 3. Should just see noise '''
        polcal(vis=datacopy, caltable=outcal, field='0', spw='3', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Dflls',
               smodel=[1.0, 0.08, 0.06, 0.0])

        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpath)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 30:40]
        tb.close()

        meanDifference = np.mean(calresult) - np.mean(refresult)

        self.assertTrue(np.all(np.isclose(meanDifference, 0, atol=4e-3)))

    def test_poltypeXf(self):

        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Xf',
               smodel=[1.0, 0.08, 0.06, 0.0], gaintable=[calpath])

        calresult = getparam(outcal)

        self.assertTrue(np.all(np.isclose(calresult, [1 + 0j])))

        shutil.rmtree(outcal)
        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Xf',
               smodel=[1.0, -0.08, -0.06, 0.0], gaintable=[calpath])

        calresult = getparam(outcal)

        self.assertTrue(np.all(np.isclose(calresult, [-1 + 0j])))


    def test_poltypeXfParangQU(self):

        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='Xfparang+QU',
               smodel=[1.0, 0.08, 0.06, 0.0], gaintable=[calpath])

        # Try flipping the signs here too for smodel
        calresult = getparam(outcal)

        print(calresult)

        self.assertTrue(np.all(np.isclose(calresult, [1+0j])))

    def test_posAngLinear(self):

        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='PosAng',
               smodel=[1.0, 0.08, 0.06, 0.0], gaintable=[calpath])

        tb.open(outcal)
        calresult = tb.getcol('FPARAM')
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, 0)))

        #with more interesting PosAng
        polcal(vis=datacopy, caltable=outcal, field='1', spw='1', solint='inf',
               combine='scan', preavg=101.0, minsnr=0.0, poltype='PosAng',
               smodel=[1.0, 0.06, 0.08, 0.0], gaintable=[calpath])

        tb.open(outcal)
        calresult = tb.getcol('FPARAM')
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, -0.1418972)))

    # Add the Xf then D test case

    ### Circular Test cases ###
    # only use outcal prefix?
    @unittest.skip("Fix in another ticket")
    def test_unpolarizedDfCirc(self):
        ''' Test unpolarized calibration Q=U=0 for Df mode'''
        polcal(vis=datacopyCirc, caltable=outcal, field='0', spw='1,2,3', preavg=101.0,
               minsnr=0.0, poltype='Df', smodel=[1.0, 0, 0, 0], refant='5')
        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpathCirc)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 10:40]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    @unittest.skip("Fix in another ticket")
    def test_unpolarizedDfQU(self):
        ''' Test unpolarized calibraion Q=U=0 for Df+QU mode '''
        polcal(vis=datacopyCirc, caltable=outcal, field='0', spw='1,2,3', preavg=101.0,
               minsnr=0.0, poltype='Df+QU', smodel=[1.0, 0, 0, 0], refant='5')
        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpathCirc)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, 10:40]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    @unittest.skip("Fix in another ticket")
    def test_unknownPolDfQU(self):
        ''' Test polarized calibration with unknown polarization for Df+QU mode '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', preavg=101.0,
               minsnr=0.0, poltype='Df+QU', smodel=[1.0, 0, 0, 0], refant='5')
        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpathCirc)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, :]
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, refresult)))

    def test_knownPolDflls(self):
        ''' Test polarized calibration with known polarization for Dflls mode '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', preavg=101.0,
               minsnr=0.0, poltype='Dflls', smodel=[1.0, 0.08, 0.06, 0], refant='5')
        # get the CPARAM of the polcal run
        calresult = getparam(outcal)
        # get CPARAM values from ref cal
        tb.open(calpathCirc)
        # Each 10 rows is a spectral window
        refresult = tb.getcol('CPARAM')[:, :, :]
        tb.close()

        print(refresult - calresult)

        self.assertTrue(np.all(np.isclose(calresult, refresult, atol=8e-4)))

    ### Xf ###

    def test_XfCorrectQU(self):
        ''' Test poltype Xf and assume the correct Q, U '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='Xf', smodel=[1.0, 0.08, 0.06, 0], refant='5', gaintable=[datapathCirc+'.Dtrue'])

        calresult = getparam(outcal)

        print(calresult)

        self.assertTrue(np.all(np.isclose(calresult, [1 + 0j], atol=8e-4)))

    def test_XfAlternateQU(self):
        ''' Test poltype Xf and assume alternate Q, U values '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='Xf', smodel=[1.0, 0.1, 0.0, 0], refant='5', gaintable=[datapathCirc + '.Dtrue'])

        calresult = getparam(outcal)

        self.assertTrue(np.all(np.isclose(calresult, [0.8 + 0.6j], atol=8e-4)))

    def test_XfNegatedQU(self):
        ''' Test poltype Xf and assume Q, U with flipped signs '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='Xf', smodel=[1.0, -0.08, -0.06, 0], refant='5', gaintable=[datapathCirc + '.Dtrue'])

        calresult = getparam(outcal)

        self.assertTrue(np.all(np.isclose(calresult, [-1 + 0j], atol=8e-4)))

    ### PosAng ###

    def test_PosAngCorrectQU(self):
        ''' Test poltype PosAng and assume the correct Q, U '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='PosAng', smodel=[1.0, 0.08, 0.06, 0], gaintable=[datapathCirc + '.Dtrue'])

        tb.open(outcal)
        calresult = tb.getcol('FPARAM')
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, 0, rtol=1e-4, atol=8e-4)))

    def test_PosAngAlternateQU(self):
        ''' Test poltype PosAng and assume alternate Q, U values '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='PosAng', smodel=[1.0, 0.1, 0.0, 0], gaintable=[datapathCirc + '.Dtrue'])

        tb.open(outcal)
        calresult = tb.getcol('FPARAM')
        tb.close()

        print(np.mean(np.degrees(calresult)))
        # Answer is ~ 0.32193175 radians
        self.assertTrue((np.isclose(np.degrees(np.mean(calresult)), 18.43, atol=1e-3, rtol=1e-3)))

    def test_PosAngNegatedQU(self):
        ''' Test poltype PosAng and assume Q, U with flipped signs '''
        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='PosAng', smodel=[1.0, -0.08, -0.06, 0], gaintable=[datapathCirc + '.Dtrue'])

        tb.open(outcal)
        calresult = tb.getcol('FPARAM')
        tb.close()

        calresult = np.mean(abs(np.degrees(calresult)))
        # ans ~ -1.5706 radians
        self.assertTrue(np.isclose(calresult, 90, atol=1e-3, rtol=1e-4))
    def test_PosAngNegatedQUApply(self):
        ''' Test applying the negated table to a new polcal call '''
        polcal(vis=datacopyCirc, caltable=outcal+'.PA-rel', field='1', spw='', solint='inf',
               minsnr=0.0, poltype='PosAng', smodel=[1.0, -0.08, -0.06, 0], gaintable=[datapathCirc + '.Dtrue'])

        polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
               minsnr=0.0, poltype='PosAng', smodel=[1.0, -0.08, -0.06, 0], gaintable=[datapathCirc+'.Dtrue',
                                                                                       outcal+'.PA-rel'])

        tb.open(outcal)
        calresult = tb.getcol('FPARAM')
        tb.close()

        self.assertTrue(np.all(np.isclose(calresult, 0, atol=1e-4)))

    ### Xfparang+QU

    def test_XfParangQUCorrectQU(self):
        ''' Test poltype XfParang+QU and assume the correct Q, U '''
        P = polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
                   minsnr=0.0, poltype='Xfparang+QU', smodel=[1.0, 0.08, 0.06, 0])

        calresult = getparam(outcal)

        print(calresult)

        self.assertTrue(np.all(np.isclose(calresult, [1 + 0j], atol=8e-4)))

    def test_XfParangQUAlternateQU(self):
        ''' Test poltype XfParang+QU and assume the correct Q, U '''
        P = polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
                   minsnr=0.0, poltype='Xfparang+QU', smodel=[1.0, 0.1, 0.0, 0])

        calresult = getparam(outcal)

        self.assertTrue(np.all(np.isclose(calresult, [0.8 + 0.6j], atol=8e-4)))

    def test_XfParangQUNegatedQU(self):
        ''' Test poltype XfParang+QU  and assume Q, U with flipped signs '''
        P = polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
                   minsnr=0.0, poltype='Xfparang+QU', smodel=[1.0, -0.08, -0.06, 0])

        calresult = getparam(outcal)

        self.assertTrue(np.all(np.isclose(calresult, [-1 + 0j], atol=8e-4)))

    def test_XfParangQUNegatedQUApply(self):
        ''' Test applying the negated table to a new polcal call '''
        polcal(vis=datacopyCirc, caltable=outcal+'.XfpaQU-rel', field='1', spw='', solint='inf',
               minsnr=0.0, poltype='Xfparang+QU', smodel=[1.0, -0.08, -0.06, 0])

        P = polcal(vis=datacopyCirc, caltable=outcal, field='1', spw='', solint='inf',
                   minsnr=0.0, poltype='Xfparang+QU', smodel=[1.0, -0.08, -0.06, 0],
                   gaintable=[outcal + '.XfpaQU-rel', datapathCirc + '.Dtrue'])

        calresult = getparam(outcal)
        print(calresult)

        self.assertTrue(np.all(np.isclose(calresult, [1.0 + 0j], atol=1e-4)))

def suite(self):
    return [polcal_test]


if __name__ == '__main__':
    unittest.main()

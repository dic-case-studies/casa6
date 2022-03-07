##########################################################################
# test_tool_agentflagger.py
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.exportuvfits.html
#
##########################################################################

import os
import shutil
import unittest
import math
import numpy as np
import numbers

from casatools import ctsys, table, msmetadata
from casatasks import exportuvfits, importuvfits, applycal, split, casalog
_tb = table()
ctsys_resolve = ctsys.resolve
_msmd = msmetadata()

datapath = ctsys_resolve('unittest/exportuvfits/')
testdata = 'gaincalcopy.ms'
output = 'uvfitstest.uvfits'
reimport = 'uvfitsreimport.ms'
gaincaltable = os.path.join(datapath, 'gaincaltest2.ms.G0')
splitdata = 'splitdata.ms'
testlog = 'testlog.log'

class exportuvfits_test(unittest.TestCase):

    def setUp(self):
        shutil.copytree(os.path.join(datapath, 'gaincaltest2.ms'), testdata)
    
    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)

        _tb.close()
        if os.path.exists(output):
            os.remove(output)
        if os.path.exists(testdata):
            shutil.rmtree(testdata)
        if os.path.exists(reimport):
            shutil.rmtree(reimport)
        if os.path.exists(splitdata):
            shutil.rmtree(splitdata)
        if os.path.exists(testdata+'.flagversions'):
            shutil.rmtree(testdata+'.flagversions')

        if os.path.exists('imported_no_restfreqs.ms'):
            shutil.rmtree('imported_no_restfreqs.ms')
        if os.path.exists('no_rest_freqs.uvfits'):
            os.remove('no_rest_freqs.uvfits')
        if os.path.exists('no_source_table.ms'):
            shutil.rmtree('no_source_table.ms')
        if os.path.exists('no_source_table.uvfits'):
            os.remove('no_source_table.uvfits')
        if os.path.exists('uvfits_test.ms'):
            shutil.rmtree('uvfits_test.ms')
        if os.path.exists('rest_freq_test.ms'):
            shutil.rmtree('rest_freq_test.ms')
        if os.path.exists('CAS-5492.uvfits'):
            os.remove('CAS-5492.uvfits')
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')

    def test_exportOverwrite(self):
        """CAS-5492: test the overwrite parameter when exporting MSes to uvfits"""
        msname = "uvfits_test.ms" 
        shutil.copytree(os.path.join(datapath, msname), msname)
        fitsname = "CAS-5492.uvfits"
        res = exportuvfits(vis=msname, fitsfile=fitsname)
        # Not sure why all of a sudden CASA6 is returning None for tasks
        self.assertTrue(res == None, "Failed exportuvfits")

        self.assertRaises(
            Exception, exportuvfits, vis=msname, fitsfile=fitsname, overwrite=False,
            msg="exportuvfits succeeded but should have failed because "
            + "overwrite=False"
        )

        # succeed because overwrite=True
        res = exportuvfits(vis=msname, fitsfile=fitsname, overwrite=True)

        self.assertTrue(
            res == None,
            "exportuvfits failed but should have succeeded because "
            + "overwrite=True"
        )
        

    def test_noRestFreqs(self):
        """CAS-11514: test exporting an MS with no rest frequencies in the SOURCE table"""
        msname = "rest_freq_test.ms"
        shutil.copytree(os.path.join(datapath, msname), msname)
        fitsname = "no_rest_freqs.uvfits"
        res = exportuvfits(vis=msname, fitsfile=fitsname)
        self.assertTrue(res == None, "Failed exportuvfits with no rest freqs")
        # import and check the rest freqs
        # importuvfits doesn't return anything, so we cannot test the
        # return value for success
        msname = "imported_no_restfreqs.ms"
        importuvfits(fitsfile=fitsname, vis=msname)
        _msmd.open(msname)
        restfreqs = _msmd.restfreqs()
        _msmd.done()
        expec = {
            '0': {
                'type': 'frequency', 'm0': {'value': 0.0, 'unit': 'Hz'},
                'refer': 'LSRK'
            }
        }
        self.assertEqual(
            restfreqs, expec, "Got wrong restfreqs from re-imported dataset"
        )

    def test_no_source_table(self):
        """CAS-11514: test exporting an MS with no rest frequencies in the SOURCE table"""
        msname = "no_source_table.ms"
        shutil.copytree(os.path.join(datapath, msname), msname)
        fitsname = "no_source_table.uvfits"
        res = exportuvfits(vis=msname, fitsfile=fitsname)
        self.assertTrue(res == None, "Failed exportuvfits with no SOURCE table")

    def test_basicExport(self):
        '''Check that a ms can be exported as a uvfits with no additional parameters'''
        exportuvfits(vis=testdata, fitsfile=output)

        self.assertTrue(os.path.exists(output))

    def test_fieldSelection(self):
        '''Check that field selection properly selects a subset of the data'''
        # Export with field selection
        # Reimport and check that the fields were selected
        exportuvfits(vis=testdata, fitsfile=output, field='0')
        importuvfits(fitsfile=output, vis=reimport)

        expected_fields = {0}

        _tb.open(reimport)
        fields = {i for i in _tb.getcol('FIELD_ID')}
        _tb.close()

        self.assertTrue(fields == expected_fields)

    def test_spwSelection(self):
        '''Check that spw selection properly selects a subset of the data'''
        exportuvfits(vis=testdata, fitsfile=output, spw='0')
        importuvfits(fitsfile=output, vis=reimport)

        expected_spws = 1

        _tb.open(reimport+'/SPECTRAL_WINDOW')
        spws = len(_tb.getcol('NUM_CHAN'))
        _tb.close()

        self.assertTrue(spws == expected_spws)

    def test_antennaSelection(self):
        '''Check that the antenna parameter selects a subset of the data'''
        exportuvfits(vis=testdata, fitsfile=output, antenna='0')
        importuvfits(fitsfile=output, vis=reimport)

        expected_antennas = {0}

        _tb.open(reimport)
        antennas = {i for i in _tb.getcol('ANTENNA1')}
        _tb.close()

        print(antennas, expected_antennas)
        self.assertTrue(antennas == expected_antennas)

    def test_timerangeSelection(self):
        '''Check that the timerange parameter selects a subset of the data'''
        exportuvfits(vis=testdata, fitsfile=output, timerange='>5000')
        importuvfits(fitsfile=output, vis=reimport)

        _tb.open(reimport)
        selected = len(_tb.getcol('TIME'))
        _tb.close()

        self.assertTrue(selected == 95940)

    def test_columnSelection(self):
        '''Check that the data input column can be selected with the datacolumn parameter'''
        # Corrected data
        applycal(vis=testdata, gaintable=[gaincaltable])

        exportuvfits(vis=testdata, fitsfile=output, datacolumn='data')
        importuvfits(fitsfile=output, vis=reimport)
        _tb.open(reimport)
        res1 = np.mean(_tb.getcol('DATA'))
        _tb.close()

        os.remove(output)
        shutil.rmtree(reimport)

        exportuvfits(vis=testdata, fitsfile=output, datacolumn='corrected')
        importuvfits(fitsfile=output, vis=reimport)
        _tb.open(reimport)
        res2 = np.mean(_tb.getcol('DATA'))
        _tb.close()

        self.assertFalse(res1 == res2)

    def test_multiSourceImage(self):
        '''Check multisource is not overwritten and one source is present'''
        log = casalog.logfile()
        casalog.setlogfile(testlog)
        exportuvfits(vis=testdata, fitsfile=output, field='1', multisource=False)
        casalog.setlogfile(log)
        importuvfits(fitsfile=output, vis=reimport)

        # Check that it is not overwritten
        overwritten = False
        with open(testlog) as fout:
            for line in fout:
                if 'Multiple sources are present, thus written as a multi-source FITS file' in line:
                    overwritten = True
        self.assertFalse(overwritten)

        # Check that there is only one source
        _tb.open(reimport+'/SOURCE')
        ids = _tb.getcol('SOURCE_ID')
        _tb.close()

        self.assertTrue(np.all(ids == 0))


    def test_multiSourceOverwrite(self):
        '''Check that if mutisource is False when multiple sources are present then it will use multisource'''
        # switch to temp log
        log = casalog.logfile()
        casalog.setlogfile(testlog)
        exportuvfits(vis=testdata, fitsfile=output, multisource=False)
        # switch back to main log and look for the overwrite in the log
        casalog.setlogfile(log)

        overwritten = False
        with open(testlog) as fout:
            for line in fout:
                if 'Multiple sources are present, thus written as a multi-source FITS file' in line:
                    overwritten = True

        self.assertTrue(overwritten)

    def test_combineSpws(self):
        '''Check that combine spw combines all spws into one frequency group'''
        exportuvfits(vis=testdata, fitsfile=output, combinespw=False)
        importuvfits(fitsfile=output, vis=reimport)

        _tb.open(reimport+'/SPECTRAL_WINDOW')
        freqGroup = _tb.getcol('FREQ_GROUP')
        _tb.close()

        self.assertTrue(np.all(freqGroup == [0, 1, 2, 3]))

    def test_useStationName(self):
        '''Check that the station names are written with writestation=True'''
        # With station written
        exportuvfits(vis=testdata, fitsfile=output, writestation=True)
        importuvfits(fitsfile=output, vis=reimport)

        _tb.open(reimport+'/ANTENNA')
        withStation = _tb.getcol('STATION')
        _tb.close()

        # Without write station
        os.remove(output)
        shutil.rmtree(reimport)
        exportuvfits(vis=testdata, fitsfile=output, writestation=False)
        importuvfits(fitsfile=output, vis=reimport)

        _tb.open(reimport + '/ANTENNA')
        withoutStation = _tb.getcol('STATION')
        _tb.close()

        self.assertFalse(np.all(withoutStation == withStation))

    def test_padWithFlags(self):
        '''Check that missisng data is filled with flags'''
        # first remove spw 0 data
        _tb.open(testdata, nomodify=False)
        flag = _tb.getcol('FLAG')

        for i in range(len(flag[0][0])):
            for j in range(8):
                flag[0, j, i] = True
        _tb.putcol('FLAG', flag)
        _tb.close()

        split(vis=testdata, outputvis=splitdata, datacolumn='data', keepflags=False)

        # export and reimport split data
        exportuvfits(vis=splitdata, fitsfile=output, padwithflags=True)
        importuvfits(fitsfile=output, vis=reimport)

        # Check that the missing data has been added with flags
        _tb.open(reimport)
        res = _tb.getcol('FLAG')
        _tb.close()

        self.assertFalse(np.all(res == False))

    def test_missingWeights(self):
        '''Check that a WEIGHT_SPECTRUM column is created and filled if one does not exist'''
        exportuvfits(vis=testdata, fitsfile=output)
        importuvfits(fitsfile=output, vis=reimport)

        _tb.open(testdata)
        has_weight = 'WEIGHT_SPECTRUM' in _tb.colnames()
        _tb.close()
        self.assertFalse(has_weight)

        _tb.open(reimport)
        has_weight = 'WEIGHT_SPECTRUM' in _tb.colnames()
        _tb.close()
        self.assertTrue(has_weight)


if __name__ == '__main__':
    unittest.main()
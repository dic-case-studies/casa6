from __future__ import absolute_import
from __future__ import print_function
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


datapath = 'unittest/exportuvfits/'


class exportuvfits_test(unittest.TestCase):

    def setUp(self):
        self.datapath = ctsys_resolve(datapath)
        self.testdata = 'gaincalcopy.ms'
        self.output = 'uvfitstest.uvfits'
        self.reimport = 'uvfitsreimport.ms'
        self.gaincaltable = os.path.join(self.datapath, 'gaincaltest2.ms.G0')
        self.splitdata = 'splitdata.ms'
        self.testlog = 'testlog.log'

        shutil.copytree(os.path.join(self.datapath, 'gaincaltest2.ms'), self.testdata)
    
    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)
        # make sure directory is clean as per verification test requirement
        #cwd = os.getcwd()
        #for filename in os.listdir(cwd):
        #    file_path = os.path.join(cwd, filename)
        #    try:
        #        if os.path.isfile(file_path) or os.path.islink(file_path):
        #            os.unlink(file_path)
        #        elif os.path.isdir(file_path):
        #            # CASA 5 tests need this directory
        #            if filename != 'xml':
        #                shutil.rmtree(file_path)
        #    except Exception as e:
        #        print('Failed to delete %s. Reason: %s' % (file_path, e))

        _tb.close()
        if os.path.exists(self.output):
            os.remove(self.output)
        if os.path.exists(self.testdata):
            shutil.rmtree(self.testdata)
        if os.path.exists(self.reimport):
            shutil.rmtree(self.reimport)
        if os.path.exists(self.splitdata):
            shutil.rmtree(self.splitdata)
        if os.path.exists(self.testdata+'.flagversions'):
            shutil.rmtree(self.testdata+'.flagversions')

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

    def test_export_overwrite(self):
        """CAS-5492: test the overwrite parameter when exporting MSes to uvfits"""
        msname = "uvfits_test.ms" 
        shutil.copytree(os.path.join(self.datapath, msname), msname)
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
        

    def test_no_rest_freqs(self):
        """CAS-11514: test exporting an MS with no rest frequencies in the SOURCE table"""
        msname = "rest_freq_test.ms"
        shutil.copytree(os.path.join(self.datapath, msname), msname)
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
        shutil.copytree(os.path.join(self.datapath, msname), msname)
        fitsname = "no_source_table.uvfits"
        res = exportuvfits(vis=msname, fitsfile=fitsname)
        self.assertTrue(res == None, "Failed exportuvfits with no SOURCE table")

    def test_basicExport(self):
        '''Check that a ms can be exported as a uvfits with no additional parameters'''
        exportuvfits(vis=self.testdata, fitsfile=self.output)

        self.assertTrue(os.path.exists(self.output))

    def test_fieldSelection(self):
        '''Check that field selection properly selects a subset of the data'''
        # Export with field selection
        # Reimport and check that the fields were selected
        exportuvfits(vis=self.testdata, fitsfile=self.output, field='0')
        importuvfits(fitsfile=self.output, vis=self.reimport)

        fields = set()
        expected_fields = {0}

        _tb.open(self.reimport)
        for i in _tb.getcol('FIELD_ID'):
            fields.add(i)
        _tb.close()

        self.assertTrue(fields == expected_fields)

    def test_spwSelection(self):
        '''Check that spw selection properly selects a subset of the data'''
        exportuvfits(vis=self.testdata, fitsfile=self.output, spw='0')
        importuvfits(fitsfile=self.output, vis=self.reimport)

        expected_spws = 1

        _tb.open(self.reimport+'/SPECTRAL_WINDOW')
        spws = len(_tb.getcol('NUM_CHAN'))
        _tb.close()

        self.assertTrue(spws == expected_spws)

    def test_antennaSelection(self):
        '''Check that the antenna parameter selects a subset of the data'''
        exportuvfits(vis=self.testdata, fitsfile=self.output, antenna='0')
        importuvfits(fitsfile=self.output, vis=self.reimport)

        antennas = set()
        expected_antennas = {0}

        _tb.open(self.reimport)
        for i in _tb.getcol('ANTENNA1'):
            antennas.add(i)
        _tb.close()

        print(antennas, expected_antennas)
        self.assertTrue(antennas == expected_antennas)

    def test_timerangeSelection(self):
        '''Check that the timerange parameter selects a subset of the data'''
        exportuvfits(vis=self.testdata, fitsfile=self.output, timerange='>5000')
        importuvfits(fitsfile=self.output, vis=self.reimport)

        _tb.open(self.reimport)
        selected = len(_tb.getcol('TIME'))
        _tb.close()

        self.assertTrue(selected == 95940)

    def test_columnSelection(self):
        # Corrected data
        applycal(vis=self.testdata, gaintable=[self.gaincaltable])

        exportuvfits(vis=self.testdata, fitsfile=self.output, datacolumn='data')
        importuvfits(fitsfile=self.output, vis=self.reimport)
        _tb.open(self.reimport)
        res1 = np.mean(_tb.getcol('DATA'))
        _tb.close()

        os.remove(self.output)
        shutil.rmtree(self.reimport)

        exportuvfits(vis=self.testdata, fitsfile=self.output, datacolumn='corrected')
        importuvfits(fitsfile=self.output, vis=self.reimport)
        _tb.open(self.reimport)
        res2 = np.mean(_tb.getcol('DATA'))
        _tb.close()

        self.assertFalse(res1 == res2)

    def test_multiSourceImage(self):
        '''Check multisource is not overwritten and one source is present'''
        log = casalog.logfile()
        casalog.setlogfile(self.testlog)
        exportuvfits(vis=self.testdata, fitsfile=self.output, field='1', multisource=False)
        casalog.setlogfile(log)
        importuvfits(fitsfile=self.output, vis=self.reimport)

        # Check that it is not overwritten
        overwritten = False
        with open(self.testlog) as fout:
            for line in fout:
                if 'Multiple sources are present, thus written as a multi-source FITS file' in line:
                    overwritten = True
        self.assertFalse(overwritten)

        # Check that there is only one source
        _tb.open(self.reimport+'/SOURCE')
        ids = _tb.getcol('SOURCE_ID')
        _tb.close()

        self.assertTrue(np.all(ids == 0))


    def test_multiSourceOverwrite(self):
        '''Check that if mutisource is False when multiple sources are present then it will use multisource'''
        # switch to temp log
        log = casalog.logfile()
        casalog.setlogfile(self.testlog)
        exportuvfits(vis=self.testdata, fitsfile=self.output, multisource=False)
        # switch back to main log and look for the overwrite in the log
        casalog.setlogfile(log)

        overwritten = False
        with open(self.testlog) as fout:
            for line in fout:
                if 'Multiple sources are present, thus written as a multi-source FITS file' in line:
                    overwritten = True

        self.assertTrue(overwritten)

    def test_combineSpws(self):
        '''Check that combine spw combines all spws into one frequency group'''
        exportuvfits(vis=self.testdata, fitsfile=self.output, combinespw=False)
        importuvfits(fitsfile=self.output, vis=self.reimport)

        _tb.open(self.reimport+'/SPECTRAL_WINDOW')
        freqGroup = _tb.getcol('FREQ_GROUP')
        _tb.close()

        self.assertTrue(np.all(freqGroup == [0, 1, 2, 3]))

    def test_useStationName(self):
        '''Check that the station names are written with writestation=True'''
        # With station written
        exportuvfits(vis=self.testdata, fitsfile=self.output, writestation=True)
        importuvfits(fitsfile=self.output, vis=self.reimport)

        _tb.open(self.reimport+'/ANTENNA')
        withStation = _tb.getcol('STATION')
        _tb.close()

        # Without write station
        os.remove(self.output)
        shutil.rmtree(self.reimport)
        exportuvfits(vis=self.testdata, fitsfile=self.output, writestation=False)
        importuvfits(fitsfile=self.output, vis=self.reimport)

        _tb.open(self.reimport + '/ANTENNA')
        withoutStation = _tb.getcol('STATION')
        _tb.close()

        self.assertFalse(np.all(withoutStation == withStation))

    def test_padWithFlags(self):
        # first remove spw 0 data
        _tb.open(self.testdata, nomodify=False)
        flag = _tb.getcol('FLAG')

        for i in range(len(flag[0][0])):
            for j in range(8):
                flag[0, j, i] = True
        _tb.putcol('FLAG', flag)
        _tb.close()

        split(vis=self.testdata, outputvis=self.splitdata, datacolumn='data', keepflags=False)

        # export and reimport split data
        exportuvfits(vis=self.splitdata, fitsfile=self.output, padwithflags=True)
        importuvfits(fitsfile=self.output, vis=self.reimport)

        # Check that the missing data has been added with flags
        _tb.open(self.reimport)
        res = _tb.getcol('FLAG')
        _tb.close()

        self.assertFalse(np.all(res == False))

    def test_missingWeights(self):
        exportuvfits(vis=self.testdata, fitsfile=self.output)
        importuvfits(fitsfile=self.output, vis=self.reimport)

        _tb.open(self.testdata)
        has_weight = 'WEIGHT_SPECTRUM' in _tb.colnames()
        _tb.close()
        self.assertFalse(has_weight)

        _tb.open(self.reimport)
        has_weight = 'WEIGHT_SPECTRUM' in _tb.colnames()
        _tb.close()
        self.assertTrue(has_weight)

            
def suite():
    return [exportuvfits_test]        


if __name__ == '__main__':
    unittest.main()
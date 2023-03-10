##########################################################################
# test_task_listobs.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.information.listobs.html
#
# Each requirement is broken up into 4 versions testing it on a MS, MMS, time-averaged MS, and time-averaged MMS
#
# test_logread checks that the information is output in a form that follows logger guidelines
# test_file checks that the listfile parameter operates as intended by creating a file
# test_input checks that listobs can take a MS, MMS, time-averaged MS, or time-averaged MMS in the input file
# test_logfile checks that the logfile is populated when new listobs commands are run
# test_indivrow checks that one row is displayed per scan
# test_novis checks to be sure that listobs gives no visibility information
# test_scan checks that the scan information exists and the parameter accepts proper inputs
# test_field checks that the field information exists and the parameter accepts proper inputs
# test_spw checks that the spectral window information exists and the parameter accepts proper inputs
# test_corrs checks that the correlation information exists and the parameter accepts proper inputs
# test_ant checks that the antenna information exists and the parameter accepts proper inputs
# test_selectdata checks that the selectdata parameter works as intended and properly handles inputs
# test_uvrange checks that the parameter properly handles inputs
# test_timerange checks that the parameter properly handles inputs
# test_intent checks that intent information exists and inputs are properly handled
# test_array checks for the existence of the array parameter and that inputs are properly handled
# test_observation checks for the existence of the parameter and that inputs are properly handled
# test_verbose checks that verbose outputs are longer than non-verbose
# test_overwrite checks that the overwritten file is the same as the original
# test_CAS_6733 checks for an infinite loop bug
# test_avgInterval checks for the existence of the average int information
# test_listunfl checks that unflagged information is displayed by listobs
#
###########################################################################
import string
import sys
import os
import unittest
import hashlib
import subprocess
import shutil
import copy
import numpy as np
from casatestutils.compare import compare_dictionaries

import casatools
from casatasks import partition, split, listobs, casalog

ms = casatools.ms()

from casatestutils import listing as lt

    # Generate the test data
datapath = casatools.ctsys.resolve('unittest/listobs/')

# This is for tests that check what the parameter validator does when parameters are
# given wrong types - these don't exercise the task but the parameter validator!
validator_exc_type = AssertionError

# Input data
mesSet = os.path.join(datapath,'uid___X02_X3d737_X1_01_small.ms')
# Data for old test
msfile1Orig = os.path.join(datapath,'ngc5921_ut.ms')
msfile2Orig = os.path.join(datapath,'uid___X02_X3d737_X1_01_small.ms')
nep = os.path.join(datapath,'nep2-shrunk.ms')
msfile1 = 'ngc5921_ut.ms'
msfile2 = 'uid___X02_X3d737_X1_01_small.ms'
msfile3 = os.path.join(datapath, 'CAS-6733.ms')
msfile4 = os.path.join(datapath, 'lofar_small_dysco.ms')


outvis = 'genmms.mms'
if not os.path.exists(outvis):
    partition(vis=mesSet, outputvis=outvis, createmms=True)
multiMesSet = outvis

outvis = 'gentimeavgms.ms'
if not os.path.exists(outvis):
    split(vis=mesSet, outputvis=outvis, datacolumn='DATA', timebin='1s')
timeavg_ms = outvis

outvis = 'gentimeavgmms.mms'
if not os.path.exists(outvis):
    split(vis=multiMesSet, outputvis=outvis, datacolumn='DATA', timebin='1s')
timeavg_mms = outvis

logpath = casalog.logfile()
# nep = 'nep2-shrunk.ms'
# Old reffiles
reffile = os.path.join(datapath, 'listobs_reference/reflistobs')


def _sha1it(filename):
    blocksize = 65536
    hasher = hashlib.sha1()
    with open(filename, 'rb') as afile:
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)
    return hasher.hexdigest()


class listobs_test_base(unittest.TestCase):

    ref_return = {}
    ref_return_base =\
            {'BeginTime': 55248.126073333326,
             'EndTime': 55248.130800000006,
             'IntegrationTime': 408.38400173187256,
             'field_0': {'code': 'none',
                         'direction': {'m0': {'unit': 'rad', 'value': 1.4433872913993107},
                                       'm1': {'unit': 'rad', 'value': 0.2361430477948328},
                                       'refer': 'J2000',
                                       'type': 'direction'},
                         'name': 'J0530+135'},
             'field_1': {'code': 'none',
                         'direction': {'m0': {'unit': 'rad', 'value': 0.0},
                                       'm1': {'unit': 'rad', 'value': 0.0},
                                       'refer': 'J2000',
                                       'type': 'direction'},
                         'name': 'Mars'},
             'field_2': {'code': 'none',
                         'direction': {'m0': {'unit': 'rad', 'value': 1.6056810952191334},
                                       'm1': {'unit': 'rad', 'value': -0.14975884913199947},
                                       'refer': 'J2000',
                                       'type': 'direction'},
                         'name': 'J0607-085'},
             'nfields': 3,
             'numrecords': 1080,
             'scan_1': {'0': {'BeginTime': 55248.126073333326,
                              'EndTime': 55248.12846666667,
                              'FieldId': 0,
                              'FieldName': 'J0530+135',
                              'IntegrationTime': 3.0240000000000187,
                              'SpwIds': np.array([0, 1], dtype=np.int32),
                              'StateId': 0,
                              'nRow': 600,
                              'scanId': 1}},
             'scan_2': {'0': {'BeginTime': 55248.12877055556,
                              'EndTime': 55248.13014111111,
                              'FieldId': 1,
                              'FieldName': 'Mars',
                              'IntegrationTime': 3.023999999999993,
                              'SpwIds': np.array([0, 1], dtype=np.int32),
                              'StateId': 5,
                              'nRow': 360,
                              'scanId': 2}},
             'scan_3': {'0': {'BeginTime': 55248.130450000004,
                              'EndTime': 55248.130800000006,
                              'FieldId': 2,
                              'FieldName': 'J0607-085',
                              'IntegrationTime': 3.024000000000007,
                              'SpwIds': np.array([0, 1], dtype=np.int32),
                              'StateId': 8,
                              'nRow': 120,
                              'scanId': 3}},
             'timeref': 'UTC'}
    ref_return['uid___X02_X3d737_X1_01_small.ms'] = ref_return_base
    # Numerical differences in MMS
    ref_return['genmms.mms'] = copy.deepcopy(ref_return_base)
    ref_return['genmms.mms']['scan_1']['0']['IntegrationTime'] = 3.0240000000000475
    ref_return['genmms.mms']['scan_2']['0']['IntegrationTime'] = 3.023999999999984
    ref_return['gentimeavgms.ms'] = ref_return_base
    ref_return['gentimeavgmms.mms'] = ref_return['genmms.mms']

    def setUp(self):
        pass

    def tearDown(self):
        # remove files and temp logs
        os.system('rm -rf ' + 'listobs*.txt')
        os.system('rm -rf testlog.log')
        casalog.setlogfile(str(logpath))

    def check_return_dict(self, result, dataset):
        """
        Ensure the dictionary returned from listobs is as expected.
        @param result: result dictionary from a listob execution
        @param dataset: dataset/vis name
        """

        filename = os.path.basename(dataset)
        # See CAS-13170. Ref values come from RHEL. There are differences ~10-16 on Mac
        self.assertTrue(compare_dictionaries(result, self.ref_return[filename], rtol=1e-13,
                                             atol=1e-13))

    def check_file_plus_dict(self, dataset, filename):
        res = listobs(vis=dataset, listfile=filename)
        self.assertTrue(os.path.isfile('listobs.txt'))
        self.check_return_dict(res, dataset)

    def check_logread_plus_dict(self, dataset):
        casalog.setlogfile('testlog.log')
        res = listobs(vis=dataset)
        casalog.setlogfile(logpath)

        print('VERSION', ' ', sys.version_info)
        # Check that the file can be read in python session default encoding
        with open('testlog.log', 'r', encoding=sys.getdefaultencoding( )) as fout:
            fout.readlines( )

        self.check_return_dict(res, dataset)

    def check_logfilecontain_plus_dict(self, dataset):
        casalog.setlogfile('testlog.log')
        res = listobs(vis=dataset)
        self.assertTrue('listobs' in open('testlog.log').read(), msg='logfile not populated by listobs command on a MS')
        self.check_return_dict(res, dataset)

    def novis(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertFalse('Amp' in open('listobs.txt').read())

    def scancheck(self, dataset):
        # Check that listobs runs and that the scan column exists
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Scan' in open('listobs.txt').read(), msg='scan does not exist in output')
        # Check to see if a valid input and invalid input is returns the proper value (pass/fail)
        try:
            listobs(vis=dataset, scan='1')
        except Exception:
            self.fail('Scan fails to select')
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, scan=1)

        # Make temp log
        casalog.setlogfile('testlog.log')

        # test for certain warnings appearing in the log
        listobs(vis=dataset, scan='1,2')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='A warning is raised for multiple scans')
        listobs(vis=dataset, scan=['1', '2'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='fails to label incorrect data type')
        listobs(vis=dataset, scan='abc')
        self.assertTrue('Parse error' in open('testlog.log').read(), msg='fails to recognize improper string')

    def fieldcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, field=1)

        self.assertTrue('FldId' in open('listobs.txt').read(), msg='Field Id does not exist in a MS')
        # section that should raise no warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, field='1')
        listobs(vis=dataset, field='0~2')
        listobs(vis=dataset, field='0,2')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='not accepting proper input for field in a MS')
        # section that should be raising warnings
        listobs(vis=dataset, field='0-2')
        self.assertTrue('No match found for name "0-2"' in open('testlog.log').read(),
                        msg='Failed to identify improper delimiter')
        listobs(vis=dataset, field='abc')
        self.assertTrue('No match found for name "abc"' in open('testlog.log').read(),
                        msg='Failed to identify improper string')
        # return to default log file

    def spwcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, spw=1)

        self.assertTrue('SpwID' in open('listobs.txt').read(), msg='Spw does not exist in a MS')
        # section that should raise no warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, spw='0')
        listobs(vis=dataset, spw='0,1')
        listobs(vis=dataset, spw='0~1')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='not accepting proper input for spw in a MS')
        # section that should raise warnings
        listobs(vis=dataset, spw='3')
        self.assertTrue('No match found for 3' in open('testlog.log').read(),
                        msg='fails to recognize out of range values')
        listobs(vis=dataset, spw='0')
        self.assertTrue('-1' in open('testlog.log').read(), msg='Fails to recognize improper delimiter')
        listobs(vis=dataset, spw='abc')
        self.assertTrue('No match found for "abc"' in open('testlog.log').read(),
                        msg='Fails to recognize improper string')
        # return to default log file

    def corrcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Corrs' in open('listobs.txt').read(), msg='Corrs does not exist in a MS')
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, correlation=1)

        # section that should not raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, correlation='XX')
        listobs(vis=dataset, correlation='XX,YY')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='not accepting proper input for correlation in a MS')
        # section that should raise warnings
        listobs(vis=dataset, correlation=['XX', 'YY'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(),
                        msg='No warning for using a list was given')
        listobs(vis=dataset, correlation='RR')
        self.assertTrue('named RR' in open('testlog.log').read(), msg='No warning for using a absent correlation')

    def antcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Antennas' in open('listobs.txt').read(), msg='Antennas section does not exist in MS')
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, antenna=0)

        # section that should not raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, antenna='0')
        listobs(vis=dataset, antenna='0,DV01')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='Proper inputs raised warnings')
        # section that should raise warnings
        listobs(vis=dataset, antenna='abc')
        self.assertTrue('Antenna Expression: No match found for token(s)' in open('testlog.log').read(),
                        msg='No warning raise for invalid string')
        listobs(vis=dataset, antenna='3')
        self.assertTrue('No match found for the antenna specificion [ID(s): [3]]' in open('testlog.log').read(),
                        msg='No warning for ID out of range')
        # This one is marked as correct by the documentation, but CASA disagrees
        listobs(vis=dataset, antenna=['0,DV01'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(),
                        msg='Failed to recognize list as incorrect data type')
        # return to default log file

    def selectcheck(self, dataset):
        # selectdata should fail for all invalid inputs
        try:
            listobs(vis=dataset, selectdata=False)
        except Exception:
            self.fail('Passing False to select data fails on a MS')
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, selectdata=1)
        with self.assertRaises(validator_exc_type):
            listobs(vis=dataset, selectdata='str')

    def uvrangecheck(self, dataset):
        try:
            listobs(vis=dataset, uvrange='0~100klambda')
        except Exception:
            self.fail('fails to read valid input for uvrange in a MS')
        with self.assertRaises(AssertionError):
            listobs(vis=dataset, uvrange=0)
        with self.assertRaises(AssertionError):
            listobs(vis=dataset, uvrange=[1, 2])

        # Use temp log
        casalog.setlogfile('testlog.log')
        # shouldn't raise Warning
        listobs(vis=dataset, uvrange='0~100')
        listobs(vis=dataset, uvrange='0~100klambda')
        listobs(vis=dataset, uvrange='0~50,60~100')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='Warnings are raised for valid inputs')
        # should raise warnings
        listobs(vis=dataset, uvrange=['0~50', '60~100'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(),
                        msg='Fails to raise warning for wrong data type')
        listobs(vis=dataset, uvrange='0-100')
        self.assertTrue('near char. 2 in string "0-100"' in open('testlog.log').read(),
                        msg='Fails to raise warning for wrong delimiter')
        listobs(vis=dataset, uvrange='abc')
        self.assertTrue('near char. 1 in string "abc"' in open('testlog.log').read(),
                        msg='Fails to raise warning for improper string')
        # restore default log path

    def timerangecheck(self, dataset):
        # check valid entry
        try:
            listobs(vis=dataset, timerange='03:00:00~04:00:00')
        except Exception:
            self.fail()
        # create temp log and check inputs that raise no warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, timerange='3:0:0~4:0:0,4:0:0~5:0:0')
        self.assertFalse('WARN' in open('testlog.log').read())
        # check that specific warnings are raised
        listobs(vis=dataset, timerange='abc')
        self.assertTrue('Parse error at or near ' in open('testlog.log').read())
        listobs(vis=dataset, timerange=[])
        self.assertTrue('incorrect data type' in open('testlog.log').read())
        listobs(vis=dataset, timerange='03:00:00-04:00:00')
        self.assertTrue('near char. 9 in string "03:00:00-04:00:00"' in open('testlog.log').read())
        listobs(vis=dataset, timerange='3~4')
        self.assertTrue('MSSelectionNullSelection' in open('testlog.log').read())
        # Check that passing an int fails
        with self.assertRaises(AssertionError):
            listobs(vis=dataset, timerange=4)

    def intentcheck(self, dataset):
        # Returns true with a valid input and false with an int
        try:
            listobs(vis=dataset,
                    intent='CALIBRATE_PHASE.ON_SOURCE,CALIBRATE_POINTING.ON_SOURCE,CALIBRATE_WVR.ON_SOURCE')
        except Exception:
            self.fail('Fails with valid input')
        with self.assertRaises(AssertionError):
            listobs(vis=dataset, intent=1)

        # Test for the existence of the column scan intent
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('ScanIntent' in open('listobs.txt').read(), msg='There is no ScanIntent information for a MS')
        # These shouldn't raise any warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, intent='CALIBRATE_PHASE.ON_SOURCE')
        listobs(vis=dataset, intent='CALIBRATE_PHASE.ON_SOURCE,CALIBRATE_POINTING.ON_SOURCE,CALIBRATE_WVR.ON_SOURCE')
        self.assertFalse('WARN' in open('testlog.log').read(),
                         msg='There are warnings for inputs that should raise none')
        # These should raise a warning
        listobs(vis=dataset, intent=[])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='Incorrect data type list accepted')
        listobs(vis=dataset, intent='abc')
        self.assertTrue('No match found for "abc"' in open('testlog.log').read(),
                        msg='Invalid string accepted without warning')
        # Set log path back to default

    def arraycheck(self, dataset):
        try:
            listobs(vis=dataset, array='0')
        except Exception:
            self.fail('Listobs fails to recognize valid array in a MS')
        with self.assertRaises(AssertionError):
            listobs(vis=dataset, array=0)

        self.res = listobs(vis=dataset, listfile='listobs.txt')
        # These should raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, array='abc')
        self.assertTrue('Parse error' in open('testlog.log').read(),
                        msg='Listobs fails to recognize invalid string in a MS')
        listobs(vis=dataset, array='10')
        self.assertTrue('The selected table has zero rows' in open('testlog.log').read(),
                        msg='Listobs fails to recognize empty table from a MS')

        self.assertTrue('ArrayID' in open('listobs.txt').read(), msg='There is no Array information for a MS')

    def obscheck(self, dataset):
        try:
            listobs(vis=dataset, observation=0)
        except Exception:
            self.fail('Observation fails to accept Int for a MS')
        try:
            listobs(vis=dataset, observation='0')
        except Exception:
            self.fail('Observation fails to accept proper string for a MS')
        # These should raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, observation='abc')
        self.assertTrue('Parse error' in open('testlog.log').read(), msg='Listobs fails to identify improper string')
        listobs(vis=dataset, observation='10')
        self.assertTrue(('The selected table has zero rows') in open('testlog.log').read(),
                        msg='Listobs fails to identify an empty table')

        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('ObservationID' in open('listobs.txt').read(), msg='There is no Observation information')

    def verbosecheck(self, dataset):
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, verbose=False)
        nonverb = ['name', 'station']
        self.assertTrue(all(x in open('testlog.log').read() for x in nonverb),
                        msg='non-verbose not showing proper info')
        listobs(vis=dataset, verbose=True)
        items = ['Name', 'Station', r'Diam.', r'Long.', r'Lat.', 'Offset', 'ITRF', 'Scan', 'FieldName', 'SpwIds']
        self.assertTrue(all(x in open('testlog.log').read() for x in items))

    def overwritecheck(self, dataset):
        listfile = "listobs.txt"
        try:
            listobs(vis=dataset, listfile=listfile)
        except Exception:
            self.fail()
        # test default value is overwrite=False
        with self.assertRaises(RuntimeError):
            listobs(vis=dataset, listfile=listfile)
        with self.assertRaises(RuntimeError):
            listobs(vis=dataset, listfile=listfile, overwrite=False)
        expec = _sha1it(listfile)
        try:
            listobs(vis=dataset, listfile=listfile, overwrite=True)
        except Exception:
            self.fail()
        got = _sha1it(listfile)
        self.assertTrue(got == expec)

    def avgIntervalcheck(self, dataset):
        listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Average Interval' in open('listobs.txt').read(),
                        msg='There is no average interval column in a MS')

    def unfcheck(self, dataset):
        try:
            listobs(vis=dataset, listunfl=True)
        except Exception:
            self.fail()
        listobs(vis=dataset, listfile='listobs1.txt', listunfl=True)
        self.assertTrue('nUnflRows' in open('listobs1.txt').read())

        listobs(vis=dataset, listfile='listobs2.txt', listunfl=False)
        self.assertFalse('nUnflRows' in open('listobs2.txt').read())


class test_listobs(listobs_test_base):

    @classmethod
    def setUpClass(cls):

        pass

    def setUp(self):
        self.res = None

        # From merged cases
        self.ms = ms

        if (not os.path.exists(msfile1)):
            shutil.copytree(msfile1Orig, msfile1, symlinks=True)
        if (not os.path.exists(msfile2)):
            shutil.copytree(msfile2Orig, msfile2, symlinks=True)

    def tearDown(self):
        # remove files and temp logs
        os.system('rm -rf ' + 'listobs*.txt')
        os.system('rm -rf testlog.log')
        casalog.setlogfile(str(logpath))
        # From merged cases
        self.ms.done()

    @classmethod
    def tearDownClass(cls):
        # remove all the generated data
        os.system('rm -rf genmms.mms')
        os.system('rm -rf gentimeavgms.ms')
        os.system('rm -rf gentimeavgmms.mms')
        os.system('rm -rf genmms.mms.flagversions')
        # From merged tests
        shutil.rmtree(msfile1, ignore_errors=True)
        shutil.rmtree(msfile2, ignore_errors=True)
        os.system('rm -f diff*')
        os.system('rm -f newobs*')
        if os.path.exists('CAS-5203.log'): os.remove('CAS-5203.log')

    # Test that Different ms are taken
    def test_wrongInp(self):
        with self.assertRaises(AssertionError):
            listobs(vis='foo.ms')

    def test_input(self):
        '''Listobs test: Check all if listobs can take ms, mms, time averaged ms, and time averaged mms'''
        # See if listobs can take all these forms of MS
        try:
            listobs(vis=mesSet)
        except Exception:
            self.fail('Fails to take MS file')
        try:
            listobs(vis=multiMesSet)
        except Exception:
            self.fail('Fails to take MMS file')
        try:
            listobs(vis=timeavg_ms)
        except Exception:
            self.fail('Fails to take time-averaged MS file')
        try:
            listobs(vis=timeavg_mms)
        except Exception:
            self.fail('Fails to take time-averaged MMS file')

    # Test the list file
    def test_fileMS(self):
        '''Listobs test: Check to see if list file is generated from a MS'''
        self.check_file_plus_dict(mesSet, 'listobs.txt')

    def test_fileMMS(self):
        '''Listobs test: Check to see if list file is generated from a MMS'''
        self.check_file_plus_dict(multiMesSet, 'listobs.txt')

    def test_fileTimeAvgMS(self):
        '''Listobs test: Check to see if list file is generated from a time-averaged MS'''
        self.check_file_plus_dict(timeavg_ms, 'listobs.txt')

    def test_fileTimeAvgMMS(self):
        '''Listobs test: Check to see if list file is generated from a time-averaged MMS'''
        self.check_file_plus_dict(timeavg_mms, 'listobs.txt')

    # Test the log file encoding
    def test_logreadMS(self):
        '''Listobs test: Check that the log file from a MS is human readable'''
        self.check_logread_plus_dict(mesSet)

    def test_logreadMMS(self):
        '''Listobs test: Check that the log file from a MMS is human readable'''
        self.check_logread_plus_dict(multiMesSet)

    def test_logreadTimeAvgMS(self):
        '''Listobs test: Check that the log file from a time-averaged MS is human readable'''
        self.check_logread_plus_dict(timeavg_ms)

    def test_logreadTimeAvgMMS(self):
        '''Listobs test: Check that the log file from a time-averaged MMS is human readable'''
        self.check_logread_plus_dict(timeavg_mms)

    # Test log file written to

    def test_logfileMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.check_logfilecontain_plus_dict(mesSet)

    def test_logfileMMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.check_logfilecontain_plus_dict(multiMesSet)

    def test_logfileTimeAvgMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.check_logfilecontain_plus_dict(timeavg_ms)

    def test_logfileTimeAvgMMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.check_logfilecontain_plus_dict(timeavg_mms)

    # Test one row displayed per scan

    def test_indivrow(self):
        '''Listobs test: Check if one row is displayed per scan'''
        ms.open(mesSet)
        self.res = ms.summary(True, listunfl=True)
        ms.close()
        # beginning and ending times should be different
        btime = self.res['scan_1']['0']['BeginTime']
        etime = self.res['scan_1']['0']['EndTime']
        self.assertNotEqual(btime, etime, msg='Beginning and Ending times should not be equal')

    # Test that no visibility information is given

    def test_novisMS(self):
        '''Listobs test: Check if there is any visibility information when using a MS'''
        self.novis(mesSet)

    def test_novisMMS(self):
        '''Listobs test: Check if there is any visibility information when using a MMS'''
        self.novis(multiMesSet)

    def test_novisTimeAvgMS(self):
        '''Listobs test: Check if there is any visibility information when using a time-averaged MS'''
        self.novis(timeavg_ms)

    def test_novisTimeAvgMMS(self):
        '''Listobs test: Check if there is any visibility information when using a time-averaged MMS'''
        self.novis(timeavg_mms)

    # Test the scan parameter

    def test_scanMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a MS'''
        self.scancheck(mesSet)

    def test_scanMMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a MMS'''
        self.scancheck(multiMesSet)

    def test_scanTimeAvgMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a time-averaged MS'''
        self.scancheck(timeavg_ms)

    def test_scanTimeAvgMMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a time-averaged MMS'''
        self.scancheck(timeavg_mms)

    # Test the field parameter

    def test_fieldMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a MS'''
        self.fieldcheck(mesSet)

    def test_fieldMMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a MMS'''
        self.fieldcheck(multiMesSet)

    def test_fieldTimeAvgMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a time-averaged MS'''
        self.fieldcheck(timeavg_ms)

    def test_fieldTimeAvgMMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a time-averaged MMS'''
        self.fieldcheck(timeavg_mms)

    # Test the spw parameter

    def test_spwMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a MS'''
        self.spwcheck(mesSet)

    def test_spwMMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a MMS'''
        self.spwcheck(multiMesSet)

    def test_spwTimeAvgMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a time-averaged MS'''
        self.spwcheck(timeavg_ms)

    def test_spwTimeAvgMMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a time-averaged MMS'''
        self.spwcheck(timeavg_mms)

    # Test the Correlation parameter

    def test_corrsMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a MS'''
        self.corrcheck(mesSet)

    def test_corrsMMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a MMS'''
        self.corrcheck(multiMesSet)

    def test_corrsTimeAvgMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a time-averaged MS'''
        self.corrcheck(timeavg_ms)

    def test_corrsTimeAvgMMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a time-averaged MMS'''
        self.corrcheck(timeavg_mms)

    # Test the antenna parameter

    def test_antMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a MS'''
        self.antcheck(mesSet)

    def test_antMMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a MMS'''
        self.antcheck(multiMesSet)

    def test_antTimeAvgMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a time-averaged MS'''
        self.antcheck(timeavg_ms)

    def test_antTimeAvgMMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a time-averaged MMS'''
        self.antcheck(timeavg_mms)

    # Test the selectdata parameter

    def test_selectdataMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a MS'''
        self.selectcheck(mesSet)

    def test_selectdataMMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a MMS'''
        self.selectcheck(multiMesSet)

    def test_selectdataTimeAvgMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a time-averaged MS'''
        self.selectcheck(timeavg_ms)

    def test_selectdataTimeAvgMMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a time-averaged MMS'''
        self.selectcheck(timeavg_mms)

    # Test the uvrange parameter

    def test_uvrangeMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a MS'''
        self.uvrangecheck(mesSet)

    def test_uvrangeMMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a MMS'''
        self.uvrangecheck(multiMesSet)

    def test_uvrangeTimeAvgMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a time-averaged MS'''
        self.uvrangecheck(timeavg_ms)

    def test_uvrangeTimeAvgMMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a time-averaged MMS'''
        self.uvrangecheck(timeavg_mms)

    # Test time range parameter

    def test_timerangeMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a MS'''
        self.timerangecheck(mesSet)

    def test_timerangeMMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a MMS'''
        self.timerangecheck(multiMesSet)

    def test_timerangeTimeAvgMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a time-averaged MS'''
        self.timerangecheck(timeavg_ms)

    def test_timerangeTimeAvgMMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a time-averaged MMS'''
        self.timerangecheck(timeavg_mms)

    # Test the intent parameter

    def test_intentMS(self):
        '''Listobs test: Check to see that Intent info exists for a MS and accepts correct inputs'''
        self.intentcheck(mesSet)

    def test_intentMMS(self):
        '''Listobs test: Check to see that Intent info exists for a MMS and accepts correct inputs'''
        self.intentcheck(multiMesSet)

    def test_intentTimeAvgMS(self):
        '''Listobs test: Check to see that Intent info exists for a time-averaged MS and accepts correct inputs'''
        self.intentcheck(timeavg_ms)

    def test_intentTimeAvgMMS(self):
        '''Listobs test: Check to see that Intent info exists for a time-averaged MMS and accepts correct inputs'''
        self.intentcheck(timeavg_mms)

    # Test the array parameter

    def test_arrayMS(self):
        '''Listobs test: Check for the existence of the array parameter in a MS and accepts proper inputs'''
        self.arraycheck(mesSet)

    def test_arrayMMS(self):
        '''Listobs test: Check for the existence of the array parameter in a MMS and accepts proper inputs'''
        self.arraycheck(multiMesSet)

    def test_arrayTimeAvgMS(self):
        '''Listobs test: Check for the existence of the array parameter in a time-averaged MS and accepts proper inputs'''
        self.arraycheck(timeavg_ms)

    def test_arrayTimeAvgMMS(self):
        '''Listobs test: Check for the existence of the array parameter in a time-averaged MMS and accepts proper inputs'''
        self.arraycheck(timeavg_mms)

    # Test the observation parameter

    def test_observationMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a MS and check for proper inputs'''
        self.obscheck(mesSet)

    def test_observationMMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a MMS and check for proper inputs'''
        self.obscheck(multiMesSet)

    def test_observationTimeAvgMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a time-averaged MS and check for proper inputs'''
        self.obscheck(timeavg_ms)

    def test_observationTimeAvgMMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a time-averaged MMS and check for proper inputs'''
        self.obscheck(timeavg_mms)

    # Test the Verbose parameter

    def test_verboseMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a MS'''
        self.verbosecheck(mesSet)

    def test_verboseMMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a MMS'''
        self.verbosecheck(multiMesSet)

    def test_verboseTimeAvgMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a time-averaged MS'''
        self.verbosecheck(timeavg_ms)

    def test_verboseTimeAvgMMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a time-averaged MMS'''
        self.verbosecheck(timeavg_mms)

    # Test the overwrite param

    def test_overwriteMS(self):
        """Test overwrite parameter - CAS-5203: test for MS"""
        self.overwritecheck(mesSet)

    def test_overwriteMMS(self):
        """Test overwrite parameter - CAS-5203: test for MMS"""
        self.overwritecheck(multiMesSet)

    def test_overwriteTimeAvgMS(self):
        """Test overwrite parameter - CAS-5203: test for time-averaged MS"""
        self.overwritecheck(timeavg_ms)

    def test_overwriteTimeAvgMMS(self):
        """Test overwrite parameter - CAS-5203: test for time-averaged MMS"""
        self.overwritecheck(timeavg_mms)

    # Test the inf loop bug CAS-6733

    def test_CAS_6733(self):
        """Verify listobs runs to completion on data set in CAS-6733. This was an infinite loop bugfix"""
        vis = msfile3
        try:
            listobs(vis=vis)
        except Exception:
            self.fail()

    # Test average interval

    def test_avgIntervalMS(self):
        '''Listobs test: Check that the Int (s) column exists for a MS'''
        self.avgIntervalcheck(mesSet)

    def test_avgIntervalMMS(self):
        '''Listobs test: Check that the Int (s) column exists for a MMS'''
        self.avgIntervalcheck(multiMesSet)

    def test_avgIntervalTimeAvgMS(self):
        '''Listobs test: Check that the Int (s) column exists for a time-averaged MS'''
        self.avgIntervalcheck(timeavg_ms)

    def test_avgIntervalTimeAvgMMS(self):
        '''Listobs test: Check that the Int (s) column exists for a time-averaged MMS'''
        self.avgIntervalcheck(timeavg_mms)

    # Test list unflagged parameter

    def test_listunflMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a MS'''
        self.unfcheck(mesSet)

    def test_listunflMMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a MMS'''
        self.unfcheck(multiMesSet)

    def test_listunflTimeAvgMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a time-averaged MS'''
        self.unfcheck(timeavg_ms)

    def test_listunflTimeAvgMMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a time-averaged MMS'''
        self.unfcheck(timeavg_ms)

    # Test Dysco storage manager

    def test_lofarDysco(self):
        '''Verify that we can read a compressed MS that uses the Dysco storage manager'''
        try:
            listobs(vis=msfile4)
        except Exception:
            self.fail()

    # Start of merged cases from test_listobs

    def test_longFieldName(self):
        '''Listobs 2: CSV-591. Check if long field names are fully displayed'''
        self.ms.open(msfile1)
        res = self.ms.summary(True, listunfl=True)
        self.ms.close()
        name = res['field_0']['name']
        self.assertFalse(name.__contains__('*'), "Field name contains a *")
        name = res['scan_7']['0']['FieldName']
        self.assertFalse(name.__contains__('*'), "Field name contains a *")

    def test_almaOneRow(self):
        '''Listobs 3: CAS-2751. Check that ALMA MS displays one row per scan'''
        self.ms.open(msfile2)
        res = self.ms.summary(True, listunfl=True)
        self.ms.close()
        # Begin and end times should be different
        btime = res['scan_1']['0']['BeginTime']
        etime = res['scan_1']['0']['EndTime']
        self.assertNotEqual(btime, etime, "Begin and End times of scan=1 should not be equal")

        # Only one row of scan=1 should be printed
        output = 'listobs4.txt'
        out = "newobs4.txt"
        reference = reffile + '_almaOneRow'
        diff = "difflistobs4"

        listobs(vis=msfile2, verbose=True, listfile=output, listunfl=True)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff))

    def test_nonVerboseFileSave(self):
        '''Listobs 4: Save on a file, verbose=False'''
        output = 'listobs5.txt'
        out = "newobs5.txt"
        reference = reffile + '_nonVerboseFileSave'
        diff1 = "diff1listobs5"
        diff2 = "diff2listobs5"

        #        # Run it twice to check for the precision change
        self.res = listobs(vis=msfile1, verbose=False, listfile=output, listunfl=True)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff1)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different in first run. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff1))

        os.system('rm -rf ' + output + " " + out)
        self.res = listobs(vis=msfile1, verbose=False, listfile=output, listunfl=True)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff2)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different in second run. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff2))

    def test_verboseFileSave(self):
        '''Listobs 5: Save on a file, verbose=True'''
        output = 'listobs6.txt'
        out = "newobs6.txt"
        diff = "difflistobs6"
        reference = reffile + '_verboseFileSave'
        self.res = listobs(vis=msfile1, listfile=output, verbose=True, listunfl=True)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff))

    def test_scanSelectionParam(self):
        '''Listobs 6: test scan selection parameters'''
        output = "listobs7.txt"
        out = "newobs7.txt"
        diff = "difflistobs7"
        reference = reffile + '_scanSelectionParam'
        self.res = listobs(vis=msfile1, scan='2', listfile=output, verbose=True, listunfl=True)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff))

    def test_antennaSelectionParam(self):
        '''Listobs 7: test antenna selection parameters'''
        output = "listobs8.txt"
        out = "newobs8.txt"
        diff = "difflistobs8"
        reference = reffile + '_antennaSelectionParam'
        self.res = listobs(vis=msfile1, antenna='3&&4', listfile=output, verbose=True, listunfl=True)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff))

    def test_ephem(self):
        '''ephemeris objects'''
        output = "listobs9.txt"
        out = "newobs9.txt"
        diff = "difflistobs9"
        reference = reffile + '_ephem'
        self.res = listobs(vis=os.path.join(datapath, nep), listfile=output, verbose=True, listunfl=False)
        #        # Remove the name of the MS from output before comparison
        os.system("sed '1,3d' " + output + ' > ' + out)
        os.system("diff " + reference + " " + out + " > " + diff)
        self.assertTrue(lt.compare(out, reference),
                        'New and reference files are different. %s != %s. '
                        'See the diff file %s.' % (out, reference, diff))

    def test_overwrite(self):
        """Test overwrite parameter - CAS-5203"""
        listfile = "CAS-5203.log"
        listobs(vis=msfile1, listfile=listfile)

        # test default value is overwrite=False
        with self.assertRaises(RuntimeError):
            listobs(vis=msfile1, listfile=listfile)
        with self.assertRaises(RuntimeError):
            listobs(vis=msfile1, listfile=listfile, overwrite=False)

        expec = _sha1it(listfile)
        listobs(vis=msfile1, listfile=listfile, overwrite=True)
        got = _sha1it(listfile)
        self.assertTrue(got == expec)

if __name__ == '__main__':
    unittest.main()

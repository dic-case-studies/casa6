##########################################################################
#
# Copyright (C) 2019 ESO (in the framework of the ALMA collaboration)
# Copyright (C) 2019 Associated Universities, Inc. Washington DC, USA.
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
# CAS-12834 - Very minimal test to have some coverage of task listcal.
#             This is more a TODO than a test.
#
##########################################################################
CASA6=False
try:
    import casatools
    from casatasks import listcal #, casalog
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

import os
import shutil
import unittest

reg_unittest_datap = 'unittest/listcal/'
if CASA6:
    datapath = casatools.ctsys.resolve(reg_unittest_datap)
else:
    datapath = os.path.join(os.path.join(os.environ.get('CASAPATH').split()[0],
                                         'casatestdata'), reg_unittest_datap)

# This is for tests that check what the parameter validator does when parameters are
# given wrong types - these don't exercise the task but the parameter validator!
if CASA6:
    validator_exc_type = AssertionError
else:
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
    validator_exc_type = RuntimeError

def searchInFile(filename, searched):
    """
    Search for a substring in a file
    """
    with open(filename) as searched_file:
        for line in searched_file:
            if searched in line:
                return True
        return False


def getLine(filename, lineNum):

    counter = 0
    with open(filename) as fout:
        for line in fout:
            if counter == lineNum:
                return line
            else:
                counter += 1


class test_listcal_minimal(unittest.TestCase):
    """
    Most basic ways of calling listcal, with an existing small dataset from the
    'unittest' input data that is used in the bandpass and gaincal tests.
    This is just a start, to have some tests for the task listcal.
    """

    @classmethod
    def setUpClass(cls):
        if not CASA6:
            default(listcal)

        vis_name = 'ngc5921.ms'
        caltable_name = 'ngc5921.ref1a.gcal'
        shutil.copytree(os.path.join(datapath,vis_name),vis_name)
        shutil.copytree(os.path.join(datapath,caltable_name),caltable_name)
        cls._vis = vis_name
        cls._caltable = caltable_name
        cls._listfile = 'listcal_listfile.txt'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._vis)
        shutil.rmtree(cls._caltable)

    def tearDown(self):
        if os.path.isfile(self._listfile):
            os.remove(self._listfile)

    def test_small_gcal_nosel(self):
        """
        Test listcal with MS + gcal caltable, without any selection.
        """
        listcal(vis=self._vis, caltable=self._caltable, listfile=self._listfile)
        self.assertTrue(os.path.isfile(self._listfile))
        size = os.path.getsize(self._listfile)
        self.assertGreaterEqual(size, 10000)
        self.assertLessEqual(size, 12000)

    def test_small_gcal_sel_ant(self):
        """
        Test to have some additional tests with selections for the moment. Too basic.
        """
        listcal(vis=self._vis, caltable=self._caltable, antenna='VA02',
                listfile=self._listfile)
        self.assertTrue(os.path.isfile(self._listfile))
        size = os.path.getsize(self._listfile)
        self.assertGreaterEqual(size, 1000)
        self.assertLessEqual(size, 1200)

    def test_small_gcal_forget_caltable(self):
        """
        Test proper error when listcal with MS + forget caltable.
        """
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(validator_exc_type):
                listcal(vis=self._vis, listfile=self._listfile)
        else:
                listcal(vis=self._vis, listfile=self._listfile)

        self.assertFalse(os.path.isfile(self._listfile))

    def test_small_gcal_forget_ms(self):
        """
        Test proper error when listcal caltable + forget MS.
        """
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(validator_exc_type):
                listcal(caltable=self._caltable, listfile=self._listfile)
        else:
            listcal(caltable=self._caltable, listfile=self._listfile)                

        self.assertFalse(os.path.isfile(self._listfile))

    def test_small_gcal_wrong_antenna(self):
        """
        Test that gives a wrong selection (antenna).
        """
        with self.assertRaises(RuntimeError):
            listcal(vis=self._vis, caltable=self._caltable, antenna='inexistent-no-no')
        self.assertFalse(os.path.isfile(self._listfile))


class listcal_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not CASA6:
            default(listcal)

        vis_name = 'ngc5921.ms'
        caltable_name = 'ngc5921.ref1a.gcal'
        reffile_name = 'reflistcal.txt'

        gainvis_name = 'gaincaltest2.ms'
        gaincaltable_name = 'gaincaltest2.ms.G0'
        shutil.copytree(os.path.join(datapath,vis_name),vis_name)
        shutil.copytree(os.path.join(datapath,caltable_name),caltable_name)
        shutil.copytree(os.path.join(datapath, gainvis_name), gainvis_name)
        shutil.copytree(os.path.join(datapath, gaincaltable_name), gaincaltable_name)
        shutil.copyfile(os.path.join(datapath, reffile_name), reffile_name)
        cls._vis = vis_name
        cls._caltable = caltable_name
        cls._gainvis = gainvis_name
        cls._gaincaltable = gaincaltable_name
        cls._listfile = 'listcal_listfile.txt'
        cls._reffile = reffile_name

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._vis)
        shutil.rmtree(cls._caltable)
        shutil.rmtree(cls._gainvis)
        shutil.rmtree(cls._gaincaltable)
        os.remove(cls._reffile)

    def tearDown(self):
        if os.path.isfile(self._listfile):
            os.remove(self._listfile)

    def test_fieldSelection(self):
        """ Checks that specific fields can be selected """
        listcal(vis=self._vis, caltable=self._caltable, field='N5921_2', listfile=self._listfile)

        self.assertFalse(searchInFile(self._listfile, '1331+30500002_0'), msg='Found fields that are not selected')
        self.assertTrue(searchInFile(self._listfile, 'N5921_2'), msg='did not find selected field')

    def test_antennaSelection(self):
        """ Checks that specific antenna can be selected """
        listcal(vis=self._vis, caltable=self._caltable, antenna='VA01', listfile=self._listfile)

        self.assertFalse(searchInFile(self._listfile, 'VA02'), msg='Found anntennas that were not selected')
        self.assertTrue(searchInFile(self._listfile, 'VA01'), msg='The selected antenna was not found')

    def test_spwSelection(self):
        """ Check that specific spectral windows can be selected """
        listcal(vis=self._gainvis, caltable=self._gaincaltable, spw='1', listfile=self._listfile)

        self.assertFalse(searchInFile(self._listfile, 'SpwID = 0'), msg='Found spw that was not selected')
        self.assertTrue(searchInFile(self._listfile, 'SpwID = 1'), msg='Did not find the selected spw')

    def test_headerLayout(self):
        """ Check that the header is layed out in the form described in the documentation """
        listcal(vis=self._gainvis, caltable=self._gaincaltable, listfile=self._listfile)

        header = getLine(self._listfile, 0)

        self.assertTrue("SpwID" in header)
        self.assertTrue("Date" in header)
        self.assertTrue("CalTable" in header)
        self.assertTrue("MS name" in header)

    def test_spectralWindowHeader(self):
        """ Check that there is spectral window, date, caltable name, and ms name in the header """
        listcal(vis=self._gainvis, caltable=self._gaincaltable, listfile=self._listfile)

        header = getLine(self._listfile, 14)

        self.assertTrue("SpwID" in header)
        self.assertTrue("Date" in header)
        self.assertTrue("CalTable" in header)
        self.assertTrue("MS name" in header)

    def test_headerColumnLabels(self):
        """ Check that the data columns are labeled as described """
        listcal(vis=self._gainvis, caltable=self._gaincaltable, listfile=self._listfile)

        header = getLine(self._listfile, 3)
        self.assertTrue(header.startswith("Time       Field      Chn| Amp    Phs  F  Amp    Phs  F| Amp    Phs  F  Amp    Phs  F| Amp    Phs  F  Amp    Phs  F| Amp    Phs  F  Amp    Phs  F|"))

    def test_compareValues(self):
        """ Compare the values in the output to a reference output for the same table """
        listcal(vis=self._gainvis, caltable=self._gaincaltable, listfile=self._listfile)

        with open(self._listfile) as file1:
            with open(self._reffile) as file2:
                diff = set(file1).difference(file2)

        # There are differences in the file path when running the command
        self.assertTrue(len(diff) == 4)

def suite():
    return [test_listcal_minimal, listcal_test]

if __name__ == '__main__':
    unittest.main()

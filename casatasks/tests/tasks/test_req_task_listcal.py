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
    # import casatools
    from casatasks import listcal #, casalog
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

import os
import shutil
import unittest

reg_unittest_datap = 'regression/unittest'
if CASA6:
    datapath = casatools.ctsys.resolve(reg_unittest_datap)
else:
    datapath = os.path.join(os.path.join(os.environ.get('CASAPATH').split()[0],
                                         'data'), reg_unittest_datap)

class test_listcal_minimal(unittest.TestCase):
    """
    Most basic ways of calling listcal, with an existing small dataset from the
    'regression/unittest' input data that is used in the bandpass and gaincal tests.
    This is just a start, to have some tests for the task listcal.
    """

    @classmethod
    def setUpClass(cls):
        if not CASA6:
            default(listcal)

        vis_name = 'ngc5921.ms'
        caltable_name = 'ngc5921.ref1a.gcal'
        os.system('cp -RL {} .'.format(os.path.join(datapath, 'gaincal/{}'.
                                                    format(vis_name))))
        os.system('cp -RL {} .'.format(os.path.join(datapath, 'gaincal/{}'.
                                                    format(caltable_name))))
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
        if CASA6:
            with self.assertRaises(AssertionError):
                listcal(vis=self._vis, listfile=self._listfile)
        else:
                listcal(vis=self._vis, listfile=self._listfile)

        self.assertFalse(os.path.isfile(self._listfile))

    def test_small_gcal_forget_ms(self):
        """
        Test proper error when listcal caltable + forget MS.
        """
        if CASA6:
            with self.assertRaises(AssertionError):
                listcal(caltable=self._caltable, listfile=self._listfile)
        else:
            listcal(caltable=self._caltable, listfile=self._listfile)                

        self.assertFalse(os.path.isfile(self._listfile))

    def test_small_gcal_wrong_antenna(self):
        """
        Test that gives a wrong selection (antenna).
        """
        listcal(vis=self._vis, caltable=self._caltable, antenna='inexistent-no-no')
        self.assertFalse(os.path.isfile(self._listfile))


def suite():
    return [test_listcal_minimal]

if __name__ == '__main__':
    unittest.main()

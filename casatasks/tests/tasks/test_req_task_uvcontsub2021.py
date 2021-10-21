##########################################################################
# test_req_task_phaseshift.py
#
# Copyright (C) 2021 European Southern Obervatory, ALMA partnership
# Copyright (C) 2021 Associated Universities, Inc. Washington DC, USA.
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
# Test initially added in CAS-13631
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.manipulation.uvcontsub2021.html
# And the documents found in the "Development"/"Visibility Manipulation" folder of the CASA
# google drive:
# https://drive.google.com/drive/u/1/folders/1ttYI8Xcgfa-e1Dk0f8kzrv1bz7fIqylm
#
##########################################################################
import numpy as np
import os
import shutil
import unittest

CASA6 = False
try:
    from casatools import table, ctsys
    from casatasks import uvcontsub2021
    ctsys_resolve = ctsys.resolve
    CASA6 = True
except ImportError:
    from tasks import uvcontsub2021
    from taskinit import tbtool as table

    dataroot = os.path.join(
        os.environ.get('CASAPATH').split()[0], 'casatestdata'
    )

    def ctsys_resolve(mypath):
        return os.path.join(dataroot, mypath)

datadir = os.path.join('unittest', 'uvcontsub')
ms_simple = 'known0.ms'
datapath_simple = ctsys_resolve(os.path.join(datadir, ms_simple))

ms_alma = 'uid___X02_X3d737_X1_01_small.ms'
datapath_alma = ctsys_resolve(os.path.join('measurementset', 'alma', ms_alma))


class uvcontsub2021_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        shutil.copytree(datapath_simple, ms_simple)
        shutil.copytree(datapath_alma, ms_alma)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(ms_simple)

    def setUp(self):
        # Input MS is always strictly read-only, one copy in setUpClass is enough
        # Default output name for simple tests
        self.output = 'test_uvcs_output.ms'

    def tearDown(self):
        if os.path.exists(self.output):
            shutil.rmtree(self.output)

    def _check_rows(self, vis, col_name, expected_rows, expected_val=None):
        tbt = table()
        try:
            tbt.open(vis)
            col = tbt.getcol(col_name)
            nrows = tbt.nrows()  # or col.shape[-1]
            self.assertEqual(nrows, expected_rows, 'Number of rows different from expected')
            if expected_val:
                self.assertTrue(np.all(col == expected_val), "Column '{}' values different "
                                "from expected. Expected: {}. Column values: {}".
                                format(col_name, expected_val, col))
        finally:
            tbt.done()

    def _check_return(self, res):
        self.assertEqual(res, {})

    def test_makes_output_ms_data(self):
        """
        Check that in a simple command the input MS is taken and an output MS
        is created and has a data column
        """

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output)
        self._check_return(res)
        self.assertTrue(os.path.exists(self.output))
        self._check_rows(self.output, 'DATA', 340)

        # check also no-overwrite of existing MS
        with self.assertRaises(ValueError):
            uvcontsub2021(vis=ms_simple, outputvis=self.output)

    def test_select_field(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='2')
        self._check_return(res)
        self._check_rows(self.output, 'FIELD_ID', 120, 2)

    def test_select_spw(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, spw='1')
        self._check_return(res)
        self._check_rows(self.output, 'DATA_DESC_ID', 810, 1)

    def test_select_scan(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, scan='2')
        self._check_return(res)
        self._check_rows(self.output, 'SCAN_NUMBER', 360, 2)

    def test_select_intent(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, intent='*AMPLI*')
        self._check_return(res)
        self._check_rows(self.output, 'SCAN_NUMBER', 360, 2)

    def test_select_array(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, array='0')
        self._check_return(res)
        self._check_rows(self.output, 'ARRAY_ID', 1080, 0)

    def test_select_observation(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, observation='0')
        self._check_return(res)
        self._check_rows(self.output, 'OBSERVATION_ID', 1080, 0)

    def test_datacolumn(self):
        """Check different datacolumn choices and that results make sense
        depending oninput MS"""

        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='CORRECTED')

        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='MODEL')

        if CASA6:
            with self.assertRaises(AssertionError):
                res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='bogus')
        else:
            with self.assertRaises(RuntimeError):
                res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='bogus')

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='DATA')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)

    def test_fitspw_empty(self):
        """Check that fitspw works. When empty, fit all channels in all SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspw='')
        self._check_return(res)
        # TODO: better checks and most likely a better input ms (keep small)
        self._check_rows(self.output, 'DATA', 340)

    def test_fitspw_spws(self):
        """Check that fitspw works. When selecting some SPWs, fit all channels
        in those SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspw='0')
        self._check_return(res)
        # TODO: better checks and most likely a better input ms (keep small)
        self._check_rows(self.output, 'DATA', 340)

    def test_fitswp_channels(self):
        """Check that fitspw works. When selecting some channels in some SPWs,
        fit those channels in those SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspw='0:5~19')
        self._check_return(res)
        # TODO: better checks and most likely a better input ms (keep small)
        self._check_rows(self.output, 'DATA', 340)

    def test_fitspw_multifield(self):
        """Check that fitspw works. Different fitspw strings for different fields"""

        with self.assertRaises(RuntimeError):
            # TODO: support this properly
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=['2', '0:100~500;600~910;1215~1678;1810~1903'])
            self._check_return(res)
            # TODO: better checks and most likely a better input ms (keep small)
            self._check_rows(self.output, 'DATA', 1080)

    def test_fitmethod_gsl(self):
        """Check that methods work - gsl"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitmethod='gsl')
        self._check_return(res)
        # TODO: better checks / overlaps with numerical tests
        self._check_rows(self.output, 'DATA', 340)

    def test_fitmethod_casacore(self):
        """Check that methods work - casacore"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitmethod='casacore')
        self._check_return(res)
        # TODO: better checks / overlaps with numerical tests
        self._check_rows(self.output, 'DATA', 340)

    def test_fitorder1(self):
        """ Check different fit orders (0, 1, 2) work"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=1)
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)

    def test_fitorder2(self):
        """ Check different fit orders (0, 1, 2) work"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=2)
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)

    def test_writemodel(self):
        """ Check the model column is added to the output MS and its values maatch"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, writemodel=True)
        self._check_return(res)
        # TODO: get model column, check input/DATA ~= output/DATA+MODEL
        with self.assertRaises(RuntimeError):
            self._check_rows(self.output, 'MODEL', 1080)


class uvcontsub2021_numerical_verification_test(unittest.TestCase):
    """ To be defined - CAS-13632 """

    def test_sim_single_source_snr_better(self):
        pass

    def test_sim_single_source_snr_worse(self):
        pass


def suite():
    return [uvcontsub2021_test,
            uvcontsub2021_numerical_verification_test]


if __name__ == '__main__':
    unittest.main()

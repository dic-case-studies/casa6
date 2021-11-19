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

from casatools import table, ctsys
from casatasks import uvcontsub2021
ctsys_resolve = ctsys.resolve

datadir = os.path.join('unittest', 'uvcontsub')
ms_simple = 'known0.ms'
datapath_simple = ctsys.resolve(os.path.join(datadir, ms_simple))

ms_alma = 'uid___X02_X3d737_X1_01_small.ms'
datapath_alma = ctsys.resolve(os.path.join('measurementset', 'alma', ms_alma))


class uvcontsub2021_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        shutil.copytree(datapath_simple, ms_simple)
        shutil.copytree(datapath_alma, ms_alma)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(ms_simple)
        shutil.rmtree(ms_alma)

    def setUp(self):
        # Input MS is always strictly read-only, one copy in setUpClass is enough
        # Default output name for simple tests
        self.output = 'test_uvcs_output.ms'

    def tearDown(self):
        if os.path.exists(self.output):
            shutil.rmtree(self.output)

    def _check_rows(self, vis, col_name, expected_rows, expected_val=None):
        """
        Meant to check rows of a column from an output MS produced by
        uvcontsub2021. Uses unittest asserts to verify conditions.

        :param vis: MS to check
        :param col_name: name of column to check (MODEL_DATA, FIELD, etc.)
        :param expected_rows: number of rows the column must have
        :param expected_vals: for simple cases where the same value is expected in every row,
                              ensure all rows have this value
        """
        # TODO: add simple checks of data column values (mean or similar). Overlaps partially
        # with numerical tests. Like a parameter 'expected_mean=None' to compare approx.
        # Perhaps a different function specific to data column checks
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

    def _check_data_stats(self, vis, exp_mean, exp_median, exp_min, exp_max,
                          col_name='DATA'):
        """
        WIP - check basic stats of data column, for now just to prevent uninteded changes
        """
        tbt = table()
        try:
            tbt.open(vis)
            col = tbt.getcol(col_name)

            nans = np.isnan(np.sum(col))
            self.assertFalse(nans)
            # ALMA test datasets have large numbers of 0s
            # zeros_count = np.count_nonzero(col==0)
            # self.assertEqual(0, zeros_count)
            places = 5
            if exp_mean is not None:
                dmean = np.mean(col)
                print('Mean: {}'.format(dmean))
                self.assertAlmostEqual(dmean, exp_mean, places=places)
            if exp_median is not None:
                dmedian = np.median(col)
                print('Median: {}'.format(dmedian))
                self.assertAlmostEqual(dmedian, exp_median, places=places)
            if exp_min is not None:
                dmin = col.min()
                print('Min: {}'.format(dmin))
                self.assertAlmostEqual(dmin, exp_min, places=places)
            if exp_max is not None:
                dmax = col.max()
                print('Max: {}'.format(dmax))
                self.assertAlmostEqual(dmax, exp_max, places=places)
        finally:
            tbt.done()

    def test_makes_output_ms_data(self):
        """
        Check that in a simple command the input MS is taken and an output MS
        is created and has a data column
        """

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output)
        self._check_return(res)
        self.assertTrue(os.path.exists(self.output))
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

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

        with self.assertRaises(AssertionError):
            res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='bogus')

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='DATA')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))
        # TODO need a 'datacolumn' test using CORRECTED - use another ms_alma or someth else

    def test_fitspw_empty(self):
        """Check that fitspw works. When empty, fit all channels in all SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspw='')
        self._check_return(res)
        # TODO: better checks and most likely a better input ms (prob ~pl-unittest)
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitspw_spws(self):
        """Check that fitspw works. When selecting some SPWs, fit all channels
        in those SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspw='0')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitspw_channels(self):
        """Check that fitspw works. When selecting some channels in some SPWs,
        fit those channels in those SPWs (like example 2 from task page)"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspw='0:5~19')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, (-18.5249996-7.05000019j),
                               (-19.7999992-13j), (-68-17.6000004j), (28+26.3999996j))

    def test_fitspw_multifield(self):
        """Check that fitspw works. Different fitspw strings for different fields
        (like example 4 from task page)"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                            fitspw=[
                                ['1', 'NONE'],
                                ['2', '0:100~500;600~910;1215~1678;1810~1903'],
                                ['3', '0:100~1903']
                            ])
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 1080)
        self._check_data_stats(self.output, (-0.0123391838-6.38635834e-06j),
                               0, (-0.655234873+0j), (2.09603309+0j))

    def test_fitspw_multifield_blocks(self):
        """Check that fitspw works. Different fitspw strings for different fields
        but giving list of fields for a same fitspw"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                            fitspw=[
                                ['1, 2', '0:100~500;600~910;1215~1678;1810~1903'],
                                ['3', '0:100~1903']
                            ])
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 1080)
        self._check_data_stats(self.output, (-0.0123391838-6.38635834e-06j),
                               0, (-0.655234873+0j), (2.09603309+0j))

    def test_fitspw_multifield_wrong_field(self):
        """Check that wrong fitspw lists produce an exception"""

        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=[
                                    ['99', '0:100~500;600~910;1215~1678;1810~1903']
                                ])
            self._check_return(res)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=[
                                    ['-2', '0:100~500']
                                ])
            self._check_return(res)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=[
                                    ['4', '0,1']
                                ])
            self._check_return(res)

    def test_fitspw_multifield_wrong_format(self):
        """Check that fitspw works. Different fitspw strings for different fields
        (like example 4 from task page)"""

        with self.assertRaises(RuntimeError):
            # Wrong number of elements (not a list of pairs) - in third line
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=[
                                    ['1', 'NONE'],
                                    ['2', '0:100~500;600~910;1215~1678;1810~1903'],
                                    ['3', '4', '0:100~1903']
                                ])

            self._check_return(res)
            self._check_rows(self.output, 'DATA', 1080)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            # Wrong field
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=[
                                    ['1', 'NONE'],
                                    ['bla_fail', '0:600~910;1215~1678'],
                                    ['3', '0:100~1903']
                                ])

            self._check_return(res)
            self._check_rows(self.output, 'DATA', 1080)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            # Repeated indices
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspw=[
                                    ['1, 2', 'NONE'],
                                    ['2', '0:100~500;600~910;1215~1678;1810~1903'],
                                ])
            self._check_return(res)
            self._check_rows(self.output, 'DATA', 1080)

    def test_fitspw_separate_fields(self):
        """Check that fitspw works. Different fitspw strings for different
        fields, and each field to a different output MS (like example 3 from
        task page)"""

        res_f1 = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='1',
                               fitspw='0:100~500;600~910;1215~1678;1810~1903')
        self._check_return(res_f1)
        self._check_rows(self.output, 'DATA', 360)
        self._check_data_stats(self.output, (-0.013026415-2.0328553e-05j),
                               0j, (-0.633247495+0j), (2.01062560+0j))

        shutil.rmtree(self.output)
        res_f2 = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='2',
                               fitspw='0:100~1303')
        self._check_return(res_f2)
        self._check_rows(self.output, 'DATA', 120)
        self._check_data_stats(self.output, (-0.0140258259+5.83529241e-06j),
                               0j, (-0.588652134+0j), (1.83768487+0j))

    def test_fitspw_spws_sel(self):
        """Check the use of spw selection and fitspw together"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, spw='0', fitspw='0')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitspw_spws_sel_wrong(self):
        """Check that if the spw selection and fitspw are not compatible, an exception
        is produced"""

        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_simple, outputvis=self.output, spw='3', fitspw='0')
            self._check_return(res)

    def test_fitmethod_gsl(self):
        """Check that methods work - gsl"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitmethod='gsl')
        self._check_return(res)
        # TODO: better checks / but overlaps with numerical tests
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitmethod_casacore(self):
        """Check that methods work - casacore"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitmethod='casacore')
        self._check_return(res)
        # TODO: better checks / but overlaps with numerical tests
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitorder1(self):
        """ Check different fit orders (0, 1, 2) work (like example 1 from task page)"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=1,
                            fitspw='0:2~20')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)

    def test_fitorder2(self):
        """ Check different fit orders (0, 1, 2) work"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=2,
                            fitspw='0:2~20')
        self._check_return(res)
        self._check_rows(self.output, 'DATA', 340)

    def test_writemodel(self):
        """ Check the model column is added to the output MS and its values match
        (like example 5 from task page)"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, writemodel=True)
        self._check_return(res)
        # TODO: get model column, check input/DATA ~= output/DATA+MODEL
        with self.assertRaises(RuntimeError):
            self._check_rows(self.output, 'MODEL_DATA', 340)


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

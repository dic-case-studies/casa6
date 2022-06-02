##########################################################################
# test_task_uvcontsub2021.py
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
# See also requirements as listed here:
# https://open-confluence.nrao.edu/display/CASA/uvcontsub2021
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

# SPW 1 of this dataset has 1 channel
ms_alma = 'uid___X02_X3d737_X1_01_small.ms'
datapath_alma = ctsys.resolve(os.path.join(datadir, ms_alma))

# MS for tests that use CORRECTED_DATA
# Beware: this is all flagged!
ms_corr = 'uid___A002_X71a45c_X1d24.ms.split'
datapath_corr = ctsys.resolve(os.path.join(datadir, ms_corr))

# Another MS for tests that use CORRECTED_DATA
ms_papersky = 'papersky_standard.ms'
datapath_papersky = ctsys.resolve(os.path.join(datadir, ms_papersky))

# Mixed polarizations, from CAS-12283. This MS has ~60 SPWs with very mixed pols
ms_mixed_pols = 'split_ddid_mixedpol_CAS-12283.ms'
datapath_mixed_pols = ctsys.resolve(os.path.join(datadir, ms_mixed_pols))


class uvcontsub2021_test_base(unittest.TestCase):
    """
    Base class to share utility functions between the test classes below.
    """

    def _check_data_stats(self, vis, exp_mean, exp_median, exp_min, exp_max,
                          col_name='DATA', places=5, fields=None, ddis=None):
        """
        Check basic stats of data column. Can work on the whole column or
        on slices of it (fields, ddis)

        This implements simple checks of expected values for the mean,
        median, min and max across all visibilities (all polarization,
        channels, baselins and times).  Overlaps partially with
        numerical tests (class below) and could be relaxed to apply
        less strict asserts.

        :param vis: MS to check
        :param exp_mean: expected mean value, across all visibilities
        :param exp_median: expected median value, across all visibilities
        :param exp_min: expected minimum value, across all visibilities
        :param exp_max: expected maximum value, across all visibilities
        :param col_name: name of data column with visibilities
        :param places: places/digits to check for "almost" equality
        :param fields: field or list of fields (FIELD_IDs) to include in stats
        :param ddis: DDI or list of DDIs to include in stats

        """
        def to_csv(ids):
            """List of IDs to string of comma-separated IDs"""
            if isinstance(ids, list):
                return ','.join([str(item) for item in ids])
            else:
                return ids

        def get_col(vis, col_name, fields, ddis):
            """ Returns the column data in an array"""
            tbt = table()
            try:
                tbt.open(vis)

                if fields is None and ddis is None:
                    col = tbt.getcol(col_name)
                else:
                    query_str = ''
                    if fields is not None:
                        query_str += 'FIELD_ID in [{}]'.format(to_csv(fields))
                    if ddis is not None:
                        if not query_str:
                            query_str = 'DATA_DESC_ID in [{}]'.format(to_csv(ddis))
                        else:
                            query_str += ' AND DATA_DESC_ID in [{}]'.format(to_csv(ddis))
                    try:
                        query_col = tbt.query(query_str, columns=col_name, style='python')
                        col = query_col.getcol(col_name)
                    finally:
                        query_col.done()

                if len(col) == 0:
                    raise RuntimeError('Unexpected empty column or query, check test setup')

            finally:
                tbt.done()

            return col

        def check_stats(self, exp_mean, exp_median, places):
            nans = np.isnan(np.sum(col))
            self.assertFalse(nans)
            # ALMA test datasets have large numbers of 0s
            # zeros_count = np.count_nonzero(col==0)
            # self.assertEqual(0, zeros_count)
            if exp_mean is not None:
                dmean = np.mean(col)
                self.assertAlmostEqual(dmean, exp_mean, places=places)
            if exp_median is not None:
                dmedian = np.median(col)
                self.assertAlmostEqual(dmedian, exp_median, places=places)
            if exp_min is not None:
                dmin = col.min()
                self.assertAlmostEqual(dmin, exp_min, places=places)
            if exp_max is not None:
                dmax = col.max()
                self.assertAlmostEqual(dmax, exp_max, places=places)

            verbose = False
            if verbose:
                print(f'Mean, median, min, max: {dmean} {dmedian} {dmin} {dmax}')

        col = get_col(vis, col_name, fields, ddis)
        check_stats(self, exp_mean, exp_median, places)

    def _check_rows(self, vis, col_name, expected_rows, expected_val=None):
        """
        Meant to check rows of a column from an output MS produced by
        uvcontsub2021. Uses unittest asserts to verify conditions.

        :param vis: MS to check
        :param col_name: name of column to check (MODEL_DATA, FIELD, etc.)
        :param expected_rows: number of rows the column must have
        :param expected_vals: for simple cases where the same value is expected in every
                              row, ensure all rows have this value
        """

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

    def _check_task_return(self, res, fields=None):
        """
        Checks consistency of the uvcontsub task return dictionary

        :param fields: results are expected for these fields (and no other fields)
        """

        import pprint
        verbose = False
        if verbose:
            pprint.pprint(res)

        self.assertTrue('description' in res)
        self.assertTrue('goodness_of_fit' in res)
        self.assertTrue('field' in res['goodness_of_fit'])
        gof_field = res['goodness_of_fit']['field']
        if fields:
            for fid in fields:
                self.assertTrue(str(fid) in gof_field, f'Expected field {fid} is not found'
                                'in returned value')

        for fid in gof_field:
            self.assertFalse(fields and str(fid) in fields, f'Field {fid} is found in the '
                             'returned value but it was not expected')
            self.assertTrue('scan' in gof_field[fid])
            scans = gof_field[fid]['scan']
            for sid in scans:
                self.assertTrue('spw' in scans[sid])
                spws = scans[sid]['spw']
                for spw_id in spws:
                    self.assertTrue('polarization' in spws[spw_id])
                    pols = spws[spw_id]['polarization']
                    for pid in pols:
                        self.assertTrue('chi_squared') in pols[pid]
                        stats = pols[pid]['chi_squared']
                        for metric in ['average', 'min', 'max']:
                            self.assertTrue(metric in stats)
                            self.assertEqual(stats[metric].keys(), {'real', 'imag'})
                            self.assertGreaterEqual(stats[metric]['real'], 0)
                            self.assertGreaterEqual(stats[metric]['imag'], 0)

    def _check_input_output_model(self, vis, outputvis, in_col_name):
        """
        Checks (with unittest assert) that INPUT/colname = OUTPUT/DATA + OUTPUT/MODEL_DATA

        :param outputvis: name of input MS with DATA or CORRECTED_DATA
        :param outputvis: name of output with DATA and MODEL_DATA
        :param in_col_name: name of input data column to compare (DATA or CORRECTED_DATA)
        """
        tbt = table()
        try:
            tbt.open(vis)
            data_orig = tbt.getcol(in_col_name)
        finally:
            tbt.close()

        try:
            tbt.open(outputvis)
            data_sub = tbt.getcol('DATA')
            model_sub = tbt.getcol('MODEL_DATA')
        finally:
            tbt.close()

        self.assertTrue(np.allclose(data_orig, data_sub+model_sub, rtol=1e-6),
                        f'Output DATA and MODEL_DATA do not add up to input data column, '
                        f'Output DATA: {data_sub}\n'
                        f'Output MODEL_DATA: {model_sub}\n'
                        f'Input {in_col_name}: {data_orig}\n'
                        )


class uvcontsub2021_test(uvcontsub2021_test_base):
    """
    Main verification test for uvcontsub2021
    """

    @classmethod
    def setUpClass(cls):
        shutil.copytree(datapath_simple, ms_simple)
        shutil.copytree(datapath_alma, ms_alma)
        shutil.copytree(datapath_corr, ms_corr)
        shutil.copytree(datapath_papersky, ms_papersky)
        shutil.copytree(datapath_mixed_pols, ms_mixed_pols)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(ms_simple)
        shutil.rmtree(ms_alma)
        shutil.rmtree(ms_corr)
        shutil.rmtree(ms_papersky)
        shutil.rmtree(ms_mixed_pols)

    def setUp(self):
        # Input MS is always strictly read-only, one copy in setUpClass is enough
        # Default output name for simple tests
        self.output = 'test_uvcs_output.ms'

    def tearDown(self):
        if os.path.exists(self.output):
            shutil.rmtree(self.output)

    def test_makes_output_ms_data(self):
        """
        Check that in a simple command the input MS is taken and an output MS
        is created and has a data column
        """

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output)
        self._check_task_return(res, fields=[0])
        self.assertTrue(os.path.exists(self.output))
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

        # check also no-overwrite of existing MS
        with self.assertRaises(ValueError):
            uvcontsub2021(vis=ms_simple, outputvis=self.output)

    def test_select_field(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='2')
        self._check_task_return(res, fields=[2])
        self._check_rows(self.output, 'FIELD_ID', 120, 2)

    def test_select_spw(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, spw='1')
        self._check_task_return(res, fields=[0, 1, 2])
        self._check_rows(self.output, 'DATA_DESC_ID', 810, 1)

    def test_select_scan(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, scan='2')
        self._check_task_return(res, fields=[1])
        self._check_rows(self.output, 'SCAN_NUMBER', 360, 2)

    def test_select_intent(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, intent='*AMPLI*')
        self._check_task_return(res, fields=[1])  # fields 0,2 excluded by intent selection
        self._check_rows(self.output, 'SCAN_NUMBER', 360, 2)

    def test_select_array(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, array='0')
        self._check_task_return(res, fields=[0, 1, 2])
        self._check_rows(self.output, 'ARRAY_ID', 1080, 0)

    def test_select_observation(self):
        """ Check field selection works"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, observation='0')
        self._check_task_return(res, fields=[0, 1, 2])
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

        # 'datacolumn' test using DATA:
        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, datacolumn='DATA')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

        # 'datacolumn' test using CORRECTED:
        shutil.rmtree(self.output)
        res = uvcontsub2021(vis=ms_corr, outputvis=self.output, datacolumn='CORRECTED')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 2)
        self._check_data_stats(ms_corr, (1419.16761+3.06690944e-06j),
                               (1400.83813+1.36621107e-12j), (357.430084+0j),
                               (4006.86426+0.000186866833j), col_name='CORRECTED_DATA')
        self._check_data_stats(self.output, (1419.16761+3.06690944e-06j),
                               (1400.83813+1.36621107e-12j), (357.430084+0j),
                               (4006.86426+0.000186866833j))

    def test_fitspec_empty(self):
        """Check that fitspec works. When empty, fit all channels in all SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspec='')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitspec_spws(self):
        """Check that fitspec works. When selecting some SPWs, fit all channels
        in those SPWs"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspec='0')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j),
                               fields=0, ddis=0)

    def test_fitspec_spw_one_chan(self):
        """Check fitspec when giving one spw with 1 channel (perfect fit if order 0)"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='0', fitspec='1')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 600)
        self._check_data_stats(self.output, (0.115759703+7.48776972e-08j), 0j,
                               (-0.00129960873+0.000193971457j), (1.41059673+0j))
        self._check_data_stats(self.output, 0j, 0j, 0j, 0j, fields=0, ddis=1)

    def test_fitspec_dict_fitspec_one_chan(self):
        """Check fitspec when giving in a dict one spw with 1 channel (perfect fit if
        order 0)"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='0',
                            fitspec=[['0', '1', 0]])
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 600)
        self._check_data_stats(self.output, (0.115759703+7.48776972e-08j), 0j,
                               (-0.00129960873+0.000193971457j), (1.41059673+0j))

    def test_fitspec_dict_fitspec_one_chan_order1(self):
        """Check fitspec when giving in a dict one spw with 1 channel (perfect fit if
        order 0)"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='0',
                            fitspec=[['0', '1', 1]])
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 600)
        self._check_data_stats(self.output, (0.115759703+7.48776972e-08j), 0j,
                               (-0.00129960873+0.000193971457j), (1.41059673+0j))

    def test_fitspec_channels(self):
        """Check that fitspec works. When selecting some channels in some SPWs,
        fit those channels in those SPWs (like example 2 from task page)"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitspec='0:5~19')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, (-18.5249996-7.05000019j),
                               (-19.7999992-13j), (-68-17.6000004j), (28+26.3999996j),
                               fields=0, ddis=0)


    def test_fitspec_multifield(self):
        """Check that fitspec works. Different spw:chan strings for different fields
        (like example 4 from task page)"""

        res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                            fitspec=[
                                ['0', '0:100~500;600~910;1215~1678;1810~1903', 0],
                                ['1', 'NONE', 0],
                                ['2', '0:100~1903', 0]
                            ])
        self._check_task_return(res, fields=[0, 2])
        self._check_rows(self.output, 'DATA', 1080)
        self._check_data_stats(self.output, (0.0286860488-2.65735951e-06j), 0j,
                               (-0.655080259+0j), (2.09603309+0j))

    def test_fitspec_multifield_blocks(self):
        """Check that fitspec works. Different spw:chan strings for different fields
        but giving list of fields for a same fitspec"""

        # Give some fields grouped
        res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                            fitspec=[
                                ['0, 1', '0:100~500;600~900;1200~1900', 0],
                                ['2', '0:100~1903', 0]
                            ])
        self._check_task_return(res, fields=[0, 1, 2])
        self._check_rows(self.output, 'DATA', 1080)
        expected_vals = [(-0.0125788084-5.98476360e-06j), 0j,
                         (-0.664994299+0j), (2.09603309+0j)]
        self._check_data_stats(self.output, *expected_vals)

        # Giving the fields one at a time should be equivalent
        shutil.rmtree(self.output)
        res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                            fitspec=[
                                ['0', '0:100~500;600~900;1200~1900', 0],
                                ['1', '0:100~500;600~900;1200~1900', 0],
                                ['2', '0:100~1903', 0]
                            ])
        self._check_task_return(res, fields=[0, 1, 2])
        self._check_rows(self.output, 'DATA', 1080)
        self._check_data_stats(self.output, *expected_vals)

    def test_fitspec_multifield_multispw_diff_fitorder(self):
        """Check different fitorder values for different fields and spws, in
        addition to different spw:chan strings"""
        res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                            fitspec=[
                                ['0', '0:100~500;600~910;1215~1678;1810~1903', 1],
                                ['0', '1:0', 2],
                                ['1', 'NONE', 1],
                                ['2', '0:100~1903', 2]
                            ])
        self._check_task_return(res, fields=[0, 2])
        self._check_rows(self.output, 'DATA', 1080)
        self._check_data_stats(self.output, (0.0361723897-2.72414129e-06j), 0j,
                               (-1.22063279+0j), (2.01062560+0j))

    def test_fitspec_multifield_fields_sel(self):
        """Check the use of field selection (multiple) and fitspec (multiple/ per field
        list) together"""
        res = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='1,2',
                            fitspec=[
                                ['1', '0:100~500;600~900', 0],
                                ['2', '0:100~1903', 0]
                            ])
        self._check_task_return(res, fields=[1, 2])
        self._check_rows(self.output, 'DATA', 480)
        self._check_data_stats(self.output, (-0.0165902558-1.55947462e-05j),
                               0j, (-0.702126324+0j), (2.01062560+0j))

    def test_fitspec_multifield_wrong_field(self):
        """Check that wrong fitspec lists produce an exception"""

        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspec=[
                                    ['99', '0:100~500;600~910;1215~1678;1810~1903', 3]
                                ])
            self._check_task_return(res)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspec=[
                                    ['-2', '0:100~500', 2]
                                ])
            self._check_task_return(res)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspec=[
                                    ['4', '0,1', 1]
                                ])
            self._check_task_return(res)

    def test_fitspec_multifield_wrong_format(self):
        """Check that fitspec works. Different fitspec strings for different fields
        (like example 4 from task page)"""

        with self.assertRaises(RuntimeError):
            # Wrong number of elements (not a list of pairs) - in third line
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspec=[
                                    ['1', 'NONE', 0],
                                    ['2', '0:100~500;600~910;1215~1678;1810~1903', 0],
                                    ['3', '4', '0:100~1903', 0]
                                ])

            self._check_task_return(res)
            self._check_rows(self.output, 'DATA', 1080)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            # Wrong field
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspec=[
                                    ['1', 'NONE', 3],
                                    ['bla_fail', '0:600~910;1215~1678', 2],
                                    ['3', '0:100~1903', 1]
                                ])

            self._check_task_return(res)
            self._check_rows(self.output, 'DATA', 1080)

        shutil.rmtree(self.output)
        with self.assertRaises(RuntimeError):
            # Repeated indices
            res = uvcontsub2021(vis=ms_alma, outputvis=self.output,
                                fitspec=[
                                    ['1, 2', 'NONE', 1],
                                    ['2', '0:100~500;600~910;1215~1678;1810~1903', 1],
                                ])
            self._check_task_return(res)
            self._check_rows(self.output, 'DATA', 1080)

    def test_fitspec_separate_fields(self):
        """Check that fitspec works. Different spw:chan strings for different
        fields, and each field to a different output MS (like example 3 from
        task page)"""

        res_f1 = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='1',
                               fitspec='0:100~500;600~910;1215~1678;1810~1903')
        self._check_task_return(res_f1, fields=[1])
        self._check_rows(self.output, 'DATA', 360)
        self._check_data_stats(self.output, (-0.013026415-2.0328553e-05j),
                               0j, (-0.633247495+0j), (2.01062560+0j))

        shutil.rmtree(self.output)
        res_f2 = uvcontsub2021(vis=ms_alma, outputvis=self.output, field='2',
                               fitspec='0:100~1303')
        self._check_task_return(res_f2, fields=[2])
        self._check_rows(self.output, 'DATA', 120)
        self._check_data_stats(self.output, (-0.0140258259+5.83529241e-06j),
                               0j, (-0.588652134+0j), (1.83768487+0j))

    def test_fitspec_spws_sel(self):
        """Check the use of spw selection and fitspec together"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, spw='0', fitspec='0')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitspec_spws_sel_wrong(self):
        """Check that if the spw selection and fitspec are not compatible, an exception
        is produced"""

        with self.assertRaises(RuntimeError):
            res = uvcontsub2021(vis=ms_simple, outputvis=self.output, spw='3', fitspec='0')
            self._check_task_return(res)

    def test_fitmethod_gsl(self):
        """Check that methods work - gsl"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitmethod='gsl')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitmethod_casacore(self):
        """Check that methods work - casacore"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitmethod='casacore')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j))

    def test_fitorder1(self):
        """ Check different fit orders (0, 1, 2) work (like example 1 from task page)"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=1,
                            fitspec='0:2~20')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, (-8.56096448-3.18684196j),
                               (-4.05263185-2.67017555j), (-60.7157898-9.26315689j),
                               (45.24561309+28.1754398j))

    def test_fitspec_dict_fitorder1(self):
        """Check different fit orders (0, 1, 2) work, when given in field/spw
        dict, and produce same results as simple command

        """

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=1,
                            fitspec=[['0', '0:2~20', 1]])
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, (-8.56096448-3.18684196j),
                               (-4.05263185-2.67017555j), (-60.7157898-9.26315689j),
                               (45.24561309+28.1754398j))

    def test_fitorder2(self):
        """ Check different fit orders (0, 1, 2) work"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=2,
                            fitspec='0:2~20')
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, (8.891941398+3.44454119j),
                               (2.90232849+1.83334923j), (-42.0941238-14.0698833j),
                               (92.6959686+31.4091110j))

    def test_fitspec_dict_fitorder2(self):
        """Check different fit orders (0, 1, 2) work, when given in field/spw
        dict, and produce same results as simple command

        """

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, fitorder=2,
                            fitspec=[[ '0', '0:2~20', 2]])
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_data_stats(self.output, (8.891941398+3.44454119j),
                               (2.90232849+1.83334923j), (-42.0941238-14.0698833j),
                               (92.6959686+31.4091110j))
    def test_writemodel(self):
        """ Check the model column is added to the output MS and its values match
        (like example 5 from task page)"""

        res = uvcontsub2021(vis=ms_simple, outputvis=self.output, writemodel=True)
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 340)
        self._check_rows(self.output, 'MODEL_DATA', 340)
        self._check_data_stats(ms_simple, (33.875+15.75j), (15+12j), (2+3j), (98+47j),
                               col_name='DATA')
        self._check_data_stats(self.output, 0j, (-8.25-5.5j), (-42.5-11j), (53.5+33j),
                               col_name='DATA')
        self._check_data_stats(self.output, (33.875+15.75j), (33.875+15.75j), (23.25+17.5j),
                               (44.5+14j),
                               col_name='MODEL_DATA')
        # Besides checks on columns one at a time, check input/DATA == output/DATA+MODEL
        self._check_input_output_model(ms_simple, self.output, in_col_name='DATA')

    def test_writemodel_from_corrected(self):
        """ Check the model column, like test_writemodel, but taking input from
        CORRECTED_DATA"""
        res = uvcontsub2021(vis=ms_papersky, outputvis=self.output, datacolumn='CORRECTED',
                            writemodel=True)
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 6318)
        self._check_rows(self.output, 'MODEL_DATA', 6318)
        self._check_data_stats(ms_papersky, (-0.162619720-0.0525208352j),
                               (-0.0567376502-0.756820410j), (-47.9221878-9.60179615j),
                               (8.73174095-5.35791397j),
                               col_name='CORRECTED_DATA')
        self._check_data_stats(ms_papersky, (-0.162619720-0.0525208352j),
                               (-0.0567376502-0.756820410j), (-47.9221878-9.60179615j),
                               (8.73174095-5.35791397j),
                               col_name='DATA')
        self._check_data_stats(self.output, (-2.19514804e-10+7.80831479e-11j),
                               (-0.00495260954-0.259171128j), (-30.7899017-11.9215565j),
                               (18.17750740-0.306599379j),
                               col_name='DATA')
        self._check_data_stats(self.output, (-0.162619720-0.0525208355j),
                               (-0.0806645453-0.0697010457j), (-17.1322861+2.31976032j),
                               (3.02633667+0.301503152j),
                               col_name='MODEL_DATA')

        self._check_input_output_model(ms_papersky, self.output,
                                       in_col_name='CORRECTED_DATA')

    def test_writemodel_from_corrected_all_flagged(self):
        """ Check the model column, like test_writemodel, but taking input from
        CORRECTED_DATA, with all data points flagged (fit doesn't iterate, model=0) """
        res = uvcontsub2021(vis=ms_corr, outputvis=self.output, datacolumn='CORRECTED',
                            writemodel=True)
        # perhaps TODO: check all chisq == inf, and all count == 1
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 2)
        self._check_rows(self.output, 'MODEL_DATA', 2)
        self._check_data_stats(ms_corr, (1419.16761+3.06690944e-06j),
                               (1400.83813+1.36621107e-12j), (357.430084+0j),
                               (4006.86426+0.000186866833j),
                               col_name='CORRECTED_DATA')
        self._check_data_stats(ms_corr, (58.7170534+8.56816769e-08j), (58.2515029+0j),
                               (14.4238396+0j), (163.592163+7.62939453e-06j),
                               col_name='DATA')
        self._check_data_stats(self.output, (1419.16761+3.06690944e-06j),
                               (1400.83813+1.36621107e-12j), (357.430084+0j),
                               (4006.86426+0.000186866833j),
                               col_name='DATA')
        # This MS is all flagged -> fits produce chisq inf and null fit
        self._check_data_stats(self.output, 0, 0, 0, 0,
                               col_name='MODEL_DATA')

        self._check_input_output_model(ms_corr, self.output, in_col_name='CORRECTED_DATA')

    def test_mixed_pols(self):
        """ Check normal functioning with an MS with mixed polarizations in its (many)
        SPWs"""
        res = uvcontsub2021(vis=ms_mixed_pols, outputvis=self.output)

        # Check all SPWs have been taken.
        def count_entries(d):
            return sum([count_entries(v) if isinstance(v, dict) else 1 for v in d.values()])
        self.assertEqual(count_entries(res), 1331)
        self._check_task_return(res, fields=[0])
        self._check_rows(self.output, 'DATA', 1550)
        self._check_data_stats(self.output, (-0.000782105923+0.000300504051j),
                               0j, (-66.9297485+55.1413803j), (57.7461967+16.6702766j))

    def test_fitspec_spws_diff_fitorder(self):
        """Check that fitorder is applied correctly per field, per SPW, following
        the fitspec list with different fit orders (between 0 and 3) for different spws"""
        # Uses ms_mixed_pols: many SPWs (and with mixed pols), but 1 field

        # A: fit using order 0,1,2,3, without per field/spw fitspec
        res_order_0 = uvcontsub2021(vis=ms_mixed_pols, outputvis=self.output, fitorder=0,
                                    fitspec='3')
        self._check_task_return(res_order_0, fields=[0])
        stats_field_0_spw_3_order_0 = (
            (-0.000772132455+0.000219751841j), (-0.000659341924-0.0187813938j),
            (-1.42448926+0.441552132j), (0.793870270-0.328885317j))
        self._check_data_stats(self.output, *stats_field_0_spw_3_order_0, fields=0, ddis=3)
        shutil.rmtree(self.output)

        res_order_1 = uvcontsub2021(vis=ms_mixed_pols, outputvis=self.output, fitorder=1,
                                    fitspec='2:100~500')
        self._check_task_return(res_order_1, fields=[0])
        stats_field_0_spw_2_order_1 = (
            (-0.00390083172-0.00195605258j), (0.00536000915-0.0672568083j),
            (-11.8622923-2.65349817j), (14.1239767+4.88118219j))
        self._check_data_stats(self.output, *stats_field_0_spw_2_order_1, fields=0, ddis=2)
        shutil.rmtree(self.output)

        res_order_2 = uvcontsub2021(vis=ms_mixed_pols, outputvis=self.output, fitorder=2,
                                    fitspec='63:2~50;55~60')
        self._check_task_return(res_order_2, fields=[0])
        stats_field_0_spw_63_order_2 = (
            (-0.00150726591+0.00118314047j), (0.00288047199+0.00183670223j),
            (-2.88742280-0.969364464j), (2.69186521+1.11569095j))
        self._check_data_stats(self.output, *stats_field_0_spw_63_order_2, fields=0, ddis=63)
        shutil.rmtree(self.output)

        res_order_3 = uvcontsub2021(vis=ms_mixed_pols, outputvis=self.output, fitorder=3,
                                    fitspec='4:50~503;600~850;900~1002')
        self._check_task_return(res_order_3, fields=[0])
        stats_field_0_spw_4_order_3 = (
            (-0.00245836405+0.00548351408j), (0.00843574991+1.61623479j),
            (-10.5505962+6.32199335j), (9.71259212+5.48828125j))
        self._check_data_stats(self.output, *stats_field_0_spw_4_order_3, fields=0, ddis=4)
        shutil.rmtree(self.output)

        # B: cross-check A against per field/spw fitspec with corresponding fit orders
        res_order0 = uvcontsub2021(vis=ms_mixed_pols, outputvis=self.output,
                                   fitspec =[
                                       ['0', '1:0', 0],
                                       ['0', '2:100~500', 1],  # 2 pol, 512 chan
                                       ['0', '3', 0],  # 4 pol, 64 chan
                                       ['0', '4:50~503;600~850;900~1002', 3],
                                       # 2 pol, 1024 chan
                                       ['0', '5:100~495', 0],  # 2 pol, 512 chan
                                       ['0', '63:2~50;55~60', 2]  # 4 pol, 64 chan
                                   ])
        self._check_data_stats(self.output, *stats_field_0_spw_2_order_1, fields=0, ddis=2)
        self._check_data_stats(self.output, *stats_field_0_spw_3_order_0, fields=0, ddis=3)
        self._check_data_stats(self.output, *stats_field_0_spw_4_order_3, fields=0, ddis=4)
        self._check_data_stats(self.output, *stats_field_0_spw_63_order_2, fields=0, ddis=63)


class uvcontsub2021_numerical_sim_test(uvcontsub2021_test_base):
    """
    Tests of numerical behavior based on simulated datasets. To be refined - CAS-13632
    """

    @classmethod
    def setUpClass(cls):
        cls.ms_cont_nonoise_order_0 = 'sim_alma_cont_poly_order_0_nonoise.ms'
        cls.ms_cont_noise_order_0 = 'sim_alma_cont_poly_order_0_noise.ms'
        cls.ms_cont_nonoise_order_1 = 'sim_alma_cont_poly_order_1_nonoise.ms'
        cls.ms_cont_noise_order_1 = 'sim_alma_cont_poly_order_1_noise.ms'
        cls.sim_mss = [cls.ms_cont_nonoise_order_0, cls.ms_cont_noise_order_0,
                       cls.ms_cont_nonoise_order_1, cls.ms_cont_noise_order_1]

        for sim in cls.sim_mss:
            datapath_sim = ctsys.resolve(os.path.join(datadir, sim))
            shutil.copytree(datapath_sim, sim)

    @classmethod
    def tearDownClass(cls):
        for sim in cls.sim_mss:
            shutil.rmtree(sim)

    def setUp(self):
        # Input MS is always strictly read-only, one copy in setUpClass is enough
        # Default output name for simple tests
        self.output = 'test_numerical_uvcs_output.ms'

        # A few parameters shared by different simulated MSs. These should match
        # the parameters (channels, polynomial coefficients, etc.) used in the notebook
        # that produced the simulated datasets.
        # fitspec to exclude a simulated spectral line, as added in the simulation notebook
        self.fitspec = '0:0~59;86~127'

        # For order 0 cont, this constant cont is added to each visibility
        pol_coeffs_order_0 = 0.025
        # For order 1 cont, polynomial coefficients (on x=chan_number)
        pol_coeffs_order_1 = [-0.1, 0.75]
        nchan = 128

        # Constant flux of a component added to the simulated datasets
        source_comp_flux = 0.5+0j

        # Synthetic polynomial continuum from the simulations
        # (the simulations used here add the polynomial to both real and imaginary parts)
        self.exp_cont_order_0 = source_comp_flux + (1+1j) * pol_coeffs_order_0
        x_pol_1 = np.linspace(0, 1, nchan)
        self.exp_cont_order_1 = source_comp_flux + (1+1j) * np.polyval(pol_coeffs_order_1,
                                                                       x_pol_1)

        # Thesholds on fit chi-squres for the simulations without / with noise
        self.chi_square_thresholds_order_0 = {'nonoise': {'real': 1.0e-28,
                                                          'imag': 1.0e-30},
                                              'noise': {'0': {'real': 0.98,
                                                              'imag': 1.07},
                                                        '1': {'real': 1.06,
                                                              'imag': 1.07}
                                                        }
                                              }
        self.chi_square_thresholds_order_1 = {'nonoise': {'real': 1.3e-13,
                                                          'imag': 3e-14},
                                              'noise': {'0': {'real': 0.97,
                                                              'imag': 1.06},
                                                        '1': {'real': 1.05,
                                                              'imag': 1.04}
                                                        }
                                              }

    def tearDown(self):
        if os.path.exists(self.output):
            shutil.rmtree(self.output)

    def _check_diffs(self, vis, outputvis, exp_cont):
        """
        Compare numerical differences between input visibilities and
        (output visibilities + expected_continuum), where the expected
        continuum is known from simulations.

        :param vis: input MS
        :param outputvis: output, cont subtracted MS to compare against (input MS - exp_cont)
        :param exp_cont: expected continuum values (known from simulation parameters)

        :returns: a tuple with the 25th, 50th/median, and 75th percentiles of the differences
        in absolute values of the visibilities across all channels, for all polarizations,
        baselines, and times. The difference (in absolute visibility values) is calculated
        as a percentage relative to the expected continuum value at each channel.
        """
        tbt = table()
        try:
            col_name = 'DATA'
            tbt.open(vis)
            col_in = tbt.getcol(col_name)

            tbt.close()
            tbt.open(outputvis)
            col_sub = tbt.getcol(col_name)

            nans = np.isnan(np.sum(col_sub))
            self.assertFalse(nans)

            if not np.isscalar(exp_cont):
                # broadcast to pol/time
                broadcast_exp_cont = exp_cont.reshape(1, len(exp_cont), 1)
            else:
                broadcast_exp_cont = exp_cont
            diff = (col_in) - (col_sub + broadcast_exp_cont)
            dmax = np.max(diff)
            dmin = np.min(diff)
            dmedian = np.median(diff)
            print(f'Diff median: {dmedian}, min: {dmin}, max: {dmax}')

            def pc_relative_diff(diff, ref):
                return 100.0 * np.absolute(diff / ref)

            # Print these percentiles
            pc_levels = [25, 50, 75]
            # Differences in absolute val of visibilities
            diff_abs = np.absolute(diff)
            amedian = np.median(diff_abs)
            amin = np.min(diff_abs)
            amax = np.max(diff_abs)
            rdiff_abs = pc_relative_diff(diff_abs, np.absolute(broadcast_exp_cont))
            ra25, ramedian, ra75 = np.percentile(rdiff_abs, pc_levels)
            ramin = np.min(rdiff_abs)
            ramax = np.max(rdiff_abs)
            print(f' Diff in absolute values. Median: {amedian}, min: {amin}, max: {amax}.'
                  f'\n   Relative to cont, 25pc: {ra25}, median: {ramedian} %, 75pc: {ra75},'
                  f' min: {ramin} %, max: {ramax} %')

            diff_real = np.absolute(diff.real)
            rmedian = np.median(diff_real)
            rmin = np.min(diff_real)
            rmax = np.max(diff_real)
            rdiff_real = pc_relative_diff(diff.real, broadcast_exp_cont.real)
            rr25, rrmedian, rr75 = np.percentile(rdiff_real, pc_levels)
            rrmin = np.min(rdiff_real)
            rrmax = np.max(rdiff_real)
            print(f' Diff in real part. Median: {rmedian}, min: {rmin}, max: {rmax}'
                  f'\n   Relative to cont, 25pc: {rr25}, median: {rrmedian} %, 75pc: {rr75},'
                  f' min: {rrmin} %, max: {rrmax} %')

            diff_imag = np.absolute(diff.imag)
            imedian = np.median(diff_imag)
            imin = np.min(diff_imag)
            imax = np.max(diff_imag)
            rdiff_imag = pc_relative_diff(diff.imag, broadcast_exp_cont.imag)
            ri25, rimedian, ri75 = np.percentile(rdiff_imag, pc_levels)
            rimin = np.min(rdiff_imag)
            rimax = np.max(rdiff_imag)
            print(f' Diff in imag part. Median: {imedian}, min: {imin}, max: {imax}'
                  f'\n   Relative to cont, 25pc: {ri25}, median: {rimedian} %, 75pc: {ri75},'
                  f' min: {rimin} %, max: {rimax} %')

            return ra25, ramedian, ra75
        finally:
            tbt.done()

    def _assert_chi_sq_values_nonoise(self, res, chi_sq_thresholds):
        """
        Compare and enforce assert on chi_squared thresholds, for MSs
        without noise. In these MSs, avg/min/max chi_sq values are the same,
        and both polarizations are also the same.

        :param res: result dictionary from uvcontsub2021
        :param chi_sq_thresholds: thresholds to compare (resuls should be same or
                                  better==smaller chi_sq)

        """
        chi_sq_imag = chi_sq_thresholds['nonoise']['imag']
        chi_sq_real = chi_sq_thresholds['nonoise']['real']
        for agg in ['average', 'max', 'min']:
            for pol in ['0', '1']:
                self.assertLessEqual(res['goodness_of_fit']['field']['0']['scan']['1']['spw']
                                     ['0']['polarization'][pol]['chi_squared'][agg]['imag'],
                                     chi_sq_imag)
                self.assertLessEqual(res['goodness_of_fit']['field']['0']['scan']['1']['spw']
                                     ['0']['polarization'][pol]['chi_squared'][agg]['real'],
                                     chi_sq_real)

    def _assert_chi_sq_values_with_noise(self, res, chi_sq_thresholds):
        """
        Compare and enforce assert on chi_squared thresholds, for MSs
        with noise. For simplicity uses only average chi_sq per polarization.

        :param res: result dictionary from uvcontsub2021
        :param chi_sq_thresholds: thresholds on avg chi_sq to compare (resuls should
                                  be same or better==smaller chi_sq)
        """
        for pol in ['0', '1']:
            chi_sq_real = self.chi_square_thresholds_order_0['noise'][pol]['real']
            chi_sq_imag = self.chi_square_thresholds_order_0['noise'][pol]['imag']

            self.assertLessEqual(res['goodness_of_fit']['field']['0']['scan']['1']['spw']
                                 ['0']['polarization'][pol]['chi_squared']['average']
                                 ['real'],
                                 chi_sq_real)
            self.assertLessEqual(res['goodness_of_fit']['field']['0']['scan']['1']['spw']
                                 ['0']['polarization'][pol]['chi_squared']['average']
                                 ['imag'],
                                 chi_sq_imag)

    def test_sim_specline_nonoise_pol_0(self):
        """ Check fitting of continuum as polynomial order 0"""
        res = uvcontsub2021(vis=self.ms_cont_nonoise_order_0, outputvis=self.output,
                            fitorder=0, fitspec=self.fitspec)
        self._check_task_return(res, fields=[0])

        # Expected cont form simulations. Added to visibilities as a polynomial on channels
        exp_cont = self.exp_cont_order_0
        print(f'Checking numerical differences for MS {self.ms_cont_nonoise_order_0}')
        diff25, diff50, diff75 = self._check_diffs(vis=self.ms_cont_nonoise_order_0,
                                                   outputvis=self.output,
                                                   exp_cont=exp_cont)
        # These percentiles are ~4.5e-6
        self.assertLessEqual(diff25, 1e-5)
        self.assertLessEqual(diff50, 1e-5)
        self.assertLessEqual(diff50, 1e-5)

        # Check also residuals from the viewpoint of the returned dict (chi_squared)
        self._assert_chi_sq_values_nonoise(res, self.chi_square_thresholds_order_0)

    def test_sim_specline_noise_pol_0(self):
        """ Check fitting of continuum as polynomial order 0"""
        res = uvcontsub2021(vis=self.ms_cont_noise_order_0, outputvis=self.output,
                            fitorder=0, fitspec=self.fitspec)
        self._check_task_return(res, fields=[0])

        exp_cont = self.exp_cont_order_0
        print(f'Checking numerical differences for MS {self.ms_cont_noise_order_0}')
        diff25, diff50, diff75 = self._check_diffs(vis=self.ms_cont_noise_order_0,
                                                   outputvis=self.output,
                                                   exp_cont=exp_cont)
        # Values from sim dataset + small (<0.05%) tolerance
        self.assertLessEqual(diff25, 1.45)
        self.assertLessEqual(diff50, 2.25)
        self.assertLessEqual(diff50, 3.2)

        self._assert_chi_sq_values_with_noise(res, self.chi_square_thresholds_order_0)

    def test_sim_specline_nonoise_pol_1(self):
        """ Check fitting of continuum as polynomial order 1"""
        res = uvcontsub2021(vis=self.ms_cont_nonoise_order_1, outputvis=self.output,
                            fitorder=1, fitspec=self.fitspec)
        self._check_task_return(res, fields=[0])

        # Values added to visibilities as a polynomial on channels
        exp_cont = self.exp_cont_order_1
        print(f'Checking numerical differences for MS {self.ms_cont_nonoise_order_1}')
        diff25, diff50, diff75 = self._check_diffs(vis=self.ms_cont_nonoise_order_1,
                                                   outputvis=self.output,
                                                   exp_cont=exp_cont)
        # These percentiles are ~[2, 6]e-6
        self.assertLessEqual(diff25, 1e-5)
        self.assertLessEqual(diff50, 1e-5)
        self.assertLessEqual(diff50, 1e-5)

        self._assert_chi_sq_values_nonoise(res, self.chi_square_thresholds_order_1)

    def test_sim_specline_noise_pol_1(self):
        """ Check fitting of continuum as polynomial order 1. Gaussian noise included"""
        res = uvcontsub2021(vis=self.ms_cont_noise_order_1, outputvis=self.output,
                            fitorder=1, fitspec=self.fitspec)
        self._check_task_return(res, fields=[0])

        exp_cont = self.exp_cont_order_1
        print(f'Checking numerical differences for MS {self.ms_cont_noise_order_1}')
        diff25, diff50, diff75 = self._check_diffs(vis=self.ms_cont_nonoise_order_1,
                                                   outputvis=self.output,
                                                   exp_cont=exp_cont)
        # Values from sim dataset + small (<0.05%) tolerance
        self.assertLessEqual(diff25, 5.47)
        self.assertLessEqual(diff50, 8.47)
        self.assertLessEqual(diff75, 11.97)

        self._assert_chi_sq_values_with_noise(res, self.chi_square_thresholds_order_1)


def suite():
    return [uvcontsub2021_test,
            uvcontsub2021_numerical_verification_test]


if __name__ == '__main__':
    unittest.main()

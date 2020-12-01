import os
import numpy as np
import shutil
import sys
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import sdatmcor
    import casatasks.private.task_sdatmcor as sdatmcor_impl
    import casatasks.private.sdutil as sdutil
    # default isn't used in casatasks

    def default(atask):
        pass
    # for testhelper import
    sys.path.append(
        os.path.dirname(
            os.path.abspath(
                os.path.dirname(__file__))))
    from casatools import ctsys
    from casatools import calibrater
    from casatools import ms as mstool
    from casatasks import gencal, applycal

    ctsys_resolve = ctsys.resolve

else:
    from tasks import sdatmcor
    import task_sdatmcor as sdatmcor_impl

    from tasks import gencal, applycal
    from __main__ import default
    from taskinit import cbtool as calibrater
    from taskinit import mstool
    import sdutil

    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0], 'data')
        subdir_hints = ['', 'casa-data-req']
        for subdir in subdir_hints:
            path = os.path.join(dataPath, subdir, apath)
            if os.path.exists(path):
                return path


def smart_remove(name):
    if os.path.exists(name):
        if os.path.isdir(name):
            shutil.rmtree(name)
        elif os.path.isfile(name):
            os.remove(name)
        else:
            # could be a symlink
            os.remove(name)


def std_with_clip(arr, num_iter=1, factor=3):
    assert num_iter > 0
    std = np.abs(arr).max()
    for i in range(num_iter):
        std = arr[np.abs(arr) < std * factor].std()
    return std


def read_table(name, spw, cols=['STATE_ID', 'DATA']):
    ms = mstool()
    idx = ms.msseltoindex(name, spw=[int(spw)])
    ddid = idx['dd']
    with sdutil.tbmanager(name) as tb:
        tsel = tb.query('DATA_DESC_ID IN [{}]'.format(','.join([str(i) for i in ddid])))
        try:
            result_dict = dict((k, tsel.getcol(k)) for k in cols if k in tb.colnames())
        finally:
            tsel.close()
    return result_dict


def apply_gainfactor(name, spw, factor):
    ms = mstool()
    idx = ms.msseltoindex(name, spw=[int(spw)])
    ddid = idx['dd'][0]
    with sdutil.tbmanager(name, nomodify=False) as tb:
        colnames = tb.colnames()
        tsel = tb.query('DATA_DESC_ID=={}'.format(ddid))
        try:
            for colname in ['DATA', 'FLOAT_DATA', 'CORRECTED_DATA']:
                if colname in colnames:
                    data = tsel.getcol(colname)
                    tsel.putcol(colname, data * factor)
        finally:
            tsel.close()


class test_sdatmcor(unittest.TestCase):
    datapath = ctsys_resolve('visibilities/almasd')
    infile = 'X320b_sel2.ms'
    outfile = infile + '.atmcor'
    caltable = infile + '.k2jycal'

    local_unit_test = False

    def setUp(self):
        default(sdatmcor)

        # default Args
        self.args = {
            'infile': self.infile,
            'datacolumn': 'data',
            'outfile': self.outfile,
        }

        smart_remove(self.infile)
        smart_remove(self.outfile)
        smart_remove(self.caltable)
        shutil.copytree(os.path.join(self.datapath, self.infile), self.infile)

    def tearDown(self):
        print("tearDown::deleting MSs.")
        smart_remove(self.infile)
        smart_remove(self.outfile)
        smart_remove(self.caltable)

    def _check_result_spw(self, spw, is_selected, is_processed, on_source_only):

        contents_after = read_table(self.outfile, spw, ['STATE_ID', 'DATA', 'CORRECTED_DATA'])
        stateids_after = contents_after['STATE_ID']
        if 'CORRECTED_DATA' in contents_after:
            data_after = contents_after['CORRECTED_DATA'].real
        else:
            data_after = contents_after['DATA'].real

        # if is_selected is False, selected data should be empty
        if not is_selected:
            self.assertEqual(len(stateids_after), 0)
            return

        contents_before = read_table(self.infile, spw, ['STATE_ID', 'DATA', 'CORRECTED_DATA'])
        stateids_before = contents_before['STATE_ID']
        if 'CORRECTED_DATA' in contents_before:
            data_before = contents_before['CORRECTED_DATA'].real
        else:
            data_before = contents_before['DATA'].real

        # check if spw exist in the data
        self.assertGreater(len(stateids_after), 0)

        if is_processed:
            if spw == 23:
                # examine averaged spectrum
                #   - exclude edge channels (4 channels)
                #   - take average along pol and time axes
                #   - subtract zero order baseline
                #   - compute std with clipping (3 iterations => 2 clipping)
                #   - set threshold to 7.5 times std
                #   - check if any data exceeds threshold for data before correction
                #   - check if no data exceeds threshold for data after correction
                edge = 4
                # data before correction
                mask0 = np.logical_or(stateids_before == 14, stateids_before == 84)
                data_on = data_before[:, edge:-edge, mask0]
                data_mean = data_on.mean(axis=(0, 2))
                data_sub = data_mean - np.median(data_mean)
                threshold = 7.5 * std_with_clip(data_sub, num_iter=3)
                self.assertFalse(np.all(np.abs(data_sub) < threshold))
                # data after correction
                mask1 = np.logical_or(stateids_after == 14, stateids_after == 84)
                data_on = data_after[:, edge:-edge, mask1]
                data_mean = data_on.mean(axis=(0, 2))
                data_sub = data_mean - np.median(data_mean)
                threshold = 7.5 * std_with_clip(data_sub, num_iter=3)
                self.assertTrue(np.all(np.abs(data_sub) < threshold))
            elif spw == 19:
                # examine averaged spectrum
                #   - take average along pol and time axes
                #   - take difference of data before and after correction
                #   - check if diff.mean() > 0.95 * data_before_correction.mean()
                #     which means that most of the continuum-like component is
                #     originated by residual of atmospheric emission
                #   - check if std of difference is less than 30% of
                #     (diff.max() - diff.min())
                mask0 = np.logical_or(stateids_before == 14, stateids_before == 84)
                data_on_before = data_before[:, :, mask0]
                data_mean_before = data_on_before.mean(axis=(0, 2))
                mask1 = np.logical_or(stateids_after == 14, stateids_after == 84)
                data_on_after = data_after[:, :, mask1]
                data_mean_after = data_on_after.mean(axis=(0, 2))
                diff = data_mean_before - data_mean_after
                diff_mean = diff.mean()
                diff_std = diff.std()
                diff_std_norm = diff_std / (diff.max() - diff.min())
                self.assertGreater(diff_mean, data_mean_after.mean() * 0.95)
                self.assertLess(diff_std_norm, 0.3)

            if not on_source_only:
                # OFF_SOURCE data should not be touched
                mask0 = np.logical_and(stateids_before != 14, stateids_before != 84)
                data_off_before = data_before[:, :, mask0]
                mask1 = np.logical_and(stateids_after != 14, stateids_after != 84)
                data_off_after = data_after[:, :, mask1]
                self.assertTrue(np.all(data_off_before == data_off_after))
        else:
            self.assertTrue(np.all(data_after == data_before))

    def check_result(self, spwprocess, on_source_only=False):
        """Check Result

        Args:
            spwprocess (dict): key is spw id, value is whether or not the spw is processed
        """
        # outfile should exist
        self.assertTrue(os.path.exists(self.outfile))

        for spw in [19, 23]:
            is_selected = spw in spwprocess
            is_processed = spwprocess.get(spw, False)
            self._check_result_spw(spw, is_selected, is_processed, on_source_only)

    def test_sdatmcor_normal(self):
        '''test normal usage of sdatmcor'''
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data')
        self.check_result({19: True, 23: True})

    def test_sdatmcor_explicit_atmtype(self):
        '''test specifying atmtype explicitly'''
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', atmtype=2)
        self.check_result({19: True, 23: True})

    def test_sdatmcor_overwrite(self):
        '''test overwriting existing outfile'''
        os.mkdir(self.outfile)
        self.assertTrue(os.path.exists(self.outfile))
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', overwrite=True)
        self.check_result({19: True, 23: True})

    def test_sdatmcor_no_overwrite(self):
        '''test to avoid overwriting existing outfile'''
        os.mkdir(self.outfile)
        self.assertTrue(os.path.exists(self.outfile))
        with self.assertRaises(Exception):
            sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', overwrite=False)

    def test_sdatmcor_wrong_datacolumn(self):
        '''test wrong datacolumn'''
        wrong_colnames = ['corrected', 'float_data']
        for colname in wrong_colnames:
            with self.assertRaises(Exception):
                sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn=colname)

    def test_sdatmcor_corrected(self):
        '''test if CORRECTED_DATA column is handled properly'''
        # add CORRECTED_DATA column
        cb = calibrater()
        cb.open(self.infile, addcorr=True, addmodel=False)
        cb.close()
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='corrected')
        self.check_result({19: True, 23: True})

    def test_sdatmcor_spw_select_23(self):
        '''test data selection: select spw 23'''
        sdatmcor(infile=self.infile, outputspw='23', outfile=self.outfile, datacolumn='data')
        self.check_result({23: True})

    def test_sdatmcor_spw_select_all(self):
        '''test data selection: select 19 and 23 explicitly'''
        sdatmcor(infile=self.infile, outputspw='19,23', outfile=self.outfile, datacolumn='data')
        self.check_result({19: True, 23: True})

    def test_sdatmcor_spw_process_19(self):
        '''test data selection: process only spw 19'''
        sdatmcor(infile=self.infile, spw='19', outfile=self.outfile, datacolumn='data')
        self.check_result({19: True, 23: False})

    def test_sdatmcor_spw_process_all(self):
        '''test data selection: declare to process 19 and 23 explicitly'''
        sdatmcor(infile=self.infile, spw='19,23', outfile=self.outfile, datacolumn='data')
        self.check_result({19: True, 23: True})

    def test_sdatmcor_spw_process_23_select_23(self):
        '''test data selection: select and process spw 23'''
        sdatmcor(infile=self.infile, spw='23', outputspw='23', outfile=self.outfile, datacolumn='data')
        self.check_result({23: True})

    def test_sdatmcor_spw_process_all_select_19(self):
        '''test data selection: process spw 19 and 23 but output only spw 19'''
        sdatmcor(infile=self.infile, spw='19,23', outputspw='19', outfile=self.outfile, datacolumn='data')
        self.check_result({19: True})

    def test_sdatmcor_spw_process_23_select_all(self):
        '''test data selection: process only spw 23 but output both 19 and 23'''
        sdatmcor(infile=self.infile, spw='23', outputspw='19,23', outfile=self.outfile, datacolumn='data')
        self.check_result({19: False, 23: True})

    def test_sdatmcor_spw_process_99_select_all(self):
        '''test data selection: specify invalid spw to process'''
        with self.assertRaises(Exception):
            sdatmcor(infile=self.infile, spw='19,99', outputspw='19,23', outfile=self.outfile, datacolumn='data')

    def test_sdatmcor_intent_selection(self):
        '''test intent selection: test if selection of ON_SOURCE data (i.e. excluding OFF_SOURCE data) still works'''
        sdatmcor(infile=self.infile, outfile=self.outfile, intent='OBSERVE_TARGET#ON_SOURCE*', datacolumn='data')
        self.check_result({19: True, 23: True}, on_source_only=True)

    def test_sdatmcor_spw_process_less_than_20_select_all(self):
        '''test data selection: specify invalid spw to process'''
        sdatmcor(infile=self.infile, spw='<20', outputspw='', outfile=self.outfile, datacolumn='data')
        self.check_result({19: True, 23: False})

    def test_sdatmcor_scan_selection(self):
        """test data selection: select only one scan (scan 5)"""
        # just to confirm the task completes without error
        try:
            sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', scan='5')
        except Exception:
            self.fail('sdatmcor should not raise any Exception')

    def test_sdatmcor_antenna_selection(self):
        """Test antenna selection"""
        sdatmcor(infile=self.infile, antenna='PM02', outfile=self.outfile)
        self.check_result({19: True, 23: True})

    def test_sdatmcor_msselect(self):
        """Test msselect"""
        sdatmcor(infile=self.infile, msselect='ANTENNA1 == 1', outfile=self.outfile)
        self.check_result({19: True, 23: True})

        with self.assertRaises(Exception):
            sdatmcor(infile=self.infile, msselect='ANTENNA1 == 2', outfile=self.outfile)

    def test_sdatmcor_gainfactor_float(self):
        """test gainfactor: float input"""
        gainfactor = 10.0
        apply_gainfactor(self.infile, 19, gainfactor)
        apply_gainfactor(self.infile, 23, gainfactor)
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', gainfactor=gainfactor)
        self.check_result({19: True, 23: True})

    def test_sdatmcor_gainfactor_dict(self):
        """test gainfactor: dict input"""
        gainfactor = {'19': 10.0, '23': 45.0}
        apply_gainfactor(self.infile, 19, gainfactor['19'])
        apply_gainfactor(self.infile, 23, gainfactor['23'])
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', gainfactor=gainfactor)
        self.check_result({19: True, 23: True})

    def test_sdatmcor_gainfactor_caltable(self):
        """test gainfactor: caltable input"""
        gainfactor = {'19': 10.0, '23': 45.0}
        for k, v in gainfactor.items():
            p = 1 / np.sqrt(v)
            gencal(vis=self.infile, caltable=self.caltable, caltype='amp', spw=k, parameter=[p])
        applycal(vis=self.infile, gaintable=self.caltable, flagbackup=False)
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='corrected', gainfactor=gainfactor)
        self.check_result({19: True, 23: True})

    def test_parse_gainfactor_exception(self):
        """test exception raised in parse_gainfactor"""
        with self.assertRaises(RuntimeError):
            sdatmcor_impl.parse_gainfactor(self.infile)

    def test_parse_spw(self):
        """test utility functio, parse_spw"""
        test_cases = [
            ('', [19, 23]),
            ('*', [19, 23]),
            ('19,23', [19, 23]),
            ('19', [19]),
            ('<20', [19]),
            ('17,19,21', [19]),
            ('>20', [23]),
            ('0', []),
            ('1000', RuntimeError),
            ([19], TypeError),
        ]
        for spw, expected in test_cases:
            if isinstance(expected, list):
                actual = sdatmcor_impl.parse_spw(self.infile, spw)
                self.assertEqual(
                    actual,
                    expected,
                    msg='failed: input {} result {} expected {}'.format(
                        spw, actual, expected
                    )
                )
            else:
                # error cases
                with self.assertRaises(expected):
                    actual = sdatmcor_impl.parse_spw(self.infile, spw)

    def test_tweak_antenna_selection(self):
        """Test tweak of antenna selection"""
        # common test cases
        test_cases0 = [
            ('', ''),
            ('PM02', 'PM02&&&'),
            ('PM02&', 'PM02&&&'),
            ('PM02&&', 'PM02&&'),
            ('PM02&&&', 'PM02&&&'),
            ('0', '0&&&'),
            ('0&', '0&&&'),
            ('0&&', '0&&'),
            ('0&&&', '0&&&'),
            ('0&1', '0&&&;1&&&'),
            ('0&&1', '0&&&;1&&&'),
            ('0;1', '0&&&;1&&&'),
        ]
        # specific to convert_antenna_spec_autocorr
        test_cases1 = [
            ('0&;1&&', '0&&&;1&&'),
            ('0&&&;1', '0&&&;1&&&')
        ]
        for antenna, expected in test_cases0 + test_cases1:
            actual = sdutil.convert_antenna_spec_autocorr(antenna)
            self.assertEqual(actual, expected)
        # specific to get_antenna_selection_include_autocorr
        test_cases2 = [
            ('0&;1&&', '0&;1&&'),
            ('0&&&;1', '0&&&;1'),
        ]
        for antenna, expected in test_cases0 + test_cases2:
            actual = sdutil.get_antenna_selection_include_autocorr(self.infile, antenna)
            self.assertEqual(actual, expected)

    def test_get_default_antenna(self):
        """test get_default_antenna and relevant function"""
        duration_ref = {
            'PM01': 1024.320,
            'PM02': 1019.808,
            'PM03': 1027.728,
            'PM04': 1046.112,
        }
        count_ref = {'PM01': 474, 'PM02': 475, 'PM03': 474, 'PM04': 474}
        counts, durations = sdatmcor_impl.inspect_flag_cmd(self.infile)
        for k, v in count_ref.items():
            self.assertEqual(v, counts[k])
        for k, v in duration_ref.items():
            self.assertAlmostEqual(v, durations[k], places=4)

        # default antenna should be PM02 (ID 1)
        default_antenna = sdatmcor_impl.get_default_antenna(self.infile)
        self.assertEqual(default_antenna, 1)

    def test_default_antenna_with_selection(self):
        """Test default antenna determination with data selection excluding the best one"""
        # check if the "best" antenna is not in the MAIN table
        with sdutil.tbmanager(os.path.join(self.infile, 'ANTENNA'), nomodify=False) as tb:
            # replace antenna name to simulate the data selection excluding PM02
            # PM01 <-> PM02
            tb.putcell('NAME', 0, 'PM02')
            tb.putcell('NAME', 1, 'PM01')
        sdatmcor(infile=self.infile, outfile=self.outfile, spw='19', datacolumn='data')
        self.check_result({19: True, 23: False})

    def test_default_antenna_with_no_flag_commands(self):
        """Test default antenna determination for empty FLAG_CMD table"""
        # check if the task can handle empty FLAG_CMD table
        with sdutil.tbmanager(os.path.join(self.infile, 'FLAG_CMD'), nomodify=False) as tb:
            tb.removerows(np.fromiter(range(tb.nrows()), dtype=int))
            self.assertEqual(tb.nrows(), 0)
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data')
        self.check_result({19: True, 23: True})


def suite():
    return [test_sdatmcor]


if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

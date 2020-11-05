import glob
import os
import numpy as np
import shutil
import sys
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import sdatmcor
    # default isn't used in casatasks

    def default(atask):
        pass
    # for testhelper import
    sys.path.append(
        os.path.dirname(
            os.path.abspath(
                os.path.dirname(__file__))))
    from casatasks.private.sdutil import tbmanager
    from casatools import ctsys
    from casatools import calibrater
    from casatools import ms as mstool
    from casatasks import gencal, applycal

    datapath = ctsys.resolve('')

    ctsys_resolve = ctsys.resolve

else:
    from tasks import sdatmcor
    from tasks import gencal, applycal
    from __main__ import default
    from sdutil import tbmanager
    from taskinit import cbtool as calibrater
    from taskinit import mstool

    # Define the root for the data files
    datapath = os.environ.get('CASAPATH').split()[0] + ''

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
    with tbmanager(name) as tb:
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
    with tbmanager(name, nomodify=False) as tb:
        colnames = tb.colnames()
        tsel = tb.query('DATA_DESC_ID=={}'.format(ddid))
        try:
            for colname in ['DATA', 'FLOAT_DATA', 'CORRECTED_DATA']:
                if colname in colnames:
                    data = tsel.getcol(colname)
                    tsel.putcol(colname, data * factor)
        finally:
            tsel.close()


#
# Test-MS
#

# template MS (to copy from)
defInputMs = "uid___A002_Xe770d7_X320b.ms"

# testing MS (set up basic conditions)
defWorkMsBasic = "uid___A002_Xe770d7_X320b-t.ms"

# output MS #


##############
# Test Entry
##############


class test_sdatmcor(unittest.TestCase):
    datapath = ctsys_resolve('visibilities/almasd')
    infile = 'X320b_sel2.ms'
    outfile = infile + '.atmcor'
    caltable = infile + '.k2jycal'

    local_unit_test = False

    def __get_temporary_files(self):
        return glob.glob('_AtmCor_Temp*,ms')

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
        for tmp in self.__get_temporary_files():
            smart_remove(tmp)

# private function #

    def _copy_remote_file(self, infile, outfile):
        print("Copying temp MS from Remote to Local.")

    def _if_exist(self, msname):
        _filePath = os.path.join("./", msname)
        if os.path.exists(_filePath):
            return True
        else:
            return False

#
# Class Method
#
    @classmethod
    def setUpClass(cls):
        print("setUpClass::deleting existing work-MS.")

    @classmethod
    def tearDownClass(cls):
        print("tearDownClass::deleting work-MS.")

##############
# Run Task
##############
    def _run_task(self, auxArgs=None):
        print("_run_task::starts")

        if auxArgs is not None:
            for k in auxArgs:
                self.args[k] = auxArgs[k]

        if not self.local_unit_test:   # FORCE to change
            self.args = {
                'debug': 'skipTaskExec',
                'infile' : infile,
                'outfile': outfile,
                'overwrite':  True }

        # Execution.
        #  if success, returns True
        #  if any error, returns False.
        try:
            return sdatmcor(**self.args)
        except Exception:
            return False

#################
# Check Result
#################


######################
# check time
######################


######################
# check Output
######################


##################################
# Read Data from Specified MS
##################################
    # MAIN #
    def _get_main(self, msname):
        # get MAIN table data
        with tbmanager(msname) as tb:
            # Key data
            self.tm = tb.getcol('TIME')
            self.a1 = tb.getcol('ANTENNA1')
            self.a2 = tb.getcol('ANTENNA2')
            self.sc = tb.getcol('SCAN_NUMBER')
            self.fd = tb.getcol('FIELD_ID')

    # DATA (spectra) #
    def _get_spectra(self, msname, row):
        with tbmanager(msname) as tb:
            # Spectra Data
            self.data = tb.getcell('FLOAT_DATA', row)
            self.wgt = tb.getcell('WEIGHT', row)
            self.sgm = tb.getcell('SIGMA', row)

        return self.data

#####################################
# Generate Data on FLOAT_DATA column
#####################################
    def _generate_data(self, msName, stateOption=False):
        pass

##############
# MISC
##############


############################
# TEST FIXTURE
############################


    # def test_param1(self):
    #     '''sdatmcor::1:: file1  '''

    #     infile = 'uid___A002_Xe770d7_X320b-t.ms'
    #     outfile = 'sdatmcor-out1.ms'
    #     prm = {'debug': 'skipTaskExec',
    #            'infile' : infile,
    #            'outfile': outfile,
    #            'overwrite':  True }
    #     # Run Task and check
    #     self.assertTrue(self._run_task(prm))

    def _check_result_spw(self, spw, is_selected, is_processed, on_source_only, factor):

        contents_after = read_table(self.outfile, spw, ['STATE_ID', 'DATA', 'CORRECTED_DATA'])
        stateids_after = contents_after['STATE_ID']
        if 'CORRECTED_DATA' in contents_after:
            data_after = contents_after['CORRECTED_DATA']
        else:
            data_after = contents_after['DATA'].real

        # if is_selected is False, selected data should be empty
        if not is_selected:
            self.assertEqual(len(stateids_after), 0)
            return

        contents_before = read_table(self.infile, spw, ['STATE_ID', 'DATA', 'CORRECTED_DATA'])
        stateids_before = contents_before['STATE_ID']
        if 'CORRECTED_DATA' in contents_before:
            data_before = contents_before['CORRECTED_DATA']
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
                #   - check if mean of difference is around 0.34
                #   - check if std of difference is small enough
                mask0 = np.logical_or(stateids_before == 14, stateids_before == 84)
                data_on_before = data_before[:, :, mask0]
                data_mean_before = data_on_before.mean(axis=(0, 2))
                mask1 = np.logical_or(stateids_after == 14, stateids_after == 84)
                data_on_after = data_after[:, :, mask1]
                data_mean_after = data_on_after.mean(axis=(0, 2))
                diff = data_mean_before - data_mean_after
                diff_mean = diff.mean()
                diff_std = diff.std()
                self.assertAlmostEqual(diff_mean, 0.307 * factor, places=2)
                self.assertLess(diff_std, 0.0003 * factor)

            if not on_source_only:
                # OFF_SOURCE data should not be touched
                mask0 = np.logical_and(stateids_before != 14, stateids_before != 84)
                data_off_before = data_before[:, :, mask0]
                mask1 = np.logical_and(stateids_after != 14, stateids_after != 84)
                data_off_after = data_after[:, :, mask1]
                self.assertTrue(np.all(data_off_before == data_off_after))
        else:
            self.assertTrue(np.all(data_after == data_before))

    def check_result(self, spwprocess, on_source_only=False, factor={}):
        """Check Result

        Args:
            spwprocess (dict): key is spw id, value is whether or not the spw is processed
        """
        # outfile should exist
        self.assertTrue(os.path.exists(self.outfile))

        # no temporary files exist
        self.assertEqual(len(self.__get_temporary_files()), 0)

        for spw in [19, 23]:
            is_selected = spw in spwprocess
            is_processed = spwprocess.get(spw, False)
            gainfactor = factor.get(str(spw), 1.0)
            self._check_result_spw(spw, is_selected, is_processed, on_source_only, gainfactor)

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

    def test_sdatmcor_intent_selection(self):
        '''test intent selection: test if selection of ON_SOURCE data (i.e. excluding OFF_SOURCE data) still works'''
        sdatmcor(infile=self.infile, outfile=self.outfile, intent='OBSERVE_TARGET#ON_SOURCE*', datacolumn='data')
        self.check_result({19: True, 23: True}, on_source_only=True)

    def test_sdatmcor_gainfactor_float(self):
        """test gainfactor: float input"""
        gainfactor = 10.0
        apply_gainfactor(self.infile, 19, gainfactor)
        apply_gainfactor(self.infile, 23, gainfactor)
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', gainfactor=gainfactor)
        self.check_result(
            {19: True, 23: True},
            factor={'19': gainfactor, '23': gainfactor}
        )

    def test_sdatmcor_gainfactor_dict(self):
        """test gainfactor: dict input"""
        gainfactor = {'19': 10.0, '23': 45.0}
        apply_gainfactor(self.infile, 19, gainfactor['19'])
        apply_gainfactor(self.infile, 23, gainfactor['23'])
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', gainfactor=gainfactor)
        self.check_result(
            {19: True, 23: True},
            factor=gainfactor
        )

    def test_sdatmcor_gainfactor_caltable(self):
        """test gainfactor: caltable input"""
        gainfactor = {'19': 10.0, '23': 45.0}
        for k, v in gainfactor.items():
            p = 1 / np.sqrt(v)
            gencal(vis=self.infile, caltable=self.caltable, caltype='amp', spw=k, parameter=[p])
        applycal(vis=self.infile, gaintable=self.caltable, flagbackup=False)
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='corrected', gainfactor=gainfactor)
        self.check_result(
            {19: True, 23: True},
            factor=gainfactor
        )


def suite():
    return [test_sdatmcor]


if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
"""
History

3-AUG-2020   Revised for Unit-Test until the Task body is completed.
             - added 'skipTaskExec'  on debug argument.

17-SEP-2020  Totally simplified, by minimum functions.
             - This calls the Task with skipTaskExec.
               simply called and returns soon.
"""

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

    datapath = ctsys.resolve('')

    ctsys_resolve = ctsys.resolve

else:
    from tasks import sdatmcor
    from __main__ import default
    from sdutil import tbmanager

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
    infile = 'X320b_sel.ms'
    outfile = infile + '.atmcor'

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
        shutil.copytree(os.path.join(self.datapath, self.infile), self.infile)

    def tearDown(self):
        print("tearDown::deleting MSs.")
        smart_remove(self.infile)
        smart_remove(self.outfile)
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

    def check_result(self):
        # outfile should exist
        self.assertTrue(os.path.exists(self.outfile))

        # no temporary files exist
        self.assertEqual(len(self.__get_temporary_files()), 0)

        # examine averaged spectrum
        #   - exclude edge channels (4 channels)
        #   - take average along pol and time axes
        #   - subtract zero order baseline
        #   - compute std with clipping (3 iterations => 2 clipping)
        #   - set threshold to 7.5 times std
        #   - check if no data exceeds threshold
        with tbmanager(self.outfile) as tb:
            stateids = tb.getcol('STATE_ID')
            data = tb.getcol('DATA').real

        edge = 4
        data_on = data[:, edge:-edge, stateids == 14]
        data_mean = data_on.mean(axis=(0, 2))
        data_sub = data_mean - np.median(data_mean)
        threshold = 7.5 * std_with_clip(data_sub, num_iter=3)
        self.assertTrue(np.all(np.abs(data_sub) < threshold))

    def test_sdatmcor_normal(self):
        '''test normal usage of sdatmcor'''
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data')
        self.check_result()

    def test_sdatmcor_explicit_atmtype(self):
        '''test specifying atmtype explicitly'''
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', atmtype=2)
        self.check_result()

    def test_sdatmcor_overwrite(self):
        '''test overwriting existing outfile'''
        os.mkdir(self.outfile)
        self.assertTrue(os.path.exists(self.outfile))
        sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', overwrite=True)
        self.check_result()

    def test_sdatmcor_no_overwrite(self):
        '''test to avoid overwriting existing outfile'''
        os.mkdir(self.outfile)
        self.assertTrue(os.path.exists(self.outfile))
        # TODO: sdatmcor should throw exception
        self.assertFalse(
            sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn='data', overwrite=False)
        )

    def test_sdatmcor_wrong_datacolumn(self):
        '''test wrong datacolumn'''
        wrong_colnames = ['corrected', 'float_data']
        for colname in wrong_colnames:
            # TODO: sdatmcor should throw exception
            self.assertFalse(
                sdatmcor(infile=self.infile, outfile=self.outfile, datacolumn=colname)
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

import unittest
import os
# import numpy
import sys
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

else:
    from tasks import sdatmcor
    from __main__ import default
    from sdutil import tbmanager

    # Define the root for the data files
    datapath = os.environ.get('CASAPATH').split()[0] + ''

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

    local_unit_test = True

    def setUp(self):

        default(sdatmcor)

        # default Args
        self.args = {'infile': defInputMs,
                     'datacolumn': 'float_data'
                     }
        
    def tearDown(self):
        print("tearDown::deleting MSs.")


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


    def test_param1(self):
        '''sdatmcor::1:: file1  '''

        infile = 'uid___A002_Xe770d7_X320b-t.ms'
        outfile = 'sdatmcor-out1.ms'
        prm = {'debug': 'skipTaskExec',
               'infile' : infile,
               'outfile': outfile,
               'overwrite':  True }
        # Run Task and check
        self.assertTrue(self._run_task(prm))

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

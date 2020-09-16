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
    # datapath = ctsys.resolve('atmcor/')
    datapath = "/work-dev/nishie/ATM-devWork/ATM-data"

else:
    from tasks import sdatmcor
    from __main__ import default
    from sdutil import tbmanager

    # Define the root for the data files
    if False:
        datapath = os.environ.get('CASAPATH').split()[0] + "/atmcor/"
    else:
        datapath = "/work-dev/nishie/ATM-devWork/ATM-data"
#
# Test-MS
#

# template MS (to copy from)
defInputMs = "uid___A002_Xe770d7_X320b.ms"

# testing MS (set up basic conditions)
defWorkMsBasic = "uid___A002_Xe770d7_X320b-t.ms"

# output MS #

defOutputMs = "sdave.ms"        # (internal) output MS
defPrivateMs = "sdave-*.ms"           # (private) Debug output MS template
defPrivateMsForm = 'sdave-{}-{}.ms'   # Debug output MS form


##############
# Test Entry
##############


class test_sdatmcor(unittest.TestCase):
    def setUp(self):
        default(sdatmcor)

        # Check Environment 
        msFile = os.path.join(datapath, 'uid___A002_Xe770d7_X320b.ms')
        if os.path.exists(msFile):
            print( "Test MS available on 'datapath'.")
        else:
            print( "File not Found on this Environment.")
            return 

        # copy template
        self._copy_remote_file(defInputMs, defWorkMsBasic)

        # default Args (based on CASR-424)
        self.args = {'infile': defInputMs,
                     'outfile': defOutputMs,
                     'datacolumn': 'float_data'
                     }
        

    def tearDown(self):
        print("tearDown::deleting MSs.")

        """
        # delete copied in-MS and out-MS
        os.system('rm -rf ' + defInputMs)
        os.system('rm -rf ' + defOutputMs)  # Comment out , for DEBUG ##

        # CAS-13088 TENTATIVE #
        os.system('rm -rf ' + "atm*.ms")
        """
# private function #

    def _copy_remote_file(self, infile, outfile):
        print("Copying temp MS from Remote to Local.")
        os.system('cp -RL ' + os.path.join(datapath, infile) + ' ' + outfile)

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
        """
        os.system('rm -rf ' + defWorkMsBasic)  # in case, the MS already exist.
        os.system('rm -rf ' + defWorkMsTimeSpan)
        os.system('rm -rf ' + defWorkMs3NRO)
        os.system('rm -rf ' + defWorkMs3ALMA)
        os.system('rm -rf ' + defPrivateMs)
        os.system('rm -rf ' + "TEST-*.ms")
        """

    @classmethod
    def tearDownClass(cls):
        print("tearDownClass::deleting work-MS.")
        #
        # Comment Out if you reserve MS.
        #
        """
        os.system('rm -rf ' + defWorkMsBasic)
        os.system('rm -rf ' + defWorkMsTimeSpan)
        os.system('rm -rf ' + defWorkMs3NRO)
        os.system('rm -rf ' + defWorkMs3ALMA)
        os.system('rm -rf ' + defPrivateMs)
        os.system('rm -rf ' + "TEST-*.ms")
        """

##############
# Run Task
##############
    def _run_task(self, auxArgs=None):
        print("_run_task::starts")

        if auxArgs is not None:
            for k in auxArgs:
                self.args[k] = auxArgs[k]

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
        None

##############
# MISC
##############


############################
# TEST FIXTURE
############################


#
# SPW
#

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

    def test_param2(self):
        '''sdatmcor::2:: Execute minimum '''

        infile = 'uid___A002_Xe770d7_X320b-t.ms'
        outfile = 'sdatmcor-out2.ms'
        prm = {'debug'    : 'interruptCorrection,showAtmProfile',
               'infile'   : infile,
               'outfile'  : outfile,
               'overwrite':  True,
               'antenna'  : 'PM02',
               'scan'     : '3,5' }
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

"""

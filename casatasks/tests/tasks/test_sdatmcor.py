import datetime
import re
import unittest
import os
import numpy
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
    datapath = ctsys.resolve('regression/unittest/sdimaging')

else:
    from tasks import sdatmcor
    from __main__ import default
    from sdutil import tbmanager

    # Define the root for the data files
    datapath = os.environ.get('CASAPATH').split(
    )[0] + "/data/regression/unittest/sdimaging/"
"""
sdatmcor begins
"""

#
# Test-MS
#

# template MS (to copy from)
defInputMs = "sdimaging.ms"

# testing MS (set up basic conditions)
defWorkMsBasic = "sdimaging-t.ms"

# testing MS(2) (for 'state','scan' in TimeSpan)
defWorkMsTimeSpan = "sdimaging-t2.ms"

# template-MS(3) (for 'float_data' and 'data') (to copy from)
defInputMs3 = "Uranus1.cal.Ant0.spw34.ms"

# testing MS ('float_data' is used)
defWorkMs3NRO = "Uranus1.cal.Ant0.spw34-nobeyama.ms"

# testing MS ('data' is used instead of 'float_data' )
defWorkMs3ALMA = "Uranus1.cal.Ant0.spw34-ALMA.ms"

# output MS #

defOutputMs = "sdave.ms"        # (internal) output MS
defPrivateMs = "sdave-*.ms"           # (private) Debug output MS template
defPrivateMsForm = 'sdave-{}-{}.ms'   # Debug output MS form

# Test Condition #

numTune = 0                  # must be in {12,24,36...}  and 0(=no operation)
nInScan = 63                 # number of scan (CONST)
nReduce = nInScan * numTune  # nReduce MUST BE even number
nRowOrg = 3843
nRow = nRowOrg - nReduce     # Final Size

# 'scan' and 'state' condition

numOfState = 3                    # test-MS2 (for timespan)
numOfScan = int(nRow / nInScan)   # test-MS2 (for timespan), in sdtimeimaging.ms,  numOfScan=61

# Numerical Error Limit

errLimit = 1.0e-9   # numerical error Limit of ZeroSum. NonZeroSum
errLimit2 = 1.0e-9  # numerical error Limit of Sigma and Weight
testInterval = 1.0   # fundamental INTERVAL in TEST-MS (tunable)

##############
# Test Entry
##############


class test_sdatmcor(unittest.TestCase):
    def setUp(self):
        default(sdatmcor)

        # parameter on self.
        self.inpMs = defInputMs
        self.interval = testInterval

        # copy template
        #### self._copy_remote_file(defInputMs, defInputMs)

        # default Args (based on CASR-424)
        self.args = {'infile': defInputMs,
                     ## 'outfile': defOutputMs,
                     'datacolumn': 'float_data'
                     }
        """
        #
        # create TEST-MS only for the first time.
        #
        if not self._if_exist(defWorkMsBasic):
            print("- TestMS.[{}] being created on current dir.".format(defWorkMsBasic))

            # Copy template and generate Test-MS
            self._copy_remote_file(defInputMs, defWorkMsBasic)

            # Data Generation #
            self. _generate_data(defWorkMsBasic, stateOption=False)

        # create TEST-MS only for the first time.
        #  (for TimeSpan test)
        if not self._if_exist(defWorkMsTimeSpan):
            print("- TestMS(for TimeSpan.[{}] being created on current dir.".format(defWorkMsTimeSpan))

            # Copy template #
            self._copy_remote_file(defInputMs, defWorkMsTimeSpan)

            # Data Generation #
            self. _generate_data(defWorkMsTimeSpan, stateOption=True)

        # create TEST-MS "Only for the first time".
        #  ( using 'data' column, instead of 'float_data' )
        #  ( These MSs are to check for ALMA specific work in mstransform )
        if not self._if_exist(defWorkMs3NRO):
            print("- TestMS(for data/float_data.[{}] being created.".format(defWorkMs3NRO))

            # Copy template #
            self._copy_remote_file(defInputMs3, defWorkMs3NRO)
            self._copy_remote_file(defInputMs3, defWorkMs3ALMA)

            # No Data Generation, but only set TelescopeName
            self._set_telescopename(defWorkMs3NRO, "Nobeyama")  # change name
            self._set_telescopename(defWorkMs3ALMA, "ALMAtest")  # change name
        """

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

    def _checkZero(self, data):
        '''
          check all the results must be zero.
            See _generate_data() to see Test Data.
        '''
        print("-- checking Zero --")

        check = numpy.all(numpy.abs(data) < errLimit)
        self.assertTrue(check, msg='## Zero check Failed ##\n{}'.format(data))

    def _checkNonZero(self, data):
        '''
          check sum of each averaged result
            check non Zero.
        '''
        print("-- checking Non Zero --")

        check = numpy.all(numpy.abs(data) >= errLimit)
        self.assertTrue(
            check,
            msg='## Non Zero check Failed (ref={})  ##\n{}'.format(
                errLimit,
                data))

    def _checkZeroSum(self, data1, data2):
        '''
          check sum of each averaged result
            check Zero.
        '''
        print("-- checking ZeroSum of data1 and data2.")

        zSumData = numpy.abs(numpy.array(data1) + numpy.array(data2))
        check = numpy.all(zSumData < errLimit)
        self.assertTrue(
            check,
            msg='## Zero Sum check Failed (ref={})\n{}'.format(
                errLimit,
                zSumData))

######################
# check time
######################

    def _checkTime(self, msName, row, refTime):
        '''
          check time of specified row,
              compare value with the reference.
        '''
        print("-- checking Time --")

        # get time and inspection.
        self. _get_main(msName)
        # one output
        time = self.tm[row]
        # check Time
        check = (time == refTime)
        self.assertTrue(
            check,
            msg='## Time is Invalid.##\n val={} ref={}'.format(
                time,
                refTime))

######################
# check Output
######################

    def _checkOutputRec(self, msName, refNRow):
        '''
          check calculated result record count,
            compare value with the expected count.
        '''
        print("-- checking Output Record Count --")

        # get time
        self. _get_main(msName)
        # count output rows
        nrow = len(self.tm)
        # check
        check = (nrow == refNRow)
        self.assertTrue(
            check,
            msg='## Row Count in Output is Invalid.##\n val={} ref={}'.format(
                nrow,
                refNRow))

######################
# check scan
######################

    def _check_scan(self, out_msname, refValue):
        '''
          check 'scan'
            number of output must 1,
            compare value with expected value.
        '''
        print("-- checking scan selection --")
        # get table
        self. _get_main(out_msname)
        # get one value from row=0
        scan = self.sc[0]
        # check scan ID
        self.assertTrue(len(self.sc) == 1,
                        msg='## unexpected number of output. ##\n {}'.format(len(self.sc)))
        self.assertTrue(scan == refValue,
                        msg='## unexpected scan no. in output. ##\n {}'.format(scan))

########################
# check Weight/Spectra
########################

    # Check Wait and Sigma
    def _checkWeightSigma(self, msName, row, weight_ref):
        '''
          check Sigma and Weight
            compare 'weight' with expected value.
            Sigma is mathematically inspected.
        '''
        print("-- checking Weight and Sigma --")

        self._get_spectra(msName, row)

        print("Weight Ref :{0}".format(weight_ref))
        print("Weight     :{0}".format(self.wgt))
        print("Sigma      :{0}".format(self.sgm))

        # Check (based on Formula about Sigma and Weight) #
        check1 = (self.wgt[0] == weight_ref)
        check2 = (self.wgt[1] == weight_ref)
        check3 = ((1.0 / self.wgt[0]) -
                  (self.sgm[0] * self.sgm[0]) < errLimit2)
        check4 = ((1.0 / self.wgt[1]) -
                  (self.sgm[1] * self.sgm[1]) < errLimit2)

        # Assert
        self.assertTrue(
            check1, msg='## Weight[0] is unexpected. ##\n {}/{}'.format(self.wgt[0], weight_ref))
        self.assertTrue(
            check2, msg='## Weight[1] is unexpected. ##\n {}/{}'.format(self.wgt[1], weight_ref))
        self.assertTrue(
            check3,
            msg='## Sigma [0] is unexpected. ##\n sigma={}, weight={}'.format(
                self.sgm[0], self.wgt[0]))
        self.assertTrue(
            check4,
            msg='## Sigma [1] is unexpected. ##\n sigma={}, weight={}'.format(
                self.sgm[1], self.wgt[1]))

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
        print("----- Generating MS.")

        # Column Name
        dataColumnName = 'FLOAT_DATA'

        # if non-zero value is set on 'offset', an Intentional Fail is raised.
        offset = 0.0

        # Value parameters
        slope = 0.5     # (tunable, but error threshold subject to be changed.)
        baseTime = 0.0  # gives start time in JD.

        # Table Access (with numPy array operation)
        with tbmanager(msName, nomodify=False) as tb:

            # reduce MS row size if reduce size is specified.
            if (nReduce != 0):
                print("----- reducing rows count in Test-MS.")
                rows = list(range(nReduce))
                tb.removerows(rows)

            # show nRow in MS.
            NN = tb.nrows()  # NN MUST BE same as nRow
            print("Nrow = {}".format(NN))

            # initialize STATE_ID (option)
            #  (ex)  state =0,1,...numOfState-1 ,0,1,.....
            if stateOption:
                print("------ stateOption Active, putting three STATE_IDs on the MS. ")
                arrayState = numpy.mod(numpy.arange(0, NN), numOfState)
                tb.putcol("STATE_ID", arrayState)

            # get array shape of Spectra Data, by getcolshapestring(),
            # returned string is like:"[2,1024]", via list.
            tk = re.split(r",|\[|\]",   # specify delimiter as; [ , ]
                          tb.getcolshapestring(
                              dataColumnName,
                              nrow=1)[0])
            # separating to :: <zero> [<1st> ,<2nd*> ]<3rd>
            nChan = int(tk[2])
            # separating to :: <zero> [<1st*> ,<2nd> ]<3rd>
            # nPol = int(tk[1]) (not used, reserved)

            # create array (time, interval)
            arrayTime = testInterval * \
                numpy.arange(0, NN, dtype=numpy.float64) + baseTime
            arrayInterval = numpy.full(NN, testInterval, dtype=numpy.float64)

            # put to column (from numpy array)
            print("------ Putting Time,INTERVAL.")
            tb.putcol("TIME", arrayTime)
            tb.putcol("INTERVAL", arrayInterval)

            # create Test-Data
            print("------ Calculating Curve.")
            NN1 = (NN - 1) / 2
            L = numpy.linspace(-NN1, NN1, NN) * slope + offset
            VAL = numpy.tile(L, [nChan, 1])
            arrayData3 = numpy.array([VAL, VAL])

            # write to the column at once
            print("------ Putting Curve.")
            tb.putcol(dataColumnName, arrayData3)
        print("------ Done.")

    # set telescope name on MS.
    def _set_telescopename(self, msName, telName):
        print("------ changing Telscope Name. ")
        msObservation = msName + '/OBSERVATION'
        with tbmanager(msObservation, nomodify=False) as tb:
            tb.putcell('TELESCOPE_NAME', 0, telName)
            # tb.resync()

##################################
# sub function for TEST FIXTURE
##################################

    def _check_averaged_result_N1(self, outMsName):
        '''
        This function inspects the Averaged result-MS.
        All the spectral data will be averaged. One averaged result remains, which must be Zero.
        '''
        # get the result and inspect #
        fData = self._get_spectra(outMsName, 0)     # use row=0 from RESULT
        self._checkZero(fData)

        # Ref time
        refTime = (nRow - 1) / 2 * testInterval
        self._checkTime(outMsName, 0, refTime)

        # Weight, Sigma
        self._checkWeightSigma(outMsName, 0, nRow)

    def _check_averaged_result_N3(self, outMsName):
        '''
        This function inspects the Averaged result-MS by 3 averaged results.
        Three sections are Averaged. 1st.result and 3rd.result are different sign and sum =0
        The 2nd section makes Zero Sum.
        Note: In Test-MS, the test data is designed by 'Odd functional' curve.
        '''
        #
        self._get_main(outMsName)
        check = (len(self.tm) == 3)
        self.assertTrue(
            check, msg='## Unexpected Result Count  ##\n count={}'.format(len(self.tm)))

        # get the result  #
        fData0 = self._get_spectra(outMsName, 0)        # result on row=0
        fData1 = self._get_spectra(outMsName, 1)        # row=1
        fData2 = self._get_spectra(outMsName, 2)        # row=2

        # Inspection #
        # Following two sections must be different sign.
        self._checkZeroSum(fData0, fData2)
        self._checkZero(fData1)              # must be zero

        # Ref Time in 3 sections (becomes centre of the section)
        Tref0 = testInterval * (nRow / 3 - 1.0) / 2
        Tref1 = Tref0 + testInterval * (nRow / 3)
        Tref2 = Tref1 + testInterval * (nRow / 3)

        # check Time
        self._checkTime(outMsName, 0, Tref0)
        self._checkTime(outMsName, 1, Tref1)
        self._checkTime(outMsName, 2, Tref2)

        # check Weight, Sigma
        self._checkWeightSigma(outMsName, 0, (nRow / 3))
        self._checkWeightSigma(outMsName, 1, (nRow / 3))
        self._checkWeightSigma(outMsName, 2, (nRow / 3))

    def _check_averaged_result_N3TimeSpan(self, outMsName):
        '''
        This is for TimeSpan (when number of state=3)
        '''
        # check result record count. (must be same as state count)
        self._get_main(outMsName)
        check = (len(self.tm) == 3)
        self.assertTrue(
            check, msg='## Unexpected Result Count. ##\n count={}'.format(len(self.tm)))

        # get the result  #
        fData0 = self._get_spectra(outMsName, 0)        # result on row=0
        fData1 = self._get_spectra(outMsName, 1)        # row=1
        fData2 = self._get_spectra(outMsName, 2)        # row=2

        # Inspection
        # The sum of three section data must be zero (particular in TimeSpan test)
        self._checkZero(fData0 + fData1 + fData2)

        # check Weight, Sigma
        self._checkWeightSigma(outMsName, 0, (nRow / 3))
        self._checkWeightSigma(outMsName, 1, (nRow / 3))
        self._checkWeightSigma(outMsName, 2, (nRow / 3))

    def _check_averaged_result_N61(self, outMsName):
        '''
        This is for TiimeSpan (when scan=state, 61 results are inspected.)
         see numOfScan
        '''
        print("outfile ={} specified.".format(outMsName))
        # check Zero Sum
        for n in range(numOfScan):
            # symmetricaly get the data. These sum must be Zero #
            fData_1 = self._get_spectra(outMsName, n)
            fData_2 = self._get_spectra(outMsName, (numOfScan - 1) - n)
            self._checkZeroSum(fData_1, fData_2)

##############
# MISC
##############

    def _set_outfile_timebin(self, testNo, numRec):

        strTimeBin = '{}s'.format(numRec * testInterval)
        outFile = defPrivateMsForm.format(testNo, numRec)
        return outFile, strTimeBin

############################
# TEST FIXTURE
############################


#
# SPW
#

    def test_param1(self):
        '''sdatmcor::1:: file1  '''

        infile = 'uid___A002_Xe770d7_X320bi.ms'
        outfile = 'sdatmcor-out1.ms'
        prm = {'debug': 'skipTaskExec',
               'infile' : infile,
               'outfile': outfile,
               'overwrite':  True }
        # Run Task and check
        self.assertTrue(self._run_task(prm))

    def test_param2(self):
        '''sdatmcor::2:: file2   '''

        infile = 'uid___A002_Xe77cf7_X2ed5.ms'
        outfile = 'sdatmcor-out2.ms'
        prm = {'debug': 'skipTaskExec',
               'infile' : infile, 
               'outfile': outfile, 
               'overwrite':  True }
        # Run Task and check
        self.assertTrue(self._run_task(prm))

    def test_param3(self):
        '''sdatmcor::3:: file3   '''

        infile = 'uid___A002_Xe78c20_X2dc1.ms'
        outfile = 'sdatmcor-out3.ms'
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

"""

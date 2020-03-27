import datetime
import re       # Regular Expr.
import unittest
import os
import numpy
import sys
from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import sdtimeaverage
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
    from tasks import sdtimeaverage
    from __main__ import default
    from sdutil import tbmanager

    # Define the root for the data files
    datapath = os.environ.get('CASAPATH').split(
    )[0] + "/data/regression/unittest/sdimaging/"
'''''''''''''''''''''''''''
sdtimeaverage begins
'''''''''''''''''''''''''''
# local import

# MS name for basic test
defInputMs = "sdimaging.ms"     # template MS
defWorkMsBasic = "sdimaging-t.ms"   # testing MS (set up basic conditions)
# testing MS (modified 'state' for TimeSpan)
defWorkMsTimeSpan = "sdimaging-t2.ms"

# MS name for 'float_data' and 'data'
defInputMs3 = "Uranus1.cal.Ant0.spw34.ms"        # copied original.
# testing MS (using 'data' instead of 'float_data' )
defWorkMs3NRO = "Uranus1.cal.Ant0.spw34-nobeyama.ms"
# testing MS (using 'data' instead of 'float_data' )
defWorkMs3ALMA = "Uranus1.cal.Ant0.spw34-ALMA.ms"

defOutputMs = "sdave.ms"         # (internal) output MS
defPrivateMs = "sdave-*.ms"       # private  output MS form of each test
defPrivateMsForm = 'sdave-{}-{}.ms'   # debug file format of each test

# Test Condition
numTune = 0                  # must be in {12,24,36...}  and 0(=no operation)
nInScan = 63
nReduce = nInScan * numTune    # nReduce MUST BE even number
nRow = 3843 - nReduce     # Final Size

# 'scan' and 'state' condition
numOfState = 3                     # test-MS2 (for timespan)
numOfScan = int(nRow / nInScan)   # test-MS2 (for timespan) std=61

# Numerical Error Limit
errLimit = 5.0e-08   # numerical error Limit of ZeroSum
errLimit2 = 2.0e-08   # numerical error Limit of Sigma and Weight
testInterval = 1.0    # fundamental INTERVAL in TEST-MS (tunable)

##############
# Test Entry
#############


class test_sdtimeaverage(unittest.TestCase):
    def setUp(self):
        default(sdtimeaverage)

        # copy from master
        self.inpMs = defInputMs
        os.system(
            'cp -RL ' +
            os.path.join(
                datapath,
                self.inpMs) +
            ' ' +
            self.inpMs)

        # output MS
        '''  none '''

        # params
        self.interval = testInterval

        # default Args (based on R424)
        self.args = {'infile': defInputMs,
                     'outfile': defOutputMs,
                     'datacolumn': 'float_data',   # revised CASR-474
                     'timespan': 'scan'
                     }

        # +
        # create TEST-MS only for the first time.
        # -
        filePath = os.path.join("./", defWorkMsBasic)
        if not os.path.exists(filePath):
            print("- TestMS.[{}] being created.".format(filePath))

            # Copy template and generate Test-MS
            os.system(
                'cp -RL ' +
                os.path.join(
                    datapath,
                    defInputMs) +
                ' ' +
                defWorkMsBasic)
            self. generate_data(defWorkMsBasic, stateOption=False)

        # +
        # create TEST-MS only for the first time.
        #  (for TimeSpan test)
        # -
        filePath = os.path.join("./", defWorkMsTimeSpan)
        if not os.path.exists(filePath):
            print("- TestMS(for TimeSpan.[{}] being created.".format(filePath))

            # Copy template and generate Test-MS
            os.system(
                'cp -RL ' +
                os.path.join(
                    datapath,
                    defInputMs) +
                ' ' +
                defWorkMsTimeSpan)
            self. generate_data(defWorkMsTimeSpan, stateOption=True)

        # +
        # create TEST-MS only for the first time.
        #  ( using 'data' column, instead of 'float_data' )
        #  ( this MSs are check for ALMA special, in mstransform )
        # -
        filePath = os.path.join("./", defWorkMs3NRO)
        if not os.path.exists(filePath):
            print(
                "- TestMS(for data/float_data.[{}] being created.".format(filePath))

            # Copy template and generate Test-MS
            os.system(
                'cp -RL ' +
                os.path.join(
                    datapath,
                    defInputMs3) +
                ' ' +
                defWorkMs3NRO)  # ALMA
            os.system(
                'cp -RL ' +
                os.path.join(
                    datapath,
                    defInputMs3) +
                ' ' +
                defWorkMs3ALMA)  # nobeyama

            # No Data Generation, but set TelescopeName
            self.set_telescopename(defWorkMs3NRO, "Nobeyama")  # change name
            self.set_telescopename(defWorkMs3ALMA, "ALMAtest")  # change name

    def tearDown(self):

        # delete copied in-MS and out-MS
        print("tearDown::deleting MSs.")

        os.system('rm -rf ' + self.inpMs)
        os.system('rm -rf ' + defOutputMs)  # Comment out , for DEBUG ##

#
# Class Method
#
    @classmethod
    def setUpClass(cls):
        print("setUpClass::deleting existing work-MS.")
        os.system('rm -rf ' + defWorkMsBasic)  # in case, the MS already exist.
        os.system('rm -rf ' + defWorkMsTimeSpan)
        os.system('rm -rf ' + defWorkMs3NRO)
        os.system('rm -rf ' + defWorkMs3ALMA)

        os.system('rm -rf ' + defPrivateMs)

        os.system('rm -rf ' + "TEST-*.ms")

    @classmethod
    def tearDownClass(cls):
        print("tearDownClass::deleting work-MS.")
        #
        # Comment Out if you reserve MS.
        #
        os.system('rm -rf ' + defWorkMsBasic)  # delete ALL.
        os.system('rm -rf ' + defWorkMsTimeSpan)
        os.system('rm -rf ' + defWorkMs3NRO)
        os.system('rm -rf ' + defWorkMs3ALMA)

        os.system('rm -rf ' + defPrivateMs)

        os.system('rm -rf ' + "TEST-*.ms")

##############
# Run Task
##############
    def run_task(self, auxArgs=None):
        print("run_task::starts")

        if auxArgs is not None:
            for k in auxArgs:
                self.args[k] = auxArgs[k]

        # execution ## open
        try:
            return sdtimeaverage(**self.args)
        except Exception:
            return False

#################
# Check Result
#################

    def checkZero(self, data):
        print("-- checking Zero --")

        check = numpy.all(numpy.abs(data) < errLimit)
        self.assertTrue(check, msg='## Zero check Failed ##\n{}'.format(data))

    def checkNonZero(self, data):
        print("-- checking Not Zero --")

        check = numpy.all(numpy.abs(data) >= errLimit)
        self.assertTrue(
            check,
            msg='## Non Zero check Failed (ref={})  ##\n{}'.format(
                errLimit,
                data))

    def checkZeroSum(self, data1, data2):
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

    def checkTime(self, msName, row, refTime):
        print("-- checking Time --")

        # get time and inspection.
        self. get_main(msName)
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

    def checkOutputRec(self, msName, refNRow):
        print("-- checking Output Record Count --")

        # get time and inspection.
        self. get_main(msName)
        # one output
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

    def check_scan(self, outMsName, refValue):
        print("-- checking scan selection --")
        # get table and inspection.
        self. get_main(defOutputMs)
        # one output
        scan = self.sc[0]
        # check scan ID
        self.assertTrue(len(
            self.sc) == 1, msg='## unexpected number of output##\n {}'.format(len(self.sc)))
        self.assertTrue(
            scan == refValue,
            msg='## unexpected scan no. in output.##\n {}'.format(scan))

########################
# check Weight/Spectra
########################

    # Check Wait and Sigma
    def checkWeightSigma(self, msName, row, weight_ref):
        print("-- checking Weight and Sigma --")

        self.get_spectra(msName, row)

        print("Weight Ref :{0}".format(weight_ref))
        print("Weight     :{0}".format(self.wgt))
        print("Sigma      :{0}".format(self.sgm))

        # Check (based on fomula about Sigma and Waight) #
        check1 = (self.wgt[0] == weight_ref)
        check2 = (self.wgt[1] == weight_ref)
        check3 = ((1.0 / self.wgt[0]) -
                  (self.sgm[0] * self.sgm[0]) < errLimit2)
        check4 = ((1.0 / self.wgt[1]) -
                  (self.sgm[1] * self.sgm[1]) < errLimit2)

        self.assertTrue(
            check1, msg='## Weight[0] is unexpected.##\n {}/{}'.format(self.wgt[0], weight_ref))
        self.assertTrue(
            check2, msg='## Weight[1] is unexpected.##\n {}/{}'.format(self.wgt[1], weight_ref))
        self.assertTrue(
            check3,
            msg='## Sigma [0] is unexpected.##\n {}'.format(
                self.sgm[0]))
        self.assertTrue(
            check4,
            msg='## Sigma [1] is unexpected.##\n {}'.format(
                self.sgm[1]))

##################################
# Read Data from Specified MS
##################################
    # MAIN #
    def get_main(self, msName):
        # get MAIN table data
        with tbmanager(msName) as tb:

            # Key data
            self.tm = tb.getcol('TIME')
            self.a1 = tb.getcol('ANTENNA1')
            self.a2 = tb.getcol('ANTENNA2')
            self.sc = tb.getcol('SCAN_NUMBER')
            self.fd = tb.getcol('FIELD_ID')

    # DATA (spectra) #
    def get_spectra(self, msName, row):
        with tbmanager(msName) as tb:
            # Spectra
            self.data = tb.getcell('FLOAT_DATA', row)
            self.wgt = tb.getcell('WEIGHT', row)
            self.sgm = tb.getcell('SIGMA', row)

        return self.data

################################
# Generate Data on FLOAT_DATA
################################
    def generate_data(self, msName, stateOption=False):
        print("----- Generating MS.")

        # Column Name
        dataColumnName = 'FLOAT_DATA'
        # Test Slope
        # if specified non-zero, an Intensive Fail will be activated.
        offset = 0.0
        slope = 0.5     # (tunable, but error threshold subject to be changed.)
        # (comment)
        baseTime = 0.0
        # Table Access (with numPy array operation)
        with tbmanager(msName, nomodify=False) as tb:
            # +
            # reduce MS row size if reduce size is specified.
            # -
            if nReduce != 0:
                print("----- reducing rows in Test-MS.")
                rows = list(range(nReduce))
                tb.removerows(rows)

            # Confirm nRow
            NN = tb.nrows()  # NN MUST BE same as  nRow
            print("Nrow = {}".format(NN))
            # -----

            # change STATE_ID (option)
            if stateOption:
                print("------ stateOption Active, using three STATE_IDs on the MS. ")
                arrayState = numpy.mod(numpy.arange(0, NN), numOfState)
                tb.putcol("STATE_ID", arrayState)

            # get array shape of Spectra data, by getcolshapestring(),
            # returning string like:"[2,1024]" via 'list'.
            tk = re.split(
                r",|\[|\]",
                tb.getcolshapestring(
                    dataColumnName,
                    nrow=1)[0])  # delimter = [ , ]
            nChan = int(tk[2])  # separating to :: <zero> [<1st> ,<2nd*> ]<3rd>
            # (not used ) separating to :: <zero> [<1st*> ,<2nd> ]<3rd>
            # nPol = int(tk[1]) (not used, reserved)

            # create array (time, interval)
            arrayTime = testInterval * \
                numpy.arange(0, NN, dtype=numpy.float64) + baseTime
            arrayInterval = numpy.full(NN, testInterval, dtype=numpy.float64)

            # put to column
            print("------ Putting Time,INTERVAL.")
            tb.putcol("TIME", arrayTime)
            tb.putcol("INTERVAL", arrayInterval)

            # create Test-Data [use numpy.array]
            print("------ Calculating Curve.")
            NN1 = (NN - 1) / 2
            L = numpy.linspace(-NN1, NN1, NN) * slope + offset
            VAL = numpy.tile(L, [nChan, 1])
            arrayData3 = numpy.array([VAL, VAL])

            # write to the column at once
            print("------ Putting Curve.")
            tb.putcol(dataColumnName, arrayData3)
        print("------ Done.")

    # TELESCOP NAME (7-Feb-2020)

    def set_telescopename(self, msName, telName):
        print("------ changing Telscope Name. ")
        msObservation = msName + '/OBSERVATION'
        with tbmanager(msObservation, nomodify=False) as tb:
            tb.putcell('TELESCOPE_NAME', 0, telName)
            # tb.resync()

##################################
# sub function for TEST FIXTURE
##################################

    def check_averaged_result_N1(self, outMsName):
        '''
        This function inspects the Averaged result-MS
        All the spectral data is averaged. one averaged result remains, which must be Zero.
        '''
        # get the result and inspect #
        fData = self.get_spectra(outMsName, 0)     # use row=0 from RESULT
        self.checkZero(fData)

        # Ref time
        refTime = (nRow - 1) / 2 * testInterval
        self.checkTime(outMsName, 0, refTime)

        # Weight, Sigma
        self.checkWeightSigma(outMsName, 0, nRow)

    def check_averaged_result_N3(self, outMsName):
        '''
        This function inspects the Averaged result-MS by 3 averaged results.
        Three sections are Averaged. 1st.result and 3rd.result are different sign and sum =0
        The 2nd. section makes Zero Sum.
        Note: In Test-MS test data is designed by odd functional curve.
        '''
        #
        self.get_main(outMsName)
        check = (len(self.tm) == 3)
        self.assertTrue(
            check, msg='## Unexpected Result Count  ##\n count={}'.format(len(self.tm)))

        # get the result  #
        fData0 = self.get_spectra(outMsName, 0)        # result on row=0
        fData1 = self.get_spectra(outMsName, 1)        # row=1
        fData2 = self.get_spectra(outMsName, 2)        # row=2

        # inspection #
        # must be different sign = Zero Sum.
        self.checkZeroSum(fData0, fData2)
        self.checkZero(fData1)              # must be zero

        # Ref Time in 3 sections (becomes centre of the section)
        print("DBG in check_averaged_result_N3:: nRow = {} ".format(nRow))
        Tref0 = testInterval * (nRow / 3 - 1.0) / 2
        Tref1 = Tref0 + testInterval * (nRow / 3)
        Tref2 = Tref1 + testInterval * (nRow / 3)

        # check Time
        self.checkTime(outMsName, 0, Tref0)
        self.checkTime(outMsName, 1, Tref1)
        self.checkTime(outMsName, 2, Tref2)

        # check Weight, Sigma
        self.checkWeightSigma(outMsName, 0, (nRow / 3))
        self.checkWeightSigma(outMsName, 1, (nRow / 3))
        self.checkWeightSigma(outMsName, 2, (nRow / 3))

    def check_averaged_result_N3TimeSpan(self, outMsName):
        '''
        This is for TimeSpan (num of state=3)
        '''
        #
        self.get_main(outMsName)
        check = (len(self.tm) == 3)
        self.assertTrue(
            check, msg='## Unexpected Result Count  ##\n count={}'.format(len(self.tm)))

        # get the result  #
        fData0 = self.get_spectra(outMsName, 0)        # result on row=0
        fData1 = self.get_spectra(outMsName, 1)        # row=1
        fData2 = self.get_spectra(outMsName, 2)        # row=2

        # inspection #
        # must be zero (TimeSpan test particular)
        self.checkZero(fData0 + fData1 + fData2)

        # check Time (No action)
        pass

        # check Weight, Sigma
        self.checkWeightSigma(outMsName, 0, (nRow / 3))
        self.checkWeightSigma(outMsName, 1, (nRow / 3))
        self.checkWeightSigma(outMsName, 2, (nRow / 3))

    def check_averaged_result_N61(self, outMsName):
        '''
        This is for TiimeSpa (scan=state, 61 results are inspected.
        '''
        print("outfile ={} specified.".format(outMsName))
        # check Zero Sum
        for n in range(numOfScan):
            # symmetricaly get the data. These sum must be Zero #
            fData_1 = self.get_spectra(outMsName, n)
            fData_2 = self.get_spectra(outMsName, (numOfScan - 1) - n)
            self.checkZeroSum(fData_1, fData_2)

##############
# MISC
##############

    def setOutfile_Timebin(self, testNo, numRec):

        strTimeBin = '{}s'.format(numRec * testInterval)
        outFile = defPrivateMsForm.format(testNo, numRec)
        return outFile, strTimeBin

############################
# TEST FIXTURE
############################

# TIME RANGE #
    def test_param00(self):
        '''sdtimeagerage::00:: timerange = 00:00:00~01:04:03 NORMAL (3843s same as in MS)'''

        # set timebin string and private outputMS name.
        privateOutfile, dmy = self.setOutfile_Timebin(0, 3844)
        # Time string (justify to interval and row counts)
        td = datetime.timedelta(seconds=testInterval * nRow)
        timerangeString = '00:00:00~' + str(td)  # '01:04:03' =3848 #
        # Run Task
        prm = {'timerange': timerangeString,
               'timebin': '',
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}  # Specify Full-Range #
        self.run_task(prm)
        # Check Result (zerosum check)
        self.get_spectra(privateOutfile, 0)   # row=0
        self.checkZero(self.data)
        self.checkOutputRec(privateOutfile, 1)

    def test_param01E(self):
        '''sdtimeagerage::01E:: timerange = 00:00:00~01:00:00 ERROR case(3600s INSUFFICIENT)'''

        # set timebin string and private outputMS name.
        privateOutfile, dmy = self.setOutfile_Timebin(1, nRow)
        # Time string (justify to interval and row counts)
        td = datetime.timedelta(
            seconds=nRow *
            testInterval /
            10)   # relatively very short
        timerangeString = '00:00:00~' + str(td)
        # Run Task
        prm = {'timerange': timerangeString,  # orig set up = '00:00:00~00:4:00',
               'timebin': '',
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}  # Specify Full-Range #
        self.run_task(prm)
        # Check Result (zerosum check)
        self.get_spectra(privateOutfile, 0)   # row=0
        self.checkNonZero(self.data)
        self.checkOutputRec(privateOutfile, 1)

    def test_param02(self):
        '''sdtimeagerage::02":: timerange = ""   (dafault) '''

        # set timebin string and private outputMS name.
        privateOutfile, dmy = self.setOutfile_Timebin(2, nRow + 1)
        # Run Task
        prm = {'timerange': '',
               'timebin': '',
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}  # Specify Full-Range #
        self.run_task(prm)
        # Check Result (zerosum check)
        self.get_spectra(privateOutfile, 0)   # row=0
        self.checkZero(self.data)
        self.checkOutputRec(privateOutfile, 1)

# SCAN #
    def test_para10(self):
        '''sdtimeagerage::10:: scan=2 (Within the range)'''

        # Run Task
        prm = {'timebin': '',
               'scan': '2'}
        self.run_task(prm)
        # check scan
        self.check_scan(defOutputMs, 2)
        self.checkOutputRec(defOutputMs, 1)

    def test_param11a(self):
        '''sdtimeagerage::11a:: scan=61 (Within the range)'''

        # Run Task
        prm = {'timebin': '',
               'scan': '61'}  # Normal. In range. #
        self.assertTrue(self.run_task(prm))
        # check scan
        self.check_scan(defOutputMs, 61)
        self.checkOutputRec(defOutputMs, 1)

    def test_param11b(self):
        '''sdtimeagerage::11b:: scan=62 (Error Out of range) '''

        # Run Task
        prm = {'timebin': '',
               'scan': '62'}  # ERROR : out of range in MS #
        self.assertFalse(self.run_task(prm))
        # check scan
        #    no check ...

    def test_param12(self):
        '''sdtimeagerage::12:: scan='' (no number) Default action. '''

        # set timebin string and private outputMS name.
        privateOutfile, timebin_str = self.setOutfile_Timebin(12, nRow + 3)

        prm = {'timebin': '',
               'scan': '',
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}
        # Run Task
        self.run_task(prm)

        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile)
        self.checkOutputRec(privateOutfile, 1)

# FIELD #
    def test_param20(self):
        '''sdtimeaverage::20:: field = 'FLS3a*' (Exact NAME)'''

        prm = {'field': 'FLS3a*'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param21E(self):
        '''sdtimeaverage::21E:: field = 'hoge*' (Error :Bad NAME)'''

        prm = {'field': 'hoge'}
        # Run Task and check
        self.assertFalse(self.run_task(prm))

    def test_param22(self):
        '''sdtimeaverage::22:: field = '' (OK :use default)'''

        prm = {'field': '*'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

# SPW #
    def test_param30(self):
        '''sdtimeaverage::30:: spw = '1' (exist)'''

        prm = {'spw': '0'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param31E(self):
        '''sdtimeaverage::31E:: spw = '9' (Error: Not Exist)'''

        prm = {'spw': '9'}
        # Run Task and check
        self.assertFalse(self.run_task(prm))

    def test_param32(self):
        '''sdtimeaverage:: spw = '' (default)'''

        prm = {'spw': ''}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param33(self):
        '''sdtimeaverage::33:: spw = '*' (OK: Wildcard)'''

        prm = {'spw': ''}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

# ANTENNA #
    def test_param40(self):
        '''sdtimeaverage::40:: antenna = 'GBT' (Exact name)'''

        prm = {'antenna': 'GBT'}
        # Run Task and check
        self.assertTrue(self.run_task(prm), msg="Error in test_param1")
        self.checkOutputRec(defOutputMs, 1)

    def test_param41E(self):
        '''sdtimeaverage::41E antenna = 'gBT' (Error: Bad name) '''

        prm = {'antenna': 'gBT'}
        # Run Task and check
        self.assertFalse(self.run_task(prm))  # must be false

# TIMEBIN ,make Average #
    def test_param100(self):
        '''sdtimeaverage::100:: timebin=1282(N=3)  '''
        # set timebin string and private outputMS name.
        privateOutfile, timebin_str = self.setOutfile_Timebin(
            100, nRow / 3 + 0.5)  # if reduced cond. 0.0 is needed #

        prm = {'timebin': timebin_str,
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}
        # Run Task
        self.run_task(prm)

        # Check Result (zerosum check)
        self.check_averaged_result_N3(privateOutfile)
        self.checkOutputRec(privateOutfile, 3)

    def test_param101(self):
        '''sdtimeaverage::101: timebin=3846(N=1), timebin=''  '''

        # set timebin string and private outputMS name.
        privateOutfile, timebin_str = self.setOutfile_Timebin(101, nRow + 3)

        prm = {'timebin': timebin_str,             # Immediate Value ,
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}
        # Run Task
        self.run_task(prm)

        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile)
        self.checkOutputRec(privateOutfile, 1)

    def test_param103(self):
        '''sdtimeaverage::103: timebin=3846(N=1), timebin='all'  '''

        # set timebin string and private outputMS name.
        privateOutfile, dmy = self.setOutfile_Timebin(103, nRow + 3)

        prm = {'timebin': 'all',                # default = all is applied.
               'infile': defWorkMsBasic,
               'outfile': privateOutfile}
        # Run Task
        self.run_task(prm)

        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile)
        self.checkOutputRec(privateOutfile, 1)

# TIMEBIN #
    def test_param110(self):
        '''sdtimeagerage::110:: timebin='all' '''

        prm = {'timebin': 'all'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param111(self):
        '''sdtimeagerage::111:: timebin='ALL' '''

        # Run Task
        prm = {'timebin': 'ALL'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param112(self):
        '''sdtimeagerage::112:: timebin='' (default)    '''

        # Run Task
        prm = {'timebin': ''}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param113E(self):
        '''sdtimeagerage::113E:: timebin='Alles' (ERROR: Bad keyword)    '''

        # Run Task
        prm = {'timebin': 'Alles'}
        # Run Task and check
        self.assertFalse(self.run_task(prm))  # Error expected #

    def test_param114(self):
        '''sdtimeagerage::114:: timebin='aLL' (OK: Upper/Lower case mixed)    '''

        # Run Task
        prm = {'timebin': 'aLL'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.checkOutputRec(defOutputMs, 1)

    def test_param115(self):
        '''sdtimeagerage::115:: timebin='0' (No averaging, not an Error)    '''

        # Run Task
        prm = {'timebin': '0'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))

    def test_param116(self):
        '''sdtimeagerage::115:: timebin='-1' (Error. Not acceptable)    '''

        # Run Task
        prm = {'timebin': '-1'}
        # Run Task and check
        self.assertFalse(self.run_task(prm))

# +
# DATACOLUMN (Testing alternative collumn select )
# -

    def test_param50(self):
        '''sdtimeaverage::50:: MS= 'float_data'   arg = 'float_data' (NORMAL)   '''

        prm = {'infile': defWorkMs3NRO,
               'outfile': "TEST-50.ms",
               'datacolumn': 'float_data'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))

    def test_param51(self):
        '''sdtimeaverage::51:: MS= 'float_data'    arg = 'data' Column Swich.  '''

        prm = {'infile': defWorkMs3NRO,
               'outfile': "TEST-51.ms",
               'datacolumn': 'data'}

        # Run Task and check
        self.assertTrue(self.run_task(prm))  # Alternative Use will work.

    def test_param52(self):
        '''sdtimeaverage::52:: MS= 'data'    arg = 'float_data' Column Swich '''

        prm = {'infile': defWorkMs3ALMA,
               'outfile': "TEST-52.ms",
               'datacolumn': 'float_data'}

        # Run Task and check
        self.assertTrue(self.run_task(prm))  # must be false

    def test_param53(self):
        '''sdtimeaverage::53:: MS= 'data'    arg = 'data'   (NORMAL)  '''

        prm = {'infile': defWorkMs3ALMA,
               'outfile': "TEST-53.ms",
               'datacolumn': 'data'}

        # Run Task and check
        self.assertTrue(self.run_task(prm))  # must be false

# +
# ALMA option in mstransform
# -
    def test_param60Nobeyama(self):
        '''sdtimeaverage::60 Nobeyama:: simply 'scan' is applied  '''

        privateOutfile = "TEST-60-Nobeyama.ms"
        prm = {'infile': defWorkMs3NRO,
               'outfile': privateOutfile,
               'timebin': 'all',
               'timespan': 'scan'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        # see directly about infile for detail.
        self.checkOutputRec(privateOutfile, 121)

    def test_param60ALMA(self):
        '''sdtimeaverage::61:: ALMA  'scan, state' is applied. '''

        privateOutfile = "TEST-61-ALMA.ms"
        prm = {'infile': defWorkMs3ALMA,
               'outfile': privateOutfile,
               'timebin': 'all',
               'timespan': 'scan'}
        # Run Task and check
        self.assertTrue(self.run_task(prm))
        # 'scan, state' is applied. number of result = 1
        self.checkOutputRec(privateOutfile, 1)

# +
# TimeSpan
# -

    # ordinary behavior #
    def test_param70(self):
        '''sdtimeaverage::70:: timespan="scan"  '''

        privateOutfile, timebin_str = self.setOutfile_Timebin(70, nRow + 3)
        prm = {'infile': defWorkMsTimeSpan,
               'timespan': 'scan',
               'outfile': privateOutfile}

        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.check_averaged_result_N3TimeSpan(privateOutfile)
        # Averaged by each State={0,1,2} . see generate_data()
        self.checkOutputRec(privateOutfile, numOfState)

    # ordinary behavior #
    def test_param71(self):
        '''sdtimeaverage::71:: timespan="state"  '''

        privateOutfile, timebin_str = self.setOutfile_Timebin(71, nRow + 3)
        prm = {'infile': defWorkMsTimeSpan,
               'timespan': 'state',
               'outfile': privateOutfile}

        # Run Task and check
        self.assertTrue(self.run_task(prm))

        # Check Result (zerosum check)
        self.check_averaged_result_N61(privateOutfile)
        self.checkOutputRec(privateOutfile, numOfScan)

    # ordinary behavior #
    def test_param72(self):
        '''sdtimeaverage::72:: timespan="scan,state" (WARN) in NRO '''

        privateOutfile, timebin_str = self.setOutfile_Timebin(72, nRow + 3)
        prm = {'infile': defWorkMsTimeSpan,
               'timespan': 'scan,state',  # HERE is the point ##
               'outfile': privateOutfile}

        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.check_averaged_result_N1(privateOutfile)
        self.checkOutputRec(privateOutfile, 1)

    # ordinary behavior #
    def test_param73(self):
        '''sdtimeaverage::73:: timespan="state,scan" (WARN) in NRO '''

        privateOutfile, timebin_str = self.setOutfile_Timebin(73, nRow + 3)
        prm = {'infile': defWorkMsTimeSpan,
               'timespan': 'state,scan',  # HERE is the point ##
               'outfile': privateOutfile}

        # Run Task and check
        self.assertTrue(self.run_task(prm))
        self.check_averaged_result_N1(privateOutfile)
        self.checkOutputRec(privateOutfile, 1)

    # Error in sytax #
    def test_param74E(self):
        '''sdtimeaverage::79E:: timespan="hoge"  '''

        privateOutfile, timebin_str = self.setOutfile_Timebin(79, nRow + 3)
        prm = {'infile': defWorkMsTimeSpan,
               'timespan': 'hoge',
               'outfile': privateOutfile}

        # Run Task and check
        self.assertTrue(self.run_task(prm))


"""
5-Feb-2020   built a function to reduce the TEST-MS size to shorten execution of mstransform.
14-Feb-2020  editing 'data' column switch Test.
20-Feb-2020  Pushed. Merged master. Re-Buid.
23-Mar-2020  Reformed by autopep8 and manyally formed some.
"""
# Control #


def suite():
    return [test_sdtimeaverage]


if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

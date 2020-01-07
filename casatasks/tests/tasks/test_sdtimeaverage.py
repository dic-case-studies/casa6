import shutil
import unittest
import os
import numpy
import math
import sys
import filecmp
import glob
from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
#   from casatasks import nrobeamaverage
    from casatasks import sdtimeaverage
    from casatools import ms
    from casatools import table

    # default isn't used in casatasks
    def default(atask):
        pass

    ### for testhelper import
    sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
    import testhelper as th
    from casatasks.private.sdutil import tbmanager
    from casatools import ctsys
    datapath=ctsys.resolve('regression/unittest/sdimaging')

else:
#   from tasks import nrobeamaverage
    from tasks import sdtimeaverage
    from taskinit import mstool as ms
    from taskinit import tbtool as table

    from __main__ import default
    import testhelper as th
    from sdutil import tbmanager

    # Define the root for the data files
    datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/sdimaging/"
'''''''''''''''''''''''''''
sdtimeaverage begins
'''''''''''''''''''''''''''
# MS name for this test
defInputMs  = "sdimaging.ms"
defWorkMs   = "sdimaging-t.ms"
defOutputMs = "bave.ms"

# Compare err limit , ideally vector(1024 x 2) is the best
nRow     = 3843  ## DO NOT CHANGE ## 
errLimit  = 2.0e-08   # numerical error Limit of ZeroSum
errLimit2 = 5.0e-08   # numerical error Limit of Sigma and Weight
testInterval = 1.0      # fundamental INTERVAL in TEST-MS

##############
# Test Entry
#############
class test_sdtimeaverage(unittest.TestCase):
    def setUp(self):
        default(sdtimeaverage)

        # copy from master
        self.inpMs = defInputMs
        os.system('cp -RL '+ os.path.join(datapath,self.inpMs) +' '+ self.inpMs)

        # output MS 
        '''  none '''
  
        # params
        self.interval = testInterval

        # default Args (minimum)
        self.args = {'infile'     :  defInputMs,
                     'outfile'    :  defOutputMs,
                     'datacolumn' :  'float_data'    # CASR-474 (float ->data) 
                    }

        # create TEST-MS only for the first time.  
        filePath=os.path.join( "./", defWorkMs)
        if  not os.path.exists(filePath):
            print ( "- TestMS.[{}] being created.".format(filePath)    ) 

            # Copy template and generate Test-MS
            os.system('cp -RL '+ os.path.join(datapath, defInputMs) +' '+ defWorkMs)
            self. generate_data( defWorkMs )

        else:
            print( "- TestMS already exists.")

    def tearDown(self):

        # delete copied in-MS and out-MS
        print( "tearDown::deleting MSs.")

        os.system('rm -rf ' + self.inpMs )
        os.system('rm -rf ' + defOutputMs )   ## Comment out , for DEBUG ##
        os.system('rm -rf ' + "bave*.ms" )
        return

#
# Class Method
#
    @classmethod
    def setUpClass(cls):
        print( "setUpClass::deleting existing work-MS.")
        os.system('rm -rf ' + defWorkMs ) # in case, the MS already exist.

    @classmethod
    def tearDownClass(cls):
        print( "tearDownClass::deleting work-MS.")
        os.system('rm -rf ' + defWorkMs )

##############
# Run Task
##############
    def run_task(self, auxArgs=None):
        print( "run_task::starts" )

        if auxArgs is not None:
            for k in auxArgs: self.args[k] = auxArgs[k]

        ## execution ## 
        try:
            return sdtimeaverage(**self.args)
        except Exception as instance:
            return False


#################
# Check Result
#################

    def checkZero(self,data):
        print("-- checking Zero --")

        check = numpy.all(numpy.abs(data) < errLimit)
        self.assertTrue(check, msg='## Zero check Failed ##\n{}'.format(data)   )
        return True;

    def checkNonZero(self,data):
        print("-- checking Not Zero --")

        check = numpy.all(numpy.abs(data) >= errLimit)
        self.assertTrue(check, msg='## Non Zero check Failed ##\n{}'.format(data)   )
        return True;

    def checkZeroSum(self, data1,data2):
        print("-- checking ZeroSum of data1 and data2." )

        zSumData = numpy.abs(numpy.array(data1) + numpy.array(data2))
        check = numpy.all(zSumData < errLimit)
        self.assertTrue(check, msg='## Zero Sum check Failed ##\n{}'.format(zSumData)   )
        return True

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
        check = ( time == refTime )
        self.assertTrue(check, msg='## Time is not equal with expect.##\n{}'.format(time) )

######################
# check scan
######################

    def check_scan (self, outMsName, refValue ):

        # get table and inspection. 
        self. get_main(defOutputMs)
        # one output
        scan = self.sc[0]
        # check scan ID
        self.assertTrue (len(self.sc)==1 , msg='## unexpected number of output##\n {}'.format(len(self.sc)) )
        self.assertTrue (scan == refValue, msg='## unexpected scan no. in output.##\n {}'.format(scan) )

########################
# check Weight/Spectra
########################

    # Chck Wait and Sigma
    def checkWeightSigma(self, msName, row, weight_ref ):
        print( "-- checking Weight and Sigma --")

        self.get_spectra(msName, row )

        print( "Weight Ref :{0}".format(weight_ref) )
        print( "Weight     :{0}".format(self.wgt)   )
        print( "Sigma      :{0}".format(self.sgm)  )

        # check #
        check1 =  (self.wgt[0] == weight_ref)
        check2 =  (self.wgt[1] == weight_ref)
        check3 =  ( (1.0/self.wgt[0])  - (self.sgm[0] * self.sgm[0])  < errLimit2 )
        check4 =  ( (1.0/self.wgt[1])  - (self.sgm[1] * self.sgm[1])  < errLimit2 )

        self.assertTrue(check1, msg='## Weight[0] is unexpected.##\n {}'.format(self.wgt[0]) )
        self.assertTrue(check2, msg='## Weight[1] is unexpected.##\n {}'.format(self.wgt[1]) )
        self.assertTrue(check3, msg='## Sigma [0] is unexpected.##\n {}'.format(self.sgm[0]) )
        self.assertTrue(check4, msg='## Sigma [1] is unexpected.##\n {}'.format(self.sgm[1]) )

#########################
# Generating Test Data 
#########################

#+
# Read Data from Specified MS
#-
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
    def get_spectra(self,msName, row ):
        with tbmanager(msName) as tb:
            # Spectra 
            self.data = tb.getcell('FLOAT_DATA',row)
            self.wgt  = tb.getcell('WEIGHT', row)
            self.sgm  = tb.getcell('SIGMA', row)           
        
        return self.data 

################################
# Generate DATa on FLOAT_DATA
################################
    def generate_data( self, msName ):
        print( "-- Generating MS." )
        self. get_main(defInputMs )
        # Test Slope
        offset = 0.0        # if specified non-zero, intensive fail can be cauesed.
        slope  = 0.0001
        # Time
        baseTime   = 0.0
        # Table Access
        with tbmanager(msName,nomodify=False) as tb:
            # Data Buffer
            arrayData2 = tb.getcol('FLOAT_DATA')
            # create array (time, interval)
            NN = len(self.tm)

            ''' Revised to Array operation '''
            arrayTime = testInterval * numpy.arange(0,NN,dtype=numpy.float64) + baseTime
            arrayInterval = testInterval * numpy.ones(NN,dtype=numpy.float64)
        
            # put to column
            tb.putcol("TIME",       arrayTime  )
            tb.putcol("INTERVAL",   arrayInterval  )

            # create Test-Data
            for row in range(NN):
                for n in range(1024):
                    # values
                    x = row - numpy.floor(nRow/2)
                    val = offset + slope * x
                    # set to Buffer 
                    arrayData2[0][n][row] = val
                    arrayData2[1][n][row] = val
            # write to the column at once
            tb.putcol("FLOAT_DATA",   arrayData2  )
          
        return True          

##################################
# sub function for TEST FIXTURE
##################################

    def check_averaged_result_N1(self, outMsName):
        '''
        this function inspect the Averaged result-MS
        '''
        # get the result and inspect #
        fData = self.get_spectra(outMsName, 0 )     # use row=0 from RESULT 
        self.checkZero( fData )    

        # Ref time 
        refTime = (nRow -1)/2 * testInterval
        self.checkTime(outMsName, 0, refTime)

        # Weight, Sigma 
        self.checkWeightSigma(outMsName, 0, nRow )

    def check_averaged_result_N3(self, outMsName):

        # get the result and inspect #
        fData0 = self.get_spectra(outMsName, 0 )        # result on row=0
        fData1 = self.get_spectra(outMsName, 1 )        # 
        fData2 = self.get_spectra(outMsName, 2 )        # 

        # inspection # 
        self.checkZeroSum( fData0, fData2 )   # must be different sign
        self.checkZero( fData1 )              # must be zero

        # Ref Time in 3 sections
        Tref0 = testInterval * (nRow /3 -1.0)/2
        Tref1 = Tref0 + testInterval * (nRow/3)
        Tref2 = Tref1 + testInterval * (nRow/3)

        # check Time
        self.checkTime(outMsName, 0, Tref0)
        self.checkTime(outMsName, 1, Tref1)
        self.checkTime(outMsName, 2, Tref2)

        # Weight, Sigma 
        self.checkWeightSigma(outMsName, 0, (nRow/3) )
        self.checkWeightSigma(outMsName, 1, (nRow/3) )
        self.checkWeightSigma(outMsName, 2, (nRow/3) )

##############
# MISC
##############

    def setOutfile_Timebin(self, testNo, numRec ):
        
        strTimeBin = '{}s'.format(numRec * testInterval)
        outFile    = 'bave-{}-{}.ms'.format(testNo, numRec)
        return outFile, strTimeBin;

############################
# TEST FIXTURE
############################

## TIMESPAN ###
    def test_param00(self):
        '''sdtimeagerage::00:: timerange = 00:00:00~01:04:03 NORMAL (3843s same as in MS)'''

        # set timebin string and private outputMS name.
        privateOutfile, dmy  = self.setOutfile_Timebin( 0, 3844 )
        # Run Task
        prm =  {'timerange' : '00:00:00~01:04:03',
                'timebin'   : '',
                'infile'    : defWorkMs,
                'outfile'   : privateOutfile    } # Specify Full-Range #
        self.run_task( prm )
        # Check Result (zerosum check)
        self.get_spectra(privateOutfile, 0 )   # row=0
        self.checkZero (self.data )

    def test_param01E(self):
        '''sdtimeagerage::01E:: timerange = 00:00:00~01:00:00 ERROR case(3600s INSUFFICIENT)'''

        # set timebin string and private outputMS name.
        privateOutfile, dmy  = self.setOutfile_Timebin( 1, 3600 )
        # Run Task
        prm =  {'timerange' : '00:00:00~01:00:00',
                'timebin'   : '',
                'infile'    : defWorkMs,
                'outfile'   : privateOutfile    } # Specify Full-Range #
        self.run_task( prm )
        # Check Result (zerosum check)
        self.get_spectra(privateOutfile, 0 )   # row=0
        self.checkNonZero (self.data )

    def test_param02(self):
        '''sdtimeagerage::02":: timerange = ""   (dafault) '''

        # set timebin string and private outputMS name.
        privateOutfile, dmy  = self.setOutfile_Timebin( 2, 3844 )
        # Run Task
        prm =  {'timerange' : '',
                'timebin'   : '',
                'infile'    : defWorkMs,
                'outfile'   : privateOutfile    } # Specify Full-Range #
        self.run_task( prm )
        # Check Result (zerosum check)
        self.get_spectra(privateOutfile, 0 )   # row=0
        self.checkZero (self.data )

## SCAN ###
    def test_para10(self):
        '''sdtimeagerage::10:: scan=2 '''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '2'   }
        self.run_task( prm )
        # check scan 
        self.check_scan (defOutputMs, 2)

    def test_param11a(self):
        '''sdtimeagerage::11a:: scan=61  OK, In Range'''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '61'   } # Normal. In range. #
        self.assertTrue(self.run_task( prm ))
        # check scan 
        self.check_scan (defOutputMs, 61)

    def test_param11b(self):
        '''sdtimeagerage::11b:: scan=62 Error out Of Range  '''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '62'   } # ERROR : out of range in MS #
        self.assertFalse(self.run_task( prm ))
        # check scan 
        #    no check ...

    def test_param12(self):
        '''sdtimeagerage::12:: scan='' N=1 OK  '''

        # set timebin string and private outputMS name.
        privateOutfile, timebin_str  = self.setOutfile_Timebin( 42, 3846 )

        prm =  {'timebin' : '' ,
                'scan'    : '' ,
                'infile'  : defWorkMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile)

## FIELD ###
    def test_param20(self):
        '''sdtimeaverage::20:: field = 'FLS3a*' (EXACT NAME)'''

        prm =  {'field' : 'FLS3a*'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param21E(self):
        '''sdtimeaverage::21E:: field = 'hoge*' (NG NAME)'''

        prm =  {'field' : 'hoge'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm ))

    def test_param22(self):
        '''sdtimeaverage::22:: field = '' (OK :use default)'''

        prm =  {'field' : '*'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

## SPW ###
    def test_param30(self):
        '''sdtimeaverage::30:: spw = '1' (exist)'''

        prm =  {'spw' : '0'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param31E(self):
        '''sdtimeaverage::31E:: spw = '9' (Not Exist)'''

        prm =  {'spw' : '9'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm ))

    def test_param32(self):
        '''sdtimeaverage:: spw = '' (default)'''

        prm =  {'spw' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param33(self):
        '''sdtimeaverage::33:: spw = '*' (Wildcard)'''

        prm =  {'spw' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

## ANTENNA ###
    def test_param40(self): 
        '''sdtimeaverage::40:: antenna = 'GBT' ''' 

        prm =  {'antenna' : 'GBT'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm), msg="Error in test_param1")

    def test_param41E(self):
        '''sdtimeaverage::41E antenna = 'gBT' (Error) '''

        prm =  {'antenna'    : 'gBT' }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # must be false

## TIMEBIN ,make Average ###
    def test_param100(self):
        '''sdtimeaverage::100:: timebin=1282(N=3)  '''
        # set timebin string and private outputMS name.
        privateOutfile, timebin_str  = self.setOutfile_Timebin( 20, 1282 )
 
        prm =  {'timebin' : timebin_str,          
                'infile'  : defWorkMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # Check Result (zerosum check)
        self.check_averaged_result_N3(privateOutfile)

    def test_param101(self):
        '''sdtimeaverage::101: timebin=3846(N=1), timebin=''  '''

        # set timebin string and private outputMS name.
        privateOutfile, timebin_str  = self.setOutfile_Timebin( 21, 3846 )

        prm =  {'timebin' : timebin_str,             # Immediate Value ,
                'infile'  : defWorkMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )
     
        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile) 

    def test_param103(self):
        '''sdtimeaverage::103: timebin=3846(N=1), timebin='all'  '''

        # set timebin string and private outputMS name.
        privateOutfile, dmy  = self.setOutfile_Timebin( 22, 3846 )

        prm =  {'timebin' : 'all',                # default = all is applied.
                'infile'  : defWorkMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile)

## TIMEBIN ###
    def test_param110(self): 
        '''sdtimeagerage::110:: timebin='all' ''' 

        prm =  {'timebin' : 'all'   }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
 
    def test_param111(self):
        '''sdtimeagerage::111:: timebin='ALL' '''

        # Run Task
        prm =  {'timebin' : 'ALL'   }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param112(self):
        '''sdtimeagerage::112:: timebin='' (default)    '''

        # Run Task
        prm =  {'timebin' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param113E(self):
        '''sdtimeagerage::113E:: timebin='Alles' (ERROR)    '''

        # Run Task
        prm =  {'timebin' : 'Alles'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # Error expected #  

    def test_param114(self):
        '''sdtimeagerage::114:: timebin='aLL' (Capital mixed)=OK    '''

        # Run Task
        prm =  {'timebin' : 'aLL'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

## DATACOLUMN ###
    def test_param50(self):
        '''sdtimeaverage::50:: datacolumn = 'float_data' '''

        prm =  {'datacolumn' : 'float_data'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param51E(self):
        '''sdtimeaverage::51E:: datacolumn = 'data' (Error) '''

        prm =  {'datacolumn' : 'data' }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # must be false

    def test_param52E(self):
        '''sdtimeaverage::52E:: datacolumn = '' default=data is applied, only in this test,makes Error. '''

        prm =  {'datacolumn' : '' }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # must be false

#### Control ######

def suite():
    return [test_sdtimeaverage]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

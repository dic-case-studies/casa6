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

# MS name for this test
defInputMs  = "sdimaging.ms"
defWorkMs   = "sdimaging-t.ms"
defOutputMs = "bave.ms"

# Compare err limit , ideally vector(1024 x 2) is the best
nRow     = 3843  ## DO NOT CHANGE ## 
errLimit  = 5e-08   # numerical error Limit of ZeroSum
errLimit2 = 5e-08   # numerical error Limit of Sigma and Weight
interval_0 = 1.0    # fundamental INTERVAL in TEST-MS

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
        self.interval = interval_0

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

#        os.system('rm -rf ' + defWorkMs )
        os.system('rm -rf ' + "bave*.ms" )
        return

#
# Class Method
#
    @classmethod
    def setUpClass(cls):
        print( "setUpClass::deleting work-MS.")
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

    def checkZeroSum(self, data1,data2):
        print("-- checking ZeroSum." )

        asum_data = numpy.abs(numpy.array(data1) + numpy.array(data2))
        check = numpy.all(asum_data < errLimit)
        self.assertTrue(check, msg='## Zero Sum check Failed ##\n{}'.format(asum_data)   )
        return True

######################
# check time
######################

    def checkTime(self, msName, row, refTime):
        print("-- checking Time --")

        # get time and inspection. 
        self. get_main(msName)
        # one output  
        T0 = self.tm[row]
        # check Time
        check = ( T0 == refTime )
        self.assertTrue(check, msg='## Time is not equal with expect.##\n{}'.format(T0) )

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

#+
# Generate DATa on FLOAT_DATA
#-
    def generate_data( self, msName ):
        print( "-- Generating MS." )
        self. get_main(defInputMs )
        # Test Slope
        offset = 0.0        # if specified non-zero, intensive fail can be cauesed.
        slope  = 0.0001
        # Time
        baseTime   = 0
        # Table Access
        with tbmanager(msName,nomodify=False) as tb:

            # create array (time, interval)
            NN = len(self.tm)
            arrayTime     = list(range(NN))
            arrayInterval = list(range(NN))
            for row in range(NN): 
                arrayTime    [row] = baseTime + (interval_0 * row)  
                arrayInterval[row] = interval_0

            # Put Time and Interval from the array  
            tb.putcol("TIME",       arrayTime  )
            tb.putcol("INTERVAL",   arrayInterval  )

            # Put DATA 
            # Float_data 
            arrayData = [list(range(1024)), list(range(1024))]
            for row in range(NN):
                for n in range(1024):
                    # values
                    x = row - numpy.floor(nRow/2)
                    val = offset + slope * x
                    arrayData[0][n] = val
                    arrayData[1][n] = val
                # write as an Array[2,1024] 
                tb.putcell("FLOAT_DATA", row,   arrayData  )
            #
          
        return True          

#================================
# sub function for TEST FIXTURE
#================================

    def check_averaged_result_N1(self, outMsName):
        '''
        this function inspect the Averaged result-MS
        '''
        # get the result and inspect #
        f_data = self.get_spectra(outMsName, 0 )     # use row=0 from RESULT 
        self.checkZero( f_data )    

        # Ref time 
        t_ref = (nRow -1)/2 * interval_0
        self.checkTime(outMsName, 0, t_ref)

        # Weight, Sigma 
        self.checkWeightSigma(outMsName, 0, nRow )

    def check_averaged_result_N3(self, outMsName):

        # get the result and inspect #
        f_data0 = self.get_spectra(outMsName, 0 )        # result on row=0
        f_data1 = self.get_spectra(outMsName, 1 )        # 
        f_data2 = self.get_spectra(outMsName, 2 )        # 

        # inspection # 
        self.checkZeroSum( f_data0, f_data2 )   # must be different sign)
        self.checkZero( f_data1 )               # must be zero)

        # Ref Time in 3 sections
        Tref_0 = interval_0 * (nRow /3 -1.0)/2
        Tref_1 = Tref_0 + interval_0 * (nRow/3)
        Tref_2 = Tref_1 + interval_0 * (nRow/3)

        # check Time
        self.checkTime(outMsName, 0, Tref_0)
        self.checkTime(outMsName, 1, Tref_1)
        self.checkTime(outMsName, 2, Tref_2)

        # Weight, Sigma 
        self.checkWeightSigma(outMsName, 0, (nRow/3) )
        self.checkWeightSigma(outMsName, 1, (nRow/3) )
        self.checkWeightSigma(outMsName, 2, (nRow/3) )


    def setOutfile_Timebin(self, testNo, numRec ):
        
        strTimebin = '{}s'.format(numRec * interval_0)
        outFile    = 'bave-{}-{}.ms'.format(testNo, numRec)
        return outFile, strTimebin;

#=================================================
# TEST FIXTURE
#==================================================

## ANTENNA ###
    def test_param1(self): 
        '''sdtimeaverage:: antenna = 'GBT' ''' 

        prm =  {'antenna' : 'GBT'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm), msg="Error in test_param1")

    def test_param1E(self):
        '''sdtimeaverage:: antenna = 'gBT' (Error) '''

        prm =  {'antenna'    : 'gBT' }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # must be false

## TIMEBIN ,make Average ###
    def test_param20(self):
        '''sdtimeaverage::20:: timebin=1282(N=3)  '''
        # set timebin string and private outputMS name.
        privateOutfile, timebin_str  = self.setOutfile_Timebin( 20, 1282 )
 
        prm =  {'timebin' : timebin_str,          
                'infile'  : defWorkMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # Check Result (zerosum check)
        self.check_averaged_result_N3(privateOutfile)

    def test_param21(self):
        '''sdtimeaverage::21: timebin=3846(N=1), timebin=''  '''

        # set timebin string and private outputMS name.
        privateOutfile, timebin_str  = self.setOutfile_Timebin( 21, 3846 )

        prm =  {'timebin' : timebin_str,             # Immediate Value ,
                'infile'  : defWorkMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )
     
        # Check Result (zerosum check)
        self.check_averaged_result_N1(privateOutfile) 

    def test_param22(self):
        '''sdtimeaverage::22: timebin=3846(N=1), timebin='all'  '''

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
    def test_param30(self): 
        '''sdtimeagerage::30:: timebin='all' ''' 

        prm =  {'timebin' : 'all'   }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
 
    def test_param31(self):
        '''sdtimeagerage::30:: timebin='ALL' '''

        # Run Task
        prm =  {'timebin' : 'ALL'   }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param32(self):
        '''sdtimeagerage::30:: timebin=''    '''

        # Run Task
        prm =  {'timebin' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param33E(self):
        '''sdtimeagerage::30:: timebin='Alles' (ERROR)    '''

        # Run Task
        prm =  {'timebin' : 'Alles'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # Error expected #  

## SCAN ###
    def test_param40(self):
        '''sdtimeagerage::40:: scan=2 '''

        # Run Task
        prm =  {'timebin' : '', 
                'scan'    : '2'   }
        self.run_task( prm )
        # check scan 
        self.check_scan (defOutputMs, 2)

    def test_param41a(self):
        '''sdtimeagerage::41a:: scan=61  '''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '61'   }i # Normal. In range. #
        self.assertTrue(self.run_task( prm ))
        # check scan 
        self.check_scan (defOutputMs, 61)

    def test_param41b(self):
        '''sdtimeagerage::41b:: scan=62  '''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '62'   } # ERROR : out of range in MS #
        self.assertFalse(self.run_task( prm ))
        # check scan 
        #    no check ...

    def test_param42(self):
        '''sdtimeagerage::42:: scan='' N=1  '''

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

## DATACOLUMN ###
    def test_param50(self):
        '''sdtimeaverage:: datacolumn = 'float_data' '''

        prm =  {'datacolumn' : 'float_data'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param51E(self):
        '''sdtimeaverage:: datacolumn = 'data' (Error) '''

        prm =  {'datacolumn' : 'data' }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # must be false

    def test_param52E(self):
        '''sdtimeaverage:: datacolumn = '' default=data is applied, makes Error. '''

        prm =  {'datacolumn' : '' }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # must be false

## FIELD ###
    def test_param60(self):
        '''sdtimeaverage:: field = 'FLS3a*' (EXACT NAME)'''

        prm =  {'field' : 'FLS3a*'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

    def test_param61E(self):
        '''sdtimeaverage:: field = 'hoge*' (NG NAME)'''

        prm =  {'field' : 'hoge'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm ))

    def test_param62(self):
        '''sdtimeaverage:: field = '' (OK :use default)'''

        prm =  {'field' : '*'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))

#### Control ######

def suite():
    return [test_sdtimeaverage]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

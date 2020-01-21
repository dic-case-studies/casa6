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
# local import
import re       # Regular Expr.
import time	# Measure used time.
# MS name for this test
defInputMs  = "sdimaging.ms"     # template MS
defWorkMs   = "sdimaging-t.ms"   # testing MS (modified here)
defOutputMs = "sdave.ms"         # (internal) output MS

# Test Conditio , Numerical error limit.
nRow     = 3843  ## DO NOT CHANGE ## 
errLimit  = 2.0e-08   # numerical error Limit of ZeroSum
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
# check Output
######################

    def checkOutputRec(self, msName, refNRow):
        print("-- checking Output Record Count --")

        # get time and inspection. 
        self. get_main(msName)
        # one output  
        nrow = len(self.tm)
        # check
        check = ( nrow == refNRow )
        self.assertTrue(check, msg='## Unexpected Row Count in Output.##\n{}'.format(nrow) )

######################
# check scan
######################

    def check_scan (self, outMsName, refValue ):
        print("-- checking scan selection --")
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

    # Check Wait and Sigma
    def checkWeightSigma(self, msName, row, weight_ref ):
        print( "-- checking Weight and Sigma --")

        self.get_spectra(msName, row )

        print( "Weight Ref :{0}".format(weight_ref) )
        print( "Weight     :{0}".format(self.wgt)   )
        print( "Sigma      :{0}".format(self.sgm)  )

        # Check (based on fomula about Sigma and Waight) #
        check1 =  (self.wgt[0] == weight_ref)
        check2 =  (self.wgt[1] == weight_ref)
        check3 =  ( (1.0/self.wgt[0])  - (self.sgm[0] * self.sgm[0])  < errLimit2 )
        check4 =  ( (1.0/self.wgt[1])  - (self.sgm[1] * self.sgm[1])  < errLimit2 )

        self.assertTrue(check1, msg='## Weight[0] is unexpected.##\n {}'.format(self.wgt[0]) )
        self.assertTrue(check2, msg='## Weight[1] is unexpected.##\n {}'.format(self.wgt[1]) )
        self.assertTrue(check3, msg='## Sigma [0] is unexpected.##\n {}'.format(self.sgm[0]) )
        self.assertTrue(check4, msg='## Sigma [1] is unexpected.##\n {}'.format(self.sgm[1]) )

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
    def get_spectra(self,msName, row ):
        with tbmanager(msName) as tb:
            # Spectra 
            self.data = tb.getcell('FLOAT_DATA',row)
            self.wgt  = tb.getcell('WEIGHT', row)
            self.sgm  = tb.getcell('SIGMA', row)           
        
        return self.data 

################################
# Generate Data on FLOAT_DATA
################################
    def generate_data( self, msName ):
        print( "----- Generating MS." )
        self. get_main(defInputMs )
        # Test Slope
        offset = 0.0        # if specified non-zero, an Intensive Fail will be activated.
        slope  = 0.0001     # (tunable) 
        # Time ,any  Date(Julian Day） by seconds  can be set here.
        baseTime   = 0.0
        # Table Access (with numPy array operation)
        with tbmanager(msName,nomodify=False) as tb:
            # get Rec size
            NN = len(self.tm)

            # get array shape of Spectra data, by getcolshapestring(), returning string like:"[2,1024]" via 'list'. 
            tk = re.split(",|\[|\]",tb.getcolshapestring('FLOAT_DATA', nrow=1)[0] ) # delimter = [ , ]
            nChan = int(tk[2]) #             separating to :: <zero> [<1st> ,<2nd*> ]<3rd>
            nPol  = int(tk[1]) # (not used ) separating to :: <zero> [<1st*> ,<2nd> ]<3rd> 

            # create array (time, interval)
            arrayTime = testInterval * numpy.arange(0,NN,dtype=numpy.float64) + baseTime
            arrayInterval =  numpy.full(NN,testInterval, dtype=numpy.float64)
        
            # put to column
            print( "----- Putting Time,INTERVAL." )
            tb.putcol("TIME",       arrayTime  )
            tb.putcol("INTERVAL",   arrayInterval  )

            # create Test-Data [use numpy.array]
            print( "----- Calculating Curve." )
            NN1 = (NN-1)/2
            L = numpy.linspace(-NN1,NN1, NN)* slope + offset 
            VAL = numpy.tile(L, [nChan, 1])
            arrayData3 = numpy.array( [VAL,VAL] )

            # write to the column at once
            print( "----- Putting Curve."  )
            tb.putcol("FLOAT_DATA",   arrayData3  )
        print( "----- Done."  )  
        return True          

##################################
# sub function for TEST FIXTURE
##################################

    def check_averaged_result_N1(self, outMsName):
        '''
        This function inspects the Averaged result-MS
        All the spectral data is averaged. one averaged result remains, which must be Zero.
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
        '''
        This function inspects the Averaged result-MS by 3 averaged results.
        Three sections are Averaged. 1st.result and 3rd.result are different sign and sum =0
        The 2nd. section makes Zero Sum.
        Note: In Test-MS test data is designed by odd functional curve. 
        '''
        # get the result  #
        fData0 = self.get_spectra(outMsName, 0 )        # result on row=0
        fData1 = self.get_spectra(outMsName, 1 )        # row=1
        fData2 = self.get_spectra(outMsName, 2 )        # row=2

        # inspection # 
        self.checkZeroSum( fData0, fData2 )   # must be different sign = Zero Sum.
        self.checkZero( fData1 )              # must be zero

        # Ref Time in 3 sections (becomes centre of the section)
        Tref0 = testInterval * (nRow /3 -1.0)/2
        Tref1 = Tref0 + testInterval * (nRow/3)
        Tref2 = Tref1 + testInterval * (nRow/3)

        # check Time
        self.checkTime(outMsName, 0, Tref0)
        self.checkTime(outMsName, 1, Tref1)
        self.checkTime(outMsName, 2, Tref2)

        # check Weight, Sigma 
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
        self.checkOutputRec(privateOutfile, 1 )

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
        self.checkOutputRec(privateOutfile, 1 )

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
        self.checkOutputRec(privateOutfile, 1 )

## SCAN ###
    def test_para10(self):
        '''sdtimeagerage::10:: scan=2 (Within the range)'''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '2'   }
        self.run_task( prm )
        # check scan 
        self.check_scan (defOutputMs, 2)
        self.checkOutputRec(defOutputMs, 1 )

    def test_param11a(self):
        '''sdtimeagerage::11a:: scan=61 (Within the range)'''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '61'   } # Normal. In range. #
        self.assertTrue(self.run_task( prm ))
        # check scan 
        self.check_scan (defOutputMs, 61)
        self.checkOutputRec(defOutputMs, 1 )

    def test_param11b(self):
        '''sdtimeagerage::11b:: scan=62 (Error Out of range) '''

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '62'   } # ERROR : out of range in MS #
        self.assertFalse(self.run_task( prm ))
        # check scan 
        #    no check ...

    def test_param12(self):
        '''sdtimeagerage::12:: scan='' (no number) Default action. '''

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
        self.checkOutputRec(privateOutfile, 1 )

## FIELD ###
    def test_param20(self):
        '''sdtimeaverage::20:: field = 'FLS3a*' (Exact NAME)'''

        prm =  {'field' : 'FLS3a*'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

    def test_param21E(self):
        '''sdtimeaverage::21E:: field = 'hoge*' (Error :Bad NAME)'''

        prm =  {'field' : 'hoge'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm ))

    def test_param22(self):
        '''sdtimeaverage::22:: field = '' (OK :use default)'''

        prm =  {'field' : '*'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

## SPW ###
    def test_param30(self):
        '''sdtimeaverage::30:: spw = '1' (exist)'''

        prm =  {'spw' : '0'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

    def test_param31E(self):
        '''sdtimeaverage::31E:: spw = '9' (Error: Not Exist)'''

        prm =  {'spw' : '9'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm ))

    def test_param32(self):
        '''sdtimeaverage:: spw = '' (default)'''

        prm =  {'spw' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

    def test_param33(self):
        '''sdtimeaverage::33:: spw = '*' (OK: Wildcard)'''

        prm =  {'spw' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

## ANTENNA ###
    def test_param40(self): 
        '''sdtimeaverage::40:: antenna = 'GBT' (Exact name)''' 

        prm =  {'antenna' : 'GBT'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm), msg="Error in test_param1")
        self.checkOutputRec(defOutputMs, 1 )

    def test_param41E(self):
        '''sdtimeaverage::41E antenna = 'gBT' (Error: Bad name) '''

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
        self.checkOutputRec(privateOutfile, 3 )

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
        self.checkOutputRec(privateOutfile, 1 )

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
        self.checkOutputRec(privateOutfile, 1 )

## TIMEBIN ###
    def test_param110(self): 
        '''sdtimeagerage::110:: timebin='all' ''' 

        prm =  {'timebin' : 'all'   }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )
 
    def test_param111(self):
        '''sdtimeagerage::111:: timebin='ALL' '''

        # Run Task
        prm =  {'timebin' : 'ALL'   }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

    def test_param112(self):
        '''sdtimeagerage::112:: timebin='' (default)    '''

        # Run Task
        prm =  {'timebin' : ''  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

    def test_param113E(self):
        '''sdtimeagerage::113E:: timebin='Alles' (ERROR: Bad keyword)    '''

        # Run Task
        prm =  {'timebin' : 'Alles'  }
        # Run Task and check
        self.assertFalse(self.run_task( prm )) # Error expected #  

    def test_param114(self):
        '''sdtimeagerage::114:: timebin='aLL' (OK: Upper/Lower case mixed)    '''

        # Run Task
        prm =  {'timebin' : 'aLL'  }
        # Run Task and check
        self.assertTrue(self.run_task( prm ))
        self.checkOutputRec(defOutputMs, 1 )

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

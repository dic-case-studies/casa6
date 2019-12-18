import shutil
import unittest
import os
import numpy
import math
import sys
import exceptions
import filecmp
import glob
# from tasks import nrobeamaverage
from tasks import sdtimeaverage
from taskinit import mstool, tbtool
from __main__ import default
import testhelper as th
from sdutil import tbmanager, toolmanager, table_selector

# Define the root for the data files
datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/sdimaging/"

# MS name for this test
def_inputMs  = "sdimaging.ms"
def_workMs   = "sdimaging-t.ms"
def_outputMs = "bave.ms"

# Compare err limit , ideally vector(1024 x 2) is the best
nRow     = 3843  ## DO NOT CHANGE ## 
errLimit  = 5e-6
errLimit2 = 1e-08
zeroData = [[errLimit]*1024, [errLimit]*1024]
interval_0 = 1.0

##############
# Test Entry
#############
class test_sdtimeaverage(unittest.TestCase):
    def setUp(self):
        default(sdtimeaverage)

        # copy from master
        self.i_ms = def_inputMs
        os.system('cp -RL '+ os.path.join(datapath,self.i_ms) +' '+ self.i_ms)
        os.system('cp -RL '+ os.path.join(datapath,self.i_ms) +' '+ def_workMs)

        # output MS 
        self.o_ms = ""

        # params
        self.interval = interval_0
        self.tol = errLimit

        # default Args (minimum)
        self.args = {'infile'     :  def_inputMs,
                     'outfile'    :  def_outputMs,
                     'datacolumn' :  'float_data'    # CASR-474 (float ->data) 
                    }

        self. get_main( def_inputMs )
        self. generate_data( def_workMs )

    def tearDown(self):

        # delete copied in-MS and out-MS
        print( "tearDown::deleting MSs")

        os.system('rm -rf ' + self.i_ms )
        os.system('rm -rf ' + def_outputMs )   ## Comment out , for DEBUG ##

        os.system('rm -rf ' + def_workMs )
        os.system('rm -rf ' + "bave*.ms" )
        return

##############
# Run Task
##############
    def run_task(self, aux_args=None):
        print( "run_task::stars" )

        if aux_args is not None:
            for k in aux_args: self.args[k] = aux_args[k]
        sdtimeaverage(**self.args)

#################
# Check Result
#################

    def checkZero(self,data):
        print("-- checking Zero --")

        check = numpy.abs(data) < zeroData
        if check.all()==False:
            print ( "## Zero check Failed ##" )
            print ( data )
            return False
        print("-- checking Zero OK.--")
        return  True

    def checkZeroSum(self, data1,data2):
        print("-- checking ZeroSum." )

        asum_data = numpy.abs(numpy.array(data1) + numpy.array(data2))
        print(asum_data)

        check = asum_data < zeroData
        if check.all()==False:
            print ( "## Zero Sum check Failed ##" )
            return False
        print("-- checking ZeroSum. OK" )
        return True

######################
# check time
######################

    def checkTime(self, msname, row, refTime):
        print("-- checking Time --")

        # get time and inspection. 
        self. get_main(msname)
        # one output  
        T0 = self.tm[row]
        # check Time
        check = ( T0 == refTime )
        if check:
            print("-- checking Time OK --")

        self.assertTrue( check )

#########################
# Generating Test Data 
#########################

#+
# Read Data from Specified MS
#-
    # MAIN #
    def get_main(self, MsName):
        # get MAIN table data
        with tbmanager(MsName) as tb:

            # Key data 
            self.tm = tb.getcol('TIME')
            self.a1 = tb.getcol('ANTENNA1')
            self.a2 = tb.getcol('ANTENNA2')
            self.dd = tb.getcol('DATA_DESC_ID')
            self.sc = tb.getcol('SCAN_NUMBER')
            self.st = tb.getcol('STATE_ID')

    # DATA (spectra) #
    def get_spectra(self,MsName, row ):
        with tbmanager(MsName) as tb:
            # Spectra 
            self.data = tb.getcell('FLOAT_DATA',row)
            self.wgt  = tb.getcell('WEIGHT', row)
            self.sgm  = tb.getcell('SIGMA', row)

        return self.data

    # Chck Wait and Sigma
    def checkWeightSigma(self, msname, row, weight_ref ):
        print( "-- checking Weight and Sigma --")

        self.get_spectra(msname, row )

        print "Weight Ref", weight_ref
        print "Weight ",self.wgt 
        print "Gigma  ",self.sgm 

        # check #
        check =  (self.wgt[0] == weight_ref) and \
                 (self.wgt[1] == weight_ref) and \
                 ( (1.0/self.wgt[0])  - (self.sgm[0] * self.sgm[0])  < errLimit2 ) and \
                 ( (1.0/self.wgt[1])  - (self.sgm[1] * self.sgm[1])  < errLimit2 )
        if check:
            print( "-- checking Weight and Sigma OK --")

        self.assertTrue(check)
#+
# Generate DATa on FLOAT_DATA
#-
    def generate_data( self, MsName ):

        self. get_main( def_inputMs )

        # Test Slope
        offset = 0.0     # if specified non-zero, intensive fail can be cauesed.
        slope  = 0.0001
        # Time
        baseTime   = 0
        # Table Access
        with tbmanager(MsName,nomodify=False) as tb:
            # write to cell in each row #
            for row in range(len(self.tm) ):
                data_array = [list(range(1024)), list(range(1024))]
                # make const #i
                N = len(self.tm)
                for n in range(1024):
                    x = row - numpy.floor(N/2)
                    data_array[0][n] =  offset + slope * x
                    data_array[1][n] =  offset + slope * x
                # write as an Array[2,1024] 
                tb.putcell("FLOAT_DATA", row,  data_array  )
                # Time and Interval  
                tb.putcell("TIME",       row,  baseTime + (interval_0 * row)  )
                tb.putcell("INTERVAL",   row,  interval_0  )
            #endfor
        return


#=================================================
# TEST FIXTURE
#==================================================

# Generating TestMS only 
#    def test_param0(self):
#        print( "test_param0:: generating Test MS. ")
#        # test MS generation
#        self. get_main( def_inputMs )
#        self. generate_data( def_workMs )

# 'all' + Antenna Name
    def test_param1(self):
        print( "XXXXXXXX test_param(1: timebin=all, antenna=GBT ) XXXXXXXX")

        prm =  {'infile'  : def_inputMs,
                'timebin' : 'all',
                'antenna' : 'GBT'  }
        # Run Task
        self.run_task( prm )


# N=3 averaged Out.
    def test_param20(self):
        print( "XXXXXXXX test_param(20: timebin=1282) N=3 XXXXXXXX")

        timebin_str = str(1282 * interval_0)+'s'
        privateOutfile = 'bave-20-1282.ms'
        prm =  {'timebin' : timebin_str,
                'infile'  : def_workMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # get the result and inspect #
        f_data0 = self.get_spectra(privateOutfile, 0 )        # result on row=0
        f_data1 = self.get_spectra(privateOutfile, 1 )        # 
        f_data2 = self.get_spectra(privateOutfile, 2 )        # 

        # inspection # 
        self.assertTrue(self.checkZeroSum( f_data0, f_data2 ))   # assert (must be different sign)
        self.assertTrue(self.checkZero( f_data1 ))               # assert (must be zero)

        # Ref Time in 3 sections
        Tref_0 = interval_0 * (nRow /3 -1.0)/2
        Tref_1 = Tref_0 + interval_0 * (nRow/3)
        Tref_2 = Tref_1 + interval_0 * (nRow/3)

        # check Time
        self.checkTime(privateOutfile, 0, Tref_0)
        self.checkTime(privateOutfile, 1, Tref_1)
        self.checkTime(privateOutfile, 2, Tref_2)

        # Weight, Sigma 
        self.checkWeightSigma(privateOutfile, 0, (nRow/3) )
        self.checkWeightSigma(privateOutfile, 1, (nRow/3) )
        self.checkWeightSigma(privateOutfile, 2, (nRow/3) )

# N=1 (timebin = actual time)
    def test_param21(self):
        print( "XXXXXXXX test_param(21: timebin=3846) N=1 XXXXXXXX")
        timebin_str = str(3846 * interval_0)+'s'
        privateOutfile = 'bave-21-3846.ms'
        prm =  {'timebin' : timebin_str,             # Immediate Value ,
                'infile'  : def_workMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # get the result and inspect #
        f_data = self.get_spectra(privateOutfile, 0 )   # use row=0 from RESULT
        self.assertTrue(self.checkZero( f_data ))       # assert

        # Ref time 
        Tref = (nRow -1)/2 * interval_0
        # check Time
        self.checkTime(privateOutfile, 0, Tref)

        # Weight, Sigma 
        self.checkWeightSigma(privateOutfile, 0, nRow )

# N=1 (timebin=all) 
    def test_param22(self):
        print( "XXXXXXXX test_param(22: timebin=(all)) N=1 XXXXXXXX")

        privateOutfile = 'bave-22-3846.ms'
        prm =  {'timebin' : 'all',                # default = all is applied.
                'infile'  : def_workMs,
                'outfile' : privateOutfile  }
        # Run Task
        self.run_task( prm )

        # get the result and inspect #
        f_data = self.get_spectra(privateOutfile, 0 )   # use row=0 from RESULT
        self.assertTrue(self.checkZero( f_data ))          # assert

        # Ref time 
        Tref = (nRow -1)/2 * interval_0
        # check Time
        self.checkTime(privateOutfile, 0, Tref)

        # Weight, Sigma 
        self.checkWeightSigma(privateOutfile, 0, nRow )

    def test_param30(self):
        print( "XXXXXXXX test_param(30: timebin=all) XXXXXXXX")

        prm =  {'timebin' : 'all'   }
        self.run_task( prm )

    def test_param31(self):
        print( "XXXXXXXX test_param(31: timebin = ALL (Capital) XXXXXXXX")

        # Run Task
        prm =  {'timebin' : 'ALL'   }
        self.run_task( prm )

    def test_param32(self):
        print( "XXXXXXXX test_param(32: timebin='' ) XXXXXXXX")

        # Run Task
        prm =  {'timebin' : ''  }
        self.run_task( prm )

    def test_param40(self):
        print( "XXXXXXXX test_param(40: timebin='', scan=1 ) XXXXXXXX")

        # Run Task
        prm =  {'timebin' : '',
                'scan'    : '2'   }
        self.run_task( prm )

        # get table and inspection. 
        self. get_main(def_outputMs)
        # one output
        scan = self.sc[0]
        # check scan ID
        self.assertTrue(len(self.sc)==1 )
        self.assertTrue (scan == 2 )
#+
# The Last Fixture
#  - delete temporary file(s).
#-
#    def test_param99(self):
#        print( "XXXXXXXX test_param99:: deleting Test MS. XXXXXXXX")
#        if True:
#            os.system('rm -rf ' + def_workMs )
#            os.system('rm -rf ' + "bave*.ms" )

# CASA5 featre 
def suite():
    return [test_sdtimeaverage]

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
def_inputMs  = "sdimaging.ms"
def_workMs   = "sdimaging-t.ms"
def_outputMs = "bave.ms"

# Compare err limit , ideally vector(1024 x 2) is the best
nRow     = 3843
errLimit = 5e-8
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

    def tearDown(self):

        # delete copied in-MS and out-MS
        print( "tearDown::deleting MSs")

        os.system('rm -rf ' + self.i_ms )
        os.system('rm -rf ' + def_outputMs )   ## Comment out , for DEBUG ##

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

    def checkZero(self, data):
        print ("-- checking Zero.")
        for n in range(len(data)):
            if (data[0][n] >  errLimit) or (data[1][n] >  errLimit):
                print ( "## zero check Failed ##" )
                print ( data[0][n] )
                print ( data[1][n] )
                return False
        return True

    def checkZeroSum(self, data1,data2):
        print("-- checking Equal." )
        for n in range(len(data1)):
            if abs(data1[0][n] + data2[0][n]) >  errLimit:
                return False
            if abs(data1[1][n] + data2[1][n]) >  errLimit:
                return False
        return True


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
            self.sig  = tb.getcell('SIGMA', row)           
        
        return self.data 

#+
# Generate DATa on FLOAT_DATA
#-
    def generate_data( self, MsName ):

        self. get_main( def_inputMs )

        # Test Slope
        offset = 0.0
        slope  = 0.0001
        # Time
        baseTime   = 0
        # Table Access
        with tbmanager(MsName,nomodify=False) as tb:
            # write to cells #
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
                tb.putcell("TIME",       row,  baseTime + row  )
                tb.putcell("INTERVAL",   row,  interval_0  )
            #endfor
        return          
          
#=================================================
# TEST FIXTURE
#==================================================

# TestMS 
    def test_param0(self):
        print( "test_param0:: generating Test MS. ") 
        # test MS generation
        self. get_main( def_inputMs )
        self. generate_data( def_workMs )
# 'all'    
    def test_param1(self): 
        print( "XXXXXXXX test_param(1: antenna=GBT ) XXXXXXXX")

        prm =  {'infile'  : def_inputMs,
                'timebin' : 'all', 
                'antenna':'GBT'  }
        self.run_task( prm )


# N=3
    def test_param20(self):
        print( "XXXXXXXX test_param(20: timebin=1282) N=3 XXXXXXXX")

        privateOutfile = 'bave-20-1282.ms'
        prm =  {'timebin' : '1282s',
                'infile'  : def_workMs,
                'outfile' : privateOutfile  }
        self.run_task( prm )

        # get the result and inspect #
        f_data0 = self.get_spectra(privateOutfile, 0 )        # result on row=0
        f_data1 = self.get_spectra(privateOutfile, 1 )        # 
        f_data2 = self.get_spectra(privateOutfile, 2 )        # 

        self.assertTrue(self.checkZero( f_data1 ))               # assert (must be zero)
        self.assertTrue(self.checkZeroSum( f_data0, f_data2 ))   # assert (must be different sign)
  
        # get time and inspection. 
        self. get_main(privateOutfile)
        Tref_0 = (nRow /3 -1.0)/2
        Tref_1 = Tref_0 + (nRow/3)
        Tref_2 = Tref_1 + (nRow/3)
        self.assertTrue( self.tm[0]== Tref_0 )
        self.assertTrue( self.tm[1]== Tref_1 )
        self.assertTrue( self.tm[2]== Tref_2 )
 
# N=1 (timebin = actual time)
    def test_param21(self):
        print( "XXXXXXXX test_param(21: timebin=3846) N=1 XXXXXXXX")

        privateOutfile = 'bave-21-3846.ms'
        prm =  {'timebin' : '3846s',
                'infile'  : def_workMs,
                'outfile' : privateOutfile  }
        self.run_task( prm )

        # get the result and inspect #
        f_data = self.get_spectra(privateOutfile, 0 )   # use row=0 from RESULT
        self.assertTrue(self.checkZero( f_data ))          # assert
 
        # get time and inspection. 
        self. get_main(privateOutfile) 
        T0 = self.tm[0]
        Tref = (nRow -1)/2 * interval_0
        self.assertTrue( T0 == Tref ) 

# N=1 (timebin=all) 
    def test_param22(self):
        print( "XXXXXXXX test_param(21: timebin=3846) N=1 XXXXXXXX")

        privateOutfile = 'bave-22-3846.ms'
        prm =  {'timebin' : 'all',
                'infile'  : def_workMs,
                'outfile' : privateOutfile  }
        self.run_task( prm )

        # get the result and inspect #
        f_data = self.get_spectra(privateOutfile, 0 )   # use row=0 from RESULT
        self.assertTrue(self.checkZero( f_data ))          # assert
 
        # get time and inspection. 
        self. get_main(privateOutfile) 
        T0 = self.tm[0]
        Tref = (nRow -1)/2 * interval_0
        self.assertTrue( T0 == Tref )

    def test_param30(self): 
        print( "XXXXXXXX test_param(30: timebin=all) XXXXXXXX")

        prm =  {'timebin' : 'all'   }
        self.run_task( prm )
   
    def test_param31(self):
        print( "XXXXXXXX test_param(31: timebin = ALL (Capital) XXXXXXXX")

        prm =  {'timebin' : 'ALL'   }
        self.run_task( prm )
  
    def test_param32(self):
        print( "XXXXXXXX test_param(32: timebin='' ) XXXXXXXX")
        prm =  {'timebin' : ''  }
        self.run_task( prm )

 
    def test_param99(self):
        print( "test_param99:: deleting Test MS. ")
##        os.system('rm -rf ' + def_workMs )
##        os.system('rm -rf ' + "bave*.ms" )


#### Controled ######

def suite():
    return [test_sdtimeaverage]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

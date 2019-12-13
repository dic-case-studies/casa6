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

# Important parameter
interval = 2.99827  # same as the MAIN in sdimaging.ms
num_ave  = 3684
errLimit = 1e-05



def check_eq(val, expval, errLimit=None):
    """Checks that val matches expval within tol."""

    pass
    return


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
        self.interval = interval
        self.tol = errLimit

        # default Args (minimum)
        self.args = {'infile'     :  def_inputMs,
                     'outfile'    :  def_outputMs,
                     'datacolumn' :  'float_data'    # CASR-474 (float ->data) 
                    }
    def tearDown(self):

        # delete copied in-MS and out-MS
        print( "tearDown::deleting MSs")

#       os.system('rm -rf ' + self.i_ms )
#       os.system('rm -rf ' + def_workMs )
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


#+
# Get number of data (record count)
#-
    def check_num_data(self, num_ave=1):
        pass     
        return 

    def get_num_data(self, stcol=None):
        pass
        return

#+
# Check Values
#-
    def check_values(self, num_ave):

        # time

        # antenna ID

        # spectrum

        # weight and sigma

        pass
        return


    def _do_check_values(self, iidx, oidx):
        pass
        return

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
# Write Data (Scalar)
#-
    def put_scalar_data(self, MsName):
        with tbmanager(MsName,nomodify=False) as tb:
            # write value 
            val = 999.8888
            for irow in range(10):
                tb.putcell("TIME", irow,  val )

        return
#+
# Dump float_data
#-

    def dump_data( self, MsName):
       self. get_main( def_inputMs )

       # get (float)data from each row 
       for row in range(len(self.tm) ):
           print( "row=",row)
           data_tmp = self. get_spectra( def_inputMs, row )

           for pos in range(1024):
               print (  row, pos, data_tmp[0][pos], data_tmp[1][pos] )

    def generate_data( self, MsName ):

        self. get_main( def_inputMs )
        data_row0 = self. get_spectra( MsName, 0 )  # data (row 0)

        # fill const #

        for i in range(2):
            for j in range(1024):
                data_row0[i][j] = 0.1

        for row in range(len(self.tm) ):
            #print( "- generating rec=", row);
            with tbmanager(MsName,nomodify=False) as tb:
                # write as an Array[2,1024] 
                tb.putcell("FLOAT_DATA", row,  data_row0 ) 
        return   
          
#+
# Prepare variable on self context.
#-
    def declar(self):
        self.name = None     # test var.
        self.address = None  # test var.
 
#=================================================
# TEST FIXTURE
#==================================================

    def test_param1(self): 
        print( "XXXXXXXX test_param(1) XXXXXXXX")

        print ( "- reading MAIN table" )
        self. get_main( def_inputMs )

        print ( "- generationg DATA " )
        self. generate_data( def_workMs )

        prm =  {'infile'  : def_inputMs,
                'timebin' : 'all', 
                'antenna':'GBT'  }
        self.run_task( prm )


    def test_param20(self):
        print( "XXXXXXXX test_param(20: timebin=all) XXXXXXXX")

        prm =  {'timebin' : '6128s',
                'infile'  : def_workMs,
                'outfile' : 'bave-20-6128.ms'  }
        self.run_task( prm )

    def test_param21(self):
        print( "XXXXXXXX test_param(20: timebin=all) XXXXXXXX")

        prm =  {'timebin' : '12256s',
                'infile'  : def_workMs,
                'outfile' : 'bave-20-12256.ms'  }
        self.run_task( prm )



    '''

    def test_param30(self): 
        print( "XXXXXXXX test_param(30: timebin=all) XXXXXXXX")

        prm =  {'timebin' : 'all', 
                'antenna' : 'GBT',
                'outfile' : 'bave-30.ms'  }
        self.run_task( prm )
   
     
    def test_param31(self):
        print( "XXXXXXXX test_param(31: timebin Nospecified) XXXXXXXX")

        prm =  {'outfile' : 'bave-31.ms'  }
        self.run_task( prm )


    def test_param32(self):
        print( "XXXXXXXX test_param(3-2: timebin='' ) XXXXXXXX")
        prm =  {'timebin' : '',
                'outfile' : 'bave-32.ms'  }
        self.run_task( prm )

   '''

    #
    # ORIGINAL
    #

    '''
    def test_default(self): # no time averaging(timebin='0s'), rewriting beam IDs only
        print( "XXXXXXXX test_default XXXXXXXX")
        self.run_task()
        self.check_num_data()
        self.check_values()
   
    def test_time_averaging(self): # every two on-spectra are averaged into one specrum
        print( "XXXXXXXX test_time_averaging XXXXXXXX")
        self.run_task({'timebin': self.get_timebin(num_ave)})

        self.check_num_data(num_ave)
        self.check_values(num_ave=num_ave) # for the first data with state=on-source, spw=0

    def test_time_averaging2(self): # every two on-spectra are averaged into one specrum
        print( "XXXXXXXX test_time_averaging2 XXXXXXXX")
        self.run_task({'timebin': 'all'} )

        self.check_num_data(num_ave)
        self.check_values(num_ave=num_ave) # for the first data with state=on-source, spw=0

    def test_time_averaging3(self): # every two on-spectra are averaged into one specrum
        print( "XXXXXXXX test_time_averaging3 XXXXXXXX")
        self.run_task({'timebin': '100000s'} )

        self.check_num_data(num_ave)
        self.check_values(num_ave=num_ave) # for the first data with state=on-source, spw=0
    '''

#### Controled ######

def suite():
    return [test_sdtimeaverage]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

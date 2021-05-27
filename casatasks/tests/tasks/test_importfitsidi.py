#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for the FITS-IDI import to MS                   #
#    
#                                                                           #
# Rationale for Inclusion:                                                  #
#    The conversion of FITS-IDI to MS needs to be verified.                 #
#                                                                           # 
# Features tested:                                                          #
#    1) Is the import performed without raising exceptions                  #
#    2) Do all expected tables exist                                        #
#    3) Can the MS be opened                                                #
#    4) Do the tables contain expected values                               #
#    5) Can several FITS-IDI files be read into one MS                      #
#                                                                           #
# Input data:                                                               #
#                                                                           #
#############################################################################
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import numpy
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, ms, table
    from casatasks import importfitsidi

    _ms = ms( )
    _tb = table( )

    datapath = ctsys.resolve('unittest/importfitsidi/')
else:
    from __main__ import default
    from tasks import *
    from taskinit import *

    # not local tools
    _ms = ms
    _tb = tb

    datapath=os.environ.get('CASAPATH').split()[0]+'/casatestdata/unittest/importfitsidi/'

myname = 'importfitsidi-unit-test'

# default dataset name
my_dataset_names = ['n09q2_1_1-shortened.IDI1',
                    'n09q2_1_1-shortened-part1.IDI1',
                    'n09q2_1_1-shortened-part2.IDI1',
                    'emerlin_multiuv.IDI1',
                    '1331_3030_C-Band_5GHz__64.000_128.fits',
                    'VLBA_TL015A_tl015arecor_BIN0_SRC0_1_201020T164655.idifits']

# name of the resulting MS
msname = my_dataset_names[0]+'.ms'

def checktable(thename, theexpectation):
    global msname, myname
    _tb.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print(myname, ": comparing ", mycell)
        value = _tb.getcell(mycell[0], mycell[1])
        # see if value is array
        try:
            isarray = value.__len__
        except:
            # it's not an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement = (value == mycell[2])
            else:
                in_agreement = ( abs(value - mycell[2]) < mycell[3]) 
        else:
            if is_CASA6:
                stype = str
            else:
                stype = basestring
            if isinstance(value, stype):
                in_agreement = value == mycell[2]
            else:
                # it's an array
                # zero tolerance?
                if mycell[3] == 0:
                    in_agreement =  (value == mycell[2]).all() 
                else:
                    try:
                        in_agreement = (abs(value - mycell[2]) < mycell[3]).all()
                    except:
                        in_agreement = False
        if not in_agreement:
            print(myname, ":  Error in MS subtable", thename, ":")
            print("     column ", mycell[0], " row ", mycell[1], " contains ", value)
            print("     expected value is ", mycell[2])
            _tb.close()
            return False
    _tb.close()
    print(myname, ": table ", thename, " as expected.")
    return True


###########################
# beginning of actual test 

class test_importfitsidi(unittest.TestCase):
    
    def setUp(self):
        res = None

        for fname in my_dataset_names:
            if(os.path.exists(fname)):
                os.remove(fname)
            datasetPath = os.path.join(datapath,fname)
            if is_CASA6:
                datasetPath = ctsys.resolve(datasetPath)
            shutil.copy(datasetPath, fname)

        if not is_CASA6:
            default(importfitsidi)
        
    def tearDown(self):
        for fname in my_dataset_names:
            os.remove(fname)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        
    def test1(self):
        '''fitsidi-import: Test good input'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importfitsidi(my_dataset_names[0], msname,  scanreindexgap_s=100., constobsid=True)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            _ms.close()
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 2.0, 1E-8],
                         ['DATA',      42,[[ 1.06945515 +7.91624188e-09j,
                                             0.98315531 +9.31322575e-10j,
                                             1.05244470 +5.77396486e-09j,
                                             0.90496856 -0.00000000e+00j,
                                             0.93005872 -6.71682887e-09j,
                                             0.80769897 -0.00000000e+00j,
                                             0.93059886 -6.97121116e-09j,
                                             0.77081358 -0.00000000e+00j,
                                             0.93020338 -7.45058060e-09j,
                                             0.83353537 -0.00000000e+00j,
                                             0.91982168 -5.54113422e-09j,
                                             0.88411278 -4.65661287e-10j,
                                             1.02857709 +5.78550630e-09j,
                                             0.93398595 -0.00000000e+00j,
                                             1.13884020 +1.01289768e-08j,
                                             2.49237108 -0.00000000e+00j ]], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
            expected = [
                         ['UVW',       638, [171529.37575288, -786712.70341456, 210321.20978818], 1E-8],
                         ['EXPOSURE',  638,  2.0, 1E-8],
                         ['DATA',      638, [[-0.00224198+0.00067056j,
                                              -0.00475123+0.0024323j,
                                              -0.00416393+0.00212671j,
                                              -0.00565350+0.00340364j,
                                              -0.00527357+0.00011977j,
                                              -0.00292699+0.00131954j,
                                              -0.00429945+0.00035823j,
                                              -0.00545671-0.00033945j,
                                              -0.00646004+0.00037293j,
                                              -0.00419376-0.00115011j,
                                              -0.00508117+0.00045939j,
                                              -0.00501660-0.00047975j,
                                              -0.00444734-0.00101535j,
                                              -0.00384988-0.00102731j,
                                              -0.00551326+0.00101364j,
                                              -0.00337701+0.00080481j]], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [ 3370605.8469,  711917.6732,  5349830.8438], 0.0001],
                         ['DISH_DIAMETER',1, 0.0, 0.0] # the EVN default value
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "SPECTRAL_WINDOW"
            expected = [ ['NUM_CHAN',        7, 16, 0],
                         ['TOTAL_BANDWIDTH', 7, 8E6, 0],
                         ['CHAN_WIDTH',      7, [ 500000.,  500000.,  500000.,  500000.,  500000.,  500000.,
                                                  500000.,  500000.,  500000.,  500000.,  500000.,  500000.,
                                                  500000.,  500000.,  500000.,  500000.], 1E-8],
                         ['CHAN_FREQ',       7, [  4.32184900e+10,   4.32189900e+10,   4.32194900e+10,   4.32199900e+10,
                                                   4.32204900e+10,   4.32209900e+10,   4.32214900e+10,   4.32219900e+10,
                                                   4.32224900e+10,   4.32229900e+10,   4.32234900e+10,   4.32239900e+10,
                                                   4.32244900e+10,   4.32249900e+10,   4.32254900e+10,   4.32259900e+10], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test2(self):
        '''fitsidi-import: Test good input, list of two input files'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importfitsidi([my_dataset_names[1],my_dataset_names[2]], msname)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            _ms.close()
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 2.0, 1E-8],
                         ['DATA',      42,[[ 1.06945515 +7.91624188e-09j,
                                             0.98315531 +9.31322575e-10j,
                                             1.05244470 +5.77396486e-09j,
                                             0.90496856 -0.00000000e+00j,
                                             0.93005872 -6.71682887e-09j,
                                             0.80769897 -0.00000000e+00j,
                                             0.93059886 -6.97121116e-09j,
                                             0.77081358 -0.00000000e+00j,
                                             0.93020338 -7.45058060e-09j,
                                             0.83353537 -0.00000000e+00j,
                                             0.91982168 -5.54113422e-09j,
                                             0.88411278 -4.65661287e-10j,
                                             1.02857709 +5.78550630e-09j,
                                             0.93398595 -0.00000000e+00j,
                                             1.13884020 +1.01289768e-08j,
                                             2.49237108 -0.00000000e+00j ]], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
            expected = [
                         ['UVW',       638, [171529.37575288, -786712.70341456, 210321.20978818], 1E-8],
                         ['EXPOSURE',  638,  2.0, 1E-8],
                         ['DATA',      638, [[-0.00224198+0.00067056j,
                                              -0.00475123+0.0024323j,
                                              -0.00416393+0.00212671j,
                                              -0.00565350+0.00340364j,
                                              -0.00527357+0.00011977j,
                                              -0.00292699+0.00131954j,
                                              -0.00429945+0.00035823j,
                                              -0.00545671-0.00033945j,
                                              -0.00646004+0.00037293j,
                                              -0.00419376-0.00115011j,
                                              -0.00508117+0.00045939j,
                                              -0.00501660-0.00047975j,
                                              -0.00444734-0.00101535j,
                                              -0.00384988-0.00102731j,
                                              -0.00551326+0.00101364j,
                                              -0.00337701+0.00080481j]], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [ 3370605.8469,  711917.6732,  5349830.8438], 0.0001],
                         ['DISH_DIAMETER',1, 0.0, 0.0] # the EVN default value
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "SPECTRAL_WINDOW"
            expected = [ ['NUM_CHAN',        7, 16, 0],
                         ['TOTAL_BANDWIDTH', 7, 8E6, 0],
                         ['CHAN_WIDTH',      7, [ 500000.,  500000.,  500000.,  500000.,  500000.,  500000.,
                                                  500000.,  500000.,  500000.,  500000.,  500000.,  500000.,
                                                  500000.,  500000.,  500000.,  500000.], 1E-8],
                         ['CHAN_FREQ',       7, [  4.32184900e+10,   4.32189900e+10,   4.32194900e+10,   4.32199900e+10,
                                                   4.32204900e+10,   4.32209900e+10,   4.32214900e+10,   4.32219900e+10,
                                                   4.32224900e+10,   4.32229900e+10,   4.32234900e+10,   4.32239900e+10,
                                                   4.32244900e+10,   4.32249900e+10,   4.32254900e+10,   4.32259900e+10], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test3(self):
        '''fitsidi-import: Test antenna name and multi uvtable'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importfitsidi(my_dataset_names[3], msname)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            _ms.close()
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True
    
            # check main table first
            name = ""
            
            #             col name, row number, expected value, tolerance
            # this is the first row of the second uv table
            expected = [
                         ['UVW',       253, [ -209524.112917, 52443.4688607, 25501.833085 ], 1E-6],
                         ['EXPOSURE',  253, 1.0, 1E-8],
                         ['WEIGHT',    253, 500000.0, 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [ 3859711.503, -201995.077, 5056134.251], 0.001],
                         ['DISH_DIAMETER',1, 15.0, 0.001],
                         ['NAME',         1, 'Kn', 0]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "SPECTRAL_WINDOW"
            expected = [ ['NUM_CHAN',        3, 512, 0],
                         ['TOTAL_BANDWIDTH', 3, 128e+6, 0],
                         ['CHAN_WIDTH',      3, [250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,
  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000.,  250000. ], 1E-8],
                         ['CHAN_FREQ',       3, [   5.87200000e+09,   5.87225000e+09,   5.87250000e+09,   5.87275000e+09,
   5.87300000e+09,   5.87325000e+09,   5.87350000e+09,   5.87375000e+09,
   5.87400000e+09,   5.87425000e+09,   5.87450000e+09,   5.87475000e+09,
   5.87500000e+09,   5.87525000e+09,   5.87550000e+09,   5.87575000e+09,
   5.87600000e+09,   5.87625000e+09,   5.87650000e+09,   5.87675000e+09,
   5.87700000e+09,   5.87725000e+09,   5.87750000e+09,   5.87775000e+09,
   5.87800000e+09,   5.87825000e+09,   5.87850000e+09,   5.87875000e+09,
   5.87900000e+09,   5.87925000e+09,   5.87950000e+09,   5.87975000e+09,
   5.88000000e+09,   5.88025000e+09,   5.88050000e+09,   5.88075000e+09,
   5.88100000e+09,   5.88125000e+09,   5.88150000e+09,   5.88175000e+09,
   5.88200000e+09,   5.88225000e+09,   5.88250000e+09,   5.88275000e+09,
   5.88300000e+09,   5.88325000e+09,   5.88350000e+09,   5.88375000e+09,
   5.88400000e+09,   5.88425000e+09,   5.88450000e+09,   5.88475000e+09,
   5.88500000e+09,   5.88525000e+09,   5.88550000e+09,   5.88575000e+09,
   5.88600000e+09,   5.88625000e+09,   5.88650000e+09,   5.88675000e+09,
   5.88700000e+09,   5.88725000e+09,   5.88750000e+09,   5.88775000e+09,
   5.88800000e+09,   5.88825000e+09,   5.88850000e+09,   5.88875000e+09,
   5.88900000e+09,   5.88925000e+09,   5.88950000e+09,   5.88975000e+09,
   5.89000000e+09,   5.89025000e+09,   5.89050000e+09,   5.89075000e+09,
   5.89100000e+09,   5.89125000e+09,   5.89150000e+09,   5.89175000e+09,
   5.89200000e+09,   5.89225000e+09,   5.89250000e+09,   5.89275000e+09,
   5.89300000e+09,   5.89325000e+09,   5.89350000e+09,   5.89375000e+09,
   5.89400000e+09,   5.89425000e+09,   5.89450000e+09,   5.89475000e+09,
   5.89500000e+09,   5.89525000e+09,   5.89550000e+09,   5.89575000e+09,
   5.89600000e+09,   5.89625000e+09,   5.89650000e+09,   5.89675000e+09,
   5.89700000e+09,   5.89725000e+09,   5.89750000e+09,   5.89775000e+09,
   5.89800000e+09,   5.89825000e+09,   5.89850000e+09,   5.89875000e+09,
   5.89900000e+09,   5.89925000e+09,   5.89950000e+09,   5.89975000e+09,
   5.90000000e+09,   5.90025000e+09,   5.90050000e+09,   5.90075000e+09,
   5.90100000e+09,   5.90125000e+09,   5.90150000e+09,   5.90175000e+09,
   5.90200000e+09,   5.90225000e+09,   5.90250000e+09,   5.90275000e+09,
   5.90300000e+09,   5.90325000e+09,   5.90350000e+09,   5.90375000e+09,
   5.90400000e+09,   5.90425000e+09,   5.90450000e+09,   5.90475000e+09,
   5.90500000e+09,   5.90525000e+09,   5.90550000e+09,   5.90575000e+09,
   5.90600000e+09,   5.90625000e+09,   5.90650000e+09,   5.90675000e+09,
   5.90700000e+09,   5.90725000e+09,   5.90750000e+09,   5.90775000e+09,
   5.90800000e+09,   5.90825000e+09,   5.90850000e+09,   5.90875000e+09,
   5.90900000e+09,   5.90925000e+09,   5.90950000e+09,   5.90975000e+09,
   5.91000000e+09,   5.91025000e+09,   5.91050000e+09,   5.91075000e+09,
   5.91100000e+09,   5.91125000e+09,   5.91150000e+09,   5.91175000e+09,
   5.91200000e+09,   5.91225000e+09,   5.91250000e+09,   5.91275000e+09,
   5.91300000e+09,   5.91325000e+09,   5.91350000e+09,   5.91375000e+09,
   5.91400000e+09,   5.91425000e+09,   5.91450000e+09,   5.91475000e+09,
   5.91500000e+09,   5.91525000e+09,   5.91550000e+09,   5.91575000e+09,
   5.91600000e+09,   5.91625000e+09,   5.91650000e+09,   5.91675000e+09,
   5.91700000e+09,   5.91725000e+09,   5.91750000e+09,   5.91775000e+09,
   5.91800000e+09,   5.91825000e+09,   5.91850000e+09,   5.91875000e+09,
   5.91900000e+09,   5.91925000e+09,   5.91950000e+09,   5.91975000e+09,
   5.92000000e+09,   5.92025000e+09,   5.92050000e+09,   5.92075000e+09,
   5.92100000e+09,   5.92125000e+09,   5.92150000e+09,   5.92175000e+09,
   5.92200000e+09,   5.92225000e+09,   5.92250000e+09,   5.92275000e+09,
   5.92300000e+09,   5.92325000e+09,   5.92350000e+09,   5.92375000e+09,
   5.92400000e+09,   5.92425000e+09,   5.92450000e+09,   5.92475000e+09,
   5.92500000e+09,   5.92525000e+09,   5.92550000e+09,   5.92575000e+09,
   5.92600000e+09,   5.92625000e+09,   5.92650000e+09,   5.92675000e+09,
   5.92700000e+09,   5.92725000e+09,   5.92750000e+09,   5.92775000e+09,
   5.92800000e+09,   5.92825000e+09,   5.92850000e+09,   5.92875000e+09,
   5.92900000e+09,   5.92925000e+09,   5.92950000e+09,   5.92975000e+09,
   5.93000000e+09,   5.93025000e+09,   5.93050000e+09,   5.93075000e+09,
   5.93100000e+09,   5.93125000e+09,   5.93150000e+09,   5.93175000e+09,
   5.93200000e+09,   5.93225000e+09,   5.93250000e+09,   5.93275000e+09,
   5.93300000e+09,   5.93325000e+09,   5.93350000e+09,   5.93375000e+09,
   5.93400000e+09,   5.93425000e+09,   5.93450000e+09,   5.93475000e+09,
   5.93500000e+09,   5.93525000e+09,   5.93550000e+09,   5.93575000e+09,
   5.93600000e+09,   5.93625000e+09,   5.93650000e+09,   5.93675000e+09,
   5.93700000e+09,   5.93725000e+09,   5.93750000e+09,   5.93775000e+09,
   5.93800000e+09,   5.93825000e+09,   5.93850000e+09,   5.93875000e+09,
   5.93900000e+09,   5.93925000e+09,   5.93950000e+09,   5.93975000e+09,
   5.94000000e+09,   5.94025000e+09,   5.94050000e+09,   5.94075000e+09,
   5.94100000e+09,   5.94125000e+09,   5.94150000e+09,   5.94175000e+09,
   5.94200000e+09,   5.94225000e+09,   5.94250000e+09,   5.94275000e+09,
   5.94300000e+09,   5.94325000e+09,   5.94350000e+09,   5.94375000e+09,
   5.94400000e+09,   5.94425000e+09,   5.94450000e+09,   5.94475000e+09,
   5.94500000e+09,   5.94525000e+09,   5.94550000e+09,   5.94575000e+09,
   5.94600000e+09,   5.94625000e+09,   5.94650000e+09,   5.94675000e+09,
   5.94700000e+09,   5.94725000e+09,   5.94750000e+09,   5.94775000e+09,
   5.94800000e+09,   5.94825000e+09,   5.94850000e+09,   5.94875000e+09,
   5.94900000e+09,   5.94925000e+09,   5.94950000e+09,   5.94975000e+09,
   5.95000000e+09,   5.95025000e+09,   5.95050000e+09,   5.95075000e+09,
   5.95100000e+09,   5.95125000e+09,   5.95150000e+09,   5.95175000e+09,
   5.95200000e+09,   5.95225000e+09,   5.95250000e+09,   5.95275000e+09,
   5.95300000e+09,   5.95325000e+09,   5.95350000e+09,   5.95375000e+09,
   5.95400000e+09,   5.95425000e+09,   5.95450000e+09,   5.95475000e+09,
   5.95500000e+09,   5.95525000e+09,   5.95550000e+09,   5.95575000e+09,
   5.95600000e+09,   5.95625000e+09,   5.95650000e+09,   5.95675000e+09,
   5.95700000e+09,   5.95725000e+09,   5.95750000e+09,   5.95775000e+09,
   5.95800000e+09,   5.95825000e+09,   5.95850000e+09,   5.95875000e+09,
   5.95900000e+09,   5.95925000e+09,   5.95950000e+09,   5.95975000e+09,
   5.96000000e+09,   5.96025000e+09,   5.96050000e+09,   5.96075000e+09,
   5.96100000e+09,   5.96125000e+09,   5.96150000e+09,   5.96175000e+09,
   5.96200000e+09,   5.96225000e+09,   5.96250000e+09,   5.96275000e+09,
   5.96300000e+09,   5.96325000e+09,   5.96350000e+09,   5.96375000e+09,
   5.96400000e+09,   5.96425000e+09,   5.96450000e+09,   5.96475000e+09,
   5.96500000e+09,   5.96525000e+09,   5.96550000e+09,   5.96575000e+09,
   5.96600000e+09,   5.96625000e+09,   5.96650000e+09,   5.96675000e+09,
   5.96700000e+09,   5.96725000e+09,   5.96750000e+09,   5.96775000e+09,
   5.96800000e+09,   5.96825000e+09,   5.96850000e+09,   5.96875000e+09,
   5.96900000e+09,   5.96925000e+09,   5.96950000e+09,   5.96975000e+09,
   5.97000000e+09,   5.97025000e+09,   5.97050000e+09,   5.97075000e+09,
   5.97100000e+09,   5.97125000e+09,   5.97150000e+09,   5.97175000e+09,
   5.97200000e+09,   5.97225000e+09,   5.97250000e+09,   5.97275000e+09,
   5.97300000e+09,   5.97325000e+09,   5.97350000e+09,   5.97375000e+09,
   5.97400000e+09,   5.97425000e+09,   5.97450000e+09,   5.97475000e+09,
   5.97500000e+09,   5.97525000e+09,   5.97550000e+09,   5.97575000e+09,
   5.97600000e+09,   5.97625000e+09,   5.97650000e+09,   5.97675000e+09,
   5.97700000e+09,   5.97725000e+09,   5.97750000e+09,   5.97775000e+09,
   5.97800000e+09,   5.97825000e+09,   5.97850000e+09,   5.97875000e+09,
   5.97900000e+09,   5.97925000e+09,   5.97950000e+09,   5.97975000e+09,
   5.98000000e+09,   5.98025000e+09,   5.98050000e+09,   5.98075000e+09,
   5.98100000e+09,   5.98125000e+09,   5.98150000e+09,   5.98175000e+09,
   5.98200000e+09,   5.98225000e+09,   5.98250000e+09,   5.98275000e+09,
   5.98300000e+09,   5.98325000e+09,   5.98350000e+09,   5.98375000e+09,
   5.98400000e+09,   5.98425000e+09,   5.98450000e+09,   5.98475000e+09,
   5.98500000e+09,   5.98525000e+09,   5.98550000e+09,   5.98575000e+09,
   5.98600000e+09,   5.98625000e+09,   5.98650000e+09,   5.98675000e+09,
   5.98700000e+09,   5.98725000e+09,   5.98750000e+09,   5.98775000e+09,
   5.98800000e+09,   5.98825000e+09,   5.98850000e+09,   5.98875000e+09,
   5.98900000e+09,   5.98925000e+09,   5.98950000e+09,   5.98975000e+09,
   5.99000000e+09,   5.99025000e+09,   5.99050000e+09,   5.99075000e+09,
   5.99100000e+09,   5.99125000e+09,   5.99150000e+09,   5.99175000e+09,
   5.99200000e+09,   5.99225000e+09,   5.99250000e+09,   5.99275000e+09,
   5.99300000e+09,   5.99325000e+09,   5.99350000e+09,   5.99375000e+09,
   5.99400000e+09,   5.99425000e+09,   5.99450000e+09,   5.99475000e+09,
   5.99500000e+09,   5.99525000e+09,   5.99550000e+09,   5.99575000e+09,
   5.99600000e+09,   5.99625000e+09,   5.99650000e+09,   5.99675000e+09,
   5.99700000e+09,   5.99725000e+09,   5.99750000e+09,   5.99775000e+09,
   5.99800000e+09,   5.99825000e+09,   5.99850000e+09,   5.99875000e+09,
   5.99900000e+09,   5.99925000e+09,   5.99950000e+09,   5.99975000e+09,], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test4(self):
        '''fitsidi-import: Test good input, list of two input files, constobsid and scanreindexing'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importfitsidi([my_dataset_names[1],my_dataset_names[2]], msname, 
                                 constobsid=True, scanreindexgap_s=1.5)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            _ms.close()
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 2.0, 1E-8],
                         ['DATA',      42,[[ 1.06945515 +7.91624188e-09j,
                                             0.98315531 +9.31322575e-10j,
                                             1.05244470 +5.77396486e-09j,
                                             0.90496856 -0.00000000e+00j,
                                             0.93005872 -6.71682887e-09j,
                                             0.80769897 -0.00000000e+00j,
                                             0.93059886 -6.97121116e-09j,
                                             0.77081358 -0.00000000e+00j,
                                             0.93020338 -7.45058060e-09j,
                                             0.83353537 -0.00000000e+00j,
                                             0.91982168 -5.54113422e-09j,
                                             0.88411278 -4.65661287e-10j,
                                             1.02857709 +5.78550630e-09j,
                                             0.93398595 -0.00000000e+00j,
                                             1.13884020 +1.01289768e-08j,
                                             2.49237108 -0.00000000e+00j ]], 1E-8],
                         ['OBSERVATION_ID', 42, 0, 0],
                         ['SCAN_NUMBER', 42, 1, 0]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
            expected = [
                         ['UVW',       638, [171529.37575288, -786712.70341456, 210321.20978818], 1E-8],
                         ['EXPOSURE',  638,  2.0, 1E-8],
                         ['DATA',      638, [[-0.00224198+0.00067056j,
                                              -0.00475123+0.0024323j,
                                              -0.00416393+0.00212671j,
                                              -0.00565350+0.00340364j,
                                              -0.00527357+0.00011977j,
                                              -0.00292699+0.00131954j,
                                              -0.00429945+0.00035823j,
                                              -0.00545671-0.00033945j,
                                              -0.00646004+0.00037293j,
                                              -0.00419376-0.00115011j,
                                              -0.00508117+0.00045939j,
                                              -0.00501660-0.00047975j,
                                              -0.00444734-0.00101535j,
                                              -0.00384988-0.00102731j,
                                              -0.00551326+0.00101364j,
                                              -0.00337701+0.00080481j]], 1E-8],
                         ['OBSERVATION_ID', 638, 0, 0],
                         ['SCAN_NUMBER', 638, 14, 0]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            _tb.open(msname+'/OBSERVATION')
            nr = _tb.nrows()
            _tb.close()
            if not nr==1:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table OBSERVATION failed'
                
        self.assertTrue(retValue['success'])
<<<<<<< HEAD
    

    def test5(self):
        '''fitsidi-import: Test e-MERLIN polarization swapping'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        # FITS-IDI files from the e-MERLIN correlator have some
        # baselines in non-canonical order.  Make sure that
        # importfitsidi correctly swaps the cross-polarisation
        # products when it swaps the antennas.
        #
        self.res = importfitsidi(my_dataset_names[4], msname)
        print(myname, ": Success! Now checking output ...")

        # The Mk2-De baseline will be swapped.  These are the expected
        # phases for that baseline.
        expected = numpy.array([[ 3.13964248,  2.86829948,  2.58296013,  2.3439014 ,  2.10988331,
                                  1.97421145,  1.7791388 ,  1.5775888 ,  1.40460491,  1.22053266,
                                  1.01471496,  0.82222307,  0.59186804,  0.38555288,  0.16937262,
                                  -0.05033521, -0.25431359, -0.47900817, -0.69958377, -0.87468606,
                                  -1.09604585, -1.34909344, -1.62127054, -1.87543309, -2.12334824,
                                  -2.34573412, -2.53386259, -2.72510219, -2.89187169, -3.06771994,
                                  3.0699141 ,  2.88419271],
                                [-1.84054375, -0.83255982, -1.02057993, -1.24484885, -1.32794607,
                                  -0.82715571, -0.79155159, -1.19584942, -0.73648691, -0.59470183,
                                  -0.83499312, -0.48545197, -0.57062495, -0.34914687,  0.01485419,
                                  -0.17965424,  0.01576512,  0.64190328,  0.36306599,  0.6468659 ,
                                  1.15446103,  1.02624571,  0.89803523,  1.0682807 ,  0.83425844,
                                  0.8939116 ,  0.96772492,  1.20877922,  1.51828313,  1.73502862,
                                  1.48054171,  2.27795124],
                                [-1.72998369, -2.10472918, -2.28158116, -2.53716612, -2.63328958,
                                  -2.69052339, -2.89023638, -3.02556944,  3.04760623,  2.90053821,
                                  2.76279259,  2.49925351,  2.26630735,  1.77548301,  1.423805  ,
                                  1.0532515 ,  0.66610861,  0.37606233,  0.14707069,  0.03789629,
                                  -0.18920378, -0.37956986, -0.56126845, -0.69712436, -0.82414556,
                                  -1.03143954, -1.22291136, -1.41123331, -1.65069032, -1.79867935,
                                  -1.80994236, -2.02587509],
                                [ 0.1903225 ,  0.28886324,  0.27794465,  0.29231161,  0.298327  ,
                                  0.34109077,  0.45782065,  0.60621154,  0.77383608,  0.90957648,
                                  1.06189084,  1.24068403,  1.39028943,  1.51727164,  1.60704148,
                                  1.63228989,  1.67657077,  1.71906674,  1.80779397,  1.9613409 ,
                                  2.11035109,  2.254354  ,  2.42930079,  2.55850363,  2.67686605,
                                  2.78933764,  2.83647609,  2.92243385, -3.14000082, -2.84540629,
                                  -2.57244778, -2.36642051]])

        try:
            _ms.open(msname)
            _ms.select({'antenna1': [0], 'antenna2': [2]})
            _ms.selectchannel(nchan=32, start=0, width=16, inc=16)
            rec = _ms.getdata(['phase'], average=True)
            _ms.close()
        except:
            print(myname, ": Error  Cannot open MS")
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS'
        else:
            results = numpy.isclose(rec['phase'], expected, rtol=8e-7, atol=1e-8).all()
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Crosspol check for Mk2-De baseline failed'

        # The De-Pi baseline will not be swapped.  These are the expected
        # phases for that baseline.
        expected = numpy.array([[ 2.18404031,  2.52688265,  2.88025069, -3.06777024, -2.7234993 ,
                                  -2.4769671 , -2.12538004, -1.74727786, -1.40869737, -1.12853301,
                                  -0.83723795, -0.57250589, -0.23436357,  0.10187802,  0.42698029,
                                  0.77069479,  1.06355178,  1.38657355,  1.72898543,  2.03018403,
                                  2.38947916,  2.717098  ,  3.01598072, -2.97815442, -2.67935085,
                                  -2.41194654, -2.12304139, -1.84744442, -1.55729675, -1.25337708,
                                  -1.00634491, -0.74762005],
                                [ 2.4460454 ,  3.00888586, -3.08539438, -3.0410893 , -2.99218869,
                                  -2.94937682, -2.57229424, -2.26054573, -1.91930437, -1.58734751,
                                  -1.2868315 , -0.96962148, -0.59307194, -0.25247654,  0.13383511,
                                  0.34873137,  0.54140246,  0.77145398,  1.05496573,  1.18506956,
                                  1.47813749,  1.70011413,  1.94997513,  2.24114943,  2.3705833 ,
                                  2.51808047,  2.82756138,  3.0124886 , -2.97727132, -2.61935139,
                                  -2.24914908, -1.8698504 ],
                                [ 1.28013647,  1.47956204,  1.30087411,  1.58066428,  1.52697527,
                                  1.57440555,  1.64748108,  1.613729  ,  1.63689208,  1.7030493 ,
                                  1.54547369,  1.63134134,  1.5600487 ,  1.57901645,  1.56169891,
                                  1.65976775,  1.65306389,  1.659922  ,  1.46294129,  1.4702158 ,
                                  1.38168573,  1.39496899,  1.37612033,  1.25856972,  1.31672955,
                                  1.34671378,  1.25240135,  1.17149973,  1.25372207,  0.98874217,
                                  0.85313851,  0.67320693],
                                [-0.2974771 , -0.27914691, -0.29463151, -0.35890678, -0.41856381,
                                  -0.50046909, -0.56969333, -0.62361556, -0.73453844, -0.78436178,
                                  -0.84025633, -0.90566015, -0.98629397, -1.05542088, -1.12567806,
                                  -1.16428077, -1.23636246, -1.32543719, -1.40451121, -1.51268327,
                                  -1.60201359, -1.69172263, -1.74908674, -1.83386433, -1.8914839 ,
                                  -1.954831  , -2.02388644, -2.12116575, -2.21479368, -2.31664228,
                                  -2.38661528, -2.48180223]])

        try:
            _ms.open(msname)
            _ms.select({'antenna1': [2], 'antenna2': [3]})
            _ms.selectchannel(nchan=32, start=0, width=16, inc=16)
            rec = _ms.getdata(['phase'], average=True)
            _ms.close()
        except:
            print(myname, ": Error  Cannot open MS")
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS'
        else:
            results = numpy.isclose(rec['phase'], expected, rtol=8e-7, atol=1e-8).all()
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Crosspol check for De-Pi baseline failed'

        self.assertTrue(retValue['success'])

=======
>>>>>>> master
    
    def test5(self):
        '''fitsidi-import: Test e-MERLIN polarization swapping'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        # FITS-IDI files from the e-MERLIN correlator have some
        # baselines in non-canonical order.  Make sure that
        # importfitsidi correctly swaps the cross-polarisation
        # products when it swaps the antennas.
        #
        self.res = importfitsidi(my_dataset_names[4], msname)
        print(myname, ": Success! Now checking output ...")

        # The Mk2-De baseline will be swapped.  These are the expected
        # phases for that baseline.
        expected = numpy.array([[ 3.13964248,  2.86829948,  2.58296013,  2.3439014 ,  2.10988331,
                                  1.97421145,  1.7791388 ,  1.5775888 ,  1.40460491,  1.22053266,
                                  1.01471496,  0.82222307,  0.59186804,  0.38555288,  0.16937262,
                                  -0.05033521, -0.25431359, -0.47900817, -0.69958377, -0.87468606,
                                  -1.09604585, -1.34909344, -1.62127054, -1.87543309, -2.12334824,
                                  -2.34573412, -2.53386259, -2.72510219, -2.89187169, -3.06771994,
                                  3.0699141 ,  2.88419271],
                                [-1.84054375, -0.83255982, -1.02057993, -1.24484885, -1.32794607,
                                  -0.82715571, -0.79155159, -1.19584942, -0.73648691, -0.59470183,
                                  -0.83499312, -0.48545197, -0.57062495, -0.34914687,  0.01485419,
                                  -0.17965424,  0.01576512,  0.64190328,  0.36306599,  0.6468659 ,
                                  1.15446103,  1.02624571,  0.89803523,  1.0682807 ,  0.83425844,
                                  0.8939116 ,  0.96772492,  1.20877922,  1.51828313,  1.73502862,
                                  1.48054171,  2.27795124],
                                [-1.72998369, -2.10472918, -2.28158116, -2.53716612, -2.63328958,
                                  -2.69052339, -2.89023638, -3.02556944,  3.04760623,  2.90053821,
                                  2.76279259,  2.49925351,  2.26630735,  1.77548301,  1.423805  ,
                                  1.0532515 ,  0.66610861,  0.37606233,  0.14707069,  0.03789629,
                                  -0.18920378, -0.37956986, -0.56126845, -0.69712436, -0.82414556,
                                  -1.03143954, -1.22291136, -1.41123331, -1.65069032, -1.79867935,
                                  -1.80994236, -2.02587509],
                                [ 0.1903225 ,  0.28886324,  0.27794465,  0.29231161,  0.298327  ,
                                  0.34109077,  0.45782065,  0.60621154,  0.77383608,  0.90957648,
                                  1.06189084,  1.24068403,  1.39028943,  1.51727164,  1.60704148,
                                  1.63228989,  1.67657077,  1.71906674,  1.80779397,  1.9613409 ,
                                  2.11035109,  2.254354  ,  2.42930079,  2.55850363,  2.67686605,
                                  2.78933764,  2.83647609,  2.92243385, -3.14000082, -2.84540629,
                                  -2.57244778, -2.36642051]])

        try:
            _ms.open(msname)
            _ms.select({'antenna1': [0], 'antenna2': [2]})
            _ms.selectchannel(nchan=32, start=0, width=16, inc=16)
            rec = _ms.getdata(['phase'], average=True)
            _ms.close()
        except:
            print(myname, ": Error  Cannot open MS")
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS'
        else:
            results = numpy.isclose(rec['phase'], expected, rtol=8e-7, atol=1e-8).all()
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Crosspol check for Mk2-De baseline failed'

        # The De-Pi baseline will not be swapped.  These are the expected
        # phases for that baseline.
        expected = numpy.array([[ 2.18404031,  2.52688265,  2.88025069, -3.06777024, -2.7234993 ,
                                  -2.4769671 , -2.12538004, -1.74727786, -1.40869737, -1.12853301,
                                  -0.83723795, -0.57250589, -0.23436357,  0.10187802,  0.42698029,
                                  0.77069479,  1.06355178,  1.38657355,  1.72898543,  2.03018403,
                                  2.38947916,  2.717098  ,  3.01598072, -2.97815442, -2.67935085,
                                  -2.41194654, -2.12304139, -1.84744442, -1.55729675, -1.25337708,
                                  -1.00634491, -0.74762005],
                                [ 2.4460454 ,  3.00888586, -3.08539438, -3.0410893 , -2.99218869,
                                  -2.94937682, -2.57229424, -2.26054573, -1.91930437, -1.58734751,
                                  -1.2868315 , -0.96962148, -0.59307194, -0.25247654,  0.13383511,
                                  0.34873137,  0.54140246,  0.77145398,  1.05496573,  1.18506956,
                                  1.47813749,  1.70011413,  1.94997513,  2.24114943,  2.3705833 ,
                                  2.51808047,  2.82756138,  3.0124886 , -2.97727132, -2.61935139,
                                  -2.24914908, -1.8698504 ],
                                [ 1.28013647,  1.47956204,  1.30087411,  1.58066428,  1.52697527,
                                  1.57440555,  1.64748108,  1.613729  ,  1.63689208,  1.7030493 ,
                                  1.54547369,  1.63134134,  1.5600487 ,  1.57901645,  1.56169891,
                                  1.65976775,  1.65306389,  1.659922  ,  1.46294129,  1.4702158 ,
                                  1.38168573,  1.39496899,  1.37612033,  1.25856972,  1.31672955,
                                  1.34671378,  1.25240135,  1.17149973,  1.25372207,  0.98874217,
                                  0.85313851,  0.67320693],
                                [-0.2974771 , -0.27914691, -0.29463151, -0.35890678, -0.41856381,
                                  -0.50046909, -0.56969333, -0.62361556, -0.73453844, -0.78436178,
                                  -0.84025633, -0.90566015, -0.98629397, -1.05542088, -1.12567806,
                                  -1.16428077, -1.23636246, -1.32543719, -1.40451121, -1.51268327,
                                  -1.60201359, -1.69172263, -1.74908674, -1.83386433, -1.8914839 ,
                                  -1.954831  , -2.02388644, -2.12116575, -2.21479368, -2.31664228,
                                  -2.38661528, -2.48180223]])

        try:
            _ms.open(msname)
            _ms.select({'antenna1': [2], 'antenna2': [3]})
            _ms.selectchannel(nchan=32, start=0, width=16, inc=16)
            rec = _ms.getdata(['phase'], average=True)
            _ms.close()
        except:
            print(myname, ": Error  Cannot open MS")
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS'
        else:
            results = numpy.isclose(rec['phase'], expected, rtol=8e-7, atol=1e-8).all()
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Crosspol check for De-Pi baseline failed'

        self.assertTrue(retValue['success'])

    def test6(self):
        '''fitsidi-import: Test import of gain curves'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importfitsidi(my_dataset_names[5], msname, 
                                 constobsid=True, scanreindexgap_s=5)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "GAIN_CURVE/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            _ms.close()
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

            name = "GAIN_CURVE"
            expected = [ ['TYPE',       8, 'POWER(ZA)', 0],
                         ['NUM_POLY',   8, 3, 0],
                         ['GAIN',       8, [[0.80699998, 0.01596000,
                                             -0.00020470],
                                            [0.80699998, 0.01596000,
                                             -0.00020470]], 1E-8],
                         ['SENSITIVITY', 8, [0.07699999, 0.06499999], 1E-8],
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

        self.assertTrue(retValue['success'])

def suite():
    return [test_importfitsidi]
    
if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

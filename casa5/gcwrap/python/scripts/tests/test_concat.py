#############################################################################
# $Id:$
# Test Name:                                                                #
#    Unit Test Script for the concat task
#    
#                                                                           #
#############################################################################
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import glob
import unittest
from math import sqrt

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, calibrater
    from casatools import table as tbtool
    from casatools import ms as mstool
    from casatasks import split, concat

    cb = calibrater( )
    tb = tbtool( )
    ms = mstool( )

    ctsys_resolve = ctsys.resolve
else:
    from __main__ import default
    from tasks import *
    from taskinit import *

    cb = cbtool( )

    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata')
        return os.path.join(dataPath,apath)

myname = 'test_concat'

# name of the resulting MS
msname = 'concatenated.ms'

testmms=False

def checktable(thename, theexpectation, multims=False):
    global msname, myname
    if multims:        
        tb.open(msname+"/SUBMSS/"+thename)
    else:
        tb.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print(myname, ": comparing ", mycell)
        value = tb.getcell(mycell[0], mycell[1])
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
            # it's an array
            if type(value) == str: # it's actually a string
                try:
                    # zero tolerance?
                    if mycell[3] == 0:
                        in_agreement = (value == mycell[2])
                    else:
                        in_agreement = (mycell[2] in value ) # i.e. the expectation is contained
                except:
                    in_agreement = False
            else:
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
            tb.close()
            return False
    tb.close()
    print(myname, ": table ", thename, " as expected.")
    return True


###########################
# beginning of actual test 

class test_concat(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         '''A class method called before tests in an individual class run'''
#         # MS used for all tests of this class.
#         cls.msfile = 'vlass1.1_rowtest_6fields.ms'
#         os.symlink(data_path+cls.msfile, cls.msfile)
        
    @classmethod
    def setUpClass(cls):
        '''Copy input files only once for the class'''
        global testmms
        res = None

        datapath=ctsys_resolve('unittest/concat/')
        datapathmms = ''
        # Pick up alternative data directory to run tests on MMSs
        testmms = False
        if 'TEST_DATADIR' in os.environ:
            testmms = True
            DATADIR = str(os.environ.get('TEST_DATADIR'))
            if os.path.isdir(DATADIR):
                datapathmms = DATADIR+'/concat/input/'


        cpath = os.path.abspath(os.curdir)
        filespresent = sorted(glob.glob("*.ms"))
        if(datapathmms!=''): 
            print("\nTesting on MMSs ...\n")

            nonmmsinput = ['A2256LC2_4.5s-1.ms', 'part1.ms', 'part2-mod2.ms', 'shortpart1.ms', 'sim7.ms', 
                           'uid___A002_Xab8dc1_X95a-shrunk.ms', 'uid___A002_Xab8dc1_Xf13-shrunk.ms',
                           'nep1-shrunk.ms', 'nep2-shrunk.ms', 'nep3-shrunk.ms',
                           'X1127a.smaller.ms', 'Xd21.smaller.ms']
            os.chdir(datapathmms)
            myinputmslist = sorted(glob.glob("*.ms"))
            os.chdir(cpath)

            for mymsname in myinputmslist:
                if not ((mymsname in filespresent) or (mymsname in nonmmsinput)):
                    print("Copying MMS", mymsname)
                    rval = os.system('cp -R '+datapathmms+'/'+mymsname+' .')
                    if rval!=0:
                        raise Exception('Error while copying input data.')

            for mymsname in nonmmsinput:
                if not mymsname in filespresent:
                    print("Copying non-MMS ", mymsname)
                    rval = os.system('cp -R '+datapath+'/'+mymsname+' .')
                    if rval!=0:
                        raise Exception('Error while copying input data.')
                    
        else:
            os.chdir(datapath)
            myinputmslist = sorted(glob.glob("*.ms"))
            print(myinputmslist)
            os.chdir(cpath)
            for mymsname in myinputmslist:
                if not mymsname in filespresent:
                    print("Copying ", mymsname)
                    rval = os.system('cp -RL '+os.path.join(datapath,mymsname)+' .')
                    if rval!=0:
                        raise Exception('Error while copying input data.')

        os.chdir(cpath)

        # create MSs with ephemeris use
        
        if not 'xy1.ms' in filespresent:
            print("Creating MSs with ephemeris table use ...")

            split(vis='part1.ms', outputvis='xy1.ms', scan="1~30", datacolumn='data')

            tb.open('xy1.ms', nomodify=False)
            a = tb.getcol('TIME')
            delta = (54709.*86400-a[0])
            a = a + delta
            tb.putcol('TIME', a)
            a = tb.getcol('TIME_CENTROID')
            a = a + delta
            tb.putcol('TIME_CENTROID', a)
            tb.close()
            tb.open('xy1.ms/FIELD', nomodify=False)
            a = tb.getcol('TIME')
            a = a + delta
            tb.putcol('TIME', a)
            tb.close()
            tb.open('xy1.ms/OBSERVATION', nomodify=False)
            a = tb.getcol('TIME_RANGE')
            a = a + delta
            tb.putcol('TIME_RANGE', a)
            tb.close()

            shutil.copytree('xy1.ms', 'xy1-noephem.ms')

            ms.open('xy1.ms', nomodify=False)
            ms.addephemeris(0,ctsys_resolve('ephemerides/JPL-Horizons/Uranus_54708-55437dUTC.tab'),
                            'Uranus_54708-55437dUTC', '1908-201')  # this field is not really Uranus but for a test this doesn't matter
            ms.addephemeris(1,ctsys_resolve('ephemerides/JPL-Horizons/Jupiter_54708-55437dUTC.tab'),
                            'Jupiter_54708-55437dUTC', 0)
            ms.close()
            myinputmslist.append('xy1.ms')
            myinputmslist.append('xy1-noephem.ms')

        if not 'xy2.ms' in filespresent:
            split(vis='part1.ms', outputvis='xy2.ms', scan="58~65", datacolumn='data')

            tb.open('xy2.ms', nomodify=False)
            a = tb.getcol('TIME')
            delta = (54719.*86400-a[0]) # ten days later than xy1.ms !
            a = a + delta
            tb.putcol('TIME', a)
            a = tb.getcol('TIME_CENTROID')
            a = a + delta
            tb.putcol('TIME_CENTROID', a)
            tb.close()
            tb.open('xy2.ms/FIELD', nomodify=False)
            a = tb.getcol('TIME')
            a = a + delta
            tb.putcol('TIME', a)
            tb.close()
            tb.open('xy2.ms/OBSERVATION', nomodify=False)
            a = tb.getcol('TIME_RANGE')
            a = a + delta
            tb.putcol('TIME_RANGE', a)
            tb.close()

            ms.open('xy2.ms', nomodify=False)
            ms.addephemeris(0,ctsys_resolve('ephemerides/JPL-Horizons/Uranus_54708-55437dUTC.tab'),
                            'Uranus_54708-55437dUTC', '1908-201')
            ms.close()
            myinputmslist.append('xy2.ms')

        if not 'xy2late.ms' in filespresent:
            split(vis='part1.ms', outputvis='xy2late.ms', scan="58~65", datacolumn='data')

            tb.open('xy2late.ms', nomodify=False)
            a = tb.getcol('TIME')
            delta = (55438.*86400-a[0]) # much later than xy1.ms, beyond the end of Uranus_54708-55437dUTC.tab !
            a = a + delta
            tb.putcol('TIME', a)
            a = tb.getcol('TIME_CENTROID')
            a = a + delta
            tb.putcol('TIME_CENTROID', a)
            tb.close()
            tb.open('xy2late.ms/FIELD', nomodify=False)
            a = tb.getcol('TIME')
            a = a + delta
            tb.putcol('TIME', a)
            tb.close()
            tb.open('xy2late.ms/OBSERVATION', nomodify=False)
            a = tb.getcol('TIME_RANGE')
            a = a + delta
            tb.putcol('TIME_RANGE', a)
            tb.close()

            ms.open('xy2late.ms', nomodify=False)
            ms.addephemeris(0,ctsys_resolve('ephemerides/JPL-Horizons/Uranus_55437-56293dUTC.tab'),
                            'Uranus_55437-56293dUTC', '1908-201')
            ms.close()
            myinputmslist.append('xy2late.ms')

        if not 'xy2-jup-ur.ms' in filespresent:
            split(vis='xy2.ms', outputvis='xy2-jup-ur.ms', field = 'jupiter, 1908-201', datacolumn='data')
            myinputmslist.append('xy2-jup-ur.ms')

        if not 'xya.ms' in filespresent:
            split(vis='xy1.ms', outputvis='xya.ms', spw='0:0~63', datacolumn='data')
            myinputmslist.append('xya.ms')
            
        if not 'xyb.ms' in filespresent:
            split(vis='xy1.ms', outputvis='xyb.ms', spw='0:64~127', datacolumn='data')
            myinputmslist.append('xyb.ms')

        cls.listoffiles = myinputmslist

        if not is_CASA6:
            default(concat)
        return True

        
    def tearDown(self):
        shutil.rmtree(msname,ignore_errors=True)

    @classmethod
    def tearDownClass(cls):
        '''A class method called after tests in an individual class have run'''
        # remove all the generated data
        for ifile in cls.listoffiles:
            print(ifile)
            os.system('rm -rf '+ifile)
            shutil.rmtree(ifile,ignore_errors=True)

    def test1(self):
        '''Concat 1: 4 parts, same sources but different spws'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    
        
        concat(vis=['part1.ms','part2.ms','part3.ms','part4.ms'],concatvis=msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test1.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test1.ms',ignore_errors=True)
#             shutil.copytree(msname,'test1.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

            # check source table
            name = "SOURCE"
            #             col name, row number, expected value, tolerance
            expected = [
                ['SOURCE_ID',           55, 13, 0],
                ['SPECTRAL_WINDOW_ID',  55, 3, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            # check spw table
            name = "SPECTRAL_WINDOW"
            #             col name, row number, expected value, tolerance
            expected = [
                ['NUM_CHAN',           3, 128, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

        self.assertTrue(retValue['success'])


    def test2(self):
        '''Concat 2: 3 parts, different sources, different spws,  copypointing=False, visweightscale=[3.,2.,1.]'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        concat(vis=['part1.ms','part2-mod.ms','part3.ms'],concatvis=msname,
               copypointing=False, visweightscale=[3.,2.,1.])
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test2.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test2.ms',ignore_errors=True)
#             shutil.copytree(msname,'test2.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

            # check source table
            name = "SOURCE"
            #             col name, row number, expected value, tolerance
            expected = [
                ['SOURCE_ID',           41, 13, 0],
                ['SPECTRAL_WINDOW_ID',  41, 2, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            # check spw table
            name = "SPECTRAL_WINDOW"
            #             col name, row number, expected value, tolerance
            expected = [
                ['NUM_CHAN',           2, 128, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

            # collecting parameters for subsequent test of MAIN table
            msnrows = []
            oldweightbeg = []
            oldweightend = []
            oldsigmabeg = []
            oldsigmaend = []
            ii = 0
            for myms in ['part1.ms','part2-mod.ms','part3.ms']:
                tb.open(myms)
                msnrows.append(tb.nrows())
                oldweightbeg.append(tb.getcell('WEIGHT',0))
                oldweightend.append(tb.getcell('WEIGHT',tb.nrows()-1))
                oldsigmabeg.append(tb.getcell('SIGMA',0))
                oldsigmaend.append(tb.getcell('SIGMA',tb.nrows()-1))
                tb.close()


            name = "" # i.e. Main
            #             col name, row number, expected value, tolerance
            expected = [
                    ['WEIGHT', 0, 3.*oldweightbeg[0], 1E-6], # scaling uses float precision
                    ['WEIGHT', msnrows[0]-1, 3.*oldweightend[0], 1E-6],
                    ['WEIGHT', msnrows[0], 2.*oldweightbeg[1], 1E-6],
                    ['WEIGHT', msnrows[0]+msnrows[1]-1, 2.*oldweightend[1], 1E-6],
                    ['WEIGHT', msnrows[0]+msnrows[1], oldweightbeg[2], 1E-6],
                    ['WEIGHT', msnrows[0]+msnrows[1]+msnrows[2]-1, oldweightend[2], 1E-6]
                ]

            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

            expected2 = [
                    ['SIGMA', 0, 1./sqrt(3.)*oldsigmabeg[0], 1E-6], # scaling uses float precision
                    ['SIGMA', msnrows[0]-1, 1./sqrt(3.)*oldsigmaend[0], 1E-6],
                    ['SIGMA', msnrows[0], 1./sqrt(2.)*oldsigmabeg[1], 1E-6],
                    ['SIGMA', msnrows[0]+msnrows[1]-1, 1/sqrt(2.)*oldsigmaend[1], 1E-6],
                    ['SIGMA', msnrows[0]+msnrows[1], oldsigmabeg[2], 1E-6],
                    ['SIGMA', msnrows[0]+msnrows[1]+msnrows[2]-1, oldsigmaend[2], 1E-6]
                ]

            results2 = checktable(name, expected2)
            if not results2:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

        self.assertTrue(retValue['success'])


    def test3(self):
        '''Concat 3: 3 parts, different sources, same spws'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    
        concat(vis=['part1.ms','part2-mod2.ms','part3.ms'],concatvis=msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test3.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test3.ms',ignore_errors=True)
#             shutil.copytree(msname,'test3.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

            # check source table
            name = "SOURCE"
            #             col name, row number, expected value, tolerance
            expected = [
                ['SOURCE_ID',           28, 13, 0],
                ['SPECTRAL_WINDOW_ID',  28, 1, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            # check spw table
            name = "SPECTRAL_WINDOW"
            #             col name, row number, expected value, tolerance
            expected = [
                ['NUM_CHAN',           1, 128, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

        self.assertTrue(retValue['success'])


    def test4(self):
        '''Concat 4: five MSs with identical sources but different time/intervals on them (CSV-268)'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis = ['shortpart1.ms', 'shortpart2.ms', 'shortpart3.ms', 'shortpart4.ms',
                      'shortpart5.ms'],
               concatvis = msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test4.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test4.ms',ignore_errors=True)
#             shutil.copytree(msname,'test4.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True
        
    
            # check source table
            name = "SOURCE"
            #             col name, row number, expected value, tolerance
            expected = [
                ['SOURCE_ID',           0, 0, 0],
                ['SPECTRAL_WINDOW_ID',  0, 0, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            expected = [
                ['SOURCE_ID',           7, 0, 0],
                ['SPECTRAL_WINDOW_ID',  7, 7, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            expected = [
                ['SOURCE_ID',           8, 1, 0],
                ['SPECTRAL_WINDOW_ID',  8, 0, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            expected = [
                ['SOURCE_ID',           15, 1, 0],
                ['SPECTRAL_WINDOW_ID',  15, 7, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            expected = [
                ['SOURCE_ID',           16, 0, 100000],
                ['SPECTRAL_WINDOW_ID',  16, 0, 100000]
                ]
            print("The following should fail: SOURCE row 16 should not exist")
            try:
                results = checktable(name, expected)
            except:
                print("Expected error.")
                results = False
            if results: 
                retValue['success']=False
                retValue['error_msgs']='SOURCE row 16 should not existCheck of table '+name+' failed'
            # check spw table
            name = "SPECTRAL_WINDOW"
            #             col name, row number, expected value, tolerance
            expected = [
                ['NUM_CHAN',           8, 4, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
                
        self.assertTrue(retValue['success'])
        
    def test5(self):
        '''Concat 5: two MSs with different state table (CAS-2601)'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis = ['A2256LC2_4.5s-1.ms','A2256LC2_4.5s-2.ms'], concatvis = msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test5.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test5.ms',ignore_errors=True)
#             shutil.copytree(msname,'test5.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check state table
            name = "STATE"
            #             col name, row number, expected value, tolerance
            expected = [
                ['CAL',  0, 0, 0],
                ['SIG',  0, 1, 0],
                ['SUB_SCAN',  2, 1, 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test6(self):
        '''Concat 6: two MSs with different state table and feed table'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis = ['A2256LC2_4.5s-1.ms','A2256LC2_4.5s-2b.ms'], concatvis = msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test6.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test6.ms',ignore_errors=True)
#             shutil.copytree(msname,'test6.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check FEED table
            name = "FEED"
            #             col name, row number, expected value, tolerance
            expected = [
                ['SPECTRAL_WINDOW_ID',  53, 1, 0],
                ['SPECTRAL_WINDOW_ID',  54, 2, 0],
                ['SPECTRAL_WINDOW_ID',  107, 3, 0],
                ['RECEPTOR_ANGLE',  54, [-1,0], 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test7(self):
        '''Concat 7: two MSs with different antenna table such that baseline label reversal becomes necessary'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis = ['sim7.ms','sim8.ms'], concatvis = msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test7.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test7.ms',ignore_errors=True)
#             shutil.copytree(msname,'test7.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check Main table
#            tb.open('test7.ms')
            tb.open(msname)
            ant1 = tb.getcol('ANTENNA1')
            ant2 = tb.getcol('ANTENNA2')
            tb.close()
            result = True
            print(myname, ": OK. Checking baseline labels ...")
            for i in range(0,len(ant1)):
                if(ant1[i]>ant2[i]):
                    print("Found incorrectly ordered baseline label in row ", i, ": ", ant1, " ", ant2)
                    result = False
                    break

            if not result:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test8(self):
        '''Concat 8: two MSs with different antenna tables, copypointing = False'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis = ['sim7.ms','sim8.ms'], concatvis = msname, copypointing=False)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()                
#             if 'test8.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test8.ms',ignore_errors=True)
#             shutil.copytree(msname,'test8.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check Main table
#            tb.open('test8.ms')
            tb.open(msname)
            ant1 = tb.getcol('ANTENNA1')
            ant2 = tb.getcol('ANTENNA2')
            tb.close()
            result = True
            print(myname, ": OK. Checking baseline labels ...")
            for i in range(0,len(ant1)):
                if(ant1[i]>ant2[i]):
                    print("Found incorrectly ordered baseline label in row ", i, ": ", ant1, " ", ant2)
                    result = False
                    break

            if not result:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

        self.assertTrue(retValue['success'])



    def test9(self):
        '''Concat 9: 3 parts, different sources, same spws, different scratch columns: no, yes, no'''
        global testmms
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        shutil.rmtree('part2-mod2-wscratch.ms',ignore_errors=True)
        shutil.copytree('part2-mod2.ms', 'part2-mod2-wscratch.ms')
        print('creating scratch columns in part2-mod2-wscratch.ms')
        cb.open('part2-mod2-wscratch.ms') # calibrator-open creates scratch columns
        cb.close()
        
        if testmms:
            print("Expecting an Error ... ")

        concat(vis=['part1.ms','part2-mod2-wscratch.ms','part3.ms'],concatvis=msname)
        shutil.rmtree('part2-mod2-wscratch.ms',ignore_errors=True)

        if not testmms:

            print(myname, ": Now checking output ...")
            mscomponents = set(["table.dat",
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
                                "SOURCE/table.dat",
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
                                "SOURCE/table.f0",
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
            self.assertTrue(retValue['success'])
            print(myname, ": MS exists. All tables present. Try opening as MS ...")
            try:
                ms.open(msname)
            except:
                print(myname, ": Error  Cannot open MS table", tablename)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
            else:
                ms.close()
#                 if 'test9.ms' in glob.glob("*.ms"):
#                     shutil.rmtree('test9.ms',ignore_errors=True)
#                 shutil.copytree(msname,'test9.ms')
                print(myname, ": OK. Checking tables in detail ...")
                retValue['success']=True

                # check source table
                name = "SOURCE"
                #             col name, row number, expected value, tolerance
                expected = [
                    ['SOURCE_ID',           28, 13, 0],
                    ['SPECTRAL_WINDOW_ID',  28, 1, 0]
                    ]
                results = checktable(name, expected)
                if not results:
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                # check spw table
                name = "SPECTRAL_WINDOW"
                #             col name, row number, expected value, tolerance
                expected = [
                    ['NUM_CHAN',           1, 128, 0]
                    ]
                results = checktable(name, expected)
                if not results:
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

            self.assertTrue(retValue['success'])
        # endif not testmms
        
    def test10(self):
        '''Concat 10: 3 parts, different sources, same spws, different scratch columns: yes, no, no'''
        global testmms
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        if testmms:
            print("Skipping this test as concat cannot create scratch columns in MMSs.")
        else:
            shutil.rmtree('part1-wscratch.ms',ignore_errors=True)
            shutil.copytree('part1.ms', 'part1-wscratch.ms')
            print('creating scratch columns in part1-wscratch.ms')
            cb.open('part1-wscratch.ms') # calibrator-open creates scratch columns
            cb.close()

            concat(vis=['part1-wscratch.ms','part2-mod2.ms','part3.ms'],concatvis=msname)
            shutil.rmtree('part1-wscratch.ms',ignore_errors=True)

            print(myname, ": Now checking output ...")
            mscomponents = set(["table.dat",
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
                                "SOURCE/table.dat",
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
                                "SOURCE/table.f0",
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
            self.assertTrue(retValue['success'])
            print(myname, ": MS exists. All tables present. Try opening as MS ...")
            try:
                ms.open(msname)
            except:
                print(myname, ": Error  Cannot open MS table", tablename)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
            else:
                ms.close()
#                 if 'test10.ms' in glob.glob("*.ms"):
#                     shutil.rmtree('test10.ms',ignore_errors=True)
#                 shutil.copytree(msname,'test10.ms')
                print(myname, ": OK. Checking tables in detail ...")
                retValue['success']=True

                # check source table
                name = "SOURCE"
                #             col name, row number, expected value, tolerance
                expected = [
                    ['SOURCE_ID',           28, 13, 0],
                    ['SPECTRAL_WINDOW_ID',  28, 1, 0]
                    ]
                results = checktable(name, expected)
                if not results:
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                # check spw table
                name = "SPECTRAL_WINDOW"
                #             col name, row number, expected value, tolerance
                expected = [
                    ['NUM_CHAN',           1, 128, 0]
                    ]
                results = checktable(name, expected)
                if not results:
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

            self.assertTrue(retValue['success'])
        # endif

    def test11(self):
        '''Concat 11: 2 parts of same MS split in time,  use of ephemerides, first ephemeris covers both MS time ranges'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['xy1.ms','xy2.ms'],concatvis=msname, copypointing=False)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Uranus_54708-55437dUTC.tab",
                            "FIELD/EPHEM1_Jupiter_54708-55437dUTC.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test11.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test11.ms',ignore_errors=True)
#             shutil.copytree(msname,'test11.ms')
            print(myname, ": OK. Checking tables in detail ...")
#            tb.open('test11.ms/FIELD')
            tb.open(msname+'/FIELD')
            a = list(tb.getcol('NAME'))
            tb.close()
            tb.open('xy1.ms/FIELD')
            compa = list(tb.getcol('NAME'))
            compa.append('jupiter') # jupiter should occur a second time because in xy1.ms it is an ephemeris object, in xy2 it is not
            tb.close()
            retValue['success']=True
            if not (len(a)==len(compa) and a==compa):
                print("FIELD table of test11.ms has unexpected NAME column:")
                print("           ", a)
                print(" expected: ", compa)
                retValue['success']=False

        self.assertTrue(retValue['success'])

    def test12(self):
        '''Concat 12: 2 parts of same MS split in spw,  use of ephemerides'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['xya.ms','xyb.ms'],concatvis=msname, copypointing=False)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Uranus_54708-55437dUTC.tab",
                            "FIELD/EPHEM1_Jupiter_54708-55437dUTC.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        nonames = set(["FIELD/EPHEM2_Uranus_54708-55437dUTC.tab",
                       "FIELD/EPHEM3_Uranus_54708-55437dUTC.tab",
                       "FIELD/EPHEM2_Jupiter_54708-55437dUTC.tab",
                       "FIELD/EPHEM3_Jupiter_54708-55437dUTC.tab"])
        for name in nonames:
            if os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "should not exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' should not exist'
            else:
                print(myname, ": ", name, " not present as expected.")
        self.assertTrue(retValue['success'])

        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test12.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test12.ms',ignore_errors=True)
#             shutil.copytree(msname,'test12.ms')
            #print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

        self.assertTrue(retValue['success'])


    def test13(self):
        '''Concat 13: 2 parts of same MS split in time,  use of ephemerides only in second part'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['xy1-noephem.ms','xy2.ms'],concatvis=msname, copypointing=False)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Uranus_54708-55437dUTC.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test13.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test13.ms',ignore_errors=True)
#             shutil.copytree(msname,'test13.ms')
            #print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

        self.assertTrue(retValue['success'])

    def test14(self):
        '''Concat 14: 2 parts of same MS split in time,  use of ephemerides, first ephemeris covers both MS time ranges, not chronologically ordered'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['xy2.ms', 'xy1.ms'],concatvis=msname, copypointing=False)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Uranus_54708-55437dUTC.tab",
                            "FIELD/EPHEM1_Jupiter_54708-55437dUTC.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test14.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test14.ms',ignore_errors=True)
#             shutil.copytree(msname,'test14.ms')
            print(myname, ": OK. Checking tables in detail ...")
#            tb.open('test14.ms/FIELD')
            tb.open(msname+'/FIELD')
            a = list(tb.getcol('NAME'))
            tb.close()
            tb.open('xy1.ms/FIELD')
            compa = list(tb.getcol('NAME'))
            compa.append('jupiter') # jupiter should occur a second time because in xy1.ms it is an ephemeris object, in xy2 it is not
            tb.close()
            retValue['success']=True
            if not (len(a)==len(compa) and a==compa):
                print("FIELD table of test14.ms has unexpected NAME column:")
                print("           ", a)
                print(" expected: ", compa)
                retValue['success']=False

        self.assertTrue(retValue['success'])


    def test15(self):
        '''Concat 15: 2 parts of same MS split in time,  use of ephemerides, first ephemeris does not cover both MS time ranges'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['xy1.ms','xy2late.ms'],concatvis=msname, copypointing=False)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Uranus_54708-55437dUTC.tab",
                            "FIELD/EPHEM1_Jupiter_54708-55437dUTC.tab",
                            "FIELD/EPHEM2_Uranus_55437-56293dUTC.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test15.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test15.ms',ignore_errors=True)
#             shutil.copytree(msname,'test15.ms')
            print(myname, ": OK. Checking tables in detail ...")
#            tb.open('test15.ms/FIELD')
            tb.open(msname+'/FIELD')
            a = list(tb.getcol('NAME'))
            tb.close()
            tb.open('xy1.ms/FIELD')
            compa = list(tb.getcol('NAME'))
            compa.append('jupiter') # jupiter should occur a second time because in xy1.ms it is an ephemeris object, in xy2late it is not
            compa.append('1908-201') # 1908-201 should occur a second time because the ephemeris in xy1.ms does not cover the time range of xy2late 
            tb.close()
            retValue['success']=True
            if not (len(a)==len(compa) and a==compa):
                print("FIELD table of test15.ms has unexpected NAME column:")
                print("           ", a)
                print(" expected: ", compa)
                retValue['success']=False

        self.assertTrue(retValue['success'])

    def test16(self):
        '''Concat 16: 3 parts, SD data, one non-concurrent, two concurrent (CAS-5316)'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    
        concat(vis=['X39a.pm03.scan3.ms', 'X425.pm03.scan4.ms', 'X425.pm04.scan4.ms'],concatvis=msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
                            "table.f0",
                            "table.f0i",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test16.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test16.ms',ignore_errors=True)
#             shutil.copytree(msname,'test16.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

#            tb.open('test16.ms')
            tb.open(msname)
            a = tb.getcol('SCAN_NUMBER')
            tb.close()
            if not (a[0]==3 and a[59]==3 and a[60]==4 and a[len(a)-1]==4):
                print("Scan numbers not as expected. Should be == 3 up to index 59, then 4 thereafter.")
                retValue['success']=False

        self.assertTrue(retValue['success'])

    def test17(self):
        '''Concat 17: 2 completely different MSs, use of ephemerides'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['uid___A002_Xab8dc1_X95a-shrunk.ms','uid___A002_Xab8dc1_Xf13-shrunk.ms'],
               concatvis=msname, copypointing=False)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Neptune_57303.95200000.tab",
                            "FIELD/EPHEM1_Uranus_57304.04900000.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
#             if 'test17.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test17.ms',ignore_errors=True)
#             shutil.copytree(msname,'test17.ms')
            #print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

        self.assertTrue(retValue['success'])


    def test18(self):
        '''Concat 18: 3 MSs, use of ephemerides with parameter forcesingleephemfield'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['nep2-shrunk.ms','nep1-shrunk.ms','nep3-shrunk.ms'],concatvis=msname, 
               copypointing=False, forcesingleephemfield=['Neptune'])
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_Titan_57494.36900000.tab",
                            "FIELD/EPHEM1_Neptune_57494.36900000.tab",
                            "FIELD/EPHEM2_Pallas_57501.36900000.tab",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
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
                
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        
            tb.open(msname+'/FIELD/EPHEM1_Neptune_57494.36900000.tab')
            if not tb.nrows()==535:
                 print(myname, ": Error  ephemeris concatenation did not work", tablename)
                 retValue['success']=False
                 retValue['error_msgs']=retValue['error_msgs']+' ephemeris concatenation did not work'
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
            tb.close()
#             if 'test18.ms' in glob.glob("*.ms"):
#                 shutil.rmtree('test18.ms',ignore_errors=True)
#             shutil.copytree(msname,'test18.ms')
            #print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

        self.assertTrue(retValue['success'])

    def test19(self):
        '''Concat 19: two MSs with different PROCESSOR table (CAS-12937)'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis = ['X1127a.smaller.ms','Xd21.smaller.ms'], concatvis = msname)

        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                            "SOURCE/table.dat",
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
                            "SOURCE/table.f0",
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
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
            if 'test19.ms' in glob.glob("*.ms"):
                shutil.rmtree('test19.ms',ignore_errors=True)
            shutil.copytree(msname,'test19.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check processor table
            name = "PROCESSOR"
            #             col name, row number, expected value, tolerance
            expected = [
                ['FLAG_ROW',  0, False, 0],
                ['FLAG_ROW',  1, False, 0],
                ['FLAG_ROW',  2, False, 0],
                ['FLAG_ROW',  3, False, 0],
                ['FLAG_ROW',  4, False, 0],
                ['FLAG_ROW',  5, False, 0],
                ['MODE_ID',  0, 3, 0],
                ['MODE_ID',  1, 0, 0],
                ['MODE_ID',  2, 1, 0],
                ['MODE_ID',  3, 0, 0],
                ['MODE_ID',  4, 0, 0],
                ['MODE_ID',  5, 1, 0],
                ['TYPE',  0, 'RADIOMETER', 0],
                ['TYPE',  1, 'CORRELATOR', 0],
                ['TYPE',  2, 'CORRELATOR', 0],
                ['TYPE',  3, 'CORRELATOR', 0],
                ['TYPE',  4, 'RADIOMETER', 0],
                ['TYPE',  5, 'CORRELATOR', 0],
                ['TYPE_ID',  0, -1, 0],
                ['TYPE_ID',  1, -1, 0],
                ['TYPE_ID',  2, -1, 0],
                ['TYPE_ID',  3, -1, 0],
                ['TYPE_ID',  4, -1, 0],
                ['TYPE_ID',  5, -1, 0],
                ['SUB_TYPE',  0, 'SQUARE_LAW_DETECTOR', 0],
                ['SUB_TYPE',  1, 'ALMA_ACA', 0],
                ['SUB_TYPE',  2, 'ALMA_ACA', 0],
                ['SUB_TYPE',  3, 'ALMA_BASELINE', 0],
                ['SUB_TYPE',  4, 'ALMA_RADIOMETER', 0],
                ['SUB_TYPE',  5, 'ALMA_BASELINE', 0]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

            tb.open('test19.ms')
            if not (tb.getcell('PROCESSOR_ID',187503) == 5):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Reindexing of main table not correct.'
            tb.close()

        self.assertTrue(retValue['success'])

    def test20(self):
        '''Concat 20: two MSs with different antenna table in full polarization (CAS-12872)'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        concat(vis=['n08c1_swap1.ms','n08c1_swap2_swapped.ms'],concatvis=msname)
        
        print(myname, ": Now checking output ...")
        mscomponents = set(["table.dat",
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
                
        self.assertTrue(retValue['success'])
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
            tb.close()
            if 'test20.ms' in glob.glob("*.ms"):
                shutil.rmtree('test20.ms',ignore_errors=True)
            shutil.copytree(msname,'test20.ms')
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

            # check Main table
            nrows = 12
            msnrows = []
            for myms in ['n08c1_swap1.ms','n08c1_swap2_swapped.ms']:
                tb.open(myms)
                msnrows.append(tb.nrows())
                tb.close()
            tb.open('n08c1_swap2.ms')
            olddata = tb.getcol('DATA', 0, nrows)
            oldflag = tb.getcol('FLAG', 0, nrows)
            oldweight = tb.getcol('WEIGHT', 0, nrows)
            oldwtspec = tb.getcol('WEIGHT_SPECTRUM', 0, nrows)
            oldsigma = tb.getcol('SIGMA', 0, nrows)
            oldsgspec = tb.getcol('SIGMA_SPECTRUM', 0, nrows)
            olduvw = tb.getcol('UVW', 0, nrows)
            tb.close()

            ms.open('test20.ms',nomodify=False)
            ms.sort(columns=['TIME','ANTENNA1','ANTENNA2'])
            ms.done()
            tb.open('test20.ms')
            newdata = tb.getcol('DATA', msnrows[0], nrows)
            newflag = tb.getcol('FLAG', msnrows[0], nrows)
            newweight = tb.getcol('WEIGHT', msnrows[0], nrows)
            newwtspec = tb.getcol('WEIGHT_SPECTRUM', msnrows[0], nrows)
            newsigma = tb.getcol('SIGMA', msnrows[0], nrows)
            newsgspec = tb.getcol('SIGMA_SPECTRUM', msnrows[0], nrows)
            newuvw = tb.getcol('UVW', msnrows[0], nrows)
            tb.close()
            result = True
            if not (olddata == newdata).all():
                print('Mismatch in DATA column')
                result = False
            if not (oldflag == newflag).all():
                print('Mismatch in FLAG column')
                result = False
            if not (oldweight == newweight).all():
                print('Mismatch in WEIGHT column')
                result = False
            if not (oldwtspec == newwtspec).all():
                print('Mismatch in WEIGHT_SPECTRUM column')
                result = False
            if not (oldsigma == newsigma).all():
                print('Mismatch in SIGMA column')
                result = False
            if not (oldsgspec == newsgspec).all():
                print('Mismatch in SIGMA_SPECTRUM column')
                result = False
            if not (olduvw == newuvw).all():
                print('Mismatch in UVW column')
                result = False

            if not result:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of MAIN table failed'

        shutil.rmtree('test20.ms',ignore_errors=True)
        self.assertTrue(retValue['success'])


class concat_cleanup(unittest.TestCase):           
    def setUp(self):
        pass
    
    def tearDown(self):
        os.system('rm -rf *.ms')

    def testrun(self):
        '''Concat: Cleanup'''
        pass
    
def suite():
    return [test_concat,concat_cleanup]        

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

#########################################################################
# test_task_virtualconcat.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.virtualconcat.html
#
#
##########################################################################
import os
import sys
import shutil
import glob
import time
import unittest

from casatools import ctsys, calibrater, table, ms
from casatasks import virtualconcat, concat, listobs

cb = calibrater( )
tb = table( )
_ms = ms( )

myname = 'test_task_virtualconcat'

# name of input MS
datapath = ctsys.resolve('unittest/virtualconcat/')
inputmslist = ['A2256LC2_4.5s-1.ms','A2256LC2_4.5s-2.ms','part2-mod2.ms','part2.ms','part4.ms','shortpart2.ms',
               'shortpart4.ms','sim7.ms','X39a.pm03.scan3.ms','X425.pm04.scan4.ms','A2256LC2_4.5s-2b.ms','part1.ms',
               'part2-mod.ms','part3.ms','shortpart1.ms','shortpart3.ms','shortpart5.ms','sim8.ms','X425.pm03.scan4.ms']

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

class test_virtualconcat(unittest.TestCase):

    def setUp(self):
        global testmms
        res = None

        # Pick up alternative data directory to run tests on MMSs
        testmms = False
        #if 'TEST_DATADIR' in os.environ:
            #testmms = True
            #print("\nTesting on MMSs ...\n")
            #DATADIR = str(os.environ.get('TEST_DATADIR'))
            #if os.path.isdir(DATADIR):
             #   datapath = DATADIR + '/concat/input/'

        for mymsname in inputmslist:
            # Only copy if the MS doesn't exist in the working directory
            if not os.path.exists(mymsname):
                print("Copying ", mymsname)
                inpms = os.path.join(datapath,mymsname)
                print(inpms)
                shutil.copytree(inpms,mymsname)

    def tearDown(self):
        shutil.rmtree(msname,ignore_errors=True)
        if os.path.exists(self.tempname): shutil.rmtree(self.tempname)

    @classmethod
    def tearDownClass(cls):
        for inpms in inputmslist:
            if os.path.exists(inpms): shutil.rmtree(inpms)

        if os.path.exists('allparts.ms'): shutil.rmtree('allparts.ms')
        if os.path.exists('refconcatenated.ms'): shutil.rmtree('refconcatenated.ms')

        if os.path.exists('ms.txt'): os.remove('ms.txt')
        if os.path.exists('mms.txt'): os.remove('mms.txt')
        if os.path.exists('diff.txt'): os.remove('diff.txt')
        if os.path.exists('ms.txt'): os.remove('ms.txt')

    def test1(self):
        '''Virtualconcat 1: 4 parts, same sources but different spws'''
        self.tempname = self._testMethodName + '.ms'
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    
        self.res = virtualconcat(vis=['part1.ms','part2.ms','part3.ms','part4.ms'],concatvis=msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output MS ", msname)
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
        '''Virtualconcat 2: 3 parts, different sources, different spws, visweightscale=[3.,2.,1.], keepcopy=True'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        self.res = virtualconcat(vis=['part1.ms','part2-mod.ms','part3.ms'],concatvis=msname, visweightscale=[3.,2.,1.], keepcopy=True)
        self.assertEqual(self.res,None)
        
        print(myname, ": Now checking output MS ", msname)
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
            ii = 0
            for myms in ['part1.ms','part2-mod.ms','part3.ms']:
                tb.open(myms)
                msnrows.append(tb.nrows())
                oldweightbeg.append(tb.getcell('WEIGHT',0))
                oldweightend.append(tb.getcell('WEIGHT',tb.nrows()-1))
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

        self.assertTrue(retValue['success'])


    def test3(self):
        '''Virtualconcat 3: 3 parts, different sources, same spws'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    
        self.res = virtualconcat(vis=['part1.ms','part2-mod2.ms','part3.ms'],concatvis=msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output MS ", msname)
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
        '''Virtualconcat 4: five MSs with identical sources but different time/intervals on them (CSV-268)'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        self.res = virtualconcat(vis = ['shortpart1.ms', 'shortpart2.ms', 'shortpart3.ms', 'shortpart4.ms', 'shortpart5.ms'],
                          concatvis = msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
        '''Virtualconcat 5: two MSs with different state table (CAS-2601)'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        self.res = virtualconcat(vis = ['A2256LC2_4.5s-1.ms','A2256LC2_4.5s-2.ms'],
                          concatvis = msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])

    def test6(self):
        '''Virtualconcat 6: two MSs with different state table and feed table'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        self.res = virtualconcat(vis = ['A2256LC2_4.5s-1.ms','A2256LC2_4.5s-2b.ms'],
                          concatvis = msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
        '''Virtualconcat 7: two MSs with different antenna table such that baseline label reversal becomes necessary'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        self.res = virtualconcat(vis = ['sim7.ms','sim8.ms'],
                          concatvis = msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check Main table
            tb.open(self.tempname)
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
                retValue['error_msgs']=retValue['error_msgs']+'Check of table test7.ms failed'
                retValue['error_msgs']=retValue['error_msgs']+'Check of table test7.ms failed'
                
        self.assertTrue(retValue['success'])

    def test8(self):
        '''Virtualconcat 8: two MSs with different antenna tables, copypointing=False'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        os.system('rm -rf ref'+msname)
        concat(vis = ['sim7.ms','sim8.ms'],
               concatvis = "ref"+msname, copypointing=False)
        
        self.res = virtualconcat(vis = ['sim7.ms','sim8.ms'],
                          concatvis = msname, copypointing=False)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        
    
            # check Main table
            tb.open("ref"+msname)
            ant1ref = tb.getcol('ANTENNA1')
            ant2ref = tb.getcol('ANTENNA2')
            tb.close()            
            
            tb.open(self.tempname)
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

                if(ant1[i]!=ant1ref[i]):
                    print("Found disagreement in ANTENNA1 in row ", i, ": ", ant1, " ", ant1ref)
                    result = False
                    break
                    
                if(ant2[i]!=ant2ref[i]):
                    print("Found disagreement in ANTENNA2 in row ", i, ": ", ant2, " ", ant2ref)
                    result = False
                    break
                
            if result:
                print(myname, ": OK. Checking pointing table ...")

            tb.open('test8.ms/POINTING')
            pointingrows = tb.nrows()
            tb.close()
            if pointingrows>0:
                print("Pointing table should be empty!")
                result = False

            if not result:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of tables main and/or pointing failed'
                
        self.assertTrue(retValue['success'])



    def test9(self):
        '''Virtualconcat 9: 3 parts, different sources, same spws, different scratch columns: no, yes, no'''
        self.tempname = self._testMethodName + '.ms'

        global testmms
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        shutil.rmtree('part2-mod2-wscratch.ms',ignore_errors=True)
        shutil.copytree('part2-mod2.ms', 'part2-mod2-wscratch.ms', True)
        print('creating scratch columns in part2-mod2-wscratch.ms')
        if testmms:
            _ms.open('part2-mod2-wscratch.ms')
            mses = _ms.getreferencedtables()
            _ms.close()
            mses.sort()
            for mname in mses:
                cb.open(mname)
                cb.close()
        else:
            cb.open('part2-mod2-wscratch.ms') # calibrator-open creates scratch columns
            cb.close()

        self.res = virtualconcat(vis=['part1.ms','part2-mod2-wscratch.ms','part3.ms'],concatvis=msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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
        
    def test10(self):
        '''Virtualconcat 10: 3 parts, different sources, same spws, different scratch columns: yes, no, no'''
        self.tempname = self._testMethodName + '.ms'

        global testmms
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        shutil.rmtree('part1-wscratch.ms',ignore_errors=True)
        shutil.copytree('part1.ms', 'part1-wscratch.ms', True)
        print('creating scratch columns in part1-wscratch.ms')
        if testmms:
            _ms.open('part1-wscratch.ms')
            mses = _ms.getreferencedtables()
            _ms.close()
            mses.sort()
            for mname in mses:
                cb.open(mname)
                cb.close()
        else:
            cb.open('part1-wscratch.ms') # calibrator-open creates scratch columns
            cb.close()

        self.res = virtualconcat(vis=['part1-wscratch.ms','part2-mod2.ms','part3.ms'],concatvis=msname)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
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

    def test11(self):
        '''Virtualconcat 11: comparison to concat'''
        self.tempname = self._testMethodName + '.ms'

        global testmms
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        if testmms:
            print("Skipping this test as concat will not work with an MMS.")
        else:
            shutil.rmtree('allparts.ms', ignore_errors=True)
            shutil.rmtree('allparts.mms', ignore_errors=True)
            os.system('rm -f ms.txt mms.txt')
        
            thebeginning = time.time()
            concat(vis=['part1.ms','part2.ms','part3.ms','part4.ms'], concatvis='allparts.ms')
            theend = time.time()
            print("duration using concat (s) = ", theend-thebeginning)

            thebeginning = time.time()
            virtualconcat(vis=['part1.ms','part2.ms','part3.ms','part4.ms'], concatvis='allparts.mms')
            theend = time.time()
            print("duration using virtualconcat (s) =", theend-thebeginning)
        
            listobs(vis='allparts.ms', listfile='ms.txt')
            shutil.rmtree('allparts.ms')
            shutil.move('allparts.mms', 'allparts.ms') # to get same file name
            listobs(vis='allparts.ms', listfile='mms.txt')
            os.system('diff ms.txt mms.txt > diff.txt')
            os.system('cat diff.txt')
            retValue['success'] = (os.path.getsize('diff.txt') == 0)

        self.assertTrue(retValue['success'])

    def test12(self):
        '''Virtualconcat 12: two MSs with different antenna tables, copypointing=True (default)'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        
        self.res = virtualconcat(vis = ['sim7.ms','sim8.ms'],
                          concatvis = msname, copypointing=True)
        self.assertEqual(self.res,None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            _ms.close()
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname, True)
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True        

            result = True
            tb.open('test12.ms/POINTING')
            pointingrows = tb.nrows()
            tb.close()
            if pointingrows==0:
                result = False

            if not result:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of pointing table failed'
                
        self.assertTrue(retValue['success'])

    def test13(self):
        '''Virtualconcat 13: 3 parts, SD data, one non-concurrent, two concurrent (CAS-5316)'''
        self.tempname = self._testMethodName + '.ms'

        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    
        self.res = virtualconcat(vis=['X39a.pm03.scan3.ms', 'X425.pm03.scan4.ms', 'X425.pm04.scan4.ms'],concatvis=msname)
        self.assertEqual(self.res, None)

        print(myname, ": Now checking output ...")
        try:
            _ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            _ms.close()
            if self.tempname in glob.glob("*.ms"):
                shutil.rmtree(self.tempname,ignore_errors=True)
            shutil.copytree(msname,self.tempname)
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True

            tb.open(self.tempname)
            a = tb.getcol('SCAN_NUMBER')
            tb.close()
            if not (a[0]==3 and a[59]==3 and a[60]==4 and a[len(a)-1]==4):
                print("Scan numbers not as expected. Should be == 3 up to index 59, then 4 thereafter.")
                retValue['success']=False

        self.assertTrue(retValue['success'])

if __name__ == '__main__':
    unittest.main()

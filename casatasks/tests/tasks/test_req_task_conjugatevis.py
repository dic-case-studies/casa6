########################################################################
# test_task_conjugatevis.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.calibration.conjugatevis.html
#
#
##########################################################################
import sys
import os
import unittest
import shutil
import glob
import numpy
import time

import casatools
from casatasks import conjugatevis, casalog, listobs
conjtb = casatools.table()
origtb = casatools.table()
tb = casatools.table()
ms = casatools.ms()

# Define paths to sample data files used for tests.
datapath = casatools.ctsys.resolve('unittest/conjugatevis/')

logpath = casalog.logfile()
msfile = 'gaincaltest2.ms'

### Funtions from merged test test_conjugatevis ###
myname = 'test_conjugatevis'

# name of the resulting MS
msname = 'conjugated.ms'

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:
    testmms = True
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/concat/input/'
    if os.path.isdir(DATADIR):
        datapath = DATADIR
    print('conjugatevis tests will use data from %s' % datapath)


def checktable(thename, theexpectation):
    global msname, myname
    tb.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print("%s: comparing %s"% (myname,mycell))
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
            print("%s:  Error in MS subtable %s:" % (myname,thename))
            print ("     column %s row %s contains %s" % (mycell[0],mycell[1],value))
            print("     expected value is %s" % mycell[2])
            tb.close()
            return False
    tb.close()
    print("%s: table %s as expected." %  (myname, thename))
    return True

class conjugatevis_test(unittest.TestCase):

    def setUp(self):
        if not os.path.exists(msfile):
            shutil.copytree(os.path.join(datapath, msfile),msfile)
        res = None

        cpath = os.path.abspath(os.curdir)
        filespresent = sorted(glob.glob("*.ms"))
        os.chdir(datapath)
        mymsname = 'shortpart1.ms'
        if not mymsname in filespresent:
            print("Copying %s" % mymsname)
            shutil.copytree(mymsname, cpath+'/'+mymsname)
            os.chdir(cpath)

    # Remove files created during tests
    def tearDown(self):
        casalog.setlogfile(logpath)
        shutil.rmtree(msfile, ignore_errors=True)
        shutil.rmtree(msname, ignore_errors=True)
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
        if os.path.exists('gaincal2-conj.ms'):
            shutil.rmtree('gaincal2-conj.ms', ignore_errors=True)

        if os.path.exists('conjugated_copygaincal2.ms'):
            shutil.rmtree('conjugated_copygaincal2.ms', ignore_errors=True)

        if os.path.exists('copygaincal2.ms'):
            shutil.rmtree('copygaincal2.ms', ignore_errors=True)

        if os.path.exists('shortpart1.ms'):
            shutil.rmtree('shortpart1.ms')

    def test_takesMeasurementSet(self):
        ''' 1. test_takesMeasurementSet: Check that conjugatevis opens a MeasurementSet file'''
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)
        # Check that conjugatevis created an output file.
        self.assertTrue('gaincal2-conj.ms' in os.listdir('.'))

    def test_changeSignOneSpectralWindow(self):
        ''' 2. test_changeSignOnSpectralWindow: Check that conjugatevis changed the phase sign for a single spectral window'''

        # Open the original MS file and get the column of phase information for the data
        tb.open(msfile)
        originalData = tb.getcol('DATA')
        originalPhase = numpy.imag(originalData[:,:,0][0][:])
        tb.close()

        conjugatevis(vis=msfile, spwlist = 0, outputvis='gaincal2-conj.ms', overwrite=True)

        # Open the conjugated version of the data to obtain the phase information for the data
        tb.open('gaincal2-conj.ms')
        conjData = tb.getcol('DATA')
        conjPhase = numpy.imag(conjData[:,:,0][0][:])
        tb.close()

        # The conjugated phase values should be -1 * the original phase values
        self.assertTrue((conjPhase == -1 * originalPhase).all())

    def test_changeMultipleSpectralWindows(self):
        ''' 3. test_changeMultipleSpectralWindows: Check that conjugatevis changed the phase sign for multiple spectral windows'''

        # Open the original MS file to get phase information for the data
        tb.open(msfile)
        originalData = tb.getcol('DATA')
        originalPhase = numpy.imag(originalData[:,:,0][0][:])
        tb.close()

        # Open the conjugated version of the data to obtain phase information for data '''
        conjugatevis(vis=msfile, spwlist = [0,1], outputvis='gaincal2-conj.ms', overwrite=True)
        tb.open('gaincal2-conj.ms')
        conjData = tb.getcol('DATA')
        conjPhase = numpy.imag(conjData[:,:,0][0][:])
        tb.close()

        # The conjugated version of the data should be -1 * the original phase values
        self.assertTrue((conjPhase == -1 * originalPhase).all())

    def test_defaultSpectralWindowList(self):
        ''' 4. test_defaultSpectralWindows: Check that conjugatevis changes phase sign for all spectral windows for default spwlist'''
        
        # Open the original MS file to get the phase information
        tb.open(msfile)
        originalData = tb.getcol('DATA')
        originalPhase = numpy.imag(originalData[:,:,0][0][:])
        tb.close()
        
        # Open the conjugated version of the data to get phase information
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)
        tb.open('gaincal2-conj.ms')
        conjData = tb.getcol('DATA')
        conjPhase = numpy.imag(conjData[:,:,0][0][:])
        tb.close()

        # The conjugated version of the data should be -1 * the original phase values
        self.assertTrue((conjPhase == -1 * originalPhase).all())

    def test_specifiedOutputFileName(self):
        ''' 5. test_specifiedOutputFileName: Check that conjugatevis writes the output ms file with the specified file name'''
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)
        self.assertTrue('gaincal2-conj.ms' in os.listdir('.'))

    def test_defaultOutputFileName(self):
        ''' 6. test_defaultOutputFileName: Check that conjugatevis writes the output ms file as 'conjugated_'+vis'''
        
        # Copy the original file to the working directory '''
        shutil.copytree(msfile, 'copygaincal2.ms')

        # Run conjugatevis on the copied file, which is necessary because otherwise it will attempt to name the output file with the full path of the sample file
        conjugatevis(vis='copygaincal2.ms')

        # The default output file should be conjugated_copygaincal2.ms
        newfilename = 'conjugated_copygaincal2.ms'

        # Copied file is no longer needed 
        shutil.rmtree('copygaincal2.ms')
        self.assertTrue(os.path.exists(newfilename))

    def test_overwriteTrue(self):
        ''' 7. test_overwriteTrue: Check that conjugatevis will overwrite an existing file if overwrite=True'''

        # Run conjugatevis to make sure gaincal2-conj.ms exists
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)

        # Get modification time for this file
        modificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))
        
        # Run conjugatevis again to overwrite the file
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)

        # Get modifcation time for new version
        afterModificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))

        # The modifications times should be different since gaincal2-conj.ms was overwritten.
        self.assertFalse(modificationTime == afterModificationTime)

    #TODO Needs work
    def test_overwriteFalse(self):
        ''' 8. test_overwriteFalse: Check that conjugatevis will not overwrite an existing file if overwrite=False'''
        
        # Run conjugatevis to make sure gaincal2-conj.ms 
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)
        modificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))
        
        try:
            conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=False)
            self.fail()
        except Exception:
            self.assertTrue(True)

    #TODO Needs work
    def test_overwriteDefault(self):
        ''' 9. test_overwriteDefault: Check that the default setting of overwrite is False '''

        # Run conjugatevis to make sure gaincal2-conj.ms exists
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)
        modificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))

        # Run again and expect an exception to be raised
        try:
            conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=False)
            self.fail()
        except Exception:
            self.assertTrue(True)

    # MERGED TEST CASES #
    def test1(self):
        '''Conjugatevis 1: '''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        self.res = conjugatevis(vis='shortpart1.ms', spwlist=[5,7], outputvis=msname)
        self.assertEqual(self.res,None)

        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["table.dat",
                            "table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,msname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. Try opening as MS ..." %  myname)
        try:
            ms.open(msname)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,tablename))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            retValue['success']=True

            # check main table
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                ['DATA',           4000,
                 [[-0.00426177-0.00387163j],
                  [ 0.00058119+0.00283016j]],
                 0.00000001]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            expected = [
                ['DATA',           3000,
                 [[ 0.00347826-0.00406267j],
                  [ 0.00458098-0.00508398j]],
                 0.00000001]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'

# Main #
if __name__ == '__main__':
    unittest.main()


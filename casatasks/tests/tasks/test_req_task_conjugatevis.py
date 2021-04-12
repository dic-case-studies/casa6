########################################################################
# test_req_task_conjugatevis.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_conjugatevis/about
#
#
##########################################################################
''' This function is run using the command from the bin directory of the 
casa prerelease directory to be tested:
./casa -c ../lib/python2.7/runUnitTest.py test_req_task_conjugatevis '''

CASA6 = False
try:
    import casatools
    from casatasks import conjugatevis, casalog, listobs
    CASA6 = True
    conjtb = casatools.table()
    origtb = casatools.table()
    tb = casatools.table()

except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from taskinit import tbtool

    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)

import sys
import os
import unittest
import shutil
import numpy
import time

# Define paths to sample data files used for tests.
if CASA6:
    datapath = casatools.ctsys.resolve('unittest/conjugatevis/')

else:
#    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
#        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/gaincaltest2.ms'
#    else:
    datapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/conjugatevis/'

logpath = casalog.logfile()
msfile = 'gaincaltest2.ms'

class conjugatevis_test(unittest.TestCase):

    def setUp(self):
        if not os.path.exists(msfile):
            shutil.copytree(os.path.join(datapath, msfile),msfile)
        if not CASA6:
            default(conjugatevis)

    # Remove files created during tests
    def tearDown(self):
        casalog.setlogfile(logpath)
        shutil.rmtree(msfile, ignore_errors=True)
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
        if os.path.exists('gaincal2-conj.ms'):
            shutil.rmtree('gaincal2-conj.ms')

        if os.path.exists('conjugated_copygaincal2.ms'):
            shutil.rmtree('conjugated_copygaincal2.ms')

        if os.path.exists('copygaincal2.ms'):
            shutil.rmtree('copygaincal2.ms')

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
        
        if CASA6 or casa_stack_rethrow:
            # Run again and expect an exception to be raised
            try:
                conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=False)
                self.fail()
            except Exception:
                self.assertTrue(True)
        else:
            # Run conjugatevis without overwriting
            conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=False)
            afterModificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))

            # Since gaincal2-conj.ms was not overwritten, the modification times should be the same.
            self.assertTrue(modificationTime == afterModificationTime)

    #TODO Needs work
    def test_overwriteDefault(self):
        ''' 9. test_overwriteDefault: Check that the default setting of overwrite is False '''

        # Run conjugatevis to make sure gaincal2-conj.ms exists
        conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=True)
        modificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))

        if CASA6 or casa_stack_rethrow:
            # Run again and expect an exception to be raised
            try:
                conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=False)
                self.fail()
            except Exception:
                self.assertTrue(True)
        else:
            # Run conjugatevis without overwriting
            conjugatevis(vis=msfile, outputvis='gaincal2-conj.ms', overwrite=False)
            afterModificationTime = time.ctime(os.path.getmtime('gaincal2-conj.ms'))

            # Since gaincal2-conj.ms was not overwritten, the modification times should be the same.
            self.assertTrue(modificationTime == afterModificationTime)

def suite():
    return[conjugatevis_test]

# Main #
if __name__ == '__main__':
    unittest.main()


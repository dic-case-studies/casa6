########################################################################
# test_req_task_imhistory.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_imhistory/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import imhistory, casalog
    import sys
    import os
    myia = casatools.image()
    tb = casatools.table()
    sys.path.append(os.path.abspath(os.path.dirname(__file__)))
    CASA6 = True
except ImportError:
    import sys
    import os
    from __main__ import default
    from tasks import *
    from taskinit import *
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
    myia = iatool()

import unittest
import shutil

if CASA6:
    casaimagepath = casatools.ctsys.resolve('image/ngc5921.clean.image')
    fitspath = casatools.ctsys.resolve('fits/1904-66_AIR.fits')
    #miriadpath = casatools.ctsys.resolve('visibilities/other/compact.vis')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        casaimagepath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image'
        fitspath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/fits/1904-66_AIR.fits'
    else:
        casaimagepath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image'
        fitspath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/fits/1904-66_AIR.fits'
        
def change_perms(path):
    os.chmod(path, 0o777)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root,d), 0o777)
        for f in files:
            os.chmod(os.path.join(root,f), 0o777)
        
logpath = casalog.logfile()
imagecopy = 'imagecopy.image'

class imhistory_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        shutil.copytree(casaimagepath, imagecopy)
        change_perms(imagecopy)
    
    def setUp(self):
        if not CASA6:
            default(imhistory)

    def tearDown(self):
        myia.done()
        tb.done()
        casalog.setlogfile(logpath)
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
        if os.path.exists('basic'):
            shutil.rmtree('basic')
            
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(imagecopy)

    def test_takesCASAImage(self):
        ''' 1. test_takesCASAImage: Check that imhistory takes a CASA image file (*.image)'''
        messages = imhistory(imagecopy, mode='list', verbose=False)
        self.assertTrue(messages)

    def test_takesFITS(self):
        ''' 2. test_takesFITS: Check that imhistory takes a FITS file '''
        messages = imhistory(fitspath, mode='list', verbose=False)
        self.assertTrue(messages)

    def test_listModeVerbose(self):
        ''' 3. test_listModeVerbose: Check that the list mode with verbose on outputs to log file and outputs an array of strings '''
        casalog.setlogfile('testlog.log')
        historyMessages = imhistory(imagecopy, mode='list', verbose=True)
        self.assertTrue(len(historyMessages) > 0 and 'HISTORY' in open('testlog.log').read())

    def test_listModeNoVerbose(self):
        ''' 4. test_listModeNoVerbose: Check that the list mode with verbose off outputs an array of strings and does not output to the log file '''
        casalog.setlogfile('testlog.log')
        historyMessages = imhistory(imagecopy, mode='list', verbose=False)
        
        if CASA6:
            self.assertFalse(os.path.getsize("testlog.log") > 41336)
        else:
            self.assertFalse('HISTORY' in open('testlog.log').read())

    def test_appendModeNoDefaults(self):
        '''5. test_appendModeNoDefaults: Check that the append mode adds a string to the image history without use of default settings for message or origin '''
        casalog.setlogfile('testlog.log')
        success = imhistory(imagecopy, mode='append', message='TESTMESSAGEtest5', origin='TESTORIGINtest5')
        # Run imhistory again to output the history messages to the log to check if the message was added
        imhistory(imagecopy, mode='list', verbose=True)
        self.assertTrue('TESTMESSAGEtest5' in open('testlog.log').read() and 'TESTORIGINtest5' in open('testlog.log').read()) 

    def test_appendModeDefaultOrigin(self):
        ''' 6. test_appendModeDefaultOrigin: Check that append mode adds a string to the image history with the default origin setting '''
        casalog.setlogfile('testlog.log')
        #default(imhistory)
        success = imhistory(imagecopy, mode='append', message='TESTMESSAGEtest6')
        # Run imhistory again to output the history messages to the log to check if the message was added.
        if not CASA6:
            default(imhistory)
        imhistory(imagecopy, mode='list', verbose=True)
        self.assertTrue('imhistory' in open('testlog.log').read() and 'TESTMESSAGEtest6' in open('testlog.log').read())

    def test_correctReturnedParameters(self):
        ''' 7. test_correctReturnedParameters: Check that imhistory returns the expected parameters by looking for FILLM and BPASS '''
        casalog.setlogfile('testlog.log')
        historyMessages = imhistory(imagecopy, mode='list')
        self.assertTrue(('FILLM' in s for s in historyMessages) and ('BPASS' in n for n in historyMessages))
        
    def test_noExistingMode(self):
        ''' 8. test_noExistingMode: Check that an exception is raised when a non-valid mode is given '''
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(Exception):
                imhistory(imagecopy, mode='fakeMode')
        else:
            casalog.setlogfile('testlog.log')
            imhistory(imagecopy, mode='fakemode')
            self.assertTrue('SEVERE' in open('testlog.log').read())
            
            
    # merged imhistory tests start here
    # ---------------------------------------------
            
    def test_imhistory(self):
        """Test general functionality"""
        shape = [2,3,4]
        imagename = "basic"
        myia.fromshape(imagename, shape)
        myia.done()
        h = imhistory(imagename, mode="list")
        self.assertTrue(len(h) == 3, "Incorrect history length")
        for hh in h[1:2]:
            print(hh)
            self.assertTrue("fromshape" in hh, "Incorrect message")
        msg = "fred"
        self.assertTrue(
            imhistory(imagename, mode="append", message=msg),
            "Error appending message"
        )
        h = imhistory(imagename, mode="list")
        self.assertTrue(len(h) == 4, "Incorrect history length")
        for hh in h[1:2]:
            self.assertTrue("fromshape" in hh, "Incorrect message")
        self.assertTrue(msg in h[3], "Incorrect appended message")
        

def suite():
    return[imhistory_test]

# Main #
if __name__ == '__main__':
    unittest.main()



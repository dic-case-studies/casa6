##########################################################################
# test_task_listfits.py
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
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
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.information.listfits.html
#
##########################################################################
import os
import sys
import shutil
import unittest

from casatools import ctsys
from casatasks import casalog, listfits

class listfits_test(unittest.TestCase):
    
    # Input and output names
    fitsfile = 'ngc5921.fits'
    res = None

    def setUp(self):
        self.res = None
        
        if(os.path.exists(self.fitsfile)):
            os.system('rm -rf ' + self.fitsfile)

        shutil.copyfile(ctsys.resolve(os.path.join('unittest/listfits/',self.fitsfile)), self.fitsfile)
    
    def tearDown(self):
        if (os.path.exists(self.fitsfile)):
            os.system('rm -rf ' + self.fitsfile)
        
    def test1(self):
        '''Test 1: Empty input should return False'''
        # CASA5 tasks return False, casatasks throw exceptions
        fitsfile = ''
        self.assertRaises(Exception, listfits, fitsfile)
        
    def test2(self):
        '''Test 2: Good input should return None'''
        self.res=listfits(self.fitsfile)
        self.assertEqual(self.res,None)
        
    def test3(self):
        '''Test 3: list the fits into a private logfile'''
        thelogfile=casalog.logfile()
        logfile= "mylistfits.log"
        try:
           open(logfile,"w").close()
           casalog.setlogfile(logfile)
           listfits(self.fitsfile)
           casalog.setlogfile(thelogfile)
        except:
           print('could not open "%s" for writing' % logfile)
        
if __name__ == '__main__':
    unittest.main()

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys
    from casatasks import casalog, listfits
else:
    import commands
    from __main__ import default
    from tasks import *
    from taskinit import *

class listfits_test(unittest.TestCase):
    
    # Input and output names
    fitsfile = 'ngc5921.fits'
    res = None

    def setUp(self):
        self.res = None
        
        if(os.path.exists(self.fitsfile)):
            os.system('rm -rf ' + self.fitsfile)

        if is_CASA6:
            shutil.copyfile(ctsys.resolve(os.path.join('regression/ngc5921',self.fitsfile)), self.fitsfile)
        else:
            default(listfits)
            shutil.copyfile(os.environ.get('CASAPATH').split()[0] +\
                            '/data/regression/ngc5921/'+self.fitsfile, self.fitsfile)
    
    def tearDown(self):
        if (os.path.exists(self.fitsfile)):
            os.system('rm -rf ' + self.fitsfile)
        
    def test1(self):
        '''Test 1: Empty input should return False'''
        # CASA5 tasks return False, casatasks throw exceptions
        fitsfile = ''
        if is_CASA6:
            self.assertRaises(Exception, listfits, fitsfile)
        else:
            self.res = listfits(fitsfile)
            self.assertFalse(self.res)
        
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
        
def suite():
    return [listfits_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

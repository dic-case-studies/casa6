from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    ### for testhelper import
    import sys
    sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
    import testhelper as th

    from casatools import ctsys
    from casatasks import gaincal

    # default is not used in CASA6
    def default(atask):
        pass
    
    datapath = ctsys.resolve('regression/unittest/gaincal')
else:
    import testhelper as th
    from __main__ import default
    from tasks import gaincal
    from taskinit import *

    datapath = os.environ.get('CASAPATH').split()[0] +\
        '/data/regression/unittest/gaincal/'

''' Python unit tests for the gaincal task

These tests will only verify if the gain calibration
tables created for an MS and an MMS agree. These are
not full unit tests for the gaincal task.
'''

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:   
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/gaincal/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print('WARN: directory '+DATADIR+' does not exist')

print('gaincal tests will use data from %s' % datapath)


# Base class which defines setUp functions
# for importing different data sets
class test_base(unittest.TestCase):
    
    def cleanUp(self):
        shutil.rmtree(self.msfile, ignore_errors=True)
        os.system('rm -rf '+self.msfile+'.gcal')
    
    def setUp_ngc5921(self):
        
        # Input names
        prefix = 'ngc5921'
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'
            
        self.reffile = os.path.join(datapath,prefix)
        self.cleanUp()
        
        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):        
            shutil.copytree(fpath, self.msfile)
        else:
            self.fail('Data does not exist -> '+fpath)

        default('gaincal')
               
        
    def setUp_ngc4826(self):
        
        # Input names
        prefix = 'ngc4826'
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'
            
        self.reffile = os.path.join(datapath,prefix)
        self.cleanUp()

        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):
            shutil.copytree(fpath, self.msfile)
        else:
            self.fail('Data does not exist -> '+fpath)

        default('gaincal')


class gaincal1_test(test_base):

    def setUp(self):
        self.setUp_ngc5921()

    def tearDown(self):
        if os.path.lexists(self.msfile):
            shutil.rmtree(self.msfile)
        
        os.system('rm -rf ngc5921*.gcal')
        
    def test1a(self):
        '''Gaincal 1a: Default values to create a gain table'''
        msgcal = self.msfile + '.gcal'
        reference = self.reffile + '.ref1a.gcal'
        gaincal(vis=self.msfile, caltable=msgcal,uvrange='>0.0')
        self.assertTrue(os.path.exists(msgcal))
        
        # Compare the calibration table with a reference
        self.assertTrue(th.compTables(msgcal, reference, ['WEIGHT']))


    def test2a(self):
        '''Gaincal 2a: Create a gain table using field selection'''
        
        msgcal = self.msfile + '.field0.gcal'
        reference = self.reffile + '.ref2a.gcal'
        gaincal(vis=self.msfile, caltable=msgcal, 
                uvrange='>0.0',field='0', gaintype='G',solint='int',
                combine='',refant='VA02')
        self.assertTrue(os.path.exists(msgcal))

        # Compare the calibration tables
        self.assertTrue(th.compTables(msgcal, reference, ['WEIGHT']))        


class gaincal2_test(test_base):

    def setUp(self):
        self.setUp_ngc4826()
           
            
    def tearDown(self):
        if os.path.lexists(self.msfile):
            shutil.rmtree(self.msfile)

        os.system('rm -rf ngc4826*.gcal')
        
        
    def test1b(self):
        '''Gaincal 1b: Create a gain table for an MS with many spws'''
        msgcal = self.msfile + '.gcal'
        reference = self.reffile + '.ref1b.gcal'
        gaincal(vis=self.msfile, caltable=msgcal, 
                uvrange='>0.0',field='0,1',spw='0', gaintype='G',minsnr=2.0,
                refant='ANT5', solint='inf',combine='')
        self.assertTrue(os.path.exists(msgcal))

        # Compare the calibration tables
        self.assertTrue(th.compTables(msgcal, reference, ['WEIGHT']))


def suite():
    return [gaincal1_test, gaincal2_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

from __future__ import absolute_import
from __future__ import print_function
import os
import shutil

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    ### for testhelper import
    #import sys
    #sys.path.append(os.path.abspath(os.path.dirname(__file__)))
    #import testhelper as th
    from casatools import ctsys
    from casatasks import bandpass
else:
    #import testhelper as th
    from __main__ import default
    from tasks import bandpass
    from taskinit import *

from casatestutils import testhelper as th
import unittest


''' Python unit tests for the bandpass task

These tests will only verify if the bandpass calibration
tables created for an MS and an MMS agree. These are
not full unit tests for the bandpass task.
'''

if is_CASA6:
    datapath = ctsys.resolve('unittest/bandpass/')
else:
    datapath = os.environ.get('CASAPATH').split()[0] +\
               '/casatestdata/unittest/bandpass/'

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/bandpass/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print('WARN: directory '+DATADIR+' does not exist')

print('bandpass tests will use data from '+datapath)

# Base class which defines setUp functions
# for importing different data sets
class test_base(unittest.TestCase):

    def cleanUp(self):
        shutil.rmtree(self.msfile, ignore_errors=True)
        os.system('rm -rf '+self.msfile+'.bcal')
    
    def setUp_ngc5921(self):
        
        # Input names
        prefix = 'ngc5921'
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'

        self.reffile = os.path.join(datapath, prefix)
        self.cleanUp()
        
        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):        
            shutil.copytree(fpath, self.msfile)
        else:
            self.fail('Data does not exist -> '+fpath)

        if not is_CASA6:
            default('bandpass')
               
        
    def setUp_ngc4826(self):
        
        # Input names
        prefix = 'ngc4826'
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'

        self.reffile = os.path.join(datapath, prefix)
        self.cleanUp()
        
        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):
            shutil.copytree(fpath, self.msfile)
        else:
            self.fail('Data does not exist -> '+fpath)

        if not is_CASA6:
            default('bandpass')


class bandpass1_test(test_base):

    def setUp(self):
        self.setUp_ngc5921()

    def tearDown(self):
        if os.path.lexists(self.msfile):
            shutil.rmtree(self.msfile)

        os.system('rm -rf ngc5921*.bcal')
        
    def test1a(self):
        '''Bandpass 1a: Create bandpass table using field=0'''
        msbcal = self.msfile + '.bcal'
        reference = self.reffile + '.ref1a.bcal'
        bandpass(vis=self.msfile, caltable=msbcal, field='0',uvrange='>0.0',
                 bandtype='B',solint='inf',combine='scan',refant='VA15')
        self.assertTrue(os.path.exists(msbcal))

        # Compare the calibration tables
        self.assertTrue(th.compTables(msbcal, reference, ['WEIGHT']))


class bandpass2_test(test_base):

    def setUp(self):
        self.setUp_ngc4826()
                       
    def tearDown(self):
        if os.path.lexists(self.msfile):
            shutil.rmtree(self.msfile)


        os.system('rm -rf ngc4826*.bcal')
        
        
    def test1b(self):
        '''Bandpass 1b: Create cal tables for the MS and MMS split by spws'''
        ### Forced Failure
        msbcal = self.msfile + '.bcal'
        reference = self.reffile + '.ref1b.bcal'
        #bandpass(vis=self.msfile, caltable=msbcal, uvrange='>0.0',
        #         field='0',spw='0',bandtype='B',
        #         solint='inf',combine='scan',refant='ANT5')
        bandpass(vis=self.msfile, caltable=msbcal, uvrange='>0.0',
                 field='0',spw='1',bandtype='B',
                 solint='inf',combine='scan',refant='ANT5')
        self.assertTrue(os.path.exists(msbcal))

        # Compare the calibration tables
        self.assertTrue(th.compTables(msbcal, reference, ['WEIGHT']))


def suite():
    return [bandpass1_test, bandpass2_test]
        
if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

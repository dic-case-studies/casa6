import os
import sys
import shutil
import unittest
from CASAtools import ctsys, image, table
from CASAtasks import clearstat

_ia = image( )
_tb = table( )

'''
Unit tests of task clearstat. It tests the following parameters:
    clears read lock on table,
    clears write lock on table,
    clears read lock on image,
    clears write lock on image,
    clears all locks
'''
datapath = 'regression/exportasdm/input'

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:
    testmms = True
    DATADIR = os.path.join(str(os.environ.get('TEST_DATADIR')),'clearstat')
    if os.path.isdir(DATADIR):
        datapath = DATADIR
    print('clearstat tests will use data from %s' % datapath)

class clearstat_test(unittest.TestCase):
    
    # Input names
    msfile = 'Itziar.ms'
    res = None
    img = 'n4826_tmom1.im'
    
    def setUp(self):
        self.res = None
        if(os.path.exists(self.msfile)):
            os.system('rm -rf ' + self.msfile)
        if(os.path.exists(self.img)):
            os.system('rm -rf ' + self.img)
            
        shutil.copytree(ctsys.resolve(os.path.join(datapath,self.msfile)), self.msfile)
        shutil.copytree(ctsys.resolve(os.path.join('regression/ngc4826redux/reference', self.img)), self.img)
    
    def tearDown(self):
        os.system('rm -rf ' + self.msfile)
        os.system('rm -rf ' + self.img)

        _tb.close()
        if(_ia.isopen == True):
            _ia.close()
            
        
    def test1(self):
        '''Test 1: Clear table read lock'''
        _tb.open(self.msfile)
        lock = _tb.haslock(write=False)
        self.assertTrue(lock,'Cannot acquire read lock on table')
        clearstat( )
        lock = _tb.haslock(write=False)
        _tb.close( )
        self.assertFalse(lock,'Failed to clear table read lock')

    def test2(self):
        '''Test 2: Clear table write lock'''
        _tb.open(self.msfile)
        _tb.lock()
        lock = _tb.haslock(write=True)
        self.assertTrue(lock,'Cannot acquire write lock on table')
        clearstat( )
        lock = _tb.haslock(write=True)
        _tb.close()
        self.assertFalse(lock,'Failed to clear table write lock')

    def test3(self):
        '''Test 3: Clear image read lock'''
        _ia.open(self.img)
        lock = _ia.haslock()
        self.assertTrue(lock[0]==True and lock[1]==False,'Cannot acquire read lock on image')
        clearstat()
        lock = _ia.haslock()
        _ia.close()
        self.assertTrue(lock[0]==False and lock[1]==False,'Failed to clear read lock on image')

    def test4(self):
        '''Test 4: Clear image write lock'''
        _ia.open(self.img)
        _ia.lock(writelock=True)
        lock = _ia.haslock()
        self.assertTrue(lock[0]==True and lock[1]==True,'Cannot acquire write lock on image')
        clearstat()
        lock = _ia.haslock()
        _ia.close()
        self.assertTrue(lock[0]==False and lock[1]==False,'Failed to clear write lock on image')

    def test5(self):
        '''Test 5: Clear all locks'''
        _tb.open(self.msfile)
        tbreadlock = _tb.haslock(write=False)
        _tb.lock()
        tbwritelock = _tb.haslock(write=True)
        _ia.open(self.img)
        _ia.lock(writelock=True)
        lock = _ia.haslock()
        self.assertTrue(tbreadlock==True and tbwritelock==True and lock[0]==True and lock[1]==True,
                        'Cannot acquire locks on table and/or image')
        clearstat()
        tbreadlock = _tb.haslock(write=False)
        tbwritelock = _tb.haslock(write=True)
        lock = _ia.haslock()
        _tb.close()
        _ia.close()

        self.assertTrue(tbreadlock==False and tbwritelock==False and lock[0]==False and lock[1]==False,
                        'Failed to clear locks on table and/or image')

def suite():
    return [clearstat_test]

if __name__ == '__main__':
    unittest.main()

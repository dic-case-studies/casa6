import os
import io
import sys
import shutil
import subprocess
import unittest

from casatools import ctsys
from casatools.platform import bytes2str
from casatasks import casalog, listhistory

datapath = ctsys.resolve('regression/unittest/listhistory')

testmms = False
if 'TEST_DATADIR' in os.environ:
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/listhistory/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print('WARN: directory %s does not exist' % DATADIR)

print('listhistory tests will use data from %s' % datapath)

class listhistory_test(unittest.TestCase):

    # Input and output names
    msfile = 'Itziar.ms'
    itismms = testmms

    def setUp(self):
        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):
            os.symlink(fpath, self.msfile)
        else:
            self.fail('Data does not exist -> '+fpath)
            
    def tearDown(self):
        if os.path.lexists(self.msfile):
            os.unlink(self.msfile)
        
    def test1(self):
        '''Test 1: Empty input should return False'''
        myms = ''
        self.assertRaises(Exception,listhistory,myms)
        
    def test2(self):
        '''Test 2: Good input should return None'''
        res = listhistory(self.msfile)
        self.assertEqual(res,None)
        
    def test3(self):
        '''Test 3: Compare length of reference and new lists'''
        logfile= "mylisth.log"

        open(logfile,"w").close( )
        casalog.setlogfile(logfile)
        res = listhistory(self.msfile)

        # Get the number of lines in file
        refnum=10
        if self.itismms:
            refnum = 36

        cmd=['wc', '-l', logfile]
        print(cmd)
        output=bytes2str(subprocess.check_output(cmd))
        num = int(output.split()[0])
        self.assertEqual(refnum,num)


class listhistory_cleanup(unittest.TestCase):
    
    def tearDown(self):
        os.system('rm -rf *Itziar.*')

    def test_cleanup(self):
        '''listhistory: Cleanup'''
        pass
        
def suite():
    return [listhistory_test, listhistory_cleanup]

if __name__ == '__main__':
    unittest.main()

from __future__ import absolute_import
from __future__ import print_function
import os
import io
import sys
import shutil
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys
    from casatools.platform import bytes2str
    from casatasks import casalog, listhistory
    import subprocess

    datapath = ctsys.resolve('regression/unittest/listhistory')
else:
    import commands
    from __main__ import default
    from tasks import *
    from taskinit import *

    datapath = os.environ.get('CASAPATH').split()[0] +\
        '/data/regression/unittest/listhistory/'

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
        # CASA5 tasks return False, casatasks throw exceptions
        myms = ''
        if is_CASA6:
            self.assertRaises(Exception,listhistory,myms)
        else:
            res = listhistory(myms)
            self.assertFalse(res)
            
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
        # the number of expected lines differs
        if is_CASA6:
            refnum=16
            if self.itismms:
                # this is a guess, not tested
                refnum = 40
        else:
            # for CASA5, get only the relevant lines in the logfile
            newfile= "newlisth.log"
            cmd="sed -n \"/Begin Task/,/End Task/p\" %s > %s " %(logfile,newfile)
            print(cmd)
            os.system(cmd)
            logfile = newfile

            refnum = 13
            if self.itismms:
                refnum = 36

        if is_CASA6:
            cmd=['wc', '-l', logfile]
            print(cmd)
            output = bytes2str(subprocess.check_output(cmd))
        else:
            cmd="wc -l %s" % logfile
            print(cmd)
            output=commands.getoutput(cmd)

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

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

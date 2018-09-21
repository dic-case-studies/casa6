import shutil
import unittest
import os
import sys
from casatasks import flagdata, flagmanager
from casatools import ctsys, agentflagger


# Path for data
datapath = ctsys.resolve('regression/unittest/flagdata')

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/flagdata/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR

print('flagmanager tests will use data from '+datapath)

# jagonzal (CAS-4287): Add a cluster-less mode to by-pass parallel processing for MMSs as requested 
if 'BYPASS_PARALLEL_PROCESSING' in os.environ:
    ParallelTaskHelper.bypassParallelProcessing(1)

# Local copy of the agentflagger tool
aflocal = agentflagger( )

# Base class which defines setUp functions for importing different data sets
class test_base(unittest.TestCase):
    
    def setUp_flagdatatest(self):
        '''VLA data set, scan=2500~2600 spw=0 1 chan, RR,LL'''
        self.vis = "flagdatatest.ms"
        if testmms:
            self.vis = "flagdatatest.mms"

        if os.path.exists(self.vis):
            print("The MS is already around, just unflag")
        else:
            print("Moving data...")
            os.system('cp -r ' + os.path.join(datapath,self.vis) + ' '+ self.vis)

        os.system('rm -rf ' + self.vis + '.flagversions')

    def setUp_bpass_case(self):
        self.vis = "cal.fewscans.bpass"

        if os.path.exists(self.vis):
            print("The CalTable is already around, just unflag")
        else:
            print("Moving data...")
            os.system('cp -r ' + os.path.join(datapath,self.vis) + ' ' + self.vis)

        os.system('rm -rf ' + self.vis + '.flagversions')        

    def unflag_ms(self):
        aflocal.open(self.vis)
        aflocal.selectdata()
        agentUnflag={'apply':True,'mode':'unflag'}
        aflocal.parseagentparameters(agentUnflag)
        aflocal.init()
        aflocal.run(writeflags=True)
        aflocal.done()


class test_flagmanager1(test_base):
    
    def setUp(self):
        os.system("rm -rf flagdatatest.ms*") # test1 needs a clean start
        self.setUp_flagdatatest()
        
    def test1m(self):
        '''flagmanager test1m: mode=list, flagbackup=True/False'''
        self.unflag_ms()        

        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 2)
        aflocal.done()


        flagdata(vis=self.vis, mode='unflag', flagbackup=False)
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 2)
        aflocal.done()

        flagdata(vis=self.vis, mode='unflag', flagbackup=True)
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 3)
        aflocal.done()
        
        newname = 'Ha! The best version ever!'

        flagmanager(vis=self.vis, mode='rename', oldname='flagdata_1', versionname=newname, 
                    comment='This is a *much* better name')
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 3)
        aflocal.done()
        
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversion file does not exist: flags.'+newname)
        
        # Specific for MMS
        if testmms:
            areg = self.vis+'/SUBMSS/*flagversions*'
            import glob
            print('Check for .flagversions in the wrong place.')
            self.assertEqual(glob.glob(areg), [], 'There should not be any .flagversions in the'
                                                ' SUBMSS directory')
            

    def test2m(self):
        """flagmanager test2m: Create, then restore autoflag"""

        self.unflag_ms()        

        flagdata(vis=self.vis, mode='summary')
        flagmanager(vis=self.vis)
        
        flagdata(vis=self.vis, mode='manual', antenna="2", flagbackup=True)
        
        flagmanager(vis=self.vis)
        ant2 = flagdata(vis=self.vis, mode='summary')['flagged']

        print("After flagging antenna 2 there were", ant2, "flags")

        # Change flags, then restore
        flagdata(vis=self.vis, mode='manual', antenna="3", flagbackup=True)
        flagmanager(vis = self.vis)
        ant3 = flagdata(vis=self.vis, mode='summary')['flagged']

        print("After flagging antenna 2 and 3 there were", ant3, "flags")

        flagmanager(vis=self.vis, mode='restore', versionname='flagdata_2')
        restore2 = flagdata(vis=self.vis, mode='summary')['flagged']

        print("After restoring pre-antenna 2 flagging, there are", restore2, "flags; should be", ant2)

        self.assertEqual(restore2, ant2)

    def test_CAS2701(self):
        """flagmanager: Do not allow flagversions with empty versionname''"""
                  
        try:
            flagmanager(vis = self.vis,mode = "save",versionname = "")     
        except IOError as e:
            print('Expected exception: %s'%e)
        
    def test_rename(self):
        '''flagmanager: do not overwrite an existing versionname'''
        
        # Create a flagbackup
        flagdata(vis=self.vis, mode='manual', antenna="2", flagbackup=True)
        fname = 'flagdata_1'
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+fname),
                        'Flagversions file does not exist: flags.'+fname)
        
        # Rename it
        newname = 'Rename_Me'
        flagmanager(vis=self.vis, mode='rename', oldname=fname, versionname=newname, 
                    comment='CAS-3080')
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversions file does not exist: flags.'+newname)
        
        self.assertFalse(os.path.exists(self.vis+'.flagversions/flags.'+fname),
                        'Flagversions file shuold not exist: flags.'+fname)
       
    def test_caltable_flagbackup(self):
        '''Flagmanager:: cal table mode=list, flagbackup=True/False'''
        # Need a fresh start
        os.system('rm -rf cal.fewscans.bpass*')
        self.setUp_bpass_case()
        self.unflag_ms()        
       
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 2)
        aflocal.done()

        flagdata(vis=self.vis, mode='unflag', flagbackup=False)
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 2)
        aflocal.done()

        flagdata(vis=self.vis, mode='unflag', flagbackup=True)
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 3)
        aflocal.done()
        
        newname = 'Ha! The best version ever!'

        flagmanager(vis=self.vis, mode='rename', oldname='flagdata_1', versionname=newname, 
                    comment='This is a *much* better name')
        flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 3)
        aflocal.done()
        
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversion file does not exist: flags.'+newname)        

    def test_save(self):
        '''flagmanager: CAS-3080, do not overwrite an existing versionname'''
        
        # Create a flagbackup
        flagdata(vis=self.vis, mode='manual', antenna="2", flagbackup=True)
        fname = 'flagdata_1'
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+fname),
                        'Flagversions file does not exist: flags.'+fname)
        
        # Rename
        newname = 'Do_Not_Overwrite_Me'
        print('Rename versionname to Do_Not_Overwrite_Me')
        flagmanager(vis=self.vis, mode='save', versionname=newname)
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversions file does not exist: flags.'+newname)
        
        print('Move existing versionname to temporary name')
        flagmanager(vis=self.vis, mode='save', versionname=newname)
        flagmanager(vis=self.vis, mode='list')
        lf = os.listdir(self.vis+'.flagversions')
        self.assertTrue([s for s in lf if '.old.' in s])
        self.assertEqual(len(lf), 4)
        

# Cleanup class 
class cleanup(test_base):
    
    def tearDown(self):
        os.system('rm -rf flagdatatest.*ms*')

    def test_runTest(self):
        '''flagdata: Cleanup'''
        pass


def suite():
    return [test_flagmanager1, cleanup]
    
if __name__ == '__main__':
    unittest.main()

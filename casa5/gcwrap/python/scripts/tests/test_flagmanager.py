from __future__ import absolute_import
from __future__ import print_function
import shutil
import unittest
import os
import sys

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import flagdata, flagmanager
    from casatools import ctsys, agentflagger
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper

    # CASA6 doesn't need default
    def default(atask):
        pass

    # Path for data
    datapath = ctsys.resolve('unittest/flagmanager/')
else:
    from tasks import flagdata, flagmanager
    from taskinit import aftool as agentflagger
    from __main__ import default
    from parallel.parallel_task_helper import ParallelTaskHelper
    from casa_stack_manip import stack_frame_find

    # to rethrow exception - not relevant in casatasks
    import inspect
    g = stack_frame_find( )
    g['__rethrow_casa_exceptions'] = True

    # Path for data
    datapath = os.environ.get('CASAPATH').split()[0] + "/casatestdata/unittest/flagmanager/"

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
aflocal = agentflagger()

# Base class which defines setUp functions for importing different data sets
class test_base(unittest.TestCase):
    
    def setUp_flagdatatest(self):
        '''VLA data set, scan=2500~2600 spw=0 1 chan, RR,LL'''
        self.vis = "flagdatatest.ms"

        if os.path.exists(self.vis):
            print("The MS is already around, just unflag")
        else:
            print("Moving data...")
            os.system('cp -RH ' + os.path.join(datapath,self.vis) + ' '+ self.vis)

        os.system('rm -rf ' + self.vis + '.flagversions')
        
        default(flagdata)

    def setUp_bpass_case(self):
        self.vis = "cal.fewscans.bpass"

        if os.path.exists(self.vis):
            print("The CalTable is already around, just unflag")
        else:
            print("Moving data...")
            os.system('cp -RH ' + os.path.join(datapath,self.vis) + ' ' + self.vis)

        os.system('rm -rf ' + self.vis + '.flagversions')        
        default(flagdata)

    def unflag_ms(self):
        aflocal.open(self.vis)
        aflocal.selectdata()
        agentUnflag={'apply':True,'mode':'unflag'}
        aflocal.parseagentparameters(agentUnflag)
        aflocal.init()
        aflocal.run(writeflags=True)
        aflocal.done()


class test_flagmanager_ms(test_base):
    
    def setUp(self):
#        os.system("rm -rf flagdatatest.ms*") # test1 needs a clean start
        # the .flagversions of the MS is always deleted before running a new test
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
        except ValueError as e:
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

    def test_save(self):
        '''flagmanager: CAS-3080, do not overwrite an existing versionname'''
        
        # Create a flagbackup
        flagdata(vis=self.vis, mode='manual', antenna="2", flagbackup=True)
        fname = 'flagdata_1'
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+fname),
                        'Flagversions file does not exist: flags.'+fname)
        
        # Rename
        newname = 'Do_Not_Overwrite_Me'
        print ('Rename versionname to Do_Not_Overwrite_Me')
        flagmanager(vis=self.vis, mode='save', versionname=newname)
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversions file does not exist: flags.'+newname)
        
        print('Move existing versionname to temporary name')
        flagmanager(vis=self.vis, mode='save', versionname=newname)
        flagmanager(vis=self.vis, mode='list')
        lf = os.listdir(self.vis+'.flagversions')
        self.assertTrue([s for s in lf if '.old.' in s])
        self.assertEqual(len(lf), 4)


    def test_dictionary1(self):
        '''flagmanager: Test the returned dictionary in list mode'''
        self.unflag_ms()        

        # MS should start only with the standard info such as:
        # MS : /diska/scastro/work/cas11462/nosedir/cal.fewscans.bpass
        # main : working copy in main table
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 2)
        aflocal.done()

        flagsdict1 = flagmanager(vis=self.vis, mode='list')

        self.assertIsInstance(flagsdict1, dict, 'mode=list should return a dictionary of the flag versions')

        flagdata(vis=self.vis, mode='unflag', flagbackup=False)
        flagsdict2 = flagmanager(vis=self.vis, mode='list')
        self.assertDictEqual(flagsdict1, flagsdict2, 'Dictionaries should contain the same flag versions')

        flagdata(vis=self.vis, mode='unflag', flagbackup=True)
        flagsdict3 = flagmanager(vis=self.vis, mode='list')
        self.assertEqual(flagsdict3[0]['name'],'flagdata_1')
        
        # Rename the flagdata_1 version
        newname = 'A new name'

        flagmanager(vis=self.vis, mode='rename', oldname='flagdata_1', versionname=newname, 
                    comment='A new versions name')
        flagsdict4 = flagmanager(vis=self.vis, mode='list')        
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversion file does not exist: flags.'+newname)
        self.assertEqual(flagsdict4[0]['name'], newname, 'Dictionary should show new name')
        
        # Specific for MMS
        if testmms:
            areg = self.vis+'/SUBMSS/*flagversions*'
            import glob
            print('Check for .flagversions in the wrong place.')
            self.assertEqual(glob.glob(areg), [], 'There should not be any .flagversions in the'
                                                ' SUBMSS directory')


class test_flagmanager_caltables(test_base):
                                            
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
        
        newname = 'A much better caltable version !'

        flagmanager(vis=self.vis, mode='rename', oldname='flagdata_1', versionname=newname, 
                    comment='This is a *much* better name')
        flagsdict = flagmanager(vis=self.vis, mode='list')
        aflocal.open(self.vis)
        self.assertEqual(len(aflocal.getflagversionlist()), 3)
        aflocal.done()
        
        self.assertTrue(os.path.exists(self.vis+'.flagversions/flags.'+newname),
                        'Flagversion file does not exist: flags.'+newname)        

        self.assertEqual(flagsdict[0]['name'], newname, 'Dictionary should show new name of flag versions')
        print(flagsdict)

# Cleanup class 
class cleanup(test_base):
    
    def tearDown(self):
        os.system('rm -rf flagdatatest.*ms*')

    def test_runTest(self):
        '''flagdata: Cleanup'''
        pass


def suite():
    return [test_flagmanager_ms, 
            test_flagmanager_caltables,
            cleanup]

if is_CASA6:    
    if __name__ == '__main__':
        unittest.main()

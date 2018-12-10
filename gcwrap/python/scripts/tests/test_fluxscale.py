import os
import shutil
import testhelper as th
import numpy as np
from __main__ import default
from tasks import fluxscale
from taskinit import *
import unittest
import exceptions


''' Python unit tests for the fluxscale task

These tests will only verify if the fluxscale
tables created for an MS and an MMS agree. These are
not full unit tests for the fluxscale task.
'''

datapath = os.environ.get('CASAPATH').split()[0] +\
                            '/data/regression/unittest/fluxscale/'


# Pick up alternative data directory to run tests on MMSs
testmms = False
if os.environ.has_key('TEST_DATADIR'):   
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/fluxscale/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print 'WARN: directory '+DATADIR+' does not exist'

print 'fluxscale tests will use data from '+datapath         

    
class fluxscale1_test(unittest.TestCase):

    def setUp(self):
        # Input names
        self.prefix = 'ngc5921'
        self.msfile = self.prefix + '.ms'
#        if testmms:
#            self.msfile = self.prefix + '.mms'
                        
        self.gtable = self.prefix + '.ref1a.gcal'
        self.reffile = self.prefix + '.def.fcal'
        self.reffile2 = self.prefix + '.def.inc.fcal'
        self.tearDown()
        
        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):        
            shutil.copytree(fpath, self.msfile, symlinks=True)
            fpath = os.path.join(datapath,self.gtable)
            shutil.copytree(fpath, self.gtable, symlinks=True)
            fpath = os.path.join(datapath,self.reffile)
            shutil.copytree(fpath, self.reffile, symlinks=True)
            fpath = os.path.join(datapath,self.reffile2)
            shutil.copytree(fpath, self.reffile2, symlinks=True)
        else:
            self.fail('Data does not exist -> '+fpath)

        default('fluxscale')

    def tearDown(self):
        #pass
        shutil.rmtree(self.msfile, ignore_errors=True)        
        os.system('rm -rf ngc5921*.fcal')
        os.system('rm -rf ngc5921*.gcal')

        
    def test_default(self):
        '''Fluxscale test 1.1: Create a flux table using field=0 as reference'''
        
        # Input
        gtable = self.gtable
        
        # Output
        outtable = self.msfile + '.fcal'
        
        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*', 
                  transfer='1445*')
        self.assertTrue(os.path.exists(outtable))
        
        # File to compare with
        reference = self.reffile
        
        # Compare the calibration table with a reference
        self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

        # compare some determined values returned in the dict

        #refdict={'1': {'spidxerr': np.array([ 0.,  0.,  0.]), 'spidx': np.array([ 0.,  0.,  0.]), \
        #         'fluxdErr': np.array([0.00055571]), \
        #         'fieldName': '1445+09900002_0', 'numSol': np.array([54]), \
                 #'fluxd': np.array([0.16825763])}, \
                 # flux density seems changed a bit. Updated - 2013.01.29 TT
		 # for linux on current trunk 22670
		 # for OSX 10.6 got the previous value 
        #         'fluxd': np.array([0.16825765])}, \
        #         'freq': np.array([1.41266507e+09]), \
        #         'spwName': np.array(['none'], dtype='|S5'), \
        #         'spwID': np.array([0])} 
      
        # new returned dictionary (2013.09.12 TT)
        refdict={'1': {'fitRefFreq': 0.0, 
                       'spidxerr': np.array([ 0.,  0.,  0.]), 
                       'spidx': np.array([ 0.,  0.,  0.]), 
                       '0': {'fluxdErr': np.array([ 0.00055574,  0.        ,  0.        ,  0.        ]), 
                             'numSol': np.array([ 54.,   0.,   0.,   0.]), 
                             'fluxd': np.array([ 0.16825768,  0.        ,  0.        ,  0.        ])}, 
                       'fitFluxd': 0.0, 
                       'fieldName': '1445+09900002_0', 
                       'fitFluxdErr': 0.0}, 
                       'freq': np.array([  1.41266507e+09]), 
                 'spwName': np.array(['none'],dtype='|S5'), 
                 'spwID': np.array([0], dtype=np.int32)}

        
        diff_fluxd=abs(refdict['1']['0']['fluxd'][0]-thisdict['1']['0']['fluxd'][0])/refdict['1']['0']['fluxd'][0]
        #self.assertTrue(diff_fluxd<1.5e-8)
	# increase the tolerance level
        self.assertTrue(diff_fluxd<1e-5)
        
            
    def test_incremental(self): 
        '''Fluxscale test 1.2: Create an incremental flux table using field=0 as reference'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.msfile + '.inc.fcal'

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*',
                  transfer='1445*', incremental=True)
        self.assertTrue(os.path.exists(outtable))

        # File to compare with
        reference = self.reffile2

        # Compare the calibration table with a reference
        self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

    def test_gainthreshold(self):
        '''Fluxscale test 1.3: gainthreshold parameter test'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.msfile + '.thres.fcal'

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*',
                  transfer='1445*', gainthreshold=0.05,incremental=True)
        self.assertTrue(os.path.exists(outtable))

        # File to compare with
        #reference = self.reffile2

        # Compare the calibration table with a reference
        #self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

    def test_antennasel(self):
        '''Fluxscale test 1.4: antenna de-selection test'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.msfile + '.antsel.fcal'

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*',
                  transfer='1445*', antenna='!24',incremental=True)
        self.assertTrue(os.path.exists(outtable))

        # File to compare with
        #reference = self.reffile2

        # Compare the calibration table with a reference
        #self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

    def test_antennaselwithtime(self):
        '''Fluxscale test 1.5: empy selection case: antenna with time selection test'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.msfile + '.antsel.fcal'

        # This time selection deselect all the data for the reference source and would raise an exception.
        try:
          thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*',
                  transfer='1445*', antenna='!24', timerange='>1995/04/13/09:38:00', incremental=True)
        except exceptions.RuntimeError, instance:
          print "Expected exception raised:",instance

    def test_antennaselwithscan(self):
        '''Fluxscale test 1.6: antenna selection with scan selection test'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.msfile + '.antsel.fcal'

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*',
                  transfer='1445*', antenna='!24', scan='1~5', incremental=True)
        self.assertTrue(os.path.exists(outtable))

    def test_refintransfer(self):
        '''Fluxscale test 1.7: test CAS-10227 fix: reference field included in transfer fields'''

        #input
        gtable = self.gtable

        # Output
        outtable = self.msfile + '.test1.7.fcal'

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1331*',
                  transfer='1445*,1331*', incremental=True)
        self.assertTrue(os.path.exists(outtable))
        self.assertFalse(thisdict.has_key('0'))

class fluxscale2_test(unittest.TestCase):

    def setUp(self):
        # Input names
        prefix = 'ngc4826'
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'
            
        self.gtable = prefix + '.spw.gcal'
        self.reffile = prefix + '.spw.fcal'

        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):
            shutil.copytree(fpath, self.msfile, symlinks=True)
            fpath = os.path.join(datapath,self.gtable)
            shutil.copytree(fpath, self.gtable, symlinks=True)
            fpath = os.path.join(datapath,self.reffile)
            shutil.copytree(fpath, self.reffile, symlinks=True)
        else:
            self.fail('Data does not exist -> '+fpath)

        default('fluxscale')
           
            
    def tearDown(self):
        shutil.rmtree(self.msfile, ignore_errors=True)
        os.system('rm -rf ngc4826*.gcal')
        os.system('rm -rf ngc4826*.fcal')
        
    def test_spws(self):
        '''Fluxscale 2: Create a fluxscale table for an MS with many spws'''
        # Input
        gtable = self.gtable
        
        # Output
        outtable = self.msfile + '.fcal'
        
        # torelance for value tests
        tol = 1.e-5

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='3C273-F0',
                             transfer=['1310+323-F0'],refspwmap=[0,0])
        self.assertTrue(os.path.exists(outtable))
        
        # File to compare with
        reference = self.reffile
        
        # Compare the calibration table with a reference
        self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))
        
        # compare some determined values returned in the dict

        #refdict={'1': {'spidxerr': np.array([ 0.,  0.,  0.]), 'spidx': np.array([ 0.,  0.,  0.]), \
        #         'fluxdErr': np.array([-1.        ,  0.04080052, -1.        , -1.        , -1.        , -1.        ]), \
        #         'fieldName': '1310+323-F0', 'numSol': np.array([-1,  8, -1, -1, -1, -1], dtype=np.int32), \
        #         'fluxd': np.array([-1.        ,  1.44578847, -1.        , -1.        , -1.        , -1.        ])}, \
        #         'freq': np.array([  1.15138579e+11,   1.15217017e+11,  -1.00000000e+00, -1.00000000e+00,  -1.00000000e+00, \
        #                         -1.00000000e+00]), 'spwName': np.array(['', '', '', '', '', ''], dtype='|S1'), \
        #         'spwID': np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)} 

        # updated for new returned dictionary format (2013.09.12 TT)
        refdict= {'1': {'fitRefFreq': 0.0, 
                        'spidxerr': np.array([ 0.,  0.,  0.]), 
                        'fitFluxd': 0.0, 
                        'spidx': np.array([ 0.,  0.,  0.]), 
                        '1': {'fluxdErr': np.array([ 0.04080052,  0.        ,  0.        ,  0.        ]), 
                              'numSol': np.array([ 8.,  0.,  0.,  0.]), 
                              'fluxd': np.array([ 1.44578847,  0.        ,  0.        ,  0.        ])}, 
                        '0': {'fluxdErr': np.array([-1., -1., -1., -1.]), 
                              'numSol': np.array([-1., -1., -1., -1.]), 
                              'fluxd': np.array([-1., -1., -1., -1.])}, 
                        '3': {'fluxdErr': np.array([-1., -1., -1., -1.]), 
                              'numSol': np.array([-1., -1., -1., -1.]), 
                              'fluxd': np.array([-1., -1., -1., -1.])}, 
                        '2': {'fluxdErr': np.array([-1., -1., -1., -1.]), 
                              'numSol': np.array([-1., -1., -1., -1.]), 
                              'fluxd': np.array([-1., -1., -1., -1.])}, 
                        '5': {'fluxdErr': np.array([-1., -1., -1., -1.]), 
                              'numSol': np.array([-1., -1., -1., -1.]), 
                              'fluxd': np.array([-1., -1., -1., -1.])}, 
                        '4': {'fluxdErr': np.array([-1., -1., -1., -1.]), 
                              'numSol': np.array([-1., -1., -1., -1.]), 
                              'fluxd': np.array([-1., -1., -1., -1.])}, 
                        'fieldName': '1310+323-F0', 
                        'fitFluxdErr': 0.0}, 
                  'freq': np.array([  1.15138579e+11,   1.15217017e+11,  -1.00000000e+00, -1.00000000e+00,  
                                      -1.00000000e+00,  -1.00000000e+00]), 
                  'spwName': np.array(['', '', '', '', '', ''], dtype='|S1'), 
                  'spwID': np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)}

        diff_fluxd=abs(refdict['1']['1']['fluxd'][0]-thisdict['1']['1']['fluxd'][0])/refdict['1']['1']['fluxd'][0]

        self.assertTrue(diff_fluxd<tol)


class fluxscale3_test(unittest.TestCase):

    def setUp(self):
        # Input names
        # ngc2403.tutorial.ms reduced in size
        prefix = 'n2403.short'
        self.prefix = prefix
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'

        self.gtable = prefix + '.flagFld1.gcal'
        self.reffile = prefix + '.flagFld1.ref.fcal'
        self.gtable2 = prefix + '.flagPartFld3.gcal'
        self.reffile2 = prefix + '.flagPartFld3.ref.fcal'

        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):
            shutil.copytree(fpath, self.msfile, symlinks=True)
            fpath = os.path.join(datapath,self.gtable)
            shutil.copytree(fpath, self.gtable, symlinks=True)
            fpath = os.path.join(datapath,self.reffile)
            shutil.copytree(fpath, self.reffile, symlinks=True)
            fpath = os.path.join(datapath,self.gtable2)
            shutil.copytree(fpath, self.gtable2, symlinks=True)
            fpath = os.path.join(datapath,self.reffile2)
            shutil.copytree(fpath, self.reffile2, symlinks=True)
        else:
            self.fail('Data does not exist -> '+fpath)

        default('fluxscale')

    def tearDown(self):
        shutil.rmtree(self.msfile, ignore_errors=True)
        os.system('rm -rf n2403*.gcal')
        os.system('rm -rf n2403*.fcal')

    def test_flaggedref1(self):
        '''Fluxscale test3: Ref field 1 in caltable is all flagged'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.prefix + '.flagFld1.fcal'

        # torelance for value test
        tol = 1.e-5

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1,3,4',
                             transfer='2')
        self.assertTrue(os.path.exists(outtable))

        # File to compare with
        reference = self.reffile

        # Compare the calibration table with a reference
        self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

        # compare some determined values returned in the dict (tested on RHEL6)
        refdict={'freq': np.array([  1.41825202e+09]), 
                 '2': {'fitRefFreq': 0.0, 
                       'spidxerr': np.array([ 0.,  0.,  0.]),
                       'spidx': np.array([ 0.,  0.,  0.]),
                       '0': {'fluxdErr': np.array([ 0.00189683,  0.,  0.,  0.]), 
                             'numSol': np.array([ 54.,   0.,   0.,   0.]),
                             'fluxd': np.array([ 3.20832802,  0.,  0.,  0.])}, 
                 'fitFluxd': 0.0, 'fieldName': '0841+708', 'fitFluxdErr': 0.0}, 
                 'spwName': np.array(['127*24.4 kHz channels @ 1.42 GHz (BARY)'],dtype='|S40'),
                 'spwID': np.array([0], dtype=np.int32)} 

        diff_fluxd=abs(refdict['2']['0']['fluxd'][0]-thisdict['2']['0']['fluxd'][0])/refdict['2']['0']['fluxd'][0]

        self.assertTrue(diff_fluxd<tol)

 
    def test_flaggedref2(self):
        '''Fluxscale test3: Ref field 3 in caltable is partially flagged'''

        # Input
        gtable = self.gtable2

        # Output
        outtable = self.prefix + '.flagPartFld3.fcal'

        # torelance for value test
        tol = 1.e-5

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1,3,4',
                             transfer='2')
        self.assertTrue(os.path.exists(outtable))

        # File to compare with
        reference = self.reffile2

        # Compare the calibration table with a reference
        self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

        # compare some determined values returned in the dict (tested on RHEL6)
        refdict={'freq': np.array([  1.41825202e+09]), 
                 '2': {'fitRefFreq': 0.0, 'spidxerr': np.array([ 0.,  0.,  0.]), 
                       'spidx': np.array([ 0.,  0.,  0.]), 
                       '0': {'fluxdErr': np.array([ 0.0022236,  0.,  0.,  0.]), 
                             'numSol': np.array([ 54.,   0.,   0.,   0.]), 
                             'fluxd': np.array([ 3.19455455,  0.,  0.,  0.])}, 
                             'fitFluxd': 0.0, 'fieldName': '0841+708', 'fitFluxdErr': 0.0}, 
                       'spwName': np.array(['127*24.4 kHz channels @ 1.42 GHz (BARY)'], dtype='|S40'), 
                       'spwID': np.array([0], dtype=np.int32)} 


        diff_fluxd=abs(refdict['2']['0']['fluxd'][0]-thisdict['2']['0']['fluxd'][0])/refdict['2']['0']['fluxd'][0]

        self.assertTrue(diff_fluxd<tol)

class fluxscale_fit_test(unittest.TestCase):

    def setUp(self):
        # Input names
        # shorten CasaGuide data set (TW Hya)
        prefix = 'twhya.short'
        self.prefix = prefix
        self.msfile = prefix + '.ms'
#        if testmms:
#            self.msfile = prefix + '.mms'

        self.gtable = prefix + '.fittest.gcal'
        self.reffile = prefix + '.fittest.ref.fcal'

        fpath = os.path.join(datapath,self.msfile)
        if os.path.lexists(fpath):
            shutil.copytree(fpath, self.msfile, symlinks=True)
            fpath = os.path.join(datapath,self.gtable)
            shutil.copytree(fpath, self.gtable, symlinks=True)
            #fpath = os.path.join(datapath,self.reffile)
            #shutil.copytree(fpath, self.reffile, symlinks=True)
        else:
            self.fail('Data does not exist -> '+fpath)

        default('fluxscale')

    def tearDown(self):
        shutil.rmtree(self.msfile, ignore_errors=True)
        os.system('rm -rf twhya*.gcal')
        os.system('rm -rf twhya*.fcal')

    def test_spectralfit1(self):
        '''Fluxscale spectral fit test1'''
        # Input
        gtable = self.gtable

        # Output
        outtable = self.prefix + '.fittest.fcal'

        # torelance for value test
        tol = 1.e-5

        thisdict = fluxscale(vis=self.msfile, caltable=gtable, fluxtable=outtable, reference='1',
                             refspwmap=[0,1,3,3], listfile='twhya.fittest1.out.txt')

        self.assertTrue(os.path.exists(outtable))

        # File to compare with
        reference = self.reffile

        # Compare the calibration table with a reference
        #self.assertTrue(th.compTables(outtable, reference, ['WEIGHT']))

        # compare some determined values returned in the dict (tested on RHEL6)
        refdict={'0': {'0': {'fluxd': np.array([ 10.40151801,   0.        ,   0.        ,   0.        ]),
                             'fluxdErr': np.array([ 0.01823489,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 16.,   0.,   0.,   0.])},
                       '1': {'fluxd': np.array([ 10.49918965,   0.        ,   0.        ,   0.        ]),
                             'fluxdErr': np.array([ 0.02617547,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 16.,   0.,   0.,   0.])},
                       '2': {'fluxd': np.array([ 10.40252265,   0.        ,   0.        ,   0.        ]),
                             'fluxdErr': np.array([ 0.04820387,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 14.,   0.,   0.,   0.])},
                       '3': {'fluxd': np.array([ 10.06407263,   0.        ,   0.        ,   0.        ]),
                             'fluxdErr': np.array([ 0.02031524,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 16.,   0.,   0.,   0.])},
                 'covarMat': np.array([[  1.32836567e-06,  -2.99318676e-05],
                                   [ -2.99318676e-05,   2.06267774e-02]]),
                 'fieldName': '3c279',
                 'fitFluxd': 10.284275366907444,
                 'fitFluxdErr': 0.049681178272988555,
                 'fitRefFreq': 350998123030.1949,
                 'spidx': np.array([ 1.0121737 ,  0.86003842,  0.        ]),
                 'spidxerr': np.array([ 0.00209799,  0.26143238,  0.        ])},
                 '3': {'0': {'fluxd': np.array([ 1.02073715,  0.        ,  0.        ,  0.        ]),
                             'fluxdErr': np.array([ 0.01662615,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 16.,   0.,   0.,   0.])},
                       '1': {'fluxd': np.array([ 0.82954399,  0.        ,  0.        ,  0.        ]),
                             'fluxdErr': np.array([ 0.01892758,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 16.,   0.,   0.,   0.])},
                       '2': {'fluxd': np.array([ 1.00187618,  0.        ,  0.        ,  0.        ]),
                             'fluxdErr': np.array([ 0.02459043,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 14.,   0.,   0.,   0.])},
                       '3': {'fluxd': np.array([ 1.30754902,  0.        ,  0.        ,  0.        ]),
                             'fluxdErr': np.array([ 0.01957324,  0.        ,  0.        ,  0.        ]),
                             'numSol': np.array([ 16.,   0.,   0.,   0.])},
               'covarMat': np.array([[  8.53253971e-05,   9.38851637e-04],
                                 [  9.38851637e-04,   1.31604054e+00]]),
               'fieldName': 'J1147-382=QSO',
               'fitFluxd': 1.0679456164935712,
               'fitFluxdErr': 0.07095309114437398,
               'fitRefFreq': 350998123030.1949,
               'spidx': np.array([ 0.02854914, -7.30077968,  0.        ]),
               'spidxerr': np.array([ 0.02885403,  3.58345512,  0.        ])},
               '4': {'0': {'fluxd': np.array([ 0.95582976,  0.        ,  0.        ,  0.        ]),
                           'fluxdErr': np.array([ 0.00889469,  0.        ,  0.        ,  0.        ]),
                           'numSol': np.array([ 16.,   0.,   0.,   0.])},
                     '1': {'fluxd': np.array([ 0.74538553,  0.        ,  0.        ,  0.        ]),
                           'fluxdErr': np.array([ 0.01522706,  0.        ,  0.        ,  0.        ]),
                           'numSol': np.array([ 16.,   0.,   0.,   0.])},
                     '2': {'fluxd': np.array([ 0.96176398,  0.        ,  0.        ,  0.        ]),
                           'fluxdErr': np.array([ 0.01663815,  0.        ,  0.        ,  0.        ]),
                           'numSol': np.array([ 14.,   0.,   0.,   0.])},
                     '3': {'fluxd': np.array([ 1.23699681,  0.        ,  0.        ,  0.        ]),
                           'fluxdErr': np.array([ 0.01410021,  0.        ,  0.        ,  0.        ]),
                           'numSol': np.array([ 16.,   0.,   0.,   0.])},
                  'covarMat': np.array([[  4.01358110e-05,  -2.64936568e-04],
                                     [ -2.64936568e-04,   6.57883034e-01]]),
                  'fieldName': 'J1037-295=QSO',
                  'fitFluxd': 1.0177140843936607,
                  'fitFluxdErr': 0.06550518651825477,
                  'fitRefFreq': 350998123030.1949,
                  'spidx': np.array([ 0.00762578, -6.81332936,  0.        ]),
                  'spidxerr': np.array([ 0.02795337,  3.57884208,  0.        ])},
                  'freq': np.array([  3.56732311e+11,   3.57968689e+11,   3.45799939e+11, 3.43721561e+11]),
                  'spwID': np.array([0, 1, 2, 3], dtype=np.int32),
                  'spwName': np.array(['', '', '', ''], 
                 dtype='|S1')}


        diff_fluxd=abs(refdict['0']['0']['fluxd'][0]-thisdict['0']['0']['fluxd'][0])/refdict['0']['0']['fluxd'][0]
        self.assertTrue(diff_fluxd<tol)
       
        diff_covar=abs(refdict['3']['covarMat'][0][0]-thisdict['3']['covarMat'][0][0])/refdict['3']['covarMat'][0][0]
        self.assertTrue(diff_covar<tol)


def suite():
    return [fluxscale1_test, fluxscale2_test, fluxscale3_test, fluxscale_fit_test]


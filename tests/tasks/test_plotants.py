from __future__ import absolute_import
import os
import string
import sys
import shutil
import unittest
from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import plotants
    from casatools import ctsys
    ctsys_resolve = ctsys.resolve
else:
    from __main__ import default
    from tasks import *
    #from taskinit import *

    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'data')
        return os.path.join(dataPath,apath)

'''
Unit tests for task plotants. It tests the following parameters:
    vis:           wrong and correct values
    figfile:       if output is created
'''
class plotants_test(unittest.TestCase):
    # Input and output names
    msfile = 'ic2233_1.ms'
    res = None
    fig = 'plotantstest.png'

    def setUp(self):
        self.res = None
        if not is_CASA6: default(plotants)

        # It is not necessary to copy it for all tests
        if (not os.path.exists(self.msfile)):
            datapath = ctsys_resolve('regression/ic2233')
            shutil.copytree(os.path.join(datapath,self.msfile), self.msfile)

    def tearDown(self):
        if (os.path.exists(self.msfile)):
            os.system('rm -rf ' + self.msfile)

        os.system('rm -rf ' + self.fig)

    @unittest.skipIf(is_CASA6,"failure, data not found")
    def test1(self):
       '''Test 1: Default parameters'''
       self.res = plotants()
       self.assertFalse(self.res)

    @unittest.skipIf(is_CASA6,"failure, data not found")
    def test2(self):
        '''Test 2: Bad input file'''
        msfile = 'badfile'
        self.res = plotants(vis=msfile)
        self.assertFalse(self.res)

    def test3(self):
        '''Test 3: Good input file and output exists'''
        self.res = plotants(vis=self.msfile, figfile=self.fig)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test4(self):
        '''Test 4: Label antenna IDs'''
        self.res = plotants(vis=self.msfile, figfile=self.fig, antindex=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test5(self):
        '''Test 5: Logarithmic antenna positions'''
        self.res = plotants(vis=self.msfile, figfile=self.fig, logpos=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    @unittest.skipIf(is_CASA6,"failure, unknown reasons")
    def test6(self):
        '''Test 6: Exclude antenna positions'''
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            exclude='1,5,19,14,10,13')
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test7(self):
        '''Test 7: checkbaselines'''
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            checkbaselines=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test8(self):
        '''Test 8: exclude checkbaselines'''
        # antenna (name) 11 is already excluded by checkbaselines
        # (warning)
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            exclude='11', checkbaselines=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test9(self):
        '''Test 9: Title'''
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            title='IC2233')
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test10(self):
        '''Test 10: All arguments'''
        self.res = plotants(self.msfile, self.fig, True, True, '1,3,5,7,9',
            True, "IC2233")
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

def suite():
    return [plotants_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

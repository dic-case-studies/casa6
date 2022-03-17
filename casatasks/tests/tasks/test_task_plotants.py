##########################################################################
# test_task_plotants.py
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.visualization.plotants.html
#
##########################################################################
import os
import string
import sys
import shutil
import unittest

from casatasks import plotants
from casatools import ctsys

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

        # It is not necessary to copy it for all tests
        if (not os.path.exists(self.msfile)):
            datapath = ctsys.resolve('unittest/plotants/')
            shutil.copytree(os.path.join(datapath,self.msfile), self.msfile)

    def tearDown(self):
        if (os.path.exists(self.msfile)):
            os.system('rm -rf ' + self.msfile)

        os.system('rm -rf ' + self.fig)

    def test1(self):
       '''Test 1: Default parameters'''
       self.assertRaises(Exception,plotants)

    def test2(self):
        '''Test 2: Bad input file'''
        msfile = 'badfile'
        self.assertRaises(Exception,plotants,vis=msfile)

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

if __name__ == '__main__':
    unittest.main()

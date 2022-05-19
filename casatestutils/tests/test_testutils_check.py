##########################################################################
import csv
import fnmatch
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import sys
import unittest

from casatools import ctsys
import casatools
tb = casatools.table()
from casaplotms import plotms
from casatestutils import check
from casaplotms import plotmstool as pm

# Paths for data
datapath = ctsys.resolve("unittest/plotms/")
tinyim = ctsys.resolve("unittest/imsmooth/tiny.im")


class plotms_test_base(unittest.TestCase):

    outputDir="/tmp/" + str(os.getpid()) + "/"
    plotfile_jpg = "/tmp/myplot.jpg"
    display = os.environ.get("DISPLAY")
    testms  = "pm_ngc5921.ms"
    ms = os.path.join(outputDir, testms)


# ------------------------------------------------------------------------------

    def checkPlotfile(self, plotfileName, minSize, maxSize=None):
        self.assertTrue(os.path.isfile(plotfileName), "Plot was not created")
        plotSize = os.path.getsize(plotfileName)
        print(plotfileName, 'file size is', plotSize)
        self.assertGreater(plotSize, minSize)
        if maxSize:
            self.assertLess(plotSize, maxSize)

    def checkDisplay(self):
        self.assertGreater(len(self.display), 0, 'DISPLAY not set, cannot run test')

    def removePlotfile(self, plotfile=None):
        try:
            if not plotfile:
                plotfile = self.plotfile_jpg
            os.remove(plotfile)
        except OSError:  # "No such file or directory"
            pass

    def cleanUp(self):
        if os.path.exists(self.outputDir):
            shutil.rmtree(self.outputDir)

    def setUpData(self):
        if not os.path.exists(self.ms):
            shutil.copytree(os.path.join(datapath, self.testms),
                self.ms, symlinks=True)

    def tearDownData(self):
        self.cleanUp()
        pm.setPlotMSFilename("")
# ------------------------------------------------------------------------------
 
class test_check(plotms_test_base):
    ''' Test check class from casatestutils'''

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        self.ch = check.Check()

    def tearDown(self):
        self.tearDownData()

    def test_check_plotfile_minsize(self):
        '''Test Check Basic plot with minsize'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic01.jpg")
        self.removePlotfile()
        plotms(vis=self.ms, plotfile=self.plotfile_jpg, showgui=False, highres=True)
        self.assertTrue(self.ch.check_plotfile(self.plotfile_jpg, 190000))

    def test_check_plotfile_maxsize(self):
        '''Test Check Basic plot with maxsize'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic01.jpg")
        self.removePlotfile()
        plotms(vis=self.ms, plotfile=self.plotfile_jpg, showgui=False, highres=True)
        #plotSize = os.path.getsize(self.plotfile_jpg)
        #print(self.plotfile_jpg, 'file size is', plotSize)
        self.assertTrue(self.ch.check_plotfile(self.plotfile_jpg, 0,max_size=220000))

    def test_check_plotfile_bad_maxsize(self):
        '''Test Check Basic plot with bad maxsize '''
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic01.jpg")
        self.removePlotfile()
        plotms(vis=self.ms, plotfile=self.plotfile_jpg, showgui=False, highres=True)
        #plotSize = os.path.getsize(self.plotfile_jpg)
        #print(self.plotfile_jpg, 'file size is', plotSize)
        self.assertFalse(self.ch.check_plotfile(self.plotfile_jpg, 0,max_size=1))

    def test_check_pixels(self):
        '''Test Check Pixels from single image slice'''
        # Open Image, get first slice, return first slice matches itself
        tb.open(tinyim)
        image = tb.getcol('map')
        tb.close()

        ref = image[0]

        self.assertTrue(self.ch.check_pixels(imagename=tinyim, loc='0', refval=ref, rtol=1e-05, atol=1e-08))

    def test_check_pixels_multi(self):
        '''Test Check Pixels from single image multiple slices'''
        # Open Image, get first slice, return first slice matches itself
        tb.open(tinyim)
        image = tb.getcol('map')
        tb.close()

        ref = image[tuple([0,1])]

        self.assertTrue(self.ch.check_pixels(imagename=tinyim, loc='0,1', refval=ref, rtol=1e-05, atol=1e-08))

    def test_check_pixels_slice(self):
        '''Test Check Pixels from single image slice range'''
        # Open Image, get first slice, return first slice matches itself
        tb.open(tinyim)
        image = tb.getcol('map')
        tb.close()

        ref = image[tuple([slice(0,2)])]

        self.assertTrue(self.ch.check_pixels(imagename=tinyim, loc='0:2', refval=ref, rtol=1e-05, atol=1e-08))

    def test_check_pixels_exception_noref(self):
        '''Test Check Pixels throws exception for no ref'''
        with self.assertRaises(Exception):
            self.ch.check_pixels(imagename=tinyim, loc='0', rtol=1e-05, atol=1e-08)


if __name__ == '__main__':
    unittest.main()

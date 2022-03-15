########################################################################
# test_task_specflux.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.specflux.html
#
#
##########################################################################
import os
import shutil
import unittest

from casatools import ctsys
from casatasks import specflux

datapath = ctsys.resolve('unittest/specflux/')
refpath = ctsys.resolve('unittest/specflux/specflux_reference/')

im1 = os.path.join(datapath,"specflux1.im")
im2 = os.path.join(datapath,"specflux2.im")

class specflux_test(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def tearDown(self):
        for i in range(1, 11):
            f = 'log' + str(i)
            if os.path.exists(f):
                os.remove(f)
    
    def test_default(self):
        """Test default settings"""
        logfile = "log1"
        specflux(im1, logfile=logfile)
        self._compare(logfile)
        logfile = "log2"
        specflux(im2, logfile=logfile)
        self._compare(logfile)
        
    def test_box_chans(self):
        """test setting box and channel range"""
        logfile = "log3"
        specflux(im1, box="4,4,15,15", chans="5~20", logfile=logfile)
        self._compare(logfile)
        logfile = "log4"
        specflux(im2, box="10,10,19,19", chans="30~35", logfile=logfile)
        self._compare(logfile)
        
    def test_box_chans_mask(self):
        """test setting box and chans with a mask"""
        logfile = "log5"
        specflux(
            im1, box="4,4,15,15", chans="5~20",
            mask="'" + im1 + "'" + ">0", logfile=logfile
        )
        self._compare(logfile)
        logfile = "log6"
        specflux(
            im2, box="10,10,19,19", chans="30~35",
            mask="'" + im2 + "'" + "<0", logfile=logfile
        )
        self._compare(logfile)
        
    def test_unit(self):
        """test setting spectral unit"""
        logfile = "log7"
        specflux(im1, unit="GHz", logfile=logfile)
        self._compare(logfile)
        logfile = "log8"
        specflux(im2, unit="kHz", logfile=logfile)
        self._compare(logfile)
        
    def test_beam(self):
        """test setting beam info"""
        logfile = "log9"
        specflux(im1, major="8arcmin", minor="6arcmin", logfile=logfile)
        self._compare(logfile)
        logfile = "log10"
        specflux(im2, major="8arcmin", minor="6arcmin", logfile=logfile)
        # brightness unit is K, so setting beam should have no effect
        self._compare(logfile, "log2")
        
    def _compare(self, gfile, efile=""):
        if not efile:
            efile = gfile
        with open (gfile) as f:
            got = f.readlines()
        with open (os.path.join(refpath,efile)) as f:
            expec = f.readlines()
        self.assertTrue(len(got) == len(expec))
        # skip first element because paths will in general be different
        for (myg, mye) in zip(got[2:], expec[2:]):
            self.assertTrue(
                myg == mye,
                gfile + " mismatch: got: " + myg + " exp: " + mye
            )
 
if __name__ == '__main__':
    unittest.main()

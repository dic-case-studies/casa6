##########################################################################
# test_tool_coordsys.py
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.coordsys.html
#
#
##########################################################################
import shutil
import unittest
import os

from casatools import image as iatool
from casatools import coordsys
from casatools import ctsys
from casatools import table

myim = "center_0.fits"

datapath='unittest/coordsys/'

class coordsys_test(unittest.TestCase):
    
    def setUp(self):
        shutil.copy(ctsys.resolve(datapath + myim), myim)
        
    def tearDown(self):
        os.remove(myim)
        tb = table( )
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done( )
        
    def test_CAS_2724(self):
        myia = iatool()
        myia.open(myim)
        mycsys = myia.coordsys()
        self.assertRaises(Exception, mycsys.toworldmany, [1])

    def test_constructor(self):
        """Test constructors"""
        x = coordsys( )
        
    def test_findaxisbyname(self):
        myia = iatool()
        myia.fromshape("", [4, 4, 4, 4])
        csys = myia.coordsys( )
        myia.done()

        self.assertTrue(csys.findaxisbyname("ri") == 0)
        self.assertRaises(Exception, csys.findaxisbyname, "ra", False)
        self.assertTrue(csys.findaxisbyname("ra", True) == 0)
        
    def test_transpose(self):
        """Test the transpose() method"""
        myia = iatool()
        myia.fromshape("", [4, 4, 4, 4])
        csys = myia.coordsys()
        myia.done()
        
        orig = csys.names()
        self.assertTrue(csys.transpose([3,2,1,0]))
        new = csys.names()
        self.assertTrue(orig[0] == new[3])
        self.assertTrue(orig[1] == new[2])
        self.assertTrue(orig[2] == new[1])
        self.assertTrue(orig[3] == new[0])
        
    def test_findaxis(self):
        """Test the findaxis() method"""
        myia = iatool()
        myia.fromshape("", [20, 20, 4, 20])
        csys = myia.coordsys()
        self.assertRaises(Exception, csys.findaxis, True, -1)
        self.assertRaises(Exception, csys.findaxis, True, 4)
        got = csys.findaxis(True, 0)
        self.assertTrue(got == {'axisincoordinate': 0, 'coordinate': 0})
        got = csys.findaxis(True, 1)
        self.assertTrue(got == {'axisincoordinate': 1, 'coordinate': 0})
        got = csys.findaxis(True, 2)
        self.assertTrue(got == {'axisincoordinate': 0, 'coordinate': 1})
        got = csys.findaxis(True, 3)
        self.assertTrue(got == {'axisincoordinate': 0, 'coordinate': 2})
        myia.done()
        
    def test_findcoordinate(self):
        """Test the findcoordinate() method"""
        myia = iatool()
        myia.fromshape("", [20, 20, 4, 20])
        csys = myia.coordsys()
        myia.done()
        self.assertFalse(csys.findcoordinate('linear')['return'])
        self.assertFalse(
            csys.findcoordinate('direction', 2)['return']
        )
        got = csys.findcoordinate("direction")
        self.assertTrue(
            got['return'] and (got['pixel'] == [0,1]).all()
            and (got['world'] == [0,1]).all()
        )
        got = csys.findcoordinate("spectral")
        self.assertTrue(
            got['return'] and (got['pixel'] == [3]).all()
            and (got['world'] == [3]).all()
        )
        got = csys.findcoordinate("stokes")
        self.assertTrue(
            got['return'] and (got['pixel'] == [2]).all()
            and (got['world'] == [2]).all()
        )

if __name__ == '__main__':
    unittest.main()

##########################################################################
# test_me_shift.py
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
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.measures.html
#
##########################################################################

import unittest
import copy
import math
from casatools import measures
from casatools import quanta

class me_shift_test(unittest.TestCase):

    def setUp(self):
        self.me = measures( )
        self.qa = quanta( )
        pass
    
    def tearDown(self):
        self.me.done( )
        pass

    def test_shift(self):
        """Test me.shift"""
        v = self.me.direction("J2000", "4h20m30s", "+30.20.30")
        got = self.me.shift(v, "20arcmin", "0deg")
        expec = copy.deepcopy(v)
        expec['m1'] = self.qa.add(expec['m1'], "20arcmin")
        self.assertTrue(got == expec)
        got = self.me.shift(v, "20arcmin", "90deg")
        expec = 1.1433867531223854
        self.assertTrue(abs(got['m0']['value']/expec - 1) < 1e-7)
        expec = 0.5295520783025025
        self.assertTrue(abs(got['m1']['value']/expec - 1) < 1e-7)
        got = self.me.shift(v, "20arcmin", "180deg")
        self.assertTrue(got['m0']['value'] == v['m0']['value'])
        expec = self.qa.sub(v['m1'], '20arcmin')
        self.assertTrue(abs(got['m1']['value']/expec['value'] - 1) < 1e-7)

if __name__ == '__main__':
    unittest.main()

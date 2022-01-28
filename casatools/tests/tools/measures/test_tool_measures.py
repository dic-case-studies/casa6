##########################################################################
# test_measures.py
#
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
# Unit tests for the measures tool.
#
# Features tested:
#  1. me.cometdist
#  2. me.cometangdiam
#  3. me.shift
##########################################################################

import copy
import unittest

from casatools import  ctsys, measures, quanta
ctsys_resolve = ctsys.resolve

def check_eq(val, expval, tol=None):
    """Checks that val matches expval within tol."""
    if type(val) == dict:
        for k in val:
            check_eq(val[k], expval[k], tol)
    else:
        try:
            if tol and hasattr(val, '__rsub__'):
                are_eq = abs(val - expval) < tol
            else:
                are_eq = val == expval
            if hasattr(are_eq, 'all'):
                are_eq = are_eq.all()
            if not are_eq:
                raise ValueError('!=')
        except ValueError:
            errmsg = "%r != %r" % (val, expval)
            if (len(errmsg) > 66): # 66 = 78 - len('ValueError: ')
                errmsg = "\n%r\n!=\n%r" % (val, expval)
            raise ValueError(errmsg)
        except Exception as e:
            print("Error comparing %s to %s" % (val,expval))
            raise e

class Ganymede(unittest.TestCase):
    """
    Base class for Ganymede tests.
    """
    def setUp(self):
        self.me = measures( )
        self.qa = quanta()
        cometdir = ctsys_resolve("ephemerides/JPL-Horizons/")
        self.me.framecomet(cometdir + "Ganymede_55437-56293dUTC.tab")
        self.me.doframe(self.me.epoch("utc", "2011/01/03/17:00:00"))
        self.me.doframe(self.me.observatory("ALMA"))

    def tearDown(self):
        self.me.done( )
        pass

class me_test_cometdist(Ganymede):
    def test_cometdist(self):
        """Is the distance to Ganymede just right?"""
        check_eq(self.me.cometdist(), {'unit': 'AU', 'value': 5.1241}, 0.001)
            
class me_test_cometangdiam(Ganymede):
    def test_cometangdiam(self):
        """Is Ganymede's angular diameter correct?"""
        check_eq(self.me.cometangdiam(), {'unit': 'rad', 'value': 6.868e-06},
                 1.0e-9)

class me_test_shift(Ganymede):

    def test_shift(self):
        """Test me.shift"""
        v = self.me.direction("J2000", "4h20m30s", "+30.20.30")
        got = self.me.shift(v, "20arcmin", "0deg")
        expec = copy.deepcopy(v)
        expec['m1'] = self.qa.add(expec['m1'], "20arcmin")
        self.assertTrue(got == expec)
        got = self.me.shift(v, "20arcmin", "90deg")
        expec = 1.1433867531223854
        self.assertTrue(abs(got['m0']['value'] / expec - 1) < 1e-7)
        expec = 0.5295520783025025
        self.assertTrue(abs(got['m1']['value'] / expec - 1) < 1e-7)
        got = self.me.shift(v, "20arcmin", "180deg")
        self.assertTrue(got['m0']['value'] == v['m0']['value'])
        expec = self.qa.sub(v['m1'], '20arcmin')
        self.assertTrue(abs(got['m1']['value'] / expec['value'] - 1) < 1e-7)

if __name__ == '__main__':
    unittest.main()

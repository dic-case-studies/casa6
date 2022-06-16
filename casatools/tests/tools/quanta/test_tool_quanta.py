##########################################################################
# test_tool_quanta.py
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
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.quanta.html
#
#
##########################################################################
import shutil
import unittest
import os

from casatools import quanta
import numpy as np
import re
import math

qa = quanta()

class quanta_test(unittest.TestCase):
   

    def setUp(self):
       pass


    def tearDown(self):
        qa.done()
        

# Tests for quanta.convert
class quanta_convert_test(quanta_test):


    def test_exceptions(self):
        """Test various exception cases"""
        def __test_exception(method_parms, expected_msg):
            with self.assertRaises(RuntimeError) as cm:
                res = qa.convert(**method_parms)
            got_exception = cm.exception
            pos = str(got_exception).find(expected_msg)
            self.assertNotEqual(
                pos, -1, msg=f'{got_exception}'
            )


       # bogus unit
        parms = {}
        parms['v'] = qa.quantity('5kg')
        parms['outunit'] = 'bogus'
        __test_exception(
            parms,
            'Error in QuantumHolder::fromString with input string "bogus": '
            'Illegal input units or format\n in converting quantity'
        )
 

    def test_value_as_string(self):
        """Test specifying value as a string"""
        for ks in (True, False):
            res = qa.quantity('5kg', keepshape=ks)
            self.assertEqual(res, {'value': 5, 'unit': 'kg'})
            # test unitname is ignored in this case
            res = qa.quantity('5kg', 'pc', keepshape=ks)
            self.assertEqual(res, {'value': 5, 'unit': 'kg'})


    def test_valid(self):
        """Test specifying the quantity as a scalar"""
        # v is string
        res = qa.convert('5kg', 'g')
        self.assertEqual(
            res, {'value': 5000, 'unit': 'g'}, 'incorrect conversion'
        )
        # v is quantity with scalar value
        res = qa.convert(qa.quantity('5kg'), 'g')
        self.assertEqual(
            res, {'value': 5000, 'unit': 'g'}, 'incorrect conversion'
        )
        # v is quantity with vector of values
        x = [5, 9.25, 8]
        res = qa.convert(qa.quantity(x, 'kg'), 'g')
        self.assertEqual(res['unit'], 'g', 'incorrect conversion')
        self.assertTrue((res['value']/1000 == x).all(), 'incorrect conversion') 
        # v is a multi-dimensional array
        x = np.random.random([5, 2, 8, 9])
        res = qa.convert(qa.quantity(x, 'kg', keepshape=True), 'g')
        self.assertEqual(res['unit'], 'g', 'incorrect conversion')
        self.assertTrue(
            np.isclose(res['value']/1000, x).all(), 'incorrect conversion'
        )
        res = qa.convert(5, 'm')
        self.assertTrue(
            res == {'unit': 'm.m-1', 'value': 5}, 'incorrect conversion'
        )
        """
        res = qa.convert[4, 8], 'm')
        self.assertTrue(
            res == {'unit': 'm.m-1', 'value': 5}, 'incorrect conversion'
        )
        """


# Tests for quanta.quantity
class quanta_quantity_test(quanta_test):


    def test_exceptions(self):
        """Test various exception cases"""
        def __test_exception(method_parms, expected_msg):
            with self.assertRaises(RuntimeError) as cm:
                res = qa.quantity(**method_parms)
            got_exception = cm.exception
            pos = str(got_exception).find(expected_msg)
            self.assertNotEqual(
                pos, -1, msg=f'{got_exception}'
            )
            
        parms = {}
        # string is bad quantity
        parms['v'] = 'blah'
        __test_exception(
            parms, 'Error in QuantumHolder::fromString with input string '
            '"blah": Illegal input units or format\n in converting quantity'
        )
        # bad unit
        parms = {}
        parms['v'] = 5
        parms['unitname'] = 'zz'
        __test_exception(parms, "Unit::check Illegal unit string 'zz'")


    def test_value_as_string(self):
        """Test specifying value as a string"""
        for ks in (True, False):
            res = qa.quantity('5kg', keepshape=ks)
            self.assertEqual(res, {'value': 5, 'unit': 'kg'})
            # test unitname is ignored in this case
            res = qa.quantity('5kg', 'pc', keepshape=ks)
            self.assertEqual(res, {'value': 5, 'unit': 'kg'})


    def test_scalar_value(self):
        """Test specifying the quantity as a scalar"""
        for ks in (True, False):
            res = qa.quantity(5, 'kg', keepshape=ks)
            self.assertEqual(res, {'value': 5, 'unit': 'kg'})


    def test_vector_value(self):
        """Test specifying the quantity as a vector"""
        v = [5, 8, 7.5]
        for ks in (True, False):
            res = qa.quantity(v, 'kg', keepshape=ks)
            self.assertEqual(res['unit'], 'kg', 'incorrect unit')
            self.assertTrue((res['value'] == v).all(), 'incorrect vector value(s)')


    def test_multidimensional_arra_value(self):
        """Test specifying the quantity as a multidimensional array"""
        v = np.random.random((4,5))
        # keepshape = False returns a 1D array
        res = qa.quantity(v, 'kg', keepshape=False)
        self.assertEqual(res['unit'], 'kg', 'incorrect unit')
        self.assertTrue(
            np.allclose(res['value'], v.ravel('F')), 'incorrect vector value(s)'
        )
        # keepshape = True preserves input array shape
        res = qa.quantity(v, 'kg', keepshape=True)
        self.assertEqual(res['unit'], 'kg', 'incorrect unit')
        self.assertTrue(
            np.allclose(res['value'], v), 'incorrect vector value(s)'
        )


    def test_value_as_record(self):
        """
        Test value as record (eg already a quantity).
        This is essentially a reflection operator.
        """
        v = {'unit': 'arcsec', 'value': 840}
        # scalar
        self.assertEqual(qa.quantity(v), v, 'scalar reflection failed')
        # vector
        v = {'unit': 'arcsec', 'value': [840., 840., 840., 840.]}
        r = qa.quantity(v)
        self.assertEqual(
            r['unit'], v['unit'], 'vector reflection failed on unit test'
        )
        self.assertTrue(
            (r['value'] == v['value']).all(), 'vector reflection failed'
        )
        # array
        v = {'unit': 'arcsec', 'value': [[840., 840., 840., 840.]]}
        r = qa.quantity(v)
        self.assertEqual(
            r['unit'], v['unit'], 'array reflection failed on unit test'
        )
        self.assertTrue(
            (r['value'] == v['value']).all(), 'array reflection failed'
        )


if __name__ == '__main__':
    unittest.main()

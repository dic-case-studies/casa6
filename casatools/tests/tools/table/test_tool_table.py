##########################################################################
# test_tool_table.py
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
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.table.html
#
# Methods tested in this script
# getcellslice
##########################################################################
import shutil
import unittest
import os
import numpy as np
import math

from casatools import image, table

class TableBase(unittest.TestCase):

    def setUp(self):
        self.ia = image()
        self.tb = table()

    def tearDown(self):
        self.ia.done()
        self.tb.done()

# Tests for table.getcellslice
class table_getcellslice_test(TableBase):
    
    myim = 'tmp.im'
    arr = np.array([
        [1,  2,  3,  4],
        [5,  6,  7,  8],
        [9, 12, 11, 12]
    ])

    def tearDown(self):
        super().tearDown()
        for f in [self.myim]:
            if os.path.exists(f) and os.path.isdir(f):
               shutil.rmtree(f)


    def test_exceptions(self):
        """Test various exception cases"""
        def __test_exception(method_parms, expected_msg):
            with self.assertRaises(RuntimeError) as cm:
                res = self.tb.getcellslice(**method_parms)
            got_exception = cm.exception
            pos = str(got_exception).find(expected_msg)
            self.assertNotEqual(
                pos, -1, msg=f'Unexpected exception was thrown: {got_exception}'
            )
            
        self.ia.fromarray(self.myim, self.arr)
        self.ia.done()
        self.tb.open(self.myim)
        # bad column name
        parms = {}
        parms['columnname'] = 'bogus'
        parms['rownr'] = 0
        parms['blc'] = -1
        parms['trc'] = -1
        __test_exception(
            parms, f'Table column {parms["columnname"]} is unknown'
        ) 
        # bad row number
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 1
        parms['blc'] = -1
        parms['trc'] = -1
        __test_exception(
            parms, f'TableColumn: row number {parms["rownr"]} exceeds #rows 1 '
            f'in table {os.path.dirname(os.path.abspath(self.myim))}'
        )
        # bad blc
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = [0, 0, 0]
        parms['trc'] = -1
        __test_exception(parms, 'blc must have length of 2')
        # blc too low
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = [0, -1]
        parms['trc'] = -1
        __test_exception(
            parms, 'All elements of blc must be greater than or equal to 0'
        )
        # blc too high
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = [0, 4]
        parms['trc'] = -1
        __test_exception(parms, 'Element 1 of blc must be less than 4')
        # bad trc
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = -1
        parms['trc'] = [5, 5, 5]
        __test_exception(parms, 'trc must have length of 2')
        # trc negative
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = -1
        parms['trc'] = [2, -1]
        __test_exception(
            parms, 'All elements of trc must be greater than or equal to 0'
        )
        # trc too large
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = -1
        parms['trc'] = [4, 4]
        __test_exception(parms, 'Element 0 of trc must be less than 3')
        # trc less than blc
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = [1, 2]
        parms['trc'] = [2, 0]
        __test_exception(
            parms,
            'All elements of trc must be greater than their corresponding blc '
            'elements'
        )
        # trc equal to blc
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = [1, 2]
        parms['trc'] = [2, 2]
        __test_exception(
            parms,
            'All elements of trc must be greater than their corresponding blc '
            'elements'
        )
        # too many values for inc
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = -1
        parms['trc'] = -1
        parms['incr'] = [1,1,1]
        __test_exception(
            parms, 'incr must have length of 2'
        )
        # incr negative
        parms = {}
        parms['columnname'] = 'map'
        parms['rownr'] = 0
        parms['blc'] = -1
        parms['trc'] = -1
        parms['incr'] = [1, 0]
        __test_exception(
            parms, 'All elements of incr must be greater than 0'
        )


    def test_getcellslice(self):
        """tests for valid inputs"""
        self.ia.fromarray(self.myim, self.arr)
        self.ia.done()
        self.tb.open(self.myim)
        # get the entire cell
        z = self.tb.getcellslice('map', 0, -1, -1)
        self.assertTrue((z == self.arr).all(), 'getting entire array failed')
        z = self.tb.getcellslice('map', 0, [0, 0], -1)
        self.assertTrue((z == self.arr).all(), 'getting entire array failed')
        z = self.tb.getcellslice('map', 0, -1, [2 ,3])
        self.assertTrue((z == self.arr).all(), 'getting entire array failed')
        z = self.tb.getcellslice('map', 0, -1, -1, [1, 1])
        self.assertTrue((z == self.arr).all(), 'getting entire array failed')
        z = self.tb.getcellslice('map', 0, [1, 1], -1)
        self.assertTrue((
            z == self.arr[1:, 1:]).all(), 'setting blc to non-zero failed'
        )
        z = self.tb.getcellslice('map', 0, -1, [1, 1])
        self.assertTrue((
            z == self.arr[:2, :2]).all(), 'setting trc to not trc of array '
            'failed'
        )
        z = self.tb.getcellslice('map', 0, -1, -1, [2, 2])
        self.assertTrue((
            z == self.arr[::2, ::2]).all(), 'setting incr to larger than 1 '
            'failed'
        )
        self.tb.done()


if __name__ == '__main__':
    unittest.main()

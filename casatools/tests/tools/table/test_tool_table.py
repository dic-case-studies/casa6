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

    def test_getcellslice(self):
        """all tests"""
        self.ia.fromarray(self.myim, self.arr)
        self.ia.done()
        self.tb.open(self.myim)
        # bad cell name
        z = self.tb.getcellslice('map_bogus', 0, [0,0,0], [2,2,2], [1,1,1])

        # get the entire cell
        z = self.tb.getcellslice('map', 0, [0,0,0], [2,2,2], [1,1,1])
        z = self.tb.getcellslice('map', 0, -1, -1)
        self.assertTrue((z == self.arr).all(), 'getting entire array failed')
        self.tb.done()


if __name__ == '__main__':
    unittest.main()

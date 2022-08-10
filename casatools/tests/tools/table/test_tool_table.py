##########################################################################
# test_tool_table.py
#
# Copyright (C) 2022
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
import itertools
import os, stat
import numpy as np
from uuid import uuid4
from casatools import table, ctsys, image
import math


##########################################################################
ms_name = 'n08c1_swap1.ms'
orig_ms_path = ctsys.resolve( f'unittest/table/{ms_name}' )
print( f'table tool tests will use {orig_ms_path}' )


class TableBase(unittest.TestCase):
    "setup common to all tests"

    @classmethod
    def setUpClass(cls):
        cls.scratch_path = str(uuid4( ))
        cls.ms_path = os.path.join(cls.scratch_path,ms_name)

    @staticmethod
    def remove_readonly(func, path, _):
        "Clear the readonly bit and reattempt the removal"
        os.chmod(path, stat.S_IWRITE)
        func(path)

    def setUp(self):
        self.tb = table( )
        self.ia = image()
        if os.path.exists(self.scratch_path):
            shutil.rmtree( self.scratch_path, onerror=self.remove_readonly )
        if not os.path.exists(self.scratch_path):
            os.makedirs(self.scratch_path)
        shutil.copytree( orig_ms_path, self.ms_path )
        self.tb.open(self.ms_path,nomodify=False)
        self.rows = self.tb.row( )

    def tearDown(self):
        self.rows.done( )
        self.tb.close( )
        self.tb.done( )
        self.ia.done()
        tb = table()
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done()
        if os.path.exists(self.scratch_path):
            shutil.rmtree( self.scratch_path, onerror=self.remove_readonly )


class TableRowTest(TableBase):
    def test_get(self):
        """Test get function"""
        ### fetch single row
        row = self.rows.get(21)
        ### check values
        self.assertTrue(np.isclose(np.abs(np.sum(row['DATA'])),0.8696193426581154))

    def test_shape( self ):
        """Test for valid data shape"""
        ### fetch single row
        row = self.rows.get(21)
        ### check shape
        self.assertTrue(row['DATA'].shape == (4, 32))
        
    def test_put(self):
        """Test put function"""

        ### fetch single row
        row = self.rows.get(22)
        ### modify DATA field
        row['DATA'].fill( 1.0/row['DATA'].size )
        ### store modified DATA
        self.rows.put( 22, { 'DATA': row['DATA'] } )
        ### check to see if the new values are read back
        tb = table( )
        tb.open(self.tb.name( ))
        rows = tb.row( )
        storedrow = self.rows.get(22)
        rows.done( )
        tb.close( )
        tb.done( )
        self.assertTrue(np.isclose(np.abs(np.sum(storedrow['DATA'])),1.0))

    def test_all_rows(self):
        ### read all rows from table
        allrows = self.rows[:]
        ### check number of rows & sum of data
        self.assertTrue(len(allrows) == 720 and np.isclose( sum([np.abs(np.sum(x['DATA'])) for x in allrows]), 22031.1419634223 ))

    def test_some_rows(self):
        ### read slice from table
        rows = self.rows[:15]
        ### check number of rows & sum of data
        self.assertTrue(len(rows) == 15 and np.isclose( sum([np.abs(np.sum(x['DATA'])) for x in rows]), 498.3275412917137 ))

    def test_columnnames_include(self):
        self.assertTrue( list(itertools.chain(*[ x.keys() for x in self.tb.row(columnnames=['DATA'])[:5] ])) == ['DATA', 'DATA', 'DATA', 'DATA', 'DATA'] )

    def test_columnnames_exclude(self):
        self.assertTrue( set(self.tb.row(columnnames=['DATA'],exclude=True).get(0).keys( )) == set( [ 'ANTENNA1', 'ANTENNA2', 'ARRAY_ID', 'DATA_DESC_ID',
                                                                                                       'EXPOSURE', 'FEED1', 'FEED2', 'FIELD_ID', 'FLAG', 
                                                                                                       'FLAG_CATEGORY', 'FLAG_ROW', 'INTERVAL', 'OBSERVATION_ID',
                                                                                                       'PROCESSOR_ID', 'SCAN_NUMBER', 'SIGMA', 'SIGMA_SPECTRUM', 
                                                                                                       'STATE_ID', 'TIME', 'TIME_CENTROID', 'UVW', 'WEIGHT', 
                                                                                                       'WEIGHT_SPECTRUM' ] ) )

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
            'All elements of trc must be greater than or equal to their '
            'corresponding blc elements'
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

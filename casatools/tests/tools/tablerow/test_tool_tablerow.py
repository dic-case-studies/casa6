##########################################################################
# test_tool_tablerow.py
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
##########################################################################
import shutil
import unittest
import os, stat
import numpy as np
from uuid import uuid4
from casatools import table, ctsys

ms_name = 'ngc5921_ut.ms'
orig_ms_path = ctsys.resolve( f'unittest/mstool/ngc5921_ut.ms')
if orig_ms_path == 'unittest/mstool':
    raise RuntimeError(f'cannot find test data {orig_ms_path}')
print( f'tablerow tool tests will use {orig_ms_path}' )



class TableRowBase(unittest.TestCase):
    "setup common to all tests"

    def setUp(self):
        self.scratch_path = str(uuid4( ))
        self.ms_path = os.path.join(self.scratch_path,ms_name)
        self.tb = table( )
        if os.path.exists(self.scratch_path):
            shutil.rmtree( self.scratch_path, onerror=remove_readonly )
        if not os.path.exists(self.scratch_path):
            os.makedirs(self.scratch_path)
        shutil.copytree( orig_ms_path, self.ms_path )
        self.tb.open(self.ms_path,nomodify=False)
        self.rows = self.tb.row( )

    @staticmethod
    def remove_readonly(func, path, _):
        "Clear the readonly bit and reattempt the removal"
        os.chmod(path, stat.S_IWRITE)
        func(path)

    def tearDown(self):
        self.rows.done( )
        self.tb.close( )
        self.tb.done( )
        tb = table()
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done()
        if os.path.exists(self.scratch_path):
            shutil.rmtree( self.scratch_path, onerror=self.remove_readonly )

class row(TableRowBase):
    def test_get(self):
        """Test get function"""
        ### fetch single row
        row = self.rows.get(21)
        ### check values
        self.assertTrue(np.isclose(np.abs(np.sum(row['DATA'])),7468.873846292496))

    def test_shape( self ):
        """Test for valid data shape"""
        ### fetch single row
        row = self.rows.get(21)
        ### check shape
        self.assertTrue(row['DATA'].shape == (2, 63))
        
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
        tb.open(self.ms_path)
        rows = tb.row( )
        storedrow = self.rows.get(22)
        rows.done( )
        tb.close( )
        tb.done( )
        self.assertTrue(np.isclose(np.abs(np.sum(storedrow['DATA'])),1.0))

class slice(TableRowBase):
    def test_all_rows(self):
        ### read all rows from table
        allrows = self.rows[:]
        ### check number of rows & sum of data
        self.assertTrue(len(allrows) == 22653 and np.isclose( sum([np.abs(np.sum(x['DATA'])) for x in allrows]), 12600276.218082037 ))

    def test_some_rows(self):
        ### read slice from table
        rows = self.rows[:15]
        ### check number of rows & sum of data
        self.assertTrue(len(rows) == 15 and np.isclose( sum([np.abs(np.sum(x['DATA'])) for x in rows]), 111601.736328125 ))

if __name__ == '__main__':
    unittest.main()

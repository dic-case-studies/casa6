########################################################################
# test_task_splattotable.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.splattotable.html
#
#
##########################################################################
import os
import shutil
import unittest

from casatools import ctsys, table, spectralline
from casatasks import splattotable, casalog

good_list = "splattotable_list1.txt"
bad_list = "splattotable_list2.txt"

def run_sttmethod(list, tab):
    mysl = spectralline()
    try:
        return mysl.splattotable(filenames=list, table=tab)
    except:
        raise
    mysl.done()

def run_stttask(list, tab):
    return splattotable(filenames=list, table=tab)


class splattotable_test(unittest.TestCase):
    
    def setUp(self):
        self._tb = table()
        datapath=ctsys.resolve('unittest/splattotable/')

        shutil.copy(os.path.join(datapath,good_list), good_list)
        shutil.copy(os.path.join(datapath,bad_list), bad_list)

    
    def tearDown(self):
        os.remove(good_list)
        os.remove(bad_list)
        self.assertTrue(len(self._tb.showcache()) == 0)

    def test_exceptions(self):
        """splattotable: Test various exception cases"""

        # these functions are used where exceptions are expected (here)
        # so that clean-up can happen if the exceptions do NOT occur (the test fails)
        # otherwise the test may leave a file open, causing other tests to appear to fail
        def check_run_sttmethod(filenames, tab):
            rst = run_sttmethod(filenames,tab)
            rst.done()
        
        def testit(filenames, tab):
            for i in [0,1]:
                if (i==0):
                    self.assertRaises(Exception, check_run_sttmethod, filenames, tab)
                else:
                    self.assertRaises(Exception, run_stttask, filenames, tab)

                # in all cases, nothing should be opened
                self.assertTrue(len(self._tb.showcache()) == 0)

        # CASA6 throws exceptions here, CASA5 does not
        # blank output table name
        try:
            testit(good_list, "")
        except:
            casalog.post("Failure in test_exceptions testing blank output table name",'SEVERE')
            raise
        
        # bad list
        try:
            testit(bad_list, "myout");
        except:
            casalog.post("Failure in test_exceptions testing bad list",'SEVERE')
            raise
        
        # unwritable table
        try:
            testit(good_list, "foo/bar/myout");
        except:
            casalog.post("Failure in test_exceptions testing unwritable table",'SEVERE')
            raise

    def test_good_list(self):
        """splattotable: Test converting a good list"""
        def testit(filenames, tab):
            mytb = table()
            for i in [0,1]:
                tab = tab + str(i)
                if (i==0):
                    newsl = run_sttmethod(filenames, tab)
                    newsl.done()
                else:
                    try:
                        run_stttask(filenames, tab)
                    except Exception:
                        self.fail()

                self.assertTrue(mytb.open(tab))
            mytb.done()
                    
        testit(good_list, "good_table")
      
if __name__ == '__main__':
    unittest.main()

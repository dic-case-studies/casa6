##########################################################################
# imfit_test.py
#
# Copyright (C) 2008, 2009
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
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be adressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
# <author>
# Dave Mehringer
# </author>
#
# <summary>
# Test suite for the CASA task splattotable
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
#   <li> <linkto class="task_splattotable.py:description">splattotable</linkto> 
# </ul>
# </prerequisite>
#
# <etymology>
# Test for the splattotable task
# </etymology>
#
# <synopsis>
# Test the splattotable task and the sl.splattotable() method upon which it is built.
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_splattotable[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the splattotable task to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
from __future__ import absolute_import
import os
import shutil
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, table, spectralline
    from casatasks import splattotable, casalog
else:
    import casac
    from tasks import *
    from taskinit import *
    from taskinit import sltool as spectralline
    from taskinit import tbtool as table
    from __main__ import *

good_list = "list1.txt"
bad_list = "list2.txt"

def run_sttmethod(list, tab):
    mysl = spectralline()
    try:
        return mysl.splattotable(filenames=list, table=tab)
    except:
        raise
    mysl.done()

def run_stttask(list, tab):
    if not is_CASA6:
        default(splattotable)
    return splattotable(filenames=list, table=tab)


class splattotable_test(unittest.TestCase):
    
    def setUp(self):
        self._tb = table()
        if is_CASA6:
            datapath=ctsys.resolve('regression/unittest/splattotable')
        else:
            datapath=os.path.join(os.environ.get('CASAPATH').split()[0],'data/regression/unittest/splattotable')

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
                    # CASA6 task raises an exception, CASA5 returns None
                    if is_CASA6:
                        self.assertRaises(Exception, run_stttask, filenames, tab)
                    else:
                        self.assertEqual(run_stttask(filenames, tab), None)
                # in all cases, nothing should be opened
                self.assertTrue(len(self._tb.showcache()) == 0)

        # CASA6 throws exceptions here, CASA5 does not
        # blank output table name
        try:
            testit("list1.txt", "")
        except:
            casalog.post("Failure in test_exceptions testing blank output table name",'SEVERE')
            raise
        
        # bad list
        try:
            testit("list2.txt", "myout");
        except:
            casalog.post("Failure in test_exceptions testing bad list",'SEVERE')
            raise
        
        # unwritable table
        try:
            testit("list1.txt", "foo/bar/myout");
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
                    
        testit("list1.txt", "good_table")
      

def suite():
    return [splattotable_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

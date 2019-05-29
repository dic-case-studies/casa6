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
    from casatasks import splattotable
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
    restool = mysl.splattotable(filenames=list, table=tab)
    mysl.close()
    return restool

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
        
        def testit(filenames, tab):
            for i in [0,1]:
                if (i==0):
                    self.assertRaises(Exception, run_sttmethod, filenames, tab)
                else:
                    self.assertEqual(run_stttask(filenames, tab), None)

        # CASA6 throws exceptions here, CASA5 does not
        # blank output table name
        try:
            OK = False
            testit("list1.txt", "")
            if not is_CASA6:
                OK = True
        except:
            if is_CASA6:
                OK = True
        self.assertEqual(OK,True)

        # bad list
        try:
            OK = False
            testit("list2.txt", "myout");
            if not is_CASA6:
                OK = True
        except:
            if is_CASA6:
                OK = True
        self.assertEqual(OK,True)

        # unwritable table
        try:
            OK = False
            testit("list1.txt", "/myout");
            if not is_CASA6:
                OK = True
        except:
            if is_CASA6:
                OK = True
        self.assertEqual(OK,True)

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
                    self.assertTrue(run_stttask(filenames, tab))
                    
                self.assertTrue(mytb.open(tab))
            mytb.done()
                    
        testit("list1.txt", "good_table")
      

def suite():
    return [splattotable_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

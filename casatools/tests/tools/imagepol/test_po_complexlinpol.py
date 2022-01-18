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
# Test suite for the CASA tool method po.complexlinpol
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
# </ul>
# </prerequisite>
#
# <etymology>
# Test for the po.complexlinpol() tool method
# </etymology>
#
# <synopsis>
# Test the po.complexlinpol() tool method
# </synopsis>
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
#
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/gcwrap/python/scripts/regressions/admin/runUnitTest.py test_po_complexlinpol[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.tofits() tool method to ensure
# coding changes do not break the associated bits
# </motivation>
#

###########################################################################
import shutil
import unittest

from casatools import imagepol as potool, table, ctsys
import os, shutil

datapath = ctsys.resolve('unittest/imagepol/')
eq_beams = "pol_eq_beams.fits"
neq_beams = "pol_neq_beams.fits"

class po_complexlinpol_test(unittest.TestCase):
    
    def setUp(self):
        self.mypo = potool()
    
    def tearDown(self):
        self.mypo.done()
        tb = table()
        self.assertEqual(len(tb.showcache()), 0)

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        shutil.copy(datapath + eq_beams, eq_beams)
        shutil.copy(datapath + neq_beams, neq_beams)
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.complexlinpol('g'))
        mypo.done()
        os.remove(eq_beams)
        shutil.rmtree('g')
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.complexlinpol, 'hh')
        mypo.done()
        os.remove(neq_beams)
        shutil.rmtree('hh')
        
def suite():
    return [po_complexlinpol_test]

if __name__ == '__main__':
    unittest.main()

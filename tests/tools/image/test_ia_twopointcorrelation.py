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
# Test suite for the CASA tool method ia.twopointcorrelation()
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
# Test for the ia.twopointcorrelation() tool method
# </etymology>
#
# <synopsis>
# Test the ia.twopointcorrelation() tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_twopointcorrelation[test1,test2,...]
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

from CASAtools import image as iatool

class ia_twopointcorrelation_test(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def tearDown(self):
        pass
    
    def test_stretch(self):
        """ ia.twopointcorrelation(): Test stretch parameter"""
        yy = iatool()
        mymask = "maskim"
        yy.fromshape(mymask, [20, 20, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [20,20,1,20]
        yy.fromshape("", shape)
        yy.addnoise()
        self.assertRaises(
            Exception,
            yy.twopointcorrelation,
            mask=mymask + ">0", stretch=False
        )
        zz = yy.twopointcorrelation(
            mask=mymask + ">0", stretch=True
        )
        self.assertTrue(zz and type(zz) == type(True))
        yy.done()

    def test_history(self):
        """verify history writing"""
        myia = iatool()
        myia.fromshape("zz",[20, 20])
        outfile = "xyz.im"
        # does not return an ia tool, just a bool
        self.assertTrue(myia.twopointcorrelation(outfile))
        myia.open(outfile)
        msgs = myia.history()
        myia.done()       
        self.assertTrue("ia.twopointcorrelation" in msgs[-2])
        self.assertTrue("ia.twopointcorrelation" in msgs[-1])
    
def suite():
    return [ia_twopointcorrelation_test]

if __name__ == '__main__':
    unittest.main()

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
# Test suite for the CASA method ia.deconvolvefrombeam()
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
#
# <etymology>
# Test for the method ia.deconvolvefrombeam()
# </etymology>
#
# <synopsis>
# Test the method ia.deconvolvefrombeam().
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_deconvolvefrombeam[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the method ia.deconvolvefrombeam() to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import shutil
import unittest

from casatools import image as iatool
from casatools import quanta

class ia_deconvolvefrombeam_test(unittest.TestCase):
    
    def setUp(self):
        self.ia = iatool( )
        self.qa = quanta( )
    
    def tearDown(self):
        self.ia.done()
        self.qa.done( )
    
    def __near(self, got, exp, tol):
        qgot = self.qa.quantity(got)
        qexp = self.qa.quantity(exp)
        return self.qa.abs(self.qa.div(self.qa.sub(qgot, qexp), qexp))["value"] < tol

    def test_multibeams(self):
        """ ia.deconvolvefrombeam(): Basic tests"""
        print("*** start ")
        myia = self.ia
        source = ["4arcsec","3arcsec", "20deg"]
        beam = [
            ["2arcsec", "1arcsec", "20deg"],
            ["3arcsec", "2arcsec", "50deg"]
        ]
        emaj = [
            self.qa.quantity({'unit': 'arcsec', 'value': 3.4641016151377548}),
            self.qa.quantity({'unit': 'arcsec', 'value': 3.0203474964295665})
        ]
        emin = [
            self.qa.quantity({'unit': 'arcsec', 'value': 2.8284271247461894}),
            self.qa.quantity({'unit': 'arcsec', 'value': 1.6963198403637358})
        ]
        epa = [
            self.qa.quantity({'unit': 'deg', 'value': 20}),
            self.qa.quantity({'unit': 'deg', 'value': -1.9489431240069859})
        ]
        tol = 1e-10
        for i in [0, 1]:
            res = myia.deconvolvefrombeam(source, beam[i])
            fit = res["fit"]
            self.assertTrue(self.__near(fit["major"], emaj[i], tol))
            self.assertTrue(self.__near(fit["minor"], emin[i], tol))
            print("*** got %s" % fit["pa"])
            print("*** exp %s" % epa[i])
            self.assertTrue(self.__near(fit["pa"], epa[i], tol))
        
def suite():
    return [ia_deconvolvefrombeam_test]

if __name__ == '__main__':
    unittest.main()

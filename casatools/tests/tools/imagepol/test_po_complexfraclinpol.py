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
# Test suite for the CASA tool method po.complexfraclinpol
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
# Test for the po.complexfraclinpol() tool method
# </etymology>
#
# <synopsis>
# Test the po.complexfraclinpol() tool method
# </synopsis> 
#
#
# <motivation>
# To provide a test standard for the ia.tofits() tool method to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import shutil
import unittest

try:
    from casatools import imagepol as potool
    from casatools import ctsys
    ctsys_resolve = ctsys.resolve
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
        return os.path.join(dataPath,apath)

datapath = ctsys_resolve('unittest/imagepol/')
eq_beams = datapath + "pol_eq_beams.fits"
neq_beams = datapath + "pol_neq_beams.fits"

class po_complexfraclinpol_test(unittest.TestCase):
    
    def setUp(self):
        self.mypo = potool()
    
    def tearDown(self):
        self.mypo.done()
    
    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.complexfraclinpol("g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.complexfraclinpol, "hh")
        
        
def suite():
    return [po_complexfraclinpol_test]

if __name__ == '__main__':
    unittest.main()

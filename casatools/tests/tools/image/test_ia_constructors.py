###########################################################################
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
# Test suite for various ia tool constructors
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_constructors[test1,test2,...]
#
# </example>

###########################################################################
import shutil
import unittest
import numpy as np
import os
from casatools import image as iatool
from casatools import table

class ia_constructors_test(unittest.TestCase):
    
    def setUp(self):
        self.tb = table( )
        self.outfile = ''
    
    def tearDown(self):
        self.tb.done( )
        if self.outfile:
            if os.path.isfile(self.outfile):
                os.unlink(self.outfile)
            else:
                shutil.rmtree(self.outfile)
    
    def test_newimagefromarray(self):
        """ test repeated call of newimagefromarray doesn't segfault, CAS-5646"""
        my_image = np.zeros([128,128,16])
        myia = iatool()
        self.outfile = "mynewimage.image"
        zz = myia.newimagefromarray(
            outfile=self.outfile , pixels=my_image,
            overwrite=True
        )
        myia.open(self.outfile )
        self.assertRaises(
            Exception, myia.newimagefromarray,
            outfile=self.outfile ,
            pixels=my_image,
            overwrite=True
        )
        myia.done()
        zz.done()
        self.assertTrue(len(self.tb.showcache()) == 0) 

    def test_history(self):
        """verify history writing"""
        myia = iatool()
        self.outfile = "zz"
        myia.fromshape(self.outfile ,[20, 20])
        csys = myia.coordsys()
        ary = myia.getchunk()
        myia = myia.newimagefromarray(pixels=ary, csys=csys.torecord())
        msgs = myia.history()
        myia.done()       
        self.assertTrue("ia.newimagefromarray" in msgs[-2])
        self.assertTrue("ia.newimagefromarray" in msgs[-1])

    def test_history1(self):
        """verify history writing"""
        myia = iatool()
        myia = myia.newimagefromshape(shape=[20,20])
        msgs = myia.history()
        myia.done()       
        self.assertTrue("ia.newimagefromshape" in msgs[-2])
        self.assertTrue("ia.newimagefromshape" in msgs[-1])

def suite():
    return [ia_constructors_test]

if __name__ == '__main__':
    unittest.main()

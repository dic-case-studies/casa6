#############################################
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
# Test suite for the CASA task specsmooth
# </summary>
#
# <motivation>
# To provide a test standard for the specsmooth task to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import os
import shutil
import unittest
import numpy
from math import sqrt

from casatools import ctsys, image, table, quanta, regionmanager
from casatasks import specsmooth

_rg = regionmanager( )
_tb = table( )

class specsmooth_test(unittest.TestCase):
    
    def setUp(self):
        pass

    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)

    def test_general(self):
        """Test general behavior"""
        myia = image()
        length = 6
        imagename = "test_gen.im"
        myia.fromshape(imagename, [1, 1, length])
        bb = myia.getchunk()
        for i in range(length):
            bb[0, 0, i] = i*i + 1
        myia.putchunk(bb)
        for i in range(length):
            reg = _rg.box([0, 0, 0], [0, 0, i])
            outfile = "out" + str(i) + ".im"
            if (i < 2):
                self.assertRaises(
                    Exception, specsmooth, imagename=imagename, outfile=outfile,
                    region=reg, function="h", axis=2
                )
            else:
                for drop in (False, True):
                    outfile = "out" + str(i) + str(drop) + ".im"
                    if drop:
                        for dmethod in ("c", "m"):
                            outfile = "out" + str(i) + str(drop) + dmethod + ".im"
                            if i==2 or i==3:
                                if dmethod=="c":
                                    expec = [2.5]
                                else:
                                    if i == 2:
                                        expec = [3.0]
                                    elif i == 3:
                                        expec = [4.0]
                            elif i==4 or i==5:
                                if dmethod=="c":
                                    expec = [2.5, 10.5]
                                else:
                                    if i == 4:
                                        expec = [4.0, 12.0]
                                    if i == 5:
                                        expec = [4.0, 14.0]
                            specsmooth(
                                imagename=imagename, outfile=outfile,
                                region=reg, function="h", axis=2, dmethod=dmethod
                            )
                            han.open(outfile)
                            got = han.getchunk().ravel()
                            self.assertTrue((got == expec).all())
                            han.done()
                    else:
                        dmethod="c"
                        if i == 2:
                            expec = [1.5, 2.5, 3.5]
                        elif i == 3:
                            expec = [1.5, 2.5, 5.5, 7.5]
                        elif i == 4:
                            expec = [1.5, 2.5, 5.5, 10.5, 13.5]
                        elif i == 5:
                            expec = [1.5, 2.5, 5.5, 10.5, 17.5, 21.5]
                        specsmooth(
                            imagename=imagename, outfile=outfile,
                            region=reg, function="h", axis=2, dmethod=""
                        )
                        han = image( )
                        han.open(outfile)
                        got = han.getchunk().ravel()
                        self.assertTrue((got == expec).all())
                        han.done()
        myia.done()

    def test_history(self):
        """Test history records are written"""
        imagename = "zz.im"
        myia = image( )
        myia.fromshape(imagename, [20,20,20])
        outfile = "zz_out.im"
        self.assertTrue(
            specsmooth(
                imagename=imagename, outfile=outfile,
                axis=2, function="h"
            )
        )
        myia.open(outfile)
        msgs = myia.history()
        myia.done()
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        teststr = "specsmooth"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
 
def suite():
    return [specsmooth_test]

if __name__ == '__main__':
    unittest.main()


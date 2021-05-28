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
# Test suite for the CASA tool method ia.modify()
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
# Test for the ia.modify() tool method
# </etymology>
#
# <synopsis>
# Test the ia.modify() tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_modify[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.modify() tool method to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import shutil
import unittest
import numpy
import os

try:
    from casatools import image as iatool
    from casatools import componentlist as cltool
    from casatools import quanta as qatool
    from casatools import measures as metool
    from casatools import ctsys
    ctsys_resolve = ctsys.resolve
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
        return os.path.join(dataPath,apath)

datapath = ctsys_resolve('unittest/ia_modify/')

class ia_modify_test(unittest.TestCase):
    
    def setUp(self):
        self.qa = qatool( )
        self.me = metool( )
        pass
    
    def tearDown(self):
        self.qa.done( )
        self.me.done( )
        pass
    
    def test_stretch(self):
        """ ia.histogram(): Test stretch parameter"""
        mycl = cltool()
        mycl.addcomponent(flux=1, dir=['J2000', '00:00:00.00', '00.00.00.0'])
        yy = iatool()
        mymask = "maskim"
        yy.fromshape(mymask, [200, 200, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [200,200,1,20]
        yy.fromshape("", shape)
        yy.addnoise()
        self.assertRaises(
            Exception,
            yy.modify, model=mycl.torecord(),
            mask=mymask + ">0", stretch=False
        )
        zz = yy.modify(
            model=mycl.torecord(), mask=mymask + ">0", stretch=True
        )
        self.assertTrue(zz and type(zz) == type(True))
        yy.done()
        mycl.done()
        
    def test_CAS5688(self):
        """verify output is the same after this performance fix"""
        myia = iatool()
        for i in [0, 1]:
            if i == 0:
                imagename = "CAS5688_1.im"
            elif i == 1:
                imagename = "CAS5688_2.im"
            myia.fromshape(imagename, [20,20,1,10])
            if i == 0:
                world = myia.toworld([4.4, 4.4])['numeric']
            elif i == 1:
                world = myia.toworld([4.51, 4.51])['numeric']
            myia.setbrightnessunit("Jy/pixel")
            v0 = self.qa.quantity(world[0], "arcmin")
            v1 = self.qa.quantity(world[1], "arcmin")
            dir = self.me.direction("J2000", v0, v1)
            mycl = cltool()
            mycl.addcomponent(
                [1,0, 0, 0], "Jy", dir=dir, shape="point",
                polarization="Stokes", spectrumtype="spectral index", index=2.5
            )
            myia.modify(model=mycl.torecord(), subtract=False)
            bb = myia.getchunk()
            myia.done()
            mycl.done()
            myia.open(datapath + os.sep + imagename)
            cc = myia.getchunk()
            myia.done()
            self.assertTrue((bb == cc).all())

    def test_history(self):
        """Test history is added"""
        mycl = cltool()
        mycl.addcomponent(flux=1, dir=['J2000', '00:00:00.00', '00.00.00.0'])
        myia = iatool()
        myia.fromshape("", [200, 200, 1, 1])
        self.assertTrue(
            myia.modify(model=mycl.torecord()), "Failed to run ia.modify"
        )
        msgs = myia.history()
        mycl.done()
        myia.done()
        self.assertTrue("ia.modify" in msgs[-2], "History not written")
        self.assertTrue("ia.modify" in msgs[-1], "History not written")

    def test_disk(self):
        """test disk gives the right flux, CAS-10887"""
        mycl = cltool()
        mycl.addcomponent(
            dir="J2000 0:00:00 0.00.00", flux=1.0, shape="disk", majoraxis="100arcmin",
            minoraxis="100arcmin", positionangle="0deg"
        )
        myia = iatool()
        myia.fromshape("", [101,101])
        myia.modify(mycl.torecord(), subtract=False)
        mycl.done()
        self.assertTrue(numpy.isclose(myia.statistics()['sum'][0], 1, 1e-2))
        myia.done()

def suite():
    return [ia_modify_test]

if __name__ == '__main__':
    unittest.main()

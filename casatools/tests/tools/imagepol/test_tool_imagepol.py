##########################################################################
# test_tool_imagepol.py
#
# Copyright (C) 2018
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.html
#
# Methods from the image tool tested in this script
# complexfraclinpol, complexlinpol, depolration, fourierrotationmeasure, fraclinpol
# fractotpol, linpolint, linpolposang, rotationmeasure, sigma, sigmadepolratio,
# sigmafraclinpol, sigmafractotpol, sigmalinpolint, sigmalinpolposang, totpolint
#
##########################################################################

import os
import shutil
import math
import unittest

from casatools import imagepol as potool
from casatools import image as iatool
from casatools import ctsys, table, constants
ctsys_resolve = ctsys.resolve

datapath = ctsys_resolve('unittest/imagepol/')
eq_beams = "pol_eq_beams.fits"
neq_beams = "pol_neq_beams.fits"

# Base class to copy and remove data only once for all classes
class Imagepol_base(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        shutil.copy(datapath + eq_beams, eq_beams)
        shutil.copy(datapath + neq_beams, neq_beams)

    @classmethod
    def tearDownClass(cls) -> None:
        os.remove(eq_beams)
        os.remove(neq_beams)

# Tests for imagepol.complexfraclinpol
class po_complexfraclinpol_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        shutil.rmtree('g')
        shutil.rmtree('hh')

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.complexfraclinpol("g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.complexfraclinpol, "hh")

# Tests for imagepol.complexlinpol
class po_complexlinpol_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        tb = table()
        self.assertEqual(len(tb.showcache()), 0)

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.complexlinpol('g'))
        mypo.done()
#        os.remove(eq_beams)
        shutil.rmtree('g')
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.complexlinpol, 'hh')
        mypo.done()
#        os.remove(neq_beams)
        shutil.rmtree('hh')

# Tests for imagepol.depolration
class po_depolratio_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()
#        shutil.copy(datapath + eq_beams, eq_beams)
#        shutil.copy(datapath + neq_beams, neq_beams)

    def tearDown(self):
        self.mypo.done()
#        os.remove(eq_beams)
#        os.remove(neq_beams)

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.depolratio(eq_beams))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.depolratio, neq_beams)
        self.assertRaises(Exception, mypo.depolratio, eq_beams)

# Tests for imagepol.fourierrotationmeasure
class po_fourierrotationmeasure_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertRaises(Exception, mypo.fourierrotationmeasure, amp="mm")
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.fourierrotationmeasure, amp="hh")

# Tests for imagepol.fraclinpol
class po_fraclinpol_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.fraclinpol(True))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.fraclinpol, "hh")

# Tests for imagepol.fractotpol
class po_fractotpol_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.fractotpol(True))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.fractotpol, "hh")

# Tests for imagepol.linpolint
class po_linpolint_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.linpolint(True))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.linpolint, "hh")

# Tests for imagepol.linpolposang
class po_linpolposang_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        shutil.rmtree('g')

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.linpolposang("g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.linpolposang, "hh")

# Tests for imagepol.rotationmeasure
class po_rotationmeasure_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        shutil.rmtree('g', ignore_errors=True)
        shutil.rmtree('hh', ignore_errors=True)
        shutil.rmtree('rm_input.im', ignore_errors=True)
        shutil.rmtree('pa0.im', ignore_errors=True)

    #        tb = table( )
    #        self.assertTrue(len(tb.showcache()) == 0)
    #        tb.done( )

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.rotationmeasure("g"))
        mypo.done()
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.rotationmeasure, "hh")
        mypo.done()

    def test_algorithm(self):
        """Test rotation measure computation algorithm"""
        myia = iatool()
        imagename = "rm_input.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        csys = myia.coordsys()
        incr = csys.increment()['numeric']
        incr[3] = 1000 * incr[3]
        csys.setincrement(incr)
        myia.setcoordsys(csys.torecord())
        pixvals = myia.getchunk()
        # U values all 1
        U = 1
        pixvals[:, :, 2, :] = U
        c = constants.c / 100
        RM = 9.6
        pa0deg = 22.5
        pa0 = pa0deg / 180 * math.pi
        for chan in range(myia.shape()[3]):
            freq = myia.toworld([0, 0, 0, chan])['numeric'][3]
            lam = c / freq
            Q = U / math.tan(2 * (pa0 + RM * lam * lam))
            pixvals[:, :, 1, chan] = Q
        myia.putchunk(pixvals)
        myia.done()
        mypo = self.mypo
        rmim = "rm.im"
        pa0im = "pa0.im"
        sigma = 10e-8
        mypo.open(imagename)
        mypo.rotationmeasure(rm=rmim, pa0=pa0im, sigma=sigma)
        mypo.done()
        myia.open(rmim)
        stats = myia.statistics(list=True, verbose=True)
        self.assertTrue((abs(stats['min'][0] - RM)) < 1e-4)
        self.assertTrue((abs(stats['max'][0] - RM)) < 1e-4)
        myia.done(remove=True)
        myia.open(pa0im)
        stats = myia.statistics(list=True, verbose=True)
        self.assertTrue((abs(stats['min'][0] - pa0deg)) < 1e-4)
        self.assertTrue((abs(stats['max'][0] - pa0deg)) < 1e-4)
        myia.done(remove=True)

# Tests for imagepol.sigma
class po_sigma_test(Imagepol_base):
    
    def setUp(self):
        self.mypo = potool()
#        shutil.copy(datapath + eq_beams, eq_beams)
#        shutil.copy(datapath + neq_beams, neq_beams)

    def tearDown(self):
        self.mypo.done()
#        os.remove(eq_beams)
#        os.remove(neq_beams)
    
    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        print(eq_beams)
        mypo.open(eq_beams)
        self.assertTrue(mypo.sigma( ))
        mypo.open(neq_beams)
        self.assertTrue(mypo.sigma( ))

# Tests for imagepol.sigmadepolratio
class po_sigmadepolratio_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.sigmadepolratio(eq_beams))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.sigmadepolratio, eq_beams)

# Tests for imagepol.sigmafraclinpol
class po_sigmafraclinpol_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        shutil.rmtree('g')

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.sigmafraclinpol(outfile="g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.sigmafraclinpol, "hh")

# Tests for imagepol.sigmafractotpol
class po_sigmafractotpol_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        shutil.rmtree("g")

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.sigmafractotpol(outfile="g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.sigmafractotpol, "hh")

# Tests for imagepol.sigmalinpolint
class po_sigmalinpolint_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.sigmalinpolint(outfile="g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.sigmalinpolint, "hh")

# Tests for imagepol.sigmalinpolposang
class po_sigmalinpolposang_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()
        shutil.rmtree('g')

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.sigmalinpolposang(outfile="g"))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.sigmalinpolposang, "hh")

# Tests for imagepol.totpolint
class po_totpolint_test(Imagepol_base):

    def setUp(self):
        self.mypo = potool()

    def tearDown(self):
        self.mypo.done()

    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.totpolint(True))
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.totpolint, "hh")

if __name__ == '__main__':
    unittest.main()

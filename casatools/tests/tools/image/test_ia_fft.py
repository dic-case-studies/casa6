##########################################################################
# test_ia_fft.py
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
# Test suite for the CASA tool method ia.fft
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
# Test for the ia.fft tool method
# </etymology>
#
# <synopsis>
# Test the ia.fft() method.
# </synopsis> 
#
# <motivation>
# To provide a test standard for the ia.fft tool method to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import shutil
import unittest

try:
    from casatools import image as iatool
    from casatools import regionmanager as rgtool
    from casatools import table, ctsys
    ctsys_resolve = ctsys.resolve
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
        return os.path.join(dataPath,apath)

datapath = ctsys_resolve('unittest/ia_fft/')

from casatools import image as iatool
from casatools import quanta as qatool
from casatools import table

class ia_fft_test(unittest.TestCase):
    
    def setUp(self):
        self.tb = table( )
        pass
    
    def tearDown(self):
        self.assertTrue(len(self.tb.showcache()) == 0)
        self.tb.done( )
    
    def test_stretch(self):
        """ ia.fft(): Test stretch parameter"""
        yy = iatool()
        mymask = "maskim"
        yy.fromshape(mymask, [20, 20, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [20,20,1,5]
        yy.fromshape("", shape)
        #yy.addnoise()
        self.assertRaises(
            Exception,
            yy.fft, real="real1.im", mask=mymask + ">0",
            stretch=False
        )
        real = 'real2.im'
        zz = yy.fft(
            real=real, mask=mymask + ">0", stretch=True
        )
        shutil.rmtree(mymask)
        shutil.rmtree(real)
        self.assertTrue(type(zz) == type(False))
        yy.done()
        
    def test_delta(self):
        """Test fft of delta function"""
        myia = iatool()
        for t in ['f', 'c']:
            myia.fromshape("", [100, 100], type=t)
            bb = myia.getchunk()
            bb[50, 50] = 1
            myia.putchunk(bb)
            real = "real.im" + t
            imag = "imag.im" + t
            amp = "amp.im" + t
            phase = "phase.im" + t
            _complex = "complex.im" + t
            myia.fft(
                real=real, imag=imag, amp=amp,
                phase=phase, complex=_complex
            )
            for im in [real, imag, amp, phase, _complex]:
                expec = 1
                if im == imag or im == phase:
                    expec = 0
                elif im == _complex:
                    expec = 1 + 0j
                myia.open(im)
                got = myia.getchunk()
                myia.done()
                self.assertTrue((got == expec).all())
                shutil.rmtree(im)

    def test_regression(self):
        """Was regression test in imagetest"""

        # Open test image (has sky coordinates)
        testname = 'unittest/ia_fft/test_image.im'
        myia = iatool()
        testim = iatool()
        testim.open(ctsys_resolve(testname))
        self.assertTrue(testim)
        testshape = testim.shape()
        self.assertTrue(len(testshape) == 3)
        rname = 'real_reg'
        iname = 'imag_reg'
        aname = 'amp_reg'
        pname = 'phase_reg'
        self.assertTrue(
            testim.fft(
                real=rname, imag=iname, phase=pname, amp=aname
            )
        )
        im1 = myia.newimage(rname)
        self.assertTrue(im1)
        im2 = myia.newimage(iname)
        self.assertTrue(im2)
        im3 = myia.newimage(aname)
        self.assertTrue(im3)
        im4 = myia.newimage(pname)
        self.assertTrue(im4)
        trc = testim.shape()
        trc[2] = 0
        a1 = im1.getchunk(trc=trc)
        a2 = im2.getchunk(trc=trc)
        a3 = im3.getchunk(trc=trc)
        a4 = im4.getchunk(trc=trc)

        from numpy.fft import fft2
        p = testim.getchunk(trc=trc)
        c = fft2(p)
        b1 = c.real
        b2 = c.imag
        b3 = abs(c)  # sqrt( real(x)^2 + imag(x)^2 )
        
        ok =im1.remove(True) and im2.remove(True) and im3.remove(True) and im4.remove(True)
        self.assertTrue(ok)
        #
        # FFT whole image
        #
        ndim = len(testim.shape())
        axes = list(range(ndim))
        ok = testim.fft(real=rname, imag=iname, phase=pname, amp=aname, axes=axes)
        self.assertTrue(ok)
        im1 = myia.newimage(rname)
        self.assertTrue(im1)
        im2 = myia.newimage(iname)
        self.assertTrue(im2)
        im3 = myia.newimage(aname)
        self.assertTrue(im3)
        im4 = myia.newimage(pname)
        self.assertTrue(im4)
        a1 = im1.getchunk()
        a2 = im2.getchunk()
        a3 = im3.getchunk()
        a4 = im4.getchunk()
    
        p = testim.getchunk()
        c = fft2(p)
        b1 = c.real
        b2 = c.imag
        b3 = abs(c)
       
        self.assertTrue(
            testim.done() and im1.done(remove=True)
            and im2.done(remove=True) and im3.done(remove=True)
            and im4.done(remove=True)
        )

    def test_history(self):
        """verify history writing"""
        myia = iatool()
        (real, imag, amp, phase, complx) = ("myreal.im", "myimag.im", "myamp.im", "myphase.im", "mycomplex.im")
        myia.fromshape("", [20,20])
        myia.fft(real=real, imag=imag, amp=amp, phase=phase, complex=complx)
        myia.done()
        for im in (real, imag, amp, phase, complx):
            myia.open(im)
            msgs = myia.history()
            myia.done()
            self.assertTrue("ia.fft" in msgs[-2])
            self.assertTrue("ia.fft" in msgs[-1])
            shutil.rmtree(im)

    def test_units(self):
        """
        CAS-13489: test output units

        The output phase image should have units of radians
        If the input image represents the image plane and has units of Jy/beam or Jy/pixel,
        then the ouptut (uv-plane) images should have units of Jy. If the input image has
        a beam, the beam should be copied to the output images.
        If the input image represents the uv-plane, then the brightness unit of the
        output (image plane) images should be Jy/beam or Jy/pixel, depending on if the
        input image has a beam or not.
        """
        qa = qatool()
        bmaj = qa.quantity('4arcmin')
        bmin = qa.quantity('3arcmin')
        bpa = qa.quantity('60deg')
        _ia = iatool()
        for bu in ('Jy/beam', 'Jy/pixel'):
            # create image-domain image
            _ia.fromshape("", [20, 20])
            self.assertTrue(_ia.setbrightnessunit('Jy/pixel'), 'Failed to set brightness unit')
            if bu == 'Jy/beam':
                self.assertTrue(
                    _ia.setrestoringbeam(major=bmaj, minor=bmin, pa=bpa),
                    'Failed to set restoring beam'
                )
            real = "real.im"
            imag = "imag.im"
            amp = "amp.im"
            phase = "phase.im"
            _complex = "complex.im"
            self.assertTrue(
                _ia.fft(
                    real=real, imag=imag, amp=amp,
                    phase=phase, complex=_complex
                ), 'ia.fft() failed'
            )
            _ia.done()
            for im in (real, imag, amp, phase, _complex):
                _ia.open(im)
                bunit = _ia.brightnessunit()
                beam = _ia.restoringbeam()
                _ia.done(remove=True)
                expec = 'Jy'
                if im == phase:
                    expec = 'rad'
                self.assertTrue(
                    bunit == expec, 'image ' + im + ' has unit ' + bunit
                    + ' but should be ' + expec
                )
                if bu == 'Jy/pixel':
                    self.assertTrue(beam == {}, 'this image should have no restoring beam')
                elif bu == 'Jy/beam':
                    self.assertTrue(beam['major'] == bmaj, 'Incorrect restoring beam')
                    self.assertTrue(beam['minor'] == bmin, 'Incorrect restoring beam')
                    self.assertTrue(beam['positionangle'] == bpa, 'Incorrect restoring beam')


def suite():
    return [ia_fft_test]

if __name__ == '__main__':
    unittest.main()

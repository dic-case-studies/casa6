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
# Test for the ia.fft() tool method
# </etymology>
#
# <synopsis>
# Test the ia.fft() method.
# </synopsis> 
#
# <example>
# This test can be run via
# PYTHONPATH=<my python path> python casatools/tests/tools/image/test_ia_fft.py
# </example>
#
# <motivation>
# To provide a test standard for the ia.fft() tool method to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import numpy as np
import shutil
import unittest

from casatools import image 
from casatools import regionmanager
from casatools import table, ctsys
datapath = ctsys.resolve('unittest/ia_fft/')

class ia_fft_test(unittest.TestCase):
    
    def setUp(self):
        self.tb = table( )
        pass
    
    def tearDown(self):
        data = [
            "amp.imc", "amp.imf", "amp_reg",
            "complex.imc", "complex.imf",
            "imag.imc", "imag.imf", "imag_reg",
            "maskim", "myamp.im", "mycomplex.im",
            "myimag.im", "myphase.im", "myreal.im",
            "phase.imc", "phase.imf", "phase_reg",
            "real2.im", "real.imc", "real.imf", "real_reg" ]
        for f in data:
            if os.path.exists(f):
                if os.path.isfile(f) or os.path.islink(f):
                    os.unlink(f)
                else:
                    shutil.rmtree(f)


        self.assertTrue(len(self.tb.showcache()) == 0)
        self.tb.done( )
    
    def test_stretch(self):
        """ ia.fft(): Test stretch parameter"""
        yy = image()
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
        zz = yy.fft(
            real="real2.im", mask=mymask + ">0", stretch=True
        )
        self.assertTrue(type(zz) == type(False))
        yy.done()
        shutil.rmtree('real2.im')
        shutil.rmtree(mymask)
        
    def test_delta(self):
        """Test fft of delta function"""
        myia = image()
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
                myia.done(remove=True)
                self.assertTrue((got == expec).all())
                
    def test_regression(self):
        """Was regression test in imagetest"""

        # Open test image (has sky coordinates)
        testname = 'unittest/ia_fft/test_image.im'
        myia = image()
        testim = image()
        testim.open(ctsys.resolve(testname))
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
       
        ok = (
                testim.done() and im1.done(remove=True)
                and im2.done(remove=True) and im3.done(remove=True)
                and im4.done(remove=True)
            )
        self.assertTrue(ok)

    def test_history(self):
        """verify history writing"""
        myia = image()
        (real, imag, amp, phase, complx) = ("myreal.im", "myimag.im", "myamp.im", "myphase.im", "mycomplex.im")
        myia.fromshape("", [20,20])
        myia.fft(real=real, imag=imag, amp=amp, phase=phase, complex=complx)
        myia.done()
        for im in (real, imag, amp, phase, complx):
            myia.open(im)
            msgs = myia.history()
            myia.done(remove=True)
            self.assertTrue("ia.fft" in msgs[-2])
            self.assertTrue("ia.fft" in msgs[-1])
    
    def test_new_inc(self):
        """verify CAS-13629 ouput cellsize fix"""
        myia = image()
        npix = 200
        myia.fromshape("", [npix, npix])
        myia.fft(real='real.im')
        myia.open('real.im')
        csys = myia.coordsys()
        myia.done(remove=True)
        cdelt = csys.increment()['numeric']
        csys.done()
        mye = 1/(npix/60*np.pi/180)
        expec = np.array([-mye, mye])
        print((cdelt - expec)/expec)
        self.assertTrue(np.allclose(cdelt, expec))

def suite():
    return [ia_fft_test]

if __name__ == '__main__':
    unittest.main()

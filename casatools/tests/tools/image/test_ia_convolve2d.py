##########################################################################
# test_ia_convolve2d.py
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
###########################################################################

from __future__ import absolute_import
from __future__ import print_function
import math
import numpy as np
import os
import shutil
import unittest

try:
    from casatools import ctsys, image, table, quanta, regionmanager
    _tb = table()
    _rg = regionmanager()
    _qa = quanta()
    ctsys_resolve = ctsys.resolve
    is_CASA6 = True
except ImportError:
    from tasks import *
    from taskinit import *
    import casac
    from __main__ import *
    _tb = tbtool()
    _rg = rgtool()
    image = iatool
    _qa = qatool()
    dataRoot = os.path.join(os.environ.get('CASAPATH').split()[0],'casatestdata/')
    def ctsys_resolve(apath):
        return os.path.join(dataRoot,apath)
    is_CASA6 = False

def _near(got, expected, tol):
    return _qa.le(
        _qa.div(_qa.abs(_qa.sub(got, expected)), expected), tol
    )

def make_gauss2d(shape, xfwhm, yfwhm):
    fac = 4*math.log(2)
    values = np.empty(shape, dtype=float)
    for i in range(shape[0]):
        x = shape[0]/2 - i
        for j in range(shape[1]):
            y = shape[1]/2 - j
            xfac = x*x*fac/(xfwhm*xfwhm)
            yfac = y*y*fac/(yfwhm*yfwhm)
            values[i, j] = math.exp(-(xfac + yfac));
    return values

def run_convolve2d(
    imagename, major, minor, pa, targetres,
    outfile, kernel="gauss", beam={}, overwrite=False
):
    myia = image()
    myia.open(imagename)
    res = myia.convolve2d(
        type=kernel,
        major=major, minor=minor, pa=pa,
        targetres=targetres, outfile=outfile,
        beam=beam, overwrite=overwrite
    )
    myia.done()
    res.done()

class ia_convolve2d_test(unittest.TestCase):
   
    def setUp(self):
        self.datapath = 'unittest/imsmooth/'
        self.imname = ''

    def tearDown(self):
        if self.imname:
            if os.path.isfile(self.imname):
                os.unlink(self.imname)
            else:
                shutil.rmtree(self.imname)

        self.assertTrue(len(_tb.showcache()) == 0, 'table cache is not empty')

    def _compare_beams(self, beam1, beam2):
        self.assertTrue(_near(beam1["major"], beam2["major"], 2e-5))
        self.assertTrue(_near(beam1["minor"], beam2["minor"], 2e-5))
        pa = []
        for b in [beam1, beam2]:
            if "positionangle" in b:
                pa.append(b["positionangle"])
            else:
                pa.append(b["pa"])

        diff = abs(
            _qa.sub(
                _qa.quantity(pa[0]), 
                _qa.quantity(pa[1])
            )["value"]
        )
        self.assertTrue(diff < 1e-5)
 
    def test_multibeam(self):
        """Test per plane beams"""
        myia = image()
        self.imname = "test_image2dconvolver_multibeam.im"
        shutil.copytree(ctsys_resolve(os.path.join(self.datapath, self.imname)), self.imname)
        myia.open(self.imname)
        major = "10arcmin"
        minor = "8arcmin"
        pa = "80deg"
        got = myia.convolve2d(axes=[0, 1], major=major, minor=minor, pa=pa)
        shape = myia.shape()
        for i in range(5):
            blc=[0, 0, i]
            trc=[shape[0]-1, shape[1]-1, i]
            reg = _rg.box(blc=blc, trc=trc)
            xx = myia.subimage(region=reg)
            exp = xx.convolve2d(axes=[0, 1], major=major, minor=minor, pa=pa)
            expbeam = exp.restoringbeam()
            gotbeam = got.restoringbeam(channel=i)
            for j in ["major", "minor", "positionangle"]:
                self.assertTrue(_near(gotbeam[j], expbeam[j], 2e-7))
            self.assertTrue(abs(got.getchunk(blc=blc, trc=trc) - exp.getchunk()).max() < 3e-5)
            exp.done()
            xx.done()
        myia.done()
        got.done()

    def test_targetres(self):
        """Test targetres parameter"""
        myia = image()
        self.imname = "tres1.im"
        myia.fromshape(self.imname, [100, 100])
        csys = myia.coordsys()
        csys.setunits(["arcsec", "arcsec"])
        csys.setincrement([-1, 1])
        myia.setcoordsys(csys.torecord())
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg")
        shape = myia.shape()
        values = make_gauss2d(shape, 3.0, 6.0)
        expected = make_gauss2d(shape, 5.0, 10.0)
        myia.putchunk(values)
        myia.done()
        emaj = _qa.quantity("10arcsec")
        emin = _qa.quantity("5arcsec")
        epa = _qa.quantity("0deg")
        for unit in ("Jy/beam", "K"):
            myia.open(self.imname)
            myia.setbrightnessunit(unit)
            myia.done()
            expected = make_gauss2d(shape, 5.0, 10.0)
            if (unit == "K"):
                expected *= 3.0*6.0/5.0/10.0
            # for code in (run_convolve2d, run_imsmooth):
            for targetres in [False, True]:
                if not targetres:
                    major = "8arcsec"
                    minor = "4arcsec"
                    pa = "0deg"
                    outfile = "tres1" + unit[0]
                else:
                    major = "10arcsec"
                    minor = "5arcsec"
                    pa = "0deg"
                    outfile = "tres2" + unit[0]
                run_convolve2d(
                    imagename=self.imname, kernel="gaussian",
                    major=major, minor=minor, pa=pa, targetres=targetres,
                    outfile=outfile
                )
                myia.open(outfile)
                gotbeam = myia.restoringbeam()
                gotvals = myia.getchunk()
                myia.done()
                shutil.rmtree(outfile)
                self._compare_beams(
                    gotbeam, {"major": emaj, "minor": emin, "pa": epa}
                )    
    
    def test_beam(self):
        """Test the beam parameter"""
        myia = image()
        self.imname = "tbeam1.im"
        myia.fromshape(self.imname, [100, 100])
        csys = myia.coordsys()
        csys.setunits(["arcsec", "arcsec"])
        csys.setincrement([1, 1])
        myia.setcoordsys(csys.torecord())
        myia.setbrightnessunit("Jy/beam")
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg")
        shape = myia.shape()
        myia.putchunk(make_gauss2d(shape, 3.0, 6.0))
        expected = make_gauss2d(shape, 5.0, 10.0)
        for beam in [
            {"major": "8arcsec", "minor": "4arcsec", "pa": "0deg"},
            {
                "major": {"unit": "arcsec", "value": 8},
                "minor": {"unit": "arcsec", "value": 4},
                "pa": {"unit": "deg", "value": 0},
            }
        ]:
            outfile = 'convolve2d'
            x = run_convolve2d(
                imagename=self.imname, major="", minor="", pa="",
                beam=beam, outfile=outfile, targetres=False,
                overwrite=True
            )
            if type(x) == type(myia):
                x.done()
            myia.open(outfile)
            maxdiff = (abs(myia.getchunk()-expected)).max()
            self.assertTrue(maxdiff < 1e-6) 
            myia.done()
            shutil.rmtree(outfile)

    def test_history(self):
        """Test that history is written"""
        myia = image()
        self.imname = "zz.im"
        myia.fromshape(self.imname, [20,20])
        major = "2arcmin"
        minor = "2arcmin"
        pa = "0deg"
        bb = myia.convolve2d("", major=major,  minor=minor, pa=pa)        
        myia.done()
        msgs = bb.history()
        bb.done()
        teststr = "ia.convolve2d"
        self.assertTrue(teststr in msgs[-4], "'" + teststr + "' not found")     
        self.assertTrue(teststr in msgs[-3], "'" + teststr + "' not found")
        
    def test_stretch(self):
        """ ia.convolve2d(): Test stretch parameter"""
        yy = image()
        self.imname = "maskim"
        yy.fromshape(self.imname, [200, 200, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [200,200,1,20]
        yy.fromshape("", shape)
        yy.addnoise()
        self.assertRaises(
            Exception,
            yy.convolve2d, "", [0,1], "gaussian", "4arcmin", 
            "4arcmin", "0deg", mask=self.imname + ">0", stretch=False
        )
        zz = yy.convolve2d(
            "", [0,1], "gaussian", "4arcmin", "4arcmin", "0deg",
            mask=self.imname + ">0", stretch=True
        )
        self.assertTrue(type(zz) == type(yy))
        yy.done()
        zz.done()
        
    def test_precision(self):
        """Test images of various precisions"""
        yy = image()
        for mytype in ['d', 'c', 'f', 'cd']:
            yy.fromshape("", [20,  20, 1], type=mytype)
            yy.addnoise()
            if mytype == 'f' or mytype == 'd':
                zz = yy.convolve2d(
                    "", [0,1], "gaussian", "4arcmin", "4arcmin", "0deg"
                )
                self.assertTrue(zz)
                zz.done()
            else:
                self.assertRaises(
                    Exception, yy.convolve2d, "", [0,1], "gaussian",
                    "4arcmin", "4arcmin", "0deg"
                )
            yy.done()

    def test_copying_of_input_mask(self):
        """CAS-12904: copy input mask to output image"""
        self.imname = 'orig.im'
        yy = image()
        yy.fromshape(self.imname, [100, 100, 3])
        pix = yy.getchunk()
        for i in range(3):
            pix[:, :, i] = i
        yy.putchunk(pix)
        subi = yy.subimage("", mask=self.imname + '>0')
        yy.done()
        for i in range(3):
            reg = _rg.box([0, 0, i], [99, 99, i])
            npts = subi.statistics(region=reg)['npts']
            expec = 0 if i == 0 else 1
            # shows mask was created correctly
            self.assertEqual(npts.size, expec, 'wrong length npts array')
            if i>0:
                self.assertEqual(npts[0], 10000, 'wrong number of pts')
        conv = subi.convolve2d(
            major='4arcmin', minor='4arcmin', pa='0deg', mask=self.imname + '<2'
        )
        subi.done()
        for i in range(3):
            reg = _rg.box([0, 0, i], [99, 99, i])
            npts = conv.statistics(region=reg)['npts']
            expec = 1 if i == 1 else 0
            # shows mask was copied correctly
            self.assertEqual(npts.size, expec, 'wrong length npts array')
            if i==1:
                self.assertEqual(npts[0], 10000, 'wrong number of pts')
        conv.done()     
 
def suite():
    return [ia_convolve2d_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()


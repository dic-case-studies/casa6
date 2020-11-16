##########################################################################
# test_rg_fromtextfile.py
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
###########################################################################
import os
import shutil
import unittest
import numpy

try:
    from casatools import regionmanager as rgtool
    from casatools import image as iatool
    from casatools import quanta
    from casatools import ctsys
    from casatools import imagemetadata
    ctsys_resolve = ctsys.resolve
    _qa = quanta()
    _imd = imagemetadata()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    _qa = qatool()
    _imd = imdtool()
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
        return os.path.join(dataPath,apath)    

image = "imregion.fits"
text1 = "goodfile1.txt"
res1 = "res1.rgn"
cas_3258t = "CAS-3258.txt"
cas_3258r = "CAS-3258.rgn"
cas_3259t = "CAS-3259.txt"
cas_3259r = "CAS-3259.rgn"
cas_3260t = "CAS-3260.txt"
cas_3260r = "CAS-3260.rgn"
icrs_image = "icrs.im"

datapath = ctsys_resolve('unittest/rgtool/')

def deep_equality(a, b):
    if (type(a) != type(b)):
        print("types don't match, a is a " + str(type(a)) + " b is a " + str(type(b)))
        return False
    if (type(a) == dict):
        if (a.keys() != b.keys()):
            print("keys don't match, a is " + str(a.keys()) + " b is " + str(b.keys()))
            return False
        for k in a.keys():
            if (
                k == "telescope" or k == "observer"
                or k == "telescopeposition"
            ):
                continue
            elif (not deep_equality(a[k], b[k])):
                print("dictionary member inequality a[" + str(k) \
                    + "] is " + str(a[k]) + " b[" + str(k) + "] is " + str(b[k]))
                return False
        return True
    if (type(a) == float):
        if not (a == b or abs((a-b)/a) <= 1e-6):
            print("float mismatch, a is " + str(a) + ", b is " + str(b))
        return a == b or abs((a-b)/a) <= 1e-6
    if (type(a) == numpy.ndarray):
        if (a.shape != b.shape):
            print("shape mismatch a is " + str(a.shape) + " b is " + str(b.shape))
            return False
        x = a.tolist()
        y = b.tolist()
        for i in range(len(x)):
            if (not deep_equality(x[i], y[i])):
                print("array element mismatch, x is " + str(x[i]) + " y is " + str(y[i]))
                return False
        return True
    return a == b

class rg_fromtextfile_test(unittest.TestCase):
    
    _fixtures = [
        image, text1, res1, cas_3258t, cas_3258r, cas_3259t, cas_3259r,
        cas_3260t, cas_3260r
    ]

    _created = [icrs_image]
    
    def setUp(self):
        for im in self._fixtures:
            shutil.copy(datapath + im, im)
        self.ia = iatool()
        self.rg = rgtool()
    
    def tearDown(self):
        for im in self._fixtures + self._created:
            if os.path.exists(im):
                if os.path.isdir(im):
                    shutil.rmtree(im)
                else:
                    os.remove(im)
        self.ia.done()
        del self.ia
        self.rg.done()
        del self.rg

    def _testit(self, text, rgn):
        csys = self.ia.coordsys().torecord()
        shape = self.ia.shape()
        got = self.rg.fromtextfile(text, shape, csys)
        expected = self.rg.fromfiletorecord(rgn)
        expected['comment'] = ""
        self.assertTrue(deep_equality(got, expected))
        
        f = open(text, 'r')
        text = f.read()
        got = self.rg.fromtext(text, shape, csys)
        self.assertTrue(deep_equality(got, expected))

    def test_exceptions(self):
        """test exception cases"""

        # bad file
        self.assertRaises(Exception, self.rg.fromtextfile, "blah", {}, [1,1])
        # coordsys not set
        self.assertRaises(Exception, self.rg.fromtextfile, text1, {}, [1,1])

    def test_read(self):
        """Read test"""
        self.ia.open(image)
        self._testit(text1, res1)
        
    def test_CAS_3258(self):
        """Verify fix to CAS-3258"""
        self.ia.fromshape("", [250,250])
        self._testit(cas_3258t, cas_3258r)
        
    def test_CAS_3259(self):
        """Verify fix to CAS-3259"""
        self.ia.fromshape("", [250,250])
        self._testit(cas_3259t, cas_3259r)
        
    def test_CAS_3260(self):
        """Verify fix to CAS-3260"""
        self.ia.fromshape("", [250,250])
        self._testit(cas_3260t, cas_3260r)
        
    def test_CAS_4415(self):
        """Verify CAS-4415 (parser did not properly handle frquency decreasing with pixel number)"""
        shape = [50, 50, 10]
        self.ia.fromshape("", shape)
        csys = self.ia.coordsys()
        increment = csys.increment()["numeric"]
        increment[2] = -increment[2]
        csys.setincrement(increment)
        self.ia.setcoordsys(csys.torecord())
        zz = self.rg.fromtext(
            "circle[[20pix,20pix],6pix],range=[1pix,3pix]",
            shape, csys.torecord()
        )
        self.assertTrue(len(zz.keys()) > 0)

    def test_CAS_4425(self):
        """ Verify CAS-4425 (pixel units now accounted for in range and no units throws exception)"""
        shape = [100, 100, 80]
        self.ia.fromshape("", shape)
        csys = self.ia.coordsys()
        zz = self.rg.fromtext("box[[30pix, 30pix], [39pix, 39pix]], range=[55pix,59pix]", shape, csys.torecord())
        self.assertTrue(self.ia.statistics(region=zz)["npts"] == 500)
        zz = self.rg.fromtext("box[[30pix, 30pix], [39pix, 39pix]], range=[59pix,55pix]", shape, csys.torecord())
        self.assertTrue(self.ia.statistics(region=zz)["npts"] == 500)
        self.assertRaises(
            Exception, self.rg.fromtext, "box[[30pix, 30pix], [39pix, 39pix]], range=[59,55]",
            shape, csys.torecord()
        )
         
    def test_rectangle_rotation(self):
        """Test rectangle region is preserved under coordinate frame switch"""
        self.ia.fromshape("",[200, 200])
        csys = self.ia.coordsys()
        # rectangular box
        xx = self.rg.fromtext(
            "box[[5834.23813221arcmin, -3676.92506701arcmin],[5729.75600494arcmin, -3545.36602909arcmin]] coord=GALACTIC",
            csys=csys.torecord(), shape=self.ia.shape()
        )
        zz = self.ia.subimage("", region=xx)
        got = zz.getchunk(getmask=True)
        self.ia.open(datapath + "rect_rot.im")
        expec = self.ia.getchunk(getmask=True)
        self.assertTrue((got == expec).all())
        zz.done()
        
        #center box
        self.ia.fromshape("",[200, 200])
        csys = self.ia.coordsys()
        yval = "-3611.1455480499999arcmin"
        xwidth = _qa.tos(_qa.mul(_qa.quantity("104.48212727000009arcmin"),_qa.cos(yval)))
        xx = self.rg.fromtext(
            "centerbox[[5781.9970685749995arcmin, " + yval + "],[" + xwidth + ", 131.55903791999981arcmin]] coord=GALACTIC",
            csys=csys.torecord(), shape=self.ia.shape()
        )
        zz = self.ia.subimage("", region=xx)
        got = zz.getchunk(getmask=True)
        self.ia.open(datapath + "rect_rot2.im")
        expec = self.ia.getchunk(getmask=True)
        self.assertTrue((got == expec).all())
        
        zz.done()
        self.ia.done()
        
    def test_rotbox(self):
        """Test rotbox when specified in pixels (CAS-5723)"""
        self.ia.fromshape("",[200,200])
        reg = self.rg.fromtext(
            "rotbox [ [ 60pix , 50pix ] , [ 30pix , 30pix ] , 30deg ]",
            csys=self.ia.coordsys().torecord(),shape=self.ia.shape()
        )
        self.assertTrue(self.ia.statistics(region=reg)['npts'] == 901)
        csys = self.ia.coordsys()
        csys.setreferencevalue([800,70*60])
        self.ia.setcoordsys(csys.torecord())
        reg = self.rg.fromtext(
            "rotbox [ [ 60pix , 50pix ] , [ 30pix , 30pix ] , 30deg ]",
            csys=self.ia.coordsys().torecord(),shape=self.ia.shape()
        )
        self.assertTrue(self.ia.statistics(region=reg)['npts'] == 901)
        
    def test_ellipse(self):
        """Test ellipse for image in GALACTIC and file in J2000"""
        self.ia.open(datapath + "gal.im")
        reg = self.rg.fromtextfile(
            datapath + "testEllipse90deg.crtf",
            csys = self.ia.coordsys().torecord(),
            shape=self.ia.shape()
        )
        subi = self.ia.subimage("", region=reg)
        self.ia.open(datapath + "galwj2kellipse.im")
        expec = self.ia.getchunk(getmask=True)
        self.ia.done()
        got = subi.getchunk(getmask=True)
        subi.done()
        self.assertTrue((got == expec).all())
        
    def test_1000(self):
        """Test a large number of regions, CAS-7405"""
        self.ia.open(datapath + "1000regtest.im")
        self.assertTrue(self.ia.statistics()['npts'][0] == 331*331)
        self.assertTrue(self.ia.statistics(region=datapath + "1000circles.txt")['npts'][0] == 13679)
        self.ia.done()
        
    def test_CAS_8072(self):
        """Verify rest frequency precision issue has been fixed"""
        self.ia.fromshape("",[20,20,200])
        self.ia.addnoise()
        reg = reg = self.rg.fromtext(
            "box [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s], restfreq=1.42040575e+09Hz",
            csys = self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        reg1 = self.rg.fromtext(
            "box [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys = self.ia.coordsys().torecord(), shape=self.ia.shape()
        )        
        # no comma delimiters should throw exception
        self.assertRaises(
            Exception, self.rg.fromtext,
            "global coord=B1950 frame=LSRK veltype=RADIO restfreq=1.42040575e+09Hz"
            + "\nbox [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys = self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        reg3 = self.rg.fromtext(
            "global coord=J2000, frame=LSRK, veltype=RADIO, restfreq=1.42040575e+09Hz"
            + "\nbox [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys = self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        # different global rest freq
        reg4 = self.rg.fromtext(
            "global coord=J2000, frame=LSRK, veltype=RADIO, restfreq=1.42050575e+09Hz"
            + "\nbox [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys = self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        stats0 = self.ia.statistics(region=reg)
        stats1 = self.ia.statistics(region=reg1)
        stats3 = self.ia.statistics(region=reg3)
        stats4 = self.ia.statistics(region=reg4)

        self.ia.done()
        for k in ('maxpos', 'minpos'):
            self.assertTrue((stats0[k] == stats1[k]).all())
            self.assertTrue((stats0[k] == stats3[k]).all())
            # different rest freq used, so the positions will not
            # be the same
            self.assertTrue((stats0[k] != stats4[k]).any())

    def test_ICRS(self):
        """
        CAS-13074, verify that coord=ICRS works correctly

        1. Create a 100x100 image using the default coordinate system provided by ia.shape()

        2. Modify the coordinate system of the image from J2000 to GALACTIC ref frame,
           since the difference between J2000 and ICRS seems to be only about 10 marcsec,
           so we need to use a coordinate system where the values in the two systems differ
           more to get a convincing test.

        3. Set all pixels to 0, except the reference pixel 50, 50 which is set to 1.

        4. Create a CRTF region string using coords='ICRS' and a box that is centered at
           GALACTIC long=0, lat=0 (the ref direction in the image) but using ICRS coords
           of course. The hardcoded ICRS coords were determined using me.measure() to
           convert from GALACTIC to ICRS.

        5. Run ia.statistics() using region=the CRTF text string previously created.

        6. Test that the max pixel value found in the region is 1. This indicates that
           the region was applied correctly.

        """
        self.ia.fromshape(icrs_image, shape=[100, 100])
        csys = self.ia.coordsys()
        csys.setconversiontype("GALACTIC")
        self.ia.setcoordsys(csys.torecord())
        csys.done()
        pix = self.ia.getchunk()
        pix[:] = 0
        pix[50, 50] = 1
        self.ia.putchunk(pix)
        self.ia.done()
        # muck with the image's coordinate system more
        self.assertTrue(
            _imd.open(icrs_image), 'Unable to open imagemetadata object'
        )
        self.assertTrue(
            _imd.set('equinox', 'GALACTIC'), 'Failed to set equinox'
        )
        self.assertTrue(
            _imd.set('ctype1', 'LONGITUDE'), 'Failed to set ctype1'
        )
        self.assertTrue(
            _imd.set('ctype2', 'LATITUDE'), 'Failed to set ctype2'
        )
        _imd.done()
        icrs = 'box[[-1.63412rad, -0.50561rad], [-1.63296rad, -0.50445rad]] coord=ICRS\n'
        rg = rgtool()
        self.assertTrue(self.ia.open(icrs_image), 'Failed to open image')
        reg = rg.fromtext(icrs, csys=self.ia.coordsys().torecord(), shape=self.ia.shape())
        rg.done()
        stats = self.ia.statistics(region=reg)
        self.ia.done() 
        self.assertEqual(stats['max'][0], 1, 'Incorrect value for max')

def suite():
    return [rg_fromtextfile_test]

if __name__ == '__main__':
    unittest.main()

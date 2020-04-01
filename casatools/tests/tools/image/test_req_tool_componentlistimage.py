##########################################################################
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
# Test suite for componentlistimage
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
# Test suite for componentlistimage
# </etymology>
#
# <synopsis>
# Test suite for componentlistimage
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_componentlistimage[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for componentlistimage support to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
from __main__ import *
import numpy

CASA6 = False
try:
    import casatools
    # from casatasks import imhistory, casalog
    import sys
    import os
    myia = casatools.image()
    mycl = casatools.componentlist()
    csys = casatools.coordsys()
    myqa = casatools.quanta()
    sys.path.append(os.path.abspath(os.path.dirname(__file__)))
    CASA6 = True
except ImportError:
    import sys
    import os
    from __main__ import default
    from tasks import *
    from taskinit import *
    myia = iatool()
    mycl = cltool()
    csys = cstool()
    myqa = qatool()

import unittest
import shutil

if CASA6:
    estimates_file = casatools.ctsys.resolve('text/2gauss_estimates.txt')
    climage = casatools.ctsys.resolve('image/simple_cl.im')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        estimates_file = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/text/2gauss_estimates.txt'
        climage = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/simple_cl.im'
    else:
        estimates_file = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/text/2gauss_estimates.txt'
        climage = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/simple_cl.im'

class componentlistimage_test(unittest.TestCase):
    
    def setUp(self):
        self._myia = myia
        self._mycl = mycl
    
    def tearDown(self):
        self._myia.done()
        self._mycl.done()
    
    def test_fromcomponentlist(self):
        """Test ia.fromcomponentlist() functionality"""
        mycl = self._mycl
        myia = self._myia
        flux = [1, 2, 3, 4]
        direction = ['J2000', '00:00:00.00', '00.00.00.0']
        pt = "point"
        mycl.addcomponent(flux=flux, dir=direction, shape=pt)
        
        shape = [5, 5]
        self.assertTrue(myia.fromcomplist("", shape=shape, cl=mycl.torecord()))
        
        shape = [5, 5, 5]
        self.assertTrue(myia.fromcomplist("", shape=shape, cl=mycl.torecord()))
        
        shape = [5, 5, 4, 5]
        self.assertTrue(myia.fromcomplist("", shape=shape, cl=mycl.torecord()))
        
        imagename = "1ptsource.im"
        self.assertTrue(myia.fromcomplist(imagename, shape=shape, cl=mycl.torecord()))
        self.assertTrue(myia.open(imagename))
        
    def test_vals(self):
        """Test valid pixel values"""
        mycl = self._mycl
        myia = self._myia
        flux = [1, 2, 3, 4]
        direction = ['J2000', '00:00:00.00', '00.00.00.0']
        pt = "point"
        mycl.addcomponent(flux=flux, dir=direction, shape=pt)
        
        shape = [5, 5, 4, 1]
        self.assertTrue(myia.fromcomplist("", shape=shape, cl=mycl.torecord()))
        vals = myia.getchunk()
        for x in range(5):
            for y in range(5):
                for s in range(4):
                    if x == 2 and y == 2:
                        expec = flux[s]
                    else:
                        expec = 0
                    self.assertEqual(vals[x, y, s, 0], expec)
        mycsys = myia.coordsys()
        myia.done()
 
        # shuffle the stokes
        mycsys.setstokes("U V Q I")
        stokestoflux = [2, 3, 1, 0]
        
        self.assertTrue(
            myia.fromcomplist(
                "", shape=shape, cl=mycl.torecord(), csys=mycsys.torecord()
            )
        )
        vals = myia.getchunk()
        for x in range(5):
            for y in range(5):
                for s in range(4):
                    if x == 2 and y == 2:
                        expec = flux[stokestoflux[s]]
                    else:
                        expec = 0
                    self.assertEqual(vals[x, y, s, 0], expec)
        myia.done()
        mycl.done()
        major = "5arcmin"
        minor = "4arcmin"
        pa = "0deg"
        gauss = "Gaussian"
        mycl.addcomponent(
            flux=flux, dir=direction, majoraxis=major,
            minoraxis=minor, positionangle=pa, shape=gauss
        )
        shape = [30, 30, 4, 1]
        mycsys.setreferencepixel([15, 15, 0, 0])
        self.assertTrue(
            myia.fromcomplist(
                "", shape=shape, cl=mycl.torecord(), csys=mycsys.torecord()
            )
        )
        stats = myia.statistics(axes=[0, 1, 3])
        myia.done()
        for i in range(4):
            got = stats['sum'][i]
            expec = flux[stokestoflux[i]]
            self.assertTrue(
                numpy.isclose(got, expec),
                'i=' + str(i) + ' got=' + str(got) + ' expec=' + str(expec)
            )
        mycl.done()
        
    def test_gaussian(self):
        """Test gaussian produces correct results"""
        mycl = self._mycl
        myia = self._myia
        flux = [1, 2, 3, 4]
        gauss = "Gaussian"
        major = "5arcmin"
        minor = "4arcmin"
        pa = "0deg"
        stokes = ["I", "Q", "U", "V"]
        mycsys = csys.newcoordsys(
            direction=True,spectral=True, stokes=stokes
        )
        mycsys.setreferencepixel([15, 15, 0, 0])
        shape = [30, 30, 4, 1]
        dir0 = ['J2000', '00:00:00.00', '00.00.00.0']
        dir1 = ['J2000', '00:00:10.00', '-00.04.18']
        expecra = [0, 10]
        expecdec = [0, -4.3]
        j = 0
        tol = 1e-6
        for mydir in [dir0, dir1]:
            mycl.addcomponent(
                flux=flux, dir=mydir, majoraxis=major,
                minoraxis=minor, positionangle=pa, shape=gauss
            )
            self.assertTrue(
                myia.fromcomplist(
                    "", shape=shape, cl=mycl.torecord(), csys=mycsys.torecord()
                )
            )
            mycl.done()
            i = 0
            for s in stokes:
                res = myia.fitcomponents(stokes=s)
                mycl.fromrecord(res['results'])
                gotdir = mycl.getrefdir(0)
                rainsec = myqa.convert(gotdir['m0'], 's')['value']
                decinamin = myqa.convert(gotdir['m1'], 'arcmin')['value']
                self.assertTrue(numpy.isclose(rainsec, expecra[j], tol))
                self.assertTrue(numpy.isclose(decinamin, expecdec[j], tol))
                self.assertEqual(gotdir['refer'], "J2000")
                self.assertTrue(numpy.isclose(mycl.getfluxvalue(0)[i], flux[i]))
                mycl.done()
                i += 1
            myia.done()
            mycl.done()
            j += 1
            
        # try two gaussians simultaneously
        for mydir in [dir0, dir1]:
            mycl.addcomponent(
                flux=flux, dir=mydir, majoraxis=major,
                minoraxis=minor, positionangle=pa, shape=gauss
            )
        self.assertEqual(mycl.length(), 2)
        self.assertTrue(
            myia.fromcomplist(
                "", shape=shape, cl=mycl.torecord(), csys=mycsys.torecord()
            )
        )
        mycl.done()
        k = 0
        atol = 1e-5
        for s in stokes:
            res = myia.fitcomponents(stokes=s, estimates=estimates_file)
            mycl.fromrecord(res['results'])
            self.assertEqual(mycl.length(), 2)
            for i in [0, 1]:
                gotdir = mycl.getrefdir(i)
                rainsec = myqa.convert(gotdir['m0'], 's')['value']
                decinamin = myqa.convert(gotdir['m1'], 'arcmin')['value']
                self.assertTrue(
                    numpy.isclose(rainsec, expecra[i], rtol=0, atol=atol),
                    "got: " + str(rainsec) + " expec: " + str(expecra[i])
                )
                self.assertTrue(numpy.isclose(decinamin, expecdec[i], rtol=0, atol=atol))
                self.assertEqual(gotdir['refer'], "J2000")
                self.assertTrue(numpy.isclose(mycl.getfluxvalue(0)[k], flux[k]))
            mycl.done()
            k += 1
        myia.done()
        
    def test_history(self):
        """verify history writing"""
        mycl = self._mycl
        myia = self._myia
        flux = [1, 0, 0, 0]
        dir0 = ['J2000', '00:00:00.00', '00.00.00.0']
        pt = "point"
        mycl.addcomponent(flux=flux, dir=dir0,shape=pt)
        shape = [20, 20]
        myia.fromcomplist("", shape=shape, cl=mycl.torecord())
        mycl.done()
        msgs = myia.history()
        myia.done()
        teststr = "ia.fromcomplist"
        self.assertTrue(teststr in msgs[-2])
        self.assertTrue(teststr in msgs[-1])
        
    def test_multi_points_same_pixel(self):
        """Test that multiple point sources at the same pixel produce the correct result"""
        mycl = self._mycl
        myia = self._myia
        flux = [1, 0, 0, 0]
        dir0 = ['J2000', '00:00:00.00', '00.00.00.0']
        pt = "point"
        mycl.addcomponent(flux=flux, dir=dir0,shape=pt)
        mycl.addcomponent(flux=flux, dir=dir0,shape=pt)
        shape = [20, 20]
        clname = "jj.cl"
        mycl.rename(clname)
        myia.fromcomplist("", shape=shape, cl=mycl.torecord())
        mycl.done()
        stats = myia.statistics()
        myia.done()
        self.assertEqual(stats['max'], 2)
        # test opening by name
        myia.fromcomplist("", shape=shape, cl=clname)
        stats = myia.statistics()
        myia.done()
        self.assertEqual(stats['max'], 2)
        
    def test_mask(self):
        """Test support for masks"""
        mycl = self._mycl
        myia = self._myia
        flux0 = [1, 0, 0, 0]
        dir0 = ['J2000', '00:00:00.00', '00.00.00.0']
        flux1 = [2, 0, 0, 0]
        dir1 = ['J2000', '00:00:00.00', '00.05.00.0']
        pt = "point"
        mycl.addcomponent(flux=flux0, dir=dir0,shape=pt)
        mycl.addcomponent(flux=flux1, dir=dir1,shape=pt)
        shape = [20, 20]
        imagename = "jk.im"
        myia.fromcomplist(outfile=imagename, shape=shape, cl=mycl.torecord())
        mycl.done()
        stats = myia.statistics()
        self.assertEqual(stats['max'], 2)
        self.assertEqual(stats['sum'], 3)
        myia.calcmask(imagename + " > 1")
        stats = myia.statistics()
        self.assertEqual(stats['max'], 2)
        self.assertEqual(stats['sum'], 2)
        self.assertEqual(stats['npts'], 1)
        myia.calcmask(imagename + " < 2")
        myia.maskhandler("set", "mask1")
        stats = myia.statistics()
        self.assertEqual(stats['max'], 1)
        self.assertEqual(stats['sum'], 1)
        self.assertEqual(stats['npts'], 399)
        myia.maskhandler("set", "")
        stats = myia.statistics()
        self.assertEqual(stats['max'], 2)
        self.assertEqual(stats['sum'], 3)
        self.assertEqual(stats['npts'], 400)
        stats = myia.statistics(mask=imagename + " > 0")
        self.assertEqual(stats['max'], 2)
        self.assertEqual(stats['sum'], 3)
        self.assertEqual(stats['npts'], 2)
        
    def test_fromimage(self):
        """Test fromimage() supports reading from a componentlist image"""
        myia = self._myia
        infile = 'simple_cl.im'
        outfile = "akd.im"
        shutil.copytree(climage, infile)
        self.assertTrue(myia.fromimage(outfile=outfile, infile=infile))
        bb = myia.getchunk()
        myia.done()
        myia.open(climage)
        cc = myia.getchunk()
        myia.done()
        self.assertTrue((bb == cc).all())
        
def suite():
    return [componentlistimage_test]

if __name__ == '__main__':
    unittest.main()


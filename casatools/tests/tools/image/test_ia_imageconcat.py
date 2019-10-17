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
# Test suite for the CASA ia.imageconcat() tool method
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
#   <li> <linkto class="image:description">ia.imageconcat()</linkto> 
# </ul>
# </prerequisite>
#
# <etymology>
# ia_imageconcat_test stands for ia.imageconcat test
# </etymology>
#
# <synopsis>
# test_ia_imageconcat.py is a Python script that tests the correctness
# of the ia.imageconcat() tool method.
# </synopsis> 
#
# <example>
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_imageconcat[test1,test2,...]
# </example>
#
# <motivation>
# To provide a test standard to the ia.imageconcat() method to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import os
import shutil
import unittest

from casatools import image as iatool
from casatools import regionmanager
from casatools import quanta

class ia_imageconcat_test(unittest.TestCase):
   
    def make_images(self):
        myia = self._myia
        myia.fromshape("", [1, 1, 5]) 
        names = []
        rg = self.rg
        for i in range(5):
            name = "chan_" + str(i)
            names.append(name)
            subi = myia.subimage(
                name, region=rg.box([0, 0, i], [0, 0,i]), overwrite=True
            )
            got = subi.toworld([0 ,0, 0])['numeric'][2]
            expec = myia.toworld([0 ,0, i])['numeric'][2]
            self.assertTrue(got == expec)
            subi.done()
        return (names, myia)
 
    def setUp(self):
        self._myia = iatool()
        self.rg = regionmanager( )
        self.qa = quanta( )

    def tearDown(self):
        self._myia.done()


    def test_multibeam(self):
        """Test concatenating images with different beams"""
        myia = self._myia
        shape = [4, 4, 20]
        myia.fromshape("", shape)
        blc1=[0, 0, 0]
        trc1=[shape[0]-1, shape[1]-1, shape[2]/2-1]
        rg1 = self.rg.box(blc=blc1, trc=trc1)
        im1 = "image1.im"
        sub1 = myia.subimage(im1, region=rg1)
        im2 = "image2.im"
        blc2 = [0, 0, trc1[2]+1]
        trc2 = [shape[0]-1, shape[1]-1, shape[2]-1]
        rg2 = self.rg.box(blc=blc2, trc=trc2)
        sub2 = myia.subimage(im2, region=rg2)
        major = self.qa.quantity("3arcmin")
        minor = self.qa.quantity("2arcmin")
        pa = self.qa.quantity("0deg")
        major2 = self.qa.quantity("4arcmin")
        minor2 = self.qa.quantity("3arcmin")
        pa2 = self.qa.quantity("10deg")
        major3 = self.qa.quantity("5arcmin")
        minor3 = self.qa.quantity("4arcmin")
        pa3 = self.qa.quantity("20deg")

        # first image has no beam while second does
        sub1.setbrightnessunit("Jy/pixel") 
        sub2.setrestoringbeam(major=major, minor=minor, pa=pa)
        sub2.setbrightnessunit("Jy/beam")
        self.assertRaises(Exception, myia.imageconcat, "", [im1, im2])
        concat = myia.imageconcat("", [im1, im2], relax=True)
        self.assertTrue((concat.shape() == shape).all())
        
        # first image has a single beam, second has per plane beams
        sub1.setbrightnessunit("Jy/beam")
        sub1.setrestoringbeam(major=major, minor=minor, pa=pa)
        sub2.setbrightnessunit("Jy/beam")
        sub2.setrestoringbeam(remove=True)
        sub2.setrestoringbeam(
            major=major2, minor=minor2, pa=pa2,
            channel=0, polarization=-1
        )
        
        sub2.setrestoringbeam(
            major=major3, minor=minor3, pa=pa3,
            channel=5, polarization=-1
        )
        concat = myia.imageconcat("", [im1, im2])
        for i in range(concat.shape()[2]):
            beam = concat.restoringbeam(channel=i)
            if i < sub1.shape()[2]:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa))
            elif i == sub1.shape()[2] + 5:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major3))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor3))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa3))
            else:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major2))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor2))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa2))
                
        # both images have a single beam which is the same
        sub1.setbrightnessunit("Jy/beam")
        sub1.setrestoringbeam(major=major, minor=minor, pa=pa)
        sub2.setrestoringbeam(remove=True)
        sub2.setbrightnessunit("Jy/beam")
        sub2.setrestoringbeam(major=major, minor=minor, pa=pa)
        concat = myia.imageconcat("", [im1, im2])
        beam = concat.restoringbeam()
        self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major))
        self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor))
        self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa))
        
        # both images have single, unequal beams
        sub1.setbrightnessunit("Jy/beam")
        sub1.setrestoringbeam(major=major, minor=minor, pa=pa)
        sub2.setrestoringbeam(remove=True)
        sub2.setbrightnessunit("Jy/beam")
        sub2.setrestoringbeam(major=major2, minor=minor2, pa=pa2)
        concat = myia.imageconcat("", [im1, im2])
        for i in range(concat.shape()[2]):
            beam = concat.restoringbeam(channel=i)
            if i < sub1.shape()[2]:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa))
            else:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major2))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor2))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa2))
        
        # first image has per plane beams, second has single beam
        sub2.setbrightnessunit("Jy/beam")
        sub2.setrestoringbeam(remove=True)
        sub2.setrestoringbeam(major=major, minor=minor, pa=pa)
        sub1.setbrightnessunit("Jy/beam")
        sub1.setrestoringbeam(remove=True)
        sub1.setrestoringbeam(
            major=major2, minor=minor2, pa=pa2,
            channel=0, polarization=-1
        )
        sub1.setrestoringbeam(
            major=major3, minor=minor3, pa=pa3,
            channel=5, polarization=-1
        )
        concat = myia.imageconcat("", [im1, im2])
        for i in range(concat.shape()[2]):
            beam = concat.restoringbeam(channel=i)
            if i == 5:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major3))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor3))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa3))
            elif i < sub1.shape()[2]:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major2))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor2))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa2))
            else:
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor))
                self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa))
                
        # both images have a single beam which is the same
        sub1.setbrightnessunit("Jy/beam")
        sub1.setrestoringbeam(remove=True)
        sub1.setrestoringbeam(major=major, minor=minor, pa=pa)
        sub2.setrestoringbeam(remove=True)
        sub2.setbrightnessunit("Jy/beam")
        sub2.setrestoringbeam(major=major, minor=minor, pa=pa)
        concat = myia.imageconcat("", [im1, im2])
        beam = concat.restoringbeam()
        self.assertTrue(self.qa.eq(self.qa.quantity(beam["major"]), major))
        self.assertTrue(self.qa.eq(self.qa.quantity(beam["minor"]), minor))
        self.assertTrue(self.qa.eq(self.qa.quantity(beam["positionangle"]), pa))
        
    def test_basic(self):
        """Test basic functionality"""
        myia = self._myia
        myia.fromshape("", [1, 1, 5])
        names = []
        for i in range(5):
            name = "chan_" + str(i)
            names.append(name) 
            subi = myia.subimage(name, region=self.rg.box([0, 0, i], [0, 0,i]))
            got = subi.toworld([0 ,0, 0])['numeric'][2]
            expec = myia.toworld([0 ,0, i])['numeric'][2]
            self.assertTrue(got == expec)
            subi.done()
        concat = myia.imageconcat(infiles=[names[0], names[1], names[2]])
        for i in range(3):
            got = concat.toworld([0 ,0, i])['numeric'][2]
            expec = myia.toworld([0 ,0, i])['numeric'][2]
            self.assertTrue(got == expec)
            
        self.assertRaises(
            Exception, myia.imageconcat, outfile="blah.im", infiles=[names[0], names[1], names[3]]
        )
        concat = myia.imageconcat(
            infiles=[names[0], names[1], names[3]], relax=True
        )
        for i in range(3):
            k = i
            if i == 2: k = 3
            got = concat.toworld([0 ,0, i])['numeric'][2]
            expec = myia.toworld([0 ,0, k])['numeric'][2]
            self.assertTrue(got == expec)
            
        concat = myia.imageconcat(
            infiles=[names[0], names[1], names[3], names[4]],
            relax=True
        )
        for i in range(4):
            k = i
            if i >= 2: k = i+1
            got = concat.toworld([0 ,0, i])['numeric'][2]
            expec = myia.toworld([0 ,0, k])['numeric'][2]
            self.assertTrue(got == expec)
        concat.done()
        myia.done()
        
    def test_reorder(self):
        """Test reorder param functionality"""
        myia = self._myia
        myia.fromshape("", [1, 1, 5])
        names = []
        for i in range(5):
            name = "reorder_" + str(i)
            names.append(name) 
            subi = myia.subimage(name, region=self.rg.box([0, 0, i], [0, 0,i]))
            got = subi.toworld([0 ,0, 0])['numeric'][2]
            expec = myia.toworld([0 ,0, i])['numeric'][2]
            self.assertTrue(got == expec)
            bb = subi.getchunk()
            bb[:] = i
            subi.putchunk(bb)
            subi.done()
            
        concat = myia.imageconcat(
            infiles=[names[1], names[3], names[2]], reorder=True
        )
        for i in range(3):
            myia.open(names[i + 1])
            got = concat.toworld([0 ,0, i])['numeric'][2]
            expec = myia.toworld([0 ,0, 0])['numeric'][2]
            self.assertTrue(got == expec)
            self.assertTrue(
                concat.getchunk()[0, 0, i] == myia.getchunk()[0, 0, 0]
            )
        myia.done()

    def test_history(self):
        """Test history writing"""
        im1 = "hist1.im"
        im2 = "hist2.im"
        myia = self._myia
        myia.fromshape(im1, [10,10,5])
        myia.fromshape(im2, [10,10,5])
        myia.done()
        zz = myia.imageconcat("", [im1, im2], relax=True)
        msgs = zz.history()
        self.assertTrue("ia.imageconcat" in msgs[-1])
        self.assertTrue("ia.imageconcat" in msgs[-2])

    def test_precision(self):
        """Test different image precisions"""
        myia = self._myia
        shape = [4, 4, 5]
        expec = {}
        expec['f'] = 'float'
        expec['c'] = 'complex'
        expec['d'] = 'double'
        expec['cd'] = 'dcomplex'
        for mytype in ['f', 'c', 'd', 'cd']:
            out0 = "c0_" + mytype + ".im"
            myia.fromshape(out0, shape, type=mytype)
            out1 = "c1_" + mytype + ".im"
            myia.fromshape(out1, shape, type=mytype)
            myia.done()
            concat = myia.imageconcat("", [out0, out1], relax=True, axis=2)
            myia.done()
            self.assertTrue(
                concat.pixeltype() == expec[mytype], "wrong type for " + mytype
            )
            concat.done()
 
    def test_mode(self):
        """Test various output formats CAS-12600"""
        # 'm' has to be the last one...
        modes = ['c', 'n', 'p', 'm']
        (names, myia) = self.make_images()
        for mode in modes:
            outname = "loop1_" + mode + ".im"
            concat = myia.imageconcat(outname, infiles=names, mode=mode)
            for i in range(3):
                got = concat.toworld([0 ,0, i])['numeric'][2]
                expec = myia.toworld([0 ,0, i])['numeric'][2]
                self.assertTrue(got == expec)
        (names, myia) = self.make_images()
        for mode in modes:
            outname = "loop2_" + mode + ".im"
            concat = myia.imageconcat(
                outname, infiles=[names[0], names[1], names[3]], relax=True, mode=mode
            )
            for i in range(3):
                k = i
                if i == 2: k = 3
                got = concat.toworld([0 ,0, i])['numeric'][2]
                expec = myia.toworld([0 ,0, k])['numeric'][2]
                self.assertTrue(got == expec)
        (names, myia) = self.make_images()
        for mode in modes:
            outname = "loop3_" + mode + ".im"
            concat = myia.imageconcat(
                outname, infiles=[names[0], names[1], names[3],
                names[4]], relax=True, mode=mode
            )
            for i in range(4):
                k = i
                if i >= 2: k = i+1
                got = concat.toworld([0 ,0, i])['numeric'][2]
                expec = myia.toworld([0 ,0, k])['numeric'][2]
                self.assertTrue(got == expec)
        concat.done()
        myia.done() 

def suite():
    return [ia_imageconcat_test]

if __name__ == '__main__':
    unittest.main()

##########################################################################
# test_req_tool_ia_pbcor.py
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
# https://open-jira.nrao.edu/browse/CAS-12986
#
##########################################################################

from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import unittest

is_CASA6 = False
try:
    from casatools import ctsys, image, table
    _ia = image()
    _tb = table()
    datapath = ctsys.resolve('unittest/ia_pbcor/')
    is_CASA6 = True
except ImportError:
    import casac
    from tasks import *
    from taskinit import *
    from __main__ import *
    image = iatool
    _ia = iatool()
    _tb = tbtool()
    datapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/ia_pbcor/'
    is_CASA6 = False

im1 = "pbtest1_im.fits"
pb1 = "pbtest1_pb.fits"
co1_1 = "pbtest1_co1.fits"
co1_2 = "pbtest1_co2.im"

im2 = "pb2_im.fits"
pb2 = "pb2_pb.fits"  
co2 = "pb2_co.im"  

pb4 = "CAS_5096template.im"

data = [im1, pb1, co1_1, co1_2, im2, pb2, co2, pb4]

def run_pbcor(
    imagename, pbimage, outfile, overwrite, region, box, chans,
    stokes, mask, mode, cutoff
):
    myia = image()
    myia.open(imagename)
    res = myia.pbcor(
        pbimage=pbimage, outfile=outfile, overwrite=overwrite,
        region=region, box=box, chans=chans, stokes=stokes,
        mask=mask, mode=mode, cutoff=cutoff
    )
    myia.close()
    myia.done()
    return res

class ia_pbcor_test(unittest.TestCase):
    
    def setUp(self):
        for f in data:
            resolved = os.path.join(datapath, f)
            if os.path.isdir(resolved):
                shutil.copytree(resolved, f)
            if os.path.isfile(resolved):
                shutil.copy(resolved,f)
    
    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)
        # make sure directory is clean as per verification test requirement
        cwd = os.getcwd()
        for filename in os.listdir(cwd):
            file_path = os.path.join(cwd, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    # CASA 5 tests need this directory
                    if filename != 'xml':
                        shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

    def checkImage(self, gotImage, expectedName, epsilon):
        expected = image()                                
        expected.open(expectedName)
        got = image()
        if type(gotImage) == str:
            got.open(gotImage)
        else:
            got = gotImage
        self.assertTrue((got.shape() == expected.shape()).all(), 'image shape failure')
        diffData = got.getchunk()/expected.getchunk() - 1
        self.assertTrue(abs(diffData).max() <= epsilon, 'pixel value failure')
        gotCsys = got.coordsys()
        expectedCsys = expected.coordsys()
        diffPixels = gotCsys.referencepixel()['numeric'] - expectedCsys.referencepixel()['numeric']
        self.assertTrue(abs(diffPixels).max() == 0)
        denom = expectedCsys.referencevalue()['numeric']
        for i in range(len(denom)):
            if (denom[i] == 0):
                denom[i] = 1
        fracDiffRef = (
            gotCsys.referencevalue()['numeric'] - expectedCsys.referencevalue()['numeric']
        )/denom;
        self.assertTrue(abs(fracDiffRef).max() <= 1.5e-6, 'refval failure')
        gotnpts = got.statistics()['npts']
        expnpts = expected.statistics()['npts']
        self.assertTrue(gotnpts == expnpts, 'npts failure')
        got.done()
        expected.done()

    def test_exceptions(self):
        """ia.pbcor: Test various exception cases"""
        def testit(
            imagename, pbimage, outfile, overwrite, region,
            box, chans, stokes, mask, mode, cutoff, wantreturn
        ):
            self.assertRaises(
                Exception, run_pbcor, imagename=imagename,
                pbimage=pbimage, outfile=outfile, overwrite=overwrite,
                region=region, box=box, chans=chans,
                stokes=stokes, mask=mask, mode=mode,
                cutoff=cutoff
            )
                       
        # no image name given
        testit(
            imagename="", pbimage=pb1, outfile="",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )
        # bad image name given
        testit(
            imagename="totally_bogus", pbimage=pb1, outfile="jj.im",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )
        # no pbimage name given
        testit(
            imagename=im1, pbimage="", outfile="mm.im",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )
        # bad pbimage name given
        testit(
            imagename=im1, pbimage="totally_bogus2", outfile="pp.im",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )
        # unwritable outfile
        testit(
            imagename=im1, pbimage=pb1, outfile="/bogusplace/bogusimage",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )
        # bogus region
        testit(
            imagename=im1, pbimage=pb1, outfile="qq.im",
            overwrite=False, region="bogus_region", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )
        # bad mode
        testit(
            imagename=im1, pbimage=pb1, outfile="rr.im",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="zz", cutoff=-1.0,
            wantreturn=True
        )
        # incompatible image and pb
        testit(
            imagename=im1, pbimage=pb2, outfile="ss.im",
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d", cutoff=-1.0,
            wantreturn=True
        )

    def _testit(
        self, expected, imagename, pbimage, overwrite, region, box,
        chans, stokes, mask, mode, cutoff, epsilon=0.0
    ):
        myia = image()
        myia.open(pbimage)
        pbpix = myia.getchunk()
        myia.done()
        del myia
        for j in [0, 1]:
            outfile = "mypb.im"
            if j == 1:
                pbimage = pbpix
            outfile = outfile
            res = run_pbcor(
                imagename=imagename, pbimage=pbimage,
                outfile=outfile, overwrite=overwrite,
                region=region, box=box, chans=chans,
                stokes=stokes, mask=mask, mode=mode,
                cutoff=cutoff
            )
            self.assertTrue(res)
            res.done()
            self.checkImage(outfile, expected, epsilon)
            shutil.rmtree(outfile)

    def test_full_image_divide(self):
        """ia.pbcor: Test full image divide"""
        self._testit(
            expected=co1_1, imagename=im1, pbimage=pb1,
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d",
            cutoff=-1.0
        )

    def test_full_image_using_cutoff(self):
        """ia.pbcor: Test full image divide with cutoff"""
        self._testit(
            expected=co1_2, imagename=im1, pbimage=pb1,
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d",
            cutoff=0.001
        )
    
    def test_4d_image_with_2d_pb(self):
        """ia.pbcor: Test full image divide with cutoff. Primary beam is 2 D, image is 4 D"""
        self._testit(
            expected=co2, imagename=im2, pbimage=pb2,
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="d",
            cutoff=0.001
        )

    def test_multiply(self):
        """ia.pbcor: Test full image multiply with cutoff. Primary beam is 2 D, image is 4 D"""
        myia = image()
        myia.open(pb2)
        pixels = myia.getchunk()
        csys = myia.coordsys().torecord()
        myia.done()
        pixels = 1/pixels
        newpb = 'mult_pb.im'
        myia.fromarray(newpb, pixels, csys=csys)
        myia.done()
        self._testit(
            expected=co2, imagename=im2, pbimage=newpb,
            overwrite=False, region="", box="",
            chans="", stokes="", mask="", mode="m",
            cutoff=1000, epsilon=2e-7
        )

    def test_stretch(self):
        """ ia.pbcor(): Test stretch parameter"""
        yy = image()
        mymask = "maskim"
        yy.fromshape("", [113, 76, 1, 1])
        yy.addnoise()
        xx = yy.transpose(mymask, "0132")
        yy.done()
        xx.done()
        yy.open(im2)
        self.assertRaises(
            RuntimeError, yy.pbcor, pbimage=pb2,
            mask=mymask + ">0", stretch=False, outfile="garbage"
        )
        yy.open(im2)
        zz = yy.pbcor(
            pbimage=pb2, outfile="blahblah", mask=mymask + ">0", stretch=True
        )
        yy.done()
        self.assertTrue(zz)
        zz.done()
        
    def test_diff_spectral_coordinate(self):
        """Verify fix that a different spectral coordinates in target and template don't matter, CAS-5096"""
        imagename = os.path.join(datapath, "CAS_5096target.im")
        template = pb4
        outfile = "mypb.im"
        myia = image()
        myia.open(imagename)
        yy = myia.pbcor(
            pbimage=template,
            outfile=outfile, mask='"' + template + '">0.21'
        )
        myia.done()
        yy.done()
        self.assertTrue(os.path.exists(outfile))

    def test_history(self):
        """Test history records are written"""
        myia = image()
        imagename = "zz.im"
        myia.fromshape(imagename, [20,20])
        gg = myia.getchunk()
        gg[:] = 1
        myia.done()
        outfile = "pb_out.im"
        myia.open(imagename)
        yy = myia.pbcor(pbimage=gg, outfile=outfile)
        myia.done()
        msgs = yy.history()
        yy.done()
        teststr = "ia.pbcor"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")    
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")

    def test_empty_region_and_stokes(self):
        """Test specifying stokes with empty region works (CAS-11708)"""
        myia = image()
        imagename = "t_in.im"
        pbimage = "t_pb_in.im"
        shape = [20, 20, 4, 20]
        for im in [imagename, pbimage]:
            myia.fromshape(im, shape)
            myia.addnoise()
            myia.done()
        outfile = "t_out.im"
        myia.open(imagename)
        yy = myia.pbcor(pbimage, outfile=outfile, region="", stokes="I")
        myia.done()
        self.assertTrue((yy.shape() == [20,20,1,20]).all(), "Incorrect shape")
        yy.done()

def suite():
    return [ia_pbcor_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

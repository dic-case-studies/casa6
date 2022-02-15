##########################################################################
# test_tool_regionmanager.py
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
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.regionmanager.html
#
#
##########################################################################

import os
import shutil
import numpy
import unittest

from casatools import ctsys, quanta, imagemetadata, table
from casatools import regionmanager as rgtool
from casatools import image as iatool

_qa = quanta()
_imd = imagemetadata()

datapath = ctsys.resolve('unittest/rgtool/')

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

# Tests for regionmanager.frombcs
class rg_frombcs_test(unittest.TestCase):
    image = "imregion.fits"
    image_nospec = "imregion_nospec.fits"
    image_dironly = "imregion_dironly.fits"

    box1 = 1.24795026
    box2 = 0.782552901
    box3 = 1.24794616
    box4 = 0.782555814
    box5 = 1.24794206
    box6 = 0.782558727
    box7 = 1.24793797
    box8 = 0.782561641
    box9 = 1.2479338718038551
    box10 = 0.78256455381109313
    box11 = 1.2479297756405987
    box12 = 0.78256746696533663
    chan0 = 4.73510000e+09
    chan4 = 6.33510000e+09
    chan15 = 1.07351000e+10
    chan19 = 1.23351000e+10

    def run_frombcs(self, imagename, box, chans, stokes, stokes_control, region=""):
        myia = iatool()
        myia.open(imagename)
        mycsys = myia.coordsys()
        myrg = rgtool()
        res = myrg.frombcs(
            mycsys.torecord(), myia.shape(), box, chans,
            stokes, stokes_control, region
        )
        myia.close()
        del myia
        del myrg
        return res

    def recToList(self, rec):
        mykeys = sorted(rec.keys())
        return [rec[k]['value'] for k in mykeys]

    def setUp(self):
        for im in [self.image, self.image_nospec, self.image_dironly]:
            shutil.copy(datapath + im, im)

    def tearDown(self):
        for im in [self.image, self.image_nospec, self.image_dironly]:
            os.remove(im)
        shutil.rmtree('xxyy.im', ignore_errors=True)

    def compLists(self, got, exp):
        epsilon = 1e-8
        print("got " + str(got))
        print("exp " + str(exp))
        for i in range(len(got)):
            fracDiff = abs((got[i] - exp[i]) / exp[i]);
            self.assertTrue(fracDiff < epsilon)

    def test_full_image(self):
        """Test default gives region of entire image"""

        stokes = ""
        chans = ""
        stokes_control = "a"
        box = ""
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 4.0]
        self.compLists(gottrc, exptrc)

    def test_single_stokes(self):
        """Test setting a single stokes"""

        stokes = "Q"
        chans = ""
        stokes_control = "a"
        box = ""
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 2.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 2.0]

    def test_continguous_stokes(self):
        """Test setting a contiguous stokes"""

        box = ""
        chans = ""
        stokes = "QU"
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 2.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 3.0]

    def test_single_channel(self):
        """Test setting a single channel"""
        box = ""
        chans = "4"
        stokes = ""
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [1.24795230, 0.782549990, self.chan4, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan4, 4.0]

    def test_contiguous_channels(self):
        """Test setting multiple continuous channels"""

        box = ""
        chans = "0~4"
        stokes = ""
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan4, 4.0]

    def test_single_box(self):
        """Test setting single box"""

        box = "1,2,3,4"
        chans = ""
        stokes = ""
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [self.box1, self.box2, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [self.box3, self.box4, self.chan19, 4.0]

    def test_region_record(self):
        """Test setting region record"""

        box = ""
        chans = ""
        stokes = ""
        stokes_control = "a"
        blahia = iatool()
        blahia.open(self.image)
        mycsys = blahia.coordsys()
        blahia.done()
        myrg = rgtool()
        mybox = myrg.wbox(
            ["1pix", "2pix", "0pix", "0pix"],
            ["3pix", "4pix", "19pix", "3pix"],
            csys=mycsys.torecord()
        )
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control, mybox
        )
        gotblc = self.recToList(myreg["blc"])
        for i in range(len(gotblc)):
            gotblc[i] = gotblc[i] - 1
        gotblc = mycsys.toworld(gotblc)['numeric']
        expblc = [self.box1, self.box2, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        for i in range(len(gottrc)):
            gottrc[i] = gottrc[i] - 1
        gottrc = mycsys.toworld(gottrc)['numeric']
        exptrc = [self.box3, self.box4, self.chan19, 4.0]
        self.compLists(gottrc, exptrc);
        del myrg

    def test_first_stokes(self):
        """Test setting first stokes"""

        box = ""
        chans = ""
        stokes = ""
        stokes_control = "f"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )
        gotblc = self.recToList(myreg["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 1.0]

    def test_multiple_boxes(self):
        """Test setting multiple boxes"""

        stokes = ""
        chans = ""
        stokesControl = "a"
        box = "1,2,3,4,5,6,7,8,9,10,11,12"
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )

        gotblc = self.recToList(myreg["regions"]["*2"]["blc"])
        expblc = [1.24793387, 0.782564554, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*2"]["trc"])
        exptrc = [1.24792978, 0.782567467, self.chan19, 4.0]
        self.compLists(gottrc, exptrc)

        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["blc"])
        expblc = [self.box1, self.box2, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["trc"])
        exptrc = [self.box3, self.box4, self.chan19, 4.0]
        self.compLists(gottrc, exptrc);

        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, self.chan19, 4.0]
        self.compLists(gottrc, exptrc);

    def test_set_multiple_stokes_ranges(self):
        """Test setting multiple stokes ranges"""

        stokes = "IUV"
        chans = ""
        stokesControl = "a"
        box = ""
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )

        gotblc = self.recToList(myreg["regions"]["*1"]["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*1"]["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 1.0]
        self.compLists(gottrc, exptrc)

        gotblc = self.recToList(myreg["regions"]["*2"]["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 3.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*2"]["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 4.0]
        self.compLists(gottrc, exptrc)

    def test_multiple_channel_ranges(self):
        """Test multiple channel ranges"""

        stokes = ""
        chans = "<5,>=15"
        stokesControl = "a"
        box = ""
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )

        gotblc = self.recToList(myreg["regions"]["*1"]["blc"])
        expblc = [1.24795230, 0.782549990, self.chan0, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*1"]["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan4, 4.0]
        self.compLists(gottrc, exptrc)

        gotblc = self.recToList(myreg["regions"]["*2"]["blc"])
        expblc = [1.24795230, 0.782549990, self.chan15, 1.0]
        self.compLists(gotblc, expblc);
        gottrc = self.recToList(myreg["regions"]["*2"]["trc"])
        exptrc = [1.24791339, 0.782577665, self.chan19, 4.0]
        self.compLists(gottrc, exptrc)

    def test_multiple_boxes_channel_ranges_stokes_ranges(self):
        """Test multiple channel ranges, multiple stokes ranges, and multiple boxes"""

        stokes = "IQV"
        chans = "<5,>=15"
        stokesControl = "a"
        box = "1,2,3,4,5,6,7,8"
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image, box, chans, stokes, stokes_control
        )

        # box="5,6,7,8", chans="15~19", stokes="V"
        gotblc = self.recToList(myreg["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, self.chan15, 4.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, self.chan19, 4.0]
        self.compLists(gottrc, exptrc)

        # box="5,6,7,8", chans="0~4", stokes="V"
        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, self.chan0, 4.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, self.chan4, 4.0]
        self.compLists(gottrc, exptrc)

        # box="5,6,7,8", chans="15-19", stokes="IQ"
        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, self.chan15, 1.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, self.chan19, 2.0]
        self.compLists(gottrc, exptrc)

        # box="5,6,7,8", chans="0~4", stokes="IQ"
        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, self.chan0, 1.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, self.chan4, 2.0]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", chans="15-19", stokes="V"
        gotblc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*2"]["blc"]
        )
        expblc = [self.box1, self.box2, self.chan15, 4.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*2"]["trc"]
        )
        exptrc = [self.box3, self.box4, self.chan19, 4.0]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", chans="0-4", stokes="V"
        gotblc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*1"]["regions"]["*2"]["blc"]
        )
        expblc = [self.box1, self.box2, self.chan0, 4.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*1"]["regions"]["*2"]["trc"]
        )
        exptrc = [self.box3, self.box4, self.chan4, 4.0]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", chans="15-19", stokes="IQ"
        gotblc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["blc"]
        )
        expblc = [self.box1, self.box2, self.chan15, 1.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["trc"]
        )
        exptrc = [self.box3, self.box4, self.chan19, 2.0]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", chans="0~4", stokes="IQ"
        gotblc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["blc"]
        )
        expblc = [self.box1, self.box2, self.chan0, 1.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(
            myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]
            ["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["trc"]
        )
        exptrc = [self.box3, self.box4, self.chan4, 2.0]
        self.compLists(gottrc, exptrc)

    def test_multiple_boxes_multiple_stokes_no_spectral_axis(self):
        """Test multiple stokes ranges, and multiple boxes on image with no spectral axis"""

        stokes = "IQV"
        chans = ""
        stokesControl = "a"
        box = "1,2,3,4,5,6,7,8"
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image_nospec, box, chans, stokes, stokes_control
        )

        # box="5,6,7,8", stokes="V"
        gotblc = self.recToList(myreg["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, 4.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, 4.0]
        self.compLists(gottrc, exptrc)

        # box="5,6,7,8", stokes="IQ"
        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6, 1.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8, 2.0]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", stokes="V"
        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["blc"])
        expblc = [self.box1, self.box2, 4.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*2"]["trc"])
        exptrc = [self.box3, self.box4, 4.0]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", stokes="IQ"
        gotblc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["blc"])
        expblc = [self.box1, self.box2, 1.0]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["regions"]["*1"]["regions"]["*1"]["trc"])
        exptrc = [self.box3, self.box4, 2.0]
        self.compLists(gottrc, exptrc)

    def test_multiple_boxes_image_with_direction_coordinate_only(self):
        """Test multiple boxes on image with direction coordinate only"""

        stokes = ""
        chans = ""
        stokesControl = "a"
        box = "1,2,3,4,5,6,7,8"
        stokes_control = "a"
        myreg = self.run_frombcs(
            self.image_dironly, box, chans, stokes, stokes_control
        )

        # box="5,6,7,8", stokes="V"
        gotblc = self.recToList(myreg["regions"]["*2"]["blc"])
        expblc = [self.box5, self.box6]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*2"]["trc"])
        exptrc = [self.box7, self.box8]
        self.compLists(gottrc, exptrc)

        # box="1,2,3,4", stokes="IQ"
        gotblc = self.recToList(myreg["regions"]["*1"]["blc"])
        expblc = [self.box1, self.box2]
        self.compLists(gotblc, expblc)
        gottrc = self.recToList(myreg["regions"]["*1"]["trc"])
        exptrc = [self.box3, self.box4]
        self.compLists(gottrc, exptrc)

    def test_region_text_string(self):
        """Test setting a region text string"""

        region = "box[[1pix,2pix],[3pix,4pix]]\nbox[[5pix,6pix],[7pix,8pix]]\nbox[[9pix,10pix],[11pix,12pix]]"
        myreg = self.run_frombcs(
            self.image, "", "", "", "a", region
        )
        myia = iatool()
        myia.open(self.image)
        subi = myia.subimage("xxyy.im", region=myreg)
        myia.done()
        self.assertTrue((subi.shape()[0:2] == [11, 11]).all())
        got = subi.toworld([0, 0])['numeric'][0:2]
        expec = [self.box1, self.box2]
        self.compLists(got, expec)
        expec = [self.box11, self.box12]
        got = subi.toworld([10, 10])['numeric'][0:2]
        self.compLists(got, expec)
        gotmask = subi.getchunk(getmask=True)[:, :, 0, 0]
        subi.done()
        expmask = gotmask.copy()
        expmask[:] = False
        expmask[0:3, 0:3] = True
        expmask[4:7, 4:7] = True
        expmask[8:11, 8:11] = True
        self.assertTrue((gotmask == expmask).all())

    def test_complement_is_last_line(self):
        """Test CAS-12978 fix, difference line is last in file works"""
        myia = iatool()
        myia.fromshape("", [500, 500])
        csys = myia.coordsys().torecord()
        shape = myia.shape()
        myia.done()
        region = "ellipse [[0:0:0, 0.0.0], [10arcmin, 5arcmin], 45deg]\n- ellipse [[0:0:0, 0.0.0], [7arcmin, 2arcmin], 45deg]"
        myrg = rgtool()
        # not segfaulting verifies the test
        reg1 = myrg.frombcs(csys=csys, shape=shape, region=region)
        # also check space between "-" and shape is no longer necessary
        region = "ellipse [[0:0:0, 0.0.0], [10arcmin, 5arcmin], 45deg]\n-ellipse [[0:0:0, 0.0.0], [7arcmin, 2arcmin], 45deg]"
        reg2 = myrg.frombcs(csys=csys, shape=shape, region=region)
        # member by member comparison cannot be easily done because some values are numpy arrays
        self.assertTrue(str(reg1) == str(reg2))


# Tests for regionmanager.fromtextfile
class rg_fromtextfile_test(unittest.TestCase):

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
    cas_12980i = 'Cir_X-1_sci.spw37.mfs.I.manual.fits'
    cas_12980t = 'mynewregion.crtf'
    cas_12980c = 'cas_12980c.im'
    _fixtures = [
        image, text1, res1, cas_3258t, cas_3258r, cas_3259t, cas_3259r,
        cas_3260t, cas_3260r]

    _created = [icrs_image, cas_12980c]

    def setUp(self):
        for im in self._fixtures:
            shutil.copy(datapath + im, im)
        shutil.copy(datapath + self.cas_12980i, self.cas_12980i)
        shutil.copy(datapath + self.cas_12980t, self.cas_12980t)
        self.ia = iatool()
        self.rg = rgtool()

    def tearDown(self):
        for im in self._fixtures + self._created + [self.cas_12980i, self.cas_12980t]:
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
        self.assertRaises(Exception, self.rg.fromtextfile, "blah", {}, [1, 1])
        # coordsys not set
        self.assertRaises(Exception, self.rg.fromtextfile, self.text1, {}, [1, 1])

    def test_read(self):
        """Read test"""
        self.ia.open(self.image)
        self._testit(self.text1, self.res1)

    def test_CAS_3258(self):
        """Verify fix to CAS-3258"""
        self.ia.fromshape("", [250, 250])
        self._testit(self.cas_3258t, self.cas_3258r)

    def test_CAS_3259(self):
        """Verify fix to CAS-3259"""
        self.ia.fromshape("", [250, 250])
        self._testit(self.cas_3259t, self.cas_3259r)

    def test_CAS_3260(self):
        """Verify fix to CAS-3260"""
        self.ia.fromshape("", [250, 250])
        self._testit(self.cas_3260t, self.cas_3260r)

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
        self.ia.fromshape("", [200, 200])
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

        # center box
        self.ia.fromshape("", [200, 200])
        csys = self.ia.coordsys()
        yval = "-3611.1455480499999arcmin"
        xwidth = _qa.tos(_qa.mul(_qa.quantity("104.48212727000009arcmin"), _qa.cos(yval)))
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
        self.ia.fromshape("", [200, 200])
        reg = self.rg.fromtext(
            "rotbox [ [ 60pix , 50pix ] , [ 30pix , 30pix ] , 30deg ]",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        self.assertTrue(self.ia.statistics(region=reg)['npts'] == 901)
        csys = self.ia.coordsys()
        csys.setreferencevalue([800, 70 * 60])
        self.ia.setcoordsys(csys.torecord())
        reg = self.rg.fromtext(
            "rotbox [ [ 60pix , 50pix ] , [ 30pix , 30pix ] , 30deg ]",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        self.assertTrue(self.ia.statistics(region=reg)['npts'] == 901)

    def test_ellipse(self):
        """Test ellipse for image in GALACTIC and file in J2000"""
        self.ia.open(datapath + "gal.im")
        reg = self.rg.fromtextfile(
            datapath + "testEllipse90deg.crtf",
            csys=self.ia.coordsys().torecord(),
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
        self.assertTrue(self.ia.statistics()['npts'][0] == 331 * 331)
        self.assertTrue(self.ia.statistics(region=datapath + "1000circles.txt")['npts'][0] == 13679)
        self.ia.done()

    def test_CAS_8072(self):
        """Verify rest frequency precision issue has been fixed"""
        self.ia.fromshape("", [20, 20, 200])
        self.ia.addnoise()
        reg = reg = self.rg.fromtext(
            "box [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s], restfreq=1.42040575e+09Hz",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        reg1 = self.rg.fromtext(
            "box [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        # no comma delimiters should throw exception
        self.assertRaises(
            Exception, self.rg.fromtext,
            "global coord=B1950 frame=LSRK veltype=RADIO restfreq=1.42040575e+09Hz"
            + "\nbox [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        reg3 = self.rg.fromtext(
            "global coord=J2000, frame=LSRK, veltype=RADIO, restfreq=1.42040575e+09Hz"
            + "\nbox [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
        )
        # different global rest freq
        reg4 = self.rg.fromtext(
            "global coord=J2000, frame=LSRK, veltype=RADIO, restfreq=1.42050575e+09Hz"
            + "\nbox [[0pix,0pix], [19pix,19pix]], range=[1140km/s, 1142km/s]",
            csys=self.ia.coordsys().torecord(), shape=self.ia.shape()
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
        self.ia.fromshape(self.icrs_image, shape=[100, 100])
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
            _imd.open(self.icrs_image), 'Unable to open imagemetadata object'
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
        self.assertTrue(self.ia.open(self.icrs_image), 'Failed to open image')
        reg = rg.fromtext(icrs, csys=self.ia.coordsys().torecord(), shape=self.ia.shape())
        rg.done()
        stats = self.ia.statistics(region=reg)
        self.ia.done()
        self.assertEqual(stats['max'][0], 1, 'Incorrect value for max')

    def test_crtf_has_multiple_diff_and_ends_with_diff_bug_fix(self):
        """
        CAS-12980, verify fix that allows supplied CRTF file to work correctly
        1. copy necessary data
        2. run ia.subimage() on image using region file
        3. confirm that the correct number of pixels are not masked
        """
        self.ia.open(self.cas_12980i)
        xx = self.ia.subimage(region=self.cas_12980t)
        self.assertEqual(xx.statistics()['npts'], 6612, 'Wrong number of pixels masked')
        xx.done()

# Tests for regionmanager.selectedchannels
class rg_selectedchannels_test(unittest.TestCase):

    def setUp(self):
        self.rg = rgtool()

    def tearDown(self):
        self.rg.done()
        tb = table()
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done()

    def test_no_spectral_axis(self):
        """ Test no spectral axis throws exception"""
        myia = iatool()
        myia.fromshape("", [4, 4, 4])
        self.rg.setcoordinates(myia.coordsys().torecord())
        shape = myia.shape()
        self.assertRaises(
            Exception, self.rg.selectedchannels,
            "range=[1415MHz, 1415.1MHz]", shape
        )
        myia.done()

    def test_no_overlap(self):
        """ Test selection outside of image"""
        myia = iatool()
        myia.fromshape("", [4, 4, 20])
        self.rg.setcoordinates(myia.coordsys().torecord())
        shape = myia.shape()
        self.assertRaises(
            Exception, self.rg.selectedchannels,
            "range=[1315MHz, 1315.1MHz]", shape
        )
        myia.done()

    def test_range(self):
        """ Test range"""
        myia = iatool()
        myia.fromshape("", [4, 4, 20])
        self.rg.setcoordinates(myia.coordsys().torecord())
        shape = myia.shape()
        myia.done()

        chans = self.rg.selectedchannels("range=[1415MHz,1415.002MHz]", shape)
        self.assertTrue((chans == [10, 11, 12]).all())
        chans = self.rg.selectedchannels("range=[1415.002MHz, 1415MHz]", shape)
        self.assertTrue((chans == [10, 11, 12]).all())

        chans = self.rg.selectedchannels("range=[1414MHz,1415.002MHz]", shape)
        self.assertTrue((chans == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]).all())

        chans = self.rg.selectedchannels("range=[1415MHz,1415.2MHz]", shape)
        self.assertTrue((chans == [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]).all())


if __name__ == '__main__':
    unittest.main()

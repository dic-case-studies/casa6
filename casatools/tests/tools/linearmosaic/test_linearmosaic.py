##########################################################################
# test_linearmosaic.py
#
# Copyright (C) 2021
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
# Initial test script for the linearmosaic tool. CAS-13404
#
# Based on the requirements as inferred from the linearmosaic tool API
#    https://casadocs.readthedocs.io/en/stable/api/tt/casatools.linearmosaic.html
#
# Comparisons are made against a reference image created from a test dataset
# with 2 pointings and 1 source between the pointings. The source intensity is
# 1 Jy with a spectral index of -0.5 and a reference frequency at channel 1
# (the middle channel in the dataset). For any imaging scheme, the channel 1
# image should match the channel 1 primary beam for it to give the correct 1 Jy
# pb corrected source value.
#
# Input images were made using tclean with the following MS from casatestdata:
#  measurementset/evla/refim_oneshiftpoint.mosaic.ms
# 
# An image for comparison was made using both pointings with tclean and
#  imagename='jointmos'
#  niter=20
#  specmode='cube'
#  spw='*'
#  imsize=1024
#  phasecenter='J2000 19h59m28.5 +40d40m01.5'
#  cell='8.0arcsec'
#  gridder='mosaic'
#  field='0,1'
#  conjbeams=False
#  wbawp=True
#  pblimit=0.1
#  deconvolve='hogbom'
#  pbcor=True
#
# The two individual pointing images were made with tclean 
#  and the same values except field='0' or field='1'.
#
# The jointmos image is not (yet?) used here so it has not been saved
# in the casatestdata repository. The individual pointing images and
# associated pb and pbcor images are found there and used here.
#
# Not tested here: 
#   independent x and y axis sizes (ny parameter)
#   independent x and y axis cell size (celly parameter)
#
##########################################################################
 
 
####    Imports     ####
import os
import shutil
import unittest

from casatools import image
from casatools import linearmosaic
 
from casatools import ctsys

####    Tests     ####
class test_linearmosaic(unittest.TestCase):

    def readpix(self, imname=''):
        '''get the center pixel value from the named image'''
        self.ia.open(imname)
        result = self.ia.pixelvalue(self.center_pix)['value']['value']
        self.ia.close()
        return result
        
    ### Set Up
    @classmethod
    def setUpClass(cls):
        '''The test images are used in all tests, make one set of copies here'''

        datapath = ctsys.resolve('unittest/linearmosaic')

        # local image tool to be used in readpix
        cls.ia = image()

        # this is the center pixel used for comparison - channel 1 is where the source = 1 Jy
        cls.center_pix = [512,617,0,1]

        # phasecenter of joint image / combined images
        cls.phasecenter ='J2000 19h59m28.5 +40d40m01.5'

        # the individual pointing image and pb to be used during combination
        # use at their location in the datapath, no need for a copy here
        # pointing 0
        cls.pnt0_im = os.path.join(datapath,'pnt0_test_lmtool.image')
        cls.pnt0_pb = os.path.join(datapath,'pnt0_test_lmtool.pb')
        cls.pnt0_im_pbcor = os.path.join(datapath,'pnt0_test_lmtool.image.pbcor')

        # pointing 1
        cls.pnt1_im = os.path.join(datapath,'pnt1_test_lmtool.image')
        cls.pnt1_pb = os.path.join(datapath,'pnt1_test_lmtool.pb')
        cls.pnt1_im_pbcor = os.path.join(datapath,'pnt1_test_lmtool.image.pbcor')

    def setUp(self):
        '''Method called to prepare the test fixture.'''
        self.output_list = []
 
    ### Teardown
    def tearDown(self):
        '''cleanup anything in self.output_list'''
        for f in self.output_list:
            if os.path.exists(f):
                shutil.rmtree(f)
 
    @classmethod
    def tearDownClass(cls):
        cls.ia.done()

    ### Test Cases
    def test_flat_noise(self):
        '''test_flat_noise (flat noise type with flat noise input)'''

        self.output_list = ['linmos.fn.image','linmos.fn.image.weight','linmos.fn.sault.image']
 
        lm = linearmosaic()
        lm.defineoutputimage(nx=1024, cellx='8arcsec', imagecenter=self.phasecenter, outputimage='linmos.fn.image')
        lm.setlinmostype(linmostype='pbweight') ## flat noise
        lm.makemosaic(images=[self.pnt0_im, self.pnt1_im], weightimages=[self.pnt0_pb, self.pnt1_pb], imageweighttype=1, weighttype=1)
        lm.saultweightimage(outputimage='linmos.fn.sault.image',fracpeak=0.3)
        del lm
        
        # center pixel values from each image

        fn_image_val = self.readpix('linmos.fn.image')
        fn_weight_val = self.readpix('linmos.fn.image.weight')
        sw_image_val = self.readpix('linmos.fn.sault.image')

        # fn_image_val / fn_image_weight is not within a few % of 1, as it probably should be
        # For now, check against expected values in channel 1 at center pixel
        self.assertAlmostEqual(fn_image_val, 0.9902307)
        self.assertAlmostEqual(fn_weight_val, 1.3866229)

        # the sault image value should be the same as the image value
        # possibly AlmostEqual should be used here eventually, but at the
        # moment they are identical.
        self.assertEqual(fn_image_val, sw_image_val)

    def test_flat_noise_flat_sky_input(self):
        '''test_flat_noise_flat_sky_input (flat noise type with flat sky input)'''

        self.output_list = ['linmos.fn2.image','linmos.fn2.image.weight','linmos.fn2.sault.image']
 
        lm = linearmosaic()
        lm.defineoutputimage(nx=1024, cellx='8arcsec', imagecenter=self.phasecenter, outputimage='linmos.fn2.image')
        lm.setlinmostype(linmostype='pbweight') ## flat noise
        # use pbcor images as input
        lm.makemosaic(images=[self.pnt0_im_pbcor, self.pnt1_im_pbcor], weightimages=[self.pnt0_pb, self.pnt1_pb], imageweighttype=0, weighttype=1)
        lm.saultweightimage(outputimage='linmos.fn2.sault.image',fracpeak=0.3)
        del lm
        
        fn_image_val = self.readpix('linmos.fn2.image')
        fn_weight_val = self.readpix('linmos.fn2.image.weight')
        sw_image_val = self.readpix('linmos.fn2.sault.image')

        # fn_image_val / fn_image_weight is not within a few % of 1, as it probably should be
        # For now, check against expected values in channel 1 at center pixel
        self.assertAlmostEqual(fn_image_val, 0.9902307)
        self.assertAlmostEqual(fn_weight_val, 1.3866229)

        # the sault image value should be the same as the image value
        # possibly AlmostEqual should be used here eventually, but at the
        # moment they are identical.
        self.assertEqual(fn_image_val, sw_image_val)

    def test_flat_sky(self):
        '''test_flat_sky (flat sky type with flat noise input'''

        self.output_list = ['linmos.fs.image','linmos.fs.image.weight','linmos.fs.sault.image']
 
        lm = linearmosaic()
        lm.defineoutputimage(nx=1024, cellx='8arcsec', imagecenter=self.phasecenter, outputimage='linmos.fs.image')
        lm.setlinmostype(linmostype='optimal') ## flat sky
        lm.makemosaic(images=[self.pnt0_im, self.pnt1_im], weightimages=[self.pnt0_pb, self.pnt1_pb], imageweighttype=1, weighttype=1)
        lm.saultweightimage(outputimage='linmos.fs.sault.image',fracpeak=0.3)
        del lm
        
        fn_image_val = self.readpix('linmos.fs.image')
        fn_weight_val = self.readpix('linmos.fs.image.weight')
        sw_image_val = self.readpix('linmos.fs.sault.image')

        # check against expected values in channel 1 at center pixel
        self.assertAlmostEqual(fn_image_val, 0.9899020)
        self.assertAlmostEqual(fn_weight_val, 1.0370520)
        # check that the ratio of these values is as expected (within 5% of 1.0)
        self.assertTrue(abs(1.0-fn_image_val/fn_weight_val) < 0.05)

        # the sault image value should be the same as the image value
        # possibly AlmostEqual should be used here eventually, but at the
        # moment they are identical.
        self.assertEqual(fn_image_val, sw_image_val)
 
    def test_flat_sky_flat_sky_input(self):
        '''test_flat_sky_flat_sky_input (flat sky type with flat sky inputs)'''

        self.output_list = ['linmos.fs2.image','linmos.fs2.image.weight','linmos.fs2.sault.image']
 
        lm = linearmosaic()
        lm.defineoutputimage(nx=1024, cellx='8arcsec', imagecenter=self.phasecenter, outputimage='linmos.fs2.image')
        lm.setlinmostype(linmostype='optimal') ## flat sky
        # use pbcor images as input
        lm.makemosaic(images=[self.pnt0_im_pbcor, self.pnt1_im_pbcor], weightimages=[self.pnt0_pb, self.pnt1_pb], imageweighttype=0, weighttype=1)
        lm.saultweightimage(outputimage='linmos.fs2.sault.image',fracpeak=0.3)
        del lm
        
        fn_image_val = self.readpix('linmos.fs2.image')
        fn_weight_val = self.readpix('linmos.fs2.image.weight')
        sw_image_val = self.readpix('linmos.fs2.sault.image')

        # check against expected values in channel 1 at center pixel
        self.assertAlmostEqual(fn_image_val, 0.9899020)
        self.assertAlmostEqual(fn_weight_val, 1.0370520)
        # check that the ratio of these values is as expected (within 5% of 1.0)
        self.assertTrue(abs(1.0-fn_image_val/fn_weight_val) < 0.05)

        # the sault image value should be the same as the image value
        # possibly AlmostEqual should be used here eventually, but at the
        # moment they are identical.
        self.assertEqual(fn_image_val, sw_image_val)
 
    def test_flat_sky_existing_image(self):
        '''test_flat_sky_existing image (flat sky type onto an existing image)'''

        self.output_list = ['linmos.step.image','linmos.step.pb','linmos.step.sault.image']

        # copy pbcorr image and pb for pnt0 to be used as starting output image and weight
        shutil.copytree(self.pnt0_im_pbcor, 'linmos.step.image')
        shutil.copytree(self.pnt0_pb, 'linmos.step.pb')
 
        lm = linearmosaic()
        lm.setoutputimage(outputimage='linmos.step.image', outputweight='linmos.step.pb',imageweighttype=1,weighttype=1)
        lm.setlinmostype(linmostype='optimal') ## flat sky
        # second image is NOT the pbcor image, pb is used as weight
        lm.makemosaic(images=[self.pnt1_im], weightimages=[self.pnt1_pb], imageweighttype=0, weighttype=1)
        lm.saultweightimage(outputimage='linmos.step.sault.image',fracpeak=0.3)
        del lm
        
        fn_image_val = self.readpix('linmos.step.image')
        fn_weight_val = self.readpix('linmos.step.pb')
        sw_image_val = self.readpix('linmos.step.sault.image')

        # check against expected values in channel 1 at center pixel
        self.assertAlmostEqual(fn_image_val, 1.1447070)
        self.assertAlmostEqual(fn_weight_val, 1.0370520)

        # the ratio of these values is > 10% away from 1.0, do not check that ratio

        # the sault image value should be the same as the image value
        # possibly AlmostEqual should be used here eventually, but at the
        # moment they are identical.
        self.assertEqual(fn_image_val, sw_image_val)
 
####    Suite: Required for CASA5     ####
def suite():
    return[linear_mosaic]
  
####    Imports     ####
if __name__ == '__main__':
    unittest.main()
 

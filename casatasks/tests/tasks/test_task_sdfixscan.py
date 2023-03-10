########################################################################
# test_task_sdfixscan.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.single.sdfixscan.html
#
#
##########################################################################
import os
import shutil
import unittest

import numpy

from casatasks import sdfixscan
from casatools import ctsys
from casatools import image as iatool
from casatools import quanta

_ia = iatool()
qa = quanta()

#
# Unit test of sdfixscan task.
#

# Utility


def drop_stokes_axis(imagename, outimagename):
    """Create 3d image from 4d.

    The function intends to drop Stokes axis and assumes that
    the third axis is Stokes.
    """
    myia = iatool()
    myia.open(imagename)
    subimage = myia.subimage(outfile=outimagename, dropdeg=True,
                             keepaxes=[0, 1, 3])
    subimage.close()
    myia.close()


def drop_deg_axes(imagename, outimagename):
    """Create 2d image from 4d.

    The function intends to drop both Stokes and Spectral axis and assumes that
    the third axis is Spectral and fourth axis is Stokes.
    """
    myia = iatool()
    myia.open(imagename)
    subimage = myia.subimage(outfile=outimagename, dropdeg=True,
                             keepaxes=[0, 1])
    subimage.close()
    myia.close()

###
# Base class for sdfixscan unit test
###


class sdfixscan_unittest_base:
    """Base class for sdfixscan unit test."""

    taskname = 'sdfixscan'
    datapath = ctsys.resolve('unittest/sdfixscan/')

    def _checkfile(self, name):
        isthere = os.path.exists(name)
        self.assertEqual(isthere, True,
                         msg='output file %s was not created because of the task failure' % (name))

    def _check_shape(self, infile, outfile):
        self._checkfile(infile)
        self._checkfile(outfile)
        _ia.open(infile)
        inshape = _ia.shape()
        inaxistypes = _ia.coordsys().axiscoordinatetypes()
        _ia.close()
        _ia.open(outfile)
        outshape = _ia.shape()
        outaxistypes = _ia.coordsys().axiscoordinatetypes()
        _ia.close()

        self.assertEqual(len(inshape), len(outshape))
        self.assertTrue(numpy.all(inshape == outshape))
        self.assertEqual(inaxistypes, outaxistypes)

    def _flux(self, csys, ref):
        # see CAS-5779, images/Images/ImageStatistics.tcc
        axis_ra = csys.findaxis(axis=csys.findaxisbyname('ra'))['axisincoordinate']
        axis_dec = csys.findaxis(axis=csys.findaxisbyname('dec'))['axisincoordinate']
        print(axis_ra, axis_dec)
        units = csys.units()
        increments = csys.increment()['numeric']
        incr_ra = qa.quantity(abs(increments[axis_ra]), units[axis_ra])
        incr_dec = qa.quantity(abs(increments[axis_dec]), units[axis_dec])
        area = qa.convert(incr_ra, 'arcsec')['value'] * qa.convert(incr_dec, 'arcsec')['value']
        return area * ref['sum']

    def _checkstats(self, name, ref):
        self._checkfile(name)
        _ia.open(name)
        stats = _ia.statistics(list=True, verbose=True)

        # set 'flux' value to ref
        if 'flux' not in ref:
            ref['flux'] = self._flux(_ia.coordsys(), ref)

        _ia.close()

        for key in stats.keys():
            # for key in self.keys:
            message = 'statistics \'%s\' does not match: %s' % (key, str(stats[key]))
            if type(stats[key]) == str:
                self.assertEqual(stats[key], ref[key],
                                 msg=message)
            else:
                # print stats[key]-ref[key]
                ret = numpy.allclose(stats[key], ref[key])
                self.assertEqual(ret, True,
                                 msg=message)

###
# Test on bad parameter settings
###


class sdfixscan_test0(unittest.TestCase, sdfixscan_unittest_base):
    """Test on bad parameter setting."""

    # Input and output names
    rawfiles = ['scan_x.im', 'scan_y.im']
    prefix = sdfixscan_unittest_base.taskname + 'Test0'
    outfile = prefix + '.im'

    def setUp(self):
        self.res = None
        for name in self.rawfiles:
            if (not os.path.exists(name)):
                shutil.copytree(os.path.join(self.datapath, name), name)

    def tearDown(self):
        for name in self.rawfiles:
            if (os.path.exists(name)):
                shutil.rmtree(name)
        os.system('rm -rf ' + self.prefix + '*')

    def test000(self):
        """Test 000: Default parameters."""
        # casatasks throw exception, CASA 5 tasks return False on failure
        self.assertRaises(Exception, sdfixscan)

    def test001(self):
        """Test 001: only 1 image is given for Basket-Weaving."""
        try:
            sdfixscan(infiles=[self.rawfiles[0]], mode='fft_mask', direction=[0.])
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos = str(e).find('infiles should be a list of input images for Basket-Weaving.')
            self.assertNotEqual(pos, -1,
                                msg='Unexpected exception was thrown: %s' % (str(e)))

    def test002(self):
        """Test 002: direction is not given for Basket-Weaving."""
        try:
            sdfixscan(infiles=self.rawfiles, mode='fft_mask')
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos = str(e).find('direction must have at least two different direction.')
            self.assertNotEqual(pos, -1,
                                msg='Unexpected exception was thrown: %s' % (str(e)))

    def test003(self):
        """Test 003: Multiple images are given for Press."""
        try:
            sdfixscan(infiles=self.rawfiles, mode='model')
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos = str(e).find('infiles allows only one input file for pressed-out method.')
            self.assertNotEqual(pos, -1,
                                msg='Unexpected exception was thrown: %s' % (str(e)))

    def test004(self):
        """Test 004: direction is not given for Press."""
        try:
            sdfixscan(infiles=[self.rawfiles[0]], mode='model')
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos = str(e).find('direction allows only one direction for pressed-out method.')
            self.assertNotEqual(pos, -1,
                                msg='Unexpected exception was thrown: %s' % (str(e)))

    def test005(self):
        """Test 005: Existing output image file."""
        shutil.copytree(os.path.join(self.datapath, self.rawfiles[0]), self.outfile)
        self.assertRaises(Exception, sdfixscan, infiles=self.rawfiles, mode='fft_mask',
                          direction=[0., 90.0], outfile=self.outfile, overwrite=False)

    def test006(self):
        """Test 006: Zero beamsize for Press."""
        self.assertRaises(Exception, sdfixscan, infiles=[self.rawfiles[0]],
                          mode='model', beamsize=0.0, direction=[0.], outfile=self.outfile,
                          overwrite=True)

###
# Test on Pressed method
###


class sdfixscan_test1(unittest.TestCase, sdfixscan_unittest_base):
    """Test on Pressed method.

    Test data, scan_x.im, is artificial data, which is

       - 128x128 in R.A. and Dec.
       - 1 polarization component (Stokes I)
       - 1 spectral channel
       - Gaussian source in the center
       - 1% random noise
       - scanning noise in horizontal direction
       - smoothed by Gaussian kernel of FWHM of 10 pixel

    """

    # Input and output names
    rawfile = 'scan_x.im'
    rawfilemod = rawfile.replace('.im', '_mod.im')
    prefix = sdfixscan_unittest_base.taskname + 'Test1'
    outfile = prefix + '.im'
    mode = 'model'

    def setUp(self):
        self.res = None
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(os.path.join(self.datapath, self.rawfile), self.rawfile)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.rawfilemod)):
            shutil.rmtree(self.rawfilemod)
        os.system('rm -rf ' + self.prefix + '*')

    def test100(self):
        """Test 100: Pressed method using whole pixels."""
        res = sdfixscan(infiles=self.rawfile, mode=self.mode, numpoly=0, beamsize=300.0,
                        smoothsize=2.0, direction=0.0, outfile=self.outfile, overwrite=True)
        self.assertEqual(res, None,
                         msg='Any error occurred during imaging')
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.76603365]),
                    'maxpos': numpy.array([63, 63, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:47.944, +01.03.00.212, I, 1.415e+09Hz',
                    'mean': numpy.array([1.77990955e-11]),
                    'medabsdevmed': numpy.array([0.00867557]),
                    'median': numpy.array([-0.0022334]),
                    'min': numpy.array([-0.17581469]),
                    'minpos': numpy.array([25, 64, 0, 0], dtype=numpy.int32),
                    'minposf': '23:58:19.982, +01.04.00.222, I, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01854331]),
                    'rms': numpy.array([0.09749392]),
                    'sigma': numpy.array([0.09749689]),
                    'sum': numpy.array([2.91620381e-07]),
                    'sumsq': numpy.array([155.73097211]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}
        self._check_shape(self.rawfile, self.outfile)
        self._checkstats(self.outfile, refstats)

    def test101(self):
        """Test 101: Pressed method with certain threshold."""
        res = sdfixscan(infiles=self.rawfile, mode=self.mode, numpoly=2, beamsize=300.0,
                        smoothsize=2.0, direction=0.0, tmax=0.1, tmin=-0.1,
                        outfile=self.outfile, overwrite=True)
        self.assertEqual(res, None,
                         msg='Any error occurred during imaging')
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.91835594]),
                    'maxpos': numpy.array([63, 63, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:47.944, +01.03.00.212, I, 1.415e+09Hz',
                    'mean': numpy.array([0.02586879]),
                    'medabsdevmed': numpy.array([0.00599332]),
                    'median': numpy.array([0.00104136]),
                    'min': numpy.array([-0.03498788]),
                    'minpos': numpy.array([55, 0, 0, 0], dtype=numpy.int32),
                    'minposf': '23:56:19.991, +00.00.00.000, I, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01204253]),
                    'rms': numpy.array([0.11026896]),
                    'sigma': numpy.array([0.10719492]),
                    'sum': numpy.array([423.83432425]),
                    'sumsq': numpy.array([199.21705095]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}
        self._check_shape(self.rawfile, self.outfile)
        self._checkstats(self.outfile, refstats)

    def test102(self):
        """Test 102: Test mask in pressed method using whole pixels."""
        # add spurious to image and mask the spurious.
        my_ia = iatool()
        my_ia.open(self.rawfile)
        try:
            maxval = my_ia.statistics()['max'][0]
            data = my_ia.getchunk()
            data[15:25, 15:25, :, :] = 100. * maxval
            my_ia.putchunk(data)
            del data
            my_ia.calcmask("mask(%s) && '%s'<%f" % (self.rawfile, self.rawfile, 10. * maxval),
                           name="sprious", asdefault=True)
            mask_in = my_ia.getchunk(getmask=True)
            mpix = numpy.where(mask_in == False)
            self.assertEqual(len(mpix[0]), 100, "Failed to set mask properly at pre-processing.")
        except Exception:
            raise
        finally:
            my_ia.close()
        del my_ia
        # Task execution
        res = sdfixscan(infiles=self.rawfile, mode=self.mode, numpoly=0, beamsize=300.0,
                        smoothsize=2.0, direction=0.0, outfile=self.outfile, overwrite=True)
        # Test results
        self.assertEqual(res, None,
                         msg='Any error occurred during imaging')
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.75250125]),
                    'maxpos': numpy.array([63, 63, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:47.944, +01.03.00.212, I, 1.415e+09Hz',
                    'mean': numpy.array([-0.00220275]),
                    'medabsdevmed': numpy.array([0.00883335]),
                    'median': numpy.array([-0.00232762]),
                    'min': numpy.array([-0.18934998]),
                    'minpos': numpy.array([25, 64, 0, 0], dtype=numpy.int32),
                    'minposf': '23:58:19.982, +01.04.00.222, I, 1.415e+09Hz',
                    'npts': numpy.array([16284.]),
                    'quartile': numpy.array([0.01897511]),
                    'rms': numpy.array([0.09790624]),
                    'sigma': numpy.array([0.09788447]),
                    'sum': numpy.array([-35.86960435]),
                    'sumsq': numpy.array([156.09243727]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}

        self._check_shape(self.rawfile, self.outfile)
        self._checkstats(self.outfile, refstats)
        _ia.open(self.outfile)
        mask_out = _ia.getchunk(getmask=True)
        _ia.close()
        self.assertTrue((mask_out == mask_in).all(), "Unexpected mask in output image.")

    def test100_3d(self):
        """Test 100_3d: Pressed method using whole pixels for 3D image."""
        drop_stokes_axis(self.rawfile, self.rawfilemod)
        self.assertTrue(os.path.exists(self.rawfilemod))
        res = sdfixscan(infiles=self.rawfilemod, mode=self.mode, numpoly=0, beamsize=300.0,
                        smoothsize=2.0, direction=0.0, outfile=self.outfile, overwrite=True)
        self.assertEqual(res, None,
                         msg='Any error occurred during imaging')
        refstats = {'blc': numpy.array([0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, 1.415e+09Hz',
                    'max': numpy.array([0.76603365]),
                    'maxpos': numpy.array([63, 63, 0], dtype=numpy.int32),
                    'maxposf': '23:55:47.944, +01.03.00.212, 1.415e+09Hz',
                    'mean': numpy.array([1.77990955e-11]),
                    'medabsdevmed': numpy.array([0.00867557]),
                    'median': numpy.array([-0.0022334]),
                    'min': numpy.array([-0.17581469]),
                    'minpos': numpy.array([25, 64, 0], dtype=numpy.int32),
                    'minposf': '23:58:19.982, +01.04.00.222, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01854331]),
                    'rms': numpy.array([0.09749392]),
                    'sigma': numpy.array([0.09749689]),
                    'sum': numpy.array([2.91620381e-07]),
                    'sumsq': numpy.array([155.73097211]),
                    'trc': numpy.array([127, 127, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, 1.415e+09Hz'}
        self._check_shape(self.rawfilemod, self.outfile)
        self._checkstats(self.outfile, refstats)

    def test100_2d(self):
        """Test 100_2d: Pressed method using whole pixels for 2D image."""
        drop_deg_axes(self.rawfile, self.rawfilemod)
        self.assertTrue(os.path.exists(self.rawfilemod))
        res = sdfixscan(infiles=self.rawfilemod, mode=self.mode, numpoly=0, beamsize=300.0,
                        smoothsize=2.0, direction=0.0, outfile=self.outfile, overwrite=True)
        self.assertEqual(res, None,
                         msg='Any error occurred during imaging')
        refstats = {'blc': numpy.array([0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000',
                    'max': numpy.array([0.76603365]),
                    'maxpos': numpy.array([63, 63], dtype=numpy.int32),
                    'maxposf': '23:55:47.944, +01.03.00.212',
                    'mean': numpy.array([1.77990955e-11]),
                    'medabsdevmed': numpy.array([0.00867557]),
                    'median': numpy.array([-0.0022334]),
                    'min': numpy.array([-0.17581469]),
                    'minpos': numpy.array([25, 64], dtype=numpy.int32),
                    'minposf': '23:58:19.982, +01.04.00.222',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01854331]),
                    'rms': numpy.array([0.09749392]),
                    'sigma': numpy.array([0.09749689]),
                    'sum': numpy.array([2.91620381e-07]),
                    'sumsq': numpy.array([155.73097211]),
                    'trc': numpy.array([127, 127], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734'}
        self._check_shape(self.rawfilemod, self.outfile)
        self._checkstats(self.outfile, refstats)

###
# Test on FFT based Basket-Weaving
###


class sdfixscan_test2(unittest.TestCase, sdfixscan_unittest_base):
    """Test on FFT based Basket-Weaving.

    Test data, scan_x.im and scan_y.im, are artificial data, which is

       - 128x128 in R.A. and Dec.
       - 1 polarization component (Stokes I)
       - 1 spectral channel
       - Gaussian source in the center
       - 1% random noise
       - scanning noise in horizontal direction (scan_x.im)
         or vertical direction (scan_y.im)
       - smoothed by Gaussian kernel of FWHM of 10 pixel

    """

    # Input and output names
    rawfiles = ['scan_x.im', 'scan_y.im']
    rawfilesmod = list(map(lambda x: x.replace('.im', '_mod.im'), rawfiles))
    prefix = sdfixscan_unittest_base.taskname + 'Test2'
    outfile = prefix + '.im'
    mode = 'fft_mask'

    def setUp(self):
        self.res = None
        for name in self.rawfiles:
            if os.path.exists(name):
                shutil.rmtree(name)
            shutil.copytree(os.path.join(self.datapath, name), name)

    def tearDown(self):
        for name in self.rawfiles:
            if (os.path.exists(name)):
                shutil.rmtree(name)
        for name in self.rawfilesmod:
            if (os.path.exists(name)):
                shutil.rmtree(name)
        os.system('rm -rf ' + self.prefix + '*')

    def test200(self):
        """Test 200: FFT based Basket-Weaving using whole pixels."""
        sdfixscan(infiles=self.rawfiles, mode=self.mode, direction=[0.0, 90.0],
                  maskwidth=20.0, outfile=self.outfile, overwrite=True)
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.92863309]),
                    'maxpos': numpy.array([64, 64, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:43.941, +01.04.00.222, I, 1.415e+09Hz',
                    'mean': numpy.array([0.02962625]),
                    'medabsdevmed': numpy.array([0.00571492]),
                    'median': numpy.array([0.00429045]),
                    'min': numpy.array([-0.02551341]),
                    'minpos': numpy.array([56, 107, 0, 0], dtype=numpy.int32),
                    'minposf': '23:56:15.881, +01.47.01.037, I, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01154788]),
                    'rms': numpy.array([0.11279197]),
                    'sigma': numpy.array([0.1088349]),
                    'sum': numpy.array([485.39648429]),
                    'sumsq': numpy.array([208.43769809]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}
        self._check_shape(self.rawfiles[0], self.outfile)
        self._checkstats(self.outfile, refstats)

    def test201(self):
        """Test 201: FFT based Basket-Weaving with certain threshold."""
        sdfixscan(infiles=self.rawfiles, mode=self.mode, direction=[0.0, 90.0],
                  maskwidth=20.0, tmax=0.5, tmin=-0.1, outfile=self.outfile, overwrite=True)
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.99608284]),
                    'maxpos': numpy.array([63, 63, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:47.944, +01.03.00.212, I, 1.415e+09Hz',
                    'mean': numpy.array([0.02962625]),
                    'medabsdevmed': numpy.array([0.00570825]),
                    'median': numpy.array([0.00428429]),
                    'min': numpy.array([-0.02518307]),
                    'minpos': numpy.array([56, 107, 0, 0], dtype=numpy.int32),
                    'minposf': '23:56:15.881, +01.47.01.037, I, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01155156]),
                    'rms': numpy.array([0.113305]),
                    'sigma': numpy.array([0.10936653]),
                    'sum': numpy.array([485.39650849]),
                    'sumsq': numpy.array([210.3381585]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}
        self._check_shape(self.rawfiles[0], self.outfile)
        self._checkstats(self.outfile, refstats)

    def test202(self):
        """Test 202: Test mask in FFT based Basket-Weaving using whole pixels."""
        # add spurious to images and mask the spurious pixels.
        spix = [10, 15]
        epix = [20, 25]
        mask_in = []
        my_ia = iatool()
        for i in range(len(self.rawfiles)):
            name = self.rawfiles[i]
            s = spix[i % len(spix)]
            e = epix[i % len(epix)]
            my_ia.open(name)
            try:
                maxval = my_ia.statistics()['max'][0]
                data = my_ia.getchunk()
                data[s:e, s:e, :, :] = 100. * maxval
                my_ia.putchunk(data)
                del data
                my_ia.calcmask("mask(%s) && '%s'<%f" % (name, name, 10. * maxval),
                               name="sprious", asdefault=True)
                mask_in.append(my_ia.getchunk(getmask=True))
            except Exception:
                raise
            finally:
                my_ia.done()
            mpix = numpy.where(mask_in[i] == False)
            self.assertEqual(len(mpix[0]), (e - s)**2,
                             "Failed to set mask properly at pre-processing.")
        del my_ia
        mask_ref = mask_in[0]
        for msk in mask_in:
            mask_ref += msk
        del mask_in
        # Task execution
        sdfixscan(infiles=self.rawfiles, mode=self.mode, direction=[0.0, 90.0],
                  maskwidth=20.0, outfile=self.outfile, overwrite=True)
        # Test results
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.92862648]),
                    'maxpos': numpy.array([64, 64, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:43.941, +01.04.00.222, I, 1.415e+09Hz',
                    'mean': numpy.array([0.02969963]),
                    'medabsdevmed': numpy.array([0.00574234]),
                    'median': numpy.array([0.00439214]),
                    'min': numpy.array([-0.02554246]),
                    'minpos': numpy.array([56, 107, 0, 0], dtype=numpy.int32),
                    'minposf': '23:56:15.881, +01.47.01.037, I, 1.415e+09Hz',
                    'npts': numpy.array([16359.]),
                    'quartile': numpy.array([0.011633]),
                    'rms': numpy.array([0.11288381]),
                    'sigma': numpy.array([0.10891011]),
                    'sum': numpy.array([485.85626113]),
                    'sumsq': numpy.array([208.45871511]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}

        # print imstat(self.outfile)
        self._check_shape(self.rawfiles[0], self.outfile)
        self._checkstats(self.outfile, refstats)
        _ia.open(self.outfile)
        mask_out = _ia.getchunk(getmask=True)
        _ia.close()
        self.assertTrue((mask_out == mask_ref).all(), "Unexpected mask in output image.")

    def test203(self):
        """Test 203: test for len(infiles) > len(direction)."""
        infiles = self.rawfiles + self.rawfiles
        sdfixscan(infiles=infiles, mode=self.mode, direction=[0.0, 90.0],
                  maskwidth=20.0, outfile=self.outfile, overwrite=True)
        refstats = {'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, I, 1.415e+09Hz',
                    'max': numpy.array([0.92863309]),
                    'maxpos': numpy.array([64, 64, 0, 0], dtype=numpy.int32),
                    'maxposf': '23:55:43.941, +01.04.00.222, I, 1.415e+09Hz',
                    'mean': numpy.array([0.02962625]),
                    'medabsdevmed': numpy.array([0.00571492]),
                    'median': numpy.array([0.00429045]),
                    'min': numpy.array([-0.02551341]),
                    'minpos': numpy.array([56, 107, 0, 0], dtype=numpy.int32),
                    'minposf': '23:56:15.881, +01.47.01.037, I, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01154788]),
                    'rms': numpy.array([0.11279196]),
                    'sigma': numpy.array([0.1088349]),
                    'sum': numpy.array([485.39648429]),
                    'sumsq': numpy.array([208.4376955]),
                    'trc': numpy.array([127, 127, 0, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, I, 1.415e+09Hz'}
        self._check_shape(self.rawfiles[0], self.outfile)
        self._checkstats(self.outfile, refstats)

    def test200_3d(self):
        """Test 200_3d: FFT based Basket-Weaving using whole pixels for 3D image."""
        for infile, outfile in zip(self.rawfiles, self.rawfilesmod):
            drop_stokes_axis(infile, outfile)
            self.assertTrue(os.path.exists(outfile))
        sdfixscan(infiles=self.rawfilesmod, mode=self.mode, direction=[0.0, 90.0],
                  maskwidth=20.0, outfile=self.outfile, overwrite=True)
        refstats = {'blc': numpy.array([0, 0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000, 1.415e+09Hz',
                    'max': numpy.array([0.92863309]),
                    'maxpos': numpy.array([64, 64, 0], dtype=numpy.int32),
                    'maxposf': '23:55:43.941, +01.04.00.222, 1.415e+09Hz',
                    'mean': numpy.array([0.02962625]),
                    'medabsdevmed': numpy.array([0.00571492]),
                    'median': numpy.array([0.00429045]),
                    'min': numpy.array([-0.02551341]),
                    'minpos': numpy.array([56, 107, 0], dtype=numpy.int32),
                    'minposf': '23:56:15.881, +01.47.01.037, 1.415e+09Hz',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01154788]),
                    'rms': numpy.array([0.11279197]),
                    'sigma': numpy.array([0.1088349]),
                    'sum': numpy.array([485.39648429]),
                    'sumsq': numpy.array([208.43769809]),
                    'trc': numpy.array([127, 127, 0], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734, 1.415e+09Hz'}
        self._check_shape(self.rawfilesmod[0], self.outfile)
        self._checkstats(self.outfile, refstats)

    def test200_2d(self):
        """Test 200_2d: FFT based Basket-Weaving using whole pixels for 2D image."""
        for infile, outfile in zip(self.rawfiles, self.rawfilesmod):
            drop_deg_axes(infile, outfile)
            self.assertTrue(os.path.exists(outfile))
        sdfixscan(infiles=self.rawfilesmod, mode=self.mode, direction=[0.0, 90.0],
                  maskwidth=20.0, outfile=self.outfile, overwrite=True)
        refstats = {'blc': numpy.array([0, 0], dtype=numpy.int32),
                    'blcf': '00:00:00.000, +00.00.00.000',
                    'max': numpy.array([0.92863309]),
                    'maxpos': numpy.array([64, 64], dtype=numpy.int32),
                    'maxposf': '23:55:43.941, +01.04.00.222',
                    'mean': numpy.array([0.02962625]),
                    'medabsdevmed': numpy.array([0.00571492]),
                    'median': numpy.array([0.00429045]),
                    'min': numpy.array([-0.02551341]),
                    'minpos': numpy.array([56, 107], dtype=numpy.int32),
                    'minposf': '23:56:15.881, +01.47.01.037',
                    'npts': numpy.array([16384.]),
                    'quartile': numpy.array([0.01154788]),
                    'rms': numpy.array([0.11279197]),
                    'sigma': numpy.array([0.1088349]),
                    'sum': numpy.array([485.39648429]),
                    'sumsq': numpy.array([208.43769809]),
                    'trc': numpy.array([127, 127], dtype=numpy.int32),
                    'trcf': '23:51:31.537, +02.07.01.734'}
        self._check_shape(self.rawfilesmod[0], self.outfile)
        self._checkstats(self.outfile, refstats)


if __name__ == '__main__':
    unittest.main()

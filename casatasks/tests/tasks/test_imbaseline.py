import os
import re
import shutil
import unittest

import functools
import numpy as np
from casatasks.private.sdutil import tool_manager
from casatasks.private.task_imbaseline \
    import (CasaImageStack, EraseableFolder, Image2MSMethods, Image2MSParams, ImageShape, ImageSubtractionMethods,
            ImsmoothMethods, ImsmoothParams, MeasurementSetStack, MS2ImageMethods, SdbaselineMethods,
            SdbaselineParams, SdsmoothMethods, SdsmoothParams, UnerasableFolder,
            eraseable_folder_register, get_image_shape, imbaseline)

from casatools import ctsys, image, table

_tb = table()
ctsys_resolve = ctsys.resolve
DATACOLUMN = 'DATA'


class test_base(unittest.TestCase):

    @staticmethod
    def invalid_argument_case(func):
        """Decorator for the test case that is intended to fail due to invalid argument."""
        @functools.wraps(func)
        def wrapper(self):
            func(self)
            self.assertFalse(self.result, msg='The task must return False')
        return wrapper

    @staticmethod
    def exception_case(exception_type, exception_pattern):
        """Decorator for the test case that is intended to throw exception.

            exception_type: type of exception
            exception_pattern: regex for inspecting exception message
                               using re.search
        """
        def wrapper(func):
            @functools.wraps(func)
            def _wrapper(self):
                self.assertTrue(len(exception_pattern) > 0, msg='Internal Error')
                with self.assertRaises(exception_type) as ctx:
                    func(self)
                    self.fail(msg='The task must throw exception')
                the_exception = ctx.exception
                message = str(the_exception)
                self.assertIsNotNone(re.search(exception_pattern, message),
                                     msg='error message \'%s\' is not expected.' % (message))
            return _wrapper
        return wrapper

    def tearDown(self):
        eraseable_folder_register.clear()
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

    def _copy_test_files(self, basename, filename):
        """Copy files for testing into current path."""
        _base = ctsys_resolve(basename)
        copy_from = os.path.join(_base, filename)
        if not os.path.exists(copy_from) or copy_from == os.path.join(os.getcwd(), filename):
            raise RuntimeError(f'Error is occured or existed on a path {copy_from} or {filename}')

        if os.path.exists(filename):
            shutil.rmtree(filename)
        os.system('cp -RH ' + os.path.join(_base, filename) + ' ' + filename)

    def _create_image(self, datapath, val=1, shape=[0, 0, 0, 0]):
        _ia = image()
        ary = _ia.makearray(v=val, shape=shape)
        _ia.fromarray(outfile=datapath, pixels=ary, overwrite=True)
        _ia.done()


class AbstractFileStack_test(test_base):
    """Test AbstractFileStack / (Un)EraseableFolder.

    1-1. Create Stack with exist file
    1-2. Create Stack with unexist file
    1-3. push() exist file
    1-4. push() unexist file
    1-5. pop() exist stuff
    1-6. pop() unexist stuff
    1-7. peak() exist stuff
    1-8. peak() unexist stuff
    1-9. subpeak() exist stuff
    1-10. subpeak() unexist stuff
    1-11. bottom() exist stuff
    1-12. bottom() unexist stuff
    1-13. clear() exist EraseableFolder file
    1-14. clear() exist UneraseableFolder file
    1-15. erase EraseableFolder file
    1-16. erase UneraseableFolder file
    """

    dummy_folder1 = 'dummy1'
    dummy_folder2 = 'dummy2'
    unexist_folder = 'unexists'

    def setUp(self):
        if os.path.exists(self.dummy_folder1):
            shutil.rmtree(self.dummy_folder1)
        if os.path.exists(self.dummy_folder2):
            shutil.rmtree(self.dummy_folder2)
        os.system(f'mkdir {self.dummy_folder1}')
        os.system(f'mkdir {self.dummy_folder2}')
        if os.path.exists(self.unexist_folder):
            shutil.rmtree(self.unexist_folder)

    def test_1_1(self):
        stack = CasaImageStack(UnerasableFolder(self.dummy_folder1))
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(ValueError, 'file unexists is not found')
    def test_1_2(self):
        CasaImageStack(UnerasableFolder(self.unexist_folder))

    def test_1_3(self):
        stack = CasaImageStack()
        stack.push(UnerasableFolder(self.dummy_folder1))
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(ValueError, 'file unexists is not found')
    def test_1_4(self):
        stack = CasaImageStack()
        stack.push(UnerasableFolder(self.unexist_folder))

    def test_1_5(self):
        stack = CasaImageStack()
        stack.push(UnerasableFolder(self.dummy_folder1))
        obj = UnerasableFolder(self.dummy_folder2)
        stack.push(obj)
        tmp = stack.pop()
        self.assertEqual(obj, tmp)
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(RuntimeError, 'the stack cannot pop')
    def test_1_6(self):
        stack = CasaImageStack()
        stack.pop()

    def test_1_7(self):
        stack = CasaImageStack()
        obj = UnerasableFolder(self.dummy_folder1)
        stack.push(obj)
        self.assertEqual(stack.peak(), obj)

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_1_8(self):
        stack = CasaImageStack()
        stack.peak()

    def test_1_9(self):
        obj = UnerasableFolder(self.dummy_folder1)
        stack = CasaImageStack(obj)
        stack.push(UnerasableFolder(self.dummy_folder2))
        self.assertEqual(stack.subpeak(), obj)

    @test_base.exception_case(RuntimeError, 'the stack has only one stuff')
    def test_1_10(self):
        stack = CasaImageStack(UnerasableFolder(self.dummy_folder1))
        stack.subpeak()

    def test_1_11(self):
        obj = UnerasableFolder(self.dummy_folder1)
        stack = CasaImageStack(obj)
        self.assertEqual(stack.bottom(), obj)

    @test_base.exception_case(RuntimeError, 'the stack has not have enough stuff')
    def test_1_12(self):
        stack = CasaImageStack()
        stack.bottom()

    def test_1_13(self):
        file = EraseableFolder(self.dummy_folder1)
        stack = CasaImageStack(file)
        stack.clear(False)
        self.assertTrue(os.path.exists(self.dummy_folder1))
        self.assertEqual(stack.height(), 0)

    def test_1_14(self):
        stack = CasaImageStack(UnerasableFolder(self.dummy_folder1))
        stack.clear(False)
        self.assertTrue(os.path.exists(self.dummy_folder1))
        self.assertEqual(stack.height(), 0)

    def test_1_15(self):
        file = EraseableFolder(self.dummy_folder1)
        file.erase(False)
        eraseable_folder_register.pop(self.dummy_folder1)
        self.assertFalse(os.path.exists(self.dummy_folder1))

    def test_1_16(self):
        file = UnerasableFolder(self.dummy_folder1)
        file.erase(False)
        self.assertTrue(os.path.exists(self.dummy_folder1))


class ImageShape_test(test_base):
    """Test ImageShape.

    2-1. successful case
    2-2. invalid im_nchan
    2-3. invalid dir_shape
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_2_1(self):
        shape = ImageShape(im_shape=np.array([100, 100, 1, 100]), axis_dir=np.array([0, 1]), axis_sp=3, axis_pol=2)
        shape.validate()
        shape = ImageShape(im_shape=np.array([100, 100, 100, 1]), axis_dir=np.array([0, 1]), axis_sp=2, axis_pol=3)
        shape.validate()
        # any exceptions are not thrown, its OK

    @test_base.exception_case(ValueError, 'nchan \\d is too few to perform baseline subtraction')
    def test_2_2(self):
        shape = ImageShape(im_shape=np.array([100, 100, 1, 1]), axis_dir=np.array([0, 1]), axis_sp=3, axis_pol=2)
        shape.validate()

    @test_base.exception_case(ValueError, 'invalid value: dir_shape \\[\\d+\\]')
    def test_2_3(self):
        shape = ImageShape(np.array([100, 100, 1, 100]), axis_dir=np.array([0]), axis_sp=3, axis_pol=2)
        shape.validate()


class imsmooth_test(test_base):
    """Test imsmooth execution.

    Tests of imsmooth rely on ones of test_imsmooth basically, so we have minimal tests in imbaseline.

    3-1. simple successful case
    3-2. simple failure case
    3-3. check parameters of ImsmoothParams
    """

    datapath = ctsys_resolve('unittest/imsmooth/')
    tiny = 'tiny.im'

    def setUp(self):
        self._copy_test_files(self.datapath, self.tiny)

    def test_3_1(self):
        major = '2.5arcsec'
        minor = '2arcsec'
        pa = '0deg'
        dirkernel = 'gaussian'
        kimage = ''
        scale = -1

        _stack = CasaImageStack(top=UnerasableFolder(self.tiny))

        ImsmoothMethods.execute(dirkernel, major, minor, pa, kimage, scale, _stack)
        self.assertTrue(os.path.exists(_stack.peak().path))

        _stack.pop()
        ImsmoothMethods.execute(dirkernel, major, minor, pa, kimage, scale, _stack)
        self.assertTrue(os.path.exists(_stack.peak().path))

        _stack.clear()
        try:
            _stack.push(UnerasableFolder('nonexists'))
        except Exception:
            pass
        self.assertFalse(_stack.height() > 0)

    @test_base.exception_case(ValueError, 'Unsupported direction smoothing kernel, foobar')
    def test_3_2(self):
        major = '2.5arcsec'
        minor = '2arcsec'
        pa = '0deg'
        dirkernel = 'foobar'
        kimage = ''
        scale = -1

        _stack = CasaImageStack(top=UnerasableFolder(self.tiny))

        ImsmoothMethods.execute(dirkernel, major, minor, pa, kimage, scale, _stack)

    def test_3_3(self):
        targetres = stretch = False
        mask = region = box = chans = stokes = ''
        beam = {}
        infile = 'infile'
        outfile = 'outfile'
        kernel = ('none', 'image', 'gaussian', 'boxcar')
        major = '2.5arcsec'
        minor = '2arcsec'
        pa = '0deg'
        kimage = self.tiny
        scale = -2.0
        logorigin = 'imbaseline'

        # none
        valid_param = dict(targetres=targetres, mask=mask, beam=beam, region=region, box=box, chans=chans, stokes=stokes,
                           stretch=stretch, overwrite=True, imagename=infile, outfile=outfile, kernel=kernel[0], major=major,
                           minor=minor, pa=pa, kimage=kimage, scale=scale, __log_origin=logorigin)
        param = ImsmoothParams(infile, outfile, kernel[0], major, minor, pa, kimage, scale)
        param.validate()
        self.assertEqual(param(), valid_param)

        # image
        minor = major = pa = ''
        valid_param = dict(targetres=targetres, mask=mask, beam=beam, region=region, box=box, chans=chans, stokes=stokes,
                           stretch=stretch, overwrite=True, imagename=infile, outfile=outfile, kernel=kernel[1], major=major,
                           minor=minor, pa=pa, kimage=kimage, scale=scale, __log_origin=logorigin)
        param = ImsmoothParams(infile, outfile, kernel[1], major, minor, pa, kimage, scale)
        param.validate()
        self.assertEqual(param(), valid_param)

        # gaussian
        major = '2.5arcsec'
        minor = '2arcsec'
        pa = '0deg'
        kimage = ''
        scale = -1.0
        valid_param = dict(targetres=targetres, mask=mask, beam=beam, region=region, box=box, chans=chans, stokes=stokes,
                           stretch=stretch, overwrite=True, imagename=infile, outfile=outfile, kernel=kernel[2], major=major,
                           minor=minor, pa=pa, kimage=kimage, scale=scale, __log_origin=logorigin)
        param = ImsmoothParams(infile, outfile, kernel[2], major, minor, pa, kimage, scale)
        param.validate()
        self.assertEqual(param(), valid_param)

        # boxcar
        valid_param = dict(targetres=targetres, mask=mask, beam=beam, region=region, box=box, chans=chans, stokes=stokes,
                           stretch=stretch, overwrite=True, imagename=infile, outfile=outfile, kernel=kernel[3], major=major,
                           minor=minor, pa=pa, kimage=kimage, scale=scale, __log_origin=logorigin)
        param = ImsmoothParams(infile, outfile, kernel[3], major, minor, pa, kimage, scale)
        param.validate()
        self.assertEqual(param(), valid_param)


class image2ms_test(test_base):
    """Test image2ms.

    4-1. simple successful case
    4-2. invalid datacolumn
    4-3. invalid image
    4-4. set empty stack
    4-5. check Image2MSParams
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected = 'expected.im'
    dummy_folder1 = 'dummy1'
    datacolumn = DATACOLUMN

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected)
        self.image_shape = get_image_shape(os.path.join(self.datapath, self.expected))
        if os.path.exists(self.dummy_folder1):
            shutil.rmtree(self.dummy_folder1)
        os.system(f'mkdir {self.dummy_folder1}')

    def test_4_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected))
        ms_stack = MeasurementSetStack()
        Image2MSMethods.execute(self.datacolumn, self.image_shape, image_stack, ms_stack)
        self.assertEqual(ms_stack.height(), 1)
        ms_path = ms_stack.peak().path
        self.assertTrue(os.path.exists(ms_path))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'table.dat')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'ANTENNA')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'DATA_DESCRIPTION')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'FEED')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'FIELD')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'FLAG_CMD')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'HISTORY')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'OBSERVATION')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'POINTING')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'POLARIZATION')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'PROCESSOR')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'SOURCE')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'SPECTRAL_WINDOW')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'STATE')))

    @test_base.exception_case(RuntimeError, 'column INVALID does not exist')
    def test_4_2(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected))
        ms_stack = MeasurementSetStack()
        Image2MSMethods.execute('INVALID', self.image_shape, image_stack, ms_stack)

    @test_base.exception_case(RuntimeError, 'Unable to open image dummy1.')
    def test_4_3(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.dummy_folder1))
        ms_stack = MeasurementSetStack()
        Image2MSMethods.execute(self.datacolumn, self.image_shape, image_stack, ms_stack)

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_4_4(self):
        image_stack = CasaImageStack()
        ms_stack = MeasurementSetStack()
        Image2MSMethods.execute(self.datacolumn, self.image_shape, image_stack, ms_stack)

    def test_4_5(self):
        outfile = 'output_4_5.ms'
        params = Image2MSParams(self.expected, outfile, self.datacolumn, self.image_shape)
        params.validate()
        self.assertEqual(params.infile, self.expected)
        self.assertEqual(params.outfile, outfile)
        self.assertTrue(np.all(params.im_shape == self.image_shape.im_shape))
        self.assertTrue(np.all(params.axis_dir == self.image_shape.axis_dir))
        self.assertTrue(np.all(params.dir_shape == self.image_shape.dir_shape))
        self.assertEqual(params.axis_sp, self.image_shape.axis_sp)
        self.assertEqual(params.axis_pol, self.image_shape.axis_pol)
        self.assertEqual(params.im_nrow, self.image_shape.im_nrow)
        self.assertEqual(params.im_nchan, self.image_shape.im_nchan)
        self.assertEqual(params.im_npol, self.image_shape.im_npol)


class sdsmooth_test(test_base):
    """Test sdsmooth execution.

    Tests of sdsmooth rely on ones of test_sdsmooth basically, so we have minimal tests in imbaseline.

    5-1. simple successful case
    5-2. invalid ms stack
    5-3. invalid image stack
    5-4. check SdsmoothParams
    """
    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = 'expected.im'
    expected_ms = 'expected.ms'
    dummy_folder1 = 'dummy1'
    datacolumn = DATACOLUMN
    spkenel = 'gaussian'
    kwidth = 5

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_ms)
        self.image_shape = get_image_shape(os.path.join(self.datapath, self.expected_im))

    def test_5_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        SdsmoothMethods.execute(self.datacolumn, self.spkenel, self.kwidth, image_stack, ms_stack, self.image_shape)
        self.assertEqual(image_stack.height(), 2)
        self.assertEqual(ms_stack.height(), 2)
        self.assertTrue(os.path.exists(os.path.join(image_stack.peak().path, 'table.dat')))

        ms_path = ms_stack.peak().path
        self.assertTrue(os.path.exists(ms_path))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'table.dat')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'ANTENNA')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'DATA_DESCRIPTION')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'FEED')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'FIELD')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'FLAG_CMD')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'HISTORY')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'OBSERVATION')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'POINTING')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'POLARIZATION')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'PROCESSOR')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'SOURCE')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'SPECTRAL_WINDOW')))
        self.assertTrue(os.path.exists(os.path.join(ms_path, 'STATE')))

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_5_2(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        SdsmoothMethods.execute(self.datacolumn, self.spkenel, self.kwidth, image_stack, ms_stack, self.image_shape)

    @test_base.exception_case(RuntimeError, 'the stack has not have enough stuff')
    def test_5_3(self):
        image_stack = CasaImageStack()
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        SdsmoothMethods.execute(self.datacolumn, self.spkenel, self.kwidth, image_stack, ms_stack, self.image_shape)

    def test_5_4(self):
        spw = field = antenna = timerange = scan = pol = intent = ''
        reindex = overwrite = True
        infile = 'infile'
        outfile = 'outfile'
        datacolumn = DATACOLUMN
        kernel = ('none', 'gaussian', 'boxcar')
        kwidth = 5
        logorigin = 'imbaseline'

        def compare_params(_kernel):
            valid_params = dict(spw=spw, field=field, antenna=antenna, timerange=timerange, scan=scan, pol=pol, intent=intent,
                                reindex=reindex, overwrite=overwrite, infile=infile, datacolumn=datacolumn, kernel=_kernel,
                                kwidth=kwidth, outfile=outfile, __log_origin=logorigin)
            params = SdsmoothParams(infile=infile, outfile=outfile, datacolumn=datacolumn, spkernel=_kernel, kwidth=kwidth)
            params.validate()
            self.assertEqual(params(), valid_params)

        [compare_params(_kernel) for _kernel in kernel]


class sdbaseline_test(test_base):
    """Test sdbaseline execution.

    Tests of sdbaseline rely on ones of test_sdbaseline basically, so we have minimal tests in imbaseline.

    6-1. simple successful case
    6-2. invalid ms stack
    6-3. invalid image stack
    6-4. check SdbaselineParams
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = 'expected.im'
    expected_ms = 'expected.ms'
    bloutput = 'test.csv'
    maskmode = 'auto'
    blparam = 'analytic_variable_blparam.txt'
    chans = ''
    thresh = 5.0
    avg_limit = 4
    minwidth = 4
    edge = [0, 0]
    blfunc = 'cspline'
    order = 5
    npiece = 1
    applyfft = True
    fftthresh = 3.0
    addwn = [0]
    rejwn = []
    clipniter = 10
    clipthresh = 2.0
    datacolumn = DATACOLUMN

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_ms)
        self._copy_test_files(self.datapath, self.blparam)
        self.image_shape = get_image_shape(os.path.join(self.datapath, self.expected_im))

    def test_6_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        SdbaselineMethods.execute(self.datacolumn, self.bloutput, self.maskmode, self.chans, self.thresh, self.avg_limit,
                                  self.minwidth, self.edge, self.blfunc, self.order, self.npiece, self.applyfft,
                                  self.fftthresh, self.addwn, self.rejwn, self.blparam, self.clipniter, self.clipthresh,
                                  image_stack, ms_stack, self.image_shape)
        self.assertTrue(os.path.exists(ms_stack.peak().path))
        self.assertTrue(os.path.exists(self.bloutput))
        self.assertTrue(os.path.exists(image_stack.peak().path))

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_6_2(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        SdbaselineMethods.execute(self.datacolumn, self.bloutput, self.maskmode, self.chans, self.thresh, self.avg_limit,
                                  self.minwidth, self.edge, self.blfunc, self.order, self.npiece, self.applyfft,
                                  self.fftthresh, self.addwn, self.rejwn, self.blparam, self.clipniter, self.clipthresh,
                                  image_stack, ms_stack, self.image_shape)

    @test_base.exception_case(RuntimeError, 'the stack has not have enough stuff')
    def test_6_3(self):
        image_stack = CasaImageStack()
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        SdbaselineMethods.execute(self.datacolumn, self.bloutput, self.maskmode, self.chans, self.thresh, self.avg_limit,
                                  self.minwidth, self.edge, self.blfunc, self.order, self.npiece, self.applyfft,
                                  self.fftthresh, self.addwn, self.rejwn, self.blparam, self.clipniter, self.clipthresh,
                                  image_stack, ms_stack, self.image_shape)

    def test_6_4(self):
        antenna = field = timerange = scan = pol = intent = bltable = ''
        reindex = dosubtract = overwrite = True
        updateweight = showprogress = verbose = False
        blmode = 'fit'
        blformat = 'csv'
        sigmavalue = 'stddev'
        minnrow = 1000
        fftmethod = 'fft'

        infile = 'infile'
        outfile = 'outfile'
        datacolumn = 'DATA'
        bloutput = 'bloutput'
        maskmode = ('list', 'auto')
        chans = ''
        thresh = 6.0
        avg_limit = 5
        minwidth = 5
        edge = [1, 1]
        blfunc = ('poly', 'chebyshev', 'cspline', 'sinusoid', 'variable')
        order = 6
        npiece = 2
        applyfft = False
        fftthresh = 4.0
        addwn = [1]
        rejwn = [1]
        blparam = self.blparam
        clipniter = 11
        clipthresh = 3.0
        logorigin = 'imbaseline'

        def compare_params(_maskmode, _blfunc):
            valid_param = dict(antenna=antenna, field=field, spw=chans, timerange=timerange, scan=scan, pol=pol, intent=intent,
                               reindex=reindex, blmode=blmode, dosubtract=dosubtract, blformat=blformat, bltable=bltable,
                               updateweight=updateweight, sigmavalue=sigmavalue, showprogress=showprogress, minnrow=minnrow,
                               fftmethod=fftmethod, verbose=verbose, overwrite=overwrite, infile=infile, datacolumn=datacolumn,
                               maskmode=_maskmode, thresh=thresh, avg_limit=avg_limit, minwidth=minwidth, edge=edge,
                               bloutput=bloutput, blfunc=_blfunc, order=order, npiece=npiece, applyfft=applyfft,
                               fftthresh=fftthresh, addwn=addwn, rejwn=rejwn, clipthresh=clipthresh, clipniter=clipniter,
                               blparam=blparam, outfile=outfile, __log_origin=logorigin)
            params = SdbaselineParams(infile=infile, outfile=outfile, datacolumn=datacolumn, bloutput=bloutput,
                                      maskmode=_maskmode, chans=chans, thresh=thresh, avg_limit=avg_limit, minwidth=minwidth,
                                      edge=edge, blfunc=_blfunc, order=order, npiece=npiece, applyfft=applyfft,
                                      fftthresh=fftthresh, addwn=addwn, rejwn=rejwn, blparam=blparam, clipniter=clipniter,
                                      clipthresh=clipthresh)
            params.validate()
            self.assertEqual(params(), valid_param)

        [compare_params(_maskmode, _blfunc) for _maskmode in maskmode for _blfunc in blfunc]


class image_subtraction_test(test_base):
    """Test image subtractions.

    7-1. successful test: output = input_image - (smoothed_image - smoothed_and_subtracted_image)
    7-2. successful test: output = subtracted_image
    7-3. unmatch shape
    7-4. unmatch shape(exception is not thrown)
    7-5. three images subtraction test
    7-6. two images subtraction test
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = 'expected.im'
    expected_imsmoothed = 'expected.imsmooth.im'
    expected_bl = 'expected.bl.im'
    testdata_01 = ('testdata_01.im', 1.5, [64, 64, 4, 128])
    testdata_02 = ('testdata_02.im', 2.0, [64, 64, 4, 128])
    testdata_03 = ('testdata_03.im', 2.5, [64, 64, 4, 128])
    testdata_err = ('testdata_err.im', 1, [65, 64, 4, 128])

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_imsmoothed)
        self._copy_test_files(self.datapath, self.expected_bl)
        self._create_image(self.testdata_01[0], self.testdata_01[1], self.testdata_01[2])
        self._create_image(self.testdata_02[0], self.testdata_02[1], self.testdata_02[2])
        self._create_image(self.testdata_03[0], self.testdata_03[1], self.testdata_03[2])
        self._create_image(self.testdata_err[0], self.testdata_err[1], self.testdata_err[2])

    def test_7_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        image_stack.push(EraseableFolder(self.expected_imsmoothed))
        image_stack.push(EraseableFolder(self.expected_bl))
        output = 'output_7_1.im'
        ImageSubtractionMethods.execute(output, image_stack)
        self.assertTrue(os.path.exists(output))

    def test_7_2(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        image_stack.push(EraseableFolder(self.expected_bl))
        output = 'output_7_2.im'
        ImageSubtractionMethods.execute(output, image_stack)
        self.assertTrue(os.path.exists(output))

    @test_base.exception_case(ValueError, 'operands could not be broadcast together with shapes')
    def test_7_3(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.testdata_01[0]))
        image_stack.push(EraseableFolder(self.testdata_02[0]))
        image_stack.push(EraseableFolder(self.testdata_err[0]))
        output = 'output_7_3.im'
        ImageSubtractionMethods.execute(output, image_stack)

    def test_7_4(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.testdata_01[0]))
        image_stack.push(EraseableFolder(self.testdata_err[0]))
        output = 'output_7_4.im'
        ImageSubtractionMethods.execute(output, image_stack)
        self.assertTrue(os.path.exists(output))
        self.assertFalse(os.path.exists(self.testdata_err[0]))

    def test_7_5(self):
        # output = testdata_01 - (testdata_02 - testdata_03)
        image_stack = CasaImageStack(top=UnerasableFolder(self.testdata_01[0]))
        image_stack.push(EraseableFolder(self.testdata_02[0]))
        image_stack.push(EraseableFolder(self.testdata_03[0]))
        output = 'output_7_5.im'
        ImageSubtractionMethods.execute(output, image_stack)
        with tool_manager(output, image) as ia:
            arr = ia.getchunk()
            self.assertTrue(np.array_equal(arr, np.full((64, 64, 4, 128), 2.0)))

    def test_7_6(self):
        # output = testdata_02
        image_stack = CasaImageStack(top=UnerasableFolder(self.testdata_01[0]))
        image_stack.push(EraseableFolder(self.testdata_02[0]))
        output = 'output_7_6.im'
        ImageSubtractionMethods.execute(output, image_stack)
        with tool_manager(output, image) as ia:
            arr = ia.getchunk()
            self.assertTrue(np.array_equal(arr, np.full((64, 64, 4, 128), 2.0)))


class MS2Image_test(test_base):
    """Test MS2Image.

    8-1. successful test
    8-2. base image error
    8-3. MS error
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = 'expected.im'
    expected_orig_im = 'expected_orig.im'
    expected_ms = 'expected.ms'
    expected_bl_ms = 'expected.bl.ms'

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_ms)
        self._copy_test_files(self.datapath, self.expected_bl_ms)
        if os.path.exists(self.expected_im):
            os.rename(self.expected_im, self.expected_orig_im)
        else:
            raise RuntimeError('some errors occured in copying files')
        if not os.path.exists(self.expected_orig_im):
            raise RuntimeError('some errors occured in copying files')
        self.image_shape = get_image_shape(self.expected_orig_im)

    def test_8_1(self):
        MS2ImageMethods.convert(base_image=self.expected_orig_im,
                                input_ms=self.expected_ms,
                                input_image_shape=self.image_shape,
                                datacolumn=DATACOLUMN)
        self.assertTrue(os.path.exists(self.expected_im))
        with tool_manager(self.expected_im, image) as ia:
            arr1 = ia.getchunk()
        with tool_manager(self.expected_orig_im, image) as ia:
            arr2 = ia.getchunk()
        self.assertTrue(np.array_equal(arr1, arr2))

    @test_base.exception_case(TypeError, 'stat: path should be string, ')
    def test_8_2(self):
        MS2ImageMethods.convert(base_image=None,
                                input_ms=self.expected_ms,
                                input_image_shape=self.image_shape,
                                datacolumn=DATACOLUMN)

    @test_base.exception_case(TypeError, 'stat: path should be string, ')
    def test_8_3(self):
        MS2ImageMethods.convert(base_image=self.expected_orig_im,
                                input_ms=self.expected_bl_ms,
                                input_image_shape=self.image_shape,
                                datacolumn=DATACOLUMN)
        converted = 'expected.bl.im'
        self.assertTrue(os.path.exists(converted))
        with tool_manager(self.expected_orig_im, image) as ia:
            arr1 = ia.getchunk()
        with tool_manager(converted, image) as ia:
            arr2 = ia.getchunk()
        self.assertFalse(np.array_equal(arr1, arr2))

        MS2ImageMethods.convert(base_image=self.expected_orig_im,
                                input_ms=None,
                                input_image_shape=self.image_shape,
                                datacolumn=DATACOLUMN)


class misc_test(test_base):
    """Test miscellaneous methods.

    9-1. get_image_shape: successful
    9-2. get_image_shape: failure
    9-3. get_image_shape: failure
    """
    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = 'expected.im'
    g192_im = 'g192_a2.image'

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.g192_im)

    def test_9_1(self):
        shape = get_image_shape(self.expected_im)
        self.assertTrue(np.array_equal(shape.im_shape, [20, 20, 100]))
        self.assertTrue(np.array_equal(shape.axis_dir, [0, 1]))
        self.assertEqual(shape.axis_sp, 2)
        self.assertEqual(shape.axis_pol, -1)
        self.assertTrue(np.array_equal(shape.dir_shape, [20, 20]))
        self.assertEqual(shape.im_nrow, 400)
        self.assertEqual(shape.im_nchan, 100)
        self.assertEqual(shape.im_npol, 1)

        shape = get_image_shape(self.g192_im)
        self.assertTrue(np.array_equal(shape.im_shape, [512, 512, 1, 40]))
        self.assertTrue(np.array_equal(shape.axis_dir, [0, 1]))
        self.assertEqual(shape.axis_sp, 3)
        self.assertEqual(shape.axis_pol, 2)
        self.assertTrue(np.array_equal(shape.dir_shape, [512, 512]))
        self.assertEqual(shape.im_nrow, 262144)
        self.assertEqual(shape.im_nchan, 40)
        self.assertEqual(shape.im_npol, 1)

    @test_base.exception_case(ValueError, 'path \'notexists\' is not found')
    def test_9_2(self):
        get_image_shape('notexists')

    @test_base.exception_case(ValueError, 'image \'testdata_01.im\' is invalid')
    def test_9_3(self):
        testimage = 'testdata_01.im'
        self._create_image(testimage, 1.0, [64, 64])
        get_image_shape(testimage)


class imbaseline_test(test_base):
    """Test full of imbaseline.

    F-1. simple successful case
    F-2. imagefile is None
    F-3. output_cont is False
    F-4. bloutput is specified
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected = 'expected.im'

    def setUp(self):
        self.ia = image()
        self._copy_test_files(self.datapath, self.expected)

    def test_f_1(self):
        imagefile = self.expected
        linefile = 'output_f_1'
        dirkernel = 'gaussian'
        spkernel = 'gaussian'
        major = '20arcsec'
        minor = '10arcsec'
        pa = '0deg'
        blfunc = 'sinusoid'
        output_cont = True

        imbaseline(imagename=imagefile,
                   linefile=linefile,
                   dirkernel=dirkernel,
                   spkernel=spkernel,
                   major=major,
                   minor=minor,
                   pa=pa,
                   blfunc=blfunc,
                   output_cont=output_cont)
        self.assertTrue(os.path.exists(linefile))

    @test_base.exception_case(ValueError, 'file  is not found.')
    def test_f_2(self):
        imagefile = ''
        linefile = 'output_f_2'
        dirkernel = 'gaussian'
        spkernel = 'gaussian'
        major = '20arcsec'
        minor = '10arcsec'
        pa = '0deg'
        blfunc = 'sinusoid'
        output_cont = True

        imbaseline(imagename=imagefile,
                   linefile=linefile,
                   dirkernel=dirkernel,
                   spkernel=spkernel,
                   major=major,
                   minor=minor,
                   pa=pa,
                   blfunc=blfunc,
                   output_cont=output_cont)

    def test_f_3(self):
        imagefile = self.expected
        linefile = 'output_f_3'
        dirkernel = 'gaussian'
        spkernel = 'gaussian'
        major = '20arcsec'
        minor = '10arcsec'
        pa = '0deg'
        blfunc = 'sinusoid'
        output_cont = False

        imbaseline(imagename=imagefile,
                   linefile=linefile,
                   dirkernel=dirkernel,
                   spkernel=spkernel,
                   major=major,
                   minor=minor,
                   pa=pa,
                   blfunc=blfunc,
                   output_cont=output_cont)
        self.assertFalse(os.path.exists(linefile + '.cont'))

    def test_f_4(self):
        imagefile = self.expected
        linefile = 'output_f_4'
        dirkernel = 'gaussian'
        spkernel = 'gaussian'
        major = '20arcsec'
        minor = '10arcsec'
        pa = '0deg'
        blfunc = 'sinusoid'
        output_cont = True
        bloutput = self.expected + 'bloutput'

        imbaseline(imagename=imagefile,
                   linefile=linefile,
                   dirkernel=dirkernel,
                   spkernel=spkernel,
                   major=major,
                   minor=minor,
                   pa=pa,
                   blfunc=blfunc,
                   output_cont=output_cont,
                   bloutput=bloutput)
        self.assertTrue(os.path.exists(bloutput))


def suite():
    return [imsmooth_test, AbstractFileStack_test, ImageShape_test, imbaseline_test, image2ms_test, sdbaseline_test,
            image_subtraction_test, misc_test]

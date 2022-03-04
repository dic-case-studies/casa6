"""
Imbaseline test class.

This module contains tests for the classes of imbaseline. Since the main feature of the task is
made up of a combination of imsmooth/sdsmooth/sdbaseline, so tests consist of simple tests for
these tasks and white tests for the unique classes of imbaseline such as image subtraction.
"""

import functools
import os
import re
import shutil
import unittest
import uuid

import numpy as np

from casatasks import casalog
from casatasks.private.sdutil import tool_manager
from casatasks.private.task_imbaseline import (_CasaImageStack,
                                               _eraseable_folder_register,
                                               _EraseableFolder,
                                               _get_image_shape,
                                               _Image2MSMethods,
                                               _Image2MSParams, _ImageShape,
                                               _ImageSubtractionMethods,
                                               _ImsmoothMethods,
                                               _ImsmoothParams,
                                               _MeasurementSetStack,
                                               _MS2ImageMethods,
                                               _SdbaselineMethods,
                                               _SdbaselineParams,
                                               _SdsmoothMethods,
                                               _SdsmoothParams,
                                               _UnerasableFolder, imbaseline)
from casatools import ctsys, image, table

_tb = table()

# https://open-bitbucket.nrao.edu/projects/CASA/repos/casatestdata/browse/unittest
DATAPATH = ctsys.resolve('unittest/imbaseline/')

DATACOLUMN = 'DATA'
UNEXISTS = 'unexists'
DUMMY_FOLDERS = ('dummy1', 'dummy2', 'dummy3')

casalog.origin('imbaseline')


class test_base(unittest.TestCase):
    """Base class of ibmaseline testing."""

    @staticmethod
    def exception_case(exception_type, exception_pattern):
        """Decorate tests intended to throw a specific exception.

        exception_type: type of exception
        exception_pattern: regex for inspecting exception message using re.search
        """
        def wrapper(func):
            @functools.wraps(func)
            def _wrapper(self):
                self.assertTrue(len(exception_pattern) > 0,
                                msg='Internal Error')
                with self.assertRaises(exception_type) as ctx:
                    func(self)
                    self.fail(msg='The task must throw an exception')
                the_exception = ctx.exception
                message = str(the_exception)
                self.assertIsNotNone(
                    re.search(
                        exception_pattern,
                        message),
                    msg=f'Expected: \'{exception_pattern}\', got: \'{message}\'')
            return _wrapper
        return wrapper

    @classmethod
    def setUpClass(cls):
        cls.workdir = Workdir.create()
        cls.workdir.chdir()

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        cls.workdir.clean()

    def tearDown(self):
        _eraseable_folder_register.clear()
        self.assertTrue(len(_tb.showcache()) == 0)
        # make sure directory is clean as per verification test requirement
        cwd = os.getcwd()
        for filename in os.listdir(cwd):
            file_path = os.path.join(cwd, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

    def _create_dummy_folders(self):
        def _setup_folder(folder):
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.mkdir(folder)

        [_setup_folder(folder) for folder in DUMMY_FOLDERS]

    def _copy_test_files(self, filename):
        """Copy files for testing into current path."""
        src = os.path.join(DATAPATH, filename)
        dst = os.path.join(os.getcwd(), filename)

        if os.path.exists(dst):
            if os.path.isfile(dst) or os.path.islink(dst):
                os.unlink(dst)
            elif os.path.isdir(dst):
                shutil.rmtree(dst)

        if os.path.isfile(src):
            shutil.copy(src, dst)
        elif os.path.isdir(src):
            shutil.copytree(src, dst, symlinks=False)

    def _create_image(self, datapath, val=1, shape=[0, 0, 0, 0]):
        _ia = image()
        ary = _ia.makearray(v=val, shape=shape)
        _ia.fromarray(outfile=datapath, pixels=ary, overwrite=True)
        _ia.done()

    def _check_ms_tables(self, path):
        self.assertTrue(os.path.exists(path))
        for table_name in (
            '',
            'ANTENNA',
            'DATA_DESCRIPTION',
            'FEED',
            'FIELD',
            'FLAG_CMD',
            'HISTORY',
            'OBSERVATION',
            'POINTING',
            'POLARIZATION',
            'PROCESSOR',
            'SOURCE',
            'SPECTRAL_WINDOW',
                'STATE'):
            table_path = os.path.join(path, table_name)
            self.assertTrue(os.path.exists(
                os.path.join(table_path, 'table.dat')))


class TestAbstractFileStack(test_base):
    """Test AbstractFileStack / (Un)_EraseableFolder.

    1-1. successful case: Create Stack with exist file
    1-2. failure case: Create Stack with unexist file, an exception raises
    1-3. successful case: push() exist file into stack
    1-4. failure case: push() unexist file into stack, an exception raises
    1-5. successful case: pop() exist stuff into stack
    1-6. failure case: pop() unexist stuff into stack, an exception raises
    1-7. successful case: peak() exist stuff into stack
    1-8. failure case: peak() unexist stuff into stack, an exception raises
    1-9. successful case: subpeak() exist stuff into stack
    1-10. failure case: subpeak() unexist stuff into stack, an exception raises
    1-11. successful case: bottom() exist stuff into stack
    1-12. failure case: bottom() unexist stuff into stack, an exception raises
    1-13. successful case: do clear() of _EraseableFolder contains a file
    1-14. successful case: do clear() of Un_EraseableFolder contains a file
    1-15. successful case: do erase() of _EraseableFolder contains a file
    1-16. successful case: do erase() Un_EraseableFolder contains a file
    1-17. successful case: check height() of stack
    1-18. failure case: pop() when height is 1, an exception raises
    """

    def setUp(self):
        self._create_dummy_folders()
        if os.path.exists(UNEXISTS):
            shutil.rmtree(UNEXISTS)

    def test_1_1(self):
        """1-1. successful case: Create Stack with exist file."""
        stack = _CasaImageStack(_UnerasableFolder(DUMMY_FOLDERS[0]))
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(ValueError, f'file {UNEXISTS} is not found')
    def test_1_2(self):
        """1-2. failure case: Create Stack with unexist file, an exception raises."""
        _CasaImageStack(_UnerasableFolder(UNEXISTS))

    def test_1_3(self):
        """1-3. successful case: push() exist file into stack."""
        stack = _CasaImageStack()
        stack.push(_UnerasableFolder(DUMMY_FOLDERS[0]))
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(ValueError, f'file {UNEXISTS} is not found')
    def test_1_4(self):
        """1-4. failure case: push() unexist file into stack, an exception raises."""
        stack = _CasaImageStack()
        stack.push(_UnerasableFolder(UNEXISTS))

    def test_1_5(self):
        """1-5. successful case: pop() exist stuff into stack."""
        stack = _CasaImageStack()
        stack.push(_UnerasableFolder(DUMMY_FOLDERS[0]))
        obj = _UnerasableFolder(DUMMY_FOLDERS[1])
        stack.push(obj)
        tmp = stack.pop()
        self.assertEqual(obj, tmp)
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(RuntimeError, 'the stack cannot pop')
    def test_1_6(self):
        """1-6. failure case: pop() unexist stuff into stack, an exception raises."""
        stack = _CasaImageStack()
        stack.pop()

    def test_1_7(self):
        """1-7. successful case: peak() exist stuff into stack."""
        stack = _CasaImageStack()
        obj1 = _UnerasableFolder(DUMMY_FOLDERS[0])
        stack.push(obj1)
        obj2 = _UnerasableFolder(DUMMY_FOLDERS[1])
        stack.push(obj2)
        obj3 = _UnerasableFolder(DUMMY_FOLDERS[2])
        stack.push(obj3)
        self.assertEqual(stack.peak(), obj3)
        self.assertEqual(stack.subpeak(), obj2)
        self.assertEqual(stack.bottom(), obj1)

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_1_8(self):
        """1-8. failure case: peak() unexist stuff into stack, an exception raises."""
        stack = _CasaImageStack()
        stack.peak()

    def test_1_9(self):
        """1-9. successful case: subpeak() exist stuff into stack."""
        stack = _CasaImageStack()
        obj1 = _UnerasableFolder(DUMMY_FOLDERS[0])
        stack.push(obj1)
        obj2 = _UnerasableFolder(DUMMY_FOLDERS[1])
        stack.push(obj2)
        self.assertEqual(stack.subpeak(), obj1)
        self.assertEqual(stack.bottom(), obj1)
        obj3 = _UnerasableFolder(DUMMY_FOLDERS[2])
        stack.push(obj3)
        self.assertEqual(stack.subpeak(), obj2)
        self.assertEqual(stack.bottom(), obj1)

    @test_base.exception_case(RuntimeError, 'the stack has only one stuff')
    def test_1_10(self):
        """1-10. failure case: subpeak() unexist stuff into stack, an exception raises."""
        stack = _CasaImageStack(_UnerasableFolder(DUMMY_FOLDERS[0]))
        stack.subpeak()

    def test_1_11(self):
        """1-11. successful case: bottom() exist stuff into stack."""
        stack = _CasaImageStack()
        obj1 = _UnerasableFolder(DUMMY_FOLDERS[0])
        stack.push(obj1)
        self.assertEqual(stack.bottom(), obj1)
        obj2 = _UnerasableFolder(DUMMY_FOLDERS[1])
        stack.push(obj2)
        self.assertEqual(stack.bottom(), obj1)
        self.assertEqual(stack.peak(), obj2)
        obj3 = _UnerasableFolder(DUMMY_FOLDERS[2])
        stack.push(obj3)
        self.assertEqual(stack.bottom(), obj1)
        self.assertEqual(stack.peak(), obj3)
        self.assertEqual(stack.subpeak(), obj2)
        stack.pop()
        self.assertEqual(stack.peak(), obj2)
        self.assertEqual(stack.bottom(), obj1)

    @test_base.exception_case(RuntimeError,
                              'the stack has not have enough stuff')
    def test_1_12(self):
        """1-12. failure case: bottom() unexist stuff into stack, an exception raises."""
        stack = _CasaImageStack()
        stack.bottom()

    def test_1_13(self):
        """1-13. successful case: do clear() of _EraseableFolder contains a file."""
        file = _EraseableFolder(DUMMY_FOLDERS[0])
        stack = _CasaImageStack(file)
        stack.clear(False)
        self.assertTrue(os.path.exists(DUMMY_FOLDERS[0]))
        self.assertEqual(stack.height(), 0)

    def test_1_14(self):
        """1-14. successful case: do clear() of Un_EraseableFolder contains a file."""
        stack = _CasaImageStack(_UnerasableFolder(DUMMY_FOLDERS[0]))
        stack.clear(False)
        self.assertTrue(os.path.exists(DUMMY_FOLDERS[0]))
        self.assertEqual(stack.height(), 0)

    def test_1_15(self):
        """1-15. successful case: do erase() of _EraseableFolder contains a file."""
        file = _EraseableFolder(DUMMY_FOLDERS[0])
        file.erase(False)
        self.assertFalse(os.path.exists(DUMMY_FOLDERS[0]))

    def test_1_16(self):
        """1-16. successful case: do erase() Un_EraseableFolder contains a file."""
        file = _UnerasableFolder(DUMMY_FOLDERS[0])
        file.erase(False)
        self.assertTrue(os.path.exists(DUMMY_FOLDERS[0]))

    def test_1_17(self):
        """1-17. successful case: check height() of stack."""
        stack = _CasaImageStack()
        self.assertEqual(stack.height(), 0)
        stack.push(_UnerasableFolder(DUMMY_FOLDERS[0]))
        self.assertEqual(stack.height(), 1)
        stack.push(_UnerasableFolder(DUMMY_FOLDERS[1]))
        self.assertEqual(stack.height(), 2)
        stack.push(_UnerasableFolder(DUMMY_FOLDERS[2]))
        self.assertEqual(stack.height(), 3)
        stack.pop()
        self.assertEqual(stack.height(), 2)
        stack.pop()
        self.assertEqual(stack.height(), 1)

    @test_base.exception_case(RuntimeError, 'the stack cannot pop')
    def test_1_18(self):
        """1-18. failure case: pop() when height is 1, an exception raises."""
        stack = _CasaImageStack()
        stack.push(_UnerasableFolder(DUMMY_FOLDERS[0]))
        stack.pop()


class Test_ImageShape(test_base):
    """Test _ImageShape.

    2-1. successful case: axis_sp = axis_pol = 2or3, axis_sp != axis_pol
    2-2. failure case: invalid im_nchan, an exception raises
    2-3. failure case: invalid dir_shape, an exception raises
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_2_1(self):
        """2-1. successful case: axis_sp = axis_pol = 2or3, axis_sp != axis_pol."""
        shape = _ImageShape(im_shape=np.array(
            [100, 100, 1, 100]), axis_dir=np.array([0, 1]), axis_sp=3, axis_pol=2)
        shape.validate()
        shape = _ImageShape(im_shape=np.array(
            [100, 100, 100, 1]), axis_dir=np.array([0, 1]), axis_sp=2, axis_pol=3)
        shape.validate()
        # any exceptions are not thrown, its OK

    @test_base.exception_case(ValueError,
                              'nchan \\d is too few to perform baseline subtraction')
    def test_2_2(self):
        """2-2. failure case: invalid im_nchan, an exception raises."""
        shape = _ImageShape(im_shape=np.array(
            [100, 100, 1, 1]), axis_dir=np.array([0, 1]), axis_sp=3, axis_pol=2)
        shape.validate()

    @test_base.exception_case(ValueError,
                              'invalid value: dir_shape \\[\\d+\\]')
    def test_2_3(self):
        """2-3. failure case: invalid dir_shape, an exception raises."""
        shape = _ImageShape(np.array([100, 100, 1, 100]), axis_dir=np.array(
            [0]), axis_sp=3, axis_pol=2)
        shape.validate()


class TestImsmooth(test_base):
    """Test imsmooth execution.

    Tests of imsmooth rely on ones of test_imsmooth basically, so we have minimal tests in imbaseline.

    3-1. successful case: call imsmooth with some parameters
    3-2. failure case: invalid dirkernel, an exception raises
    3-3. set values for _ImsmoothParams and do validate(), and compare properties of it to the correct values
    """

    tiny = 'tiny.im'

    def setUp(self):
        self._copy_test_files(self.tiny)

    def test_3_1(self):
        """3-1. successful case: call imsmooth with some parameters."""
        major = '2.5arcsec'
        minor = '2arcsec'
        pa = '0deg'
        dirkernel = 'gaussian'
        kimage = ''
        scale = -1

        stack = _CasaImageStack(top=_UnerasableFolder(self.tiny))

        _ImsmoothMethods.execute(
            dirkernel, major, minor, pa, kimage, scale, stack)
        self.assertTrue(os.path.exists(stack.peak().path))

    @test_base.exception_case(ValueError,
                              'Unsupported direction smoothing kernel, foobar')
    def test_3_2(self):
        """3-2. failure case: invalid dirkernel, an exception raises."""
        major = '2.5arcsec'
        minor = '2arcsec'
        pa = '0deg'
        dirkernel = 'foobar'
        kimage = ''
        scale = -1

        stack = _CasaImageStack(top=_UnerasableFolder(self.tiny))

        _ImsmoothMethods.execute(
            dirkernel, major, minor, pa, kimage, scale, stack)

    def test_3_3(self):
        """3-3. set values for _ImsmoothParams and do validate(), and compare properties of it to the correct values."""
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

        def compare_params(
                _kernel,
                _major=major,
                _minor=minor,
                _pa=pa,
                _kimage=kimage,
                _scale=scale):
            valid_param = dict(
                targetres=targetres,
                mask=mask,
                beam=beam,
                region=region,
                box=box,
                chans=chans,
                stokes=stokes,
                stretch=stretch,
                overwrite=True,
                imagename=infile,
                outfile=outfile,
                kernel=_kernel,
                major=_major,
                minor=_minor,
                pa=_pa,
                kimage=_kimage,
                scale=_scale,
                __log_origin=logorigin)
            param = _ImsmoothParams(
                infile, outfile, _kernel, _major, _minor, _pa, _kimage, _scale)
            param.validate()
            self.assertEqual(param(), valid_param)

        # none
        compare_params(kernel[0])

        # image
        compare_params(kernel[1], _major='', _minor='', _pa='')

        # gaussian
        compare_params(kernel[2], _kimage='', _scale=-1.0)

        # boxcar
        compare_params(kernel[3], _kimage='', _scale=-1.0)


class TestImage2MS(test_base):
    """Test image2ms.

    4-1. successful case: create MeasurementSet from a image
    4-2. failure case: execute image2ms with invalid datacolumn parameter, an exception raises
    4-3. failure case: execute image2ms with invalid image data, an exception raises
    4-4. failure case: execute image2ms with empty stack, an exception raises
    4-5. set values for _Image2MSParams and do validate(), and compare properties of it to the correct values
    """

    expected = 'expected.im'
    datacolumn = DATACOLUMN

    def setUp(self):
        self._create_dummy_folders()
        self._copy_test_files(self.expected)
        self.image_shape = _get_image_shape(self.expected)

    def test_4_1(self):
        """4-1. successful case: create MeasurementSet from a image."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(self.expected))
        ms_stack = _MeasurementSetStack()
        _Image2MSMethods.execute(
            self.datacolumn, self.image_shape, image_stack, ms_stack)
        self.assertEqual(ms_stack.height(), 1)
        self._check_ms_tables(ms_stack.peak().path)

    @test_base.exception_case(RuntimeError, 'column INVALID does not exist')
    def test_4_2(self):
        """4-2. failure case: execute image2ms with invalid datacolumn parameter, an exception raises."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(self.expected))
        ms_stack = _MeasurementSetStack()
        _Image2MSMethods.execute(
            'INVALID', self.image_shape, image_stack, ms_stack)

    @test_base.exception_case(RuntimeError, 'Unable to open image dummy1.')
    def test_4_3(self):
        """4-3. failure case: execute image2ms with invalid image data, an exception raises."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(DUMMY_FOLDERS[0]))
        ms_stack = _MeasurementSetStack()
        _Image2MSMethods.execute(
            self.datacolumn, self.image_shape, image_stack, ms_stack)

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_4_4(self):
        """4-4. failure case: execute image2ms with empty stack, an exception raises."""
        image_stack = _CasaImageStack()
        ms_stack = _MeasurementSetStack()
        _Image2MSMethods.execute(
            self.datacolumn, self.image_shape, image_stack, ms_stack)

    def test_4_5(self):
        """4-5. set values for _Image2MSParams and do validate(), and compare properties of it to the correct values."""
        outfile = 'output_4_5.ms'
        params = _Image2MSParams(
            self.expected, outfile, self.datacolumn, self.image_shape)
        params.validate()
        self.assertEqual(params.infile, self.expected)
        self.assertEqual(params.outfile, outfile)
        for attr in ('im_shape', 'axis_dir', 'dir_shape'):
            self.assertTrue(np.all(getattr(params, attr) == getattr(self.image_shape, attr)))
        for attr in ('axis_sp', 'axis_pol', 'im_nrow', 'im_nchan', 'im_npol'):
            self.assertEqual(getattr(params, attr),
                             getattr(self.image_shape, attr))


class TestSdsmooth(test_base):
    """Test sdsmooth execution.

    Tests of sdsmooth rely on ones of test_sdsmooth basically, so we have minimal tests in imbaseline.

    5-1. successful case: call sdsmooth with some parameters
    5-2. failure case: call sdsmooth with invalid ms stack, an exception raises
    5-3. failure case: call sdsmooth with invalid image stack, an exception raises
    5-4. set values for _SdsmoothParams and do validate(), and compare properties of it to the correct values
    """

    input_image = 'expected.im'
    input_ms = 'expected.ms'
    datacolumn = DATACOLUMN
    spkenel = 'gaussian'
    kwidth = 5

    def setUp(self):
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.input_ms)
        self.image_shape = _get_image_shape(self.input_image)

    def test_5_1(self):
        """5-1. successful case: call sdsmooth with some parameters."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(self.input_image))
        ms_stack = _MeasurementSetStack()
        ms_stack.push(_EraseableFolder(self.input_ms))
        _SdsmoothMethods.execute(
            self.datacolumn,
            self.spkenel,
            self.kwidth,
            image_stack,
            ms_stack,
            self.image_shape)
        self.assertEqual(image_stack.height(), 2)
        self.assertEqual(ms_stack.height(), 2)
        self.assertTrue(os.path.exists(os.path.join(
            image_stack.peak().path, 'table.dat')))
        self._check_ms_tables(ms_stack.peak().path)

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_5_2(self):
        """5-2. failure case: call sdsmooth with invalid ms stack, an exception raises."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(self.input_image))
        ms_stack = _MeasurementSetStack()
        _SdsmoothMethods.execute(
            self.datacolumn,
            self.spkenel,
            self.kwidth,
            image_stack,
            ms_stack,
            self.image_shape)

    @test_base.exception_case(RuntimeError,
                              'the stack has not have enough stuff')
    def test_5_3(self):
        """5-3. failure case: call sdsmooth with invalid image stack, an exception raises."""
        image_stack = _CasaImageStack()
        ms_stack = _MeasurementSetStack()
        ms_stack.push(_EraseableFolder(self.input_ms))
        _SdsmoothMethods.execute(
            self.datacolumn,
            self.spkenel,
            self.kwidth,
            image_stack,
            ms_stack,
            self.image_shape)

    def test_5_4(self):
        """5-4. set values for _SdsmoothParams and do validate(), and compare properties of it to the correct values."""
        spw = field = antenna = timerange = scan = pol = intent = ''
        reindex = overwrite = True
        infile = 'infile'
        outfile = 'outfile'
        datacolumn = DATACOLUMN
        kernel = ('none', 'gaussian', 'boxcar')
        kwidth = 5
        logorigin = 'imbaseline'

        def compare_params(_kernel):
            valid_params = dict(
                spw=spw,
                field=field,
                antenna=antenna,
                timerange=timerange,
                scan=scan,
                pol=pol,
                intent=intent,
                reindex=reindex,
                overwrite=overwrite,
                infile=infile,
                datacolumn=datacolumn,
                kernel=_kernel,
                kwidth=kwidth,
                outfile=outfile,
                __log_origin=logorigin)
            params = _SdsmoothParams(
                infile=infile,
                outfile=outfile,
                datacolumn=datacolumn,
                spkernel=_kernel,
                kwidth=kwidth)
            params.validate()
            self.assertEqual(params(), valid_params)

        [compare_params(_kernel) for _kernel in kernel]


class TestSdbaseline(test_base):
    """Test sdbaseline execution.

    Tests of sdbaseline rely on ones of test_sdbaseline basically, so we have minimal tests in imbaseline.

    6-1. successful case: call sdbaseline with some parameters
    6-2. failure case: call sdbaseline with invalid ms stack, an exception raises
    6-3. failure case: call sdbaseline with invalid image stack, an exception raise
    6-4. set values for _SdbaselineParams and do validate(), and compare properties of it to the correct values
    """

    input_image = 'expected.im'
    input_ms = 'expected.ms'
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
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.input_ms)
        self._copy_test_files(self.blparam)
        self.image_shape = _get_image_shape(self.input_image)

    def test_6_1(self):
        """6-1. successful case: call sdbaseline with some parameters."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(self.input_image))
        ms_stack = _MeasurementSetStack()
        ms_stack.push(_EraseableFolder(self.input_ms))
        _SdbaselineMethods.execute(
            self.datacolumn,
            self.bloutput,
            self.maskmode,
            self.chans,
            self.thresh,
            self.avg_limit,
            self.minwidth,
            self.edge,
            self.blfunc,
            self.order,
            self.npiece,
            self.applyfft,
            self.fftthresh,
            self.addwn,
            self.rejwn,
            self.blparam,
            self.clipniter,
            self.clipthresh,
            image_stack,
            ms_stack,
            self.image_shape)
        self.assertTrue(os.path.exists(ms_stack.peak().path))
        self.assertTrue(os.path.exists(self.bloutput))
        self.assertTrue(os.path.exists(image_stack.peak().path))

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_6_2(self):
        """6-2. failure case: call sdbaseline with invalid ms stack, an exception raises."""
        image_stack = _CasaImageStack(top=_UnerasableFolder(self.input_image))
        ms_stack = _MeasurementSetStack()
        _SdbaselineMethods.execute(
            self.datacolumn,
            self.bloutput,
            self.maskmode,
            self.chans,
            self.thresh,
            self.avg_limit,
            self.minwidth,
            self.edge,
            self.blfunc,
            self.order,
            self.npiece,
            self.applyfft,
            self.fftthresh,
            self.addwn,
            self.rejwn,
            self.blparam,
            self.clipniter,
            self.clipthresh,
            image_stack,
            ms_stack,
            self.image_shape)

    @test_base.exception_case(RuntimeError,
                              'the stack has not have enough stuff')
    def test_6_3(self):
        """6-3. failure case: call sdbaseline with invalid image stack, an exception raise."""
        image_stack = _CasaImageStack()
        ms_stack = _MeasurementSetStack()
        ms_stack.push(_EraseableFolder(self.input_ms))
        _SdbaselineMethods.execute(
            self.datacolumn,
            self.bloutput,
            self.maskmode,
            self.chans,
            self.thresh,
            self.avg_limit,
            self.minwidth,
            self.edge,
            self.blfunc,
            self.order,
            self.npiece,
            self.applyfft,
            self.fftthresh,
            self.addwn,
            self.rejwn,
            self.blparam,
            self.clipniter,
            self.clipthresh,
            image_stack,
            ms_stack,
            self.image_shape)

    def test_6_4(self):
        """6-4. set values for _SdbaselineParams and do validate(), and compare properties of it to the correct values."""
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
        spw = '0'
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
            valid_param = dict(
                antenna=antenna,
                field=field,
                spw=spw,
                timerange=timerange,
                scan=scan,
                pol=pol,
                intent=intent,
                reindex=reindex,
                blmode=blmode,
                dosubtract=dosubtract,
                blformat=blformat,
                bltable=bltable,
                updateweight=updateweight,
                sigmavalue=sigmavalue,
                showprogress=showprogress,
                minnrow=minnrow,
                fftmethod=fftmethod,
                verbose=verbose,
                overwrite=overwrite,
                infile=infile,
                datacolumn=datacolumn,
                maskmode=_maskmode,
                thresh=thresh,
                avg_limit=avg_limit,
                minwidth=minwidth,
                edge=edge,
                bloutput=bloutput,
                blfunc=_blfunc,
                order=order,
                npiece=npiece,
                applyfft=applyfft,
                fftthresh=fftthresh,
                addwn=addwn,
                rejwn=rejwn,
                clipthresh=clipthresh,
                clipniter=clipniter,
                blparam=blparam,
                outfile=outfile,
                __log_origin=logorigin)
            params = _SdbaselineParams(
                infile=infile,
                outfile=outfile,
                datacolumn=datacolumn,
                bloutput=bloutput,
                maskmode=_maskmode,
                chans=chans,
                thresh=thresh,
                avg_limit=avg_limit,
                minwidth=minwidth,
                edge=edge,
                blfunc=_blfunc,
                order=order,
                npiece=npiece,
                applyfft=applyfft,
                fftthresh=fftthresh,
                addwn=addwn,
                rejwn=rejwn,
                blparam=blparam,
                clipniter=clipniter,
                clipthresh=clipthresh)
            params.validate()
            self.assertEqual(params(), valid_param)

        [compare_params(_maskmode, _blfunc)
         for _maskmode in maskmode for _blfunc in blfunc]


class TestImageSubtraction(test_base):
    """Test image subtractions.

    7-1. successful test: subtracted output = input_image - (smoothed_image - smoothed_and_subtracted_image)
    7-2. successful test: subtracted output = subtracted_image
    7-3. failure case: subtract three images have unmatched shape, an exception raises
    7-4. successful test: subtract two images have unmatched shape (any exceptions do not raise)
    7-5. output data check: three images subtraction test
    7-6. output data check: two images subtraction test
    """

    existing_image = 'expected.im'
    existing_imsmoothed_image = 'expected.imsmooth.im'
    existing_baselined_image = 'expected.bl.im'
    input_image = ('input_image.im', 1.5, [64, 64, 4, 128])
    smoothed_image = ('smoothed_image.im', 2.0, [64, 64, 4, 128])
    smoothed_and_subtracted_image = (
        'smoothed_and_subtracted_image.im', 2.5, [64, 64, 4, 128])
    testdata_err = ('testdata_err.im', 1, [65, 64, 4, 128])

    def setUp(self):
        self._copy_test_files(self.existing_image)
        self._copy_test_files(self.existing_imsmoothed_image)
        self._copy_test_files(self.existing_baselined_image)
        self._create_image(
            self.input_image[0], self.input_image[1], self.input_image[2])
        self._create_image(
            self.smoothed_image[0],
            self.smoothed_image[1],
            self.smoothed_image[2])
        self._create_image(
            self.smoothed_and_subtracted_image[0],
            self.smoothed_and_subtracted_image[1],
            self.smoothed_and_subtracted_image[2])
        self._create_image(
            self.testdata_err[0], self.testdata_err[1], self.testdata_err[2])

    def test_7_1(self):
        """7-1. successful test: subtracted output = input_image - (smoothed_image - smoothed_and_subtracted_image)."""
        image_stack = _CasaImageStack(
            top=_UnerasableFolder(self.existing_image))
        image_stack.push(_EraseableFolder(self.existing_imsmoothed_image))
        image_stack.push(_EraseableFolder(self.existing_baselined_image))
        output = 'output_7_1.im'
        _ImageSubtractionMethods.execute(output, image_stack)
        self.assertTrue(os.path.exists(output))

    def test_7_2(self):
        """7-2. successful test: subtracted output = subtracted_image."""
        image_stack = _CasaImageStack(
            top=_UnerasableFolder(self.existing_image))
        image_stack.push(_EraseableFolder(self.existing_baselined_image))
        output = 'output_7_2.im'
        _ImageSubtractionMethods.execute(output, image_stack)
        self.assertTrue(os.path.exists(output))

    @test_base.exception_case(ValueError,
                              'operands could not be broadcast together with shapes')
    def test_7_3(self):
        """7-3. failure case: subtract three images have unmatched shape, an exception raises."""
        image_stack = _CasaImageStack(
            top=_UnerasableFolder(self.input_image[0]))
        image_stack.push(_EraseableFolder(self.smoothed_image[0]))
        image_stack.push(_EraseableFolder(self.testdata_err[0]))
        output = 'output_7_3.im'
        _ImageSubtractionMethods.execute(output, image_stack)

    def test_7_4(self):
        """7-4. successful test: subtract two images have unmatched shape (any exceptions do not raise)."""
        image_stack = _CasaImageStack(
            top=_UnerasableFolder(self.input_image[0]))
        image_stack.push(_EraseableFolder(self.testdata_err[0]))
        output = 'output_7_4.im'
        _ImageSubtractionMethods.execute(output, image_stack)
        self.assertTrue(os.path.exists(output))
        self.assertFalse(os.path.exists(self.testdata_err[0]))

    def test_7_5(self):
        """7-5. output data check: three images subtraction test."""
        # output = input_image - (smoothed_image - smoothed_and_subtracted_image)
        image_stack = _CasaImageStack(
            top=_UnerasableFolder(self.input_image[0]))
        image_stack.push(_EraseableFolder(self.smoothed_image[0]))
        image_stack.push(_EraseableFolder(
            self.smoothed_and_subtracted_image[0]))
        output = 'output_7_5.im'
        _ImageSubtractionMethods.execute(output, image_stack)
        with tool_manager(output, image) as ia:
            arr = ia.getchunk()
            self.assertTrue(np.array_equal(
                arr, np.full((64, 64, 4, 128), 2.0)))

    def test_7_6(self):
        """7-6. output data check: two images subtraction test."""
        # output = smoothed_image
        image_stack = _CasaImageStack(
            top=_UnerasableFolder(self.input_image[0]))
        image_stack.push(_EraseableFolder(self.smoothed_image[0]))
        output = 'output_7_6.im'
        _ImageSubtractionMethods.execute(output, image_stack)
        with tool_manager(output, image) as ia:
            arr = ia.getchunk()
            self.assertTrue(np.array_equal(
                arr, np.full((64, 64, 4, 128), 2.0)))


class TestMS2Image(test_base):
    """Test MS2Image.

    8-1. successful case: convert a MeasurementSet to image and compare it to the correct chunk
    8-2. failure case: attempt to convert a MeasurementSet without base image, an exception raises
    8-3. failure case: attempt to convert a MeasurementSet unexisted, an exception raises
    """

    input_image = 'expected.im'
    original_image = 'expected_orig.im'
    input_ms = 'expected.ms'
    baselined_ms = 'expected.bl.ms'

    def setUp(self):
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.input_ms)
        self._copy_test_files(self.baselined_ms)
        if os.path.exists(self.input_image):
            os.rename(self.input_image, self.original_image)
        else:
            raise RuntimeError('some errors occured in copying files')
        if not os.path.exists(self.original_image):
            raise RuntimeError('some errors occured in copying files')
        self.image_shape = _get_image_shape(self.original_image)

    def test_8_1(self):
        """8-1. successful case: convert a MeasurementSet to image and compare it to the correct chunk."""
        _MS2ImageMethods.convert(base_image=self.original_image,
                                 input_ms=self.input_ms,
                                 input_image_shape=self.image_shape,
                                 datacolumn=DATACOLUMN)
        self.assertTrue(os.path.exists(self.input_image))
        with tool_manager(self.input_image, image) as ia:
            arr1 = ia.getchunk()
        with tool_manager(self.original_image, image) as ia:
            arr2 = ia.getchunk()
        self.assertTrue(np.array_equal(arr1, arr2))

    @test_base.exception_case(TypeError, 'stat: path should be string, ')
    def test_8_2(self):
        """8-2. failure case: attempt to convert a MeasurementSet without base image, an exception raises."""
        _MS2ImageMethods.convert(base_image=None,
                                 input_ms=self.input_ms,
                                 input_image_shape=self.image_shape,
                                 datacolumn=DATACOLUMN)

    @test_base.exception_case(TypeError, 'stat: path should be string, ')
    def test_8_3(self):
        """8-3. failure case: attempt to convert a MeasurementSet unexisted, an exception raises."""
        _MS2ImageMethods.convert(base_image=self.original_image,
                                 input_ms=self.baselined_ms,
                                 input_image_shape=self.image_shape,
                                 datacolumn=DATACOLUMN)
        converted = 'expected.bl.im'
        self.assertTrue(os.path.exists(converted))
        with tool_manager(self.original_image, image) as ia:
            arr1 = ia.getchunk()
        with tool_manager(converted, image) as ia:
            arr2 = ia.getchunk()
        self.assertFalse(np.array_equal(arr1, arr2))

        _MS2ImageMethods.convert(base_image=self.original_image,
                                 input_ms=None,
                                 input_image_shape=self.image_shape,
                                 datacolumn=DATACOLUMN)


class TestModuleMethodsOfImbaseline(test_base):
    """Test global methods.

    The white tests of global methods of imbaseline module must be implemented in this class.

    9-1. _get_image_shape: successful case: get an image shape and check properties of it
    9-2. _get_image_shape: failure case: attempt to read an image unexisted, an exception raises
    9-3. _get_image_shape: failure case: attempt to read an image has invalid shape, an exception raises
    """

    input_image = 'expected.im'
    g192_im = 'g192_a2.image'

    def setUp(self):
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.g192_im)

    def test_9_1(self):
        """9-1. _get_image_shape: successful case: get an image shape and check properties of it."""
        shape = _get_image_shape(self.input_image)
        self.assertTrue(np.array_equal(shape.im_shape, [20, 20, 100]))
        self.assertTrue(np.array_equal(shape.axis_dir, [0, 1]))
        self.assertEqual(shape.axis_sp, 2)
        self.assertEqual(shape.axis_pol, -1)
        self.assertTrue(np.array_equal(shape.dir_shape, [20, 20]))
        self.assertEqual(shape.im_nrow, 400)
        self.assertEqual(shape.im_nchan, 100)
        self.assertEqual(shape.im_npol, 1)

        shape = _get_image_shape(self.g192_im)
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
        """9-2. _get_image_shape: failure case: attempt to read an image unexisted, an exception raises."""
        _get_image_shape('notexists')

    @test_base.exception_case(ValueError,
                              'image \'testdata_01.im\' is invalid')
    def test_9_3(self):
        """9-3. _get_image_shape: failure case: attempt to read an image has invalid shape, an exception raises."""
        testimage = 'testdata_01.im'
        self._create_image(testimage, 1.0, [64, 64])
        _get_image_shape(testimage)


class TestImbaseline(test_base):
    """imbaseline tests.

    F-2. failure case: attempt to read imagefile set None, an exception raises
    F-3. successful case: execute imbaseline with output_cont set False, cont file does not generate
    F-4. successful case: execute imbaseline with bloutput set output path, bloutput generates
    """

    input_image = 'ref_multipix.signalband'
    blparam = 'analytic_variable_blparam_spw1.txt'
    f_1_count = 1

    def setUp(self):
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.blparam)

    @test_base.exception_case(ValueError, 'Error: file  is not found.')
    def test_f_2(self):
        """F-2. failure case: attempt to read imagefile set None, an exception raises."""
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
        """F-3. successful case: execute imbaseline with output_cont set False, cont file does not generate."""
        imagefile = self.input_image
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
        """F-4. successful case: execute imbaseline with bloutput set output path, bloutput generates."""
        imagefile = self.input_image
        linefile = 'output_f_4'
        dirkernel = 'gaussian'
        spkernel = 'gaussian'
        major = '20arcsec'
        minor = '10arcsec'
        pa = '0deg'
        blfunc = 'sinusoid'
        output_cont = True
        bloutput = self.input_image + 'bloutput'

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


class TestImbaselineExecution(test_base):
    """Imbaseline execution testing.

    This test class generate tests dynamically and register them with the class
    while module initialisation.
    """

    input_image = 'ref_multipix.signalband'
    blparam = 'analytic_variable_blparam_spw1.txt'
    test_name_prefix = 'test_imbaseline_execution'
    linefile = 'output_f_1'
    output_cont = True
    bloutput = input_image + '.bloutput'
    chans = ''
    thresh = 5.0
    avg_limit = 5
    minwidth = 5
    edge = [0, 0]
    order = 5
    npiece = 3
    applyfft = True
    fftthresh = 3.0
    addwn = [0]
    rejwn = []
    clipniter = 0
    clipthresh = 3.0
    major = '20arcsec'
    minor = '10arcsec'
    pa = '0deg'
    kimage = os.path.join(DATAPATH, "bessel.im")
    scale = -1.0
    kwidth = 5
    filenames_existence_check = (linefile, bloutput)

    test_no = 1

    def setUp(self):
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.blparam)

    @staticmethod
    def generate_tests():
        maskmode = ('auto', 'list')
        blfunc = ('poly', 'chebyshev', 'cspline', 'sinusoid', 'variable')
        dirkernel = ('none', 'gaussian', 'boxcar', 'image')
        spkernel = ('none', 'gaussian', 'boxcar')

        def __register_a_test_with_the_class(
                _class, maskmode, blfunc, dirkernel, spkernel):
            test_name = f'{_class.test_name_prefix}_{maskmode}_{blfunc}_{dirkernel}_{spkernel}'
            setattr(
                _class,
                test_name,
                _class._generate_a_test(
                    maskmode,
                    blfunc,
                    dirkernel,
                    spkernel,
                    test_name))

        [__register_a_test_with_the_class(__class__, _maskmode, _blfunc, _dirkernel, _spkernel)
         for _maskmode in maskmode
         for _blfunc in blfunc
         for _dirkernel in dirkernel
         for _spkernel in spkernel]

    @staticmethod
    def _generate_a_test(maskmode, blfunc, dirkernel, spkernel, test_name):
        def test_method(self):
            """TestImbaselineExecution method No."""
            params = dict(
                imagename=self.input_image,
                linefile=self.linefile,
                output_cont=self.output_cont,
                bloutput=self.bloutput,
                maskmode=maskmode,
                chans=self.chans,
                thresh=self.thresh,
                avg_limit=self.avg_limit,
                minwidth=self.minwidth,
                edge=self.edge,
                blfunc=blfunc,
                order=self.order,
                npiece=self.npiece,
                applyfft=self.applyfft,
                fftthresh=self.fftthresh,
                addwn=self.addwn,
                rejwn=self.rejwn,
                blparam=self.blparam,
                clipniter=self.clipniter,
                clipthresh=self.clipthresh,
                dirkernel=dirkernel,
                major=self.major,
                minor=self.minor,
                pa=self.pa,
                kimage=self.kimage,
                scale=self.scale,
                spkernel=spkernel,
                kwidth=self.kwidth)
            casalog.post(
                f'{test_name} [maskmode={maskmode}, blfunc={blfunc}, '
                f'dirkernel={dirkernel}, spkernel={spkernel}]', 'INFO')
            imbaseline(**params)
            for file in self.filenames_existence_check:
                self.assertTrue(os.path.exists(file))
        test_method.__doc__ += (f'{TestImbaselineExecution.test_no:03} [maskmode={maskmode}, '
                                f'blfunc={blfunc}, dirkernel={dirkernel}, spkernel={spkernel}]')
        TestImbaselineExecution.test_no += 1
        return test_method


# generate test methods of TestImbaselineExecution dynamically
TestImbaselineExecution.generate_tests()


class TestImbaselineOutputs(test_base):
    """Imbaseline output testing.

    This test class generate tests dynamically and register them with the class
    while module initialisation.
    """

    input_image = 'ref_multipix.signalband'
    blparam = 'analytic_variable_blparam_spw1.txt'
    TEST_IMAGE_SHAPE = [128, 128, 1, 10]
    TEST_IMAGE_VALUE = 2.0
    expected_output_chunk = np.full(TEST_IMAGE_SHAPE, 0.0)
    expected_cont_chunk = np.full(TEST_IMAGE_SHAPE, TEST_IMAGE_VALUE)

    test_name_prefix = 'test_imbaseline_outputs'
    linefile = 'output.im'
    output_cont = True
    bloutput = linefile + '.bloutput'
    test_image = 'input_image.im'
    major = '20arcsec'
    minor = '10arcsec'
    pa = '0deg'

    test_no = 1

    def setUp(self):
        self._copy_test_files(self.input_image)
        self._copy_test_files(self.blparam)

    @staticmethod
    def generate_tests():
        blfunc = ('poly', 'chebyshev', 'cspline', 'sinusoid')
        dirkernel = ('none', 'gaussian', 'boxcar')
        spkernel = ('none', 'gaussian', 'boxcar')

        def __register_a_test_with_the_class(
                _class, blfunc, dirkernel, spkernel):
            test_name = f'{_class.test_name_prefix}_{blfunc}_{dirkernel}_{spkernel}'
            setattr(
                _class,
                test_name,
                _class._generate_a_test(
                    blfunc,
                    dirkernel,
                    spkernel,
                    test_name))

        [__register_a_test_with_the_class(__class__, _blfunc, _dirkernel, _spkernel)
         for _blfunc in blfunc
         for _dirkernel in dirkernel
         for _spkernel in spkernel]

    @staticmethod
    def _generate_a_test(blfunc, dirkernel, spkernel, test_name):
        def test_method(self):
            """TestImbaselineOutputs test."""
            self._create_image(
                self.test_image, self.TEST_IMAGE_VALUE, self.TEST_IMAGE_SHAPE)
            params = dict(
                imagename=self.test_image,
                linefile=self.linefile,
                output_cont=self.output_cont,
                bloutput=self.bloutput,
                blfunc=blfunc,
                dirkernel=dirkernel,
                spkernel=spkernel,
                major=self.major,
                minor=self.minor,
                pa=self.pa)
            casalog.post(
                f'{test_name} [maskmode=auto, blfunc={blfunc}, '
                f'dirkernel={dirkernel}, spkernel={spkernel}]', 'INFO')
            imbaseline(**params)
            if os.path.exists(self.linefile):
                with tool_manager(self.linefile, image) as ia:
                    chunk = ia.getchunk()
                    self._summary(False, test_name, chunk)
                    self.assertTrue(np.allclose(
                        chunk, self.expected_output_chunk, atol=2.0))
            if os.path.exists(self.test_image + '.cont'):
                with tool_manager(self.test_image + '.cont', image) as ia:
                    chunk = ia.getchunk()
                    self._summary(True, test_name, chunk)
                    self.assertTrue(np.allclose(
                        chunk, self.expected_cont_chunk, atol=2.0))

        test_method.__doc__ += (f' No.{TestImbaselineOutputs.test_no:03} [blfunc={blfunc}, '
                                f'dirkernel={dirkernel}, spkernel={spkernel}]')
        TestImbaselineOutputs.test_no += 1
        return test_method

    def _summary(self, is_cont, test_name, chunk):  # temporary method
        m = re.match(
            r'test_imbaseline_outputs_([^_]+)_([^_]+)_([^_]+)', test_name)
        prefix = 'cont' if is_cont else ' out'
        print(
            f'{prefix} blfunc:{m[1]} dirkernel:{m[2]} spkernel:{m[3]}, {np.max(chunk)}, '
            f'{np.min(chunk)}, {np.average(chunk)}, {np.median(chunk)}')


# generate test methods of TestImbaselineOutputs dynamically
TestImbaselineOutputs.generate_tests()


class Workdir():
    """Workdir manipulation class."""

    def __init__(self, path: str):
        """Initialize function."""
        self.path = path

    def clean(self, dry_run: bool = False):
        if os.path.exists(self.path) and dry_run is False:
            try:
                shutil.rmtree(self.path)
            except Exception:
                casalog.post(
                    'Some errors occured when clearing work dir', 'SEVERE')
                raise RuntimeError

    def chdir(self):
        if os.path.exists(self.path):
            os.chdir(self.path)
        else:
            casalog.post(
                'Some errors occured when chdir to work dir', 'SEVERE')
            raise RuntimeError

    @classmethod
    def create(cls, parent_path: str = None) -> 'Workdir':
        def is_valid_parent_path(path):
            return \
                os.path.exists(path) and \
                os.path.isdir(path) and \
                os.access(path, os.W_OK)

        if not parent_path:
            parent_path = os.getcwd()

        if is_valid_parent_path(parent_path):
            path = parent_path
            while True:
                path = os.path.join(parent_path, str(uuid.uuid4()))
                if not os.path.exists(path):
                    os.mkdir(path)
                    casalog.post(f'created working directory: {path}', 'WARN')
                    break
            if path == parent_path:
                raise RuntimeError(
                    'Some errors occured when creating work dir')
            return Workdir(path)

        raise RuntimeError('Some errors occured when creating work dir')


def suite():
    """Unittest suite definition."""
    return [
        TestImsmooth,
        TestAbstractFileStack,
        Test_ImageShape,
        TestImbaseline,
        TestImage2MS,
        TestSdbaseline,
        TestSdsmooth,
        TestMS2Image,
        TestImageSubtraction,
        TestModuleMethodsOfImbaseline,
        TestImbaselineExecution,
        TestImbaselineOutputs]


if __name__ == '__main__':
    unittest.main()

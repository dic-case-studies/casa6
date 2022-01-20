import os
import re
import shutil
import unittest

import numpy as np
from casatasks import casalog
from casatasks.private.sdutil import calibrater_manager, table_manager
from casatasks.private.task_imbaseline import *
from casatools import ctsys, image, quanta, regionmanager, table

_ia = image()
_rg = regionmanager()
_tb = table()
_qa = quanta()
ctsys_resolve = ctsys.resolve


class test_base(unittest.TestCase):

    @staticmethod
    def invalid_argument_case(func):
        """
        Decorator for the test case that is intended to fail
        due to invalid argument.
        """
        import functools

        @functools.wraps(func)
        def wrapper(self):
            func(self)
            self.assertFalse(self.result, msg='The task must return False')
        return wrapper

    @staticmethod
    def exception_case(exception_type, exception_pattern):
        """
        Decorator for the test case that is intended to throw
        exception.

            exception_type: type of exception
            exception_pattern: regex for inspecting exception message
                               using re.search
        """
        def wrapper(func):
            import functools

            @functools.wraps(func)
            def _wrapper(self):
                self.assertTrue(len(exception_pattern) > 0, msg='Internal Error')
                with self.assertRaises(exception_type) as ctx:
                    func(self)
                    self.fail(msg='The task must throw exception')
                the_exception = ctx.exception
                message = str(the_exception)
                self.assertIsNotNone(re.search(exception_pattern, message), msg='error message \'%s\' is not expected.' % (message))
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
        _base = ctsys_resolve(basename)
        copy_from = os.path.join(_base, filename)
        if not os.path.exists(copy_from) or copy_from == os.path.join(os.getcwd(), filename):
            raise RuntimeError(f"Error is occured or existed on a path {copy_from} or {filename}")

        if os.path.exists(filename):
            shutil.rmtree(filename)
        os.system('cp -RH ' + os.path.join(_base, filename) + ' ' + filename)


def _near(got, expected, tol):
    return _qa.le(
        _qa.div(
            _qa.abs(_qa.sub(got, expected)),
            expected
        ),
        tol
    )


class AbstractFileStack_test(test_base):
    """AbstractFileStack / (Un)EraseableFolder test

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
    """ImageShape test

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

    @test_base.exception_case(ValueError, 'nchan \d is too few to perform baseline subtraction')
    def test_2_2(self):
        shape = ImageShape(im_shape=np.array([100, 100, 1, 1]), axis_dir=np.array([0, 1]), axis_sp=3, axis_pol=2)
        shape.validate()

    @test_base.exception_case(ValueError, 'invalid value: dir_shape \[\d+\]')
    def test_2_3(self):
        shape = ImageShape(np.array([100, 100, 1, 100]), np.array([0]), 3, 2)
        shape.validate()


class imsmooth_test(test_base):
    """imsmooth test

    Tests of imsmooth rely on test_imsmooth basically, so we have minimal tests in imbaseline.

    3-1. simple successful case
    3-2. simple failure case
    """

    datapath = ctsys_resolve('unittest/imsmooth/')
    tiny = "tiny.im"

    def setUp(self):
        self._copy_test_files(self.datapath, self.tiny)

    def test_3_1(self):
        major = "2.5arcsec"
        minor = "2arcsec"
        pa = "0deg"
        dirkernel = "gaussian"
        kimage = ''
        scale = -1

        _stack = CasaImageStack(top=UnerasableFolder(self.tiny))

        execute_imsmooth(dirkernel, major, minor, pa, kimage, scale, _stack)
        self.assertTrue(os.path.exists(_stack.peak().path))

        _stack.pop()
        execute_imsmooth(dirkernel, major, minor, pa, kimage, scale, _stack)
        self.assertTrue(os.path.exists(_stack.peak().path))

        _stack.clear()
        try:
            _stack.push(UnerasableFolder("nonexists"))
        except Exception:
            pass
        self.assertFalse(_stack.height() > 0)

    @test_base.exception_case(ValueError, 'Unsupported direction smoothing kernel, foobar')
    def test_3_2(self):
        major = "2.5arcsec"
        minor = "2arcsec"
        pa = "0deg"
        dirkernel = "foobar"
        kimage = ''
        scale = -1

        _stack = CasaImageStack(top=UnerasableFolder(self.tiny))

        execute_imsmooth(dirkernel, major, minor, pa, kimage, scale, _stack)


class image2ms_test(test_base):
    """image2ms test

    4-1. simple successful case
    4-2. invalid datacolumn
    4-3. invalid image
    4-4. set empty stack
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected = "expected.im"
    dummy_folder1 = "dummy1"
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
        execute_image2ms(self.datacolumn, self.image_shape, image_stack, ms_stack)
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
        execute_image2ms('INVALID', self.image_shape, image_stack, ms_stack)

    @test_base.exception_case(RuntimeError, 'Unable to open image dummy1.')
    def test_4_3(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.dummy_folder1))
        ms_stack = MeasurementSetStack()
        execute_image2ms(self.datacolumn, self.image_shape, image_stack, ms_stack)

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_4_4(self):
        image_stack = CasaImageStack()
        ms_stack = MeasurementSetStack()
        execute_image2ms(self.datacolumn, self.image_shape, image_stack, ms_stack)


class sdsmooth_test(test_base):
    """sdsmooth test

    Tests of sdsmooth rely on test_sdsmooth basically, so we have minimal tests in imbaseline.

    5-1. simple successful case
    5-2. invalid ms stack
    5-3. invalid image stack
    """
    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = "expected.im"
    expected_ms = "expected.ms"
    dummy_folder1 = "dummy1"
    datacolumn = DATACOLUMN
    spkenel = "gaussian"
    kwidth = 5

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_ms)
        self.image_shape = get_image_shape(os.path.join(self.datapath, self.expected_im))

    def test_5_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        execute_sdsmooth(self.datacolumn, self.spkenel, self.kwidth, image_stack, ms_stack, self.image_shape)
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
        execute_sdsmooth(self.datacolumn, self.spkenel, self.kwidth, image_stack, ms_stack, self.image_shape)

    @test_base.exception_case(RuntimeError, 'the stack has not have enough stuff')
    def test_5_3(self):
        image_stack = CasaImageStack()
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        execute_sdsmooth(self.datacolumn, self.spkenel, self.kwidth, image_stack, ms_stack, self.image_shape)


class sdbaseline_test(test_base):
    """sdbaseline test

    Tests of sdbaseline rely on test_sdbaselinebasically, so we have minimal tests in imbaseline.

    6-1. simple successful case
    6-2. invalid ms stack
    6-3. invalid image stack
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = "expected.im"
    expected_ms = "expected.ms"
    bloutput = "test.csv"
    maskmode = "auto"
    chans = ""
    thresh = 5.0
    avg_limit = 4
    minwidth = 4
    edge = [0, 0]
    blfunc = "cspline"
    order = 5
    npiece = 1
    applyfft = True
    fftthresh = 3.0
    addwn = [0]
    rejwn = []
    blparam = ''
    clipniter = 10
    clipthresh = 2.0
    datacolumn = DATACOLUMN

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_ms)
        self.image_shape = get_image_shape(os.path.join(self.datapath, self.expected_im))

    def test_6_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        execute_sdbaseline(self.datacolumn, self.bloutput, self.maskmode, self.chans, self.thresh, self.avg_limit, self.minwidth,
                           self.edge, self.blfunc, self.order, self.npiece, self.applyfft, self.fftthresh, self.addwn, self.rejwn, self.blparam,
                           self.clipniter, self.clipthresh, image_stack, ms_stack, self.image_shape)
        self.assertTrue(os.path.exists(ms_stack.peak().path))
        self.assertTrue(os.path.exists(self.bloutput))
        self.assertTrue(os.path.exists(image_stack.peak().path))

    @test_base.exception_case(RuntimeError, 'the stack is empty')
    def test_6_2(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        ms_stack = MeasurementSetStack()
        execute_sdbaseline(self.datacolumn, self.bloutput, self.maskmode, self.chans, self.thresh, self.avg_limit, self.minwidth,
                           self.edge, self.blfunc, self.order, self.npiece, self.applyfft, self.fftthresh, self.addwn, self.rejwn, self.blparam,
                           self.clipniter, self.clipthresh, image_stack, ms_stack, self.image_shape)

    @test_base.exception_case(RuntimeError, 'the stack has not have enough stuff')
    def test_6_3(self):
        image_stack = CasaImageStack()
        ms_stack = MeasurementSetStack()
        ms_stack.push(EraseableFolder(self.expected_ms))
        execute_sdbaseline(self.datacolumn, self.bloutput, self.maskmode, self.chans, self.thresh, self.avg_limit, self.minwidth,
                           self.edge, self.blfunc, self.order, self.npiece, self.applyfft, self.fftthresh, self.addwn, self.rejwn, self.blparam,
                           self.clipniter, self.clipthresh, image_stack, ms_stack, self.image_shape)


class image_subtraction_test(test_base):
    """Image subtraction test
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected_im = "expected.im"
    expected_imsmoothed = "expected.imsmooth.im"
    expected_bl = "expected.bl.im"

    def setUp(self):
        self._copy_test_files(self.datapath, self.expected_im)
        self._copy_test_files(self.datapath, self.expected_imsmoothed)
        self._copy_test_files(self.datapath, self.expected_bl)

    def test_7_1(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        image_stack.push(EraseableFolder(self.expected_imsmoothed))
        image_stack.push(EraseableFolder(self.expected_bl))
        output = "output_7_1.im"
        execute_image_subtraction(output, image_stack)
        self.assertTrue(os.path.exists(output))

    def test_7_2(self):
        image_stack = CasaImageStack(top=UnerasableFolder(self.expected_im))
        image_stack.push(EraseableFolder(self.expected_bl))
        output = "output_7_2.im"
        execute_image_subtraction(output, image_stack)
        self.assertTrue(os.path.exists(output))


class imbaseline_test(test_base):
    """Full test.

    x-1. simple successful case
    """

    datapath = ctsys_resolve('unittest/imbaseline/')
    expected = "expected.im"

    def setUp(self):
        self.ia = image()
        self._copy_test_files(self.datapath, self.expected)

    def test_x_1(self):
        imagefile = self.expected
        linefile = 'output_7_1'
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


def suite():
    return [imsmooth_test, AbstractFileStack_test, ImageShape_test, imbaseline_test, image2ms_test, sdbaseline_test, image_subtraction_test]


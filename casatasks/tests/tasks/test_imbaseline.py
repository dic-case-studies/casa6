import math
import os
import random
import re
import shutil
import sys
import unittest

import numpy as np
from casatasks import casalog
from casatasks.private.sdutil import calibrater_manager, table_manager
from casatasks.private.task_imbaseline import (AbstractFileStack, CasaImageStack,
                                               EraseableFolder, ImageShape,
                                               UnerasableFolder,
                                               execute_imsmooth, imbaseline)
from casatools import (componentlist, ctsys, image, ms, quanta, regionmanager,
                       table)
from scipy import signal

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

    def tearDown(self):
        if os.path.exists(self.dummy_folder1):
            shutil.rmtree(self.dummy_folder1)
        if os.path.exists(self.dummy_folder2):
            shutil.rmtree(self.dummy_folder2)

    def test_1_1(self):
        stack = CasaImageStack(UnerasableFolder(self.dummy_folder1))
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(ValueError, 'file path unexists is not found')
    def test_1_2(self):
        CasaImageStack(UnerasableFolder(self.unexist_folder))

    def test_1_3(self):
        stack = CasaImageStack()
        stack.push(UnerasableFolder(self.dummy_folder1))
        self.assertTrue(stack.height() == 1)

    @test_base.exception_case(ValueError, 'file path unexists is not found')
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
        stack = CasaImageStack(EraseableFolder(self.dummy_folder1))
        stack.clear(False)
        self.assertFalse(os.path.exists(self.dummy_folder1))
        self.assertEqual(stack.height(), 0)

    def test_1_14(self):
        stack = CasaImageStack(UnerasableFolder(self.dummy_folder1))
        stack.clear(False)
        self.assertTrue(os.path.exists(self.dummy_folder1))
        self.assertEqual(stack.height(), 0)

    def test_1_15(self):
        file = EraseableFolder(self.dummy_folder1)
        file.erase(False)
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

    3-1. simple successful case
    3-2. simple failure case
    """

    datapath = ctsys_resolve('unittest/imsmooth/')
    targetres_im = "imsmooth_targetres.fits"
    tiny = "tiny.im"
    tiny_dummy = tiny + ".dummy"
    image_names = ['g192_a2.image', 'g192_a2.image-2.rgn']

    def setUp(self):
        if(os.path.exists(self.image_names[0])):
            for file in self.image_names:
                os.system('rm -rf ' + file)

        for file in self.image_names:
            os.system('cp -RH ' + os.path.join(self.datapath, file) + ' ' + file)
        self.ia = image()
        for f in [self.targetres_im, self.tiny]:
            if(os.path.exists(f)):
                os.system('rm -rf ' + f)
            os.system('cp -RH ' + os.path.join(self.datapath, f) + ' ' + f)
        self._create_dummy()

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

    def _create_dummy(self):
        if os.path.exists(self.tiny_dummy):
            os.system('rm -rf ' + self.tiny_dummy)
        os.system('cp -RH ' + os.path.join(self.datapath, self.tiny) + ' ' + self.tiny_dummy)

    ####################################################################
    # Incorrect inputs to parameters.  The parameters are:
    #    imagename
    #    outfile
    #    kernel
    #    major
    #    minor
    #    mask
    #    region
    #    box
    #    chans
    #    stokes
    #
    # Returns True if successful, and False if it has failed.
    ####################################################################

    def _compare_beams(self, beam1, beam2):
        self.assertTrue(_near(beam1["major"], beam2["major"], 2e-5))
        self.assertTrue(_near(beam1["minor"], beam2["minor"], 2e-5))
        pa = []
        for b in [beam1, beam2]:
            if "positionangle" in b:
                pa.append(b["positionangle"])
            else:
                pa.append(b["pa"])

        diff = abs(
            _qa.sub(
                _qa.quantity(pa[0]),
                _qa.quantity(pa[1])
            )["value"]
        )
        self.assertTrue(diff < 1e-5)

    def test_input(self):
        '''Imsmooth: Testing INPUT/OUTPUT tests'''
        casalog.post("Starting imsmooth INPUT/OUTPUT tests.", 'NORMAL2')

        #######################################################################
        # Testing the imagename parameter.
        #    1. Bad file name should throw and exception
        #    2. Good file name, a file should be
        #######################################################################
        casalog.post("The IMAGENAME parameter tests will cause errors to occur, do not be alarmed", 'WARN')
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


    ####################################################################
    # image2ms test
    #   4-1. simple successful case
    ####################################################################


    ####################################################################
    # sdsmooth test
    #   5-1. simple successful case
    #   5-2. simple failure case
    ####################################################################


    ####################################################################
    # ms2image test
    #   6-1. simple successful case
    #   6-2. simple failure case
    ####################################################################


    ####################################################################
    # full test
    #   7-1. simple successful case
    #   7-2. simple failure case
    ####################################################################

class imbaseline_test(test_base):

    datapath = ctsys_resolve('unittest/imsmooth/')
    targetres_im = "imsmooth_targetres.fits"
    tiny = "tiny.im"
    tiny_dummy = tiny + ".dummy"
    image_names = ['g192_a2.image', 'g192_a2.image-2.rgn']

    def setUp(self):
        if(os.path.exists(self.image_names[0])):
            for file in self.image_names:
                os.system('rm -rf ' + file)

        for file in self.image_names:
            os.system('cp -RH ' + os.path.join(self.datapath, file) + ' ' + file)
        self.ia = image()
        for f in [self.targetres_im, self.tiny]:
            if(os.path.exists(f)):
                os.system('rm -rf ' + f)
            os.system('cp -RH ' + os.path.join(self.datapath, f) + ' ' + f)
        self._create_dummy()

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

    def _create_dummy(self):
        if os.path.exists(self.tiny_dummy):
            os.system('rm -rf ' + self.tiny_dummy)
        os.system('cp -RH ' + os.path.join(self.datapath, self.tiny) + ' ' + self.tiny_dummy)

    def test_7_1(self):
        imagefile = self.image_names[0]
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
    return [imsmooth_test, AbstractFileStack_test, ImageShape_test, imbaseline_test]


import random
import os
import numpy
import re
import shutil
import sys
import unittest
import math
from scipy import signal

from casatools import ctsys, image, regionmanager, componentlist, table, quanta ,ms
from casatasks import casalog
from casatasks.private.sdutil import table_manager, calibrater_manager
from casatasks.private.task_imbaseline import ProcessingFileStack, Unerasable, execute_imsmooth

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
                self.assertIsNotNone(re.search(exception_pattern, message), msg='error message \'%s\' is not expected.'%(message))
            return _wrapper
        return wrapper


datapath = ctsys_resolve('unittest/imsmooth/')
targetres_im = "imsmooth_targetres.fits"
tiny = "tiny.im"
image_names = ['g192_a2.image', 'g192_a2.image-2.rgn']


def _near(got, expected, tol):
    return _qa.le(
        _qa.div(
            _qa.abs(_qa.sub(got, expected)),
            expected
        ),
        tol
    )


class imsmooth_test(unittest.TestCase):

    def setUp(self):
        if(os.path.exists(image_names[0])):
            for file in image_names:
                os.system('rm -rf ' + file)

        for file in image_names:
            os.system('cp -RH ' + os.path.join(datapath, file) + ' ' + file)
        self.ia = image()
        for f in [targetres_im, tiny]:
            if(os.path.exists(f)):
                os.system('rm -rf ' + f)
            os.system('cp -RH ' + os.path.join(datapath, f) + ' ' + f)

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

        processing_file_stack = ProcessingFileStack(top=Unerasable(tiny))

        execute_imsmooth(dirkernel, major, minor, pa, kimage, scale, processing_file_stack)
        self.assertTrue(os.path.exists(processing_file_stack.top()))

        processing_file_stack.pop()
        execute_imsmooth(dirkernel, major, minor, pa, kimage, scale, processing_file_stack)
        self.assertTrue(os.path.exists(processing_file_stack.top()))

        processing_file_stack.clear()
        try:
            processing_file_stack.push(Unerasable("nonexists"))
        except Exception:
            pass
        self.assertFalse(processing_file_stack.length() > 0)


def suite():
    return [imsmooth_test]

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
from casatasks import imbaseline, casalog, imsubimage
from casatasks.private.sdutil import table_manager, calibrater_manager

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


if __name__ == '__main__':
    os.chdir("/work/dev/shimada/casa6.13520.new/tmp")
    unittest.main()
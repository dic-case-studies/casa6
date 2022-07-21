##########################################################################
# test_tool_logsink.py
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
# [CAS-13813]
# Created new method 'getorigin' of logsink.
#
##########################################################################

import unittest

import numpy as np

import casatools


class LogsinkTest(unittest.TestCase):
    """The unittest class of logsink."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getorigin_initial_origin(self):
        """Check the initial value of logsink.origin by getorigin."""
        self.assertEqual(casatools.logsink().getorigin(), '')

    def __getorigin_subtest(self, test_cases):
        """Run subtest of getorigin with test cases.

        Parameters
        ----------
        test_cases : list of tuple of {str, numeric, None}
            List of test cases. Each test case is a tuple of origin and expected value.
        """
        casalog = casatools.logsink()
        for input_origin, expected in test_cases:
            with self.subTest(origin=input_origin):
                casalog.origin(input_origin)
                self.assertEqual(casalog.getOrigin(), expected)

    def test_getorigin_set_strings(self):
        """Check that the string values set logsink.origin can get correctly by getorigin."""
        self.__getorigin_subtest(
            [
                ('test', 'test'),
                ('test\n', 'test\n')
            ]
        )

    def test_getorigin_set_nullvalues(self):
        """Check that the null values set logsink.origin can get correctly by getorigin."""
        self.__getorigin_subtest(
            [
                ('\0', ''),
                (None, 'None'),
                ('', '')
            ]
        )

    def test_getorigin_set_num(self):
        """Check that the numeric values set logsink.origin can get correctly by getorigin."""
        self.__getorigin_subtest(
            [
                (1, '1'),
                (np.pi, str(np.pi)),
                (1 + 1j, str(1 + 1j))
            ]
        )


if __name__ == '__main__':

    unittest.main()

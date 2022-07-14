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


class logsink_test(unittest.TestCase):
    """The unittest class of logsink."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getorigin_initial_origin(self):
        """Check the initial value of logsink.origin by getorigin."""
        self.assertEqual(casatools.logsink().getorigin(), '')

    def test_getorigin_set_strings(self):
        """Check that the string values set logsink.origin can get correctly by getorigin."""

        casalog = casatools.logsink()
        casalog.origin('test')
        self.assertEqual(casalog.getorigin(), 'test')
        casalog.origin('test\n')
        self.assertEqual(casalog.getorigin(), 'test\n')
        del casalog

    def test_getorigin_set_nullvalues(self):
        """Check that the null values set logsink.origin can get correctly by getorigin."""

        casalog = casatools.logsink()
        casalog.origin('\0')
        self.assertEqual(casalog.getorigin(), '')
        casalog.origin(None)
        self.assertEqual(casalog.getorigin(), 'None')
        casalog.origin('')
        self.assertEqual(casalog.getorigin(), '')
        del casalog

    def test_getorigin_set_num(self):
        """Check that the numeric values set logsink.origin can get correctly by getorigin."""

        casalog = casatools.logsink()
        casalog.origin(1)
        self.assertEqual(casalog.getorigin(), '1')
        casalog.origin(np.pi)
        self.assertEqual(casalog.getorigin(), str(np.pi))
        casalog.origin(1 + 1j)
        self.assertEqual(casalog.getorigin(), str(1 + 1j))
        del casalog


if __name__ == '__main__':

    unittest.main()

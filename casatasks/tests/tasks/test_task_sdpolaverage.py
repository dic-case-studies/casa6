########################################################################
# test_task_sdpolaverage.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.single.sdpolaverage.html
#
#
##########################################################################
import math
import os
import sys
import unittest

from casatasks import sdpolaverage
from casatasks.private.sdutil import table_manager
from casatools import ctsys

datapath = ctsys.resolve('unittest/sdpolaverage/')


def weighToSigma(weight):
    if weight > sys.float_info.min:
        return 1.0 / math.sqrt(weight)
    else:
        return -1.0


def sigmaToWeight(sigma):
    if sigma > sys.float_info.min:
        return 1.0 / math.pow(sigma, 2)
    else:
        return 0.0


def check_eq(val, expval, tol=None):
    """Check that val matches expval within tol."""
#    print val
    if type(val) == dict:
        for k in val:
            check_eq(val[k], expval[k], tol)
    else:
        try:
            if tol and hasattr(val, '__rsub__'):
                are_eq = abs(val - expval) < tol
            else:
                are_eq = val == expval
            if hasattr(are_eq, 'all'):
                are_eq = are_eq.all()
            if not are_eq:
                raise ValueError('!=')
        except ValueError:
            errmsg = "%r != %r" % (val, expval)
            if (len(errmsg) > 66):  # 66 = 78 - len('ValueError: ')
                errmsg = "\n%r\n!=\n%r" % (val, expval)
            raise ValueError(errmsg)
        except Exception as e:
            print("Error comparing {} to {}".format(val, expval))
            raise e


class test_sdpolaverage(unittest.TestCase):
    def setUp(self):
        self.inputms = "analytic_type1.fit.ms"
        self.outputms = "polave.ms"
        os.system('cp -RH ' + datapath + self.inputms + ' ' + self.inputms)

    def tearDown(self):
        os.system('rm -rf ' + self.inputms)
        os.system('rm -rf ' + self.outputms)

    def test_default(self):
        sdpolaverage(infile=self.inputms, outfile=self.outputms, datacolumn='float_data')
        with table_manager(self.inputms) as tb:
            indata = tb.getcell('FLOAT_DATA', 0)
        with table_manager(self.outputms) as tb:
            outdata = tb.getcell('FLOAT_DATA', 0)

        self.assertEqual(len(indata), len(outdata), 'Input and output data have different shape.')
        for i in range(len(indata)):
            for j in range(len(indata[0])):
                self.assertEqual(indata[i][j], outdata[i][j], 'Input and output data unidentical.')

    def test_stokes_float_data(self):
        sdpolaverage(infile=self.inputms, outfile=self.outputms,
                     polaverage='stokes', datacolumn='float_data')
        # check data
        with table_manager(self.inputms) as tb:
            indata = tb.getcell('FLOAT_DATA', 0)
        with table_manager(self.outputms) as tb:
            outdata = tb.getcell('FLOAT_DATA', 0)

        self.assertEqual(len(outdata), 1, 'No averaging over polarization?')
        tol = 1e-5
        for i in range(len(indata[0])):
            mean = 0.5 * (indata[0][i] + indata[1][i])
            check_eq(outdata[0][i], mean, tol)

        # check polarization id (should be 1)
        with table_manager(self.outputms) as tb:
            outddesc = tb.getcell('DATA_DESC_ID', 0)
        with table_manager(self.outputms + '/DATA_DESCRIPTION') as tb:
            outpolid = tb.getcol('POLARIZATION_ID')
        with table_manager(self.outputms + '/POLARIZATION') as tb:
            outpoltype = tb.getcell('CORR_TYPE', outpolid[outddesc])

        self.assertEqual(len(outpoltype), 1, 'Polarization id is inconsistent with data.')
        self.assertEqual(outpoltype[0], 1, 'Has wrong polarization id.')

    def test_stokes_corrected_data(self):
        sdpolaverage(infile=self.inputms, outfile=self.outputms,
                     polaverage='stokes', datacolumn='corrected')
        # check data
        with table_manager(self.inputms) as tb:
            indata = tb.getcell('CORRECTED_DATA', 0)
        with table_manager(self.outputms) as tb:
            outdata = tb.getcell('DATA', 0)

        self.assertEqual(len(outdata), 1, 'No averaging over polarization?')
        tol = 1e-5
        for i in range(len(indata[0])):
            mean = 0.5 * (indata[0][i] + indata[1][i])
            check_eq(outdata[0][i].real, mean.real, tol)
            check_eq(outdata[0][i].imag, mean.imag, tol)

        # check polarization id (should be 1)
        with table_manager(self.outputms) as tb:
            outddesc = tb.getcell('DATA_DESC_ID', 0)
        with table_manager(self.outputms + '/DATA_DESCRIPTION') as tb:
            outpolid = tb.getcol('POLARIZATION_ID')
        with table_manager(self.outputms + '/POLARIZATION') as tb:
            outpoltype = tb.getcell('CORR_TYPE', outpolid[outddesc])

        self.assertEqual(len(outpoltype), 1, 'Polarization id is inconsistent with data.')
        self.assertEqual(outpoltype[0], 1, 'Has wrong polarization id.')


if __name__ == '__main__':
    unittest.main()

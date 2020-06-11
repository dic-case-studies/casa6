##########################################################################
# test_req_task_rerefant.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_rerefant/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import rerefant, casalog, fringefit

    CASA6 = True
    tb = casatools.table()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy as np

# import pyfits

## DATA ##

if CASA6:
    casavis = casatools.ctsys.resolve('visibilities/vla/ngc5921.ms/')
    casacal = casatools.ctsys.resolve('caltables/ngc5921.ref1a.gcal')
    # old test path
    # datadir = casatools.ctsys.resolve('/data/regression/evn/')
    src = casatools.ctsys.resolve('visibilities/other/n08c1.ms')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        casavis = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/ngc5921.ms/'
        casacal = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ngc5921.ref1a.gcal'
        # old test path
        src = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/other/n08c1.ms'
    else:
        casavis = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/ngc5921.ms/'
        casacal = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ngc5921.ref1a.gcal'
        # old test path
        src = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/other/n08c1.ms'

copyvis = 'vis.ms'
copycal = 'copycal.gcal'

logpath = casalog.logfile()


def file_copy(filename, perm):
    os.chmod(filename, perm)
    for root, dirs, files in os.walk(filename):
        for d in dirs:
            os.chmod(os.path.join(root, d), perm)
        for f in files:
            os.chmod(os.path.join(root, f), perm)


class rerefant_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(rerefant)
        shutil.copytree(casavis, copyvis)
        shutil.copytree(casacal, copycal)
        file_copy(copyvis, 493)
        file_copy(copycal, 493)

    def tearDown(self):
        shutil.rmtree(copyvis)
        shutil.rmtree(copycal)
        casalog.setlogfile(logpath)

        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
        if os.path.exists('out.cal'):
            shutil.rmtree('out.cal')
        if os.path.exists('n08c1_reref.ms'):
            shutil.rmtree('n08c1_reref.ms')
        if os.path.exists("fringe.cal"):
            shutil.rmtree("fringe.cal")
        if os.path.exists("reref.cal"):
            shutil.rmtree("reref.cal")

    @classmethod
    def tearDownClass(cls):
        pass

    def test_takesms(self):
        '''
            test_taskesms
            ---------------

            Test that the task takes a valid ms.

            To operate this task also required a tablein to be given as well but this test only checks valid vs invalid ms
        '''

        casalog.setlogfile('testlog.log')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01')
        self.assertTrue('SEVERE' not in open('testlog.log').read())
        shutil.rmtree('out.cal')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01')
        self.assertTrue('SEVERE' not in open('testlog.log').read())

    def test_takestablein(self):
        '''
            test_takestablein
            -------------------

            Test that the task takes a valid caltable

            Check that the invalid inputs will raise a sever error

            TODO once the multiple antenna option is fixed test further
        '''

        casalog.setlogfile('testlog.log')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01')
        self.assertTrue('SEVERE' not in open('testlog.log').read())
        shutil.rmtree('out.cal')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01')
        self.assertTrue('SEVERE' not in open('testlog.log').read())

    def test_modeflex(self):
        '''
            test_modeflex
            ------------------

            Check that the flex parameter doesn't raise any SEVERE errors

            Check that it takes multiple antenna (and can switch antennas and back)

            TODO Need to be able to check if the antenna drops out and back in but it's not working with multiple antenna right now
        '''

        casalog.setlogfile('testlog.log')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA02', refantmode='flex')
        self.assertTrue('SEVERE' not in open('testlog.log').read())

        tb.open('out.cal')
        antennacol = tb.getcol('ANTENNA2')
        self.assertTrue(np.all(antennacol == 1))
        tb.close()

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01,VA02', refantmode='flex')
        self.assertTrue('SEVERE' not in open('testlog.log').read())

        tb.open('out.cal')
        antennacol = tb.getcol('ANTENNA2')
        self.assertTrue(np.all(antennacol == 0))
        tb.close()

    def test_modestrict(self):
        '''
            test_modestrict
            -----------------

                This mode will flag all antennas if the current refant is absent for a solution

                If a list is provided only use the first element of that list

                TODO come back to this one (It may be that the muti antenna problem has to be solved before I can test this fully)
        '''

        casalog.setlogfile('testlog.log')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01', refantmode='strict')
        self.assertTrue('SEVERE' not in open('testlog.log').read())

        tb.open('out.cal')
        # print(tb.getcol('FLAG'))
        tb.close()

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01,VA02', refantmode='strict')
        self.assertTrue('SEVERE' not in open('testlog.log').read())

    def test_takescaltable(self):
        '''
            test_takescaltable
            --------------------

            Test that the tasks makes an output caltable when ran
        '''

        casalog.setlogfile('testlog.log')

        rerefant(vis=copyvis, tablein=copycal, caltable='out.cal', refant='VA01')
        self.assertTrue('SEVERE' not in open('testlog.log').read())
        self.assertTrue(os.path.exists('out.cal'))

    # old test case

    def test_fringefit(self):
        dst = "n08c1_reref.ms"
        calfile = "fringe.cal"
        recalfile = "reref.cal"

        shutil.copytree(src, dst)

        # Perform fringe fit with EF as the reference station.
        fringefit(vis=dst, caltable=calfile, refant="EF")
        self.assertTrue(os.path.exists(calfile))

        # Rereference the results using ON as the new reference station.
        rerefant(vis=dst, tablein=calfile, caltable=recalfile, refant="ON")
        self.assertTrue(os.path.exists(recalfile))

        # Check original calibration table.
        tb.open(calfile)
        ant2 = tb.getcol('ANTENNA2')
        taql = "ANTENNA1 == 3 && ANTENNA2 == 0"
        tsel = tb.query(taql)
        fparam = tsel.getcol('FPARAM')
        tb.close()
        # Reference anttena is EF, aka antenna number 0.
        self.assertTrue(len(set(ant2)) == 1)
        self.assertTrue(0 in set(ant2))

        # Check rereferenced calibration table.
        tb.open(recalfile)
        ant2 = tb.getcol('ANTENNA2')
        taql = "ANTENNA1 == 0 && ANTENNA2 == 3"
        tsel = tb.query(taql)
        refparam = tsel.getcol('FPARAM')
        tb.close()
        # Reference anttena is ON, aka antenna number 3.
        self.assertTrue(len(set(ant2)) == 1)
        self.assertTrue(3 in set(ant2))

        # Parameters on EF-ON baseline should be opposite.
        self.assertTrue(np.isclose(refparam, -fparam, 1e-15).all())

        # Rereferenced parameters for ON should all be zero.
        tb.open(recalfile)
        taql = "ANTENNA1 == 3 && ANTENNA2 == 3"
        tsel = tb.query(taql)
        refparam = tsel.getcol('FPARAM')
        tb.close()
        self.assertTrue(not refparam.any())

        # Rereferenced parameters for EF should not be zero.
        tb.open(recalfile)
        taql = "ANTENNA1 == 0 && ANTENNA2 == 3"
        tsel = tb.query(taql)
        refparam = tsel.getcol('FPARAM')
        tb.close()
        self.assertTrue(refparam.any())


def suite():
    return [rerefant_test]


if __name__ == '__main__':
    unittest.main()


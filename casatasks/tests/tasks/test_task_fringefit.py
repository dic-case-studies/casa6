#########################################################################
# test_task_fringefit.py
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
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.calibration.fringefit.html
#
##########################################################################
import os
import sys
import shutil
import unittest
import itertools
import numpy as np

# For information about parameters that are unexpectedly zero, set
# VERBOSE to true.  Currently there are none, so this is for developers
# only
VERBOSE = False

from casatools import ms, ctsys, table
from casatasks import fringefit, flagmanager, flagdata
from casatestutils import testhelper as th

tblocal = table()

datapath = ctsys.resolve('unittest/fringefit/')

class Fringefit_tests(unittest.TestCase):
    prefix = 'n08c1'
    msfile = prefix + '.ms'
    uvfile = 'gaincaltest2copy.ms'

    def setUp(self):
        shutil.copytree(os.path.join(datapath, self.msfile), self.msfile)
        shutil.copytree(os.path.join(datapath, 'gaincaltest2.ms'), self.uvfile)

    def tearDown(self):
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.prefix + '.sbdcal', True)
        shutil.rmtree(self.prefix + '-zerorates.sbdcal', True)
        shutil.rmtree(self.prefix + '.mbdcal', True)
        shutil.rmtree(self.uvfile, True)
        shutil.rmtree('uvrange_with.cal', True)

    def test_sbd(self):
        sbdcal = self.prefix + '.sbdcal'
        fringefit(vis=self.msfile, caltable=sbdcal, refant='EF')
        reference = os.path.join(datapath, sbdcal)
        self.assertTrue(th.compTables(sbdcal, reference, ['WEIGHT', 'SNR']))

    def test_mbd(self):
        sbdcal = self.prefix + '-zerorates.sbdcal'
        mbdcal = self.prefix + '.mbdcal'
        fringefit(vis=self.msfile, caltable=sbdcal, field='4C39.25',
                  refant='EF', zerorates=True)
        fringefit(vis=self.msfile, caltable=mbdcal, field='J0916+3854',
                   combine='spw', gaintable=[sbdcal], refant='EF')
        reference = os.path.join(datapath, mbdcal)
        self.assertTrue(th.compTables(mbdcal, reference, ['WEIGHT', 'SNR']))

    def test_uvrange(self):
        ''' Check that the uvrnage parameter excludes antennas '''
        # create a caltable with uvrange selection
        fringefit(vis=self.uvfile, caltable='uvrange_with.cal', spw='2', refant='0', uvrange='<1160')

        # get the subset of antennas that are used vs all
        tblocal.open('uvrange_with.cal')
        output = tblocal.getcol('FLAG')
        antennas = tblocal.getcol('ANTENNA1')
        tblocal.close()

        flagged = set()
        intended_flagged = {5,8}

        for i in range(len(antennas)):
            if np.all(output[:, :, i] == True):
                flagged.add(antennas[i])

        self.assertTrue(flagged == intended_flagged)


class Fringefit_single_tests(unittest.TestCase):
    prefix = 'n08c1-single'
    msfile = prefix + '.ms'

    def setUp(self):
        shutil.copytree(os.path.join(datapath, self.msfile), self.msfile)

    def tearDown(self):
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.prefix + '.sbdcal', True)
        shutil.rmtree(self.prefix + '-2.sbdcal', True)

    def test_single(self):
        sbdcal = self.prefix + '.sbdcal'
        fringefit(vis=self.msfile, caltable=sbdcal, refant='EF')
        tblocal.open(sbdcal)
        flag = tblocal.getcol('FLAG')
        tblocal.close()
        # CAS-12693: Check that the right parameters are flagged
        self.assertFalse(flag[0, 0, 0])
        self.assertFalse(flag[0, 0, 1])
        self.assertTrue(flag[0, 0, 2])
        self.assertFalse(flag[0, 0, 3])
        self.assertTrue(flag[0, 0, 4])
        self.assertTrue(flag[4, 0, 0])
        self.assertTrue(flag[4, 0, 1])
        self.assertTrue(flag[4, 0, 2])
        self.assertTrue(flag[4, 0, 3])
        self.assertTrue(flag[4, 0, 4])
    def test_param(self):
        sbdcal = self.prefix + '-2.sbdcal'
        # We make a triple cartesian product of booleans
        # to test all possible paramactive values
        eps = 1e-20
        refant = 0
        refant_s = str(refant) 
        for pactive in itertools.product(*3*[[False, True]]):
            fringefit(vis=self.msfile, paramactive=list(pactive), caltable=sbdcal, refant=refant_s)
            tblocal.open(sbdcal)
            fparam = tblocal.getcol('FPARAM')
            flag = tblocal.getcol('FLAG')
            tblocal.close()
            param_names = ['delay', 'rate', 'dispersivity']
            if VERBOSE: 
                print(pactive, file=sys.stderr)
            # Loop over parameters
            for i in range(2):
                # Loop over stations; it seems like station 2 is the one with non zero results
                for j in range(4):
                    if not pactive[i]:
                        self.assertTrue(abs(fparam[i+1, 0, j]) < eps )
                        self.assertTrue(abs(fparam[i+5, 0, j]) < eps)
                    else:
                        # Obviously the reference antenna show have zero values for all parameters!
                        if j==refant: continue
                        # We don't report dispersion being zero when it
                        # is included in the parameters to solve
                        # because this branch doesn't yet *solve* for
                        # dispersion
                        if i==3: continue
                        if (abs(fparam[i+1, 0, j]) < eps) and (not flag[i+1,0,j]):
                            name = param_names[i]
                            v1 = fparam[i+1, 0, j]
                            v1 = fparam[i+5, 0, j]
                            if VERBOSE: 
                                print("   Parameter {} for antenna {} is {}, {}".format(name, j, v1, v1),
                                      "when it doesn't have to be",
                                      file=sys.stderr)


class Fringefit_dispersive_tests(unittest.TestCase):
    prefix = 'n14p1'
    msfile = prefix+'.ms'
    
    def setUp(self):
        shutil.copytree(os.path.join(datapath, self.msfile), self.msfile)
        flagdata(self.prefix + '.ms', mode='manual', spw='*:0~2;29~31')

    def tearDown(self):
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.prefix + '.ms' + '.flagversions')
        shutil.rmtree(self.prefix + '.mpc', True)
        shutil.rmtree(self.prefix + '.disp', True)

    def test_manual_phase_cal(self):
        fringefit(vis="n14p1.ms", caltable="n14p1.mpc",
                  scan="1", solint="300", refant="WB",
                  minsnr=50, zerorates=True,
                  globalsolve=True, niter=100, gaintable=[],
                  parang=True)
        mpcal = self.prefix + '.mpc'
        reference = os.path.join(datapath, mpcal)
        self.assertTrue(th.compTables(mpcal, reference, ['WEIGHT', 'SNR']))
        fringefit(vis=self.prefix + '.ms',
                  caltable='n14p1.disp', refant="WB",
                  scan='1', solint='60', spw='0,1,2',
                  paramactive=[True, True, True],
                  minsnr=50,
                  niter=100,
                  gaintable=['n14p1.mpc'],
                  parang=True)
        dispcal = self.prefix + '.disp'
        reference = os.path.join(datapath, dispcal)
        self.assertTrue(th.compTables(dispcal, reference, ['WEIGHT', 'SNR']))


class Fringefit_refant_bookkeeping_tests(unittest.TestCase):
    prefix = 'n08c1-single'
    msfile = prefix + '.ms'
    sbdcal = prefix + '-book.sbdcal'

    def setUp(self):
        shutil.copytree(os.path.join(datapath, self.msfile), self.msfile)
        flagdata(self.prefix + '.ms', mode='manual', spw='*:0~2;29~31')
        flagdata(self.prefix + '.ms', mode='manual', antenna='EF')

    def tearDown(self):
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.msfile + '.flagversions')
        shutil.rmtree(self.sbdcal, True)

    def test_bookkeeping(self):
        eps = 1e-2
        refant_ind = 1
        fringefit(vis=self.msfile, caltable=self.sbdcal, refant='WB')
        tblocal.open(self.sbdcal)
        fparam = tblocal.getcol('FPARAM')
        flag = tblocal.getcol('FLAG')
        tblocal.close()
        for i in range(8):
            self.assertTrue(abs(fparam[i, 0, refant_ind]) < eps )
            self.assertTrue(abs(fparam[i, 0, refant_ind]) < eps)

class FreqMetaTests(unittest.TestCase):
    prefix = 'n08c1'
    msfile = prefix + '.ms'

    def setUp(self):
        shutil.copytree(os.path.join(datapath, self.msfile), self.msfile)

    def tearDown(self):
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.prefix + '.mbdcal', True)

    def test_metadata(self):
        sbdcal = self.prefix + '-zerorates.sbdcal'
        mbdcal = self.prefix + '.mbdcal'
        fringefit(vis=self.msfile, caltable=sbdcal, field='4C39.25',
                  refant='EF', zerorates=True)
        fringefit(vis=self.msfile, caltable=mbdcal, spw="0,1,2,3", field='J0916+3854', timerange="17:10:00~17:11:00",
                   combine='spw', gaintable=[sbdcal], refant='EF')
        tblocal.open(mbdcal + '/SPECTRAL_WINDOW')
        flagrow = tblocal.getcol('FLAG_ROW')
        # There are only 4 spws in n08c1!
        tblocal.close()
        self.assertTrue((flagrow==[False] + 3*[True]).all())
        # Only two scans in unit test ms
        # /scratch/small/UnitTest/n08c1.ms
        # 10-Mar-2008/17:06:00.0 - 17:09:00.0     1      3 4C39.25                   4320  [0,1,2,3]  [1, 1, 1, 1] 
        #             17:09:00.0 - 17:11:00.0     2      2 J0916+3854                2880  [0,1,2,3]  [1, 1, 1, 1] 
        try:
            fringefit(vis=self.msfile, caltable=mbdcal, spw="0:6~32,1,2,3", field='J0916+3854', timerange="17:10:00~17:11:00",
                      combine='spw', gaintable=[sbdcal], refant='EF', append=True)
            print("In test_metadata, a fringefit which should have thrown an exception did not!")
            self.assertTrue(False)
        except RuntimeError as e: 
            print(e)
            self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()

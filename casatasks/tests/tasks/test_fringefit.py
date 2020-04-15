from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import unittest
import itertools

# is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import ms, ctsys, table
    from casatasks import fringefit

    tblocal = table()
    ctsys_resolve = ctsys.resolve
else:
    from __main__ import default
    from tasks import *
    from taskinit import tbtool

    dataRoot = os.path.join(os.environ.get('CASAPATH').split()[0],'data')
    tblocal = tbtool()

    def ctsys_resolve(apath):
        return os.path.join(dataRoot,apath)
        
if is_python3:
    ### for testhelper import
    sys.path.append(os.path.abspath(os.path.dirname(__file__)))
    import testhelper as th
else:
    import testhelper as th

datapath = ctsys_resolve('regression/evn')

class Fringefit_tests(unittest.TestCase):
    prefix = 'n08c1'
    msfile = prefix + '.ms'

    def setUp(self):
        shutil.copytree(os.path.join(datapath, self.msfile), self.msfile)

    def tearDown(self):
        shutil.rmtree(self.msfile)
        shutil.rmtree(self.prefix + '.sbdcal', True)
        shutil.rmtree(self.prefix + '-zerorates.sbdcal', True)
        shutil.rmtree(self.prefix + '.mbdcal', True)

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
        for pactive in itertools.product(*3*[[False, True]]):
            fringefit(vis=self.msfile, paramactive=pactive, caltable=sbdcal, refant='EF')
            tblocal.open(sbdcal)
            fparam = tblocal.getcol('FPARAM')
            flag = tblocal.getcol('FLAG')
            tblocal.close()
            param_names = ['delay', 'rate', 'dispersivity']
            print(pactive, file=sys.stderr)
            # Loop over parameters
            for i in range(3):
                # Loop over stations; it seems like station 2 is the one with non zero results
                for j in range(4):
                    if not pactive[i]:
                        self.assertTrue(abs(fparam[i+1, 0, j]) < eps )
                        self.assertTrue(abs(fparam[i+5, 0, j]) < eps)
                    else:
                        if (abs(fparam[i+1, 0, j]) < eps) and (not flag[i+1,0,j]):
                            name = param_names[i]
                            v1 = fparam[i+1, 0, j]
                            v1 = fparam[i+5, 0, j]
                            print("   Parameter {} for antenna {} is {}, {}".format(name, j, v1, v1),
                                  "when it doesn't have to be",
                                  file=sys.stderr)

def suite():
    return [Fringefit_tests, Fringefit_single_tests]

if __name__ == '__main__':
    unittest.main()

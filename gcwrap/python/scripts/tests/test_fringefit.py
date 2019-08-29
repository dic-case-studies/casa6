from tasks import fringefit
import os
import shutil
import unittest

import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import testhelper as th

datapath = os.environ.get('CASAPATH').split()[0] +\
                            '/data/regression/evn/'

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
        self.assertTrue(th.compTables(sbdcal, reference, ['WEIGHT']))

    def test_mbd(self):
        sbdcal = self.prefix + '-zerorates.sbdcal'
        mbdcal = self.prefix + '.mbdcal'
        fringefit(vis=self.msfile, caltable=sbdcal, field='4C39.25',
                  refant='EF', zerorates=True)
        fringefit(vis=self.msfile, caltable=mbdcal, field='J0916+3854',
                   combine='spw', gaintable=[sbdcal], refant='EF')
        reference = os.path.join(datapath, mbdcal)
        self.assertTrue(th.compTables(mbdcal, reference, ['WEIGHT']))

def suite():
    return [Fringefit_tests]

if __name__ == '__main__':
    unittest.main()

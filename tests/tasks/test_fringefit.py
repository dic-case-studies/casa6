from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import unittest

# is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import ms, ctsys
    from casatasks import fringefit

    ctsys_resolve = ctsys.resolve
else:
    from __main__ import default
    from tasks import *
    from taskinit import *

    dataRoot = os.path.join(os.environ.get('CASAPATH').split()[0],'data')

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

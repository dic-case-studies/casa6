from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import unittest

# is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import ms, ctsys, table
    from casatasks import fringefit, flagmanager, flagdata

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

    
def suite():
    return [Fringefit_tests, Fringefit_single_tests]

if __name__ == '__main__':
    unittest.main()

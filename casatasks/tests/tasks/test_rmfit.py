import math
import shutil
import unittest

from casatasks import rmfit

from casatools import imagepol as potool
from casatools import image as iatool
from casatools import table
from casatools import constants

datapath='regression/unittest/po_tool/'
eq_beams = datapath + "pol_eq_beams.fits"
neq_beams = datapath + "pol_neq_beams.fits"

class rmfit_test(unittest.TestCase):
    
    def setUp(self):
        self.mypo = potool()
    
    def tearDown(self):
        self.mypo.done()
        tb = table( )
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done( )
    
    def test_rmfit_basics(self):
        """Sanity tests for task rmfit"""
        myia = iatool()
        outfile = "xx.im"
        myia.fromshape(outfile, [20, 20, 4, 20])
        myia.addnoise()
        myia.done()
        myrm = "rm1.im"
        self.assertTrue(rmfit(imagename=outfile, rm=myrm))
        myia.open(myrm)
        self.assertTrue((myia.shape() == [20, 20]).all())
        got1 = myia.statistics(list=True, verbose=True)['sumsq']
        myia.done()
        
        # test concatenation of images
        outfile = "yy.im"
        myia.fromshape(outfile, [20, 20, 4, 20])
        myia.addnoise()
        csys = myia.coordsys()
        refval = csys.referencevalue()['numeric']
        refval[3] = 1.5e9
        csys.setreferencevalue(refval)
        myia.setcoordsys(csys.torecord())
        myia.done()
        images = ["xx.im", "yy.im"]
        myrm = "rm2.im"
        self.assertTrue(rmfit(imagename=images, rm=myrm))
        myia.open(myrm)
        self.assertTrue((myia.shape() == [20, 20]).all())
        got2 = myia.statistics(list=True, verbose=True)['sumsq']
        myia.done()
        self.assertTrue(abs(got1 - got2) > 0.1)

    def test_algorithm(self):
        """Test rotation measure computation algorithm"""
        myia = iatool()
        imagename = "rm_input.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        csys = myia.coordsys()
        incr = csys.increment()['numeric']
        incr[3] = 1000*incr[3]
        csys.setincrement(incr)
        myia.setcoordsys(csys.torecord())
        pixvals = myia.getchunk()
        # U values all 1
        U = 1
        pixvals[:,:,2,:] = U
        c = constants.c/100
        RM = 9.6
        pa0deg = 22.5
        pa0 = pa0deg/180*math.pi
        for chan in range(myia.shape()[3]):
            freq = myia.toworld([0,0,0,chan])['numeric'][3]
            lam = c/freq
            Q = U/math.tan(2*(pa0 + RM*lam*lam))
            pixvals[:,:,1,chan] = Q
        myia.putchunk(pixvals)
        myia.done()
        mypo = self.mypo
        rmim = "rm.im"
        pa0im = "pa0.im"
        sigma = 10e-8
        rmfit(imagename=imagename, rm=rmim, pa0=pa0im, sigma=sigma)
        myia.open(rmim)
        stats = myia.statistics(list=True, verbose=True)
        self.assertTrue((abs(stats['min'][0] - RM)) < 1e-4)
        self.assertTrue((abs(stats['max'][0] - RM)) < 1e-4)
        myia.done(remove=True)
        myia.open(pa0im)
        stats = myia.statistics(list=True, verbose=True)
        self.assertTrue((abs(stats['min'][0] - pa0deg)) < 1e-4)
        self.assertTrue((abs(stats['max'][0] - pa0deg)) < 1e-4)
        myia.done(remove=True)

def suite():
    return [rmfit_test]

if __name__ == '__main__':
    unittest.main()

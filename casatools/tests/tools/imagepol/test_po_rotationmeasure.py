import shutil
import unittest
import math

try:
    from casatools import constants
    from casatools import imagepol as potool
    from casatools import image as iatool
    from casatools import table
    from casatools import ctsys
    ctsys_resolve = ctsys.resolve
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'data')
        return os.path.join(dataPath,apath)

datapath = ctsys_resolve('regression/unittest/po_tool/')
eq_beams = datapath + "pol_eq_beams.fits"
neq_beams = datapath + "pol_neq_beams.fits"

class po_rotationmeasure_test(unittest.TestCase):
    
    def setUp(self):
        self.mypo = potool()
    
    def tearDown(self):
        self.mypo.done()
        tb = table( )
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done( )
    
    def test_multibeam(self):
        """Test multibeam images for correct behavior"""
        mypo = self.mypo
        mypo.open(eq_beams)
        self.assertTrue(mypo.rotationmeasure("g"))
        mypo.done()
        mypo.open(neq_beams)
        self.assertRaises(Exception, mypo.rotationmeasure, "hh")
        mypo.done()
        
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
        mypo.open(imagename)
        mypo.rotationmeasure(rm=rmim, pa0=pa0im, sigma=sigma)
        mypo.done()
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
    return [po_rotationmeasure_test]

if __name__ == '__main__':
    unittest.main()

import unittest
import copy
import math
from casatools import measures
from casatools import quanta

class me_shift_test(unittest.TestCase):

    def setUp(self):
        self.me = measures( )
        self.qa = quanta( )
        pass
    
    def tearDown(self):
        self.me.done( )
        pass

    def test_shift(self):
        """Test me.shift"""
        v = self.me.direction("J2000", "4h20m30s", "+30.20.30")
        got = self.me.shift(v, "20arcmin", "0deg")
        expec = copy.deepcopy(v)
        expec['m1'] = self.qa.add(expec['m1'], "20arcmin")
        self.assertTrue(got == expec)
        got = self.me.shift(v, "20arcmin", "90deg")
        expec = 1.1433867531223854
        self.assertTrue(abs(got['m0']['value']/expec - 1) < 1e-7)
        expec = 0.5295520783025025
        self.assertTrue(abs(got['m1']['value']/expec - 1) < 1e-7)
        got = self.me.shift(v, "20arcmin", "180deg")
        self.assertTrue(got['m0']['value'] == v['m0']['value'])
        expec = self.qa.sub(v['m1'], '20arcmin')
        self.assertTrue(abs(got['m1']['value']/expec['value'] - 1) < 1e-7)
        
def suite():
    return [me_shift_test]

if __name__ == '__main__':
    unittest.main()

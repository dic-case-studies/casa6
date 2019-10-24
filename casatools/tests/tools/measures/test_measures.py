import unittest
from casatools import measures

'''
Unit tests for the measures tool.

Features tested:
  1. me.cometdist
  2. me.cometangdiam
'''

def check_eq(val, expval, tol=None):
    """Checks that val matches expval within tol."""
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
            if (len(errmsg) > 66): # 66 = 78 - len('ValueError: ')
                errmsg = "\n%r\n!=\n%r" % (val, expval)
            raise ValueError(errmsg)
        except Exception as e:
            print("Error comparing %s to %s" % (val,expval))
            raise e

class Ganymede(unittest.TestCase):
    """
    Base class for Ganymede tests.
    """
    def setUp(self):
        self.me = measures( )
        cometdir = "ephemerides/JPL-Horizons/"
        self.me.framecomet(cometdir + "Ganymede_55437-56293dUTC.tab")
        self.me.doframe(self.me.epoch("utc", "2011/01/03/17:00:00"))
        self.me.doframe(self.me.observatory("ALMA"))

    def tearDown(self):
        self.me.done( )
        pass

class me_test_cometdist(Ganymede):
    def test_cometdist(self):
        """Is the distance to Ganymede just right?"""
        check_eq(self.me.cometdist(), {'unit': 'AU', 'value': 5.1241}, 0.001)
            
class me_test_cometangdiam(Ganymede):
    def test_cometangdiam(self):
        """Is Ganymede's angular diameter correct?"""
        check_eq(self.me.cometangdiam(), {'unit': 'rad', 'value': 6.868e-06},
                 1.0e-9)

def suite():
    return [me_test_cometdist, me_test_cometangdiam]

if __name__ == '__main__':
    unittest.main()

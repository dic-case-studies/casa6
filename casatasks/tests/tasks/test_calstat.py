import os
import unittest
from casatools import table,ctsys
from casatasks import calstat

epsilon = 0.0001

class calstat_test(unittest.TestCase):
    def setUp(self):
        self.epsilon = 0.0001

        self.caltables = [ 'ggtau.1mm.amp.gcal',
                           'ggtau.1mm.bpoly',
                           'ggtau.1mm.ph.gcal',
                           'ggtau.1mm.ph.gcal0',
                           'ggtau.3mm.amp.gcal',
                           'ggtau.3mm.bpoly',
                           'ggtau.3mm.ph.gcal',
                           'ggtau.3mm.ph.gcal0',
                           'ggtau.co.bpoly',
                           'ggtau.hco.bpoly' ]
        self.tb = table( )

    def tearDown(self):
        self.tb.done( )

    def data(self):
        return map( lambda x: ctsys.resolve("unittest/calstat/" + x), self.caltables )

    def test_cs(self):
        expected = {'ggtau.3mm.ph.gcal0':
                    {'SPLINE_KNOTS_PHASE':{'rms': 4362063360.0,
                                           'medabsdevmed': 13056.0,
                                           'min': 4362050048.0,
                                           'max': 4362076160.0,
                                           'sum': 872412620800.0,
                                           'quartile': 26112.0,
                                           'median': 4362063104.0,
                                           'sumsq': 3.80551890468e+21,
                                           'stddev': 11866.4301499,
                                           'var': 140812164.503,
                                           'npts': 200,
                                           'mean': 4362063104.0}},
                    'ggtau.1mm.ph.gcal0':
                    {'SPLINE_KNOTS_PHASE':{'rms': 4362063360.0,
                                           'medabsdevmed': 13056.0,
                                           'min': 4362050048.0,
                                           'max': 4362076160.0,
                                           'sum': 872412620800.0,
                                           'quartile': 26112.0,
                                           'median': 4362063104.0,
                                           'sumsq': 3.80551890468e+21,
                                           'stddev': 11866.4301499,
                                           'var': 140812164.503,
                                           'npts': 200,
                                           'mean': 4362063104.0}}}

        for caltable in self.data( ):

            print("Testing with data", caltable, "...")
            name = os.path.basename(caltable)
            if name in expected:

                axis='spline_knots_phase'
                s = calstat(caltable=caltable, axis=axis)

                if s.keys() != expected[name].keys():
                    raise Exception("Wrong dictionary keys. Expected %s, got %s" % \
                                    (expected[name], s))
                            
                print("Expected =", expected[name])
                print("Got = ", s)
                if not 'SPLINE_KNOTS_PHASE' in s:
                    raise Exception("Dictionary returned from calstat does not have key SPLINE_KNOTS_PHASE")

                for e in expected[name]['SPLINE_KNOTS_PHASE'].keys():
                    print("Checking %s: %s vs %s" % \
                        (e, expected[name]['SPLINE_KNOTS_PHASE'][e], s['SPLINE_KNOTS_PHASE'][e]))
                    failed = False
                    if expected[name]['SPLINE_KNOTS_PHASE'][e] == 0:
                        if s['SPLINE_KNOTS_PHASE'][e] != 0:
                            failed = True
                    else:
                        if abs((expected[name]['SPLINE_KNOTS_PHASE'][e] - s['SPLINE_KNOTS_PHASE'][e])/expected[name]['SPLINE_KNOTS_PHASE'][e]) > epsilon:
                            failed = True

                    # Remove these 3 lines of code, once CAS-1671 is solved
                    if failed == True and e in ['var', 'stddev']:
                        print("Ignoring this known problem on 64bit!")
                        failed = False

                    
                    if failed:
                        print("test failed...")
                        raise Exception("Numbers differ, expected %s, got %s" % \
                                        (str(expected[name]['SPLINE_KNOTS_PHASE'][e]), str(s['SPLINE_KNOTS_PHASE'][e])))
                
            self.tb.open(caltable)
            cols = self.tb.colnames()
            self.tb.close()

            cplx = ['amp', 'amplitude', 'phase', 'imag', 'imaginary', 'real']
            for x in cplx:
                cols.append(x)
            print(cols)
            # remove complex columns
            cols.remove('GAIN')
            if 'SCALE_FACTOR' in cols: cols.remove('SCALE_FACTOR')
            if 'SIDEBAND_REF' in cols: cols.remove('SIDEBAND_REF')
            # don't try string columns
            cols.remove('FREQ_GROUP_NAME')
            cols.remove('FIELD_NAME')
            cols.remove('FIELD_CODE')
            cols.remove('SOURCE_NAME')
            cols.remove('SOURCE_CODE')
            if 'POLY_TYPE' in cols: cols.remove('POLY_TYPE')
            if 'POLY_MODE' in cols: cols.remove('POLY_MODE')
            if 'PHASE_UNITS' in cols: cols.remove('PHASE_UNITS')

            # empty column:
            if 'VALID_DOMAIN' in cols: cols.remove('VALID_DOMAIN')

            cols = [x.lower() for x in cols]

            print("Trying these column names", cols)

            for col in cols:
                data_cols = ['']
                if col in cplx:
                    data_cols = ['gain', 'scale_factor']
                
                for dc in data_cols:
                    print("Call with caltable =", caltable, "; axis =", col, "; datacolumn =", dc)
                    if dc != '':
                        s = calstat(caltable=caltable, axis=col, datacolumn=dc)
                    else:
                        s = calstat(caltable=caltable, axis=col)
                    if col.upper() == "FLAG_CATEGORY":
                        # The MSs used have no data in FLAG_CATEGORY, therefore
                        # calstat() should fail
                        if s != None:
                            raise Exception("Error! " + str(s))
                    elif not type(s) is dict:
                        raise Exception("Error! Return value " + str(s) + " is not a dictionary")

            self.assertTrue(True)

def suit():
    return [calstat_test]

if __name__ == '__main__':
    unittest.main()

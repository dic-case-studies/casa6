##########################################################################
# test_req_task_calstat.py
#
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_calstat/about
#
# Test_logreturn checks to make sure a logfile is generated and populated
# Test_dictreturn checks that the result is a python dict object containing keys specified in the documentation
# Test_takescal checks that a caltable is accepted and non-cal tables are rejected
# Test_axis checks that different axis vaules will provide different information
# Test_axisvals checks that the values for axis provided in the documentatin are accepted as valid values
# Test_datacolumn checks that different datacolumn values provide different information
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import casalog, calstat
    from casatools import table, ctsys
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)
    table = tb
import sys
import os
import unittest
import shutil

if CASA6:
    datapath = casatools.ctsys.resolve('caltables/ggtau.1mm.amp.gcal')
    datapath_visibilities = casatools.ctsys.resolve('visibilities/')
    #filepath = casatools.ctsys.resolve('testlog.log')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ggtau.1mm.amp.gcal'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ggtau.1mm.amp.gcal'
    datapath_visibilities = os.environ.get('CASAPATH').split()[0] + '/data/visibilities/'
    #filepath = os.environ.get('CASAPATH').split()[0] + '/bin/nosedir/testlog.log'
        
logpath = casalog.logfile()
contained = ['max', 'mean', 'medabsdevmed', 'median', 'min', 'npts', 'quartile', 'rms', 'stddev', 'sum', 'sumsq', 'var']
epsilon = 0.0001

class calstat_test(unittest.TestCase):
     
    def setUp(self):
        if not CASA6:
            default(calstat)

        self.epsilon = 0.0001
        self.caltables = ['ggtau.1mm.amp.gcal',
                          'ggtau.1mm.bpoly',
                          'ggtau.1mm.ph.gcal',
                          'ggtau.1mm.ph.gcal0',
                          'ggtau.3mm.amp.gcal',
                          'ggtau.3mm.bpoly',
                          'ggtau.3mm.ph.gcal',
                          'ggtau.3mm.ph.gcal0',
                          'ggtau.co.bpoly',
                          'ggtau.hco.bpoly']
        self.tb = table()
     
    def tearDown(self):
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
        self.tb.done()
     
    def test_logreturn(self):
        '''logreturn test: Test that a logfile is written and populated with the expected information'''
        casalog.setlogfile('testlog.log')
        # Check that a logfile is populated and contains the correct infromation
        calstat(datapath, axis='TIME')
        for item in contained:
            self.assertTrue( item in open('testlog.log').read(), msg='Fails to write required information to the log')
        
    def test_dictreturn(self):
        '''dictionary test: Test that calstat makes a python dict with the expected keys'''
        # Check that type of returned object is a python dict
        caldict = calstat(datapath)
        self.assertTrue(isinstance(caldict, dict), msg='calstat does not return a python dict')
        # Check that the dict contains the correct values
        self.assertTrue(all (k in caldict['GAIN'] for k in contained), msg='Dictionalry does not contain all the correct values')
            
    def test_takescal(self):
        '''takes cal test: Test that calstat only takes cal tables'''
        self.assertTrue(calstat(datapath), msg='calstat fails to take a caltable')

        # No type checking for CASA 6, but some for CASA 5
        if CASA6 or casa_stack_rethrow:
            with self.assertRaises(RuntimeError, msg='Fails to recognize non-caltable'):
                calstat(datapath_visibilities)
        else:
            self.assertFalse(calstat(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/'), msg='Fails to recognize non-caltable')
            self.assertFalse(calstat(1), msg='Fails to refuse int type')
            self.assertFalse(calstat([]), msg='Fails to refuse a list')
            
    def test_axis(self):
        '''axis test: Check that the axis param makes different selections'''
        timeAxis = calstat(datapath, axis='TIME')
        ampAxis = calstat(datapath, axis='amp')
        self.assertFalse(timeAxis == ampAxis, msg='different axis selections should give different values')
        
    def test_axisvals(self):
        '''axis val test: Test that the inputs listed in the documentatin function properly'''
        # Some type checking for CASA 5 - this doesn't make much sense to me, it should also
        # check that exceptions are raised in CASA6 (or if casa_stack_rethrow)
        if not CASA6 and not casa_stack_rethrow:
            self.assertFalse(calstat(datapath, axis='abc'), msg='Fails to recognize non-existing axis')
            self.assertFalse(calstat(datapath, axis=1), msg='Takes int as input when it should not')
            self.assertFalse(calstat(datapath, axis=[]), msg='Takes list as input when it should not')
        # Check that the documented inputs function
        for item in ['amp', 'amplitude', 'phase', 'real', 'imag', 'imaginary']:
            self.assertTrue(calstat(datapath, axis=item), msg='axis {} is not recognized'.format(item))
    
    def test_dataColumn(self):
        '''datacolumn test: Check that the datacolumn param makes unique selections'''
        # Some type checking for CASA 5
        if not CASA6 and not casa_stack_rethrow:
            self.assertFalse(calstat(datapath, axis='amp', datacolumn='abc'), msg='Fails to recognize non-existing datacolumn')
            self.assertFalse(calstat(datapath, axis='amp', datacolumn=1), msg='Takes an int as input when it should not')
            self.assertFalse(calstat(datapath, axis='amp', datacolumn=[]), msg='Takes list as input when it should not')
        # Check that datacolumn selects different data when it should
        timeCol = calstat(datapath, axis='amp', datacolumn='TIME')
        gainCol = calstat(datapath, axis='amp', datacolumn='GAIN')
        self.assertFalse(timeCol == gainCol, msg='different datacolumns should give different results')
        # Added for additional test coverage. These should correspond to CORRECTED_DATA and MODEL_DATA
        # These columns are not present for the current cal table
        #self.assertTrue(calstat(datapath, axis='amp', datacolumn='corrected'))
        #self.assertTrue(calstat(datapath, axis='amp', datacolumn='model'))

    # Test cases from test_calstat
    def data(self):
        return map( lambda x: ctsys.resolve("regression/calstat/" + x), self.caltables )

    def test_cs(self):
        expected = {'ggtau.3mm.ph.gcal0':
                        {'SPLINE_KNOTS_PHASE': {'rms': 4362063360.0,
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
                        {'SPLINE_KNOTS_PHASE': {'rms': 4362063360.0,
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

        for caltable in self.data():

            print("Testing with data", caltable, "...")
            name = os.path.basename(caltable)
            if name in expected:

                axis = 'spline_knots_phase'
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
                        if abs((expected[name]['SPLINE_KNOTS_PHASE'][e] - s['SPLINE_KNOTS_PHASE'][e]) /
                               expected[name]['SPLINE_KNOTS_PHASE'][e]) > epsilon:
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
        

def suite():
    return[calstat_test]

if __name__ == '__main__':
    unittest.main()

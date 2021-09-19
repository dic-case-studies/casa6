from __future__ import absolute_import
import os
import sys
import shutil
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, table
    from casatasks import accor,smoothcal

    _tb = table()
else:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)

    _tb = tbtool()

'''
Unit tests for task smoothcal. It tests the following parameters:
    vis:           wrong and correct values
    tablein:       wrong and correct values
    caltable:      existence of output
    field:         wrong field value; non-default value
    smmothtype:    unsupported value; non-default value
    smoothtime:    unsupported value; non-default values

    Other tests: check the values of column smoothed GAIN against reference.
                 check calibration tables produces by accor are smoothable.
'''
class smoothcal_test(unittest.TestCase):

    # Input and output names
    msfile = 'ngc1333_ut.ms'
    gcal = 'ngc1333_ut_nct.gcal'   # New format caltables
    ref = 'ngc1333_ut_nct.ref'
    res = None
    vlbams = 'ba123a.ms'
    accor = 'smoothcal_accor'
    out = 'smoothcal_test'

    def setUp(self):
        self.res = None
        if not is_CASA6:
            default(smoothcal)
        if is_CASA6:
            datapath = ctsys.resolve('unittest/smoothcal/')
        else:
            datapath = os.path.join(os.environ.get('CASAPATH').split()[0],'casatestdata/unittest/smoothcal/')

        shutil.copytree(os.path.join(datapath,self.msfile), self.msfile)
        shutil.copytree(os.path.join(datapath,self.gcal), self.gcal)
        shutil.copytree(os.path.join(datapath,self.ref), self.ref)
        shutil.copytree(os.path.join(datapath,self.vlbams), self.vlbams)
    
    def tearDown(self):
        if (os.path.exists(self.msfile)):
            os.system('rm -rf ' + self.msfile)
        if (os.path.exists(self.gcal)):
            os.system('rm -rf ' + self.gcal)
        if (os.path.exists(self.ref)):
            os.system('rm -rf ' + self.ref)
        if (os.path.exists(self.vlbams)):
            os.system('rm -rf ' + self.vlbams)
        if (os.path.exists(self.accor)):
            os.system('rm -rf ' + self.accor)
        if (os.path.exists(self.out)):
            os.system('rm -rf ' + self.out)
        
    def getvarcol(self,table,colname):
        '''Return the requested column'''
        _tb.open(table)
        col = _tb.getvarcol(colname)
        _tb.close()
        return col

    def test0(self):
        '''Test 0: Missing input table caught by parameter checking (exception thrown)
        '''
        # CASA5 returns False (depending on __rethrow_casa_exceptions)
        try:
            OK = False
            self.res = smoothcal()
            if not is_CASA6:
                OK = not self.res
        except:
            if is_CASA6 or casa_stack_rethrow:
                OK = True
        self.assertTrue(OK)

    def test1(self):
        """Test 1: Wrong input MS caught by parameter checking (exception thrown)
        """
        msfile = 'badmsfile'
        # CASA5 returns False (depending on __rethrow_casa_exceptions)
        try:
            OK = False
            self.res = smoothcal(vis=msfile,tablein=self.gcal,caltable=self.out)
            if not is_CASA6:
                OK = not self.res
        except:
            if is_CASA6 or casa_stack_rethrow:
                OK = True
        self.assertTrue(OK)

    def test2(self):
        """Test 2: Wrong input gcal caught by parameter checking (exception thrown)
        """
        gcal = 'badgcal'
        # CASA5 returns False
        try:
            OK = False
            self.res = smoothcal(vis=self.msfile,tablein=gcal,caltable=self.out)
            if not is_CASA6:
                OK = not self.res
        except:
            if is_CASA6 or casa_stack_rethrow:
                OK = True
        self.assertTrue(OK)

    def test3(self):
        """Test 3: Good input should return None"""
        self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.out))
        
    def test4(self):
        """Test 4: Unsupported smoothtype"""
        # CASA5 returns False
        try:
            OK = False
            self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,smoothtype='average')
            if not is_CASA6:
                OK = not self.res
        except:
            if is_CASA6 or casa_stack_rethrow:
                OK = True
        self.assertTrue(OK)

    def test5(self):
        '''Test 5: Non-default smoothtype'''
        self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,smoothtype='mean')
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.out))

    def test6(self):
        '''Test 6: Unsupported smoothtime'''
        # CASA5 returns False
        try:
            OK = False
            self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,smoothtime=-1)
            if not is_CASA6:
                OK = not self.res
        except:
            if is_CASA6 or casa_stack_rethrow:
                OK = True
        self.assertTrue(OK)

    def test7(self):
        '''Test 7: Non-default smoothtype and smoothtime'''
        self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,smoothtype='mean',
                         smoothtime=7200)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.out))

    def test8(self):
        '''Unsupported field values'''
        # CASA5 returns False
        try:
            OK = False
            self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,field='23~30')
            if not is_CASA6:
                OK = not self.res
        except:
            if is_CASA6 or casa_stack_rethrow:
                OK = True
        self.assertTrue(OK)

    def test9(self):
        '''Test 9: Non-default field selection'''
        self.res = smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,smoothtype='median',
                         smoothtime=5000.5,field='2')
#        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.out))

    def test10(self):
        '''Test10: Compare smoothed values with reference'''
        self.res=smoothcal(vis=self.msfile,tablein=self.gcal,caltable=self.out,smoothtype='mean',
                         smoothtime=7200.)
        self.assertEqual(self.res,None)
        refcol = self.getvarcol(self.ref, 'CPARAM')
        smcol = self.getvarcol(self.out, 'CPARAM')
        nrows = len(refcol)
        EPS = 1e-5;  # Logical "zero"
        # Loop over every row,pol and get the data
        for i in range(1,nrows,1) :
            row = 'r%s'%i     
            # polarization is 0-1
            for pol in range(0,2) :     
                refdata = refcol[row][pol]
                smdata = smcol[row][pol]
                self.assertTrue(abs(refdata - smdata) < EPS)

    def test11(self):
        '''Test11: Smooth accor table'''
        accor(vis=self.vlbams,caltable=self.accor,corrdepflags=True)
        self.res=smoothcal(vis=self.vlbams,tablein=self.accor,caltable=self.out)
        self.assertTrue(os.path.exists(self.out))

def suite():
    return [smoothcal_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

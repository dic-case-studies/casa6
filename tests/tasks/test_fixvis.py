import os
import sys
import shutil
import unittest

from casatools import ctsys, table, quanta
from casatasks import fixvis, imstat, tclean, split

'''
Unit tests for task fixvis.

Features tested:
  1. Do converted directions in the FIELD table have the right shape? 
  2. Does the phase center shifting result in the expected shifts?
  3. Does the distances parameter work

Note: The equinox_vis regression is a more general test of fixvis.
'''
inpms = '0420+417.ms'
outms = 'output.ms'
inpms2 = 'twocenteredpointsources.ms'
outms2 = 'testx.ms'

_tb = table( )
_qa = quanta( )

class fixvis_test1(unittest.TestCase):
    def setUp(self):
        res = None
        if not os.path.exists(inpms):
            shutil.copytree(ctsys.resolve(os.path.join('regression/0420+417',inpms)), inpms)
        if not os.path.exists(inpms2):
            shutil.copytree(ctsys.resolve(os.path.join('regression/unittest/fixvis',inpms2)), inpms2)
        shutil.rmtree(outms, ignore_errors=True)
        shutil.rmtree(outms2, ignore_errors=True)

    def tearDown(self):
        shutil.rmtree(inpms)
        shutil.rmtree(inpms2)
        shutil.rmtree(outms, ignore_errors=True)
        shutil.rmtree(outms2, ignore_errors=True)
        os.system('rm -rf test[yz]*')

    def _fixvis_and_get_stats(self, phasecent, dist=""):
        refcode = 'J2000'
        shutil.rmtree(outms2, ignore_errors=True)
        mystats = ''
        try:
            self.res = fixvis(inpms2, outms2, field='0', refcode=refcode,
                              phasecenter=phasecent, distances=dist)
            self.assertTrue(self.res)
            mystats = self._get_stats(0, 'testy')
        except Exception as e:
            print("*** Unexpected error *** %s" % e)

        return mystats

    def _get_stats(self, fld, imgname):
        os.system('rm -rf ' + imgname + '*')
        tclean( vis=outms2, imagename=imgname, field=str(fld), threshold='0.1mJy',
                deconvolver='clark', gridder='standard', mask='user', imsize=[128,128],
                cell=['0.10000008arcsec','0.10000008arcsec'], weighting='natural',
                uvtaper=[] )
        return imstat(imgname + '.image')

    def test1(self):
        '''Test1: Do converted directions in the FIELD table have the right shape?'''
        refcode = 'J2000'
        self.res = fixvis( inpms, outms, refcode=refcode )
        _tb.open(outms + '/FIELD')
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        def record_error(errmsg, retValue):
            """Helper function to print and update retValue on an error."""
            print("test_fixvis.test1: Error: %s" % errmsg)
            retValue['success'] = False
            retValue['error_msgs'] += errmsg + "\n"

        mscomponents = set( [ "table.dat", "table.f0", "table.f1", "table.f2",
                              "table.f3", "table.f4", "table.f5", "table.f6",
                              "table.f7", "table.f8", "table.f9", "table.f10",
                              "table.f11", "FIELD/table.dat", "FIELD/table.f0" ] )

        for name in mscomponents:
            if not os.access(outms + "/" + name, os.F_OK):
                record_error(outms + "/" + name + " does not exist.", retValue)

        if retValue['success']:
            exp_npoly = 0
            try:
                _tb.open(outms + '/FIELD')
                npoly = _tb.getcell('NUM_POLY', 0)
                if npoly != exp_npoly:
                    record_error('FIELD/NUM_POLY[0], ' + str(npoly) + ' != expected '
                                 + str(exp_npoly), retValue)
                exp_shape = '[2, ' + str(npoly + 1) + ']'
                for dircol in ('PHASE_DIR', 'DELAY_DIR', 'REFERENCE_DIR'):
                    ref = _tb.getcolkeywords(dircol)['MEASINFO']['Ref']
                    if ref != refcode:
                        record_error(dircol + "'s stated frame, " + ref
                                     + ', != expected ' + refcode, retValue)
                    dirshape = _tb.getcolshapestring(dircol)
                    if dirshape[0] != exp_shape:
                        record_error(dircol + "'s shape, " + dirshape
                                     + ', != expected ' + exp_shape)
            except:
                record_error('Error: Cannot get FIELD directions.', retValue)
            else:
                _tb.close()

        self.assertTrue(retValue['success'])

    def test2(self):
        '''Test2: Apply trivial phase center shift, i.e. none.'''
        refcode = 'J2000'
        shutil.rmtree(outms2, ignore_errors=True)

        mystats = self._fixvis_and_get_stats('J2000 18h00m02.3092s -29d59m29.9987s')
        self.assertTrue( mystats['maxposf'] == '18:00:02.309, -29.59.29.999, I, 2.25982e+11Hz' and
                         (mystats['maxpos'] == [64, 64, 0, 0]).all( ) )

    def test3(self):
        '''Test3: Apply positive phase center shift along DEC.'''
        refcode = 'J2000'
        shutil.rmtree(outms2, ignore_errors=True)
        mystats = self._fixvis_and_get_stats('J2000 18h00m02.3092s -29d59m26.9987s')
        self.assertTrue( mystats['maxposf'] == '18:00:02.309, -29.59.29.999, I, 2.25982e+11Hz' and
                         (mystats['maxpos'] == [64,34,0,0]).all() )

    def test4(self):
        '''Test4: Apply negative phase center shift along DEC using offset syntax.'''
        refcode = 'J2000'
        shutil.rmtree(outms2, ignore_errors=True)
        mystats = self._fixvis_and_get_stats('0h -0d0m3s')
        self.assertTrue( mystats['maxposf']=='18:00:02.309, -29.59.29.999, I, 2.25982e+11Hz' and
                         (mystats['maxpos']==[64,94,0,0]).all() )
    def test5(self):
        '''Test5: Apply positive phase center shift along RA.'''
        mystats = self._fixvis_and_get_stats('J2000 18h00m02.5401s -29d59m29.9987s')
        self.assertTrue(mystats['maxposf']=='18:00:02.309, -29.59.29.999, I, 2.25982e+11Hz' and
                        (mystats['maxpos']==[94,64,0,0]).all())

    def test6(self):
        '''Test6: Apply negative phase shift along RA using offset syntax (offset is an angle).'''
        mystats = self._fixvis_and_get_stats('-0d0m3s 0deg')
        self.assertTrue(mystats['maxposf']=='18:00:02.309, -29.59.29.999, I, 2.25982e+11Hz' and
                        (mystats['maxpos']==[34,64,0,0]).all())

    def test7(self):
        '''Test7: Apply negative phase shift along RA in field 1 (using offset syntax, offset is a time), no shift in field 0.'''
        refcode = 'J2000'
        shutil.rmtree(outms2, ignore_errors=True)
        os.system('cp -R ' + inpms2 + ' ' + outms2)

        mystats0 = ''
        mystats1 = ''
        try:
            x = _qa.div(_qa.quantity(-3./3600., 'deg'),
                       _qa.cos(_qa.quantity(31.,'deg')))['value']*24./360.*3600. # (seconds)
            phc = str(x) + 's 0deg'
            self.res = fixvis(vis=outms2, outputvis=outms2, field='1', refcode=refcode,
                              phasecenter=phc, datacolumn='all')
            self.assertTrue(self.res)
            mystats0 = self._get_stats(0, 'testy')
            mystats1 = self._get_stats(1, 'testz')
        except:
            print("*** Unexpected error ***")
            self.assertFalse(True)

        self.assertTrue(mystats0['maxposf']=='18:00:02.309, -29.59.29.999, I, 2.25982e+11Hz' and
                        (mystats0['maxpos']==[64,64,0,0]).all() and
                        mystats1['maxposf']=='18:00:02.333, -30.59.29.999, I, 2.25982e+11Hz' and
                        (mystats1['maxpos']==[34,64,0,0]).all())

def suite():
    return [fixvis_test1]

if __name__ == '__main__':
    unittest.main()

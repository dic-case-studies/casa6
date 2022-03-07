from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import unittest
import math
import numpy as np
import numbers

try:
    from casatools import ctsys, table, msmetadata
    from casatasks import exportuvfits, importuvfits
    _tb = table()
    ctsys_resolve = ctsys.resolve
    _msmd = msmetadata()
    is_CASA6 = True
except ImportError:
    from tasks import *
    from taskinit import *
    import casac
    from __main__ import *
    _tb = tbtool()
    _msmd = msmdtool()
    is_CASA6 = False
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)

    data_root = os.environ.get('CASAPATH').split()[0] + '/casatestdata/'
    def ctsys_resolve(apath):
        return os.path.join(data_root, apath)

datapath = 'unittest/exportuvfits/'

class exportuvfits_test(unittest.TestCase):

    def setUp(self):
        self.datapath = ctsys_resolve(datapath)
        self.msname = ''
        self.reimportms = ''
        self.fitsname = ''
    
    def tearDown(self):
        if os.path.exists(self.msname):
            shutil.rmtree(self.msname)
        if os.path.exists(self.reimportms):
            shutil.rmtree(self.reimportms)
        if os.path.exists(self.fitsname):
            os.remove(self.fitsname)

    def test_export_overwrite(self):
        """CAS-5492: test the overwrite parameter when exporting MSes to uvfits"""
        self.msname = "uvfits_test.ms"
        shutil.copytree(os.path.join(self.datapath, self.msname), self.msname)
        self.fitsname = "CAS-5492.uvfits"
        res = exportuvfits(vis=self.msname, fitsfile=self.fitsname)
        if is_CASA6:
            # Not sure why all of a sudden CASA6 is returning None for tasks
            self.assertTrue(res == None, "Failed exportuvfits")
        else:
            self.assertTrue(res, "Failed exportuvfits")
        # fail because overwrite=False.
        # CASA 6 throws an exception, CASA 5 returns False
        if is_CASA6 or casa_stack_rethrow:
            self.assertRaises(
                Exception, exportuvfits, vis=self.msname, fitsfile=self.fitsname, overwrite=False,
                msg="exportuvfits succeeded but should have failed because "
                + "overwrite=False"
            )
        else:
            self.assertFalse(
                exportuvfits(vis=self.msname, fitsfile=self.fitsname, overwrite=False),
                "exportuvfits succeeded but should have failed because "
                + "overwrite=False"
            )
        # succeed because overwrite=True
        res = exportuvfits(vis=self.msname, fitsfile=self.fitsname, overwrite=True)
        if is_CASA6:
            self.assertTrue(
                res == None,
                "exportuvfits failed but should have succeeded because "
                + "overwrite=True"
            )
        else:
            self.assertTrue(
                res,
                "exportuvfits failed but should have succeeded because "
                + "overwrite=True"
            )
        

    def test_no_rest_freqs(self):
        """CAS-11514: test exporting an MS with no rest frequencies in the SOURCE table"""
        self.msname = "rest_freq_test.ms"
        shutil.copytree(os.path.join(self.datapath, self.msname), self.msname)
        self.fitsname = "no_rest_freqs.uvfits"
        res = exportuvfits(vis=self.msname, fitsfile=self.fitsname)
        if is_CASA6:
            self.assertTrue(res == None, "Failed exportuvfits with no rest freqs")
        else:
            self.assertTrue(res, "Failed exportuvfits with no rest freqs")
        # import and check the rest freqs
        # importuvfits doesn't return anything, so we cannot test the
        # return value for success
        self.reimportms = "imported_no_restfreqs.ms"
        importuvfits(fitsfile=self.fitsname, vis=self.reimportms)
        _msmd.open(self.reimportms)
        restfreqs = _msmd.restfreqs()
        _msmd.done()
        expec = {
            '0': {
                'type': 'frequency', 'm0': {'value': 0.0, 'unit': 'Hz'},
                'refer': 'LSRK'
            }
        }
        self.assertEqual(
            restfreqs, expec, "Got wrong restfreqs from re-imported dataset"
        )

    def test_no_source_table(self):
        """CAS-11514: test exporting an MS with no rest frequencies in the SOURCE table"""
        self.msname = "no_source_table.ms"
        shutil.copytree(os.path.join(self.datapath, self.msname), self.msname)
        self.fitsname = "no_source_table.uvfits"
        res = exportuvfits(vis=self.msname, fitsfile=self.fitsname)
        if is_CASA6:
            self.assertTrue(res == None, "Failed exportuvfits with no SOURCE table")
        else:
            self.assertTrue(res, "Failed exportuvfits with no SOURCE table")
            
def suite():
    return [exportuvfits_test]        

if __name__ == '__main__':
    unittest.main()
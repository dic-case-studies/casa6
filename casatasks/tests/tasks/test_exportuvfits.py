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
    data_root = os.environ.get('CASAPATH').split()[0] + '/data'
    def ctsys_resolve(apath):
        return os.path.join(data_root, apath)

datapath = 'regression/unittest/uvfits'

class exportuvfits_test(unittest.TestCase):

    def setUp(self):
        self.datapath = ctsys_resolve(datapath)
    
    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)
        # make sure directory is clean as per verification test requirement
        cwd = os.getcwd()
        for filename in os.listdir(cwd):
            file_path = os.path.join(cwd, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    # CASA 5 tests need this directory
                    if filename != 'xml':
                        shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

    def test_export_overwrite(self):
        """CAS-5492: test the overwrite parameter when exporting MSes to uvfits"""
        msname = "uvfits_test.ms" 
        shutil.copytree(os.path.join(self.datapath, msname), msname)
        fitsname = "CAS-5492.uvfits"
        self.assertTrue(
            exportuvfits(vis=msname, fitsfile=fitsname),
            "Failed exportuvfits"
        ) 
        # fail because overwrite=False.
        # CASA 6 throws an exception, CASA 5 returns False
        if is_CASA6:
            self.assertRaises(
                Exception, exportuvfits, vis=msname, fitsfile=fitsname, overwrite=False,
                msg="exportuvfits succeeded but should have failed because "
                + "overwrite=False"
            )
        else:
            self.assertFalse(
                exportuvfits(vis=msname, fitsfile=fitsname, overwrite=False),
                "exportuvfits succeeded but should have failed because "
                + "overwrite=False"
            )
        # succeed because overwrite=True
        self.assertTrue(
            exportuvfits(vis=msname, fitsfile=fitsname, overwrite=True),
            "exportuvfits failed but should have succeeded because "
            + "overwrite=True"
        )

    def test_no_rest_freqs(self):
        """CAS-11514: test exporting an MS with no rest frequencies in the SOURCE table"""
        msname = "rest_freq_test.ms"
        shutil.copytree(os.path.join(self.datapath, msname), msname)
        fitsname = "no_rest_freqs.uvfits"
        self.assertTrue(
            exportuvfits(vis=msname, fitsfile=fitsname),
            "Failed exportuvfits with no rest freqs"
        ) 
        # import and check the rest freqs
        # importuvfits doesn't return anything, so we cannot test the
        # return value for success
        msname = "imported_no_restfreqs.ms"
        importuvfits(fitsfile=fitsname, vis=msname)
        _msmd.open(msname)
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
        msname = "no_source_table.ms"
        shutil.copytree(os.path.join(self.datapath, msname), msname)
        fitsname = "no_source_table.uvfits"
        self.assertTrue(
            exportuvfits(vis=msname, fitsfile=fitsname),
            "Failed exportuvfits with no SOURCE table"
        )
            
def suite():
    return [exportuvfits_test]        

if __name__ == '__main__':
    unittest.main()

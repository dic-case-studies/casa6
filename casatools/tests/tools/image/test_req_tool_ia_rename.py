########################################################################
# test_req_tool_ia_transpose.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
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
# CAS-12700
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_imtrans/about
#
#
##########################################################################

from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import unittest

try:
    from casatools import ctsys, image, table
    ctsys_resolve = ctsys.resolve
    myia = image()
    _tb = table()
    is_CASA6 = True
except ImportError:
    from tasks import *
    from taskinit import *
    import casac
    from __main__ import *
    image = iatool
    myia = iatool()
    _tb = tbtool()
    is_CASA6 = False
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        data_root = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'
    else:
        data_root = os.environ.get('CASAPATH').split()[0] + '/casa-data-req'
    def ctsys_resolve(apath):
        return os.path.join(data_root, apath)

class ia_rename_test(unittest.TestCase):
    
    def setUp(self):
        pass
    
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

    def test_rename(self):
        """verify history writing"""
        myia.fromshape("zz",[20, 20])
        newname = "xx.im"
        self.assertTrue(myia.rename(newname), "rename unsuccessful")
        got = myia.name(strippath=True)
        self.assertTrue(
            got == newname,
            "wrong name " + got + " should be " + newname
        )
        myia.done()

    def test_overwrite(self):
        newname = "kfe.im"
        myia.fromshape(newname, [20, 20])
        myia.done()
        name = "jfjd.im"
        myia.fromshape(name, [5, 5])
        try:
            myia.rename(newname, overwrite=False)
            thrown = False
        except:
            thrown = True
        self.assertTrue(thrown, "overwrite=False exception not thrown")
        myia.open(name)
        res = myia.rename(newname, overwrite=True)
        print("res", res)
        self.assertTrue(myia.rename(newname, overwrite=True), "overwrite=True unsuccessful")
        self.assertTrue(myia.name(strippath=True) == newname, "wrong name")
        myia.done()

    def test_history(self):
        """verify history writing"""
        myia.fromshape("zz",[20, 20])
        myia.rename("zy.im")
        msgs = myia.history()
        myia.done()
        self.assertTrue("ia.rename" in msgs[-3], "wrong history")
        self.assertTrue("ia.rename" in msgs[-2], "wrong history")
        
def suite():
    return [ia_rename_test]

if __name__ == '__main__':
    unittest.main()


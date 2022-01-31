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

import os
import shutil
import unittest


from casatools import ctsys, image, table
ctsys_resolve = ctsys.resolve
myia = image()
_tb = table()



class ia_rename_test(unittest.TestCase):
    
    def setUp(self):
        self.newname = ''
    
    def tearDown(self):
        if self.newname:
            if os.path.isfile(self.newname):
                os.unlink(self.newname)
            else:
                shutil.rmtree(self.newname)
        self.assertTrue(len(_tb.showcache()) == 0)

    def test_rename(self):
        """verify history writing"""
        myia.fromshape("zz",[20, 20])
        self.newname = "xx.im"
        self.assertTrue(myia.rename(self.newname), "rename unsuccessful")
        got = myia.name(strippath=True)
        self.assertTrue(
            got == self.newname,
            "wrong name " + got + " should be " + self.newname
        )
        myia.done()

    def test_overwrite(self):
        self.newname = "kfe.im"
        myia.fromshape(self.newname, [20, 20])
        myia.done()
        name = "jfjd.im"
        myia.fromshape(name, [5, 5])
        try:
            myia.rename(self.newname, overwrite=False)
            thrown = False
        except:
            thrown = True
        self.assertTrue(thrown, "overwrite=False exception not thrown")
        myia.open(name)
        res = myia.rename(self.newname, overwrite=True)
        print("res", res)
        self.assertTrue(myia.rename(self.newname, overwrite=True), "overwrite=True unsuccessful")
        self.assertTrue(myia.name(strippath=True) == self.newname, "wrong name")
        myia.done()

    def test_history(self):
        """verify history writing"""
        myia.fromshape("zz",[20, 20])
        self.newname = "zy.im"
        myia.rename(self.newname)
        msgs = myia.history()
        myia.done()
        self.assertTrue("ia.rename" in msgs[-3], "wrong history")
        self.assertTrue("ia.rename" in msgs[-2], "wrong history")


if __name__ == '__main__':
    unittest.main()


##########################################################################
# test_componentlist.py
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.componentlist.html
#
#
##########################################################################

import unittest

from casatools import componentlist as cltool

class componentlist_test(unittest.TestCase):
    
    def setUp(self):
        return

    def tearDown(self):
        return
 
    def test_summarize(self):
        """Test the cl.summarize() method"""
        # make me a list
        mycl = cltool()
        mycl.addcomponent(
            [1,0,0,0],'Jy','Stokes',['J2000', '10:30:00.00', '-20.00.00.0'],
            'gaussian','4arcsec','2arcsec','30deg'
        )
        self.assertTrue(mycl.summarize(0))
        self.assertRaises(Exception, mycl.summarize, which=1)

    def test_getfluxerror(self):
        """Test cl.getfluxerror()"""
        mycl = cltool()
        mycl.addcomponent(
            [1,0,0,0],'Jy','Stokes',['J2000', '10:30:00.00', '-20.00.00.0'],
            'gaussian','4arcsec','2arcsec','30deg'
        )
        ferror = [0.2, 0.3, 0, 0]
        mycl.setflux(0, value=mycl.getfluxvalue(0), error=ferror)
        got = mycl.getfluxerror(0)
        self.assertTrue((got == ferror).all())

if __name__ == '__main__':
    unittest.main()

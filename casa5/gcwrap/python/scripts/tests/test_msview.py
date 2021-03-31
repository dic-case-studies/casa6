##########################################################################
# test_msview.py
#
# Copyright (C) 2017
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
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be adressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
# <author>
# Darrell Schiebel
# </author>
#
# <summary>
# Minimal test of CASA viewer
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
# </ul>
# </prerequisite>
#
# <etymology>
# </etymology>
#
# <synopsis>
# </synopsis> 
#
# <example>
# </example>
#
# <motivation>
# Because the viewer has no way to specify specific dimensions and because the size
# of the bitmaps that it produces are dependent on the dimensions of the viewer, it
# is difficult to create good tests to verify that the viewer works properly. As a
# result, this test only verifies that it produces the expected file.
# </motivation>
#
###########################################################################
from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import os.path
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys
    from casaviewer import msview
    ctsys_resolve = ctsys.resolve
else:
    from taskinit import *
    from tasks import msview

    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ.get('CASAPATH').split()[0],'data')
        return os.path.join(dataPath,apath)

class msview_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        """copying data"""
        self.testms = "tms"+str(os.getpid())+".ms"
        self.outfiles = { }
        for t in ['jpg', 'pdf', 'eps', 'ps', 'png', 'xbm', 'xpm', 'ppm']:
            self.outfiles[t] = "tms"+str(os.getpid())+"."+t
        os.system('cp -R %s %s' % (ctsys_resolve('regression/fits-import-export/input/test.ms'),self.testms))

    @classmethod
    def tearDownClass(self):
        """removing test images"""
        if os.path.exists(self.testms):
            shutil.rmtree(self.testms)
        for outfileType in self.outfiles:
            thisOutfile = self.outfiles[outfileType]
            if os.path.exists(thisOutfile):
                os.system('rm -rf ' + thisOutfile)
        
    def setUp(self):
        pass

    def test_xbm(self):
        """Test production of Xbm file"""
        msview(self.testms,outfile=self.outfiles['xbm'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['xbm']),"viewer failed to produce an Xbm file")

    def test_jpg(self):
        """Test production of JPEG file"""
        msview(self.testms,outfile=self.outfiles['jpg'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['jpg']),"viewer failed to produce an JPEG file")

    def test_pdf(self):
        """Test production of PDF file"""
        msview(self.testms,outfile=self.outfiles['pdf'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['pdf']),"viewer failed to produce an PDF file")

    def test_eps(self):
        """Test production of EPS file"""
        msview(self.testms,outfile=self.outfiles['eps'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['eps']),"viewer failed to produce an EPS file")

    def test_ps(self):
        """Test production of PS file"""
        msview(self.testms,outfile=self.outfiles['ps'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['ps']),"viewer failed to produce an PS file")

    def test_xpm(self):
        """Test production of XPM file"""
        msview(self.testms,outfile=self.outfiles['xpm'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['xpm']),"viewer failed to produce an XPM file")

    def test_ppm(self):
        """Test production of PPM file"""
        msview(self.testms,outfile=self.outfiles['ppm'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['ppm']),"viewer failed to produce an PPM file")

    def test_png(self):
        """Test production of PNG file"""
        msview(self.testms,outfile=self.outfiles['png'],gui=False)
        self.assertTrue(os.path.isfile(self.outfiles['png']),"viewer failed to produce an PNG file")

def suite():
    return [msview_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

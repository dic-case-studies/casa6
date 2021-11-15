##########################################################################
#
# Copyright (C) 2019 ESO (in the framework of the ALMA collaboration)
# Copyright (C) 2019 Associated Universities, Inc. Washington DC, USA.
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
#
#
##########################################################################
CASA6=False
try:
    import casatools
    from casatasks import importfits, imhead
    tb = casatools.table()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

import os
import shutil
import unittest
import numpy as np

reg_unittest_datap = 'unittest/importfits/'
datapath = casatools.ctsys.resolve(reg_unittest_datap)
fitsdata = datapath + 'two_gaussian_model.fits'
multihdu = datapath + 'vogtstar_awt.fits'
axesdata = datapath + 'test_image2dconvolver.fits'
outpath = 'fitstestout.im'
outpath2 = 'fitstestout2.im'


class importfits_test(unittest.TestCase):

    def tearDown(self):
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        if os.path.exists(outpath2):
            shutil.rmtree(outpath2)

    def test_createsCasaImage(self):
        ''' Check that the task creates a casa image from a fits image '''

        importfits(fitsimage=fitsdata, imagename=outpath)

        self.assertTrue(os.path.exists(outpath))

    def test_whichrep(self):
        ''' Check that the whichrep parameter can change coordinate representation '''

        # Need a dataset with multiple reps
        importfits(fitsimage=fitsdata, imagename=outpath, whichrep=1)

        self.assertTrue(os.path.exists(outpath))

    def test_whichHDU(self):
        ''' Check that specific header data can be selected if the fits contains multiple images '''

        importfits(fitsimage=multihdu, imagename=outpath, whichhdu=1)

        self.assertTrue(os.path.exists(outpath))

    def test_defaultAxes(self):
        """Check that default coordinate axes are added where they are missing"""

        importfits(fitsimage=axesdata, imagename=outpath, defaultaxes=True, defaultaxesvalues=['19h30m00', '-02d30m00', '88.5GHz', 'Q'])
        importfits(fitsimage=axesdata, imagename=outpath2, defaultaxes=False)

        # Get the two different images
        tb.open(outpath)
        withaxes = tb.getcol('map')
        tb.close()

        tb.open(outpath2)
        withoutaxes = tb.getcol('map')
        tb.close()

        self.assertFalse(np.array_equal(withaxes, withoutaxes))

    def test_overwrite(self):
        """Check that the task can only overwrite a table if overwrite=True"""

        importfits(fitsimage=fitsdata, imagename=outpath)
        # Run the task again and it should fail
        try:
            importfits(fitsimage=fitsdata, imagename=outpath)
        except RuntimeError:
            noOverwrite = False
        # This should return false
        self.assertFalse(noOverwrite)

        # Running again with overwrite should give no return value
        withOverwrite = importfits(fitsimage=fitsdata, imagename=outpath, overwrite=True)
        self.assertTrue(withOverwrite == None)

    def test_beam(self):
        """Provide values to be used with the synthesized beam"""

        importfits(fitsimage=fitsdata, imagename=outpath, beam=['0.35arcsec', '0.24arcsec', '25deg'])
        summary = imhead(outpath, mode='summary')
        beamresult = summary['restoringbeam']

        # Check that the image has a beam with the specified properties
        # Major Beam
        self.assertTrue(beamresult['major']['unit'] == 'arcsec')
        self.assertTrue(beamresult['major']['value'] == 0.35)
        # Minor Beam
        self.assertTrue(beamresult['minor']['unit'] == 'arcsec')
        self.assertTrue(beamresult['minor']['value'] == 0.24)
        # Position Angle
        self.assertTrue(beamresult['positionangle']['unit'] == 'deg')
        self.assertTrue(beamresult['positionangle']['value'] == 25.0)

def suite():
    return [importfits_test]


if __name__ == '__main__':
    unittest.main()

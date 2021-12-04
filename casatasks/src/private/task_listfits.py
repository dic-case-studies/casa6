from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ms
    _ms = ms( )
else:
    from taskinit import *
    # not a local tool
    _ms = ms

def listfits(fitsfile=None):
    """

    """
    #Python script

    if ((type(fitsfile)==str) & (os.path.exists(fitsfile))):
        _ms.listfits(fitsfile);
    else:
        raise ValueError('fits file not found - please verify the name')

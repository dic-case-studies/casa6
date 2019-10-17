from __future__ import absolute_import
from __future__ import print_function
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
    try:
        if ((type(fitsfile)==str) & (os.path.exists(fitsfile))):
            _ms.listfits(fitsfile);    
        else:
            raise Exception('fits file not found - please verify the name')
    except Exception as instance:
        print('*** Error *** %s' % instance)

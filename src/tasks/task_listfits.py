import os
from casatools import ms

_ms = ms( )

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

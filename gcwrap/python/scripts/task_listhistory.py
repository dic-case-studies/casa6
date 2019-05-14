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

def listhistory(vis=None):
    """List the processing history of a dataset:
    The list of all task processing steps will be
    given in the logger.

    Keyword arguments:
    vis -- Name of input visibility file (MS)
            default: none; example: vis='ngc5921.ms'

    """
    #Python script
    try:
        if ((type(vis)==str) & (os.path.exists(vis))):
            _ms.open(vis)
        else:
            raise Exception('Visibility data set not found - please verify the name')
        _ms.listhistory()
        _ms.close()
    except Exception as instance:
        print('*** Error *** %s' % instance)

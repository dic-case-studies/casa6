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
        _ms.open(vis)
        _ms.listhistory()
    finally:
        _ms.close()

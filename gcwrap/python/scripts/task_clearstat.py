from __future__ import absolute_import
from __future__ import print_function
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import table
    from casatasks import casalog
    _tb = table( )
else:
    from taskinit import *
    _tb = tb

def clearstat():
    """Clear all read/write locks on tables. This can be used if a task has
       indicated that it is trying to get a lock on a file.

    """
    try:
        _tb.clearlocks( )
    except Exception as instance:
        if is_CASA6:
            casalog.post('*** Error *** %s' % instance)
        else:
            print('*** Error ***',instance)


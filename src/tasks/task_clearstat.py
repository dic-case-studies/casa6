import os
from CASAtools import table
from CASAtasks import casalog

_tb = table( )

def clearstat():
    """Clear all read/write locks on tables. This can be used if a task has
       indicated that it is trying to get a lock on a file.

    """
    try:
        _tb.clearlocks( )
    except Exception as instance:
        casalog.post('*** Error *** %s' % instance)

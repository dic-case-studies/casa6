# Task listvis

from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ms
    from casatasks import casalog
else:
    from taskinit import *
    ms = casac.ms

def listvis(vis, options, datacolumn, field,spw, selectdata, antenna, timerange,
            correlation, scan, feed, array, observation, uvrange, average,
            showflags, pagerows, listfile):
    """List visibilities on terminal."""
        
    casalog.origin('listvis')
    myms = ms()
    
    isInteractive=False;
    
    try:
        if ((type(vis)==str) & (os.path.exists(vis))):
            myms.open(vis)
        else:
            raise Exception('Visibility data set not found - please verify the name')
                
        myms.lister(options, datacolumn, field, spw, antenna, timerange,
                    correlation, scan, feed, array, str(observation), uvrange,
                    average, showflags, "", pagerows, listfile)

    finally:
        myms.close()

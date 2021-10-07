from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ms
    from casatasks import casalog
else:
    from taskinit import *
    ms = mstool
    
def listobs(
    vis, selectdata, spw, field, antenna, uvrange,
    timerange, correlation, scan, intent, feed,
    array, observation, verbose, listfile,
    listunfl, cachesize, overwrite
):
    
    """List data set summary in the logger:

       Lists following properties of a measurement set:
       scan list, field list, spectral window list with
       correlators, antenna locations, ms table information.

       Keyword arguments:
       vis -- Name of input visibility file
               default: none. example: vis='ngc5921.ms'
       selectdata -- select data from the MS
       verbose -- level of detail
             verbose=True: (default); scan and antenna lists
             verbose=False: less information
       listfile -- save the output in a file
             default: none. Example: listfile="mylist.txt"

       """

    casalog.origin('listobs')

    try:
        myms = ms()
        myms.open(thems=vis, check=True)
                
        sel = {}
        if (selectdata):
            sel['spw'] = spw
            sel['time'] = timerange
            sel['field'] = field
            sel['baseline'] = antenna
            sel['scan'] = scan
            sel['scanintent'] = intent
            sel['polarization'] = correlation
            sel['uvdist'] = uvrange
            sel['observation'] = str(observation)
            sel['array'] = array
            sel['feed'] = feed

        # Select the data. Only-parse is set to false.
        myms.msselect(sel, False)
        obs_dict = myms.summary(
            verbose=verbose, listfile=listfile, listunfl=listunfl,
            cachesize=cachesize, overwrite=overwrite, wantreturn=True
        )
        return obs_dict
    finally:
        myms.close()

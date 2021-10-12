from __future__ import absolute_import
import os

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import ms
    from casatasks import casalog
    from .mstools import write_history
else:
    from taskinit import *
    from mstools import write_history

    ms = mstool

def importuvfits(fitsfile, vis, antnamescheme=None):
    """

    Convert a UVFITS file to a CASA visibility data set (MS):

    Keyword arguments:
    fitsfile -- Name of input UV FITS file
        default = none; example='3C273XC1.fits'
    vis -- Name of output visibility file (MS)
        default = none; example: vis='3C273XC1.ms'
    antnamescheme -- Naming scheme for VLA/JVLA/CARMA antennas
        default = old;
        old: Antenna name is a number, '04'
             This option exists for backwards compatibility
             but can lead to ambiguous results when antenna
             indices are used for data selection.
        new: Antenna name is not a number, 'VA04' or 'EA04'
             With this scheme, data selection via
             antenna names and indices is non-ambiguous.
        default = false; do not run asychronously


    """
    myms = ms()
    try:
        casalog.origin('importuvfits')
        casalog.post("")
        myms.fromfits(vis, fitsfile, antnamescheme=antnamescheme)
        myms.close()

        # Write the args to HISTORY.
        try:
            param_names = importuvfits.__code__.co_varnames[:importuvfits.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            write_history(
                myms, vis, 'importuvfits', param_names, 
                param_vals, casalog
            )
        except Exception:
            casalog.post("Failed to updated HISTORY table", 'WARN')

    finally:
        if (myms):
            #myms.close()
            del myms 



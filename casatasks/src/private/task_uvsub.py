#
# This file was generated using xslt from its XML file
#
# Copyright 2007, Associated Universities Inc., Washington DC
#
from __future__ import absolute_import
import os

# get is_python3 and is_CASA6
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import ms as mstool
    from casatasks import casalog
    from .mstools import write_history

    _ms = mstool()
else:
    from taskinit import mstool, casalog
    from mstools import write_history

    # uses the global ms tool
    _ms = mstool()

def uvsub(vis=None,reverse=False):

    """Subtract model from the corrected visibility data
    
        uvsub(vis='ngc5921.ms', reverse=false)
        
        This function subtracts model visibility data from corrected visibility
        data leaving the residuals in the corrected data column.  If the
        parameter 'reverse' is set true, the process is reversed.
        
        Keyword arguments:
        vis -- Name of input visibility file (MS)
                default: none; example: vis='ngc5921.ms'
        reverse -- Reverse the operation (add rather than subtract)
                default: false; example: reverse=true
        
        uvsub(vis='ngc5921.ms', reverse=false)
    
    """

    #Python script
    #

    try:
        casalog.origin('uvsub')
        if ((type(vis)==str) & (os.path.exists(vis))):
            _ms.open(thems=vis,nomodify=False)
        else:
            raise ValueError('Visibility data set not found - please verify the name')
            return
        _ms.uvsub(reverse)

        # Write history to MS
        try:
            param_names = uvsub.__code__.co_varnames[:uvsub.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            write_history(mstool(), vis, 'uvsub', param_names,
                          param_vals, casalog)
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                         'WARN')

    finally:
        _ms.close()

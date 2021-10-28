from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
    from casatools import calibrater
    _cb = calibrater()
else:
    from taskinit import *
    _cb = cbtool( )

def rerefant(vis,tablein,caltable,refantmode,refant):
    """ Smooth calibration solution(s) derived from one or more sources:

    Keyword arguments:
    vis -- Name of input visibility file (MS)
            default: none; example: vis='ngc5921.ms'
    tablein -- Input calibration table (any type)
            default: none; example: tablein='ngc5921.gcal'
    caltable -- Output calibration table (re-refant-ed)
            default: none
    refantmode -- Reference antenna application mode (TBD)
    refant -- Reference antenna name(s)
    """

    #Python script
    try:
        casalog.origin('rerefant')
        if ((type(vis)==str) & (os.path.exists(vis))):
            _cb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
                   
        else:
            raise Exception('Visibility data set not found - please verify the name')

        _cb.rerefant(tablein=tablein,tableout=caltable,refantmode=refantmode,refant=refant);

    finally:
        _cb.close()


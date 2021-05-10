from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import calibrater
    from casatasks import casalog

    _cb = calibrater()
else:
    from taskinit import *

    _cb = cbtool( )

def smoothcal(vis,tablein,caltable,field,smoothtype,smoothtime):
        """ Smooth calibration solution(s) derived from one or more sources:

        Keyword arguments:
        vis -- Name of input visibility file (MS)
                default: none; example: vis='ngc5921.ms'
        tablein -- Input calibration table (any type)
                default: none; example: tablein='ngc5921.gcal'
        caltable -- Output calibration table (smoothed)
                default: none; example: caltable='ngc5921_smooth.gcal'
        field -- subset of sources to select
                default: '' means all; example: field='0319_415_1 3C286'
        smoothtype -- The smoothing filter to be used
                default: 'median'; example: smoothtype='mean'
                Options: 'median','mean'
        smoothtime -- Smoothing filter time (sec)
                default: 60.0; example: smoothtime=10.
        """

        #Python script
        try:
                casalog.origin('smoothcal')
                if ((type(vis)==str) & (os.path.exists(vis))):
                        _cb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
                               
                else:
                        raise ValueError('Visibility data set not found - please verify the name')

                _cb.smooth(tablein=tablein,tableout=caltable,field=field,smoothtype=smoothtype,smoothtime=smoothtime)

        finally:
                _cb.close()

from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
       from casatasks import casalog
       from casatools import calibrater
else:
       from taskinit import *
       calibrater = cbtool

def listcal(vis=None,caltable=None,field=None,antenna=None,spw=None,
            listfile=None,pagerows=None):
       """List calibration solutions (amp and phase)."""

       casalog.origin('listcal')

       #Python script

       try:
              mycb = calibrater()
              if ((type(vis)==str) & (os.path.exists(vis))):
                     mycb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
              else:
                     raise Exception('Visibility data set not found - please verify the name')
              mycb.listcal(caltable=caltable,field=field,antenna=antenna,spw=spw,
                         listfile=listfile,pagerows=pagerows)

       finally:
              mycb.close()

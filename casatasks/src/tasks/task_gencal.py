from __future__ import absolute_import
import sys
import os
import warnings

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import calibrater
    from casatasks import casalog
    from . import correct_ant_posns as getantposns 

    _cb = calibrater( )
else:
    from taskinit import *
    import correct_ant_posns as getantposns 
    (_cb,) = gentools(['cb'])

def gencal(vis=None,caltable=None,caltype=None,infile=None,
        spw=None,antenna=None,pol=None,
        parameter=None,uniform=None):

    """ Externally specify calibration solutions af various types
    """
    
    #Python script
    try:

        if ((type(vis)==str) & (os.path.exists(vis))):
            # don't need scr col for this
            _cb.open(filename=vis,compress=False,addcorr=False,addmodel=False)  
        else:
            raise ValueError('Visibility data set not found - please verify the name')

        if (caltable==''):
            raise ValueError('A caltable name must be specified')

        if caltype=='tecim' and not (type(infile)==str and os.path.exists(infile)):
            raise ValueError('An existing tec map must be specified in infile')

        # call a Python function to retreive ant position offsets automatically (currently EVLA only)
        if (caltype=='antpos' and antenna==''):
            casalog.post(" Determine antenna position offsets from the baseline correction database")
            # correct_ant_posns returns a list , [return_code, antennas, offsets]
            antenna_offsets=getantposns.correct_ant_posns(vis,False)
            if ((len(antenna_offsets)==3) and
                    (int(antenna_offsets[0])==0) and
                    (len(antenna_offsets[1])>0) ) :
                antenna = antenna_offsets[1]
                parameter = antenna_offsets[2] 
             else:
                #raise Exception, 'No offsets found. No caltable created.'
                warnings.simplefilter('error',UserWarning)
                warnings.warn('No offsets found. No caltable created.')

        _cb.specifycal(caltable=caltable,time="",spw=spw,antenna=antenna,pol=pol,
                caltype=caltype,parameter=parameter,infile=infile,
                uniform=uniform)

        #_cb.close()
    
    except UserWarning as instance:
        casalog.post('*** UserWarning *** %s' % instance, 'WARN')

    finally:
        _cb.close()

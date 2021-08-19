from __future__ import absolute_import

import os
import sys
import warnings

import numpy as np

from casatasks.private.casa_transition import is_CASA6

if is_CASA6:
    from casatasks import casalog
    from casatools import calibrater
    from . import correct_ant_posns as getantposns
    from .jyperk import gen_factor_via_web_api, JyPerKReader4File

    _cb = calibrater()
else:
    import correct_ant_posns as getantposns
    from taskinit import *

    (_cb,) = gentools(['cb'])


def gencal(vis=None, caltable=None, caltype=None, infile='None',
           endpoint='asdm', timeout=180, retry=3, retry_wait_time=5,
           spw=None, antenna=None, pol=None,
           parameter=None, uniform=None):
    """Externally specify calibration solutions of various types.

    Arguments:
        vis {str} -- The file path stored the visibility data.
        caltable {str} -- A file name which store the caltable.
        caltype {str} -- The calibration type.
        infile (str) -- Specifies the name of the file to read.
        subparameter of caltype='jyperk:
            endpoint (str) -- The endpoint of the Jy/K DB Web API to access.
                options are 'asdm' (default), 'model-fit', 'interpolation'.
            timeout {int} --- Maximum waiting time [sec] for the Web API access,
                defaults to 180 sec.
            retry {int} -- Number of retry when the Web API access fails,
                defaults to 3 times.
            retry_wait_time {int} -- Waiting time [sec] until next query 
                when the Web API access fails, defaults to 5 sec.
        spw {str} -- The spectral windows.
        antenna {str} --
        pol {str} --
        parameter {} --
        uniform {} --
    """

    # validate arguments
    if (caltable==''):
        raise ValueError('A caltable name must be specified')

    if caltype=='tecim' and not (type(infile)==str and os.path.exists(infile)):
        raise ValueError('An existing tec map must be specified in infile')

    if caltype == 'jyperk' and not endpoint in ['asdm', 'interpolation', 'model-fit']:
        raise ValueError('When the caltype is jyperk, endpoint must be one of asdm, interpolation or model-fit')

    #Python script
    try:
        if ((type(vis)==str) & (os.path.exists(vis))):
            # don't need scr col for this
            _cb.open(filename=vis,compress=False,addcorr=False,addmodel=False)  
        else:
            raise ValueError('Visibility data set not found - please verify the name')


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
       
    except UserWarning as instance:
        casalog.post('*** UserWarning *** %s' % instance, 'WARN')

    finally:
        _cb.close()

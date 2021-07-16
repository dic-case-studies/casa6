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

    # check arguments
    if (caltable == ''):
        raise ValueError('A caltable name must be specified')

    if caltype == 'tecim' and not (type(infile) == str and os.path.exists(infile)):
        raise ValueError('An existing tec map must be specified in infile')

    if caltype == 'jyperk' and not endpoint in ['asdm', 'interpolation', 'model-fit']:
        raise ValueError('When the caltype is jyperk, endpoint must be one of asdm, interpolation or model-fit')

    if not type(vis) == str:
        raise ValueError('Visibility data set not found - vis should be str')

    if not os.path.exists(vis):
        raise ValueError('Visibility data set not found - please verify the name')

    try:
        _cb.open(filename=vis, compress=False, addcorr=False, addmodel=False)

        # call a Python function to retreive ant position offsets automatically (currently EVLA only)
        if (caltype == 'antpos' and antenna == ''):
            antenna, parameter = __complete_antpos(vis)

        if caltype == 'jyperk':
            for selection, param in __gen_specifycal_input(vis=vis, endpoint=endpoint, infile=infile,
                                                           timeout=timeout, retry=retry,
                                                           retry_wait_time=retry_wait_time):
                _cb.specifycal(caltable='amp', time='', spw=selection['spw'],
                               antenna=selection['antenna'], pol=selection['pol'],
                               parameter=param, infile='', uniform=uniform)
        else:
            _cb.specifycal(caltable=caltable, time='', spw=spw, antenna=antenna, pol=pol,
                           caltype='antpos', parameter=parameter, infile=infile,
                           uniform=uniform)

    except UserWarning as instance:
        casalog.post('*** UserWarning *** %s' % instance, 'WARN')

    finally:
        _cb.close()

def __complete_antpos(vis):
    casalog.post(" Determine antenna position offsets from the baseline correction database")
    # correct_ant_posns returns a list , [return_code, antennas, offsets]
    antenna_offsets = getantposns.correct_ant_posns(vis, False)
    if ((len(antenna_offsets) == 3) and
            (int(antenna_offsets[0]) == 0) and
            (len(antenna_offsets[1]) > 0)):
        antenna = antenna_offsets[1]
        parameter = antenna_offsets[2]
        return antenna, parameter
    else:
        # raise Exception, 'No offsets found. No caltable created.'
        warnings.simplefilter('error', UserWarning)
        warnings.warn('No offsets found. No caltable created.')

def __gen_specifycal_input(vis=None, endpoint='asdm', infile=None,
                           timeout=180, retry=3, retry_wait_time=5):
    if not infile is None:
        f = JyPerKReader4File(infile)
        factors = f.get()

    elif infile is None:
        factors = gen_factor_via_web_api(vis, endpoint=endpoint, 
                                         timeout=timeout, retry=retry, 
                                         retry_wait_time=retry_wait_time)

    for factor in factors:
        selection = {}
        selection['antenna'] = factor[1]
        selection['spw'] = factor[2]
        selection['pol'] = factor[3]
        
        yield selection, 1/np.sqrt(float(factor[4]))
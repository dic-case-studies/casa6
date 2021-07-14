from __future__ import absolute_import

import os
import sys
import warnings

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


def gencal(vis=None, caltable=None, caltype=None, infile=None,
           spw=None, antenna=None, pol=None,
           parameter=None, uniform=None,
           timeout=180, retry=3, retry_wait_time=5):

    """Externally specify calibration solutions af various types.

    Arguments:
        vis {str} -- The file path of vis.
        timeout {int} --- Maximum waiting time when accessing the web API. Second.
        retry {int} -- Number of times to retry when the web API access fails.
        retry_wait_time {int} -- The waiting time when the web request fails. Second.
    """

    # check arguments
    if (caltable == ''):
        raise ValueError('A caltable name must be specified')

    if caltype == 'tecim' and not (type(infile) == str and os.path.exists(infile)):
        raise ValueError('An existing tec map must be specified in infile')

    # Python script
    try:
        if ((type(vis) == str) & (os.path.exists(vis))):
            # don't need scr col for this
            _cb.open(filename=vis, compress=False, addcorr=False, addmodel=False)
        else:
            raise ValueError('Visibility data set not found - please verify the name')

        # call a Python function to retreive ant position offsets automatically (currently EVLA only)
        if (caltype == 'antpos' and antenna == ''):
            casalog.post(" Determine antenna position offsets from the baseline correction database")
            # correct_ant_posns returns a list , [return_code, antennas, offsets]
            antenna_offsets = getantposns.correct_ant_posns(vis, False)
            if ((len(antenna_offsets) == 3) and
                    (int(antenna_offsets[0]) == 0) and
                    (len(antenna_offsets[1]) > 0)):
                antenna = antenna_offsets[1]
                parameter = antenna_offsets[2]
            else:
                # raise Exception, 'No offsets found. No caltable created.'
                warnings.simplefilter('error', UserWarning)
                warnings.warn('No offsets found. No caltable created.')

        if caltype in ['asdm', 'interpolation', 'model-fit'] and not infile is None:
            f = JyPerKReader4File(infile)
            caltable = f.get()

        if caltype in ['asdm', 'interpolation', 'model-fit'] and infile is None:
            caltable = gen_factor_via_web_api(vis, endpoint=caltype, 
                                              timeout=timeout, retry=retry, 
                                              retry_wait_time=retry_wait_time)

        _cb.specifycal(caltable=caltable, time="", spw=spw, antenna=antenna, pol=pol,
                       caltype=caltype, parameter=parameter, infile=infile,
                       uniform=uniform)

        # _cb.close()

    except UserWarning as instance:
        casalog.post('*** UserWarning *** %s' % instance, 'WARN')

    finally:
        _cb.close()

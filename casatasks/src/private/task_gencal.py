from __future__ import absolute_import
import os
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
        Subparameter of caltype='gc|gceff|tecim|jyperk'
            infile {str} -- Specifies the name of the file to read.
        subparameter of caltype='jyperk:
            endpoint {str} -- The endpoint of the Jy/K DB Web API to access.
                options are 'asdm' (default), 'model-fit', 'interpolation'.
            timeout {int} -- Maximum waiting time [sec] for the Web API access,
                defaults to 180 sec.
            retry {int} -- Number of retry when the Web API access fails,
                defaults to 3 times.
            retry_wait_time {int} -- Waiting time [sec] until next query
                when the Web API access fails, defaults to 5 sec.
        spw {str} -- The spectral windows.
        antenna {str} -- Select data based on antenna/baseline.
        pol {str} -- Polarization selection for specified parameters.
        parameter {doubleVec} -- The calibration values.
        uniform {bool} -- Assume uniform calibration values across the array.
    """

    # validate arguments
    if (caltable == ''):
        raise ValueError('A caltable name must be specified')

    if caltype == 'tecim' and not (type(infile) == str and os.path.exists(infile)):
        raise ValueError('An existing tec map must be specified in infile')

    if caltype == 'jyperk' and endpoint not in ['asdm', 'interpolation', 'model-fit']:
        raise ValueError('When the caltype is jyperk, endpoint must be one of asdm, interpolation or model-fit')

    if not ((type(vis) == str) and (os.path.exists(vis))):
        raise ValueError('Visibility data set not found - please verify the name')

    if caltype not in ['antpos', 'jyperk']:
        gencal_type = 'general'
    else:
        gencal_type = caltype

    # If developer want to add new gencal task, developer should develop a new gencal class and
    # add the class to __gencal_factory dictionary.
    gencal = __gencal_factory[gencal_type]
    gencal.gencal(vis=vis, caltable=caltable, caltype=caltype, infile=infile,
                  endpoint=endpoint, timeout=timeout, retry=retry, retry_wait_time=retry_wait_time,
                  spw=spw, antenna=antenna, pol=pol, parameter=parameter, uniform=uniform)


class GeneralGencal():
    @classmethod
    def gencal(cls, vis=None, caltable=None, caltype=None, infile='None',
               endpoint='asdm', timeout=180, retry=3, retry_wait_time=5,
               spw=None, antenna=None, pol=None,
               parameter=None, uniform=None):
        try:
            # don't need scr col for this
            _cb.open(filename=vis, compress=False, addcorr=False, addmodel=False)
            _cb.specifycal(caltable=caltable, time='', spw=spw, antenna=antenna, pol=pol,
                           caltype=caltype, parameter=parameter, infile=infile,
                           uniform=uniform)

        except UserWarning as instance:
            casalog.post('*** UserWarning *** %s' % instance, 'WARN')

        finally:
            _cb.close()


class AntposGencal():
    @classmethod
    def gencal(cls, vis=None, caltable=None, caltype=None, infile='None',
               endpoint='asdm', timeout=180, retry=3, retry_wait_time=5,
               spw=None, antenna=None, pol=None,
               parameter=None, uniform=None):
        try:
            # don't need scr col for this
            _cb.open(filename=vis, compress=False, addcorr=False, addmodel=False)

            # call a Python function to retreive ant position offsets automatically (currently EVLA only)
            if antenna == '':
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

            _cb.specifycal(caltable=caltable, time='', spw=spw, antenna=antenna, pol=pol,
                           caltype=caltype, parameter=parameter, infile=infile,
                           uniform=uniform)

        except UserWarning as instance:
            casalog.post('*** UserWarning *** %s' % instance, 'WARN')

        finally:
            _cb.close()


class JyperkGencal():
    """A class to generate caltable using the Jy/K DB or the factor CSV file.

    This class will be called if the caltype is 'jyperk'.
    """

    @classmethod
    def gencal(cls, vis=None, caltable=None, caltype=None, infile='None',
               endpoint='asdm', timeout=180, retry=3, retry_wait_time=5,
               spw=None, antenna=None, pol=None,
               parameter=None, uniform=None):
        """Generate calibration table."""
        try:
            # don't need scr col for this
            _cb.open(filename=vis, compress=False, addcorr=False, addmodel=False)

            for selection, param in \
                JyperkGencal.__gen_specifycal_input(vis=vis, spw=spw,
                                                    endpoint=endpoint, infile=infile,
                                                    timeout=timeout, retry=retry,
                                                    retry_wait_time=retry_wait_time):

                _cb.specifycal(caltable=caltable, time='', spw=selection['spw'],
                               caltype='amp', antenna=selection['antenna'],  # pol=selection['pol'],
                               parameter=param, infile='', uniform=uniform)

        except UserWarning as instance:
            casalog.post('*** UserWarning *** %s' % instance, 'WARN')

        finally:
            _cb.close()

    @classmethod
    def __gen_specifycal_input(cls, vis=None, spw='*',
                               endpoint='asdm', infile=None,
                               timeout=180, retry=3, retry_wait_time=5):
        # The default infile is defined 'string' in gencal.xml.
        if infile == '':
            infile = None

        if type(infile) is str:
            f = JyPerKReader4File(infile)
            factors = f.get()

        elif infile is None:
            factors = gen_factor_via_web_api(vis, spw=spw,
                                             endpoint=endpoint,
                                             timeout=timeout, retry=retry,
                                             retry_wait_time=retry_wait_time)
        else:
            raise Exception('The infile argument should be str or None.')

        factors = JyperkGencal.__extract_valid_factor(factors, os.path.basename(vis))
        if len(factors) == 0:
            raise Exception('There is no factor.')

        for factor in factors:
            selection = {}
            selection['vis'] = factor[0]
            selection['antenna'] = factor[1]
            selection['spw'] = factor[2]
            selection['pol'] = factor[3]

            yield selection, 1/np.sqrt(float(factor[4]))

    @classmethod
    def __extract_valid_factor(cls, factors, vis_key):
        valid_factors = []
        for factor in factors:
            if factor[0] == vis_key:
                valid_factors.append(factor)
        return valid_factors


__gencal_factory = {
    'general': GeneralGencal,
    'antpos': AntposGencal,
    'jyperk': JyperkGencal,
}

import collections
import contextlib
import itertools
import os
import shutil

import numpy as np
from casatasks.private.casa_transition import is_CASA6

if is_CASA6:
    from casatasks import casalog
    from casatasks.private import sdutil, simutil
    from casatools import ms as mstool
    from casatools import msmetadata, quanta, singledishms

    ut = simutil.simutil()
    qa = quanta()
    sdms = singledishms()

    class ATMParameterConfigurator(collections.abc.Iterator):
        def __init__(self, key, value, do_config=True):
            data = [(key, value)] if do_config else []
            self._iter = iter(data)

        def __next__(self):
            return next(self._iter)

else:
    import sdutil
    from simutil import simutil
    from taskinit import casalog, gentools
    from taskinit import msmdtool as msmetadata
    from taskinit import mstool, qa
    from taskinit import tbtool as table

    ut = simutil()
    (sdms,) = gentools(['sdms'])

    class ATMParameterConfigurator(collections.Iterator):
        def __init__(self, key, value, do_config=True):
            data = [(key, value)] if do_config else []
            self._iter = iter(data)

        def next(self):
            return next(self._iter)


@contextlib.contextmanager
def open_msmd(path):
    msmd = msmetadata()
    msmd.open(path)
    try:
        yield msmd
    finally:
        msmd.close()


def _ms_remove(path):
    if (os.path.exists(path)):
        if (os.path.isdir(path)):
            shutil.rmtree(path)
        else:
            os.remove(path)


def get_default_params():
    # Default constant: taken from atmcor_20200807.py (CSV-3320)
    atmtype = 2         ### atmType parameter for at (1: tropical, 2: mid lat summer, 3: mid lat winter, etc)
    maxalt = 120        ### maxAltitude parameter for at (km)
    lapserate = -5.6    ### dTem_dh parameter for at (lapse rate; K/km)
    scaleht = 2.0       ### h0 parameter for at (water scale height; km)

    dosmooth = False    ### convolve dTa* spectra with [0.25, 0.5, 0.25] to mimic Hanning spectral response;
    # set to True if spectral averaging was not employed for the spw
    dp = 10.0   ### initATMProfile DEFAULT ###
    dpm = 1.2   ### initATMProfile DEFAULT ###
    return locals()


def parse_gainfactor(gainfactor):
    """Parse gainfactor parameter

    Parse gainfactor parameter.

    Args:
        gainfactor (float, dict, str): gain factor.
            if float value is given, it applies to all spws.
            if dict is given, spw id and corresponding factor
            should be provided as key-value pair.
            if str is given, it should be the name of caltable.
            factors are derived as inverse-square of values
            stored in the caltable.

    Raises:
        FileNotFoundError: specified caltable does not exist.

    Returns:
        dictionary whose keys are spw id in string while values
        are the factors to be applied to each spw.
        dictionary is defined as collections.defaultdict that
        returns 1.0 as a default value.
    """
    gaindict = collections.defaultdict(lambda: 1.0)
    if isinstance(gainfactor, dict):
        # make sure keys are str
        d = dict((str(k), v) for k, v in gainfactor.items())
        gaindict.update(d)
    elif isinstance(gainfactor, str):
        # should be the name of caltable
        if not os.path.exists(gainfactor):
            raise FileNotFoundError('"{}" should exist.'.format(gainfactor))
        with sdutil.table_manager(gainfactor) as tb:
            if 'FPARAM' in tb.colnames():
                col = 'FPARAM'
            elif 'CPARAM' in tb.colnames():
                col = 'CPARAM'
            else:
                raise RuntimeError('{} is not a caltable'.format(gainfactor))
            spw_list = set(tb.getcol('SPECTRAL_WINDOW_ID'))
            for spw in spw_list:
                tsel = tb.query('SPECTRAL_WINDOW_ID=={}'.format(spw))
                try:
                    v = tsel.getcol(col).real
                finally:
                    tsel.close()
                factor = np.mean(1 / np.square(v))
                gaindict[str(spw)] = factor
    else:
        # should be float
        v = float(gainfactor)
        gaindict = collections.defaultdict(lambda: v)
    return gaindict


def gaindict2list(msname, gaindict):
    with sdutil.table_manager(os.path.join(msname, 'SPECTRAL_WINDOW')) as tb:
        nspw = tb.nrows()

    gainlist = np.ones(nspw, dtype=float)
    if isinstance(gaindict, collections.defaultdict) and len(gaindict.keys()) == 0:
        gainlist[:] = gaindict[0]
    else:
        for k, v in gaindict.items():
            spw = int(k)
            if 0 <= spw and spw < nspw:
                gainlist[spw] = v

    return gainlist


def get_all_spws_from_main(msname):
    """
    Extract all spectral window ids that have any
    associated data in MS MAIN table.

    Args:
        msname (str): name of MS

    Returns:
        list: list of available spectral window ids
    """
    with sdutil.table_manager(msname) as tb:
        ddids = np.unique(tb.getcol('DATA_DESC_ID'))
    with open_msmd(msname) as msmd:
        spws_all = [msmd.spwfordatadesc(ddid) for ddid in ddids]
    return spws_all


def get_selected_spws(msname, spw):
    """
    Get selected spectral window ids.

    Args:
        msname (str): name of MS
        spw (str): spectral window selection

    Raises:
        TypeError: spw is not string

    Returns:
        list: list of selected spectral window ids
    """
    if not isinstance(spw, str):
        raise TypeError('spw selection must be string')
    elif len(spw) == 0:
        # '' indicates all spws, which is equivalent to '*'
        spwsel = '*'
    else:
        spwsel = spw
    ms = mstool()
    sel = ms.msseltoindex(msname, spw=spwsel)
    return sel['spw']


def parse_spw(msname, spw=''):
    """Parse spw selection into list of spw ids

    Parse spw selection into list of spw ids that have
    associated data in the MAIN table of given MS.

    Args:
        msname (str): name of MS
        spw (str): spw selection

    Raises:
        TypeError: spw selection is not str
        RuntimeError: spw selection cause empty result

    Returns:
        list: list of selected spw ids
    """
    spws_all = get_all_spws_from_main(msname)
    spws_sel = get_selected_spws(msname, spw)
    spws = set(spws_all).intersection(set(spws_sel))
    return list(spws)


def get_mount_off_source_commands(msname):
    """Return list of flag commands whose reason is "Mount_is_off_source".

    Args:
        msname (str): name of MS

    Returns:
        np.ndarray: list of flag commands
    """
    with sdutil.table_manager(os.path.join(msname, 'FLAG_CMD')) as tb:
        if tb.nrows() > 0:
            tsel = tb.query('REASON=="Mount_is_off_source"')
            try:
                commands = tsel.getcol('COMMAND')
            finally:
                tsel.close()
        else:
            commands = []
    return commands


def get_antenna_name(antenna_selection):
    """Extract antenna name from the antenna selection string.

    Here, antenna_selection is assumed to be a string
    in the form '<ANTENNA_NAME>&&*'.

    Args:
        antenna_selection (str): antenna selection string

    Returns:
        str: antenna name
    """
    return antenna_selection.split('=')[1].strip("'&*")


def get_time_delta(time_range):
    """Convert time range string into time duration in sec.

    Here, time_range is assumed to be a string in the form
    'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'

    Args:
        time_range (str): time range string

    Returns:
        float: time duration in sec
    """
    timestrs = time_range.split('=')[1].strip("'").split('~')
    timequanta = [qa.quantity(t) for t in timestrs]
    timedelta = qa.convert(qa.sub(timequanta[1], timequanta[0]), 's')['value']
    return abs(timedelta)


def cmd_to_ant_and_time(cmd):
    """Extract antenna name and time duration from the flag command.

    Args:
        cmd (str): flag command

    Returns:
        tuple: antenna name and time duration
    """
    sels = cmd.split()
    asel = list(filter(lambda x: x.startswith('antenna'), sels))[0]
    tsel = list(filter(lambda x: x.startswith('time'), sels))[0]

    antenna_name = get_antenna_name(asel)
    time_delta = get_time_delta(tsel)

    return antenna_name, time_delta


def inspect_flag_cmd(msname):
    """Inspect FLAG_CMD table.

    Search flag commands whose reason is Mount_is_off_source and
    extract antenna name and time duration from the commands.

    Args:
        msname (str): name of MS

    Returns:
        tuple: two dictionaries containing the inspection result.
               The first one is number of command counts per antenna
               while the second one is total duration flagged by the
               commands per antenna.
    """
    commands = get_mount_off_source_commands(msname)

    cmd_counts = collections.defaultdict(lambda: 0)
    time_counts = collections.defaultdict(lambda: 0)

    for cmd in commands:
        ant, dt = cmd_to_ant_and_time(cmd)
        cmd_counts[ant] += 1
        time_counts[ant] += dt

    return cmd_counts, time_counts


# Argument parameter handling
def parse_atm_params(user_param, user_default, task_default, default_unit=''):
    """Parse ATM parameters.

    Args:
        user_param (str,int,float): User input.
        user_default (str,int,float): User default.
        task_default (str,int,float): Task default.
        default_unit (str): Default unit.

    Raises:
        ValueError: user_param is invalid.

    Returns:
        Tuple: Two-tuple, resulting value as quantity and boolean
               value indicating if the value is equal to user_default.
    """
    is_customized = user_param != user_default and user_param is not None

    try:
        if qa.isquantity(task_default):
            task_default_quanta = qa.quantity(task_default)
        else:
            task_default_quanta = qa.quantity(task_default, default_unit)
    except Exception as e:
        casalog.post('INTERNAL ERROR: {}'.format(e), priority='SEVERE')
        raise

    if not is_customized:
        param = task_default_quanta['value']
    else:
        user_param_quanta = qa.quantity(user_param)
        if user_param_quanta['unit'] == '':
            user_param_quanta = qa.quantity(
                user_param_quanta['value'],
                default_unit
            )
        else:
            user_param_quanta = qa.convert(
                user_param_quanta,
                default_unit
            )
        is_compatible = qa.compare(user_param_quanta, task_default_quanta)
        if is_compatible:
            param = user_param_quanta['value']
        else:
            raise ValueError('User input "{}" should have the unit compatible with "{}"'.format(
                user_param,
                default_unit
            ))

    return param, is_customized


def parse_atm_list_params(user_param, user_default='', task_default=[], default_unit=''):
    """Parse ATM parameters.

    Args:
        user_param (str,list): User input.
        user_default (str): User default.
        task_default (list): Task default.
        default_unit (str): Unit for output values.

    Raises:
        ValueError: user_param is invalid.

    Returns:
        Tuple: Two-tuple, resulting value as quantity and boolean
               value indicating if the value is equal to user_default.
    """
    is_customized = user_param != user_default and user_param is not None

    if not is_customized:
        return task_default, is_customized

    if isinstance(user_param, (list, np.ndarray)):
        try:
            param = [parse_atm_params(p, user_default, 0, default_unit=default_unit)[0] for p in user_param]
            param = [qa.convert(p, default_unit)['value'] for p in param]
        except Exception as e:
            casalog.post('ERROR during handling list input: {}'.format(e))
            raise ValueError('list input "{}" is invalid.'.format(user_param))
        return param, is_customized
    elif isinstance(user_param, str):
        try:
            split_param = user_param.split(',')
            param, _ = parse_atm_list_params(split_param, user_default, task_default, default_unit)
        except Exception as e:
            casalog.post('ERROR during handling comma-separated str input: {}'.format(e))
            raise ValueError('str input "{}" is invalid.'.format(user_param))
        return param, is_customized
    else:
        raise ValueError('user_param for parse_atm_list_params should be either list or str.')


def get_default_antenna(msname):
    """Determine default antenna id based on the FLAG_CMD table.

    Procedure is as follows.

      (1) extract flag commands whose reason is "Mount_is_off_source".
      (2) compile the commands into a number of commands and flagged
          time durations for each antenna.
      (3) select antenna with the shortest flagged duration.
      (4) if multiple antennas match in (3), select antenna with
          the least number of commands among them.
      (5) if multiple antennas match in (4), select the first
          antenna among them.

    Args:
        msname (str): name of MS

    Raises:
        Exception: no antenna was found in the MAIN table

    Returns:
        int: default antenna id
    """
    # get list of antenna Ids from MAIN table
    with sdutil.table_manager(msname) as tb:
        ant_list = np.unique(tb.getcol('ANTENNA1'))

    # No Available antenna
    if len(ant_list) == 0:
        raise Exception("No Antenna was found.")

    # get antenna names list by antenna Id
    with open_msmd(msname) as msmd:
        ant_name = [msmd.antennanames(i)[0] for i in ant_list]

    # dictionary to map antenna name to antenna Id
    ant_dict = dict((k, v) for k, v in zip(ant_name, ant_list))

    # determine default antenna id
    cmd_counts, flagged_durations = inspect_flag_cmd(msname)

    if len(cmd_counts) == 0:
        # No flag command exists. All the antennas should be healthy
        # so just pick up the first antenna.
        default_id = ant_list[0]
        default_name = ant_name[0]
    else:
        flagged_durations_filtered = dict((k, flagged_durations[k]) for k in ant_dict.keys())
        min_duration = min(flagged_durations_filtered.values())
        candidate_antennas = [k for k, v in flagged_durations_filtered.items() if v == min_duration]

        if len(candidate_antennas) == 1:
            default_name = candidate_antennas[0]
            default_id = ant_dict[default_name]
        else:
            _counts = [cmd_counts[a] for a in candidate_antennas]
            min_count = min(_counts)
            candidate_antennas2 = [a for i, a in enumerate(candidate_antennas) if _counts[i] == min_count]
            default_name = candidate_antennas2[0]
            default_id = ant_dict[default_name]
    casalog.post('Select {} (ID {}) as a default antenna'.format(default_name, default_id))
    return default_id


def get_default_altitude(msname, antid):
    """ Get default altitude of the antenna
    decide default value of 'Altitude' for Atm Correction.
    This requires to calculate Elevation from Antenna Position Information.
    """
    with sdutil.table_manager(os.path.join(msname, 'ANTENNA')) as tb:
        # obtain the antenna Position (Earth Center) specified by antid
        X, Y, Z = (float(i) for i in tb.getcell('POSITION', antid))

        #  xyz2long()   -- https://casa.nrao.edu/casadocs/casa-5.6.0/simulation/simutil
        #
        #  When given ITRF Earth-centered (X, Y, Z, using the parameters x, y, and z) coordinates [m] for a point,
        #  this method returns geodetic latitude and longitude [radians] and elevation [m].
        #  Elevation is measured relative to the closest point to the (latitude, longitude)
        #  on the WGS84 (World Geodetic System 1984) reference ellipsoid.

        P = ut.xyz2long(X, Y, Z, 'WGS84')   # [0]:longitude, [1]:latitude, [2]:elevation (geodetic elevation)
        geodetic_elevation = P[2]

        ref = tb.getcolkeyword('POSITION', 'MEASINFO')['Ref']

    casalog.post("Default Altitude")
    casalog.post(" - Antenna ID: %d. " % antid)
    casalog.post(" - Ref = %s. " % ref)
    casalog.post(" - Position: (%s, %s, %s)." % (X, Y, Z))
    casalog.post("   Altitude (geodetic elevation):  %f" % geodetic_elevation)

    return geodetic_elevation


class ATMScalarParameterConfigurator(ATMParameterConfigurator):
    def __init__(self, key, user_input, impl_default, default_unit, api_default='', is_mandatory=True, is_effective=True):
        value, is_customized = parse_atm_params(user_param=user_input, user_default=api_default, task_default=impl_default, default_unit=default_unit)
        do_config = is_mandatory or (is_effective and is_customized)
        # TODO: remove arguments for super, i.e. just super().__init__(...) once we completely get rid of CASA5
        super(self.__class__, self).__init__(key=key, value=value, do_config=do_config)


class ATMListParameterConfigurator(ATMParameterConfigurator):
    def __init__(self, key, user_input, impl_default, default_unit, api_default='', is_mandatory=True, is_effective=True):
        value, is_customized = parse_atm_list_params(user_param=user_input, user_default=api_default, task_default=impl_default, default_unit=default_unit)
        do_config = is_mandatory or (is_effective and is_customized)
        # TODO: remove arguments for super, i.e. just super().__init__(...) once we completely get rid of CASA5
        super(self.__class__, self).__init__(key=key, value=value, do_config=do_config)


def get_configuration_for_atmcor(infile, spw, outputspw, gainfactor, user_inputs):
    # requested list of output spws and processing spws
    # processing spws are the intersection of these
    outputspws_param = parse_spw(infile, outputspw)
    spws_param = parse_spw(infile, spw)
    all_processing_spws = np.asarray(list(set(spws_param).intersection(set(outputspws_param))))

    # generate gain factor dictionary
    gaindict = parse_gainfactor(gainfactor)
    gainlist = gaindict2list(infile, gaindict)

    # default parameter values (from Tsuyoshi's original script)
    default_params = get_default_params()

    # reference antenna_id to calculate Azimuth/Elevation
    reference_antenna = int(get_default_antenna(infile))

    # altitude of reference antenna
    default_altitude = get_default_altitude(infile, reference_antenna)
    user_altitude = user_inputs['altitude'] if user_inputs['atmdetail'] else ''

    parameters = [
        ATMParameterConfigurator(key='processspw', value=all_processing_spws),
        ATMParameterConfigurator(key='gainfactor', value=gainlist),
        ATMParameterConfigurator(key='refant', value=reference_antenna),
        ATMParameterConfigurator(key='atmType', value=user_inputs['atmtype']),
        ATMParameterConfigurator(key='maxAltitude', value=float(default_params['maxalt'])),
        ATMScalarParameterConfigurator(
            key='lapseRate', user_input=user_inputs['dtem_dh'],
            impl_default=default_params['lapserate'], default_unit='K/km',
        ),
        ATMScalarParameterConfigurator(
            key='scaleHeight', user_input=user_inputs['h0'],
            impl_default=default_params['scaleht'], default_unit='km',
        ),
        ATMScalarParameterConfigurator(
            key='pressureStep', user_input=user_inputs['dp'],
            impl_default=default_params['dp'], default_unit='mbar'
        ),
        ATMScalarParameterConfigurator(
            key='pressureStepFactor', user_input=user_inputs['dpm'],
            impl_default=default_params['dpm'], default_unit='', api_default=-1
        ),
        ATMScalarParameterConfigurator(
            key='siteAltitude', user_input=user_altitude,
            impl_default=default_altitude, default_unit='m',
            api_default='',
            is_mandatory=True, is_effective=user_inputs['atmdetail']
        ),
        ATMScalarParameterConfigurator(
            key='pressure', user_input=user_inputs['pressure'],
            impl_default=0, default_unit='mbar',
            is_mandatory=False, is_effective=user_inputs['atmdetail']
        ),
        ATMScalarParameterConfigurator(
            key='temperature', user_input=user_inputs['temperature'],
            impl_default=0, default_unit='K',
            is_mandatory=False, is_effective=user_inputs['atmdetail']
        ),
        ATMScalarParameterConfigurator(
            key='humidity', user_input=user_inputs['humidity'],
            impl_default=0, default_unit='%',
            api_default=-1,
            is_mandatory=False, is_effective=user_inputs['atmdetail']
        ),
        ATMScalarParameterConfigurator(
            key='pwv', user_input=user_inputs['pwv'],
            impl_default=0, default_unit='mm',
            is_mandatory=False, is_effective=user_inputs['atmdetail']
        ),
        ATMListParameterConfigurator(
            key='layerBoundaries', user_input=user_inputs['layerboundaries'],
            impl_default=[], default_unit='m',
            is_mandatory=False, is_effective=user_inputs['atmdetail']
        ),
        ATMListParameterConfigurator(
            key='layerTemperatures', user_input=user_inputs['layertemperature'],
            impl_default=[], default_unit='K',
            is_mandatory=False, is_effective=user_inputs['atmdetail']
        )
    ]

    config = dict(itertools.chain(*parameters))

    return config


@sdutil.sdtask_decorator
def sdatmcor(
        infile=None, datacolumn=None, outfile=None, overwrite=None,
        field=None, spw=None, scan=None, antenna=None,
        correlation=None, timerange=None, intent=None,
        observation=None, feed=None, msselect=None,
        outputspw=None,
        gainfactor=None,
        dtem_dh=None, h0=None, atmtype=None,
        atmdetail=None,
        altitude=None, temperature=None, pressure=None, humidity=None, pwv=None,
        dp=None, dpm=None,
        layerboundaries=None, layertemperature=None):

    try:
        # Input/Output error check and internal set up.
        if infile == '':
            errmsg = "infile MUST BE specified."
            raise Exception(errmsg)

        if outfile == '':
            errmsg = "outfile MUST BE specified."
            raise Exception(errmsg)

        # Protection, in case infile == outfile
        if infile == outfile:
            errmsg = "You are attempting to write the output on your input file."
            raise Exception(errmsg)

        # File Info
        casalog.post("INPUT/OUTPUT")
        casalog.post("  Input MS file   = %s " % infile)
        casalog.post("  Output MS file  = %s " % outfile)

        # infile Inaccessible
        if not os.path.exists(infile):
            errmsg = "Specified infile does not exist."
            raise Exception(errmsg)

        # outfile Protected
        if os.path.exists(outfile):
            if overwrite:
                casalog.post("Overwrite:: Overwrite specified. Once delete the existing output file. ")
                _ms_remove(outfile)
            else:
                errmsg = "Specified outfile already exist."
                raise Exception(errmsg)

        # Inspect atmtype
        atmtype_int = int(atmtype)
        if atmtype_int not in (1, 2, 3, 4, 5):
            errmsg = "atmtype (=%s) should be any one of (1, 2, 3, 4, 5)." % atmtype
            raise Exception(errmsg)

        # Inspect humidity (float). The range must be 0.0 gt. Humidity  gt. 100.0 [%]
        humidity_float = float(humidity)
        if humidity_float != -1.0 and not (0.0 <= humidity_float and humidity_float <= 100.0):
            errmsg = "humidity (=%s) should be in range 0~100" % humidity
            raise Exception(errmsg)

        # datacolumn check (by XML definition)
        datacolumn_upper = datacolumn.upper()
        if datacolumn_upper not in ['DATA', 'CORRECTED', 'FLOAT_DATA']:
            errmsg = "Specified column name (%s) Unacceptable." % datacolumn
            raise Exception(errmsg)

        # tweak antenna selection string to include autocorr data
        antenna_autocorr = sdutil.get_antenna_selection_include_autocorr(infile, antenna)

        # C++ re-implementation
        sdms.open(infile)
        sdms.set_selection(spw=outputspw, field=field,
                           antenna=antenna_autocorr,
                           timerange=timerange, scan=scan,
                           polarization=correlation, intent=intent,
                           observation=observation, feed=feed,
                           taql=msselect,
                           reindex=False)

        config = get_configuration_for_atmcor(
            infile=infile,
            spw=spw,
            outputspw=outputspw,
            gainfactor=gainfactor,
            user_inputs=locals()
        )

        sdms.atmcor(config=config, datacolumn=datacolumn, outfile=outfile)

    except Exception as err:
        casalog.post('%s' % err, priority='SEVERE')
        raise

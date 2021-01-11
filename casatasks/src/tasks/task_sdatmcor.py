import os
import pylab as pl
import shutil
import contextlib
import collections
import numpy as np

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
    from casatasks.private import sdutil
    from casatools import quanta, table, msmetadata
    from casatools import atmosphere
    from casatools import ms as mstool
    from casatasks import mstransform
    import casatasks.private.simutil as simutil

    ut = simutil.simutil()
    qa = quanta()
    at = atmosphere()

else:
    from taskinit import tbtool as table
    from taskinit import mstool, casalog, qa
    from taskinit import msmdtool as msmetadata
    from casac import casac
    from tasks import mstransform
    from simutil import simutil
    import sdutil

    ut = simutil()
    at = casac.atmosphere()

# Task name
origin = 'sdatmcor'


@contextlib.contextmanager
def open_table(path, nomodify=True):
    tb = table()
    tb.open(path, nomodify=nomodify)
    try:
        yield tb
    finally:
        tb.close()


@contextlib.contextmanager
def open_ms(path):
    ms = mstool()
    ms.open(path)
    try:
        yield ms
    finally:
        ms.close()


@contextlib.contextmanager
def open_msmd(path):
    msmd = msmetadata()
    msmd.open(path)
    try:
        yield msmd
    finally:
        msmd.close()


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
        with open_table(gainfactor) as tb:
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
    if not isinstance(spw, str):
        raise TypeError('spw selection must be string')

    with open_table(msname) as tb:
        ddids = np.unique(tb.getcol('DATA_DESC_ID'))
    with open_msmd(msname) as msmd:
        spws_all = [msmd.spwfordatadesc(ddid) for ddid in ddids]

    if len(spw) == 0:
        # '' indicates all spws, which is equivalent to '*'
        spwsel = '*'
    else:
        spwsel = spw
    ms = mstool()
    sel = ms.msseltoindex(msname, spw=spwsel)
    spws = set(spws_all).intersection(set(sel['spw']))
    return list(spws)


def get_mount_off_source_commands(msname):
    """Return list of flag commands whose reason is "Mount_is_off_source".

    Args:
        msname (str): name of MS

    Returns:
        np.ndarray: list of flag commands
    """
    with open_table(os.path.join(msname, 'FLAG_CMD')) as tb:
        tsel = tb.query('REASON=="Mount_is_off_source"')
        try:
            commands = tsel.getcol('COMMAND')
        finally:
            tsel.close()
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

    # task name
    casalog.origin(origin)

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
        _msg("INPUT/OUTPUT")
        _msg("  Input MS file   = %s " % infile)
        _msg("  Output MS file  = %s " % outfile)

        # infile Inaccessible
        if not _file_exist(infile):
            errmsg = "Specified infile does not exist."
            raise Exception(errmsg)

        # outfile Protected
        if _file_exist(outfile):
            if overwrite:
                _msg("Overwrite:: Overwrite specified. Once delete the existing output file. ")
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

        # generate gain factor dictionary
        gaindict = parse_gainfactor(gainfactor)

        # Data Selection Section
        _msg("Calling mstransform for Data Selection. Output file = %s " % outfile)

        # datacolumn check (by XML definition)
        datacolumn = datacolumn.upper()
        if datacolumn not in ['DATA', 'CORRECTED', 'FLOAT_DATA']:
            errmsg = "Specified column name (%s) Unacceptable." % datacolumn
            raise Exception(errmsg)

        # tweak antenna selection string to include autocorr data
        antenna_autocorr = sdutil.get_antenna_selection_include_autocorr(infile, antenna)

        mstransform(
            vis=infile,
            outputvis=outfile,
            datacolumn=datacolumn,
            field=field,
            spw=outputspw,   # Use 'outputspw' for Data Selection.
            scan=scan,
            antenna=antenna_autocorr,
            correlation=correlation,
            timerange=timerange,
            intent=intent,
            observation=observation,
            feed=feed,
            taql=msselect,
            reindex=False)   # Must be False

        # resume 'origin'. A strange behavior in casalog/CASA6
        casalog.origin(origin)

        # result check if output was generated
        if not _file_exist(outfile):
            errmsg = "No outfile was generated by mstransform."
            raise Exception(errmsg)

        # data column name is always 'DATA' after mstransform whatever datacolumn input MS has
        datacolumn_name = 'DATA'

        # Call main body Function
        calc_sdatmcor(
            infile, datacolumn_name, outfile,
            spw,
            gaindict,
            dtem_dh, h0, atmtype_int,
            atmdetail,
            altitude, temperature, pressure, humidity_float, pwv, dp, dpm,
            layerboundaries,
            layertemperature)

    except Exception as err:
        casalog.post('%s' % err, 'SEVERE')
        raise


def _ms_remove(path):
    if (os.path.exists(path)):
        if (os.path.isdir(path)):
            shutil.rmtree(path)
        else:
            os.remove(path)


def _file_exist(path):
    return os.path.exists(path)


# Argument parameter handling
def parse_atm_params(user_param, user_default, task_default):
    """Parse ATM parameters.

    Args:
        user_param (str,int,float): User input.
        user_default (str,int,float): User default.
        task_default (str,int,float): Task default.

    Raises:
        ValueError: user_param is invalid.

    Returns:
        Tuple: Two-tuple, resulting value as quantity and boolean
               value indicating if the value is equal to user_default.
    """
    is_customized = user_param != user_default and user_param is not None

    try:
        task_default_quanta = qa.quantity(task_default)
        default_unit = task_default_quanta['unit']
    except Exception as e:
        _msg('INTERNAL ERROR: {}'.format(e), 'SEVERE')
        raise

    if not is_customized:
        param = task_default_quanta
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
            param = user_param_quanta
        else:
            raise ValueError('User input "{}" should have the unit compatible with "{}"'.format(
                user_param,
                default_unit
            ))

    return param, is_customized


def parse_atm_list_params(user_param, user_default='', task_default=[], element_unit=''):
    """Parse ATM parameters.

    Args:
        user_param (str,list): User input.
        user_default (str): User default.
        task_default (list): Task default.
        element_unit (str): Unit for output values.

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
            param = [parse_atm_params(p, user_default, qa.quantity(0, element_unit))[0] for p in user_param]
            param = [qa.convert(p, element_unit)['value'] for p in param]
        except Exception as e:
            _msg('ERROR during handling list input: {}'.format(e))
            raise ValueError('list input "{}" is invalid.'.format(user_param))
        return param, is_customized
    elif isinstance(user_param, str):
        try:
            split_param = user_param.split(',')
            param, _ = parse_atm_list_params(split_param, user_default, task_default, element_unit)
        except Exception as e:
            _msg('ERROR during handling comma-separated str input: {}'.format(e))
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
    with open_table(msname) as tb:
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
    _msg('Select {} (ID {}) as a default antenna'.format(default_name, default_id))
    return default_id


def get_default_altitude(msname, antid):
    """ Get default altitude of the antenna
    decide default value of 'Altitude' for Atm Correction.
    This requires to calculate Elevation from Antenna Position Information.
    """
    with open_table(os.path.join(msname, 'ANTENNA')) as tb:
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

    _msg("Default Altitude")
    _msg(" - Antenna ID: %d. " % antid)
    _msg(" - Ref = %s. " % ref)
    _msg(" - Position: (%s, %s, %s)." % (X, Y, Z))
    _msg("   Altitude (geodetic elevation):  %f" % geodetic_elevation)

    return geodetic_elevation


#
# show ATM Profile
# https://casa.nrao.edu/docs/CasaRef/atmosphere-Tool.html#x995-10120004.1.1
#

def show_atm_info(atm):
    """
     Returned atm (from initAtmProfile) may have a different structure.
     In casa6, additional variable information is added in Dict type.
    """
    version = at.getAtmVersion()
    _msg("\nAtomosphere Tool:: version = %s\n" % version)

    # ATM Profile
    if is_CASA6:
        _msg(atm[0])
    else:
        _msg(atm)


def show_layer_info(at):
    p = at.getProfile()

    # From Casa Tool Documentation:
    #
    # returns a tuple of
    # 0 - string listing of layer values, and arrays of layer, 1 - thickness,
    # 2 - temperature, 3 - watermassdensity, 4 - water (number density),
    # 5 - pressure, 6 - O3 (number density), 7 - CO, 8 - N2O

    # for i in range(at.getNumLayers()):
    #     # Print atmospheric profile returned by at.getProfile():
    #     # Layer thickness (idx=1), Temperature (idx=2),
    #     # Number density of water vapor(idx=4), and Pressure (idx=5)
    #     print(p[1]['value'][i],p[2]['value'][i],p[4]['value'][i],p[5]['value'][i])

    _msg(p[0])
    return


#
# Logging CASA LOG (INFO/WARN/SEVERE)
#
def _msg(msg, priority='INFO'):
    # Information message
    # if priority == 'INFO':
    #     print(msg)
    # other Warning/Error message
    casalog.post(msg, priority=priority, origin=origin)


############################################################
# Calculation Method (Replaced to C++ in next development.)
#    originated by atmcorr_20200807.py by Sawada san.
#    formed as a casa task by CAS-13160.
#    - p_xxxx arguments are for task.
#    - param_xxxx arguments are for Atm-correction.
############################################################


def calc_sdatmcor(
        p_infile,
        p_datacolumn,
        p_outfile,
        p_spw,
        gaindict,
        param_dtem_dh,
        param_h0,
        param_atmtype,
        atmdetail,
        param_altitude,
        param_temperature,
        param_pressure,
        param_humidity,
        param_PWV,
        param_dp,
        param_dpm,
        param_layerboundaries,
        param_layertemperature):

    # Argument dump.
    #   if needed, set True.
    if True:
        _msg("***********************************")
        _msg("**   calc_sdatmcor::             **")
        _msg("***********************************")
        _msg('infile      = %s' % p_infile)
        _msg('datacolumn  = %s' % p_datacolumn)
        _msg('outfile     = %s' % p_outfile)
        _msg('spw         = %s' % p_spw)
        _msg('dtem_dh     = %s' % param_dtem_dh)
        _msg('h0          = %s' % param_h0)
        _msg('atmtype     = %s' % param_atmtype)
        _msg('atmdetail   = %s' % atmdetail)
        _msg('altitude    = %s' % param_altitude)
        _msg('temperature = %s' % param_temperature)
        _msg('pressure    = %s' % param_pressure)
        _msg('humidity    = %s' % param_humidity)
        _msg('PWV         = %s' % param_PWV)
        _msg('dp          = %s' % param_dp)
        _msg('dpm         = %s' % param_dpm)
        _msg('layerboundaries   = %s' % param_layerboundaries)
        _msg('layertemperature  = %s' % param_layertemperature)
        _msg("*****************************")

    # Internal file names
    rawms = p_infile
    calms = p_outfile
    corms = p_outfile

    # CAS-13160 Original Script starts from here
    #   - in many sections, Task Code is inserted.
    #   - comment is generally added to distinguish the original and the task code.

    ##################################################
    #   Inside Constant for ATM
    ##################################################

    # Following variables have initial parameters according to the Original script.
    #  - these values will be translated with unit to var.:'atm_xxxxxxxxx' in the Task code.

    # Default constant
    atmtype = 2         ### atmType parameter for at (1: tropical, 2: mid lat summer, 3: mid lat winter, etc)
    maxalt = 120        ### maxAltitude parameter for at (km)
    lapserate = -5.6    ### dTem_dh parameter for at (lapse rate; K/km)
    scaleht = 2.0       ### h0 parameter for at (water scale height; km)

    dosmooth = False    ### convolve dTa* spectra with [0.25, 0.5, 0.25] to mimic Hanning spectral response;
    # set to True if spectral averaging was not employed for the spw
    dp = 10.0   ### initATMProfile DEFAULT ###
    dpm = 1.2   ### initATMProfile DEFAULT ###

    # other
    nband = 1         # number of band in initSpectralWindow()

    # normalized spectral response for Hanning window, FWHM=10
    hanning = [-0.00098041, -0.00202866, -0.00265951, -0.00222265,
               0.00000000, 0.00465696, 0.01217214, 0.02260546, 0.03556241,
               0.05017949, 0.06519771, 0.07911911, 0.09042184, 0.09779662,
               0.10035898, 0.09779662, 0.09042184, 0.07911911, 0.06519771,
               0.05017949, 0.03556241, 0.02260546, 0.01217214, 0.00465696,
               0.00000000, -0.00222265, -0.00265951, -0.00202866, -0.00098041]

    # Tool Default constant (CAS-13160)
    # - Antenna and its altitude are determined.

    antenna = get_default_antenna(rawms)  # default antenna_id
    altitude = get_default_altitude(rawms, antenna)  # default altitude

    # CAS-13160
    #     From here, main part of the Original script starts.
    #     In some sections, additional statements for Task are inserted.

    ################################################################
    # Get metadata
    ################################################################
    chanfreqs = {}
    _msg("\nSDATMCOR main body starts. rawms=%s\n" % rawms)

    with open_msmd(rawms) as msmd:
        # (original)
        # tmonsource = msmd.timesforintent('OBSERVE_TARGET#ON_SOURCE')
        tmoffsource = msmd.timesforintent('OBSERVE_TARGET#OFF_SOURCE')
        fdmspws = msmd.fdmspws()
        tdmspws = msmd.tdmspws()
        intentspws = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
        spws = list(set(intentspws) & (set(fdmspws) | set(tdmspws)))
        spwnames = msmd.namesforspws(spws)
        for spwid in spws:
            chanfreqs[spwid] = msmd.chanfreqs(spw=spwid)

    bnd = (pl.diff(tmoffsource) > 1)
    w1 = pl.append([True], bnd)
    w2 = pl.append(bnd, [True])
    tmoffsource = (tmoffsource[w1] + tmoffsource[w2]) / 2.  ### midpoint of OFF subscan

    ddis = {}
    with open_msmd(calms) as msmd:
        for spwid in spws:
            ddis[spwid] = msmd.datadescids(spw=spwid)[0]

    nchanperbb = [0, 0, 0, 0]
    bbprs = {}
    for i, spwid in enumerate(spws):
        bbp = int(spwnames[i].split('#')[2][3]) - 1
        bbprs[spwid] = bbp
        nchanperbb[bbp] += len(chanfreqs[spwid])

    # show initial spw info (DEBUG)
    _msg(" - fdm    spws = %s" % fdmspws, 'DEBUG1')
    _msg(" - tdm    spws = %s" % tdmspws, 'DEBUG1')
    _msg(" - intent spws = %s" % intentspws, 'DEBUG1')
    _msg(" -        spws = %s" % spws, 'DEBUG1')

    # (Task Section)
    # OFF_SOURCE Error check
    n_tmoffsource = len(tmoffsource)
    msg = "# OFF_SOURCE: count of tmoffsource = %d" % n_tmoffsource
    _msg(msg)
    if (n_tmoffsource == 0):
        msg = "Can't find the OFF_SOURCE data."
        _msg(msg, 'SEVERE')
        raise(msg)

    # requested list of output spws
    outputspws_param = parse_spw(corms, '')

    # requested list of processing spws
    spws_param = parse_spw(rawms, p_spw)

    # processing_spws must be a subset of spws as well as outputspws_param
    all_processing_spws = set(spws).intersection(set(outputspws_param))
    if not set(spws_param).issubset(all_processing_spws):
        _msg("Some of the specified spw(s) cannot be processed. Try possible one(s)", 'WARN')
    processing_spws = list(all_processing_spws.intersection(set(spws_param)))

    # Spw:: No Target check
    if len(processing_spws) == 0:
        raise Exception("No available Spw Targets. Abort.")

    _msg("Final Determined Spws Information")
    _msg('-- Target     Spws = %s' % spws)
    _msg('-- Output     Spws = %s' % outputspws_param)
    _msg('-- Requested  Spws = %s' % spws_param)
    _msg('-- Correcting Spws = %s' % processing_spws)

    # (Original)
    # Data Query on Pointing Table.
    #   output
    #   - tmpointing, elev

    with open_table(os.path.join(calms, 'POINTING')) as tb:
        # (org.) key for ANTENNA_ID select
        querytext = 'ANTENNA_ID==%s' % antenna
        subtb = tb.query(querytext)

        # (org.) Access Table
        tmpointing = subtb.getcol('TIME')
        elev = subtb.getcol('DIRECTION').squeeze()[1]
        _msg("- reading column:'TIME' and 'DIRECTION' completed.")

        subtb.close()

    ################################################################
    # Get atmospheric parameters for ATM
    ################################################################

    # ASDM_CALWVR
    with open_table(os.path.join(rawms, 'ASDM_CALWVR')) as tb:
        # confirm
        # _msg("tmonsource: %f, %f" % (tmonsource.min(), tmonsource.max()))
        tmpwv_all = tb.getcol('startValidTime')
        pwv_all = tb.getcol('water')

    # ASDM_CALATMOSPHERE
    with open_table(os.path.join(rawms, 'ASDM_CALATMOSPHERE')) as tb:
        tmatm_all = tb.getcol('startValidTime')
        tground_all = tb.getcol('groundTemperature')
        pground_all = tb.getcol('groundPressure')
        hground_all = tb.getcol('groundRelHumidity')

    # _msg("median PWV = %fm, T = %fK, P = %fPa, H = %f%%" % (pwv, tground, pground, hground))
    pwv, tground, pground, hground = [], [], [], []
    tmatm = pl.unique(tmatm_all)
    for tt in tmatm:
        deltat = abs(tmpwv_all - tt)
        pwv.append(pl.median(pwv_all[deltat == deltat.min()]))
        tground.append(pl.median(tground_all[tmatm_all == tt]))
        pground.append(pl.median(pground_all[tmatm_all == tt]))
        hground.append(pl.median(hground_all[tmatm_all == tt]))
        _msg('PWV = %fm, T = %fK, P = %fPa, H = %f%% at %s' % (pwv[-1], tground[-1], pground[-1], hground[-1], qa.time('%fs' % tt, form='fits')[0]))

    ################################################################
    # Looping over spws
    ################################################################
    with open_table(corms, nomodify=False) as tb:
        # Note CAS-13160:
        #     spw for-loop:
        #       - Inside this loop, calculation is executed on each spwid.
        #       - The corrected result is written upon OutputFile with selected Spw.
        #       - please see the inserted block.
        #     The original script does not explicitly show outputspw/spw handling logic.
        #     When intended spw is not in 'spws', no-corrected output is performed.

        # (original)
        for spwid in processing_spws:  # (original) for spwid in spws
            _msg("\nProcessing spw %d in %s. \n" % (spwid, processing_spws))

            # gain factor
            spwkey = str(spwid)
            factor = gaindict[spwkey]
            _msg('Applying gain factor {} to spw {}'.format(factor, spwkey))

            istdm = False
            nchan = len(chanfreqs[spwid])
            # fcenter = (chanfreqs[spwid][nchan/2-1]+chanfreqs[spwid][nchan/2])/2.             # PY2
            fcenter = (chanfreqs[spwid][int(nchan / 2) - 1] + chanfreqs[spwid][int(nchan / 2)]) / 2.   # PY3
            chansep = (chanfreqs[spwid][-1] - chanfreqs[spwid][0]) / (nchan - 1)

            _msg('- checking nchanperbb[bbprs[spwid]]')
            if nchanperbb[bbprs[spwid]] in [128, 256]:

                _msg('Spw %d seems to be TDM-like. More accurate Hanning smoothing is applied. istdm=True' % spwid)

                # (orginal comment)
                # ATM model is calculated for finer spectral resolution (5x original),
                # convolved with Hanning spectral response with FWHM=10,
                # then resampled to the original channel freqs

                istdm = True
                nchan *= 5
                chansep /= 5.
                chanfreqs_high = chanfreqs[spwid][0] + chansep * pl.arange(nchan)
                # fcenter = (chanfreqs_high[nchan/2-1]+chanfreqs_high[nchan/2])/2.              # PY2
                fcenter = (chanfreqs_high[int(nchan / 2) - 1] + chanfreqs_high[int(nchan / 2)]) / 2.    # PY3

            ################################################################
            # Calculate and apply correction values
            ################################################################

            _msg("Calculate and apply correction values.")

            # original: make ddis[spwid]
            _msg("- Selecting DATA_DESC_ID == %s" % ddis[spwid])

            # Query Text (DESC_ID)
            querytext = 'DATA_DESC_ID in %s' % ddis[spwid]

            # Query
            subtb = tb.query(querytext)

            # Data
            _msg("- getting tm and data. datacolumn [%s] is used." % p_datacolumn)
            tmdata = subtb.getcol('TIME')
            data = subtb.getcol(p_datacolumn)
            npol = data.shape[0]

            # Smoothing control
            if nchanperbb[bbprs[spwid]] * npol in [256, 8192]:
                _msg('Spw %d in BB_%d (total Nchan within BB is %d, sp avg likely not applied).  dosmooth=True' %
                     (spwid, bbprs[spwid] + 1, nchanperbb[bbprs[spwid]] * npol))
                dosmooth = True
            else:
                _msg('Spw %d in BB_%d (total Nchan within BB is %d, sp avg likely applied).  dosmooth=False' %
                     (spwid, bbprs[spwid] + 1, nchanperbb[bbprs[spwid]] * npol))
                dosmooth = False

            ###########################
            # Correction Main Loop
            ###########################

            _msg("\nExecuting ATM Correction(N=%d), and writing to output MS. \n" % len(tmdata))

            prevtmatm = 0.
            cdata = data.copy()
            for i, t in enumerate(tmdata):

                dt = tmoffsource - t
                if (dt < 0).sum() == 0 or (dt > 0).sum() == 0:  ### any data before first OFF or after last OFF are disregarded
                    continue

                dtoff0 = dt[dt < 0].max()
                dtoff1 = dt[dt > 0].min()
                eon, eoff0, eoff1 = pl.interp([t, t + dtoff0, t + dtoff1], tmpointing, elev)
                eoff = (dtoff1 * eoff0 - dtoff0 * eoff1) / (dtoff1 - dtoff0)

                tt = tmatm[tmatm < t]  # timestamps of CalAtmosphere results earlier than t
                if len(tt) == 0:
                    continue
                if tt[-1] != prevtmatm:
                    ################################################################
                    # Set parameters for ATM and obtain zenith opacity
                    # caveat: median values are used for calibrating the entire execution
                    ################################################################

                    _msg("- set parameters for initATM and obtain zenith opacity")

                    # Note CAS-13160: Inserted part to the original.
                    #   Setting up parameter variables to initAtmProfile()
                    #     - Following codes are extracted from the original script, which used to be
                    #       embeded inside the initAtmProfile() arguments.
                    #     - Overwrite the 'argument parameter' to 'atm_xxxx' arguments for initAtmProfile.
                    #     - The default values are up to official Task Specification.

                    _msg('Initializing ATM profile.  Data timestamp %s, weather data timestamp %s' %
                         (qa.time('%fs' % t, form='fits')[0], qa.time('%fs' % tt[-1], form='fits')[0]))

                    idx = pl.where(tmatm == tt[-1])[0][0]
                    atm_dtem_dh = qa.quantity(lapserate, 'K/km')
                    atm_h0 = qa.quantity(scaleht, 'km')
                    atm_atmtype = atmtype                                # int
                    atm_atmdetail = atmdetail                            # Bool
                    atm_altitude = qa.quantity(altitude, 'm')
                    atm_temperature = qa.quantity(tground[idx], 'K')          # tground (original)
                    atm_pressure = qa.quantity(pground[idx] / 100.0, 'mbar')  # pground (original) in  [Pa]  convert to [mbar]
                    atm_humidity = qa.quantity(hground[idx], '%')             # hground (original) in  [%]
                    atm_pwv = qa.quantity(pwv[idx] * 1000.0, 'mm')            # pwv (original) in [m] converto to [mm]
                    atm_dp = qa.quantity(dp, 'mbar')
                    atm_dpm = qa.quantity(dpm, '')                      # float
                    atm_maxAltitude = qa.quantity(maxalt, 'km')

                    # Edit Flag (True: the parameter was given and will be applied to initAtmProfile)
                    atm_altitude_set = False
                    atm_temperature_set = False
                    atm_pressure_set = False
                    atm_humidity_set = False
                    atm_pwv_set = False
                    atm_dp_set = False
                    atm_dpm_set = False

                    # (from 'help' information)
                    # User-Defined Profile (example)
                    #    myalt = [ 5071.72200397, 6792.36546384, 15727.0776121, 42464.18192672 ] #meter
                    #    mytemp = [ 270., 264., 258., 252. ] #Kelvin
                    atm_layerboundaries = []   # initially null-List.   array: layer boundary [m]
                    atm_layertemperature = []   # initially null-list.   array: layer temperature [K]

                    # ATM fundamental parameters activation.
                    atm_dtem_dh, atm_dtem_dh_set = parse_atm_params(
                        user_param=param_dtem_dh,
                        user_default='',
                        task_default=atm_dtem_dh
                    )
                    atm_h0, atm_h0_set = parse_atm_params(
                        user_param=param_h0,
                        user_default='',
                        task_default=atm_h0
                    )
                    atm_atmtype = int(param_atmtype)
                    if atm_atmtype not in (1, 2, 3, 4, 5):
                        raise ValueError('Given atmtype "{}" is invalid. Should be 1~5.'.format(param_atmtype))

                    # Sub parameter activation
                    #  when atmdetail is True
                    if atm_atmdetail:
                        _msg("\nSub Parameter from the Arguments will be set, if specified.\n")
                        atm_altitude, atm_altitude_set = parse_atm_params(
                            user_param=param_altitude,
                            user_default='',
                            task_default=atm_altitude
                        )
                        atm_temperature, atm_temperature_set = parse_atm_params(
                            user_param=param_temperature,
                            user_default='',
                            task_default=atm_temperature)
                        atm_pressure, atm_pressure_set = parse_atm_params(
                            user_param=param_pressure,
                            user_default='',
                            task_default=atm_pressure
                        )
                        atm_humidity, atm_humidity_set = parse_atm_params(
                            user_param=param_humidity,
                            user_default=-1,
                            task_default=atm_humidity
                        )
                        atm_pwv, atm_pwv_set = parse_atm_params(
                            user_param=param_PWV,
                            user_default='',
                            task_default=atm_pwv
                        )
                        atm_dp, atm_dp_set = parse_atm_params(
                            user_param=param_dp,
                            user_default='',
                            task_default=atm_dp
                        )
                        atm_dpm, atm_dpm_set = parse_atm_params(
                            user_param=param_dpm,
                            user_default=-1,
                            task_default=atm_dpm
                        )

                        # User-Defined Profile.
                        #  when the arg is active, directly pass the arg(=param_layerXXXXX) to atm_XXXXX for initAtmProfile.
                        atm_layerboundaries, _ = parse_atm_list_params(
                            user_param=param_layerboundaries,
                            user_default='',
                            task_default=[],
                            element_unit='m'
                        )
                        atm_layertemperature, _ = parse_atm_list_params(
                            user_param=param_layertemperature,
                            user_default='',
                            task_default=[],
                            element_unit='K'
                        )
                        _msg('layerboundaries {} layertemperature: {}'.format(atm_layerboundaries, atm_layertemperature))
                        if len(atm_layerboundaries) != len(atm_layertemperature):
                            raise ValueError('list length of layerboundaries and layertemperature should be the same.')

                    else:
                        _msg("\nSub Parameters were not used, due to 'atmdetail' is not True.\n")

                    # determined parameters to give to initAtmProfile
                    _msg("====================================================================")
                    _msg("  initATMProfile Parameters     [atmdetail = %s]" % atm_atmdetail)
                    _msg("------------------+-----------------------------+-------------------")
                    _msg("  parameter       | value [unit]                | specified by arg. ")
                    _msg("------------------+-----------------------------+-------------------")
                    _msg(" atmtype          |%-18s " % at.listAtmosphereTypes()[atm_atmtype - 1])    # type =1,2,3,4,5
                    _msg(" dTem_dh          |%-18s [%5s]   | %s " % (atm_dtem_dh['value'], atm_dtem_dh['unit'], atm_dtem_dh_set))
                    _msg(" h0               |%-18s [%5s]   | %s " % (atm_h0['value'], atm_h0['unit'], atm_h0_set))
                    _msg(" altitude         |%-18s [%5s]   | %s " % (atm_altitude['value'], atm_altitude['unit'], atm_altitude_set))
                    _msg(" temperature      |%-18s [%5s]   | %s " % (atm_temperature['value'], atm_temperature['unit'], atm_temperature_set))
                    _msg(" pressure         |%-18s [%5s]   | %s " % (atm_pressure['value'], atm_pressure['unit'], atm_pressure_set))
                    _msg(" humidity         |%-18s [%5s]   | %s " % (atm_humidity['value'], atm_humidity['unit'], atm_humidity_set))
                    _msg(" pwv              |%-18s [%5s]   | %s " % (atm_pwv['value'], atm_pwv['unit'], atm_pwv_set))
                    _msg(" dp               |%-18s [%5s]   | %s " % (atm_dp['value'], atm_dp['unit'], atm_dp_set))
                    _msg(" dpm              |%-18s [%5s]   | %s " % (atm_dpm['value'], atm_dpm['unit'], atm_dpm_set))
                    _msg("*maxAltitude      |%-18s [%5s]   | (FIXED CONST) " % (atm_maxAltitude['value'], atm_maxAltitude['unit']))
                    _msg(" layerboundaries  |%-18s " % atm_layerboundaries)
                    _msg(" layertemperature |%-18s " % atm_layertemperature)
                    _msg("------------------+-----------------------------------------------")

                    ###################
                    # initATMProfile
                    ###################
                    atm = at.initAtmProfile(humidity=atm_humidity['value'],
                                            temperature=atm_temperature,
                                            altitude=atm_altitude,
                                            pressure=atm_pressure,
                                            atmType=atm_atmtype,
                                            maxAltitude=atm_maxAltitude,
                                            h0=atm_h0,
                                            dTem_dh=atm_dtem_dh,
                                            dP=atm_dp,
                                            dPm=atm_dpm['value'],
                                            layerBoundaries=atm_layerboundaries,
                                            layerTemperature=atm_layertemperature)

                    # show info determined by intAtmProfile
                    show_atm_info(atm)

                    #  Layer Information generated by initAtmProfile
                    #   - automatically shown, when the arg is specified)
                    if len(atm_layerboundaries) != 0:
                        show_layer_info(at)

                    ###################
                    # Spectral Window #
                    ###################
                    at.initSpectralWindow(nband,
                                          qa.quantity(fcenter, 'Hz'),
                                          qa.quantity(nchan * chansep, 'Hz'),
                                          qa.quantity(chansep, 'Hz'))

                    # H2O
                    at.setUserWH2O(atm_pwv)
                    prevtmatm = tt[-1]

                at.setAirMass(1 / pl.cos(pl.pi / 2. - eon))
                tskyon = at.getTrjSkySpec()[1]['value']

                at.setAirMass(1 / pl.cos(pl.pi / 2. - eoff))
                tskyoff = at.getTrjSkySpec()[1]['value']
                tau0 = at.getDryOpacitySpec()[1] + at.getWetOpacitySpec()[1]['value']
                tauoff = tau0 / pl.cos(pl.pi / 2. - eoff)

                dTa = (tskyon - tskyoff) * pl.exp(tauoff)
                if istdm:
                    dTa = pl.convolve(dTa, hanning, 'same')
                    if chanfreqs[spwid][0] < chanfreqs[spwid][1]:
                        dTa = pl.interp(chanfreqs[spwid], chanfreqs_high, dTa)
                    else:
                        dTa = pl.flipud(pl.interp(pl.flipud(chanfreqs[spwid]), pl.flipud(chanfreqs_high), pl.flipud(dTa)))
                elif dosmooth:
                    dTa = pl.convolve(dTa, [0.25, 0.5, 0.25], 'same')

                for ipol in range(npol):
                    cdata[ipol, :, i] -= dTa * factor

            subtb.putcol(p_datacolumn, cdata)
            subtb.close()

        # end for spwid in spws:
        tb.flush()

        # end of with statement

    return

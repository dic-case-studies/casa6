import os
import pylab as pl
import shutil
import contextlib
import collections
import numpy as np

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
    from casatools import quanta, table, msmetadata
    from casatools import atmosphere
    from casatasks import mstransform
    import casatasks.private.simutil as simutil

    ut = simutil.simutil()
    msmd = msmetadata()
    tb = table()
    qa = quanta()
    at = atmosphere()

else:
    from taskinit import tbtool, casalog, qa
    from taskinit import msmdtool as msmetadata
    from casac import casac
    from tasks import mstransform
    from simutil import simutil

    ut = simutil()
    msmd = msmetadata()
    tb = tbtool()
    at = casac.atmosphere()

# Task name #
origin = 'sdatmcor'

# Table/MS safe access #
@contextlib.contextmanager
def open_table(path, nomodify=True):
    # tb = table( )
    tb.open(path, nomodify=nomodify)
    # print("tentative::TB OPENED.")
    try:
        yield tb
    finally:
        # print("tentative::TB CLOSED.")
        tb.close()


# MS METADATA  safe access #
@contextlib.contextmanager
def open_msmd(path):
    msmd.open(path)
    # print("tentative::MSMD OPENED.")
    try:
        yield msmd
    finally:
        # print("tentative::MSMD CLOSED.")
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
                assert False
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


def sdatmcor(
        infile, datacolumn, outfile, overwrite,
        field, spw, scan, antenna,
        correlation, timerange, intent,
        observation, feed, msselect,
        outputspw,
        gainfactor,
        dtem_dh, h0, atmtype,
        atmdetail,
        altitude, temperature, pressure, humidity, PWV,
        dp, dpm,
        layerboundaries, layertemperature):

    # Information #
    casalog.origin(origin)
    _msg("\nSDATMCOR revision 1106-SiVasParaChile (06-Nov-2020) .\n")

    #
    # Input/Output error check and internal set up.
    #
    if infile == '':
        errmsg = "infile MUST BE specified."
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)
    if outfile == '':
        errmsg = "outfile MUST BE specified."
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)       
    # Protection. In case infile == outfile #
    if infile == outfile:
        errmsg = "You are attempting to write the output on your input file."
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)

    # File Info. #
    _msg("INPUT/OUTPUT")
    _msg("  Input MS file   = %s " % infile)
    _msg("  Output MS file  = %s " % outfile)

    # Existence #
    infile_exist = _file_exist(infile)
    outfile_exist = _file_exist(outfile)

    # infile Inaccessible #
    if not infile_exist:
        errmsg = "Specified infile does not exist."
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)

    # outfile Protected #
    if outfile_exist:
        if overwrite:
            _msg("Overwrite:: Overwrite specified. Once delete the existing output file. ")
            _ms_remove(outfile)
        else:
            errmsg = "Specified outfile already exist."
            _msg("\nERROR::%s\n" % errmsg, 'ERROR')
            raise Exception(errmsg)

    # Antenna Key word, form without &&& # 
    if antenna != '':
        antenna = antenna.replace('&','')

    #
    # Inspect arguments.
    #  - inspect Unit.
    #  - default value ARE NOT considered here.
    #  - convert values to string form.
    #  - Sub-function calc_sdatmcor() accepts args basically by string.
    #
    dtem_dh     = _check_unit_and_formToStr(dtem_dh, ['K/km'])
    h0          = _check_unit_and_formToStr(h0, ['km'])

    altitude    = _check_unit_and_formToStr(altitude, ['m'])
    temperature = _check_unit_and_formToStr(temperature, ['K'])
    pressure    = _check_unit_and_formToStr(pressure, ['mbar', 'hPa'])
    humidity    = _check_unit_and_formToStr(humidity, [''])  # through (string or float)
    PWV         = _check_unit_and_formToStr(PWV, ['mm'])
    dp          = _check_unit_and_formToStr(dp,  ['mbar', 'hPa'])
    dpm         = _check_unit_and_formToStr(dpm, [''])       # through (string or float)

    # Inspect atmtype ('str or int'). The range is checked and accept atmtype==''  #
    if not _inspect_strint_range(atmtype, 1, 5):
        errmsg = "atmtype (=%s) Out of Range or Unacceptable." % atmtype
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)

    #
    # User-Defined Profile inspection.
    # after this step, the two args changes to List
    #

    #   User-Defined-Profile parameters  conversion. An empty arg makes [] list. #
    layerboundaries = _convert_userdefinedparam_to_list(layerboundaries)
    layertemperature = _convert_userdefinedparam_to_list(layertemperature)

    #   Length of the two args must be same #
    len_1 = len(layerboundaries)
    len_2 = len(layertemperature)
    if len_1 != len_2:
        errmsg = "Data count mismatches in specified User-Defined parameter. len=[%d, %d] \n" % (len_1, len_2)
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)

    # generate gain factor dictionary
    gaindict = parse_gainfactor(gainfactor)

    #
    # Call calc Function
    #
    return calc_sdatmcor(
        infile, datacolumn, outfile, overwrite,
        field, spw, scan, antenna, correlation, timerange, intent, observation, feed, msselect,
        outputspw,
        gaindict,
        dtem_dh, h0, atmtype,
        atmdetail,
        altitude, temperature, pressure, humidity, PWV, dp, dpm,
        layerboundaries,
        layertemperature)

#
# SUBROUTINES
#  for Task Handling
#

def _ms_remove(path):
    if (os.path.exists(path)):
        if (os.path.isdir(path)):
            shutil.rmtree(path)
        else:
            os.remove(path)


def _ms_copy(src, dst):
    shutil.copytree(src, dst)


def _file_exist(path):
    return os.path.exists(path)


#
# Unit handling service
#

# inspect the input value is in the specified range #
def _inspect_strint_range(indata, minimum, maximum):
    if type(indata) is str:
        if indata == '':
            return True
        elif indata.isdigit():
            return minimum <= int(indata) <= maximum
        else:
            return False
    elif type(indata) is int:
        return minimum <= indata <= maximum
    else:
        # INTERNAL ERROR:: unexpected data type. #
        assert(False)

# inspect the input value is consistent with the Unit. #
def _check_unit_and_formToStr(data, base_unit):
    try:
        if type(data) is str:
            if (data == ''):
                # No data #
                return ''
            ext_unit = qa.getunit(data)
            if (ext_unit in base_unit):
                # With Unit #
                # _msg("Unit Conversion::Data with Unit '%s'" % data)
                return str(qa.getvalue(data)[0])
            elif (ext_unit == ''):
                # Without Unit and added  #
                # _msg("Unit Conversion::No unit specified in %s . Assumed '%s'" % (data, base_unit))
                return data
            else:
                # Mismatch (ERROR) #
                errmsg = "Unit conversion:: Unexpected Unit '%s' in %s ." % (ext_unit, data)
                _msg("ERROR::%s" % errmsg, 'ERROR')
                raise Exception(errmsg)
        elif (type(data) is int) or (type(data) is float): 
            if data == -1:
                # float no input #
                return ''
            else:
                # float specified #
                return str(data)   # available  input#
        else:
            # INTERNAL ERROR:: Arg type is not expected due to the I/F Design. #
            assert(False)
    except Exception as err:
        casalog.post('%s' % err, 'SEVERE')
        raise Exception("internal function error.")

#
# Argument parameter handling
# type : int, list, float 
#
# (action) check in_para and returns chosen value either default or argument.
#  if arg_value is available , return atm_parm_variable with being converted.
#  otherwise, returns def_para to use as a default parameter.
#
def _set_int_atmparam_from_args(arg_value, atm_parm_variable):
    if (arg_value != ''):
        return int(arg_value)
    else:
        return  atm_parm_variable


def _set_list_atmparam_from_args(set_list):
    if (len(set_list) != 0):
        return set_list
    else:
        return []


def _set_float_atmparam_from_args(arg_value, atm_parm_variable, unit):
    if type(arg_value) is str:
        if (arg_value != ''):
            return qa.quantity(float(arg_value), unit), True
        else:
            return atm_parm_variable, False

# UNDER REVISION #
def _make_list_from_separatedstring(separated_string, dType):
    # No input #
    if separated_string == '':
        return []

    try:
        if type(separated_string) is str:
            tmp_list = separated_string.split(',')
            if dType == 'str':
                out_list = [str(s) for s in tmp_list]
            elif dType == 'int':
                out_list = [int(s) for s in tmp_list]
            else:
                out_list = [s for s in tmp_list]
            return out_list
        elif type(separated_string) is list:
            return separated_string
        else:
            return []
    except Exception as err:
        _msg("Error in comma-separated string.")
        casalog.post('%s' % err, 'SEVERE')
        raise Exception("internal function error.")

# UNDER REVISION #
def _convert_userdefinedparam_to_list(in_arg):
    # No input #
    if (type(in_arg) is str) and (in_arg == ''):
        return [] 

    try:
        if  type(in_arg) is list:
            out_list = [float(s) for s in in_arg]
            return out_list
        elif type(in_arg) is str:
            tmp_list = in_arg.split(',')
            out_list = [float(s) for s in tmp_list]
            return out_list
        else:
            _msg("\nERROR::Invalid arg type. Expecting only separated string or list.\n", 'SEVERE')
            raise Exception("internal function error.")

    except Exception as err:
        _msg("Error in converting an element in the List.")
        casalog.post('%s' % err, 'SEVERE')
        raise Exception("internal function error.")

#
# decide Default Antenna ID
#
def get_default_antenna(msname, antenna):
    # Choose base-antenna from selected Antenna ID
    #  - Search Priority ID is   1 > 2 > 3 > x
    #  - if ONLY one  antenna(=x) is available, use this.
    # The Rule is defined in CASR-552
    #  - iAnt = 1 is experimentally preferable. 

    with open_msmd(msname) as msmd:
        ant_list = msmd.antennaids(antenna)
        n_ant = len(ant_list)

        # Choose One #
        if 1 in ant_list:
            i_ant = 1
        elif 2 in ant_list:
            i_ant = 2
        elif 3 in ant_list:
            i_ant = 3
        elif n_ant == 1:
            i_ant = ant_list[0]
        else:
            errmsg="Irregular antenna ID detected."
            _msg("\nERROR::%s\n" % errmsg, 'ERROR')
            # raise Exception(errmsg)

        # INFO #
        ant_name = msmd.antennanames(i_ant)[0]

        _msg("Default Antenna")
        _msg(" - totally %d antenna(s) were picked up. [query=%s]" % (n_ant, antenna))
        _msg(" - Antenna ID  = %d was chosen. Name= %s" % (i_ant, ant_name))

    return i_ant

#
# decide default value of 'Altitude'
# - This requires to calculate Elevation from Antenna Position Information.
#
def get_default_altitude(msname, antid):
    with open_table(os.path.join(msname, 'ANTENNA')) as tb:
        # obtain the antenna Position (Earth Center) specified by antid #
        pos = tb.getcell('POSITION', antid)
        X = float(pos[0])
        Y = float(pos[1])
        Z = float(pos[2])

        # 
        #  xyz2long()   -- https://casa.nrao.edu/casadocs/casa-5.6.0/simulation/simutil
        #
        #  When given ITRF Earth-centered (X, Y, Z, using the parameters x, y, and z) coordinates [m] for a point, 
        #  this method returns geodetic latitude and longitude [radians] and elevation [m]. 
        #  Elevation is measured relative to the closest point to the (latitude, longitude) 
        #  on the WGS84 (World Geodetic System 1984) reference ellipsoid.
        P = ut.xyz2long(X, Y, Z, 'WGS84')     #  P[0]=longitude, P[1]=latitude, P[2]=elevation (geodetic elevation)
        geodetic_elevation = P[2]

        # ref #
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

def showAtmInfo(atm):
    """
     Returned atm (from initAtmProfile) may have a different structure.
     In casa6, additional variable information is added in Dict type. 
    """
    version = at.getAtmVersion()
    _msg("\nAtomosphere Tool:: version = %s\n" % version)

    # ATM Profile #
    if is_CASA6:
        _msg(atm[0])
    else:
        _msg(atm)

def showLayerInfo(at):
    p = at.getProfile()

    #
    # From Casa Tool Documentation:
    #
    # returns a tuple of
    # 0 - string listing of layer values, and arrays of layer, 1 - thickness,
    # 2 - temperature, 3 - watermassdensity, 4 - water (number density),
    # 5 - pressure, 6 - O3 (number density), 7 - CO, 8 - N2O

#    for i in range(at.getNumLayers()):
#        Print atmospheric profile returned by at.getProfile():
#        Layer thickness (idx=1), Temperature (idx=2),
#        Number density of water vapor(idx=4), and Pressure (idx=5)
#        print(p[1]['value'][i],p[2]['value'][i],p[4]['value'][i],p[5]['value'][i])

    _msg(p[0])
    return


#
# Logging CASA LOG (INFO/WARN/SEVERE)
#
def _msg(msg, msgtype='INFO'):
    # Information message #
    if msgtype == 'INFO':
        print(msg)
    # other Warning/Error message #
    casalog.post(msg, msgtype, origin=origin)


#
# Data Selection
# use mstransform,  msselect being included.
#
def atmMst(
    infile, datacolumn, outfile, overwrite,
    field, spw, scan, antenna,
    correlation, timerange, intent,
    observation, feed, msselect):

    # Antenna Key word, to adjust to mstransform # 
    if antenna != '':
        antenna = antenna.replace('&','')
        antenna += '&&&'

    mstransform(
        vis=infile,
        outputvis=outfile,
        datacolumn=datacolumn,
        field=field,
        spw=spw,
        scan=scan,
        antenna=antenna,
        correlation=correlation,
        timerange=timerange,
        intent=intent,
        observation=observation,
        feed=feed,
        reindex=False)   # Must be False #

############################################################
# Calculation Method (Replaced to C++ in next development.)
#    originated by atmcorr_20200602.py by Sawada san.
#    formed as a casa task by CAS-13160.
#    - p_xxxx arguments are for task.
#    - param_xxxx arguments are for Atm-correction.
############################################################
def calc_sdatmcor(
        p_infile,
        p_datacolumn,
        p_outfile,
        p_overwrite,
        p_field,
        p_spw,
        p_scan,
        p_antenna,
        p_correlation,
        p_timerange,
        p_intent,
        p_observation,
        p_feed,
        p_msselect,
        p_outputspw,
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

    #
    # Argument dump.
    #   if need, set True.
    #
    if True:
        _msg("***********************************")
        _msg("**   calc_sdatmcor::             **")
        _msg("***********************************")
        _msg('infile      = %s' % p_infile)
        _msg('datacolumn  = %s' % p_datacolumn)
        _msg('outfile     = %s' % p_outfile)
        _msg('overwrite   = %s' % p_overwrite)
        _msg('field       = %s' % p_field)
        _msg('spw         = %s' % p_spw)
        _msg('scan        = %s' % p_scan)
        _msg('antenna     = %s' % p_antenna)
        _msg('correlation = %s' % p_correlation)
        _msg('timerange   = %s' % p_timerange)
        _msg('intent      = %s' % p_intent)
        _msg('observation = %s' % p_observation)
        _msg('feed        = %s' % p_feed)
        _msg('msselect    = %s' % p_msselect)
        _msg('outputspw   = %s' % p_outputspw)
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

    # Debug flags. #
    showCorrection = False        # show index information while Correction.
    interruptCorrection = False   # Interrupt Correction
    interruptCorrectionCnt = 200  # (limit count)

    # datacolumn (XML fills default) ,to UPPER CASE #
    datacolumn = p_datacolumn.upper()
    if (datacolumn == 'CORRECTED'):    # 'CORRECTED' means column:'CORRECTED_DATA'
        datacolumn = 'CORRECTED_DATA'

    # Internal file names  #
    rawms = p_infile
    calms = p_outfile
    corms = p_outfile

    #
    # CAS-13160
    # Original Script starts from here
    #  - in many sections, Task Code is inserted.
    #  - comment is generally added to distinguish the original and the task code.
    #

    ##################################################
    #   Inside Constant for ATM
    ##################################################

    _msg("Initializing Tool Constant.")

    #
    # Following variables have initial parameters according to the Original script.
    #  - these values will be translated with unit to var.:'atm_xxxxxxxxx' in the Task code.
    #

    # Default constant. Refer to CAS-13160 #

    atmtype = 2         ### atmType parameter for at (1: tropical, 2: mid lat summer, 3: mid lat winter, etc)
    maxalt = 120        ### maxAltitude parameter for at (km)
    lapserate = -5.6    ### dTem_dh parameter for at (lapse rate; K/km)
    scaleht = 2.0       ### h0 parameter for at (water scale height; km)

    dosmooth = False    ### convolve dTa* spectra with [0.25, 0.5, 0.25] to mimic Hanning spectral response;
                        ### set to True if spectral averaging was not employed for the spw
    dp = 10.0   ### initATMProfile DEFAULT ###
    dpm = 1.2   ### initATMProfile DEFAULT ###


    # other #
    nband = 1         # number of band in initSpectralWindow()

    # normalized spectral response for Hanning window, FWHM=10
    hanning = [-0.00098041, -0.00202866, -0.00265951, -0.00222265,
               0.00000000,  0.00465696,  0.01217214,  0.02260546,  0.03556241,
               0.05017949,  0.06519771,  0.07911911,  0.09042184,  0.09779662,
               0.10035898,  0.09779662,  0.09042184,  0.07911911,  0.06519771,
               0.05017949,  0.03556241,  0.02260546,  0.01217214,  0.00465696,
               0.00000000, -0.00222265, -0.00265951, -0.00202866, -0.00098041]

    # Data Selection
    _msg("Data Selection in progress. Output file =  %s " % corms)
    try:
        atmMst(
            infile=rawms,
            datacolumn=p_datacolumn,
            outfile=corms,      
            overwrite=p_overwrite,
            field=p_field,
            spw = p_outputspw,
            scan=p_scan,
            antenna=p_antenna,
            correlation=p_correlation,
            timerange=p_timerange,
            intent=p_intent,
            observation=p_observation,
            feed=p_feed,
            msselect=p_msselect)

    except Exception as err:
        casalog.post('%s' % err, 'SEVERE')
        errmsg = "Something is wrong in atmMst. "
        raise Exception(errmsg)

    # Resume 'origin'. A strange behavior in casalog/CASA6 #
    casalog.origin(origin)

    # Result check if output exists. #
    if not _file_exist(corms):
        errmsg = "No outfile has been generated by mstransform."
        raise Exception(errmsg)

    #
    # Tool Default constant (CAS-13160)
    # - Antenna and its altitude are determined.
    #
    antenna = get_default_antenna(rawms, p_antenna)       # default antenna_id   (UNDER RE-CONSTRUCTION)
    altitude = get_default_altitude(rawms, antenna)       # default altitude  - see inside in detail (UNDER RE-CONSTRUCTION)

    #
    # Use 'DATA' after  mstransform.
    #  only when CORRECTED_DATA is originally specified.  
    #
    if datacolumn == 'CORRECTED_DATA':
        datacolumn = 'DATA'

    #
    # CAS-13160
    #
    #   From here, main part of the Original script starts.
    #   In some sections, additional statements for Task are inserted.
    #

    ################################################################
    # Get metadata
    ################################################################
    chanfreqs = {}
    _msg("\nSDATMCOR main body starts. rawms=%s\n" % rawms)

    try:
        with open_msmd(rawms) as msmd:
            # (original)
            tmonsource = msmd.timesforintent('OBSERVE_TARGET#ON_SOURCE')
            tmoffsource = msmd.timesforintent('OBSERVE_TARGET#OFF_SOURCE')
            fdmspws = msmd.fdmspws()
            tdmspws = msmd.tdmspws()
            intentspws = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
            spws = list(set(intentspws) & (set(fdmspws) | set(tdmspws)))
            spwnames = msmd.namesforspws(spws)

            # show to user. #
            ## _msg(" - tm ON  source  = %s" % tmonsource)
            ## _msg(" - tm OFF source  = %s" % tmoffsource)
            _msg(" - fdm    spws = %s" % fdmspws)
            _msg(" - tdm    spws = %s" % tdmspws)
            _msg(" - intent spws = %s" % intentspws)
            _msg(" -        spws = %s" % spws)

            # save spws (CAS-13160) #
            rawmsSpws = spws

            # (Task Section) check count of ON/OFF SOURCE
            n_tmonsource = len(tmonsource)
            n_tmoffsource = len(tmoffsource)
            msg = "Target Information. \n"   \
                + "# ON_SOURCE: count of tmonsource   = %d\n" % n_tmonsource  \
                + "# OFF_SOURCE: count of tmoffsource = %d" % n_tmoffsource
            _msg(msg)

            # (Task Section) OFF_SOURCE check #
            if (n_tmoffsource == 0):
                msg = "Can't find the OFF_SOURCE data."
                _msg(msg, 'SEVERE')
                raise(msg)

            # (Task Section) ON_SOURCE check #
            if (n_tmoffsource == 0):
                msg = "Can't find the ON_SOURCE data."
                _msg(msg, 'SEVERE')
                raise(msg)

            #
            # (Task Section)
            #     'OutputSpw' 
            #     must be a set of rawmsSpw
            #

            # request by argument  #
            outputspws_param = _make_list_from_separatedstring(p_outputspw, dType='int')

            # Must be a subset, locate the initial set.  #
            if set(rawmsSpws) >= set(outputspws_param):
                outputspws = list(set(outputspws_param))
            else:
                _msg("Some of the specified outputspw(s) cannot be processed. Try to continue", 'WARN')
                outputspws = list(set(rawmsSpws) & set(outputspws_param))


            # If default, apply all-spws #
            if (p_outputspw == ''):
                outputspws = rawmsSpws
     
            _msg("Determined outputSpws Information")
            _msg('- rawms      Spws       = %s' % rawmsSpws)
            _msg('- requested  outputSpws = %s' % outputspws_param)
            _msg('- determined outputSpws = %s' % outputspws)

 
            #
            # (Task Section )
            #     'processing Spw'
            #      must be a set of rawmsSpw
            #

            # request by argument  #
            spws_param = _make_list_from_separatedstring(p_spw, dType='int')

            # Must be a subset, locate the initial set.  #
            if set(rawmsSpws) >= set(spws_param):
                spws = list(set(spws_param))
            else:
                _msg("Some of the specified spw(s) cannot be processed. Try to continue", 'WARN')
                spws = list(set(rawmsSpws) & set(spws_param))

            # If default, apply all-spws #
            if (p_spw == ''):
                spws = rawmsSpws

            _msg("Determined Spws Information")
            _msg('- rawms      Spws = %s' % rawmsSpws)
            _msg('- requested  Spws = %s' % spws_param)
            _msg('- determined Spws = %s' % spws)

            #
            # (Task Section )
            #     outputSpw and Spw Consistency. For example;
            #      - reject  when Spw=[17,19,21], outputSpw=[19,21]
            #      - accept  when Spw=[21], outputSpw = [19,21]
            noCorSpws = []
            if set(outputspws) >= set(spws):
                noCorSpws = list(set(outputspws) - set(spws))
            else:
                spws= list(set(outputspws) & set(spws)) 

            # outputSpw, Spw:: No Target check #
            if spws == []:
                raise Exception("No Spw Targets. Abort.")
            if outputspws == []:
                raise Exception("No outputSpw targets. Abort.")

            #
            # (Original Section)
            #
 
            # CAS-13160 Changed #
            for spwid in rawmsSpws:
            ## for spwid in spws:
                chanfreqs[spwid] = msmd.chanfreqs(spw=spwid)

            # end of with

    except Exception as err:
        casalog.post('%s' % err, 'ERROR')
        errmsg = "Error in opening rawms"
        raise Exception(errmsg)

    # (original)
    bnd = (pl.diff(tmoffsource) > 1)
    w1 = pl.append([True], bnd)
    w2 = pl.append(bnd, [True])
    tmoffsource = (tmoffsource[w1]+tmoffsource[w2])/2.  ### midpoint of OFF subscan

    ddis = {}
    with open_msmd(calms) as msmd:
        # CAS-13160 Changed #
        for spwid in rawmsSpws:
        ## for spwid in spws:
            ddis[spwid] = msmd.datadescids(spw=spwid)[0]
    print("- ddis[] = %s" % ddis)

    nchanperbb = [0, 0, 0, 0]
    bbprs = {}

    # CAS-13160 Changed #
    for i, spwid in enumerate(rawmsSpws):
    ## for i, spwid in enumerate(spws):
        bbp = int(spwnames[i].split('#')[2][3])-1
        bbprs[spwid] = bbp
        nchanperbb[bbp] += len(chanfreqs[spwid])

    # (Original)
    # Data Query on Pointing Table.
    #   output::
    #   - tmpointing, elev
    #
    try:
        with open_table(os.path.join(calms, 'POINTING')) as tb:
            # (org.) key for ANTENNA_ID select
            querytext = 'ANTENNA_ID==%s' % antenna
            subtb = tb.query(querytext)

            # (org.) Access Table
            tmpointing = subtb.getcol('TIME')
            elev = subtb.getcol('DIRECTION').squeeze()[1]
            _msg("- reading column:'TIME' and 'DIRECITON' completed.")

            subtb.close()

    except Exception as err:     
        casalog.post('%s' % err, 'ERROR')
        errmsg = "Error in opening POINTING."
        raise Exception(errmsg)

    ################################################################
    # Get atmospheric parameters for ATM
    ################################################################
    try:
        with open_table(os.path.join(rawms, 'ASDM_CALWVR')) as tb:
            # confirm
            _msg("tmonsource: %f, %f" % (tmonsource.min(), tmonsource.max()))
            pwv = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                           (tmonsource.min(), tmonsource.max())).getcol('water')

    except Exception as err:
        casalog.post('%s' % err, 'ERROR')
        errmsg = "Error in opening rawms/'ASDM_CALWVR'."
        raise Exception(errmsg)

    # pick up pwv
    pwv = pl.median(pwv)

    # ASDM_CALATMOSPHERE
    try:
        with open_table(os.path.join(rawms, 'ASDM_CALATMOSPHERE')) as tb:
            subtb = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                             (tmonsource.min(),
                              tmonsource.max()))

            tground = pl.median(subtb.getcol('groundTemperature'))
            pground = pl.median(subtb.getcol('groundPressure'))
            hground = pl.median(subtb.getcol('groundRelHumidity'))

            subtb.close()
    except Exception as err:
        casalog.post('%s' % err, 'ERROR')
        errmsg = "Error in opening rawms/'ASDM_CALATMOSPHERE'."
        raise Exception(errmsg)

    _msg("median PWV = %fm, T = %fK, P = %fPa, H = %f%%" % (pwv, tground, pground, hground))

    ################################################################
    # Looping over spws
    ################################################################
    with open_table(corms, nomodify=False) as tb:
        #
        # Note CAS-13160:
        #     spw for-loop:
        #       - Inside this loop, calculation is executed with each spwid.
        #       - The corrected result is written upon OutputFile with selected Spw.
        #       - OutputSpw occasionally intend No-Corrected value.
        #       - please see the inserted block.
        #        - No-Corrected output is controled by 'noCorSpws' list.
        # 
        #      The original script does not explicitly imply the logic.
        #      Such No-corrected output was appended with minimum changes.
  
        for spwid in outputspws: # this was originally 'in spws' #
            # Log #
            _msg("\nProcessing spw %d in %s. Correction is not applied to spw %s .\n" % (spwid, outputspws, noCorSpws))

            # gain factor
            spwkey = str(spwid)
            factor = gaindict[spwkey]
            #_msg('Applying gain factor {} to spw {}'.format(spwkey, factor))
            _msg('Applying gain factor {} to spw {}'.format(factor, spwkey))

            istdm = False
            nchan = len(chanfreqs[spwid])
            # fcenter = (chanfreqs[spwid][nchan/2-1]+chanfreqs[spwid][nchan/2])/2.             # PY2
            fcenter = (chanfreqs[spwid][int(nchan/2)-1]+chanfreqs[spwid][int(nchan/2)])/2.   # PY3
            chansep = (chanfreqs[spwid][-1]-chanfreqs[spwid][0])/(nchan-1)

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
                fcenter = (chanfreqs_high[int(nchan/2)-1]+chanfreqs_high[int(nchan/2)])/2.    # PY3

            ################################################################
            # Set parameters for ATM and obtain zenith opacity
            # caveat: median values are used for calibrating the entire execution
            ################################################################

            _msg("- set parameters for initATM and obtain zenith opacity")

            # Note CAS-13160: Revised part from original.
            #
            #    Setting up parameter variables to initAtmProfile()
            #
            #    - Following codes are extracted from the original script, which used to be
            #      embeded inside the initAtmProfile() arguments.
            #    - Overwrite the  'parameter' to 'atm' arguments of initAtmProfile.
            #    - The default values are up to Task Specification rather than original script.

            atm_dtem_dh     = qa.quantity(lapserate, 'K/km')
            atm_h0          = qa.quantity(scaleht, 'km')
            atm_atmtype     = atmtype
            atm_atmdetail   = atmdetail  # Bool
            atm_altitude    = qa.quantity(altitude, 'm')
            atm_temperature = qa.quantity(tground, 'K')              # tground (original)
            atm_pressure    = qa.quantity(pground/100.0, 'mbar')     # pground (original) in  [Pa]  convert to [mbar]
            atm_humidity    = qa.quantity(hground, '%')              # hground (original) in  [%]
            atm_pwv         = qa.quantity(pwv * 1000.0, 'mm')        # pwv (original) in [m] converto to [mm]
            atm_dp          = qa.quantity(dp, 'mbar')
            atm_dpm         = qa.quantity(dpm, '')                   # keep float
            atm_maxAltitude = qa.quantity(maxalt, 'km')

            # Edit Flag (True: the parameter is given and will be applied.) #
            atm_dtem_dh_set  = False
            atm_h0_set       = False
            atm_altitude_set = False
            atm_temperature_set  = False
            atm_pressure_set     = False
            atm_humidity_set     = False
            atm_pwv_set      = False
            atm_dp_set       = False
            atm_dpm_set      = False

            #
            # (from 'help' information)
            # User-Defined Profile (example)
            #    myalt = [ 5071.72200397, 6792.36546384, 15727.0776121, 42464.18192672 ] #meter
            #    mytemp = [ 270., 264., 258., 252. ] #Kelvin
            #
            atm_layerboundaries  = []   # initially null-List.   array: layer boundary [m]
            atm_layertemperature = []   # initially null-list.   array: layer temperature [K]

            #
            # ATM fundamental parameters activation.
            #
            atm_dtem_dh, atm_dtem_dh_set = _set_float_atmparam_from_args(param_dtem_dh, atm_dtem_dh, 'K/km')
            atm_h0,      atm_h0_set      = _set_float_atmparam_from_args(param_h0, atm_h0, 'km')
            atm_atmtype   = _set_int_atmparam_from_args(param_atmtype, atm_atmtype)

            #
            # Sub parameter activation
            #  only when atmdetail is True
            #
            if atm_atmdetail:
                _msg("\nSub Parameter from the Arguments will be set, if specified.\n")
                atm_altitude,    atm_altitude_set     = _set_float_atmparam_from_args(param_altitude, atm_altitude, 'm')
                atm_temperature, atm_temperature_set  = _set_float_atmparam_from_args(param_temperature, atm_temperature, 'K')
                atm_pressure,    atm_pressure_set     = _set_float_atmparam_from_args(param_pressure, atm_pressure, 'mbar')
                atm_humidity,    atm_humidity_set     = _set_float_atmparam_from_args(param_humidity, atm_humidity, '%')
                atm_pwv,         atm_pwv_set          = _set_float_atmparam_from_args(param_PWV, atm_pwv, 'mm')
                atm_dp,          atm_dp_set           = _set_float_atmparam_from_args(param_dp, atm_dp, 'mbar')
                atm_dpm,         atm_dpm_set          = _set_float_atmparam_from_args(param_dpm, atm_dpm, '')

                # User-Defined Profile. 
                #  when the arg is active, directly pass the arg(param_xxxx) to atm_xxxx for initAtmProfile. 
                atm_layerboundaries  = _set_list_atmparam_from_args(param_layerboundaries)
                atm_layertemperature = _set_list_atmparam_from_args(param_layertemperature)

            else:
                _msg("\nSub Parameters were not used, due to 'atmdetail' is not True.\n")

            #
            # print and log, to confirm Atm parameters.
            #
            _msg("====================================================================")
            _msg("  initATMProfile Parameters     [atmdetail = %s]" % atm_atmdetail)
            _msg("------------------+-----------------------------+-------------------")
            _msg("  parameter       | value [unit]                | specified by arg. ")
            _msg("------------------+-----------------------------+-------------------")
            _msg(" atmtype          |%-18s " % at.listAtmosphereTypes()[atm_atmtype-1])    # type =1,2,3,4,5, no Zero
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

            #
            # show intAtmProfile result
            #
            showAtmInfo(atm)

            #
            #  Layer Information
            #   - automatically showing, when the arg is specified)
            #
            if len(atm_layerboundaries) != 0:
                showLayerInfo(at)
     
            ###################
            # Spectral Window #
            ###################
            at.initSpectralWindow(nband,
                                  qa.quantity(fcenter, 'Hz'),
                                  qa.quantity(nchan*chansep, 'Hz'),
                                  qa.quantity(chansep, 'Hz'))

            # H2O
            at.setUserWH2O(atm_pwv)

            ################################################################
            # Calculate and apply correction values
             ################################################################

            _msg("\nCalculate and apply correction values.\n")

            # original: make ddis[spwid]
            _msg("- Selecting DATA_DESC_ID == %s" % ddis[spwid])

            # Query Text (DESC_ID) #
            querytext = 'DATA_DESC_ID in %s' % ddis[spwid] 

            # Applying MSSELECT if anything is given #
            if (p_msselect != ''):
                querytext += p_msselect

            # Query # 
            subtb = tb.query(querytext)
  
            # Data #
            _msg("- getting tm and data. datacolumn [%s] is used." % datacolumn)
            tmdata = subtb.getcol('TIME')
            data = subtb.getcol(datacolumn)
            npol = data.shape[0]
      
            # Smoothing control
            if nchanperbb[bbprs[spwid]]*npol in [256, 8192]:
                _msg('Spw %d in BB_%d (total Nchan within BB is %d, sp avg likely not applied).  dosmooth=True' %
                     (spwid, bbprs[spwid]+1, nchanperbb[bbprs[spwid]]*npol))
                dosmooth = True
            else:
                _msg('Spw %d in BB_%d (total Nchan within BB is %d, sp avg likely applied).  dosmooth=False' %
                     (spwid, bbprs[spwid]+1, nchanperbb[bbprs[spwid]]*npol))
                dosmooth = False

            ###########################
            # Correction Main Loop
            ###########################

            _msg("\nExecuting ATM Correction(N=%d), and writing to output MS. \n" % len(tmdata))

            cdata = data.copy()
            for i, t in enumerate(tmdata):
                # debug option, interrupt the correction loop.
                if interruptCorrection:
                    if i > interruptCorrectionCnt:
                        msg = "Correction Loop was interrupted. \n---"
                        _msg(msg)
                        break

                # (original script) #
                dt = tmoffsource - t
                if (dt < 0).sum() == 0 or (dt > 0).sum() == 0:  ### any data before first OFF or after last OFF are disregarded
                    continue

                dtoff0 = dt[dt < 0].max()
                dtoff1 = dt[dt > 0].min()
                eon, eoff0, eoff1 = pl.interp([t, t+dtoff0, t+dtoff1], tmpointing, elev)
                eoff = (dtoff1*eoff0 - dtoff0*eoff1)/(dtoff1-dtoff0)

                at.setAirMass(1/pl.cos(pl.pi/2.-eon))
                tskyon = at.getTrjSkySpec()[1]['value']

                at.setAirMass(1/pl.cos(pl.pi/2.-eoff))
                tskyoff = at.getTrjSkySpec()[1]['value']
                tau0 = at.getDryOpacitySpec()[1] + at.getWetOpacitySpec()[1]['value']
                tauoff = tau0/pl.cos(pl.pi/2.-eoff)

                dTa = (tskyon-tskyoff) * pl.exp(tauoff)
                if istdm:
                    dTa = pl.convolve(dTa, hanning, 'same')
                    dTa = pl.interp(chanfreqs[spwid], chanfreqs_high, dTa)
                elif dosmooth:
                    dTa = pl.convolve(dTa, [0.25, 0.5, 0.25], 'same')

                #
                # debug option for task.
                #
                if showCorrection:
                    if spwid not in noCorSpws:
                        _msg("spw=%3d, i=%5d, time=%15s, Max(dTa)=%19s, Min(dTa)=%19s" % (spwid, i, t, max(dTa), min(dTa)))
                    else:
                        _msg("attempt to write through. spw=%3d, i=%5d" % (spwid, i))

                # CAS-13160 changed.
                # Adjust Body 
                #  - apply (the spwid is in both outputspws and spws) 
                #  - no apply (the spwid is in outputspws, but not in spws )
                #    This may happen in ordinary use-case.
                if spwid not in noCorSpws:
                    for ipol in range(npol):
                        cdata[ipol, :, i] -= dTa * factor

            subtb.putcol(datacolumn, cdata)
            subtb.close()

        # end for spwid in spws:
        tb.flush()

        # end of with statement

    return True

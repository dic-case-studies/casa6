import os
import pylab as pl
import shutil

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

def sdatmcor(
        infile, datacolumn, outfile, overwrite,
        field, spw, scan, antenna,
        correlation, timerange, intent,
        observation, feed, msselect,
        outputspw,
        dtem_dh, h0, atmtype,
        atmdetail,
        altitude, temperature, pressure, humidity, PWV,
        dp, dpm,
        layerboundaries, layertemperature ):

    # Information #
    casalog.origin(origin)
    _msg("\nSDATMCOR revision 1019 (19-Oct-2020) .\n")

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
#
# Inspect arguments.
#  - inspect Unit.
#  - default value ARE NOT considered here.
#  - convert values to string form.
#  - Subfunction calc_sdatmcor() accepts args basically by string.
#
    dtem_dh     = _inspect_value_unit(dtem_dh, ['K/km'])
    h0          = _inspect_value_unit(h0, ['km'])

    altitude    = _inspect_value_unit(altitude, ['m'])
    temperature = _inspect_value_unit(temperature, ['K'])
    pressure    = _inspect_value_unit(pressure, ['mbar', 'hPa'])
    humidity    = _inspect_value_unit(humidity, [''])  # through (string or float)
    PWV         = _inspect_value_unit(PWV, ['mm'])
    dp          = _inspect_value_unit(dp,  ['mbar', 'hPa'])
    dpm         = _inspect_value_unit(dpm, [''])       # through (string or float)

    # Inspect  atmtype ('str or int'). The range is checked and accept atmtype==''  #
    if not _inspect_str_int(atmtype, 1, 5):
        _msg("\nERROR:: atmtype (=%s) Out of Range or Unacceptable.\n" % atmtype, 'ERROR')
        return False

#
# User-Defined Profile inspection.
# after this step, the two args changes to List
#
    # Type conversion. An empty arg makes [] list. #
    layerboundaries = _conv_to_doubleArrayList(layerboundaries)
    layertemperature = _conv_to_doubleArrayList(layertemperature)

    # inspect counts, length of the two args must be same #
    len_1 = len(layerboundaries)
    len_2 = len(layertemperature)
    if len_1 != len_2:
        errmsg = "Data count mismatches in specified User-Defined parameter. len=[%d, %d] \n" % (len_1, len_2)
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)        

#
# Call calc Function
#
    return calc_sdatmcor(
        infile, datacolumn, outfile, overwrite,
        field, spw, scan, antenna, correlation, timerange, intent, observation, feed, msselect,
        outputspw,
        dtem_dh, h0, atmtype,
        atmdetail,
        altitude, temperature, pressure, humidity, PWV, dp, dpm,
        layerboundaries,
        layertemperature)

#
# SUBROUTINES
#  for Task Handling
#

#
# file handling (use Wrapper for readability)
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

# inspect the data can be translated to an Int #
def _inspect_str_int(data, minimum, maximum):
    if type(data) is str:
        if data == '':
            return True
        elif data.isdigit():
            return minimum <= int(data) <= maximum
        else:
            return False
    elif type(data) is int:
        return minimum <= data <= maximum
    else:
        # INTERNAL ERROR:: unexpected data type. #
        assert(False)

# inspect the data is consistent with the Unit. #
def _inspect_value_unit(data, base_unit):
    try:
        if type(data) is str:
            if (data == ''):
                # No data #
                return ''
            ext_unit = qa.getunit(data)
            if (ext_unit in base_unit):
                # With Unit #
                _msg("Unit Conversion::Data with Unit '%s'" % data)
                return str(qa.getvalue(data)[0])
            elif (ext_unit == ''):
                # Without Unit and added  #
                _msg("Unit Conversion::No unit specified in %s . Assumed '%s'" % (data, base_unit))
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
                # float specifiled #
                return str(data)   # available  input#
        else:
            # INTERNAL ERROR:: Arg type is not expected due to the I/F Design. #
            assert(False)
    except Exception as err:
        casalog.post('%s' % err, 'SEVERE')
        raise Exception("internal function error.")

#
# Argument parameter handling
#
# (action) check in_para and set default argument.
#  if in_para is available , return in_para with being converted.
#  otherwise, retrurns def_para to use as a default parameter.
#
def _set_int_ifactive(set_arg, through_value):
    if (set_arg != ''):
        return int(set_arg)
    else:
        return  through_value


def _copy_list_ifactive(set_list):
    if (len(set_list) != 0):
        return set_list
    else:
        return []


def _set_floatquantity_ifactive(set_arg, through_value, unit):
    if type(set_arg) is str:
        if (set_arg != ''):
            new_val = qa.quantity(float(set_arg), unit)  # CASA5 needs cast to float  ? #
            return new_val, True
        else:
            return through_value, False

def _list_comma_string(separated_string, dType):
    # make a list by separated by comma #
    try:
        if type(separated_string) is str:
            tmp_list = separated_string.split(',')
            if dType == 'str':
                out_list = [str(s) for s in tmp_list]  # convert to list['str','str', ...]
            elif dType == 'int':
                out_list = [int(s) for s in tmp_list]  # convert to list[int,int, ...]
            else:
                out_list = [s for s in tmp_list]  # No convert list [ data, data, ...]
            return out_list
        elif type(separated_string) is list:
            return separated_string
        else:
            return []
    except Exception as err:
        _msg("Error in comma-separated string.")
        casalog.post('%s' % err, 'SEVERE')
        return False


def _conv_to_doubleArrayList(in_list):
    # convert elements in  list or separated str =>  float list. #

    # check Empty #
    if (type(in_list) is str) and (in_list == ''):
        return [] 

    # conversion (non-decimal expression string will fail.) #
    try:
        if  type(in_list) is list:
            _msg("- converting a List which contains numerical expression or int/float to Float-List.")
            out_list = [float(s) for s in in_list]  # force to convert to list[floatm, ...]
            return out_list

        elif type(in_list) is str:
            _msg("- converting a Comma-separated string which contains numerical expression or int/float to Float-List.")
            tmp_list = in_list.split(',')  # convert to List #
            out_list = [float(s) for s in tmp_list]  # force to convert to list[float, ...]
            return out_list
        else:
            _msg("\nERROR::Invalid arg type, expecting separated string or list.\n", 'SEVERE')
            return []   # invalid type

    except Exception as err:
        _msg("Error in converting an element in the List.")
        casalog.post('%s' % err, 'SEVERE')
        return False

#
# decide Default Antenna ID
#
def get_default_antenna(msname, antenna):
    # Chose base-antenna from selected Antenna ID
    #  - Search Priority ID is   1 > 2 > 3 > x
    #  - if ONLY one  antenna(=x) is available, use this.
    # The Rule is defined in CASR discussion. 
    msmd.open(msname)
    ant_list = msmd.antennaids(antenna)
    n_ant = len(ant_list)
    msmd.close()

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
        errmsg="Illegular antenna ID detected."
        _msg("\nERROR::%s\n" % errmsg, 'ERROR')
        raise Exception(errmsg)

    # INFO #
    msmd.open(msname)
    ant_name = msmd.antennanames(i_ant)[0]
    msmd.close()

    _msg("Default Antenna")
    _msg(" - totally %d antenssas were picked up. [query=%s]" % (n_ant, antenna))
    _msg(" - Antenna ID  = %d was chosen. Name= %s" % (i_ant, ant_name))

    return i_ant

#
# decide default value of 'Altitude'
# - This requires to calculate Elevation from Antenna Position Information.
#
def get_default_altitude(msname, antid):
    tb.open(os.path.join(msname, 'ANTENNA'))
    # ref #
    ref = tb.getcolkeyword('POSITION', 'MEASINFO')['Ref']
    # obtain the antenna Position (Earth Center) spified by antid #
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
    P = ut.xyz2long(X, Y, Z, 'WGS84')     #  P[0]=longitude, P[1]=latitude, P[2]=elevation
    elev = P[2]
    tb.close()

    _msg("Default Altitude")
    _msg(" - Antenna ID: %d. " % antid)
    _msg(" - Ref = %s. " % ref)
    _msg(" - Position: (%s, %s, %s)." % (X, Y, Z))
    _msg("   Relative Altitude:  %f" % elev)

    return elev


#
# show ATM Profile
# https://casa.nrao.edu/docs/CasaRef/atmosphere-Tool.html#x995-10120004.1.1
#

def showAtmInfo(atm):
    """
     Returned atm (from initAtmProfile) may have a different structure.
     In casa6, additional variable information is assed in Dict type. 
    """
    version = at.getAtmVersion()
    _msg("\nAtomosphere Tool:: version = %s\n" % version)

    # ATM Profile #
    if is_CASA6:
        for s in atm:
            if type(s) is str:
                _msg(s)
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
    # Informatin meesage #
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

    # Antenna Key word for SD #
    if antenna != '':
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
        # msselect=msselect,   ## task mstransform does not support 'msselect' in args.
        reindex=False)   # Must be False #

############################################################
# Calculation Method (Replaced to C++ in next development.)
#    originated by atmcorr_20200602.py by Sawada san.
#    formed as a casa task by CAS-13160.
#    - p_xxxx arguments are for task.
#    - a_xxxx arguments are for Atm-correction.
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
        a_dtem_dh,
        a_h0,
        a_atmtype,
        atmdetail,
        a_altitude,
        a_temperature,
        a_pressure,
        a_humidity,
        a_PWV,
        a_dp,
        a_dpm,
        a_layerboundaries,
        a_layertemperature):

    #
    # Argument dump.
    #   if need to check args,
    #   please insert here like:  print('antenna  =', p_antenna, type(p_antenna))
    #
    if True:  # flag option is reserved. #
        print("***********************************")
        print("**   calc_sdatmcor::             **")
        print("***********************************")
        print('infile      =', p_infile)
        print('datacolumn  =', p_datacolumn)
        print('outfile     =', p_outfile)
        print('overwrite   =', p_overwrite)
        print('field       =', p_field)
        print('spw         =', p_spw)
        print('scan        =', p_scan)
        print('antenna     =', p_antenna)
        print('correlation =', p_correlation)
        print('timerange   =', p_timerange)
        print('intent      =', p_intent)
        print('observation =', p_observation)
        print('feed        =', p_feed)
        print('msselect    =', p_msselect)
        print('outputspw   =', p_outputspw)
        print('dtem_dh     =', a_dtem_dh)
        print('h0          =', a_h0)
        print('atmtype     =', a_atmtype)
        print('atmdetail   =', atmdetail)

        print('altitude    =', a_altitude)
        print('temperature =', a_temperature)
        print('pressure    =', a_pressure)
        print('humidity    =', a_humidity)
        print('PWV         =', a_PWV)
        print('dp          =', a_dp)
        print('dpm         =', a_dpm)
        # reserved #
        print('layerboundaries   =', a_layerboundaries)
        print('layertemperature  =', a_layertemperature)

        print("*****************************")

    # TENTATIVE:: Following flags are deleted soon.  
    skipTaskExec = False          # skip execution at the beginning of calc_sdatmcor.

    # obsoleted debug Vars., soon deleted. #
    showCorrection = False        # show index information while Correction.
    interruptCorrection = False   # Interrupt Correction
    interruptCorrectionCnt = 200  # (limit count)


    # TENTATIVE:: skip task execution, until test-MS for UT is ready.  #
    if(skipTaskExec is True):
        _msg("-------  Task Execution is skipped.")
        return True

    #
    #  File name section
    #

    # datacolumn (XML fills default) ,to UPPER CASE #
    datacolumn = p_datacolumn.upper()
    if (datacolumn == 'CORRECTED'):    # 'CORRECTED' means column:'CORRECTED_DATA'
        datacolumn = 'CORRECTED_DATA'

    # rewrite (Under Construction) - simplified.  #

    rawms = p_infile
    calms = p_outfile
    corms = p_outfile

    # existence info.  #
    rawms_exist = _file_exist(rawms)
    calms_exist = _file_exist(calms)
    corms_exist = _file_exist(corms)

    # File Info. #
    _msg("INPUT and OUTPUT")
    _msg("  default MS file (rawms) = %s , Exist =%s" % (rawms, rawms_exist))
    _msg("  default MS file (calms) = %s , Exist =%s" % (calms, calms_exist))
    _msg("  default MS file (corms) = %s , Exist =%s" % (corms, corms_exist))

    # infile Inaccesible #
    if not rawms_exist:
        _msg("\nERROR::Specified infile does not exist.\n", 'ERROR')
        return False
  
    # outfile Protected #
    if corms_exist:
        if p_overwrite:
            _msg("Overwrite:: Overwrite specified. Once delete the existing output file. ")
            _ms_remove(corms)
        else:
            _msg("\nERROR in Overwrite:: Specified outputfile already exist.\n", 'ERROR')
            return False

    #
    # CAS-13160
    # Original Script starts from here
    #  - in some section, Task-code is inserted.
    #  - comment is also added to distinguish the original and the task code.
    #

    ##################################################
    #   Inside Constant for ATM
    ##################################################

    _msg("Initializing Tool Constant.")

    #
    # Following variables have initial parameters according to the Original script.
    #  - these values will be translated with unit to var.'t_xxxxxxxxx' in Task code.
    #

    # Default constant. refer to CAS-13160 #
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

    #
    # Data Selection
    #

    # Data Selection (internally execute mstransform #
    _msg("Data Selection in progress.")

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
        _msg("Exception:: in atmMst. ")
        casalog.post('%s' % err, 'SEVERE')
        return False

    # Resume 'origin'. A strange behavior in casalog/CASA6 #
    casalog.origin(origin)

    # Result check if output exists. #
    if not _file_exist(corms):
        _msg("\nERROR:: No outfile has been generated by mstransform.\n", 'SEVERE')
        return False

    _msg("Data Selection was applied. Output file =  %s " % corms)

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
        msmd.open(rawms)

        tmonsource = msmd.timesforintent('OBSERVE_TARGET#ON_SOURCE')
        tmoffsource = msmd.timesforintent('OBSERVE_TARGET#OFF_SOURCE')
        fdmspws = msmd.fdmspws()
        tdmspws = msmd.tdmspws()
        intentspws = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
        spws = list(set(intentspws) & (set(fdmspws) | set(tdmspws)))
        spwnames = msmd.namesforspws(spws)

        # show to user. #
        _msg(" - tm ON  source  = %s" % tmonsource)
        _msg(" - tm OFF source  = %s" % tmoffsource)
        _msg(" - fdm    spws = %s" % fdmspws)
        _msg(" - tdm    spws = %s" % tdmspws)
        _msg(" - intent spws = %s" % intentspws)
        _msg(" -        spws = %s" % spws)

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
        # (Task Section )
        # force to change 'processing Spw'
        # - If spw is pecified in the arg, set up 'spws'
        #
        _msg("Determine active spw when spw is in the arg.")
        if (p_spw != ''):
            spws_param = _list_comma_string(p_spw, dType='int')

            # Including check
            #  - all the element in spws_param MUST be in Spws
            if set(spws) >= set(spws_param):   # B is included in A #
                spws = spws_param   # update
                _msg(" - updated 'spws'=%s. by given spw." % spws)
            else:
                _msg("\nERROR:: Some of the specified 'spw' are not in the raw MS. Cannot continue.\n", "ERROR")
                _msg("  - spws (MS)     = %s" % set(spws))
                _msg("  - spws (arg)    = %s" % set(spws_param))
                ## UNDER CONSTRUCTION. This case can continue ##
                return False
 
        # (original) get chanfreqs[spwid] info.
        for spwid in spws:
            chanfreqs[spwid] = msmd.chanfreqs(spw=spwid)

    except Exception as err:
        _msg("ERROR:: opening rawms")
        casalog.post('%s' % err, 'ERROR')
        return False

    finally:
        msmd.close()
        _msg("msmd was successfully closed.")  # TENTATIVE msg (checking try-finally)

    # (original)
    bnd = (pl.diff(tmoffsource) > 1)
    w1 = pl.append([True], bnd)
    w2 = pl.append(bnd, [True])
    tmoffsource = (tmoffsource[w1]+tmoffsource[w2])/2.  ### midpoint of OFF subscan

    ddis = {}

    msmd.open(calms)
    for spwid in spws:
        ddis[spwid] = msmd.datadescids(spw=spwid)[0]
    msmd.close()

    print(" - ddis[] = %s" % ddis)

    nchanperbb = [0, 0, 0, 0]
    bbprs = {}

    for i, spwid in enumerate(spws):
        bbp = int(spwnames[i].split('#')[2][3])-1
        bbprs[spwid] = bbp
        nchanperbb[bbp] += len(chanfreqs[spwid])

    #
    # (Task Section)
    # Output Spw from Arguments. 'Spw' is already decided.
    #
    _msg("Determine outputspw")

    # set Output Spw (if no arg, use spws) #
    if (p_outputspw == ''):
        outputspws = spws         # use calculated 'spws' above
    else:
        outputspws = _list_comma_string(p_outputspw, dType='int')

    _msg("Determined Spw Information")
    _msg('- spws %s' % spws)
    _msg('- outputspws %s' % outputspws)    # expected outputspw

    #
    # (Task Section)
    # Check if specified outputspw is in the list of spws
    #

    if set(spws) >= set(outputspws):   # B is included in A #
        pass
    else:
        _msg("\nERROR:: Some of specified 'outputspw' are not in the raw MS. Cannot continue.\n", "ERROR" )
        _msg("  - spws      = %s" % spws)
        _msg("  - outputspws = %s" % outputspws)

        return False

    #
    # (Original)
    # Data Query on Pointing Table.
    #   output::
    #   - tmpointing, elev
    #
    try:
        tb.open(os.path.join(calms, 'POINTING'))  # CAS-13160:: use tempMS, (in org. using rawms)

        # (org.) key for ANTENNA_ID select
        querytext = 'ANTENNA_ID==%s' % antenna
        subtb = tb.query(querytext)

        # (org.) Access Table
        tmpointing = subtb.getcol('TIME')
        elev = subtb.getcol('DIRECTION').squeeze()[1]
        _msg("- reading column:'TIME' and 'DIRECITON' completed.")

        subtb.close()
        tb.close()

    except Exception as instance:
        _msg("ERROR:: opening POINTING.")
        casalog.post('%s' % instance, 'ERROR')
        raise

    ################################################################
    # Get atmospheric parameters for ATM
    ################################################################
    try:
        tb.open(os.path.join(rawms, 'ASDM_CALWVR'))
        # confirm
        _msg("tmonsource: %f, %f" % (tmonsource.min(), tmonsource.max()))
        pwv = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                       (tmonsource.min(), tmonsource.max())).getcol('water')
        tb.close()

    except Exception as err:
        _msg("ERROR:: opening rawms/'ASDM_CALWVR'.")
        casalog.post('%s' % err, 'ERROR')
        raise

    # pick up pwv
    pwv = pl.median(pwv)

    # ASDM_CALATMOSPHERE
    try:
        _msg("reading rawms/'ASDM_CALATMOSPHERE'.")
        tb.open(os.path.join(rawms, 'ASDM_CALATMOSPHERE'))
        subtb = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                         (tmonsource.min(),
                          tmonsource.max()))

        tground = pl.median(subtb.getcol('groundTemperature'))
        pground = pl.median(subtb.getcol('groundPressure'))
        hground = pl.median(subtb.getcol('groundRelHumidity'))

        subtb.close()
        tb.close()

    except Exception as err:
        _msg("ERROR:: opening rawms/'ASDM_CALATMOSPHERE'.")
        casalog.post('%s' % err, 'ERROR')
        raise

    _msg("median PWV = %fm, T = %fK, P = %fPa, H = %f%%" % (pwv, tground, pground, hground))

    ################################################################
    # Looping over spws
    ################################################################
    tb.open(corms, nomodify=False)

    # Note CAS-13160:
    # spw for-loop:
    # - Inside this loop, calculation is executed with each spwid
    # - The result is wriiten upon OutputFile with selected Spw, which is specified by OutputSpw
    # - List of OutputSpw must be the part of List of Spws
    # - please see the inserted block.
    #

    # for spwid in spws:
    for spwid in outputspws:
        # Log #
        _msg("Processing spw %d in %s " % (spwid, spws))

        istdm = False
        nchan = len(chanfreqs[spwid])
#       fcenter = (chanfreqs[spwid][nchan/2-1]+chanfreqs[spwid][nchan/2])/2.             # PY2
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
#           fcenter = (chanfreqs_high[nchan/2-1]+chanfreqs_high[nchan/2])/2.              # PY2
            fcenter = (chanfreqs_high[int(nchan/2)-1]+chanfreqs_high[int(nchan/2)])/2.    # PY3

        ################################################################
        # Set parameters for ATM and obtain zenith opacity
        # caveat: median values are used for calibrating the entire execution
        ################################################################

        _msg("- set parameters for initATM and obtain zenith opacity")

        #
        # Note CAS-13160: Revised part from original.
        #
        #   Setting up parameter variables.
        #   - Following codes are mainly extracted from the original script, which used to be
        #     coded inside the initAtmProfile().
        #   - Overwrite the parameter to arguments of initAtmProfile.
        #   - This section is up to Task Specification rather than original script.
        #
        t_dtem_dh     = qa.quantity(lapserate, 'K/km')
        t_h0          = qa.quantity(scaleht, 'km')
        t_atmtype     = atmtype
        t_atmdetail   = atmdetail  # Bool
        t_altitude    = qa.quantity(altitude, 'm')
        t_temperature = qa.quantity(tground, 'K')              # tground (original)
        t_pressure    = qa.quantity(pground/100.0, 'mbar')     # pground (original) in  [Pa]  convert to [mbar]
        t_humidity    = qa.quantity(hground, '%' )             # hground (original) in  [%]
        t_pwv         = qa.quantity(pwv * 1000.0, 'mm')        # pwv (original) in [m] converto to [mm]
        t_dp          = qa.quantity(dp, 'mbar')
        t_dpm         = qa.quantity(dpm, '')                   # keep float
        t_maxAltitude = qa.quantity(maxalt, 'km')

        # Edit Flag (True: the parameter is given) #
        t_dtem_dh_set  = False
        t_h0_set       = False
        t_altitude_set = False
        t_temperature_set  = False
        t_pressure_set     = False
        t_humidity_set     = False
        t_pwv_set      = False
        t_dp_set       = False
        t_dpm_set      = False

        #
        # (from 'help' infomation)
        # User-Defined Profile (example)
        #    myalt = [ 5071.72200397, 6792.36546384, 15727.0776121, 42464.18192672 ] #meter
        #    mytemp = [ 270., 264., 258., 252. ] #Kelvin
        #
        t_layerboundaries  = []   # initially null-List.   array: layer boundary [m]
        t_layertemperature = []   # initiallu null-list.   array: layer temperature [K]

        #
        # ATM fundamental parameters activation.
        #
        t_dtem_dh, t_dtem_dh_set = _set_floatquantity_ifactive(a_dtem_dh, t_dtem_dh, 'K/km')
        t_h0,      t_h0_set      = _set_floatquantity_ifactive(a_h0, t_h0, 'km')
        t_atmtype   = _set_int_ifactive(a_atmtype, t_atmtype)

        #
        # Sub parameter activation
        #  only when atmdetail is True
        #
        if t_atmdetail:
            _msg("\nSub Parameter from the Arguments will be set, if specified.\n")
            t_altitude,    t_altitude_set     = _set_floatquantity_ifactive(a_altitude, t_altitude, 'm')
            t_temperature, t_temperature_set  = _set_floatquantity_ifactive(a_temperature, t_temperature, 'K')
            t_pressure,    t_pressure_set     = _set_floatquantity_ifactive(a_pressure, t_pressure, 'mbar')
            t_humidity,    t_humidity_set     = _set_floatquantity_ifactive(a_humidity, t_humidity, '%')
            t_pwv,         t_pwv_set          = _set_floatquantity_ifactive(a_PWV, t_pwv, 'mm')
            t_dp,          t_dp_set           = _set_floatquantity_ifactive(a_dp, t_dp, 'mbar')
            t_dpm,         t_dpm_set          = _set_floatquantity_ifactive(a_dpm, t_dpm, '')

            # User-Uefined Profile. 
            #  when the arg is active, directly pass the arg(a_xxxx) to t_xxxx for initAtmProfile. 
            t_layerboundaries  = _copy_list_ifactive(a_layerboundaries)
            t_layertemperature = _copy_list_ifactive(a_layertemperature)

        else:
            _msg("\nSub Parameters were not used, due to 'atmdetail' is not True.\n")

        #
        # print and log to confirm.
        #
        _msg("====================================================================")
        _msg("  initATMProfile Parameters     [atmdetail = %s]" % t_atmdetail)
        _msg("------------------+-----------------------------+-------------------")
        _msg("  parameter       | value [unit]                | specified by arg. ")
        _msg("------------------+-----------------------------+-------------------")
        _msg(" atmtype          |%-18s " % at.listAtmosphereTypes()[t_atmtype-1] )    # type =1,2,3,4,5, no Zero
        _msg(" dTem_dh          |%-18s [%5s]   | %s " % (t_dtem_dh['value'], t_dtem_dh['unit'], t_dtem_dh_set))
        _msg(" h0               |%-18s [%5s]   | %s " % (t_h0['value'], t_h0['unit'], t_h0_set))
        _msg(" altitude         |%-18s [%5s]   | %s " % (t_altitude['value'], t_altitude['unit'], t_altitude_set))
        _msg(" temperature      |%-18s [%5s]   | %s " % (t_temperature['value'], t_temperature['unit'], t_temperature_set))
        _msg(" pressure         |%-18s [%5s]   | %s " % (t_pressure['value'], t_pressure['unit'], t_pressure_set))
        _msg(" humidity         |%-18s [%5s]   | %s " % (t_humidity['value'], t_humidity['unit'], t_humidity_set))
        _msg(" pwv              |%-18s [%5s]   | %s " % (t_pwv['value'], t_pwv['unit'], t_pwv_set))
        _msg(" dp               |%-18s [%5s]   | %s " % (t_dp['value'], t_dp['unit'], t_dp_set))
        _msg(" dpm              |%-18s [%5s]   | %s " % (t_dpm['value'], t_dpm['unit'], t_dpm_set))
        _msg("*maxAltitude      |%-18s [%5s]   | (FIXED CONST) " % (t_maxAltitude['value'], t_maxAltitude['unit']))
        _msg(" layerboundaries  |%-18s " % t_layerboundaries)
        _msg(" layertemperature |%-18s " % t_layertemperature)
        _msg("------------------+-----------------------------------------------")

        ###################
        # initATMProfile
        ###################
        atm = at.initAtmProfile(humidity=t_humidity['value'],
                                temperature=t_temperature,
                                altitude=t_altitude,
                                pressure=t_pressure,
                                atmType=t_atmtype,
                                maxAltitude=t_maxAltitude,
                                h0=t_h0,
                                dTem_dh=t_dtem_dh,
                                dP=t_dp,
                                dPm=t_dpm['value'],
                                layerBoundaries=t_layerboundaries,
                                layerTemperature=t_layertemperature)

        #
        # show intAtmProfile result
        #
        showAtmInfo(atm)

        #
        #  Layer Information
        #   - automatically showing, when the arg is specified)
        #
        if len(t_layerboundaries) != 0:
            showLayerInfo(at)
     
        ###################
        # Spectral Window #
        ###################
        at.initSpectralWindow(nband,
                              qa.quantity(fcenter, 'Hz'),
                              qa.quantity(nchan*chansep, 'Hz'),
                              qa.quantity(chansep, 'Hz'))

        # H2O
        at.setUserWH2O(t_pwv)

        ################################################################
        # Calculate and apply correction values
        ################################################################

        _msg("\nCalculate and apply correction values.\n")

        # original: make ddis[spwid]
        _msg("- Selecting DATA_DESC_ID == %s" % ddis[spwid])
        querytext = 'DATA_DESC_ID in %s' % ddis[spwid]
        subtb = tb.query(querytext)
  
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

        _msg("\nExecuting  ATM Correction(N=%d), and writing to output MS. \n" % len(tmdata))

        cdata = data.copy()
        for i, t in enumerate(tmdata):
            #
            # debug option, interrupt the correction loop.
            #
            if interruptCorrection:
                if i > interruptCorrectionCnt:
                    msg = "Correction Loop was interrupted. \n---"
                    _msg(msg)
                    break

            # (org script) #
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
                _msg("spw=%3d, i=%5d, time=%15s, Max(dTa)=%19s, Min(dTa)=%19s" % (spwid, i, t, max(dTa), min(dTa)))

            # Adjust Body (dTa is vector) #
            for ipol in range(npol):
                cdata[ipol, :, i] -= dTa

        subtb.putcol(datacolumn, cdata)
        subtb.close()

    #
    # end for spwid in spws:
    #
    tb.flush()
    tb.close()

    return True

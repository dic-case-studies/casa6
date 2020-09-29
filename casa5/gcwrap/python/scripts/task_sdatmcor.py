import os
import pylab as pl
import shutil
import datetime

# for CASA
from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
    from casatools import quanta, table, msmetadata
    from casatools import atmosphere

    # mstransform #
    from casatasks import mstransform

    msmd = msmetadata()
    tb = table()
    qa = quanta()
    at = atmosphere()

    origin = 'sdatmcor'

else:
    from taskinit import tbtool, casalog, qa
    # CAS-13088
    from taskinit import msmdtool as msmetadata
    from casac import casac

    # mstransform #
    from tasks import mstransform

    msmd = msmetadata()
    tb = tbtool()
    at = casac.atmosphere()

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
        layerboundaries, layertemperature,
        debug):

# Information
    casalog.origin(origin)
    msg = "Revision  sdatmcor 0928 for Release #12 and Code-Review."
    _msg(msg)

#
# Input/Output File Handling
#
    # infile oufile, must be specified.
    if infile == '':
        _msg("ERROR:: infile MUST BE  specified.", 'ERROR')
        return False
    if outfile == '':
        _msg("ERROR:: outfile MUST BE specified.", 'ERROR')
        return False

    # infile
    infile_without_ext = os.path.splitext(os.path.basename(infile))[0]
    infile_ext = os.path.splitext(os.path.basename(infile))[1]
    infile = [infile_without_ext, infile_ext]
    # outfile
    outfile_without_ext = os.path.splitext(os.path.basename(outfile))[0]
    outfile_ext = os.path.splitext(os.path.basename(outfile))[1]
    outfile = [outfile_without_ext, outfile_ext]

    # in case infile == outfile
    if infile == outfile:
        _msg("ERROR:: You are attempting to write on input file.", 'ERROR')
        return False

#
# Unit Conversion
#
    dtem_dh     = _form_value_unit(dtem_dh, ['K/km'])
    h0          = _form_value_unit(h0, ['km'])

    altitude    = _form_value_unit(altitude, ['m'])
    temperature = _form_value_unit(temperature, ['K'])
    pressure    = _form_value_unit(pressure, ['mbar', 'hPa'])
    humidity    = humidity  # through (string or float)
    PWV         = _form_value_unit(PWV, ['mm'])
    dp          = _form_value_unit(dp,  ['mbar', 'hPa'])
    dpm         = dpm  # through (string or float)

# User-Defined Profile (nothing =[] )

    layerboundaries = _conv_to_doubleArrayList(layerboundaries)
    layertemperature = _conv_to_doubleArrayList(layertemperature)

    if(len(layerboundaries) != len(layertemperature)):
        _msg("Data count mismatches in specified user-defined parameter.",'ERROR')
        return False

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
        layertemperature,
        debug)

##################
# Data Selection
##################
def atmMst(
    infile, datacolumn, outfile, overwrite,
    field, spw, scan, antenna,
    correlation, timerange, intent,
    observation, feed, msselect):

    # Antenna Key word #
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
        # msselect=msselect,      (Not supported, due to CAS-13160 discussion.)
        reindex=False)   # Must be False #

##########################
# Subroutines
##########################

#
# file handling
#

def _ms_remove(path):
    if (os.path.exists(path)):
        _msg("- Attempt to delete [%s]." % path)
        if (os.path.isdir(path)):
            shutil.rmtree(path)
        else:
            os.remove(path)
    else:
        _msg("- No file to delete [%s]" % path)


def _ms_copy(src, dst):
    _msg("- Copying [%s] ->[%s]." % (src, dst))
    shutil.copytree(src, dst)

def file_exist(path):
    if (os.path.exists(path)):
        return True
    else:
        return False


#
# Unit handling service
#

def _form_value_unit(data, base_unit):
    if (data == ''):
        return ''

    ext_unit = qa.getunit(data)
    if (ext_unit in base_unit):
        # With Unit #
        _msg("Unit Conversion::Data with Unit '%s'" % data)
        return qa.getvalue(data)[0]
    elif (ext_unit == ''):
        # Unit Added #
        _msg("Unit Conversion::No unit specified in %s . Assumed '%s'" % (data, base_unit), 'WARN')
        return data
    else:
        # Mismatch #
        _msg("Unit conversion:: Unexpected Unit '%s' in %s . Abort." % (ext_unit, data), 'SEVERE')
        raise Exception

#
# Argument parameter handling
# (action) check in_para and set default argument.
#  if in_para is available , return in_para with being converted.
#  otherwise, retrurns def_para to use as a default parameter.

def _set_float_param(in_arg, def_para):
    if (in_arg != -1)and(in_arg > 0):
        return in_arg
    else:
        return  def_para


def _set_int_param(in_arg, def_para):
    if (in_arg != ''):
        return int(in_arg)
    else:
        return  def_para


def _set_list_param(in_arg, list_para):
    if (len(in_arg) == 0):
        return []
    else:
        return list_para


def _set_floatquantity_param(in_arg, def_para, unit):
    if (in_arg != ''):
        return  qa.quantity(float(in_arg), unit)
    else:
        return  def_para


def _list_comma_string(separated_string, dType):
    if type(separated_string) is str:
        tmp_list = separated_string.split(',')  # convert to List #
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


def _conv_to_doubleArrayList(in_list):
    # convert elements in  list or separated str =>  float list.
    if  type(in_list) is list:
        out_list = [float(s) for s in in_list]  # force to convert to list[floatm, ...]
        return out_list

    elif type(in_list) is str:
        tmp_list = in_list.split(',')  # convert to List #
        out_list = [float(s) for s in tmp_list]  # force to convert to list[float, ...]
        return out_list
    else:
        _msg("Invalid arg type, expecting separated string or list.", 'SEVERE')
        return []   # invalid type


#
# Logging CASA LOG (INFO/WARN/SEVERE)
#
def _msg(msg, msgtype='INFO'):
    print(msg)
    casalog.post(msg, msgtype, origin=origin)


#########################
# Search ANTENNA_ID
#########################
def get_antennaId(msname, antennaName):
    tb.open(os.path.join(msname, 'ANTENNA'))
    # Query
    queryText = "NAME=='%s'" % antennaName
    subtb = tb.query(queryText)

    # Access Table
    rows = subtb.rownumbers()
    iAnt = rows[0]   # Single Hit is assumed

    # close
    subtb.close()
    tb.close()
    _msg("get_antennaId:: iAnt = %d" % iAnt)
    return iAnt
    

#
# show ATM Profile
# https://casa.nrao.edu/docs/CasaRef/atmosphere-Tool.html#x995-10120004.1.1
#

def showAtmInfo(atm):

    version = at.getAtmVersion()
    _msg("\nAtomosphere Tool:: version = %s\n" % version)
    # ATM Profile #
    if is_CASA6:  # CASA6
        for s in atm:
            if type(s) is str:
                _msg(s)
    else:  # CASA5
        _msg(atm)

def showLayerInfo(at):
    p = at.getProfile()

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
# Calculation Method (Replaced to C++ in next develpment.)
#
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
        a_outputspw,
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
        a_layertemperature,
        debug):

    # Argument dump.
    #   if need to check args,
    #   pleas insert here:  print('antenna     =', p_antenna, type(p_antenna))

    # debug flags  #
    skipTaskExec = False          # skip execution at the begining of calc_sdatmcor.
    showCorrection = False        # show index information while Correction.
    interruptCorrection = False   # Interrupt Correction
    interruptCorrectionCnt = 1000  #  (limit count)

    keepMstTemp = False       # keep tempFile(mstransform out) for debug

    if('skipTaskExec' in debug):
        skipTaskExec = True

    if('interruptCorrection' in debug):
        interruptCorrection = True

    if('showCorrection' in debug):
        showCorrection = True

    if('keepMstTemp' in debug):
        keepMstTemp = True


    # TENTATIVE::  skip task execution for pseudo unit-test.
    if(skipTaskExec is True):
        msg = "-------  Task Execution is skipped."
        _msg(msg)
        return True

    #-----------------------------------------------
    #   Inside Constant for ATM (from org. script)
    #-----------------------------------------------
    atmtype = 2         ### atmType parameter for at (1: tropical, 2: mid lat summer, 3: mid lat winter, etc)
    maxalt = 120        ### maxAltitude parameter for at (km)
    lapserate = -5.6    ### dTem_dh parameter for at (lapse rate; K/km)
    scaleht = 2.0       ### h0 parameter for at (water scale height; km)

    dosmooth = False    ### convolve dTa* spectra with [0.25, 0.5, 0.25] to mimic Hanning spectral response;
                        ### set to True if spectral averaging was not employed for the spw
    dp = 10.0   ### initATMProfile DEFAULT ###
    dpm = 1.2   ### initATMProfile DEFAULT ###

    # other constant
    antenna = 1      # default antenna_id
    altitude = 5073  # default altitude
    nband = 1         # number of band in initSpectralWindow()

    # normalized spectral response for Hanning window, FWHM=10
    hanning = [-0.00098041, -0.00202866, -0.00265951, -0.00222265,
               0.00000000,  0.00465696,  0.01217214,  0.02260546,  0.03556241,
               0.05017949,  0.06519771,  0.07911911,  0.09042184,  0.09779662,
               0.10035898,  0.09779662,  0.09042184,  0.07911911,  0.06519771,
               0.05017949,  0.03556241,  0.02260546,  0.01217214,  0.00465696,
               0.00000000, -0.00222265, -0.00265951, -0.00202866, -0.00098041]

    #
    #  File name section
    #     - 'atmtype' can be used for output file, if necessary.
    #

    # atmtype #
    atmtype_for_file     = _set_int_param(a_atmtype, atmtype)  # 'atmtype_for_file'  is NotACTIVE, reserved. #

    # datacolumn (XML fills default) ,to UPPER CASE #
    datacolumn = p_datacolumn.upper()
    if (datacolumn == 'CORRECTED'):    # 'CORRECTED' means column:'CORRECTED_DATA'
        datacolumn = 'CORRECTED_DATA'

    # infile extension #
    ebuid = p_infile[0]  # filename(with path) without extension.
    ebext = p_infile[1]  # given Extension.

    # outfile, extension #
    outfile = p_outfile[0]
    outext  = p_outfile[1]

    # default file format (original style) #
    rawms = '%s%s' % (ebuid, ebext)
    calms = '%s%s' % (ebuid, ebext)                    ### name of MS in which (standard-)calibrated spectra are stored
    corms = '%s.ms.atm%d' % (ebuid, atmtype_for_file)  ### name of MS form (based on Original)

    # outfile set (=corms), if atmtype is need in outfile, write here as follows, including atmtype.
    #     corms = '%s%s.atm%d' %(ebuid, ebext, atmtype_for_file)     # use same mane as infile #
    #     corms = '%s%s.atm%d' %(outfile, ebext, atmtype_for_file)   # use same mane as infile +'<infileext>' #
    #     corms = '%s%s.atm%d' %(outfile, outext, atmtype_for_file)  # use specified outfile name.ext + atm.n #
    
    corms = '%s%s' % (outfile, outext)     # by specified args. with no other attributes. #

    # existence check
    rawms_exist =  file_exist(rawms)
    calms_exist =  file_exist(calms)
    corms_exist =  file_exist(corms)

    # check #
    _msg("INPUT and OUTPUT")
    _msg("  default MS file (rawms) = %s , Exist =%s" % (rawms, rawms_exist))
    _msg("  default MS file (calms) = %s , Exist =%s" % (calms, calms_exist))
    _msg("  default MS file (corms) = %s , Exist =%s" % (corms, corms_exist))
    _msg("  The 'calms' will soon be swithced to mstransform outfile to use selected data.")

    # infile Inaccesible
    if not rawms_exist:
        _msg("ERROR::Specified infile does not exist..", 'ERROR')
        return False
  
    # outfile Protection
    if corms_exist:
        if p_overwrite:
            _msg("Overwrite:: Overwrite specified. Once delete the existing output file. ")
            _ms_remove(corms)
        else:
            _msg("Overwrite:: Specified outputfile already exist. Abort.",'ERROR')
            return False

    #########################
    # Data Selection
    #  (under construction)
    #########################

    # Temp File (cleaned use)
    tempName = './_AtmCor-Temp'
    tempFileName = tempName + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')+'.ms'
    _msg("- using temp file [%s] for mstransform" % tempFileName)

    # internal Mstransform
    atmMst(
        infile=rawms,           # Tentatively use 'rawms' name
        datacolumn=p_datacolumn,
        outfile=tempFileName,   # Temp name is used inside.
        overwrite=p_overwrite,  # passed by the actual ARG.
        field=p_field,
        spw=p_spw,
        scan=p_scan,
        antenna=p_antenna,
        correlation=p_correlation,
        timerange=p_timerange,
        intent=p_intent,
        observation=p_observation,
        feed=p_feed,
        msselect=p_msselect)

    # set 'origin' again. due to strange behavior in casalog/CASA6 #
    casalog.origin(origin)

    # Result check #
    if not file_exist(tempFileName):
        msg = "No outfile from  mstransform."
        _msg(msg, 'SEVERE')
        return False

    # Antenna ID
    #   only Single antenna is available.
    #   multiple antenna  TBD

    # antenna List, when multiple antennas are given
    ant_list = _list_comma_string(p_antenna, dType='str')
    _msg("- antenna_list = %s" % ant_list)

    # antenna (this is used inside the script)
    antenna = get_antennaId(rawms, p_antenna)

    # enable Data Selection
    # copy temp-out to calms for input MS.  #
    msg = "Data Selection was applied. Use %s as 'calms' " % tempFileName
    _msg(msg)

    calms = tempFileName
    _ms_copy(src=tempFileName, dst=corms)

    ################################################################
    # Get metadata
    ################################################################
    chanfreqs = {}
    _msg("\nSDATMCOR main body starts. rawms=%s\n" % rawms)

    try:
        _msg("msmd.open(rawms)", rawms)
        msmd.open(rawms)
        tmonsource = msmd.timesforintent('OBSERVE_TARGET#ON_SOURCE')
        tmoffsource = msmd.timesforintent('OBSERVE_TARGET#OFF_SOURCE')
        fdmspws = msmd.fdmspws()
        tdmspws = msmd.tdmspws()
        intentspws = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
        spws = list(set(intentspws) & (set(fdmspws) | set(tdmspws)))
        spwnames = msmd.namesforspws(spws)
 
    except Exception as err:
        _msg("ERROR:: opening rawms")
        casalog.post('%s' % err, 'ERROR')

    # check count of ON/OFF SOURCE
    n_tmonsource = len(tmonsource)
    n_tmoffsource = len(tmoffsource)
    msg = "Target Information. \n"   \
        + "# ON_SOURCE: count of tmonsource   = %d\n" % n_tmonsource  \
        + "# OFF_SOURCE: count of tmoffsource = %d" % n_tmoffsource
    _msg(msg)

    # OFF_SOURCE check
    if (n_tmoffsource == 0):
        msg = "Can't find the OFF_SOURCE data."
        _msg(msg, 'SEVERE')
        return False

    # ON_SOURCE check
    if (n_tmoffsource == 0):
        msg = "Can't find the ON_SOURCE data."
        _msg(msg, 'SEVERE')
        return False

    for spwid in spws:
        chanfreqs[spwid] = msmd.chanfreqs(spw=spwid)
    msmd.close()
    _msg("- closed msmd")

    bnd = (pl.diff(tmoffsource) > 1)
    w1 = pl.append([True], bnd)
    w2 = pl.append(bnd, [True])
    tmoffsource = (tmoffsource[w1]+tmoffsource[w2])/2.  ### midpoint of OFF subscan

    ddis = {}

    _msg("msmd.open(calms).")
    msmd.open(calms)
    for spwid in spws:
        ddis[spwid] = msmd.datadescids(spw=spwid)[0]
    msmd.close()
    _msg("closed msmd")

    nchanperbb = [0, 0, 0, 0]
    bbprs = {}

    for i, spwid in enumerate(spws):
        bbp = int(spwnames[i].split('#')[2][3])-1
        bbprs[spwid] = bbp
        nchanperbb[bbp] += len(chanfreqs[spwid])

    # Data Query (on Pointing Table)
    try:
        _msg("opening calms/POINTING.")
        tb.open(os.path.join(calms, 'POINTING'))  # CAS-13160:: use tempMS, (in org. using rawms)

        querytext = 'ANTENNA_ID==%s' % antenna
        subtb = tb.query(querytext)

        # Access Table
        _msg("- reading column:'TIME'.")
        tmpointing = subtb.getcol('TIME')
        _msg("- reading column:'DIRECITON'.")
        elev = subtb.getcol('DIRECTION').squeeze()[1]
        subtb.close()
        tb.close()
        _msg("- closed POINTING.")

    except Exception as instance:
        _msg("ERROR:: opening POINTING.")
        casalog.post('%s' % instance, 'ERROR')
        raise

    ################################################################
    # Get atmospheric parameters for ATM
    ################################################################
    try:
        _msg("opening rawms/'ASDM_CALWVR'.")
        tb.open(os.path.join(rawms, 'ASDM_CALWVR'))
        # confirm #
        _msg("tmonsource: %f, %f" % (tmonsource.min(), tmonsource.max()))
        pwv = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                       (tmonsource.min(), tmonsource.max())).getcol('water')
        tb.close()
        _msg("- closed rawms/'ASDM_CALWVR'.")

    except Exception as err:
        _msg("ERROR:: opening rawms/'ASDM_CALWVR'.")
        casalog.post('%s' % err, 'ERROR')
        raise

    # pick up pwv #
    pwv = pl.median(pwv)

    # ASDM_CALATMOSPHERE
    try:
        _msg("opening rawms/'ASDM_CALATMOSPHERE'.")
        tb.open(os.path.join(rawms, 'ASDM_CALATMOSPHERE'))
        subtb = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                         (tmonsource.min(),
                          tmonsource.max()))

        tground = pl.median(subtb.getcol('groundTemperature'))
        pground = pl.median(subtb.getcol('groundPressure'))
        hground = pl.median(subtb.getcol('groundRelHumidity'))

        subtb.close()
        tb.close()
        _msg("- closed rawms/'ASDM_CALATMOSPHERE'.")

    except Exception as err:
        _msg("ERROR:: opening rawms/'ASDM_CALATMOSPHERE'.")
        casalog.post('%s' % err, 'ERROR')
        raise

    _msg("median PWV = %fm, T = %fK, P = %fPa, H = %f%%" % (pwv, tground, pground, hground))

    #
    # prepare SPW, output SPW
    #
    _msg("Determine spws, outputspw")

    # set processing SPW #
    if (p_spw != ''):
        spws = _list_comma_string(p_spw, dType='int')

    # set Output SPW (if no arg, use spws) #
    if (a_outputspw == ''):
        outputspws = spws                    # use calculated 'spws' above
    else:
        outputspws = _list_comma_string(a_outputspw, dType='int')

    _msg('- spws %s' % spws)
    _msg('- outputspws %s' % outputspws)    # expected outputspw

    # Check if specified if outputspw is in the list of spws
    _found = False
    for spw in outputspws:
        if spw in spws:
            _found = True
            break
        else:
            _found = False
 
    if not _found:
        msg = "None of the output-spws %s is not in the processing spws %s. Abort." % (outputspws, spws)
        _msg(msg, 'ERROR')
        return False

    ################################################################
    # Looping over spws
    ################################################################
    _msg("opening corms[%s] to write ATM-Corrected Data." % corms)
    tb.open(corms, nomodify=False)
    for spwid in spws:
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
        # Default Value and Argument parameter set
        # initial tool parameters are ready. use these is needed.
        #
        t_dtem_dh     = qa.quantity(lapserate, 'K/km')
        t_h0          = qa.quantity(scaleht, 'km')
        t_atmtype     = atmtype
        t_atmdetail   = atmdetail  # Bool, directly fron the arg.
        t_altitude    = qa.quantity(altitude, 'm')
        t_temperature = qa.quantity(tground, 'K')              # tground (original)
        t_pressure    = qa.quantity(pground/100.0, 'mbar')     # pground (original) in  [Pa]  convert to [mbar]
        t_humidity    = hground  # keep float                  # hground  in [percent =%]
        t_pwv         = qa.quantity(pwv * 1000.0, 'mm')        # pwv (original) im [m] converto to [mm]
        t_dp          = qa.quantity(dp, 'mbar')
        t_dpm         = dpm      # keep float
        t_maxAltitude = qa.quantity(maxalt, 'km')

        # user-defined Profile (example)
        #    myalt = [ 5071.72200397, 6792.36546384, 15727.0776121, 42464.18192672 ] #meter
        #    mytemp = [ 270., 264., 258., 252. ] #Kelvin
        t_layerboundaries  = a_layerboundaries   # array: layer boundary [m]
        t_layertemperature = a_layerboundaries   # array: layer temerature [K]

        #
        # Overwrite Values, if the Argument is specified.
        #

        # ATM fundamental #
        t_dtem_dh     = _set_floatquantity_param(a_dtem_dh, t_dtem_dh, 'K/km')
        t_h0          = _set_floatquantity_param(a_h0, t_h0, 'km')
        t_atmtype     = _set_int_param(a_atmtype, t_atmtype)

        # Sub parameter, only when atmdetai is True #
        if t_atmdetail:
            _msg("\nSub Parameter from the Arguments will be set, if specified.\n")
            t_altitude    = _set_floatquantity_param(a_altitude, t_altitude, 'm')
            t_temperature = _set_floatquantity_param(a_temperature, t_temperature, 'K')
            t_pressure    = _set_floatquantity_param(a_pressure, t_pressure, 'mbar')
            t_humidity    = _set_float_param(a_humidity, t_humidity)
            t_pwv         = _set_floatquantity_param(a_PWV, t_pwv, 'mm')
            t_dp          = _set_floatquantity_param(a_dp, t_dp, 'mbar')
            t_dpm         = _set_float_param(a_dpm, t_dpm)
            # user-defined profile.
            t_layerboundaries  = _set_list_param(a_layerboundaries, a_layerboundaries)
            t_layertemperature = _set_list_param(a_layertemperature, a_layertemperature)

        else:
            _msg("\nSub Parameters were not used, due to 'atmdetail' is not True.\n")

        # print and log to confirm #
        _msg("==========================================================")
        _msg("  initATMProfile Parameters to set up. [atmdetail = %s]   " % t_atmdetail)
        _msg("-----------------+----------------------------------------")
        _msg(" dTem_dh         |%s" % t_dtem_dh)
        _msg(" h0              |%s" % t_h0)
        _msg(" altitude        |%s" % t_altitude)
        _msg(" temperature     |%s" % t_temperature)
        _msg(" pressure        |%s" % t_pressure)
        _msg(" humidity        |%s" % t_humidity)
        _msg(" pwv             |%s" % t_pwv)
        _msg(" dp              |%s" % t_dp)
        _msg(" dpm             |%s" % t_dpm)
        _msg("*maxAltitude     |%s" % t_maxAltitude)
        _msg(" layerboundaries |%s" % t_layerboundaries)
        _msg(" layerboundaries |%s" % t_layertemperature)
        _msg("-----------------+----------------------------------------")

        # initATMProfile #
        atm = at.initAtmProfile(humidity=t_humidity,
                                temperature=t_temperature,
                                altitude=t_altitude,
                                pressure=t_pressure,
                                atmType=t_atmtype,
                                maxAltitude=t_maxAltitude,
                                h0=t_h0,
                                dTem_dh=t_dtem_dh,
                                dP=t_dp, dPm=t_dpm,
                                layerBoundaries=t_layerboundaries,
                                layerTemperature=t_layertemperature)

        # show intAtm result Information #
        showAtmInfo(atm)

        #  Layer Info (automatically showing, when the arg is specified) #
        if len(t_layerboundaries) != 0:
            showLayerInfo(at)

        # Spectral Window #
        at.initSpectralWindow(nband,
                              qa.quantity(fcenter, 'Hz'),
                              qa.quantity(nchan*chansep, 'Hz'),
                              qa.quantity(chansep, 'Hz'))

        # H2O  #
        at.setUserWH2O(t_pwv)

        ################################################################
        # Calculate and apply correction values
        ################################################################

        _msg("- Calculate and apply correction values")
        _msg("- Selecting DATA_DESC_ID == %s" % ddis[spwid])

        # Skip outputspw
        if not spwid in outputspws:
        # if not ddis[spwid] in outputspws:
            msg = "This spw %d is skipped, due to not in the output spw list." % spwid
            _msg(msg)
            continue

        # Essential Query #
        querytext = 'DATA_DESC_ID in %s' % ddis[spwid]
        subtb = tb.query(querytext)

        # time data and numPol #
        _msg("- getting tm and data. datacolumn [%s] is used." % datacolumn)
        tmdata = subtb.getcol('TIME')
        data = subtb.getcol(datacolumn)
        npol = data.shape[0]
      
        # Smoothing control #
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

        msg = "\nWriting ATM Correction Data. (N=%d) \n" % len(tmdata)
        _msg(msg)

        cdata = data.copy()
        for i, t in enumerate(tmdata):

            # debug option, interrupt the correction loop.
            if interruptCorrection:
                if i > interruptCorrectionCnt:
                    msg = "Correction Loop was interrupted. \n---"
                    _msg(msg)
                    break

            # (org script)
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

            # debug option for task. #
            if showCorrection:
                _msg("spw=%3d, i=%5d, time=%15s, Max(dTa)=%19s, Min(dTa)=%19s" % (spwid, i, t, max(dTa), min(dTa)))

            # Adjust Body (dTa is vector) #
            for ipol in range(npol):
                cdata[ipol, :, i] -= dTa

        subtb.putcol(datacolumn, cdata)
        subtb.close()

    # end for spwid in spws:
    tb.flush()
    tb.close()

    # LOG #
    _msg("- closed MS[%s] to write." % corms)

    # delete temp file. (with Tentative debug option)
    if not keepMstTemp:
        _ms_remove(tempFileName)

    # finish #
    return True

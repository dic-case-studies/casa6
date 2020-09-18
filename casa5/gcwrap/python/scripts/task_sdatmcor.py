# common import
# from Script
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
    msg = "Revision  sdatmcor 0917-1 testing data-selection"
    print(msg)
    casalog.post(msg, 'INFO', origin=origin)

#
# File Handling
#
    # infile oufile, must be specified.
    if infile == '':
        print ("FATAL:: infile MUST BE  specified.")
        return False
    if outfile == '':
        print ("FATAL:: outfile MUST BE specified.")
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
        msg = "You are atempting to write on input file."
        print(msg)
        casalog.post(msg, 'SEVERE', origin=origin)

        raise Exception

#
# Unit Conversion
#
    dtem_dh     = form_value_unit(dtem_dh, ['K/km'])
    h0          = form_value_unit(h0, ['km'])

    altitude    = form_value_unit(altitude, ['m'])
    temperature = form_value_unit(temperature, ['K'])
    pressure    = form_value_unit(pressure, ['mbar', 'hPa'])
    humidity    = humidity  # through (string or float)
    PWV         = form_value_unit(PWV, ['mm'])
    dp          = form_value_unit(dp,  ['mbar', 'hPa'])
    dpm         = dpm  # through (string or float)

# User-Define Profile (nothing =[] )

    # if(len(layerboundaries) != len(layertemperature)):
    #     print("WARN: specified Count of Bounday and Temperature mismatch.")
    #  param = ['1.1', '2.2', '3.3' .....]  
    layerboundaries = conv_to_doubleArrayList(layerboundaries)
    layertemperature = conv_to_doubleArrayList(layertemperature)

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
        vis=infile,               ## Full file spec.
        outputvis=outfile,        ## overwrite is not allowed. 
        datacolumn=datacolumn,
        field=field, 
        spw=spw, 
        scan=scan,
        antenna=antenna,        ## ex)'PM01&&&', 
        correlation=correlation, 
        timerange=timerange, 
        intent=intent,
        observation=observation, 
        feed=feed, 
#        msselect=msselect,      ## THIS MAKES ERROR when used. ##
        reindex=False)   # CAUTION #

##########################
# Subroutines
##########################

# by  Wataru, thanks. #
def ms_remove(path):
    if (os.path.exists(path)):
        print("- Atempt to delete [%s]."%path)
        if (os.path.isdir(path)):
            shutil.rmtree(path)
        else:
            os.remove(path)
    else:
        print("- No file to delete [%s]"%path)


def ms_copy(src, dst):
    print("- Copying [%s] ->[%s]."%(src, dst))
    shutil.copytree(src, dst)

def file_exist(path):
    if (os.path.exists(path)):
        return True
    else:
        return False


def form_value_unit(data, base_unit):
    if (data == ''):
        return ''  

    ext_unit = qa.getunit(data)
    if (ext_unit in base_unit):
        # INFO #
        msg = "INFO::Data with Unit '%s'" % data
        print(msg)
        casalog.post(msg, 'INFO', origin=origin)
        return qa.getvalue(data)[0]
    elif (ext_unit == ''):
        # INFO #
        msg = "WARN::No unit specified in %s . Assumed '%s'" % (data, base_unit)
        casalog.post(msg, 'WARN', origin=origin)
        return data
    else:
        # FATAL #
        msg = "FATAL:: Unexpected Unit '%s' in %s . Aborted." % (ext_unit,data)
        casalog.post(msg, 'SEVERE', origin=origin)

        raise Exception


def set_float_param(in_arg, def_para):
    #   print( "set_float_param::", in_arg, def_para )
    if (in_arg != -1)and(in_arg > 0):
        return in_arg
    else:
        return  def_para


def set_int_param(in_arg, def_para):
    if (in_arg != ''):
        return int(in_arg)
    else:
        return  def_para


def set_any_param(in_arg, def_para):
    if (in_arg != ''):
        return in_arg
    else:
        return def_para

def set_list_param(in_arg, list_para):
    if (len(in_arg) == 0):
        return []
    else:
        return list_para


def set_floatquantity_param(in_arg, def_para, unit):
    if (in_arg != ''):
        return  qa.quantity(float(in_arg), unit)
    else:
        return  def_para


def set_antenna_param(in_arg, def_para):
    if (in_arg != ''):
        return int(in_arg)
    else:
        return  def_para


def list_comma_string( separated_string, dType):
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


def listInt_comma_string( separated_string):
    if type(separated_string) is str:
        tmp_list = separated_string.split(',')  # convert to List #
        out_list = [int(s) for s in tmp_list]  # convert to list[int]
        return out_list
    elif type(separated_string) is list:
        return separated_string
    else:
        return []

def listStr_comma_string( separated_string):
    """
      make a list from comma separated string
    """
    if type(separated_string) is str:
        tmp_list = separated_string.split(',')  # convert to List #
        out_list = [str(s) for s in tmp_list]  # convert to list[int]
        return out_list
    elif type(separated_string) is list:
        return separated_string
    else:
        return []

def conv_to_doubleArrayList( in_list ):
    """
      convert elements in a list, to double espression
    """
    if  not (type(in_list) is list):
        return []
    out_list = [float(s) for s in in_list]  # convert to list[int]
    return out_list

#########################
# Search ANTENNA_ID 
#########################
def get_antennaId( msname, antennaName):
    tb.open(os.path.join(msname, 'ANTENNA'))
    # Query
    queryText = "NAME=='%s'" % antennaName
    subtb = tb.query(queryText)

    # Access Table
    rows = subtb.rownumbers()
    iAnt=rows[0]   # Single Hit is assumed 

    # close
    subtb.close()
    tb.close()
    print('get_antennaId:: iAnt = ', iAnt)
    return iAnt
    

#
# ATM Profile
#

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

    print(p[0])
    return
#
# Sum Method (Replaced to C++ in future)
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

    if True:  # flag option is reserved. #
        print("*********************************************")
        print("**   calc_sdatmcor::                       **")
        print("*********************************************")
        print('infile      =', p_infile, type(p_infile))
        print('datacolumn  =', p_datacolumn, type(p_datacolumn))
        print('outfile     =', p_outfile, type(p_outfile))
        print('overwrite   =', p_overwrite, type(p_overwrite))
        print('field       =', p_field, type(p_field))
        print('spw         =', p_spw,type(p_spw))
        print('scan        =', p_scan, type(p_scan))
        print('antenna     =', p_antenna, type(p_antenna))
        print('correlation =', p_correlation, type(p_correlation))
        print('timerange   =', p_timerange, type(p_timerange))
        print('intent      =', p_intent, type(p_intent))
        print('observation =', p_observation, type(p_observation))
        print('feed        =', p_feed, type(p_feed))
        print('msselect    =', p_msselect,type(p_msselect))
        print('outputspw   =', a_outputspw, type(a_outputspw))
        print('dtem_dh     =', a_dtem_dh, type(a_dtem_dh))
        print('h0          =', a_h0, type(a_h0))
        print('atmtype     =', a_atmtype, type(a_atmtype))
        print('atmdetail   =', atmdetail, type(atmdetail))

        print('altitude    =', a_altitude, type(a_altitude))
        print('temperature =', a_temperature, type(a_temperature))
        print('pressure    =', a_pressure, type(a_pressure))
        print('humidity    =', a_humidity, type(a_humidity))
        print('PWV         =', a_PWV, type(a_PWV))
        print('dp          =', a_dp, type(a_dp))
        print('dpm         =', a_dpm, type(a_dpm))
        # reserved #
        print('layerboundaries   =', a_layerboundaries, type(a_layerboundaries))
        print('layertemperature  =', a_layertemperature,type(a_layertemperature))

        print('debug       =', debug)
        print("*****************************")

    # debug flags  #
    skipTaskExec = False          # skip execution in do_sdatmcor() until official test code is ready
    skipCorrection = False        # skip Correction procedure in the script, to save time.
    showCorrection = False        # show index information while Correction.
    interruptCorrection = False   # Interrupt Correction
    interruptCorrectionCnt =1000  #  (limit count)

    keepMstTemp = False       # keep tempFile(mstransform out) for debug 

    if('skipTaskExec' in debug):
        skipTaskExec = True

    if('skipCorrection' in debug):
        skipCorrection = True

    if('interruptCorrection' in debug):
        interruptCorrection = True

    if('showCorrection' in debug):
        showCorrection = True

    if('keepMstTemp' in debug):
        keepMstTemp = True


    #--------------------------------------
    #   Inside Constant for ATM (from org. script)
    #-------------------------------------
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
    nband =1         # number of band in initSpectralWindow()

    # normalized spectral response for Hanning window, FWHM=10
    hanning = [-0.00098041, -0.00202866, -0.00265951, -0.00222265,
        0.00000000,  0.00465696,  0.01217214,  0.02260546,  0.03556241,
        0.05017949,  0.06519771,  0.07911911,  0.09042184,  0.09779662,
        0.10035898,  0.09779662,  0.09042184,  0.07911911,  0.06519771,
        0.05017949,  0.03556241,  0.02260546,  0.01217214,  0.00465696,
        0.00000000, -0.00222265, -0.00265951, -0.00202866, -0.00098041]

    # TENTATIVE::  skip task execution while UnitTest is incomplete.
    if(skipTaskExec is True):
        msg = "-------  Task Execution is skipped."
        print(msg)
        casalog.post(msg)
        return True

    #
    #  File name section
    #     - 'atmtype' is needed in advance
    #

    # atmtype #
    atmtype_for_file     = set_int_param(a_atmtype,atmtype)

    # datacolumn (XML fills default) ,to UPPER CASE #
    datacolumn = p_datacolumn.upper()
    if (datacolumn == 'CORRECTED'):  # add '_DATA' if CORRECED #
        datacolumn = 'CORRECTED_DATA' 

    # infile, supprty Extension #
    ebuid = p_infile[0]  ### The Arg must have a filename WITHOUT extension###
    ebext = p_infile[1]  ### given Extension  ###

    # outfile, support Extension #
    outfile = p_outfile[0]  ### The Arg must have a filename WITHOUT extension###
    outext  = p_outfile[1]  ### given Extension  ###

    # default file format (original style) #
    rawms = '%s%s' % (ebuid, ebext)
    calms = '%s%s' % (ebuid, ebext)                    ### name of MS in which (standard-)calibrated spectra are stored
    corms = '%s.ms.atm%d' % (ebuid, atmtype_for_file)  ### name of MS form (based on Original)

    # outfile (=corms), Assist (reserved if needed)
    if False:
        if (outfile == ''):
            corms = '%s%s.atm%d' %(ebuid, ebext, atmtype_for_file)     # use same mane as infile #
        elif (outext == ''):
            corms = '%s%s.atm%d' %(outfile, ebext, atmtype_for_file)   # use same mane as infile +'<infileext>' #
        else:
            corms = '%s%s.atm%d' %(outfile, outext, atmtype_for_file)  # use specified outfile name.ext + atm.n #
    else:
        corms = '%s%s' %(outfile, outext)     # use same mane as infile #


    # existence check
    rawms_exist =  file_exist(rawms)
    calms_exist =  file_exist(calms)
    corms_exist =  file_exist(corms)

    # check #
    print("INPUT and OUTPUT")
    print("default MS file (rawms) = %s , Exist =%s" % (rawms, rawms_exist))
    print("default MS file (calms) = %s , Exist =%s" % (calms, calms_exist))
    print("default MS file (corms) = %s , Exist =%s" % (corms, corms_exist)) 

    # infile inaccesible  
    if not rawms_exist:
        print("FATAL:: Specified infile does not exist..")
        return False
  
    # Overwrite Protection 
    if corms_exist:
        if p_overwrite:
            print("INFO:: Overwriting to output, once delete the existing file. " )
            ms_remove(corms)
        else:
            print("FATAL:: Specified outputfile already exist. Abort.")
            return False

#
# TENTATIVE: pre-process (Data Selection)
#    (3-Sep-2020 ~ underconstruction)
#

    #########################
    # Data Selection
    #  (under construction)
    #########################

    # Temp File (cleaned use)
    tempName = './_AtmCor-Temp'
    tempFileName = tempName + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')+'.ms'
    print("- use temp file [%s]" % tempFileName)

    # internal Mstransform
    atmMst(
        infile=rawms,           # Tentatively use 'rawms' name
        datacolumn=p_datacolumn, 
        outfile= tempFileName,      # Temp name is used inside. Ignoring the arg. # 
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

    #
    # Antenna ID
    #   only Single antenna is available    
    #   multiple antenna  TBD
    #

    # antenna List, when muslipe antennas are given
    ant_list = list_comma_string(p_antenna, dType='str')
    print("- antenna_list =", ant_list )
    # antenna (this is used inside the script) 
    antenna = get_antennaId( rawms, p_antenna)

    #
    # TESTING :: Switch cal-MS
    #
    if False:  # (this is original)
        # Prepare MS. Ignore the tempFile.
        ms_remove(corms)
        ms_copy(src=calms, dst=corms)
    else: 
        # Here is the DEFAUT route. #
        msg = "Data Selection is applied." 
        print(msg)
        casalog.post(msg, 'INFO', origin=origin)

        # Prepare MS, hand to AtmCor
        calms = tempFileName
        ms_copy(src=tempFileName, dst=corms)

    ################################################################
    # Get metadata
    ################################################################
    chanfreqs = {}
    print("SDATMCOR main body.", rawms)

    print("msmd.open(rawms)", rawms)
    msmd.open(rawms)
    tmonsource = msmd.timesforintent('OBSERVE_TARGET#ON_SOURCE')
    tmoffsource = msmd.timesforintent('OBSERVE_TARGET#OFF_SOURCE')
    fdmspws = msmd.fdmspws()
    tdmspws = msmd.tdmspws()
    intentspws = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
    spws = list(set(intentspws) & (set(fdmspws) | set(tdmspws)))
    spwnames = msmd.namesforspws(spws)

    print("#ON_SOURCE: count of tmonsource   = %d" % len(tmonsource))
    print("#OFF_SOURCE: count of tmoffsource = %d" % len(tmoffsource))

    for spwid in spws:
        chanfreqs[spwid] = msmd.chanfreqs(spw=spwid)
    msmd.close()
    print("closed msmd")

    bnd = (pl.diff(tmoffsource) > 1)
    w1 = pl.append([True], bnd)
    w2 = pl.append(bnd, [True])
    tmoffsource = (tmoffsource[w1]+tmoffsource[w2])/2.  ### midpoint of OFF subscan

    ddis = {}

    print("msmd.open(calms).")
    msmd.open(calms)
    for spwid in spws:
        ddis[spwid] = msmd.datadescids(spw=spwid)[0]
    msmd.close()
    print("closed msmd")

    nchanperbb = [0, 0, 0, 0]
    bbprs = {}

    for i, spwid in enumerate(spws):
        bbp = int(spwnames[i].split('#')[2][3])-1
        bbprs[spwid] = bbp
        nchanperbb[bbp] += len(chanfreqs[spwid])

    # Data Query (on Pointing Table)
    if False:  # (original)
        print("Opening rawms/POINTING.")
        tb.open(os.path.join(rawms, 'POINTING'))  # use original MS.
    else:
        print("Opening calms/POINTING.") 
        tb.open(os.path.join(calms, 'POINTING'))  # use tempMS

    querytext = 'ANTENNA_ID==%s' % antenna
    subtb = tb.query(querytext)

    # Access Table
    tmpointing = subtb.getcol('TIME')
    elev = subtb.getcol('DIRECTION').squeeze()[1]
    subtb.close()
    tb.close()
    print("closeed POINTING.")

    ################################################################
    # Get atmospheric parameters for ATM
    ################################################################
    print("opening rawms/'ASDM_CALWVR'.")
    tb.open(os.path.join(rawms, 'ASDM_CALWVR'))
    # confirm #
    print("tmonsource:",tmonsource.min(),tmonsource.max() )
    pwv = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                   (tmonsource.min(), tmonsource.max())).getcol('water')
    tb.close()
    print("closed rawms/'ASDM_CALWVR'.")

    # pick up pwv #
    pwv = pl.median(pwv)

    # ASDM_CALATMOSPHERE
    print("opening rawms/'ASDM_CALATMOSPHERE'.")
    tb.open(os.path.join(rawms, 'ASDM_CALATMOSPHERE'))
    subtb = tb.query('%.3f<=startValidTime && startValidTime<=%.3f' %
                     (tmonsource.min(),
                      tmonsource.max()))

    tground = pl.median(subtb.getcol('groundTemperature'))
    pground = pl.median(subtb.getcol('groundPressure'))
    hground = pl.median(subtb.getcol('groundRelHumidity'))

    subtb.close()
    tb.close()
    print("closed rawms/'ASDM_CALATMOSPHERE'.")

    print('median PWV = %fm, T = %fK, P = %fPa, H = %f%%' % (pwv, tground, pground, hground))

    # prepare SPW

    print("Determine spws, outputspw")
    # set processing SPW, (not output_SPW)  #
    if (p_spw != ''):
        spws = list_comma_string(p_spw, dType='int')

    # set Output SPW (if ommited, use spws)
    if (a_outputspw == ''):
        outputspws = spws                    # use calculated 'spws' above
    else:
        outputspws = list_comma_string(a_outputspw, dType='int')

    print('- spws %s' % spws)
    print('- outputspws %s' %outputspws)    # This is the expected outputspw

    # Check if specified  outputspw is in spws
    _found = False
    for spw in outputspws:
        if spw in spws:
            _found = True
            break
        else:
            _found = False 
 
    if not _found:
        # LOG #
        msg = "None of the output-spws (%s) is not in the processing spws[%s]. Abort." % (outputspws,spws)
        casalog.post(msg,'SEVERE',origin=origin)
        return False 

        """
        # if continue, helped by outputspw = spws #
        outputspws = spws    
        print('- spws = %s, outputspws = %s ', (spws,outputspws))
        """
    ################################################################
    # Looping over spws
    ################################################################
    print("opening corms[%s] to write ATM-Corrected Data."% corms)
    tb.open(corms, nomodify=False)
    for spwid in spws:
        # Log #
        msg = "Processing spw %d in %s " % (spwid,spws)
        print(msg)
        casalog.post(msg,'INFO',origin=origin)

        istdm = False
        nchan = len(chanfreqs[spwid])
#       fcenter = (chanfreqs[spwid][nchan/2-1]+chanfreqs[spwid][nchan/2])/2.             # PY2
        fcenter = (chanfreqs[spwid][int(nchan/2)-1]+chanfreqs[spwid][int(nchan/2)])/2.   # PY3
        chansep = (chanfreqs[spwid][-1]-chanfreqs[spwid][0])/(nchan-1)

        print('- checking nchanperbb[bbprs[spwid]]')
        if nchanperbb[bbprs[spwid]] in [128, 256]:

            print('Spw %d seems to be TDM-like. More accurate Hanning smoothing is applied. istdm=True' % spwid)

            # ATM model is calculated for finer spectral resolution (5x original),
            # convolved with Hanning spectral response with FWHM=10,
            # then resampled to the original channel freqs

            istdm = True
            nchan *= 5
            chansep /= 5.
            chanfreqs_high = chanfreqs[spwid][0] + chansep * pl.arange(nchan)
#           fcenter = (chanfreqs_high[nchan/2-1]+chanfreqs_high[nchan/2])/2.              # PY"
            fcenter = (chanfreqs_high[int(nchan/2)-1]+chanfreqs_high[int(nchan/2)])/2.    # PY3

        ################################################################
        # Set parameters for ATM and obtain zenith opacity
        # caveat: median values are used for calibrating the entire execution
        ################################################################

        print("- set parameters for initATM and obtain zenith opacity")

        #--------------------------------------------
        # Default Value and Argument parameter set
        # initial fundamental parameter value
        #-------------------------------------------
        t_dtem_dh     = qa.quantity(lapserate, 'K/km')
        t_h0          = qa.quantity(scaleht, 'km')
        t_atmtype     = atmtype
        t_atmdetail   = atmdetail  # Bool, directly fron the arg.
        t_altitude    = qa.quantity(altitude, 'm')       
        t_temperature = qa.quantity(tground, 'K')              # tground (original)
        t_pressure    = qa.quantity(pground/100.0, 'mbar')     # pground (original) in  [Pa]  convert to [mbar]
        t_humidity    = hground  # float                       # hground  in [percent =%]
        t_pwv         = qa.quantity(pwv * 1000.0, 'mm')        # pwv (original) im [m] converto to [mm]
        t_dp          = qa.quantity(dp, 'mbar')  
        t_dpm         = dpm      # float      
        t_maxAltitude = qa.quantity(maxalt, 'km')

        t_layerboundaries  = a_layerboundaries   # array: layer boundary [m]
        t_layertemperature = a_layerboundaries   # array: layer temerature [K]

        # user-defined Profile (example)
        #    myalt = [ 5071.72200397, 6792.36546384, 15727.0776121, 42464.18192672 ] #meter
        #    mytemp = [ 270., 264., 258., 252. ] #Kelvin

        #------------------------------------------
        # Change value, if the Argument  specified 
        #------------------------------------------

        # ATM fundamental
        t_dtem_dh     = set_floatquantity_param(a_dtem_dh, t_dtem_dh, 'K/km')
        t_h0          = set_floatquantity_param(a_h0, t_h0, 'km')
        t_atmtype     = set_int_param(a_atmtype, t_atmtype)

        # sub parameter, when atmdetai is True #
        if t_atmdetail:
            print ("-- Sub Parameter from the Arguments will be set, if specified." )
            t_altitude    = set_floatquantity_param(a_altitude, t_altitude, 'm')
            t_temperature = set_floatquantity_param(a_temperature, t_temperature, 'K')
            t_pressure    = set_floatquantity_param(a_pressure, t_pressure, 'mbar')
            t_humidity    = set_float_param(a_humidity, t_humidity)
            t_pwv         = set_floatquantity_param(a_PWV, t_pwv, 'mm')
            t_dp          = set_floatquantity_param(a_dp, t_dp, 'mbar')
            t_dpm         = set_float_param(a_dpm, t_dpm)
            # user-defined profile.
            t_layerboundaries  = set_list_param(a_layerboundaries, a_layerboundaries) 
            t_layertemperature = set_list_param(a_layertemperature, a_layertemperature)

        else:
            print ("-- Sub Parameters were not used, due to 'atmdetail' is not True." )

        # print to confirm #
        print("==========================================================")
        print("  initATMProfile Parameters to set up. [atmdetail = %s]   " % t_atmdetail)
        print("-----------------+----------------------------------------")
        print(" dTem_dh         |%s" % t_dtem_dh)
        print(" h0              |%s" % t_h0)
        print(" altitude        |%s" % t_altitude)
        print(" temperature     |%s" % t_temperature)
        print(" pressure        |%s" % t_pressure)
        print(" humidity        |%s" % t_humidity)
        print(" pwv             |%s" % t_pwv)
        print(" dp              |%s" % t_dp)
        print(" dpm             |%s" % t_dpm)
        print("*maxAltitude     |%s" % t_maxAltitude)
        print(" layerboundaries |%s" % t_layerboundaries)
        print(" layerboundaries |%s" % t_layertemperature)
        print("-----------------+----------------------------------------")
        # initATMProfile #
        myAtm = at.initAtmProfile(
            humidity=t_humidity, 
            temperature=t_temperature, 
            altitude=t_altitude, 
            pressure=t_pressure, 
            atmType=t_atmtype, 
            maxAltitude=t_maxAltitude, 
            h0=t_h0, 
            dTem_dh=t_dtem_dh,
            dP=t_dp, dPm=t_dpm,
            layerBoundaries=t_layerboundaries,
            layerTemperature=t_layertemperature )

        # ATM Profile #
        if is_CASA6: 
            for s in myAtm:
                print(s)
        else:
            print(myAtm)

        #  LAYER INFO (when the arg specified) #
        if len(t_layerboundaries) !=0:
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

        print("- Calculate and apply correction values")
        print("- Selecting DATA_DESC_ID == %s"%ddis[spwid])

        # Skip outputspw
        if not spwid in outputspws:  
#       if not ddis[spwid] in outputspws:
            msg = "This spw %d is skipped, due to this is not in the output spw." % spwid
            print(msg)
            casalog.post(msg,'INFO', origin=origin)
            continue

        # Essential Query (required by org. script) #
        querytext = 'DATA_DESC_ID in %s' % ddis[spwid] 
        subtb = tb.query(querytext)

        # time data and numPol #
        print( "- datacolumn [%s] is used." % datacolumn)
        tmdata = subtb.getcol('TIME')
        data = subtb.getcol(datacolumn)
        npol = data.shape[0]
      
        # Smoothing control #
        if nchanperbb[bbprs[spwid]]*npol in [256, 8192]:
            print('Spw %d in BB_%d (total Nchan within BB is %d, sp avg likely not applied).  dosmooth=True' %
                  (spwid, bbprs[spwid]+1, nchanperbb[bbprs[spwid]]*npol))
            dosmooth = True
        else:
            print('Spw %d in BB_%d (total Nchan within BB is %d, sp avg likely applied).  dosmooth=False' %
                  (spwid, bbprs[spwid]+1, nchanperbb[bbprs[spwid]]*npol))
            dosmooth = False

        # Debug(Tentatitve) Skip the main correction loop. 
        if skipCorrection:
            msg = "Correction loop will be skipped, due to Debug option. (nRow=%d)\n---" % len(tmdata)
            print(msg)
            casalog.post(msg,'INFO',origin=origin)
            continue

        ###########################
        # Correction Main Loop
        ###########################

        # LOG #
        msg = "ATM Correction for loop (N=%d), " % len(tmdata)
        print(msg)
        casalog.post(msg,'INFO',origin=origin)

        cdata = data.copy()
        for i, t in enumerate(tmdata):

            # debug option, limit the correction loop. 
            if interruptCorrection:
                if i > interruptCorrectionCnt:
                    msg = "Correction Loop was interrupted. \n---"
                    print(msg)
                    casalog.post(msg,'INFO',origin=origin)
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
                print("spw, i, time, Max(dTa),Min(dTa)", spwid, i, t, max(dTa),min(dTa))

            # Adjust Body (dTa is vector) #
            for ipol in range(npol):
                cdata[ipol, :, i] -= dTa

        subtb.putcol(datacolumn, cdata)
        subtb.close()

    # end for spwid in spws:
    tb.flush()
    tb.close()

    # LOG #
    msg = "closed MS[%s] to write."% corms
    print(msg)
    casalog.post(msg,'INFO',origin=origin)

    # delete temp file. TENTATIVE
    if keepMstTemp:
        None
    else:
        ms_remove(tempFileName)

    # finish #
    return True




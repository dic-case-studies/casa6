from __future__ import absolute_import
import shutil

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import quanta, measures, table, ms, imager
    from casatasks import casalog
    from .mstools import write_history
    
    _myqa = quanta( )
    _myme = measures( )
else:
    from taskinit import *
    # not local tools
    _myqa = qa
    _myme = me

def fixvis(vis, outputvis='',field='', refcode='', reuse=True, phasecenter='', distances='', datacolumn='all'):
    """
    Input Parameters
    vis        -- Name of the input visibility set
    
    outputvis  -- Name of the output visibility set, default: same as vis

    field      -- field selection string

    refcode    -- Reference frame to convert to,
                  default: the refcode of PHASE_DIR in the FIELD table
                  example: 'B1950'
    
    reuse      -- base recalculation on existing UVW coordinates? default=True
                  ignored if parameter 'phasecenter' is set

    phasecenter  --  if set to a valid direction: change the phase center for the given field
                     to this value
                     example: 'J2000 9h25m00s 05d12m00s'
                     If given without the equinox, e.g. '0h01m00s 00d12m00s', the parameter
                     is interpreted as a pair of offsets in RA and DEC to the present phasecenter.

    distances -- (experimental) List of the distances (as quanta) of the fields selected by field
                 to be used for refocussing.
                 If empty, the distances of all fields are assumed to be infinity.
                 If not a list but just a single value is given, this is applied to
                 all fields.
                 default: []
                 examples: ['2E6km', '3E6km']   '15au'

    datacolumn -- when applying a phase center shift, modify visibilities only in this/these column(s)
                  default: 'all' (DATA, CORRECTED, and MODEL)
                  example: 'DATA,CORRECTED' (will not modify MODEL)

    Examples:

    fixvis('NGC3256.ms','NGC3256-fixed.ms')
          will recalculate the UVW coordinates for all fields based on the existing
          phase center information in the FIELD table.

    fixvis('0925+0512.ms','0925+0512-fixed.ms','0925+0512', '', 'J2000 9h25m00s 05d12m00s')
          will set the phase center for field '0925+0512' to the given direction and recalculate
          the UVW coordinates.
    """

    # Note: this is duplicated in task_cvel, and really needing CASA-wide harmonization
    # (CAS-12871)
    def copy_ms(src, dest):
        """ This is a MMS-safe copy of an MMS tree directory.
        :param src: path to the source MS
        :param dest: path to the destination MS
        """
        shutil.copytree(src, dest, symlinks=True)

    casalog.origin('fixvis')

    if vis == outputvis or outputvis == '':
        casalog.post('Will overwrite original MS ...', 'NORMAL')
        outputvis = vis
    else:
        casalog.post('Copying original MS to outputvis ...', 'NORMAL')

        shutil.rmtree(outputvis, ignore_errors=True)
        copy_ms(vis, outputvis)

    if is_CASA6:
        tbt = table( )
        myms = ms( )
        myim = imager( )
    else:
        tbt, myms, myim = gentools(['tb', 'ms', 'im'])

    if field == '' or isinstance(field,list) and len(field) == 0:
        field='*'

    fields = myms.msseltoindex(vis=outputvis,field=field)['field']

    if len(fields) == 0:
        casalog.post( "Field selection returned zero results.", 'WARN')
        return

    thedistances = []

    if distances == "":
        distances = []
    elif distances != []:
        if type(distances) == str and _myqa.isquantity(distances):
            thedist = _myqa.canonical(distances)
            if thedist['unit'] == 'm': # a length
                for f in fields: # put nfields copies into the list
                    thedistances.append(thedist['value'])
            else:
                msg = "Parameter distances needs to contain quanta with units of length."
                raise ValueError(msg)

        elif type(distances) == list:
            if len(fields) != len(distances):
                msg = "You selected "+str(len(fields))+" fields but gave "+str(len(distances))+" distances,"
                raise ValueError(msg)
            else:
                for d in distances:
                    if _myqa.isquantity(d):
                        thedist = _myqa.canonical(d)
                        if thedist['unit'] == 'm': # a length
                            thedistances.append(thedist['value'])
                        else:
                            msg = "Parameter distances needs to contain quanta with units of length."
                            raise ValueError(msg)
        else:
            msg = "Invalid parameter distances."
            raise ValueError(msg)


    if thedistances != []:
        casalog.post('Will refocus to the given distances: '+str(distances), 'NORMAL')

    #determine therefcode, the reference frame to be used for the output UVWs
    tbt.open(outputvis+"/FIELD")
    numfields = tbt.nrows()
    therefcode = 'J2000'
    ckwdict = tbt.getcolkeyword('PHASE_DIR', 'MEASINFO')
    tbt.close()
    if refcode == '':
        if 'TabRefTypes' in ckwdict: # we have a variable reference column
            therefcode = 'J2000' # always use "J2000"
        else: # not a variable reference column
            therefcode = ckwdict['Ref']
    else: # a refcode was given, enforce validity
        if not (type(refcode)==str):
            msg = 'Invalid refcode '+str(refcode)
            raise RuntimeError(msg)
        if 'TabRefTypes' in ckwdict: # variable ref column
            refcodelist = ckwdict['TabRefTypes'].tolist()
            ref = 0
            if not (refcode in refcodelist):
                msg = 'Invalid refcode '+refcode
                raise RuntimeError(msg)
        else: # not a variable reference column
            if not (refcode in get_validcodes()):
                msg = 'Invalid refcode '+refcode
                raise RuntimeError(msg)
        #endif
        therefcode = refcode
    #end if

    if phasecenter == '': # we are only modifying the UVW coordinates
        casalog.post('Will leave phase centers unchanged.', 'NORMAL')
        casalog.post("Recalculating the UVW coordinates ...", 'NORMAL')

        fldids = []
        for i in range(numfields):
            if (i in fields):
                fldids.append(i)

        # 
        myim.open(outputvis, usescratch=False)
        myim.calcuvw(fields=fldids, refcode=therefcode, reuse=reuse)
        myim.close()
    else: # we are modifying UVWs and visibilities
        ## if observation:
        ##     casalog.post('Modifying the phase tracking center(s) is imcompatible', 'SEVERE')
        ##     casalog.post('with operating on only a subset of the observation IDs', 'SEVERE')
        ##     return False

        if type(phasecenter) != str:
            msg = "Invalid phase center."
            raise ValueError(msg)

        theoldref, theoldrefstr, ckwdict, isvarref, flddict = get_oldref(outputvis, tbt)

        # for the case of a non-variable reference column and several selected fields 
        commonoldrefstr = ''

        # handle the datacolumn parameter
        if (not type(datacolumn)==str):
            casalog.post("Invalid parameter datacolumn", 'SEVERE')
        elif datacolumn=='' or datacolumn.lower()=='all':
            datacolumn='all'
            casalog.post("Will modify the visibilities in all columns", 'NORMAL')
        else:
            # need to check datacolumn before any part of the MS gets modified
            wantedcolumns = datacolumn.split(',')
            tbt.open(outputvis)
            thecolumns = tbt.colnames()
            tbt.close()
            for col in wantedcolumns:
                if not (col.lower() in ['data','observed','corrected', 'corrected_data','model','model_data']):
                    msg = "Invalid datacolumn: \""+col+"\""
                    raise RuntimeError(msg)
                if (col.lower()=='observed'):
                    col = 'DATA'
                if (col.lower()=='corrected'):
                    col = 'CORRECTED_DATA'
                if (col.lower()=='model'):
                    col = 'MODEL_DATA'
                if not col.upper() in thecolumns:
                    msg = "Datacolumn "+col+" not present"
                    raise RuntimeError(msg)

            casalog.post("Will only modify the visibilities in the columns "+datacolumn.upper(), 'NORMAL')

        for fld in fields:
            allselected = True
            for i in range(numfields):
                if not (i in fields):
                    allselected = False
                    break

            commonoldrefstr = modify_fld_vis(fld, outputvis, tbt, myim,
                                             commonoldrefstr, phasecenter,
                                             therefcode, reuse, numfields,
                                             ckwdict, theoldref, theoldrefstr,
                                             isvarref, flddict, datacolumn,
                                             allselected, thedistances)
            if commonoldrefstr == False:
                raise RuntimeError('Failure in modify_fld_vis)')
    #endif change phasecenter

    # Write history to output MS
    try:
        param_names = fixvis.__code__.co_varnames[:fixvis.__code__.co_argcount]
        if is_python3:
            vars = locals( )
            param_vals = [vars[p] for p in param_names]
        else:           
            param_vals = [eval(p) for p in param_names]   
        write_history(myms, outputvis, 'fixvis', param_names, param_vals,
                      casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')        

def get_validcodes(code=None):
    """Returns a list of valid refcodes."""
    if not get_validcodes._hc:    # Is it not already cached?
        # Making it once per session is often enough.
        validcodes = _myme.listcodes(_myme.direction('J2000', '0','0'))
        get_validcodes._normal = validcodes['normal'].tolist()
        get_validcodes._extra  = validcodes['extra'].tolist()
    if code == 'extra':
        return get_validcodes._extra
    elif code == 'normal':
        return get_validcodes._normal
    else:
        return get_validcodes._normal + get_validcodes._extra

get_validcodes._hc = False

def is_valid_refcode(refcode):
    """Is refcode usable?"""
    return refcode in ['J2000', 'B1950', 'B1950_VLA', 'HADEC', 'ICRS'] \
            + get_validcodes('extra')

def log_phasecenter(oldnew, refstr, ra, dec):
    """Post a phase center to the logger, along with whether it is old or new."""
    casalog.post(oldnew + ' phasecenter RA, DEC ' + refstr + ' '
                 + _myqa.time(_myqa.quantity(ra, 'rad'), 10)[0] # 10 digits precision
                 + " " + _myqa.angle(_myqa.quantity(dec, 'rad'), 10)[0], 'NORMAL')
    casalog.post('          RA, DEC (rad) ' + refstr + ' '
                 + str(ra) + " " + str(dec), 'NORMAL')

def get_oldref(outputvis, tbt):
    """Returns the original reference code, string, ckwdict, and isvarref."""
    theoldref = -1
    theoldrefstr = ''    
    tbt.open(outputvis + "/FIELD")
    ckwdict = tbt.getcolkeyword('PHASE_DIR', 'MEASINFO')
    flddict = {}
    colstoget = ['PHASE_DIR', 'NAME']
    if 'TabRefTypes' in ckwdict and 'TabRefCodes' in ckwdict:
        colstoget.append('PhaseDir_Ref')
    for c in colstoget:
        flddict[c] = tbt.getcol(c)
    if flddict['PHASE_DIR'].shape[1] > 1:
        casalog.post('The input PHASE_DIR column is of order '
                     + str(flddict['PHASE_DIR'].shape[1] - 1), 'WARN')
        casalog.post('Orders > 0 are poorly tested.', 'WARN')
    flddict['PHASE_DIR'] = flddict['PHASE_DIR'].transpose((2, 0, 1))
    tbt.close()
    if 'TabRefTypes' in ckwdict and 'TabRefCodes' in ckwdict:
        isvarref = True
    else:
        isvarref = False
        theoldrefstr = ckwdict['Ref']
    return theoldref, theoldrefstr, ckwdict, isvarref, flddict

def modify_fld_vis(fld, outputvis, tbt, myim, commonoldrefstr, phasecenter,
                   therefcode, reuse, numfields, ckwdict, theoldref,
                   theoldrefstr, isvarref, flddict, datacol, allselected,
                   thedistances):
    """Modify the UVW and visibilities of field fld."""
    viaoffset = False
    thenewra_rad = 0.
    thenewdec_rad = 0.
    thenewref = -1
    theolddir = flddict['PHASE_DIR'][fld]
    fieldname = flddict['NAME'][fld]
    if 'TabRefTypes' in ckwdict and 'TabRefCodes' in ckwdict:
        # determine string name of the phase dir reference frame
        theoldref = flddict['PhaseDir_Ref'][fld]
        refcodestrlist = ckwdict['TabRefTypes'].tolist()
        refcodelist = ckwdict['TabRefCodes'].tolist()
        if not (theoldref in refcodelist):
            msg = 'Invalid refcode in FIELD column PhaseDir_Ref: ' + str(theoldref)
            raise RuntimeError(msg)
        theoldrefstr = refcodestrlist[refcodelist.index(theoldref)]

    if not isvarref:
        if not (commonoldrefstr == ''):
            theoldrefstr = commonoldrefstr
        else:
            commonoldrefstr = theoldrefstr

    
    theoldphasecenter = theoldrefstr + ' ' + \
                        _myqa.time(_myqa.quantity(theolddir[0], 'rad'), 14)[0] + ' ' + \
                        _myqa.angle(_myqa.quantity(theolddir[1],'rad'), 14)[0]

    if not is_valid_refcode(theoldrefstr):
        casalog.post('Refcode for FIELD column PHASE_DIR is valid but not yet supported: '
                     + theoldrefstr, 'WARN')
        casalog.post('Output MS may not be valid.', 'WARN')

    casalog.post('field: ' + fieldname, 'NORMAL')
    log_phasecenter('old', theoldrefstr, theolddir[0], theolddir[1])

    if therefcode != 'J2000':
        casalog.post(
                 "When changing phase center, can only write new UVW coordinates in J2000.",
                 'WARN')
        therefcode = 'J2000'
    if reuse:
        casalog.post("When changing phase center, UVW coordinates will be recalculated.",
                     'NORMAL')
        reuse = False

    dirstr = parse_phasecenter(phasecenter, isvarref, theoldref, theoldrefstr, theolddir)
    if not dirstr:
        raise RuntimeError('Failed to parse phasecenter')

    if isvarref:
        thenewrefindex = ckwdict['TabRefTypes'].tolist().index(dirstr[0])
        thenewref = ckwdict['TabRefCodes'][thenewrefindex]
        thenewrefstr = dirstr[0]
    else: # not a variable ref col
        validcodes = get_validcodes()
        if dirstr[0] in validcodes:
            thenewref = validcodes.index(dirstr[0])
            thenewrefstr = dirstr[0]
        else:
            msg = 'Invalid refcode ' + dirstr[0]
            raise RuntimeError(msg)
        if dirstr[0] != ckwdict['Ref']:
            if numfields > 1 and not allselected:
                        msg = ("You have not selected all " + str(numfields) +
                               " fields and PHASE_DIR is not a variable reference column.\n"
                               " Please use split or provide phase dir in " + ckwdict['Ref']
                               + ".")
                        raise RuntimeError(msg)
            else:
                casalog.post(
            "The direction column reference frame in the FIELD table will be changed from "
                             + ckwdict['Ref'] + " to " + dirstr[0], 'NORMAL')
    #endif isvarref

    try:
        thedir = _myme.direction(thenewrefstr, dirstr[1], dirstr[2])
        thenewra_rad = thedir['m0']['value']
        thenewdec_rad = thedir['m1']['value']
    except Exception as instance:
        msg = "*** Error \'%s\' when interpreting parameter \'phasecenter\': " % (instance)
        raise RuntimeError(msg)

    if not is_valid_refcode(thenewrefstr):
        casalog.post('Refcode for the new phase center is valid but not yet supported: '
                     + thenewrefstr, 'WARN')
        casalog.post('Output MS may not be valid.', 'WARN')

    if theolddir[0] >= _myqa.constants('pi')['value']: # old RA uses range 0 ... 2 pi, not -pi ... pi
        while (thenewra_rad < 0.): # keep RA positive in order not to confuse the user
            thenewra_rad += 2. * _myqa.constants('pi')['value']
                
    log_phasecenter('new', thenewrefstr, thenewra_rad, thenewdec_rad)

    # modify FIELD table                
    tbt.open(outputvis + '/FIELD', nomodify=False)
    pcol = tbt.getcol('PHASE_DIR')
    pcol[0][0][fld] = thenewra_rad
    pcol[1][0][fld] = thenewdec_rad
    tbt.putcol('PHASE_DIR', pcol)
        
    casalog.post("FIELD table PHASE_DIR column changed for field " + str(fld) + ".",
                 'NORMAL')

    if thenewref != -1:
        # modify reference of the phase dir column; check also the
        # other direction columns
        theoldref2 = -1
        theoldref3 = -1
        if isvarref:
            pcol = tbt.getcol('PhaseDir_Ref')
            #theoldref was already determined further above
            #theoldref = pcol[fld]
            pcol[fld] = thenewref
            
            pcol2 = tbt.getcol('DelayDir_Ref')
            theoldref2 = pcol2[fld]
            
            pcol3 = tbt.getcol('RefDir_Ref')
            theoldref3 = pcol3[fld]

            if theoldref != thenewref:
                tbt.putcol('PhaseDir_Ref', pcol)
                casalog.post(
                    "FIELD table phase center direction reference frame for field "
                    + str(fld) + " set to " + str(thenewref) + " ("
                    + thenewrefstr + ")", 'NORMAL')
                if not (thenewref == theoldref2 and thenewref == theoldref3):
                    casalog.post(
            "*** The three FIELD table direction reference frame entries for field "
                                 + str(fld)
                                 + " will not be identical in the output data: "
                                 + str(thenewref) + ", " + str(theoldref2) + ", "
                                 + str(theoldref3), 'NORMAL')
                    if not (theoldref == theoldref2 and theoldref == theoldref3):
                        casalog.post(
            "*** The three FIELD table direction reference frame entries for field "
                                     + str(fld)
                                     + " were not identical in the input data either: "
                                     + str(theoldref) + ", " + str(theoldref2)
                                     + ", " + str(theoldref3), 'NORMAL')
            else:
                casalog.post(
            "FIELD table direction reference frame entries for field  " + str(fld)
                             + " unchanged.", 'NORMAL')

        else: # not a variable reference column
            tmprec = tbt.getcolkeyword('PHASE_DIR', 'MEASINFO')
            if theoldrefstr != thenewrefstr:
                tmprec['Ref'] = thenewrefstr
                tbt.putcolkeyword('PHASE_DIR', 'MEASINFO', tmprec) 
                casalog.post(
                    "FIELD table phase center direction reference frame changed from "
                             + theoldrefstr + " to " + thenewrefstr, 'NORMAL')
    tbt.close()
    
    fldids = []
    phdirs = []
    for i in range(numfields):
        if (i==fld):
            fldids.append(i)
            phdirs.append(theoldphasecenter)

    if thedistances==[]:
        thedistances = 0. # the default value

    # 
    myim.open(outputvis, usescratch=False)
    myim.fixvis(fields=fldids, phasedirs=phdirs, refcode=therefcode, datacolumn=datacol, distances=thedistances)
    myim.close()
    return commonoldrefstr

def parse_phasecenter(phasecenter, isvarref, ref, refstr, theolddir):
    dirstr = phasecenter.split(' ')
    if len(dirstr) == 2: # interpret phasecenter as an offset
        casalog.post("No equinox given in parameter \'phasecenter\': "
                     + phasecenter, 'NORMAL')         
        casalog.post("Interpreting it as pair of offsets in (RA,DEC) ...",
                     'NORMAL')

        if isvarref and ref > 31:
            casalog.post('*** Refcode in FIELD column PhaseDir_Ref is a solar system object: '
                         + refstr, 'NORMAL')
            casalog.post(
        '*** Will use the nominal entry in the PHASE_DIR column to calculate new phase center',
                         'NORMAL')
                    
        qra = _myqa.quantity(theolddir[0], 'rad') 
        qdec = _myqa.quantity(theolddir[1], 'rad')
        qraoffset = _myqa.quantity(dirstr[0])
        qdecoffset = _myqa.quantity(dirstr[1])
        if not _myqa.isangle(qdecoffset):
            msg = "Invalid phasecenter parameter. DEC offset must be an angle."
            raise RuntimeError(msg)
        qnewdec = _myqa.add(qdec,qdecoffset)
        qnewra = qra
        ishms = (_myqa.canonical(qraoffset)['unit'] == 'rad') and (('h' in dirstr[0] and 'm' in dirstr[0] and 's' in dirstr[0]) or (dirstr[0].count(':')==2))
        if (_myqa.canonical(qraoffset)['unit'] == 'rad') and not ishms:
            casalog.post(
            "RA offset is an angle (not a time). Will divide by cos(DEC) to compute time offset.")
            if _myqa.cos(qnewdec)['value'] == 0.0:
                casalog.post(
                   "Resulting DEC is at celestial pole. Will ignore RA offset.", 'WARN')
            else:
                qraoffset = _myqa.div(qraoffset, _myqa.cos(qnewdec))
                qnewra = _myqa.add(qnewra, qraoffset)
        else:
            if not ((_myqa.canonical(qraoffset)['unit'] == 's') or ishms):
                msg = "Invalid phasecenter parameter. RA offset must be an angle or a time."
                raise RuntimeError(msg)
            # RA offset was given as a time; apply as is
            qraoffset = _myqa.convert(qraoffset, 'deg')
            qnewra = _myqa.add(qnewra, qraoffset)

        dirstr = [refstr, _myqa.time(qnewra,12)[0], _myqa.angle(qnewdec,12)[0]]

    elif not len(dirstr) == 3:
        msg = 'Incorrectly formatted parameter \'phasecenter\': ' + phasecenter
        raise RuntimeError(msg)
    return dirstr

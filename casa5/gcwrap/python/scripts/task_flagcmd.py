from __future__ import absolute_import
import os
import copy
import numpy as np
from collections import defaultdict

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatasks import casalog
    from casatools import ms, quanta, table, agentflagger
    from .mstools import write_history
    from . import flaghelper as fh

    qalocal = quanta( )
    tblocal = table( )
else:
    from taskinit import casalog, casac, qa, tb
    from mstools import write_history
    import flaghelper as fh

    # naked tool constructors ala CASA6
    agentflagger = casac.agentflagger
    ms = casac.ms

    # not really local
    qalocal = qa
    tblocal = tb

# common function to use to get a dictionary item iterator
if is_python3:
    def lociteritems(adict):
        return adict.items()
else:
    def lociteritems(adict):
        return adict.iteritems()

def flagcmd(
    vis=None,
    inpmode=None,
    inpfile=None,
    tablerows=None,
    reason=None,
    useapplied=None,
    tbuff=None,
    ants=None,
    action=None,
    flagbackup=None,
    clearall=None,
    rowlist=None,
    plotfile=None,
    savepars=None,
    outfile=None,
    overwrite=None
    ):

    #
    # Task flagcmd
    #    Reads flag commands from file or string and applies to MS

    try:
        from xml.dom import minidom
    except Exception as exc:
        raise ImportError('Failed to load xml.dom.minidom into python: {}'.format(exc))

    casalog.origin('flagcmd')

    aflocal = agentflagger()
    mslocal = ms()
    mslocal2 = ms()           

    try:
        # Use a default ntime to open the MS. The user-set ntime will be
        # used in the tool later
        ntime = 0.0

        # Open the MS and attach it to the tool
        if (type(vis) == str) & os.path.exists(vis):
            aflocal.open(vis, ntime)
        else:
            raise ValueError( 'Visibility data set not found - please verify the name' )

        # Check if vis is a cal table:
        # typevis = 1 --> cal table
        # typevis = 0 --> MS
        # typevis = 2 --> MMS
        iscal = False
        typevis = fh.isCalTable(vis)
        if typevis == 1:
            iscal = True
            
            if action != 'apply' and action != 'list':
                raise ValueError( 'Unsupported action for cal tables. Only apply and list are supported.' )
            
            if inpmode == 'table' and isinstance(inpfile, str) and inpfile == '':
                raise ValueError( 'inpmode=\'table\' needs an MS as inpfile' )
                                 
            flagcmds = {}       
            if inpmode == 'table' and fh.isCalTable(inpfile) == 0:
                # Read flag cmds from the MS
                flagcmds = readCalCmds(vis, inpfile, [], tablerows, reason, useapplied)
                listmode = 'cmd'
                
            elif inpmode == 'list':
                # Read flag cmds from a list
                flagcmds = readCalCmds(vis, '', inpfile, [], reason, True)  
                listmode = ''              
            else:
                raise ValueError( 'Unsupported inpmode for cal tables' )
                        
            # Apply flag cmds
            if len(flagcmds.keys( )) == 0:
                raise RuntimeError( 'There are no unapplied flags in the input. '\
                                 'Set useapplied=True to also use the previously-applied flags.' )
            
            # List flags on the screen/logger
            if action == 'list':
                casalog.post('Executing action = list')
                listFlagCmd(myflags=flagcmds, myoutfile='', listmode=listmode)
                
            elif action == 'apply':
                casalog.post('Executing action = apply')
                applyCalCmds(aflocal, vis, flagcmds, tablerows, flagbackup, outfile)
                
            # Save the flag cmds to an output file
            if savepars:
                if not overwrite and os.path.exists(outfile):
                    raise ValueError( 'You have set overwrite to False. Remove %s before saving the flag commands'%outfile )
                    
                fh.writeFlagCommands(vis, flagcmds, False, '', outfile, False)
                                                                
        else:
            # Input vis is an MS

            # Get overall MS time range for later use (if needed)
            mslocal2.open(vis)
            timd = mslocal2.range(['time'])
            mslocal2.close()
            if len(timd) != 0:
                ms_startmjds = timd['time'][0]
                ms_endmjds = timd['time'][1]
                t = qalocal.quantity(ms_startmjds, 's')
                t1sdata = t['value']
                ms_starttime = qalocal.time(t, form='ymd', prec=9)[0][0]
                t = qalocal.quantity(ms_endmjds, 's')
                t2sdata = t['value']
                ms_endtime = qalocal.time(t, form='ymd', prec=9)[0]
                # NOTE: could also use values from OBSERVATION table col TIME_RANGE
                casalog.post('MS spans timerange ' + ms_starttime + ' to '
                             + ms_endtime)
    
            myflagcmd = {}
    
            if action == 'clear':
                casalog.post('Action "clear" will disregard inpmode (no reading)')
                # Clear flag commands from FLAG_CMD in vis
                msfile = vis
    
                if clearall:
                    casalog.post('Deleting all rows from FLAG_CMD in MS '
                                 + msfile)
                    clearFlagCmd(msfile, myrowlist=rowlist)
                else:
                    casalog.post('Safety Mode: you chose not to set clearall=True, no action'
                                 )
                return
            
            elif inpmode == 'table':
    
                casalog.post('Reading from FLAG_CMD table')
                # Read from FLAG_CMD table into command list
                if inpfile == '':
                    msfile = vis
                else:
                    msfile = inpfile
    
                myflagcmd = readFromTable( msfile, myflagrows=tablerows,
                                           useapplied=useapplied, myreason=reason)
    
                listmode = 'table'
            elif inpmode == 'list':
                if action == 'unapply':
                    casalog.post("The unapply action can only be used with inpmode='table'",'WARN')
                    casalog.post("Save the commands to the FLAG_CMD table before using unapply",'WARN')
                    raise ValueError( "Unsupported action='unapply' for inpmode='list'" )
    
                # ##### TO DO: take time ranges calculation into account ??????
                # Parse the input file
                try:            
                    # Input commands are given in a list
                    if isinstance(inpfile, list):
                        
                        # It is a list of input files
                        if os.path.isfile(inpfile[0]):
                            flaglist = []
                            for ifile in inpfile:
                                casalog.post('Will read commands from the file '+ifile)                    
                                flaglist = flaglist + fh.readFile(ifile)
                            
                            myflagcmd = fh.parseDictionary(flaglist, reason)
                        
                        # It is a list of strings with flag commands
                        else:
                            casalog.post('Will read commands from a Python list')
                            myflagcmd = fh.parseDictionary(inpfile, reason)
                        listmode = ''
                        
                    # Input commands are given in a file
                    elif isinstance(inpfile, str):
                        
                        if inpfile == '':
                             casalog.post('Input file is empty', 'ERROR')
                             
                        casalog.post('Will read commands from the file '+inpfile)
                        flaglist = fh.readFile(inpfile)
                        casalog.post('%s'%flaglist,'DEBUG')
                        
                        myflagcmd = fh.parseDictionary(flaglist, reason)
                        listmode = 'file'
                    
                    else:
                        casalog.post('Input type is not supported', 'ERROR')
                        
                    listmode = 'list'
                    casalog.post('%s'%myflagcmd,'DEBUG1')
                                                                                            
                except Exception as instance:
                    raise Exception( 'Error reading the input list: {} '.format(instance))
                
    
            elif inpmode == 'xml':
    
                casalog.post('Reading from Flag.xml')
                if action == 'unapply':
                    casalog.post("The unapply action can only be used with inpmode='table' or 'list';'",'WARN')
                    casalog.post("save the commands to the FLAG_CMD table before using unapply.",'WARN')
                    raise ValueError( "Unsupported action='unapply' for inpmode='xml'" )

                # Read from Flag.xml (also needs Antenna.xml)
                if inpfile == '':
                    flagtable = vis
                else:
                    flagtable = inpfile
    
                # Actually parse table. Fail if Flag.xml or Antenna.xml is not found
                try:
                    myflags = fh.parseXML(flagtable, mytbuff=tbuff)            
                except Exception as exc:
                    raise RuntimeError('Error while parsing XML: {}'.format(exc))
    
                casalog.post('%s' % myflags, 'DEBUG')
    
                # Construct flags per antenna, selecting by reason if desired
                if ants != '' or reason != 'any':
                    myflagcmd = selectXMLFlags(myflags, myantenna=ants,myreason=reason)
                else:
                    myflagcmd = myflags
    
#                listmode = 'online'
                listmode = 'xml'
                
            else:
                raise ValueError( 'Input type is not supported' )

            # Before performing any action on the flag cmds, check them! 
            vrows = list(myflagcmd.keys())
            if len(vrows) == 0:
                raise RuntimeError( 'There are no unapplied flags in the input. Set useapplied=True to also use the previously-applied flags.' )
            else:
                casalog.post('Read ' + str(len(vrows))
                             + ' lines from input')
                
            casalog.post('Flagcmd dictionary is: %s'%myflagcmd, 'DEBUG1')
              
            #
            # ACTION to perform on input file
            #
            casalog.post('Executing action = ' + action)
    
            # List the flag commands from inpfile on the screen
            # and save them or not to outfile
            if action == 'list':
    
                # List the flag cmds on the screen
                listFlagCommands(myflagcmd, listmode=listmode)

                # Save the flag cmds to the outfile
                if savepars:
                    # These cmds came from the internal FLAG_CMD, only list on the screen
                    if outfile == '':
                        if inpmode == 'table' and inpfile == '':
                            pass
                        else:
                            casalog.post('Saving commands to FLAG_CMD')
                            fh.writeFlagCommands(vis, myflagcmd, False, 
                                                 '', '', True)
                    elif not overwrite and os.path.exists(outfile):
                        raise RuntimeError( 'You have set overwrite to False. Remove %s before saving the flag commands'%outfile )

                    else:
                        casalog.post('Saving commands to ' + outfile)
                        fh.writeFlagCommands(vis, myflagcmd, False, '', outfile, False)
                        
            elif action == 'apply' or action == 'unapply':
    
                # Apply/Unapply the flag commands to the data
                apply = True

                # Select the data
    
                # Select a loose union of the data selection from the list.
                # The loose union will be calculated for field and spw only.
                # Antenna, correlation and timerange should be handled by the agent
                unionpars = {}
#                 if len(vrows) > 1:
#                     unionpars = fh.parseUnion(vis, myflagcmd)
#                     if len(unionpars.keys()) > 0:
#                         casalog.post('Pre-selecting a subset of the MS: ')
#                         casalog.post('%s' % unionpars)
#                     else:
#                         casalog.post('Iterating through the entire MS')
    
#                elif len(vrows) == 1:
                if len(vrows) == 1:
    
                # Get all the selection parameters, but set correlation to ''
                    # if the table was selected by row, we need to
                    # know what is the key number in the dictionary
                    cmd0 = myflagcmd[vrows[0]]['command']
                    unionpars = fh.parseSelectionPars(cmd0)
                    casalog.post('The selected subset of the MS will be: ')
                    casalog.post('%s' % unionpars)
    
                aflocal.selectdata(unionpars)
    
                # Parse the agents parameters
                if action == 'unapply':
                    apply = False
    
                casalog.post('Parse the parameters for the agents')
                fh.parseAgents(aflocal, myflagcmd, [], apply, True, '')
    
                # Initialize the Agents
                aflocal.init()
    
                # Backup the flags before running
                if flagbackup:
                    fh.backupFlags(aflocal, msfile='', prename='flagcmd')
    
                # Run the tool
                aflocal.run(True)
    
                aflocal.done()
    
                # Update the APPLIED column
                if savepars:
                    # These flags came from internal FLAG_CMD. Always update APPLIED
                    if outfile == '':
                        if inpmode == 'table' and inpfile == '':
                            updateTable(vis, mycol='APPLIED',
                                    myval=apply, myrowlist=vrows)
                        else:
                            # save to FLAG_CMD
                            casalog.post('Saving commands to FLAG_CMD')
                            fh.writeFlagCommands(vis, myflagcmd, apply, '', 
                                                 '', True)
                    # Save to a file
                    else:                    
                        # Still need to update APPLIED column
                        if inpmode == 'table' and inpfile == '':
                            updateTable(vis, mycol='APPLIED',
                                    myval=apply, myrowlist=vrows)

                        if not overwrite and os.path.exists(outfile):
                            raise RuntimeError( 'You have set overwrite to False. Remove %s before saving the flag commands'%outfile )

                        else:
                            casalog.post('Saving commands to file '+ outfile)
                            fh.writeFlagCommands(vis, myflagcmd, apply, '', outfile, False)
                            
                # Do not save cmds but maybe update APPLIED
                else:
                    if inpmode == 'table' and inpfile == '':
                        updateTable(vis, mycol='APPLIED', myval=apply,
                                    myrowlist=vrows)
                    
            elif action == 'plot':
    
                keylist = myflagcmd.keys()
                if len(keylist) > 0:
                    # Plot flag commands from FLAG_CMD or xml
                    casalog.post('Warning: will only reliably plot individual per-antenna flags'
                                 )
                    fns = newplotflags(myflagcmd, plotfile, t1sdata, t2sdata)
                    return {'plotfiles': fns}

                else:
                    casalog.post('Warning: empty flag dictionary, nothing to plot'
                                 )
            elif action == 'extract':
                # Make the command dictionary a string again
                outdict = copy.deepcopy(myflagcmd)
                for key in myflagcmd.keys():
                    cmddict = myflagcmd[key]['command']
                    cmdline = ""
                    for k,v in lociteritems(cmddict):
                        cmdline = cmdline + k + '=' + str(v) + ' '
                    cmdline.rstrip()
                    outdict[key]['command'] = cmdline
                    
                casalog.post('Returning extracted dictionary')
                return outdict
            
    
    finally:
        aflocal.done()
        
    # Write history only to action='apply' or 'unapply'
    # write history
    if not iscal and (action == 'apply' or action == 'unapply'):
        try:
            param_names = flagcmd.__code__.co_varnames[:flagcmd.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]

            write_history(mslocal, vis, 'flagcmd', param_names,
                          param_vals, casalog)

        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                         'WARN')


# ************************************************************************
#                    Helper Functions
# ************************************************************************


def readFromTable(
    msfile,
    myflagrows=[],
    useapplied=True,
    myreason='any',
    ):
    '''Read flag commands from rows of the FLAG_CMD table of msfile
    If useapplied=False then include only rows with APPLIED=False
    If myreason is anything other than '', then select on that'''

    #
    # Return flagcmd structure:
    #
    # The flagcmd key is the integer row number from FLAG_CMD
    #
    #   Dictionary structure:
    #   key : 'id' (string)
    #         'mode' (string)         flag mode '','clip','shadow','quack'
    #         'antenna' (string)
    #         'timerange' (string)
    #         'reason' (string)
    #         'time' (float)          in mjd seconds
    #         'interval' (float)      in mjd seconds
    #         'cmd' (string)          string (for COMMAND col in FLAG_CMD)
    #         'type' (string)         'FLAG' / 'UNFLAG'
    #         'applied' (bool)        set to True here on read-in
    #         'level' (int)           set to 0 here on read-in
    #         'severity' (int)        set to 0 here on read-in

    # Open and read columns from FLAG_CMD
    mstable = os.path.join(msfile,'FLAG_CMD')

    # Note, tb.getcol doesn't allow random row access, read all

    try:
        tblocal.open(mstable)
        f_time = tblocal.getcol('TIME')
        f_interval = tblocal.getcol('INTERVAL')
        f_type = tblocal.getcol('TYPE')
        f_reas = tblocal.getcol('REASON')
        f_level = tblocal.getcol('LEVEL')
        f_severity = tblocal.getcol('SEVERITY')
        f_applied = tblocal.getcol('APPLIED')
        f_cmd = tblocal.getcol('COMMAND')
        tblocal.close()
    except:
        casalog.post('Error reading table ' + mstable, 'ERROR')
        raise Exception

    nrows = len(f_time)

    myreaslist = []

    # Parse myreason
    if type(myreason) == str:
        if myreason != 'any':
            myreaslist.append(myreason)
    elif type(myreason) == list:
        myreaslist = myreason
    else:
        casalog.post('Cannot read reason; it contains unknown variable types'
                     , 'ERROR')
        return
    
    myflagcmd = {}

    if nrows > 0:
        nflagd = 0
        if len(myflagrows) > 0:
            rowlist = myflagrows
        else:
            rowlist = list(range(nrows))
        # Prune rows if needed
        if not useapplied:
            rowl = []
            for i in rowlist:
                if not f_applied[i]:
                    rowl.append(i)
            rowlist = rowl
        if len(myreaslist) > 0:
            rowl = []
            for i in rowlist:
                if myreaslist.count(f_reas[i]) > 0:
                    rowl.append(i)
            rowlist = rowl
            
        # Define the way to parse the strings
        myParser = fh.Parser(' ', '=')
        
        for i in rowlist:
            flagd = {}
            cmd = f_cmd[i]
            if cmd == '':
                casalog.post('Ignoring empty COMMAND string', 'WARN')
                continue

            # Extract antenna and timerange strings from cmd
            flagd['id'] = str(i)
            flagd['antenna'] = ''
            flagd['mode'] = ''
            flagd['timerange'] = ''
            flagd['time'] = f_time[i]
            flagd['interval'] = f_interval[i]
            flagd['type'] = f_type[i]
            flagd['reason'] = f_reas[i]
            flagd['level'] = f_level[i]
            flagd['severity'] = f_severity[i]
            flagd['applied'] = f_applied[i]
            
            # If shadow, remove the addantenna dictionary
            if cmd.__contains__('shadow') \
                and cmd.__contains__('addantenna'):
                i0 = cmd.rfind('addantenna')
                if cmd[i0 + 11] == '{':
                    # It is a dictionary. Remove it from line
                    i1 = cmd.rfind('}')
                    antpar = cmd[i0 + 11:i1 + 1]
                    temp = cmd[i0:i1 + 1]
                    newcmd = cmd.replace(temp, '')
                    antpardict = fh.convertStringToDict(antpar)
                    flagd['addantenna'] = antpardict
                    cmd = newcmd

            # Get a dictionary without type evaluation
            preparsing = myParser.parseNoEval(cmd)
            
            # Evaluate the types into a new dictionary
            parsed = fh.evaluateParameters(preparsing)

            flagd['command'] = parsed
            
            if 'mode' in parsed:
                flagd['mode'] = parsed['mode']
            if 'timerange' in parsed:
                flagd['timerange'] = parsed['timerange']
            if 'antenna' in parsed:
                flagd['antenna'] = parsed['antenna']
            if 'id' in parsed:
                flagd['id'] = parsed['id']

            # Keep original key index, might need this later
            myflagcmd[i] = flagd
            nflagd += 1

    else:
        casalog.post('FLAG_CMD table in %s is empty, no flags extracted'
                      % msfile, 'WARN')

    return myflagcmd

def readFromCmd(cmdlist, ms_startmjds, ms_endmjds):
    '''Read the parameters from a list of commands'''

    # Read a list of strings and return a dictionary of parameters
    myflagd = {}
    nrows = len(cmdlist)
    if nrows == 0:
        casalog.post('WARNING: empty flag command list', 'WARN')
        return myflagd

    t = qalocal.quantity(ms_startmjds, 's')
    ms_startdate = qalocal.time(t, form=['ymd', 'no_time'])[0]
    t0 = qalocal.totime(ms_startdate + '/00:00:00.0')
    # t0d = qalocal.convert(t0,'d')
    t0s = qalocal.convert(t0, 's')

    ncmds = 0
    for i in range(nrows):
        cmdstr = cmdlist[i]
        # break string into key=val sets
        keyvlist = cmdstr.split()
        if len(keyvlist) > 0:
            ant = ''
            timstr = ''
            tim = 0.5 * (ms_startmjds + ms_endmjds)
            intvl = ms_endmjds - ms_startmjds
            reas = ''
            cmd = ''
            fid = str(i)
            typ = 'FLAG'
            appl = False
            levl = 0
            sevr = 0
            fmode = ''
            for keyv in keyvlist:
                # check for comment character #
                if keyv.count('#') > 0:
                    # break out of loop parsing keyvals
                    break
                try:
                    (xkey, val) = keyv.split('=')
                except:
                    casalog.post('Not a key=val pair: ' + keyv, 'WARN')
                    break
                xval = val
                # Use eval to deal with conversion from string
                # xval = eval(val)
                # strip quotes from value (if still a string)
                if type(xval) == str:
                    if xval.count("'") > 0:
                        xval = xval.strip("'")
                    if xval.count('"') > 0:
                        xval = xval.strip('"')

                # Strip these out of command string
                if xkey == 'reason':
                    reas = xval
                elif xkey == 'applied':
                    appl = False
                    if xval == 'True':
                        appl = True
                elif xkey == 'level':
                    levl = int(xval)
                elif xkey == 'severity':
                    sevr = int(xval)
                elif xkey == 'time':
                    tim = xval
                elif xkey == 'interval':
                    intvl = xval
                else:
                    # Extract (but keep in string)
                    if xkey == 'timerange':
                        timstr = xval
                        # Extract TIME,INTERVAL
                        try:
                            (startstr, endstr) = timstr.split('~')
                        except:
                            if timstr.count('~') == 0:
                            # casalog.post('Assuming a single start time ')
                                startstr = timstr
                                endstr = timstr
                            else:
                                raise Exception("too may ~'s. Not a start~end range. Error "
                                                "parsing " + timstr )
                        t = qalocal.totime(startstr)
                        starts = qalocal.convert(t, 's')
                        if starts['value'] < 1.E6:
                            # assume a time offset from ref
                            starts = qalocal.add(t0s, starts)
                        startmjds = starts['value']
                        if endstr == '':
                            endstr = startstr
                        t = qalocal.totime(endstr)
                        ends = qalocal.convert(t, 's')
                        if ends['value'] < 1.E6:
                            # assume a time offset from ref
                            ends = qalocal.add(t0s, ends)
                        endmjds = ends['value']
                        tim = 0.5 * (startmjds + endmjds)
                        intvl = endmjds - startmjds
                    elif xkey == 'antenna':

                        ant = xval
                    elif xkey == 'id':
                        fid = xval
                    elif xkey == 'unflag':
                        if xval == 'True':
                            typ = 'UNFLAG'
                    elif xkey == 'mode':
                        fmode = xval
                    cmd = cmd + ' ' + keyv
            # Done parsing keyvals
            # Make sure there is a non-blank command string after reason/id extraction
            if cmd != '':
                flagd = {}
                flagd['id'] = fid
                flagd['mode'] = fmode
                flagd['antenna'] = ant
                flagd['timerange'] = timstr
                flagd['reason'] = reas
                flagd['command'] = cmd
                flagd['time'] = tim
                flagd['interval'] = intvl
                flagd['type'] = typ
                flagd['level'] = levl
                flagd['severity'] = sevr
                flagd['applied'] = appl
                # Insert into main dictionary
                myflagd[ncmds] = flagd
                ncmds += 1

    casalog.post('Parsed ' + str(ncmds) + ' flag command strings')

    return myflagd


def readFromFile(
    cmdlist,
    ms_startmjds,
    ms_endmjds,
    myreason='',
    ):
    '''Parse list of flag command strings and return dictionary of flagcmds
    Inputs:
       cmdlist (list,string) list of command strings (default for TIME,INTERVAL)
       ms_startmjds (float)  starting mjd (sec) of MS (default for TIME,INTERVAL)
       ms_endmjds (float)    ending mjd (sec) of MS'''

#
#   Usage: myflagcmd = getflags(cmdlist)
#
#   Dictionary structure:
#   fid : 'id' (string)
#         'mode' (string)         flag mode '','clip','shadow','quack'
#         'antenna' (string)
#         'timerange' (string)
#         'reason' (string)
#         'time' (float)          in mjd seconds
#         'interval' (float)      in mjd seconds
#         'cmd' (string)          string (for COMMAND col in FLAG_CMD)
#         'type' (string)         'FLAG' / 'UNFLAG'
#         'applied' (bool)        set to True here on read-in
#         'level' (int)           set to 0 here on read-in
#         'severity' (int)        set to 0 here on read-in
#
# v3.2 Updated STM 2010-12-03 (3.2.0) handle comments # again
# v3.2 Updated STM 2010-12-08 (3.2.0) bug fixes in flagsort use, parsing
# v3.3 Updated STM 2010-12-20 (3.2.0) bug fixes parsing errors
#

    myflagd = {}
    nrows = len(cmdlist)
    if nrows == 0:
        casalog.post('WARNING: empty flag command list', 'WARN')
        return myflagd

    # Parse the reason
    reasonlist = []
    if type(myreason) == str:
        if myreason != '':
            reasonlist.append(myreason)
    elif type(myreason) == list:
        reasonlist = myreason
    else:
        casalog.post('Cannot read reason; it contains unknown variable types'
                     , 'ERROR')
        return

    t = qalocal.quantity(ms_startmjds, 's')
    ms_startdate = qalocal.time(t, form=['ymd', 'no_time'])[0]
    t0 = qalocal.totime(ms_startdate + '/00:00:00.0')
    # t0d = qalocal.convert(t0,'d')
    t0s = qalocal.convert(t0, 's')

    rowlist = []

    # IMPLEMENT THIS LATER
    # First select by reason. Simple selection...
#    if len(reasonlist) > 0:
#        for i in range(nrows):
#            cmdstr = cmdlist[i]
#            keyvlist = cmdstr.split()
#            if len(keyvlist) > 0:
#                for keyv in keyvlist:
#                    (xkey, xval) = keyv.split('=')
#
#                    if type(xval) == str:
#                        if xval.count("'") > 0:
#                            xval = xval.strip("'")
#                        if xval.count('"') > 0:
#                            xval = xval.strip('"')
#
#                    if xkey == 'reason':
#                        if reasonlist.count(xval) > 0
#    else:
    rowlist = list(range(nrows))

    # Now read the only the commands from the file that satisfies the reason selection

    ncmds = 0
#    for i in range(nrows):
    for i in rowlist:
        cmdstr = cmdlist[i]

        # break string into key=val sets
        keyvlist = cmdstr.split()
        if len(keyvlist) > 0:
            ant = ''
            timstr = ''
            tim = 0.5 * (ms_startmjds + ms_endmjds)
            intvl = ms_endmjds - ms_startmjds
            reas = ''
            cmd = ''
            fid = str(i)
            typ = 'FLAG'
            appl = False
            levl = 0
            sevr = 0
            fmode = ''
            for keyv in keyvlist:
                # check for comment character #
                if keyv.count('#') > 0:
                    # break out of loop parsing keyvals
                    break
                try:
                    (xkey, val) = keyv.split('=')
                except:
                    casalog.post('Not a key=val pair: ' + keyv, "WARN")
                    break
                xval = val
                # Use eval to deal with conversion from string
                # xval = eval(val)
                # strip quotes from value (if still a string)
                if type(xval) == str:
                    if xval.count("'") > 0:
                        xval = xval.strip("'")
                    if xval.count('"') > 0:
                        xval = xval.strip('"')

                # Strip these out of command string
                if xkey == 'reason':
                    reas = xval
                elif xkey == 'applied':

                    appl = False
                    if xval == 'True':
                        appl = True
                elif xkey == 'level':
                    levl = int(xval)
                elif xkey == 'severity':
                    sevr = int(xval)
                elif xkey == 'time':
                    tim = xval
                elif xkey == 'interval':
                    intvl = xval
                else:
                    # Extract (but keep in string)
                    if xkey == 'timerange':
                        timstr = xval
                        # Extract TIME,INTERVAL
                        try:
                            (startstr, endstr) = timstr.split('~')
                        except:
                            if timstr.count('~') == 0:
                            # casalog.post('Assuming a single start time ')
                                startstr = timstr
                                endstr = timstr
                            else:
                                raise Exception( "too may ~'s. Not a start~end range. Error "
                                                 "parsing " + timstr )
                        t = qalocal.totime(startstr)
                        starts = qalocal.convert(t, 's')
                        if starts['value'] < 1.E6:
                            # assume a time offset from ref
                            starts = qalocal.add(t0s, starts)
                        startmjds = starts['value']
                        if endstr == '':
                            endstr = startstr
                        t = qalocal.totime(endstr)
                        ends = qalocal.convert(t, 's')
                        if ends['value'] < 1.E6:
                            # assume a time offset from ref
                            ends = qalocal.add(t0s, ends)
                        endmjds = ends['value']
                        tim = 0.5 * (startmjds + endmjds)
                        intvl = endmjds - startmjds
                    elif xkey == 'antenna':

                        ant = xval
                    elif xkey == 'id':
                        fid = xval
                    elif xkey == 'unflag':
                        if xval == 'True':
                            typ = 'UNFLAG'
                    elif xkey == 'mode':
                        fmode = xval
                    cmd = cmd + ' ' + keyv
            # Done parsing keyvals
            # Make sure there is a non-blank command string after reason/id extraction
            if cmd != '':
                flagd = {}
                flagd['id'] = fid
                flagd['mode'] = fmode
                flagd['antenna'] = ant
                flagd['timerange'] = timstr
                flagd['reason'] = reas
                flagd['command'] = cmd
                flagd['time'] = tim
                flagd['interval'] = intvl
                flagd['type'] = typ
                flagd['level'] = levl
                flagd['severity'] = sevr
                flagd['applied'] = appl
                # Insert into main dictionary
                myflagd[ncmds] = flagd
                ncmds += 1

    casalog.post('Parsed ' + str(ncmds) + ' flag command strings')

    return myflagd


def updateTable(
    msfile,
    mycol='',
    myval=None,
    myrowlist=[],
    ):
    '''Update commands in myrowlist of the FLAG_CMD table of msfile    
       Usage: updateflagcmd(msfile,myrow,mycol,myval)'''

    # Example:
    #
    #    updateflagcmd(msfile,mycol='APPLIED',myval=True)
    #       Mark all rows as APPLIED=True
    #
    #    updateflagcmd(msfile,mycol='APPLIED',myval=True,myrowlist=[0,1,2])
    #       Mark rows 0,1,2 as APPLIED=True
    #

    if mycol == '':
        casalog.post('WARNING: No column to was specified to update; doing nothing'
                     , 'WARN')
        return

    # Open and read columns from FLAG_CMD
    mstable = os.path.join(msfile,'FLAG_CMD')
    try:
        tblocal.open(mstable, nomodify=False)
    except:
        raise Exception( 'Error opening table ' + mstable )

    nrows = int(tblocal.nrows())

    # Check against allowed colnames
    colnames = tblocal.colnames()
    if colnames.count(mycol) < 1:
        casalog.post('Error: column mycol=' + mycol + ' not one of: '
                     + str(colnames))
        return

    nlist = len(myrowlist)

    if nlist > 0:
        rowlist = myrowlist
    else:

        rowlist = list(range(nrows))
        nlist = nrows

    if nlist > 0:
        try:
            tblocal.putcell(mycol, rowlist, myval)
        except:
            raise Exception( 'Error updating FLAG_CMD column ' + mycol \
                + ' to value ' + str(myval) )

        casalog.post('Updated ' + str(nlist)
                     + ' rows of FLAG_CMD table in MS')
    tblocal.close()

def listFlagCommands(myflags=None, listmode=''):
    '''List flags from MS or a file. The flags are read from
       a dictionary created by fh.parseDictionary()
    Format according to listmode:
        =''          do nothing
        ='list'      Format for flag command strings
        ='table'     Format for FLAG_CMD flags
        ='xml'    Format for online flags'''

    #
    #   Dictionary structure:
    #   fid : 'id' (string)
    #         'mode' (string)         flag mode '','clip','shadow','quack'
    #         'antenna' (string)
    #         'timerange' (string)
    #         'reason' (string)
    #         'time' (float)          in mjd seconds
    #         'interval' (float)      in mjd seconds
    #         'cmd' (string)          string (for COMMAND col in FLAG_CMD)
    #         'type' (string)         'FLAG' / 'UNFLAG'
    #         'applied' (bool)
    #         'level' (int)
    #         'severity' (int)
    #


    # list to logger and screen
    if listmode == 'table':
        phdr = '%-8s %-32s %-7s %s' % (
            'Row',
            'Reason',
            'Applied',
            'Command'
            )
        casalog.post(phdr)
        mydash=80*'-'
        casalog.post('%-80s'%mydash)
        for k in myflags.keys():
#            time = myflags[k]['TIME']
            row = myflags[k]['id']
            reason = myflags[k]['reason']
            applied = myflags[k]['applied']
            
            cmddict = myflags[k]['command']
            cmdline = ""
            for key,val in lociteritems(cmddict):
                cmdstr = ""
                if isinstance(val, str):
                    # Add quotes to string values
                    cmdstr = "'"+val+"'"
                    val = cmdstr
                cmdline = cmdline + key + '=' + str(val) + ' '
            
            # Print to logger
            pstr = '%-8s %-32s %-7s %s' % (
                row, reason,applied,cmdline)
            casalog.post(pstr)

    elif listmode == 'list':
        phdr = '%-8s %-32s %s' % ('Key', 'Reason', 'Command')
        casalog.post(phdr)
        mydash=80*'-'
        casalog.post('%-80s'%mydash)
        for k in myflags.keys():
            cmddict = myflags[k]['command']
            reason = myflags[k]['reason']
            if 'reason' in cmddict:
                cmddict.pop('reason')
                
            cmdline = ""
            for key,val in lociteritems(cmddict):
                cmdstr = ""
                if isinstance(val, str):
                    # Add quotes to string values
                    cmdstr = "'"+val+"'"
                    val = cmdstr
                cmdline = cmdline + key + '=' + str(val) + ' '
            
            # Print to logger
            pstr = '%-8s %-32s %s' % (k, reason, cmdline)
            casalog.post(pstr)
            
    elif listmode == 'xml':
        phdr = '%-8s %-8s %-48s %-32s' % ('Key', 'Antenna',
                'Timerange', 'Reason')
        casalog.post(phdr)
        mydash=80*'-'
        casalog.post('%-80s'%mydash)
        for k in myflags.keys():
            reason = ''
            antenna = ''
            timerange = ''
            cmddict = myflags[k]
            if 'reason' in cmddict:
                reason = cmddict['reason']
            if 'antenna' in cmddict:
                antenna = cmddict['antenna']
            if 'timerange' in cmddict:
                timerange = cmddict['timerange']
                
            pstr = '%-8s %-8s %-48s %-32s' % (k, antenna,timerange,reason)
            casalog.post(pstr)

    return

def listFlagCmd(
    myflags=None,
    myantenna='',
    myreason='',
    myoutfile='',
    listmode='',
    ):
    '''List flags in myflags dictionary
    
    Format according to listmode:
        =''          do nothing
        ='file'      Format for flag command strings
        ='cmd'       Format for FLAG_CMD flags
        ='online'    Format for online flags'''

    #
    #   Dictionary structure:
    #   fid : 'id' (string)
    #         'mode' (string)         flag mode '','clip','shadow','quack'
    #         'antenna' (string)
    #         'timerange' (string)
    #         'reason' (string)
    #         'time' (float)          in mjd seconds
    #         'interval' (float)      in mjd seconds
    #         'cmd' (string)          string (for COMMAND col in FLAG_CMD)
    #         'type' (string)         'FLAG' / 'UNFLAG'
    #         'applied' (bool)
    #         'level' (int)
    #         'severity' (int)
    #

    useid = False

    if myoutfile != '':
        try:
            lfout = open(myoutfile, 'w')
        except:
            raise Exception( 'Error opening list output file ' \
                + myoutfile )

    keylist = myflags.keys()
    if len(keylist) == 0:
        casalog.post('There are no flags to list', 'WARN')
        return
    # Sort keys
#    keylist.sort

    # Set up any selection
    if myantenna != '':
        casalog.post('Selecting flags by antenna="' + str(myantenna)
                     + '"')
    myantlist = myantenna.split(',')

    if myreason != '':
        casalog.post('Selecting flags by reason="' + str(myreason) + '"'
                     )
    myreaslist = myreason.split(',')

    if listmode == 'online':
        phdr = '%8s %12s %8s %32s %48s' % ('Key', 'FlagID', 'Antenna',
                'Reason', 'Timerange')
    elif listmode == 'cmd':
        phdr = '%8s %45s %32s %6s %7s %3s %3s %s' % (
            'Row',
            'Timerange',
            'Reason',
            'Type',
            'Applied',
            'Level',
            'Severity',
            'Command',
            )
    elif listmode == 'file':
        phdr = '%8s %32s %s' % ('Key', 'Reason', 'Command')
    else:
        return

    if myoutfile != '':
        # list to output file
        lfout.write(phdr + '\n')
    else:
        # list to logger
        casalog.post(phdr)

    # Loop over flags
    for key in keylist:
        fld = myflags[key]
        # Get fields
        skey = str(key)
        if 'id' in fld and useid:
            fid = fld['id']
        else:
            fid = str(key)
        if 'antenna' in fld:
            ant = fld['antenna']
        else:
            ant = 'Unset'
        if 'timerange' in fld:
            timr = fld['timerange']
        else:
            timr = 'Unset'
        if 'reason' in fld:
            reas = fld['reason']
        else:
            reas = 'Unset'
        if 'command' in fld:
            cmd = fld['command']
            # To be verified
#            if 'addantenna' in fld:
#                addantenna = fld['addantenna']
#                cmd = cmd + ' addantenna=' + str(addantenna)
#        else:

#            cmd = 'Unset'
        if 'type' in fld:
            typ = fld['type']
        else:
            typ = 'FLAG'
        if 'level' in fld:
            levl = str(fld['level'])
        else:
            levl = '0'
        if 'severity' in fld:
            sevr = str(fld['severity'])
        else:
            sevr = '0'
        if 'applied' in fld:
            appl = str(fld['applied'])
        else:
            appl = 'Unset'

        # Print out listing
        if myantenna == '' or myantlist.count(ant) > 0:
            if myreason == '' or myreaslist.count(reas) > 0:
                if listmode == 'online':
                    pstr = '%8s %12s %8s %32s %48s' % (skey, fid, ant,
                            reas, timr)
                elif listmode == 'cmd':
                    # Loop over dictionary with commands
                    cmdline = ""
                    for k,v in lociteritems(cmd):
                        cmdline = cmdline + k + '=' + str(v) + ' '
                    
                    cmdline = cmdline.rstrip()                       
                    pstr = '%8s %45s %32s %6s %7s %3s %3s %s' % (
                        skey,
                        timr,
                        reas,
                        typ,
                        appl,
                        levl,
                        sevr,
                        cmdline,
                        )
                else:
                    cmdline = ""
                    for k,v in lociteritems(cmd):
                        cmdline = cmdline + k + '=' + str(v) + ' '
                    
                    cmdline = cmdline.rstrip()                       
                    pstr = '%8s %45s %32s %6s %7s %3s %3s %s' % (
                        skey,
                        timr,
                        reas,
                        typ,
                        appl,
                        levl,
                        sevr,
                        cmdline,
                        )
                    pstr = '%8s %32s %s' % (skey, reas, cmdline)
                if myoutfile != '':
                    # list to output file
                    lfout.write(pstr + '\n')
                else:
                    # list to logger
                    casalog.post(pstr)
    if myoutfile != '':
        lfout.close()


def selectXMLFlags(
    myflags=None,
    myantenna='',
    myreason='any',
    ):

    #
    # Return dictionary of input flags using selection by antenna/reason
    # and grouped/sorted by flagsort.
    #
    #   selectFlags: Return dictionary of flags using selection by antenna/reason
    #              and grouped/sorted by flagsort.
    #      myflags (dictionary)  input flag dictionary (e.g. from readflagxml
    #      myantenna (string)    selection by antenna(s)
    #      myreason (string)     selection by reason(s)
    #
    #   Usage: myflagd = selectFlags(myflags,antenna,reason)
    #
    #   Dictionary structure:
    #   fid : 'id' (string)
    #         'mode' (string)         flag mode '','clip','shadow','quack','online'
    #         'antenna' (string)
    #         'timerange' (string)
    #         'reason' (string)
    #         'time' (float)          in mjd seconds
    #         'interval' (float)      in mjd seconds
    #         'cmd' (string)          command string (for COMMAND col in FLAG_CMD)
    #         'type' (string)         'FLAG' / 'UNFLAG'
    #         'applied' (bool)        set to True here on read-in
    #         'level' (int)           set to 0 here on read-in
    #         'severity' (int)        set to 0 here on read-in
    #
    #
    #
    # Check if any operation is needed
    if myantenna == '' and myreason == '':
        casalog.post('No selection or sorting needed - sortflags returning input dictionary')
        flagd = myflags
        return flagd
    #
    flagd = {}
    nflagd = 0
    keylist = myflags.keys()
    casalog.post('Selecting from ' + str(len(keylist)) \
                 + ' flagging commands')
        
    if len(keylist) == 0:
        casalog.post('No flags found in input dictionary')
        return myflags

    #
    # Sort by key
    #
    keylist.sort()
    #
    # Construct flag command list for selected ant,reason
    #
    casalog.post('Selecting flags by antenna="' + str(myantenna)
                 + '"')
    myantlist = myantenna.split(',')

    casalog.post('Selecting flags by reason="' + str(myreason) + '"'
                 )
    myreaslist = []
# Parse myreason
    if type(myreason) == str:
        if myreason == '':
            casalog.post('WARNING: reason= is treated as selection on a blank REASON!'
                         , 'WARN')
        if myreason != 'any':
            myreaslist.append(myreason)
    elif type(myreason) == list:
        myreaslist = myreason
    else:
        casalog.post('ERROR: reason contains unknown variable type'
                     , 'SEVERE')
        return
    if len(myreaslist) > 0:
        casalog.post('Selecting for reasons: ' + str(myreaslist))

# Note antenna and reason selection checks for inclusion not exclusivity
    doselect = myantenna != '' or len(myreaslist) > 0

# Now loop over flags, break into sorted and unsorted groups
    nunsortd = 0
    unsortd = {}
    unsortdlist = []
    
# All flags are in unsorted list
    unsortd = myflags.copy()
    unsortdlist = unsortd.keys()
    nunsortd = len(unsortdlist)

# selection on unsorted flags
    if doselect and nunsortd > 0:
        keylist = unsortd.keys()
        for key in keylist:
            myd = unsortd[key]
            ant = myd['antenna']
            antlist = ant.split(',')
            reas = myd['reason']
            reaslist = reas.split(',')
    # break this flag by antenna
            antstr = ''
            reastr = ''
            addf = False

            for a in antlist:
                if myantenna == '' or myantlist.count(a) > 0:
                    addr = False
                    if len(myreaslist) > 0:
                        for r in myreaslist:
                            if reas == r or reaslist.count(r) > 0:
                                addr = True
                # check if this is a new reason
                                rlist = reastr.split(',')
                                if reastr != '':
                                    rlist = reastr.split(',')
                                    if rlist.count(r) == 0:
                                        reastr += ',' + r
                                else:
                                    reastr = r
                    else:
                        addr = True
                        reastr = reas
                    if addr:
                        addf = True
                        if antstr != '':
                # check if this is a new antenna
                            alist = antstr.split(',')
                            if alist.count(a) == 0:
                                antstr += ',' + a
                        else:
                            antstr = a
            if addf:
                flagd[nflagd] = myd
                flagd[nflagd]['antenna'] = antstr
                flagd[nflagd]['reason'] = reastr
                nflagd += 1
        flagdlist = flagd.keys()
    elif nunsortd > 0:
    # just copy to flagd w/o selection
        flagd = unsortd.copy()
        flagdlist = flagd.keys()
        nflagd = len(flagdlist)

    if nflagd > 0:
        casalog.post('Found total of ' + str(nflagd)
                     + ' flags meeting selection criteria')
    else:
        casalog.post('No flagging commands found meeting criteria')

    return flagd


# Done


def clearFlagCmd(msfile, myrowlist=[]):
    #
    # Delete flag commands (rows) from the FLAG_CMD table of msfile
    #
    # Open and read columns from FLAG_CMD

    mstable = os.path.join(msfile,'FLAG_CMD')
    try:
        tblocal.open(mstable, nomodify=False)
    except:
        raise Exception( 'Error opening table ' + mstable )

    nrows = int(tblocal.nrows())
    casalog.post('There were ' + str(nrows) + ' rows in FLAG_CMD')
    if nrows > 0:
        if len(myrowlist) > 0:
            rowlist = myrowlist
        else:
            rowlist = list(range(nrows))
        try:
            tblocal.removerows(rowlist)
            casalog.post('Deleted ' + str(len(rowlist))
                         + ' from FLAG_CMD table in MS')
        except:
            tblocal.close()
            raise Exception( 'Error removing rows ' + str(rowlist) \
                + ' from table ' + mstable )

    else:
        casalog.post('No rows to clear')

    tblocal.close()

def newplotflags(
    myflags,
    plotname,
    t1sdata,
    t2sdata,
    ):
    #
    # Function to plot flagging dictionary
    # Adapted from J.Marvil
    # Updated STM v4.1 2011-11-02 to handle ALMA flags
    # Updated STM v4.2 2012-02-16 trim flag times to data times
    # Updated STM v4.2 2012-04-10 bug fix in trim flag times to data times
    # Updated STM v4.2 2012-04-10 messages to logger

    # After the swig converstion, it seems that the following
    # line is not needed anymore
    #    qa = casac.qa = qatool = casac.quanta()

    try:
        import pylab as pl
    except ImportError as exc:
        raise ImportError('Failed to load pylab, required by flagcmd: {}'.format(exc))

    # list of supported colors (in order)
    colorlist = [
        'red',
        'blue',
        'green',
        'black',
        'cyan',
        'magenta',
        'yellow',
        'orange',
        ]
    ncolors = len(colorlist)

    # get list of flag keys
    keylist = myflags.keys()

    # get lists of antennas and reasons
    # make plotting dictionary
    myants = []
    myreas = []
    plotflag = {}
    ipf = 0
    for key in keylist:
        antstr = myflags[key]['antenna']
        reastr = myflags[key]['reason']
        timstr = myflags[key]['timerange']
        if antstr != '':
            # flags that have antenna specified
            antlist = antstr.split(',')
            nantlist = len(antlist)
        else:
            # Special
            antlist = ['All']
            nantlist = 1
        #
        realist = reastr.split(',')
        nrealist = len(realist)
        #
        timlist = timstr.split(',')
        ntimlist = len(timlist)
        #
        # Break these into nants x ntimes flags
        # Trick is assigning multiple reasons
        # Normal cases:
        # A. One reason, single/multiple antennas x times
        # B. Multiple reasons=times, single/multiple antenna(s)
        # C. Single time, multiple antennas/reasons
        # D. Multiple reasons, no way to correspond with times
        #
        timmin = 1.0E11
        timmax = 0.0
        if nrealist == 1:
            # simplest case, single reason
            reas = realist[0]
            if reas == '':
                reas = 'Unknown'
            if myreas.count(reas) == 0:
                myreas.append(reas)
            for ia in range(nantlist):
                ant = antlist[ia]
                if myants.count(ant) == 0:
                    myants.append(ant)
                for it in range(ntimlist):
                    times = timlist[it]
                    plotflag[ipf] = {}
                    plotflag[ipf]['antenna'] = ant
                    plotflag[ipf]['reason'] = reas
                    plotflag[ipf]['timerange'] = times
                    plotflag[ipf]['show'] = True
                    ipf += 1
        elif nrealist == ntimlist:
            # corresponding reasons and times
            for ia in range(nantlist):
                ant = antlist[ia]
                if myants.count(ant) == 0:
                    myants.append(ant)
                for it in range(ntimlist):
                    times = timlist[it]
                    reas = realist[it]
                    if reas == '':
                        reas = 'Unknown'
                    if myreas.count(reas) == 0:
                        myreas.append(reas)
                    plotflag[ipf] = {}
                    plotflag[ipf]['antenna'] = ant
                    plotflag[ipf]['reason'] = reas
                    plotflag[ipf]['timerange'] = times
                    plotflag[ipf]['show'] = True
                    ipf += 1
        else:
            # no correspondence between multiple reasons and ants/times
            # assign reason 'Miscellaneous'
            reas = 'Miscellaneous'
            if myreas.count(reas) == 0:
                myreas.append(reas)
            for ia in range(nantlist):
                ant = antlist[ia]
                if myants.count(ant) == 0:
                    myants.append(ant)
                for it in range(ntimlist):
                    times = timlist[it]
                    plotflag[ipf] = {}
                    plotflag[ipf]['antenna'] = ant
                    plotflag[ipf]['reason'] = reas
                    plotflag[ipf]['timerange'] = times
                    plotflag[ipf]['show'] = True
                    ipf += 1

    myants.sort()
    nants = len(myants)
    nreas = len(myreas)
    casalog.post('Found ' + str(nreas) + ' reasons to plot for '
                 + str(nants) + ' antennas')
    npf = ipf
    casalog.post('Found ' + str(npf) + ' total flag ranges to plot')

    # sort out times
    for ipf in range(npf):
        times = plotflag[ipf]['timerange']
        if times != '':
            if times.count('~') > 0:
                t1 = times[:times.find('~')]
                t2 = times[times.find('~') + 1:]
            else:
                t1 = times
                t2 = t1
            (t1s, t2s) = (qalocal.convert(t1, 's')['value'], qalocal.convert(t2,
                          's')['value'])
            plotflag[ipf]['t1s'] = t1s
            plotflag[ipf]['t2s'] = t2s
            if t1s < timmin:
                timmin = t1s
            if t2s > timmax:
                timmax = t2s
    # min,max times
    q1 = qalocal.quantity(timmin, 's')
    time1 = qalocal.time(q1, form='ymd', prec=9)[0]
    q2 = qalocal.quantity(timmax, 's')
    time2 = qalocal.time(q2, form='ymd', prec=9)[0]
    casalog.post('Found flag times from ' + time1 + ' to ' + time2)

    # sort out blank times
    for ipf in range(npf):
        times = plotflag[ipf]['timerange']
        if times == '':
            if t2sdata >= t1sdata > 0:
                plotflag[ipf]['t1s'] = t1sdata
                plotflag[ipf]['t2s'] = t2sdata
            else:
                plotflag[ipf]['t1s'] = timmin
                plotflag[ipf]['t2s'] = timmax

    # if flag times are beyond range of data, trim them
    # Added STM 2012-02-16, fixed STM 2012-04-10
    ndropped = 0
    if t2sdata >= t1sdata > 0 and (timmin < t1sdata or timmax
                                   > t2sdata):
        # min,max data times
        q1 = qalocal.quantity(t1sdata, 's')
        tdata1 = qalocal.time(q1, form='ymd', prec=9)[0]
        q2 = qalocal.quantity(t2sdata, 's')
        tdata2 = qalocal.time(q2, form='ymd', prec=9)[0]
        casalog.post('WARNING: Trimming flag times to data limits '
                     + tdata1 + ' to ' + tdata2)

        for ipf in range(npf):
            t1s = plotflag[ipf]['t1s']
            t2s = plotflag[ipf]['t2s']
            if t1s < t1sdata:
                if t2s >= t1sdata:
                    # truncate to t1sdata
                    plotflag[ipf]['t1s'] = t1sdata
                else:
                    # entirely outside data range, do not plot
                    plotflag[ipf]['show'] = False
                    ndropped += 1

            if t2s > t2sdata:
                if t1s <= t2sdata:
                    # truncate to t2sdata
                    plotflag[ipf]['t2s'] = t2sdata
                else:
                    # entirely outside data range, do not plot
                    plotflag[ipf]['show'] = False
                    ndropped += 1

        if ndropped > 0:
            casalog.post('WARNING: Trimming dropped ' + str(ndropped)
                         + ' flags entirely')

    # make reason dictionary with mapping of colors and offsets (-0.3 to 0.3)
    readict = {}
    reakeys = []
    if nreas > ncolors:
        for i in range(nreas):
            reas = myreas[i]
            readict[reas] = {}
            if i < ncolors - 1:
                colr = colorlist[i]
                readict[reas]['color'] = colr
                readict[reas]['index'] = i
                offs = 0.3 - float(i) * 0.6 / float(ncolors - 1)
                readict[reas]['offset'] = offs
                reakeys.append(reas)
            else:
                colr = colorlist[ncolors - 1]
                readict[reas]['color'] = colr
                readict[reas]['index'] = ncolors - 1
                readict[reas]['offset'] = -0.3
        reakeys.append('Other')
        readict['Other'] = {}
        readict['Other']['color'] = colorlist[ncolors - 1]
        readict['Other']['index'] = ncolors - 1
        readict['Other']['offset'] = -0.3
    else:
        for i in range(nreas):
            reas = myreas[i]
            reakeys.append(reas)
            colr = colorlist[i]
            offs = 0.3 - float(i) * 0.6 / float(ncolors - 1)
            readict[reas] = {}
            readict[reas]['color'] = colr
            readict[reas]['index'] = i
            readict[reas]['offset'] = offs
    nlegend = len(reakeys)
    casalog.post('Will plot ' + str(nlegend) + ' reasons in legend')

    if plotname == '':
        pl.ion()
    else:
        pl.ioff()

    plotflagperant = defaultdict(list)
    for ipf, flag in lociteritems(plotflag):
        if not flag['show']:
            continue
        nflag = flag.copy()
        nflag['color'] = readict[flag['reason']]['color']
        nflag['offset'] = readict[flag['reason']]['offset']
        plotflagperant[flag['antenna']].append(nflag)

    nplotted = sum(len(x) for x in plotflagperant.values())
    casalog.post('Plotted %d flags' % nplotted)

    figs = []
    figsize = (8, 6)
    # maximum number of antennas per plot page (CAS-5187)
    antlimit = 28
    if len(myants) <= antlimit:
        figs.append(pl.figure(figsize=figsize))
        _plotants(figs[0], plotflagperant, myants, readict)
    else:
        # prefer DA on first page
        da = [x for x in myants if ('DA' in x) or ('CM' in x)]
        no_da = [x for x in myants if x not in da]
        # but limit to 28 antennas per figure
        if len(da) > antlimit:
            no_da.extend(da[antlimit:])
            da = da[:antlimit]

        if da:
            figs.append(pl.figure(figsize=figsize))
            _plotants(figs[-1], plotflagperant, da, readict)

        # stuff the rest on other figures
        while no_da:
            figs.append(pl.figure(figsize=figsize))
            _plotants(figs[-1], plotflagperant, no_da[:antlimit], readict)
            no_da = no_da[antlimit:]

    filenames = []
    if plotname == '':
        pl.draw()
    else:
        if len(figs) == 1:
            figs[0].savefig(plotname, dpi=150)
            filenames.append(plotname)
        else:
            fdirname = os.path.dirname(plotname)
            filename, ext = os.path.splitext(os.path.basename(plotname))
            for i, f in enumerate(figs):
                fn = '%s-%03d%s' % (os.path.join(fdirname, filename), i + 1, ext)
                filenames.append(fn)
                f.savefig(fn, dpi=150)
    return filenames


def _plotants(figure, plotflagperant, antlist, readict_inp):
    """
    plot flags of antennas

    Parameters
    ----------
    figure : matplotlib figure
    plotflagperant: dict
        dictionary of antennas containing list of flag plot dictionaries
    antlist: dict
        list of antenna names to plot
    readict:
        list of reasons for the legend
    """

    ax1 = figure.add_axes([.15, .1, .75, .85])
    nants = len(antlist)
    readict = dict()
    used_reasons = set()
    # These style params can be critical to produce meaningful (or not too
    # misleading) plots (CAS-13100)
    style_params = {'alpha': .7, 'marker': '.', 'markersize': 1, 'linewidth': 1}
    for antind, thisant in enumerate(antlist):
        for flag in plotflagperant[thisant]:
            thisoffset = flag['offset'] + antind + 1
            ax1.plot([flag['t1s'], flag['t2s']], [thisoffset] * 2,
                     color=flag['color'], **style_params)
            used_reasons.add(flag['reason'])

    # remove reasons that are not needed
    for k, v in readict_inp.items():
        if k in used_reasons:
            readict[k] = v

    myXlim = ax1.get_xlim()
    myXrange = myXlim[1] - myXlim[0]
    # Pad the time axis?
    PadTime = 0.050000000000000003
    if PadTime > 0:
        xPad = PadTime * myXrange
        x0 = myXlim[0] - xPad
        x1 = myXlim[1] + xPad
        ax1.set_xlim(x0, x1)
        myXrange = x1 - x0
    else:
        # casalog.post('  Rescaled x axis')
        x0 = myXlim[0]
        x1 = myXlim[1]

    legendFontSize = 12
    myYupper = nants + len(readict) + 1.5
    # place legend text
    x = x0 + 0.050000000000000003 * myXrange
    for i, reas in enumerate(readict.keys()):
        colr = readict[reas]['color']
        ax1.text(x, i + nants + 1, reas, color=colr, size=legendFontSize)
    ax1.set_ylim([0, myYupper])

    ax1.set_yticks(range(1, len(antlist) + 1))
    ax1.set_yticklabels(antlist)
    # casalog.post('  Relabled y axis')

    nxticks = 3
    ax1.set_xticks(np.linspace(myXlim[0], myXlim[1], nxticks))

    mytime = []
    myTimestr = []
    for itim in range(nxticks):
        time = myXlim[0] + (myXlim[1] - myXlim[0]) * float(itim) \
            / float(nxticks - 1)
        mytime.append(time)
        q1 = qalocal.quantity(time, 's')
        time1 = qalocal.time(q1, form='ymd', prec=9)[0]
        if itim > 0:
            time1s = time1[11:]
        else:
            time1s = time1
        myTimestr.append(time1s)

    ax1.set_xticklabels(myTimestr)

# def isModeValid(line):
#     '''Check if mode is valid based on a line
#        molinede --> line with strings
#         Returns True if mode is either one of the following:
#            '',manual,clip,quack,shadow,elevation      '''
# 
#     if line.__contains__('mode'):
#         if line.__contains__('manual') or line.__contains__('clip') \
#             or line.__contains__('quack') or line.__contains__('shadow'
#                 ) or line.__contains__('elevation'):
#             return True
#         else:
#             return False
# 
#     # No mode means manual
#     return True


#
#******************** CAL TABLE FUNCTIONS   **************************
#
def readCalCmds(caltable, msfile, flaglist, rows, reason, useapplied):
    '''Flag a cal table
    
    caltable    cal table name
    msfile      optional MS with flag cmds
    flagcmds    list with flag cmds or [] when msfile is given
    reason      select only flag cmds with this reason(s)
    useapplied  select APPLIED true or false
    '''
            
    myflagcmd = {}
    if msfile != '':   
        casalog.post('Reading flag cmds from FLAG_CMD table of MS')    
        # Read only the selected rows for action = apply and
        myflagcmd = readFromTable(msfile, myflagrows=rows, useapplied=useapplied, myreason=reason)
            
    elif flaglist != []:    
        # Parse the input file                    
        if isinstance(flaglist, list):
            casalog.post('Reading from input list')
            cmdlist = flaglist

            casalog.post('Input ' + str(len(cmdlist))
                         + ' lines from input list')
            # Make a FLAG_CMD compatible dictionary and select by reason
            myflagcmd = fh.parseDictionary(cmdlist, reason, False)

        elif isinstance(flaglist, str):

            casalog.post('Reading from input file '+flaglist)
            cmdlist = fh.readFile(flaglist)

            # Make a FLAG_CMD compatible dictionary and select by reason
            myflagcmd = fh.parseDictionary(cmdlist, reason, False)
            
        else:
            casalog.post('Unsupported inpfile type', 'ERROR')
                                            
    return myflagcmd

def applyCalCmds(aflocal, caltable, myflagcmd, tablerows, flagbackup, outfile):
    
    # Get the list of parameters
    cmdkeys = list(myflagcmd.keys())
    
    # Select the data
    selpars = {}    
    if len(cmdkeys) == 1:
        # Get all the selection parameters, but set correlation to ''
        cmd0 = myflagcmd[cmdkeys[0]]['command']
        selpars = fh.parseSelectionPars(cmd0)
        casalog.post('The selected subset of the MS will be: ')
        casalog.post('%s' % selpars)
    
    aflocal.selectdata(selpars)
        
    fh.parseAgents(aflocal, myflagcmd, [], True, True, '')
    
    # Initialize the Agents
    aflocal.init()
    
    # Backup the flags before running
    if flagbackup:
        fh.backupFlags(aflocal, msfile='', prename='flagcmd')
    
    # Run the tool
    aflocal.run(True, True)
    
    aflocal.done()

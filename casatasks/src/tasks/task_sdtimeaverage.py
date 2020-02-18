import os, re
import shutil
import string
import copy
import math
import time
import datetime
import contextlib

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from . import sdutil
    from casatasks import casalog
    from casatools import quanta, ms, table, mstransformer
    from .mstools import write_history
    from .parallel.parallel_data_helper import ParallelDataHelper
    from . import flaghelper as fh
    from .update_spw import update_spwchan
    from .callibrary import callibrary
    # CAS-12721
    from casatasks.private.sdutil import tbmanager   #, toolmanager, table_selector

else:
    import sdutil
    from taskinit import mttool as mstransformer
    from taskinit import mstool as ms
    from taskinit import tbtool as table
    from taskinit import qatool as quanta
    from mstools import write_history
    from parallel.parallel_data_helper import ParallelDataHelper
    import flaghelper as fh
    from update_spw import update_spwchan
    from callibrary import callibrary

@contextlib.contextmanager
def open_table(path, nomodify=True):
    tb = table( )
    tb.open(path, nomodify=nomodify)
    try:
        yield tb
    finally:
        tb.close()

def sdtimeaverage(
             infile,
             datacolumn,
             field,
             spw, 
             timerange, 
             scan,
             antenna,  
             timebin,
             timespan,
             outfile):
    #+
    #  defaut value (timebin=all) is to be handled.
    #-
    capTimebin = timebin.upper()
    if (capTimebin == 'ALL') or (capTimebin == ''):
        timebin =  calc_timebin(infile)+'s'
    #+
    # datacolumn alternative access 
    #  In case 'float_data' does not exists, attempt to use 'data'
    #-

    # know existence of data-column on specified MS. 
    ex_float_data, ex_data = check_column(infile) 
 
    # change datacolumn 'data' to 'float_data"
    if (datacolumn == 'float_data') :
        if (ex_float_data == False)and(ex_data == True):
            datacolumn = "data"
            msg = "No FLOAT_DATA column. DATA column will be used alternatively."
            casalog.post( msg, "INFO" )

    # change datacolumn 'float_data' to 'data' 
    if (datacolumn == 'data') :
        if (ex_float_data == True)and(ex_data == False):
            datacolumn = "float_data"
            msg = "No DATA column. FLOAT_DATA column will be used alternatively."
            casalog.post( msg, "INFO" )

    #+
    # Antanna ID (add extra &&& if needed) This is Single Dish specific 
    #-
    if (len(antenna) != 0) and (antenna.find('&') == -1):
        antenna = antenna + '&&&'

    #+
    # ALMA Warning
    #- 
    if ('scan' in timespan): 
        msg = "If the TELESCOPE NAME in specified MS contains ‘ALMA’, then mstransform automatically changes timescan = ’scan, state’ . You may receive unexpected result. "

        casalog.post( msg, "WARNING" )
   

    casalog.origin('sdtimeaverage')

    try:
        # Select Data and make Average. 
        st =do_mst(infile=infile, datacolumn=datacolumn,
                   field=field, spw=spw, timerange=timerange, scan=scan, antenna=antenna, 
                   timebin=timebin, timespan=timespan, outfile=outfile)
        # History 
        add_history(casalog=casalog, infile=infile, datacolumn=datacolumn,
                    field=field, spw=spw, timerange=timerange, scan=scan,
                    timebin=timebin, timespan=timespan, antenna=antenna, outfile=outfile)

    except Exception as e:
        casalog.post('Exception from task_sdtimeaverage : ' + str(e), "SEVERE", origin=origin)
        return False
    finally:
        pass

    return st

#  Getting Column Names
def check_column(msName):
    with open_table(msName) as tb:
        columnNames = tb.colnames()
        exist_float_data  = "FLOAT_DATA" in columnNames 
        exist_data        =  "DATA" in columnNames    
        return exist_float_data, exist_data

#  Calculation range time in input MS.  
def calc_timebin(msName):
    with open_table(msName) as tb:
        tm    = tb.getcol('TIME')
        iv    = tb.getcol('INTERVAL')

    interval = iv[0]  # interval

    time_first = min(tm)
    time_last =  max(tm)

    timebin = time_last - time_first ;
    timebin += 4.0 * interval   + 1.5  # ??? expand range. 

    return str(timebin)

# call mstransform by provided procedure #
def do_mst(infile, datacolumn, field, spw, timerange, scan, antenna, timebin, timespan, outfile):
    # followings are parameters of mstransform but not used by THIS.
    # just putting default values
  
    vis = infile             # needed for ParallelDataHelper
    outputvis = outfile      # needed for ParallelDataHelper
    tileshape = [0]

    intent = ""
    correlation = ""
    array = ""
    uvrange = ""
    observation = ""
    feed = ""

    realmodelcol = False
    usewtspectrum = False
    chanbin = 1
    mode = "channel"
    start = 0
    width = 1

    timeaverageAct = False
    maxuvwdistance = 0.0

    ddistart = -1
    reindex = True
    
    # Initialize the helper class
    pdh = ParallelDataHelper('sdtimeaverage', locals())
    pdh.bypassParallelProcessing(0)
    
    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post('%s'%instance,'ERROR')
        return False

    # Create a local copy of the MSTransform tool
    mtlocal = mstransformer()  # CASA6 changed.
    mslocal = ms( )            # CASA6 changed.
    
    try:
        # Gather all the parameters in a dictionary.
        config = {}
        
        # set config param.
        config = pdh.setupParameters(inputms=infile, outputms=outfile, field=field, 
                    spw=spw, array=array, scan=scan, antenna=antenna, correlation=correlation,
                    uvrange=uvrange, timerange=timerange, intent=intent, observation=str(observation),
                                     feed=feed, taql='')
        
        # ddistart will be used in the tool when re-indexing the spw table
        config['ddistart'] = ddistart
        
        # re-index parameter is used by the pipeline to not re-index any sub-table and the associated IDs
        config['reindex'] = reindex        
        
        config['datacolumn'] = datacolumn
        dc = datacolumn.upper()            
        # Make real a virtual MODEL column in the output MS
        if "MODEL" in dc or dc == 'ALL':
            config['realmodelcol'] = realmodelcol

        config['usewtspectrum'] = usewtspectrum
        config['tileshape'] = tileshape

        # Only parse timeaverage parameters when timebin > 0s
        if is_CASA6:
                qa = quanta( ) # CASA6 needed

        tbin = qa.convert(qa.quantity(timebin), 's')['value']
        if tbin < 0:
            raise Exception("Parameter timebin must be > '0s' to do time averaging")
  
        # set config for Averaging
        timeaverageAct = (tbin > 0)
        if timeaverageAct:
            casalog.post('Parse time averaging parameters')
            config['timeaverage'] = True
            config['timebin'] = timebin
            config['timespan'] = timespan
            config['maxuvwdistance'] = maxuvwdistance
            
        # Configure the tool and all the parameters
        casalog.post('%s'%config, 'DEBUG')
        mtlocal.config(config)
        
        # Open the MS, select the data and configure the output
        mtlocal.open()
        
        # Run the tool
        casalog.post('Apply the transformations')
        mtlocal.run()        

        mtlocal.done()
                    
    except Exception as instance:
        mtlocal.done()
        casalog.post('%s'%instance,'ERROR')
        return False

    #+
    # CAS-12721:
    # Note: Following section were written concerning with CAS-7751 or other(s)
    #       Program logic is copied and used without change. 
    #-

    # Update the FLAG_CMD sub-table to reflect any spw/channels selection
    # If the spw selection is by name or FLAG_CMD contains spw with names, skip the updating    
    
    if (spw != '') and (spw != '*'):
        isopen = False

        try:
            mytb = table( )
            mytb.open(outfile + '/FLAG_CMD', nomodify=False)
            isopen = True
            nflgcmds = mytb.nrows()
            
            if nflgcmds > 0:
                updateFlagCmd = False

                # If spw selection is by name in FLAG_CMD, do not update, CAS-7751
                mycmd = mytb.getcell('COMMAND', 0)
                cmdlist = mycmd.split()
                for cmd in cmdlist:
                    # Match only spw indices, not names
                    if cmd.__contains__('spw'):
                        cmd = cmd.strip("spw=")
                        spwstr = re.search('^[^a-zA-Z]+$', cmd)
                        if spwstr != None and spwstr.string.__len__() > 0:
                            updateFlagCmd = True
                            break                
                

                if updateFlagCmd:
                    mademod = False
                    cmds = mytb.getcol('COMMAND')
                    widths = {}
                    #print("width =", width)
                    if hasattr(chanbin, 'has_key'):
                        widths = chanbin
                    else:
                        if hasattr(chanbin, '__iter__') and len(chanbin) > 1:
                            for i in range(len(chanbin)):
                                widths[i] = chanbin[i]
                        elif chanbin != 1:
    #                        print('using ms.msseltoindex + a scalar width')
                            numspw = len(mslocal.msseltoindex(vis=infile,
                                                         spw='*')['spw'])
                            if hasattr(chanbin, '__iter__'):
                                w = chanbin[0]
                            else:
                                w = chanbin
                            for i in range(numspw):
                                widths[i] = w
    #                print('widths =', widths)
                    for rownum in range(nflgcmds):
                        # Matches a bare number or a string quoted any way.
                        spwmatch = re.search(r'spw\s*=\s*(\S+)', cmds[rownum])
                        if spwmatch:
                            sch1 = spwmatch.groups()[0]
                            sch1 = re.sub(r"[\'\"]", '', sch1)  # Dequote
                            # Provide a default in case the split selection excludes
                            # cmds[rownum].  update_spwchan() will throw an exception
                            # in that case.
                            cmd = ''
                            try:
                                #print('sch1 =', sch1)
                                sch2 = update_spwchan(infile, spw, sch1, truncate=True,
                                                      widths=widths)
                                #print('sch2 =', sch2)
                                ##print('spwmatch.group() =', spwmatch.group())
                                if sch2:
                                    repl = ''
                                    if sch2 != '*':
                                        repl = "spw='" + sch2 + "'"
                                    cmd = cmds[rownum].replace(spwmatch.group(), repl)
                            #except: # cmd[rownum] no longer applies.
                            except Exception as e:
                                casalog.post(
                                    "Error %s updating row %d of FLAG_CMD" % (e,
                                                                              rownum),
                                             'WARN')
                                casalog.post('sch1 = ' + sch1, 'DEBUG1')
                                casalog.post('cmd = ' + cmd, 'DEBUG1')
                            if cmd != cmds[rownum]:
                                mademod = True
                                cmds[rownum] = cmd
                    if mademod:
                        casalog.post('Updating FLAG_CMD', 'INFO')
                        mytb.putcol('COMMAND', cmds)

                else:
                    casalog.post('FLAG_CMD table contains spw selection by name. Will not update it!','DEBUG')
                
            mytb.close()
            
        except Exception as instance:
            if isopen:
                mytb.close()
            mslocal = None
            mytb = None
            casalog.post("*** Error \'%s\' updating FLAG_CMD" % (instance),
                         'SEVERE')
            return False

    mytb = None
    mslocal = None

    return True

# revised in CASA6 
def add_history(casalog, infile, datacolumn, field, spw, timerange, scan, timebin, timespan, antenna, outfile):
    mslocal = ms( )
    # Write history to output MS, not the input ms.
    try:
        code_object = sdtimeaverage.__code__                             #CASA6
        param_names = code_object.co_varnames[:code_object.co_argcount]  #CASA6
        local_vals = locals()                                            #CASA6
        param_vals = [local_vals.get(p, None) for p in param_names]      #CASA6
        write_history(mslocal, outfile, 'sdtimeaverage', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                     'WARN')
        return False

    mslocal = None

    return True
#END
   

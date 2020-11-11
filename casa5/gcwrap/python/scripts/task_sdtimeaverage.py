import re
import numpy

import sdutil
from taskinit import mttool, mstool, tbtool, casalog, qa
from mstools import write_history
from parallel.parallel_data_helper import ParallelDataHelper
from update_spw import update_spwchan


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
    #  When 'all'(default) or '' is specified, make timebin to cover all the data.
    if timebin.upper() in ['ALL', '']:
        timebin = set_timebin_all() + 's'

    # Switch Alternative Column if needed.
    active_datacolumn = use_alternative_column(infile, datacolumn)

    # Antenna ID (add extra &&& if needed) This is Single Dish specific.
    if (len(antenna) != 0) and (not('&' in antenna)):
        antenna = antenna + '&&&'

    # 'scan,state' Warning
    #   (unexpected result warning)
    if ('scan' in timespan) and ('state' in timespan):
        #  WARN msg, to explain NRO specific issue.
        msg = '\n'.join(["Explicitly both 'scan' and 'state' were specified in TIMESPAN.",
                         "   If 'state' distinguishes OBSERVE_TARGET#ON_SOURCE / OBSERVE_TARGET#OFF_SOURCE,",
                         "   these two states are mixed, and unexpectedly averaged results might be generated.",
                         "(Suggestion) Please specify timespan = 'scan'",
                         "             to separate OBSERVE_TARGET#ON_SOURCE and OBSERVE_TARGET#OFF_SOURCE."])
        casalog.post(msg, 'WARN')

    # (Note) When timespan='' is specified.
    #  timespan will be directly posted to mstransform.

    # Convert to check timebin
    tbin = qa.convert(qa.quantity(timebin), 's')['value']

    # Time average, once Enable.
    do_timeaverage = True

    # Check timebin
    if tbin < 0:     # Error, raise Exception.
        raise Exception(
            "Parameter timebin must be >= '0s' to do time averaging")
    elif tbin == 0:  # No averaging, when tbin == 0
        msg = 'Parameter timebin equals zero. No averaging will be performed.'
        casalog.post(msg, 'WARN')
        do_timeaverage = False

    origin = 'sdtimeaverage'
    casalog.origin(origin)

    try:
        # Select Data and make Average.
        st = do_mst(
            infile=infile,
            datacolumn=active_datacolumn,
            field=field,
            spw=spw,
            timerange=timerange,
            scan=scan,
            antenna=antenna,
            timebin=timebin,
            timespan=timespan,
            outfile=outfile,
            do_timeaverage=do_timeaverage)

        # History
        add_history(
            casalog=casalog,
            infile=infile,
            datacolumn=active_datacolumn,
            field=field,
            spw=spw,
            timerange=timerange,
            scan=scan,
            timebin=timebin,
            timespan=timespan,
            antenna=antenna,
            outfile=outfile)

    except Exception as e:
        casalog.post(
            'Exception from task_sdtimeaverage : ' +
            str(e),
            'SEVERE',
            origin=origin)
        return False

    return st


def use_alternative_column(infile, datacolumn):
    """
     Alternatively use datacolumn if specified column does not exist.
       In case 'float_data' does not exist, attempt to use 'data'
       and vice versa. (For user's convenience)
    """
    #  know existence of data-column on specified MS.
    ex_float_data, ex_data = check_column(infile)

    # alter datacolumn if available
    if (datacolumn == 'float_data'):  # Change 'float_data' to 'data'
        if (not ex_float_data) and (ex_data):
            datacolumn = 'data'
            msg = 'No FLOAT_DATA column. DATA column will be used alternatively.'
            casalog.post(msg, 'INFO')
    elif (datacolumn == 'data'):      # Change 'data' to 'float_data'
        if (ex_float_data) and (not ex_data):
            datacolumn = 'float_data'
            msg = 'No DATA column. FLOAT_DATA column will be used alternatively.'
            casalog.post(msg, 'INFO')

    return datacolumn


def check_column(msname):
    """ Check specified column if exists. """
    with sdutil.tbmanager(msname) as tb:
        columnNames = tb.colnames()
        exist_float_data = 'FLOAT_DATA' in columnNames
        exist_data = 'DATA' in columnNames
        return exist_float_data, exist_data


def set_timebin_all():
    """ Calculation range time in input MS. """
    # assign very large value to timebin
    #  to cover 'all'.
    timebin = numpy.finfo(float).max
    return str(timebin)


def do_mst(
        infile,
        datacolumn,
        field,
        spw,
        timerange,
        scan,
        antenna,
        timebin,
        timespan,
        outfile,
        do_timeaverage):
    """
      call mstransform by provided procedure.
        Followings are parameters of mstransform, but not used by THIS,
        just only putting default values
    """
    vis = infile             # needed for ParallelDataHelper
    outputvis = outfile      # needed for ParallelDataHelper
    tileshape = [0]

    intent = ''
    correlation = ''
    array = ''
    uvrange = ''
    observation = ''
    feed = ''

    realmodelcol = False
    usewtspectrum = False
    chanbin = 1
    mode = 'channel'
    start = 0
    width = 1

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
        casalog.post('%s' % instance, 'ERROR')
        return False

    # Create a local copy of the MSTransform tool
    mtlocal = mttool()  # CASA5 feature
    mslocal = mstool()  # CASA5 feature

    try:
        # Gather all the parameters in a dictionary.
        config = {}

        # set config param.
        config = pdh.setupParameters(
            inputms=infile,
            outputms=outfile,
            field=field,
            spw=spw,
            array=array,
            scan=scan,
            antenna=antenna,
            correlation=correlation,
            uvrange=uvrange,
            timerange=timerange,
            intent=intent,
            observation=str(observation),
            feed=feed,
            taql='')

        # ddistart will be used in the tool when re-indexing the spw table
        config['ddistart'] = ddistart

        # re-index parameter is used by the pipeline to not re-index any
        # sub-table and the associated IDs
        config['reindex'] = reindex

        config['datacolumn'] = datacolumn
        dc = datacolumn.upper()
        # Make real a virtual MODEL column in the output MS
        if 'MODEL' in dc or dc == 'ALL':
            config['realmodelcol'] = realmodelcol

        config['usewtspectrum'] = usewtspectrum
        config['tileshape'] = tileshape

        # set config for Averaging
        if do_timeaverage:
            casalog.post('Parse time averaging parameters')
            config['timeaverage'] = True
            config['timebin'] = timebin
            config['timespan'] = timespan
            config['maxuvwdistance'] = maxuvwdistance

        # Configure the tool and all the parameters
        casalog.post('%s' % config, 'DEBUG')
        mtlocal.config(config)

        # Open the MS, select the data and configure the output
        mtlocal.open()

        # Run the tool
        casalog.post('Apply the transformations')
        mtlocal.run()

        mtlocal.done()

    except Exception as instance:
        mtlocal.done()
        casalog.post('%s' % instance, 'ERROR')
        return False

    """
      CAS-12721:
      Note: Following section were written concerning with CAS-7751 or others.
            Program logic is copied and used without change.
    """
    # Update the FLAG_CMD sub-table to reflect any spw/channels selection
    # If the spw selection is by name or FLAG_CMD contains spw with names,
    # skip the updating

    if (spw != '') and (spw != '*'):
        isopen = False

        try:
            mytb = tbtool()
            mytb.open(outfile + '/FLAG_CMD', nomodify=False)
            isopen = True
            nflgcmds = mytb.nrows()

            if nflgcmds > 0:
                update_flag_cmd = False

                # If spw selection is by name in FLAG_CMD, do not update, CAS-7751
                mycmd = mytb.getcell('COMMAND', 0)
                cmdlist = mycmd.split()
                for cmd in cmdlist:
                    # Match only spw indices, not names
                    if cmd.__contains__('spw'):
                        cmd = cmd.strip('spw=')
                        spwstr = re.search('^[^a-zA-Z]+$', cmd)
                        if spwstr is not None and spwstr.string.__len__() > 0:
                            update_flag_cmd = True
                            break

                if update_flag_cmd:
                    mademod = False
                    cmds = mytb.getcol('COMMAND')
                    widths = {}
                    # casalog.post "width =", width
                    if hasattr(chanbin, 'has_key'):
                        widths = chanbin
                    else:
                        if hasattr(chanbin, '__iter__') and len(chanbin) > 1:
                            for i in xrange(len(chanbin)):
                                widths[i] = chanbin[i]
                        elif chanbin != 1:
                            # casalog.post('using ms.msseltoindex + a scalar width')
                            numspw = len(mslocal.msseltoindex(vis=infile,
                                                              spw='*')['spw'])
                            if hasattr(chanbin, '__iter__'):
                                w = chanbin[0]
                            else:
                                w = chanbin
                            for i in xrange(numspw):
                                widths[i] = w
                    # casalog.post('widths =', widths)
                    for rownum in xrange(nflgcmds):
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
                                # casalog.post('sch1 =', sch1)
                                sch2 = update_spwchan(
                                    infile, spw, sch1, truncate=True, widths=widths)
                                # casalog.post('sch2 =', sch2)
                                # casalog.post('spwmatch.group() =', spwmatch.group())
                                if sch2:
                                    repl = ''
                                    if sch2 != '*':
                                        repl = "spw='" + sch2 + "'"
                                    cmd = cmds[rownum].replace(
                                        spwmatch.group(), repl)
                            # except: # cmd[rownum] no longer applies.
                            except Exception as e:
                                casalog.post(
                                    "Error %s updating row %d of FLAG_CMD" %
                                    (e, rownum), 'WARN')
                                casalog.post('sch1 = ' + sch1, 'DEBUG1')
                                casalog.post('cmd = ' + cmd, 'DEBUG1')
                            if cmd != cmds[rownum]:
                                mademod = True
                                cmds[rownum] = cmd
                    if mademod:
                        casalog.post('Updating FLAG_CMD', 'INFO')
                        mytb.putcol('COMMAND', cmds)

                else:
                    casalog.post(
                        'FLAG_CMD table contains spw selection by name. Will not update it!', 'DEBUG')

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


# CASA5 / CASA6 different (see CASA6 code)
def add_history(
        casalog,
        infile,
        datacolumn,
        field,
        spw,
        timerange,
        scan,
        timebin,
        timespan,
        antenna,
        outfile):
    mslocal = mstool()
    # Write history to output MS, not the input ms.
    try:
        param_names = sdtimeaverage.func_code.co_varnames[:sdtimeaverage.func_code.co_argcount]  # CASA5
        param_vals = [eval(p) for p in param_names]                                              # CASA5
        write_history(mslocal, outfile, 'sdtimeaverage', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                     'WARN')
        return False

    mslocal = None

    return True
# END

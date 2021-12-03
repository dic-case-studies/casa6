import inspect
from types import CodeType

import numpy
from casatasks.private.casa_transition import is_CASA6

if is_CASA6:
    from casatasks import casalog
    from casatools import quanta

    from . import sdutil
else:
    import sdutil
    from taskinit import casalog
    from taskinit import qatool as quanta

qa = quanta()


@sdutil.sdtask_decorator
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
        msg = '\n'.join(["Explicitly both 'scan' and 'state' were specified in 'timespan'.",
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
        raise ValueError(
            "Parameter timebin must be >= '0s' to do time averaging")
    elif tbin == 0:  # No averaging, when tbin == 0
        msg = 'Parameter timebin equals zero. No averaging will be performed.'
        casalog.post(msg, 'WARN')
        do_timeaverage = False

    # extra parameter for do_mst
    ext_config = {'do_timeaverage': do_timeaverage}

    caller: CodeType = inspect.currentframe().f_code

    # Select Data and make Average.
    sdutil.do_mst(
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
        intent='',
        caller=caller,
        ext_config=ext_config)

    # History
    sdutil.add_history(caller, casalog, outfile)


def use_alternative_column(infile, datacolumn):
    """
     Alternatively use datacolumn if the specified column does not exist.
       In case 'float_data' does not exist, sdtimeaverage attempt to use 'data'
       and vice versa. (For user's convenience)
    """
    #  obtain the existence of data-column on specified MS.
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
    """ Check the specified column if it exists. """
    with sdutil.table_manager(msname) as tb:
        columnNames = tb.colnames()
        exist_float_data = 'FLOAT_DATA' in columnNames
        exist_data = 'DATA' in columnNames
        return exist_float_data, exist_data


def set_timebin_all():
    """
      Synthesize timebin
        assign very large value to cover 'all'.
    """
    timebin = numpy.finfo(float).max
    return str(timebin)

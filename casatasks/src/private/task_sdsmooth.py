import os

# get is_python3 and is_CASA6
from casatasks.private.casa_transition import *

if is_CASA6:
    from casatasks import casalog
    from casatools import ms as mstool
    from casatools import singledishms, table

    from . import sdutil
    from .mstools import write_history

    ms = mstool()
    sdms = singledishms()
    tb = table()
else:
    import sdutil
    from mstools import write_history
    from taskinit import casalog, gentools
    ms, sdms, tb = gentools(['ms', 'sdms', 'tb'])


@sdutil.callable_sdtask_decorator
def sdsmooth(infile=None, datacolumn=None, antenna=None,
             field=None, spw=None, timerange=None, scan=None,
             pol=None, intent=None, reindex=None,
             kernel=None, kwidth=None,
             outfile=None, overwrite=None):

    try:
        if len(outfile) == 0:
            errmsg = 'outfile is empty.'
            raise ValueError(errmsg)
        
        if (os.path.exists(outfile)) and (not overwrite):
            errmsg = outfile+' exists.'
            raise ValueError(errmsg)

        sdms.open(infile)
        sdms.set_selection(spw=spw, field=field, 
                           antenna=antenna,
                           timerange=timerange, scan=scan,
                           polarization=pol, intent=intent,
                           reindex=reindex)
        sdms.smooth(type=kernel, width=kwidth, datacolumn=datacolumn, outfile=outfile)
        
        # Write to HISTORY of outfile MS
        param_names = sdsmooth.__code__.co_varnames[:sdsmooth.__code__.co_argcount]
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        write_history(ms, outfile, 'sdsmooth', param_names,
                      param_vals, casalog)

    finally:
        sdms.close()

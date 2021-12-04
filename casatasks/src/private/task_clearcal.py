from __future__ import absolute_import
import os

from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import calibrater, table, ms
    from casatasks import casalog
    from .mstools import write_history
    from .parallel.parallel_data_helper import ParallelDataHelper
    from .parallel.parallel_task_helper import ParallelTaskHelper
else:
    from taskinit import cbtool as calibrater
    from taskinit import tbtool as table
    from taskinit import mstool as ms
    from taskinit import casalog
    from mstools import write_history
    from parallel.parallel_data_helper import ParallelDataHelper
    from parallel.parallel_task_helper import ParallelTaskHelper

def clearcal(
    vis=None,
    field=None,
    spw=None,
    intent=None,
    addmodel=None,
    ):

    casalog.origin('clearcal')

    # Do the trivial parallelization
    if ParallelDataHelper.isMMSAndNotServer(vis):
        helper = ParallelTaskHelper('clearcal', locals())
        helper.go()
        return

    # Local versions of the tools
    tblocal = table()
    cblocal = calibrater()
    mslocal = ms()

    # we will initialize scr cols only if we don't create them
    doinit = False

    if (type(vis) == str) & os.path.exists(vis):
        tblocal.open(vis)
        doinit = tblocal.colnames().count('CORRECTED_DATA') > 0
        tblocal.close()

        # We ignore selection if creating the scratch columns
        if not doinit:
            casalog.post('Need to create scratch columns; ignoring selection.')

        cblocal.setvi(old=True,quiet=False);  # Old VI for now
        cblocal.open(vis, addmodel=addmodel)
    else:
        raise Exception('Visibility data set not found - please verify the name')

    # If necessary (scr col not just created), initialize scr cols
    if doinit:
        cblocal.selectvis(field=field, spw=spw, intent=intent)
        cblocal.initcalset(1)
    cblocal.close()

    # Write history to the MS
    if is_python3:
        vars = locals( )
        param_names = clearcal.__code__.co_varnames[:clearcal.__code__.co_argcount]
        param_vals = [vars[p] for p in param_names]
    else:
        param_names = clearcal.__code__.co_varnames[:clearcal.__code__.co_argcount]
        param_vals = [eval(p) for p in param_names]

    casalog.post('Updating the history in the output', 'DEBUG1')
    write_history(mslocal, vis, 'clearcal', param_names,
                  param_vals, casalog)


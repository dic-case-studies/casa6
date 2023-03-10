from __future__ import absolute_import
import os

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import calibrater, ms
    from casatasks import casalog

    from .mstools import write_history
    from .parallel.parallel_data_helper import ParallelDataHelper
    from .parallel.parallel_task_helper import ParallelTaskHelper
else:
    from taskinit import *
    from mstools import write_history
    from parallel.parallel_data_helper import ParallelDataHelper
    from parallel.parallel_task_helper import ParallelTaskHelper

    calibrater = cbtool
    ms = mstool

def initweights(vis=None,wtmode=None,tsystable=None,gainfield=None,interp=None,spwmap=None,dowtsp=None):

    casalog.origin('initweights')

    # Do the trivial parallelization
    if ParallelTaskHelper.isMPIEnabled() and ParallelDataHelper.isMMSAndNotServer(vis):
        tsystable = ParallelTaskHelper.findAbsPath(tsystable)
        helper = ParallelTaskHelper('initweights', locals())
        helper.go()
        # Write history to MS.
        try:
            param_names = initweights.__code__.co_varnames[:initweights.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            casalog.post('Updating the history in the output', 'DEBUG1')
            write_history(ms, vis, 'initweights', param_names,
                          param_vals, casalog)
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                         'WARN')
        
        return


    #Python script
    mycb=calibrater()

    # only if vis exists...
    if ((type(vis)==str) & (os.path.exists(vis))):
        if wtmode.upper().find("TSYS") > -1:
            if not os.path.exists(tsystable):
                raise Exception('Tsys calibration table %s not found' % tsystable)
            if len(spwmap)==0:
                spwmap=[-1]
            if interp=="":
                interp="linear"
        # ... and we are asked to do something...
        # open without adding anything!
        mycb.open(vis,compress=False,addcorr=False,addmodel=False)
        mycb.initweights(wtmode=wtmode,dowtsp=dowtsp,tsystable=tsystable,gainfield=gainfield,interp=interp,spwmap=spwmap)
        mycb.close()
    else:
        raise ValueError('Visibility data set not found - please verify the name')

    # Write history to MS.
    # When running in parallel, history will be written in the parallel section above
    # normal MSs should write the history here
    if ParallelTaskHelper.isMPIClient():
        myms=ms()
        try:
            param_names = initweights.__code__.co_varnames[:initweights.__code__.co_argcount]
            if is_python3:
                vars = locals()
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            casalog.post('Updating the history in the output', 'DEBUG1')
            write_history(myms, vis, 'initweights', param_names,
                          param_vals, casalog)
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % instance,'WARN')


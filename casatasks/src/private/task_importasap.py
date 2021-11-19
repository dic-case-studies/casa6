import os
import re

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *

if is_CASA6:
    from casatasks import casalog
    from casatools import agentflagger, calibrater, ms, singledishms
    from . import sdutil

    from .mstools import write_history

    mysdms = singledishms()
    mycb = calibrater()
    myms = ms()
else:
    from casac import casac
    from mstools import write_history
    from taskinit import *
    import sdutil

    agentflagger = casac.agentflagger

    mysdms, mycb, myms = gentools(['sdms', 'cb', 'ms'])

@sdutil.sdtask_decorator
def importasap(infile=None, outputvis=None, flagbackup=None, overwrite=None, parallel=None):
    """
    """

    try:
        if infile is None or outputvis is None:
            raise RuntimeError('Error: infile and outputvis must be specified.')

        # default value
        if flagbackup is None:
            flagbackup = True

        if overwrite is None:
            overwrite = False

        # basic check
        if os.path.exists(outputvis) and not overwrite:
            raise RuntimeError('%s exists.'%(outputvis))
        
        if not _is_scantable(infile):
            raise RuntimeError('%s is not a valid Scantable.'%(infile))

        # import
        status = mysdms.importasap(infile, outputvis, parallel)

        if status:
            # flagversions file must be deleted 
            flagversions = outputvis.rstrip('/') + '.flagversions'
            if os.path.exists(flagversions):
                os.system('rm -rf %s'%(flagversions))

            # initialize weights using cb tool
            mycb.open(outputvis, compress=False, addcorr=False, addmodel=False)
            mycb.initweights(wtmode='nyq')

            # create flagbackup file if user requests it
            if flagbackup:
                aflocal = agentflagger()
                aflocal.open(outputvis)
                aflocal.saveflagversion('Original',
                                        comment='Original flags at import into CASA using importasap',
                                        merge='save')
                aflocal.done()
        else:
            raise RuntimeError('Failure in importasap')

        # Write history to output MS
        param_names = importasap.__code__.co_varnames[:importasap.__code__.co_argcount]
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        write_history(myms, outputvis, 'importasap', param_names,
                      param_vals, casalog)

    finally:
        mycb.close()


def _is_scantable(filename):
    """
    Check if given data is Scantable or not
    """
    ret = False
    if os.path.isdir(filename) and os.path.exists(filename+'/table.info') \
        and os.path.exists(filename+'/table.dat'):
        with open(filename+'/table.info') as f:
            l=f.readline()
            f.close()
        match_pattern = '^Type = (Scantable)? *$'
        ret = re.match(match_pattern, l) is not None
    return ret

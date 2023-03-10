import datetime
import os
import sys
import shutil

from casatasks import casalog
from casatools import calibrater, ms, singledishms

from . import sdutil
from .mstools import write_history

mysdms = singledishms()
mycb = calibrater()
myms = ms()


@sdutil.sdtask_decorator
def importnro(infile=None, outputvis=None, overwrite=None, parallel=None):
    """
    """
    status = True

    try:
        outputvis_temp = outputvis + '-backup-' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')

        if os.path.exists(outputvis):
            if overwrite:
                os.rename(outputvis, outputvis_temp)
            else:
                raise RuntimeError('%s exists.' % (outputvis))

        if not _is_nostar(infile):
            raise RuntimeError('%s is not a valid NOSTAR data.' % (infile))

        status = mysdms.importnro(infile, outputvis, parallel)

        if status:
            # initialize weights using cb tool
            mycb.open(outputvis, compress=False, addcorr=False, addmodel=False)
            mycb.initweights(wtmode='nyq')
            if os.path.exists(outputvis_temp):
                shutil.rmtree(outputvis_temp)
        else:
            if os.path.exists(outputvis):
                shutil.rmtree(outputvis)
            if os.path.exists(outputvis_temp):
                os.rename(outputvis_temp, outputvis)
            raise RuntimeError('import failed.')

        # Write parameters to HISTORY table of MS
        param_names = importnro.__code__.co_varnames[:importnro.__code__.co_argcount]
        vars = locals()
        param_vals = [vars[p] for p in param_names]
        write_history(myms, outputvis, 'importnro', param_names,
                      param_vals, casalog)

    finally:
        if status:
            mycb.close()


def _is_nostar(filename):
    """Check if given data is NOSTAR or not."""
    ret = False
    if os.path.getsize(filename) >= 15136:  # size of observation header
        with open(filename, 'rb') as f:
            if f.read(8).replace(b'\x00', b'') == b'RW':
                ret = (f.read(15136 - 8 + 4)[-4:].replace(b'\x00', b'').decode( ) == 'LS')
            f.close()
    return ret

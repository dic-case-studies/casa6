from __future__ import absolute_import
import datetime
import os
import re
import shutil

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import singledishms, calibrater, ms
    from casatools.platform import bytes2str
    from casatasks import casalog
    from .mstools import write_history

    mysdms = singledishms( )
    mycb = calibrater( )
    myms = ms( )
else:
    from taskinit import *
    from mstools import write_history

    mysdms, mycb, myms = gentools(['sdms', 'cb', 'ms'])

def importnro(infile=None, outputvis=None, overwrite=None, parallel=None):
    """
    """
    casalog.origin('importnro')
    status = True

    try:
        outputvis_temp = outputvis + '-backup-' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        
        if os.path.exists(outputvis):
            if overwrite:
                os.rename(outputvis, outputvis_temp)
            else:
                raise RuntimeError('%s exists.'%(outputvis))
        
        if not _is_nostar(infile):
            raise RuntimeError('%s is not a valid NOSTAR data.'%(infile))

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
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        write_history(myms, outputvis, 'importnro', param_names, 
                          param_vals, casalog) 
        

    finally:
        if status: 
            mycb.close()


def _is_nostar(filename):
    """
    Check if given data is NOSTAR or not
    """
    ret = False
    if os.path.getsize(filename) >= 15136: # size of observation header
        with open(filename, 'rb') as f:
            if is_CASA6:
                if bytes2str(f.read(8).replace(b'\x00', b'')) == 'RW':
                    ret = (bytes2str(f.read(15136-8+4)[-4:].replace(b'\x00', b'')) == 'LS')
            else:
               if f.read(8).replace('\x00', '') == 'RW':
                    ret = f.read(15136-8+4)[-4:].replace('\x00', '') == 'LS'
            f.close()

    return ret

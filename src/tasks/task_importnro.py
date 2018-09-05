import datetime
import os
import re
import shutil

from CASAtools import singledishms, calibrater
from CASAtools.platform import bytes2str
from CASAtasks import casalog

mysdms = singledishms( )
mycb = calibrater( )

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
        
        return status
    except Exception as instance:
        casalog.post('*** Error *** %s' % instance, 'SEVERE')
        raise
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
            if bytes2str(f.read(8).replace(b'\x00', b'')) == 'RW':
                ret = (bytes2str(f.read(15136-8+4)[-4:].replace(b'\x00', b'')) == 'LS')
            f.close()

    return ret

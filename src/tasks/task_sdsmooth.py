import numpy
import os

from CASAtools import singledishms, table
from CASAtools import ms as mstool
from CASAtasks import casalog
from . import sdutil

ms = mstool( )
sdms = singledishms( )
tb = table( )

def sdsmooth(infile=None, datacolumn=None, antenna=None, 
              field=None, spw=None, timerange=None, scan=None, 
              pol=None, intent=None, reindex=None,
              kernel=None, kwidth=None,
              outfile=None, overwrite=None):

    casalog.origin('sdsmooth')

    try:
        if len(outfile) == 0:
            errmsg = 'outfile is empty.'
            raise_exception(errmsg)
        
        if (os.path.exists(outfile)) and (not overwrite):
            errmsg = outfile+' exists.'
            raise_exception(errmsg)

        sdms.open(infile)
        sdms.set_selection(spw=spw, field=field, 
                           antenna=antenna,
                           timerange=timerange, scan=scan,
                           polarization=pol, intent=intent,
                           reindex=reindex)
        sdms.smooth(type=kernel, width=kwidth, datacolumn=datacolumn, outfile=outfile)
    except Exception:
        raise
    finally:
        sdms.close()

def raise_exception(errmsg):
    casalog.post(errmsg, priority='SEVERE')
    raise Exception(errmsg)

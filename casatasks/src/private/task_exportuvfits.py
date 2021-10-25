import os
from casatools import ms

from casatasks import casalog

def exportuvfits( vis, fitsfile, datacolumn, field, spw, antenna, time,
                  writesyscal, multisource, combinespw,
                  writestation, padwithflags, overwrite ):

    casalog.origin('exportuvfits')

    try:
        myms = ms( )
        if ((type(vis)==str) & (os.path.exists(vis))):
            myms.open( vis, lock=True )
        else:
            raise ValueError('Visibility data set not found - please verify the name')
        writesyscal=False #until ms syscal table defined
        res = myms.tofits( fitsfile=fitsfile, column=datacolumn, field=field, spw=spw,
                           baseline=antenna, time=time, writesyscal=writesyscal,
                           multisource=multisource, combinespw=combinespw,
                           writestation=writestation, padwithflags=padwithflags,
                           overwrite=overwrite )
        if not res:
            raise Exception("exportuvfits failed")
    finally:
        if myms:
            myms.done()


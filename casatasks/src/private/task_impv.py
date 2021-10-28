from __future__ import absolute_import
import sys

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image
    from casatasks import casalog
    from .ialib import write_image_history
else:
    from taskinit import *
    from ialib import write_image_history
    image = iatool

def impv(
    imagename, outfile, mode, start, end, center, length, pa, width,
    unit, overwrite, region, chans, stokes, mask, stretch
):
    casalog.origin('impv')
    try:
        if len(outfile) == 0:
            raise Exception("outfile must be specified.")
        mymode = mode.lower()
        if mymode.startswith('c'):
            if len(start) == 0 or len(end) == 0:
                raise ValueError("When mode='coords', start and end must both be specified.")
            center = ""
            length = ""
            pa = ""
        elif mymode.startswith('l'):
            if (
                len(center) == 0 
                or (
                    not isinstance(length, (int, float))
                    and len(length) == 0
                )
                or len(pa) == 0
            ):
                raise ValueError("When mode='length', center, length, and pa must all be specified.")
            start = ""
            end = ""
        else:
            raise ValueError("Unsupported value for mode.")
        myia = image()
        outia = None
        myia.dohistory(False)
        if (not myia.open(imagename)):
            raise RuntimeError("Cannot create image analysis tool using %s" % imagename)
        outia = myia.pv(
            outfile=outfile, start=start, end=end, center=center,
            length=length, pa=pa, width=width, unit=unit,
            overwrite=overwrite, region=region, chans=chans,
            stokes=stokes, mask=mask, stretch=stretch, wantreturn=True
        )

        try:
            param_names = impv.__code__.co_varnames[:impv.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]   
            write_image_history(
                outia, sys._getframe().f_code.co_name,
                param_names, param_vals, casalog
            )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')

    finally:
        if (myia):
            myia.done()
        if (outia):
            outia.done()

from __future__ import absolute_import
import sys

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image, coordsys, regionmanager
    from casatasks import casalog

    from .ialib import write_image_history
else:
    from taskinit import *

    from ialib import write_image_history

    image = iatool
    regionmanager = rgtool
    coordsys = cstool

def imdev(
    imagename, outfile, region, box, chans,
    stokes, mask, overwrite, stretch,
    grid, anchor, xlength, ylength, interp, stattype, statalg,
    zscore, maxiter
):
    _myia = image()
    _myrg = regionmanager()
    _mycs = coordsys()
    try:
        casalog.origin('imdev')
        _myia.open(imagename)
        _mycs = _myia.coordsys()
        csrec = _mycs.torecord()
        shape =  _myia.shape()
        reg = _myrg.frombcs(
            csrec, shape,
            box, chans, stokes, "a", region
        )
        outia = _myia.deviation(
            outfile=outfile, region=reg, mask=mask,
            overwrite=overwrite, stretch=stretch, grid=grid,
            anchor=anchor, xlength=xlength, ylength=ylength,
            interp=interp, stattype=stattype, statalg=statalg,
            zscore=zscore, maxiter=maxiter
        )
        try:
            param_names = imdev.__code__.co_varnames[:imdev.__code__.co_argcount]
            if is_python3:
                vars = locals()
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
        _myia.done()
        _myrg.done()
        _mycs.done()
        if outia:
            outia.done()

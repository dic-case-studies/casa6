from __future__ import absolute_import
import sys

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import regionmanager, image
    from casatasks import casalog
    from .ialib import write_image_history

    _rg = regionmanager( )
else:
    from taskinit import *
    from ialib import write_image_history

    _rg = rgtool()
    image = iatool

def imsubimage(
    imagename, outfile, box, region, chans, stokes, mask, dropdeg,
    overwrite, verbose, stretch, keepaxes
):
    casalog.origin('imsubimage')
    myia = image()
    myia.dohistory(False)
    outia = None
    tmp_csys = None
    try:
        if (not myia.open(imagename)):
            raise RuntimeError("Cannot create image analysis tool using " + imagename)
        if (len(outfile) == 0):
            raise ValueError("outfile must be specified.")
        xregion = region
        if (type(region) != type({})):
            tmp_csys = myia.coordsys()
            xregion = _rg.frombcs(
                csys=tmp_csys.torecord(), shape=myia.shape(), box=box,
                chans=chans, stokes=stokes, stokescontrol="a", region=region
            )
            tmp_csys.done()
        outia = myia.subimage(
            outfile=outfile, region=xregion, mask=mask, dropdeg=dropdeg,
            overwrite=overwrite, list=verbose, stretch=stretch, keepaxes=keepaxes
        )
        try:
            param_names = imsubimage.__code__.co_varnames[:imsubimage.__code__.co_argcount]
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
        myia.done()
        _rg.done()
        if outia:
            outia.done()
        if tmp_csys:
            tmp_csys.done()
        

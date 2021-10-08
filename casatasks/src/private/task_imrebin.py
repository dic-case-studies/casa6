from __future__ import absolute_import
import sys

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image, regionmanager
    from casatasks import casalog
    from .ialib import write_image_history
else:
    from taskinit import *
    from ialib import write_image_history

    image = iatool
    regionmanager = rgtool

def imrebin(
    imagename, outfile, factor, region, box, chans, stokes, mask,
    dropdeg, overwrite, stretch, crop
):
    casalog.origin('imrebin')
    valid = True
    # because there is a bug in the tasking layer that allows float
    # arrays through when the spec is for intArray
    for x in factor:
        if x != int(x):
            valid = False
            break
    if not valid:
        for i in range(len(factor)):
            factor[i] = int(factor[i])   
        casalog.post(
            "factor is not an int array, it will be adjusted to "
                + str(factor),
            'WARN'
        )
    myia = image( )
    myia.dohistory(False)
    outia = None
    try:
        if (not myia.open(imagename)):
            raise RuntimeError("Cannot create image analysis tool using %s" % imagename)
        if (len(outfile) == 0):
            raise ValueError("outfile must be specified.")
        if (type(region) != type({})):
            myrg = regionmanager( )
            reg = myrg.frombcs(
                csys=myia.coordsys().torecord(), shape=myia.shape(), box=box,
                chans=chans, stokes=stokes, stokescontrol="a", region=region
            )
        else:
            reg = region
        outia = myia.rebin(
            outfile=outfile, bin=factor, region=reg, mask=mask, dropdeg=dropdeg,
            overwrite=overwrite, stretch=stretch, crop=crop
        )

        try:
            param_names = imrebin.__code__.co_varnames[:imrebin.__code__.co_argcount]
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
        if myia:
            myia.done()
        if outia:
            outia.done()
        if myrg:
            myrg.done()

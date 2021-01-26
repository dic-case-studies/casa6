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

def imfit(
    imagename, box, region, chans, stokes,
    mask, includepix, excludepix, residual,
    model, estimates, logfile, append,
    newestimates, complist, overwrite, dooff,
    offset, fixoffset, stretch, rms, noisefwhm,
    summary
):
    casalog.origin('imfit')
    myia = image()
    try:
        myia.dohistory(False)
        if (not myia.open(imagename)):
            raise Exception("Cannot create image analysis tool using %s" % imagename)
        result_dict = myia.fitcomponents(
            box=box, region=region, chans=chans, stokes=stokes,
            mask=mask, includepix=includepix,
            excludepix=excludepix, residual=residual,
            model=model, estimates=estimates, logfile=logfile,
            append=append, newestimates=newestimates,
            complist=complist, overwrite=overwrite, dooff=dooff,
            offset=offset, fixoffset=fixoffset, stretch=stretch,
            rms=rms, noisefwhm=noisefwhm, summary=summary
        )
        try:
            param_names = imfit.__code__.co_varnames[:imfit.__code__.co_argcount]
            if is_python3:
                vars = locals()
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]  
            for im in [residual, model]: 
                write_image_history(
                    im, sys._getframe().f_code.co_name,
                    param_names, param_vals, casalog
                )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')

        return result_dict

    finally:
        myia.done()

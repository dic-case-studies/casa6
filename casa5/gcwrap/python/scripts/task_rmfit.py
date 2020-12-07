import tempfile
import shutil

from casatasks.private.casa_transition import *
if is_CASA6:
    from casatasks import casalog
    from casatools import image, imagepol
    from .ialib import write_image_history
else:
    from taskinit import *
    from ialib import write_image_history
    image = iatool
    imagepol = potool

def rmfit(
    imagename, rm, rmerr, pa0, pa0err, nturns, chisq,
    sigma, rmfg, rmmax, maxpaerr
):
    casalog.origin('prom')
    myia = image()
    myia.dohistory(False)
    mypo = imagepol()
    tmpim = ""
    try:
        if len(imagename) == 0:
            raise ValueError("imagename must be specified.")
        if type(imagename) == type(['s']):
            # negative axis value means concatenate along spectral axis
            tmpim = tempfile.mkdtemp(suffix=".im", prefix="_rmfit_concat")
            myia = myia.imageconcat(
                outfile=tmpim, infiles=imagename, relax=True,
                axis=-1, overwrite=True
            )
            if not myia:
                raise RuntimeError("Unable to concatenate images.")
            myia.done()
            mypo.open(tmpim)
        else:
            if not mypo.open(imagename):
                raise RuntimeError("Cannot create image analysis tool using " + imagename)
        mypo.rotationmeasure(
            rm=rm, rmerr=rmerr, pa0=pa0, pa0err=pa0err, nturns=nturns, chisq=chisq,
            sigma=sigma, rmfg=rmfg, rmmax=rmmax, maxpaerr=maxpaerr
        )
        try:
            param_names = rmfit.func_code.co_varnames[:rmfit.func_code.co_argcount]
            param_vals = [eval(p) for p in param_names]
            for im in [rm, rmerr, pa0, pa0err, nturns, chisq]:
                write_image_history(
                    im, sys._getframe().f_code.co_name,
                    param_names, param_vals, casalog
                )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        # tasks no longer return bools
    finally:
        if (myia):
            myia.done()
        if (mypo):
            mypo.done()
        if len(tmpim) > 0:
            try:
                shutil.rmtree(tmpim)
            except Exception as exc:
                casalog.post("Could not remove " + tmpim + " because " + str(exc))

from __future__ import absolute_import

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import image, regionmanager, coordsys
    from casatasks import casalog
else:
    from taskinit import *

    image = iatool
    regionmanager = rgtool
    coordsys = cstool
    
def imstat(
    imagename, axes, region, box, chans,
    stokes, listit, verbose, mask, stretch,
    logfile, append, algorithm, fence, center,
    lside, zscore, maxiter, clmethod, niter
):
    _myia = image()
    _myrg = regionmanager()
    _mycs = coordsys()
    try:
        casalog.origin('imstat')
        _myia.open(imagename)
        _mycs = _myia.coordsys()
        csrec = _mycs.torecord()
        shape =  _myia.shape()
        reg = _myrg.frombcs(
            csrec, shape,
            box, chans, stokes, "a", region
        )
        return _myia.statistics(
            axes=axes, region=reg, list=listit,
            verbose=verbose, robust=True, mask=mask,
            stretch=stretch, logfile=logfile, append=append,
            algorithm=algorithm, fence=fence, center=center,
            lside=lside, zscore=zscore, maxiter=maxiter,
            clmethod=clmethod, niter=niter
        )

    finally:
        _myia.done()
        _myrg.done()
        _mycs.done()

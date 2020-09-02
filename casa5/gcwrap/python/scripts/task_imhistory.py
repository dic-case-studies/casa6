from __future__ import absolute_import

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import image
    from .. import casalog
else:
    from taskinit import *
    image = iatool

def imhistory(
    imagename, mode, verbose, origin, message
):
    _myia = image()
    try:
        casalog.origin('imhistory')
        _myia.open(imagename)
        if mode.startswith("l") or mode.startswith("L"):
            return _myia.history(verbose)
        elif mode.startswith("a") or mode.startswith("A"):
            return _myia.sethistory(origin=origin, history=message)
        raise ValueError("Unsopported mode " + mode)

    finally:
        _myia.done()

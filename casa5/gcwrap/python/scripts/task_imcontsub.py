from __future__ import absolute_import
import os
import re
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

def imcontsub(
    imagename, linefile, contfile, fitorder,
    region, box, chans, stokes
):
    casalog.origin('imcontsub')

    if ( len( linefile ) > 0 ):
        if ( os.path.exists( linefile ) ):
            raise ValueError('Error: file ' + linefile
                             +' already exists, please delete before continuing.',\
                             'SEVERE' )
    else:
        casalog.post("The linefile parameter is empty, consequently the"
                      +" spectral line image will NOT be\nsaved on disk.", \
                      'WARN')
            
    if ( len( contfile ) > 0 ):
            if ( os.path.exists( contfile ) ):
                raise ValueError('Error: Unable to continue file '+contfile\
                                 +' already exists, please delete before continuing.')
    else:
        casalog.post("The contfile parameter is empty, consequently the"
                      +" continuum image will NOT be\nsaved on disk.", \
                      'WARN')
    
    _myia = image()
    _myia.dohistory(False)
    _myia.open(imagename)
    mycsys = _myia.coordsys()
    if isinstance(box, list):
        box = ', '.join([str(b) for b in box])

    # Don't mix chans up with reg!  reg selects a subset for output, and chans
    # selects a subset to define the line-free channels.
    myrg = regionmanager()
    reg = myrg.frombcs(
        csys=mycsys.torecord(), shape=_myia.shape(),
        box=box, stokes=stokes, stokescontrol="f",
        region=region
    )
    channels = []
    if chans != None and len(chans) > 0:
        channels = myrg.selectedchannels(chans, _myia.shape())
    
    try:
        # Now do the continuum subtraction.
        lineim = _myia.continuumsub(
            outline=linefile, outcont=contfile, region=reg,
            channels=channels, fitorder=fitorder,
            overwrite=False
        )
        if not lineim:
            raise Exception("ia.continuumsub did not complete successfully")
        try:
            param_names = imcontsub.__code__.co_varnames[:imcontsub.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            for x in [lineim, contfile]:
                write_image_history(
                    x, sys._getframe().f_code.co_name,
                    param_names, param_vals, casalog
                )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        lineim.done()

    finally:
        _myia.done()
        if (lineim):
            lineim.done()
        if ( reg != None ):
            del reg

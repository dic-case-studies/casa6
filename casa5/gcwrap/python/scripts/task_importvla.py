from __future__ import absolute_import
import os

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import vlafiller, ms, agentflagger, table
    from casatasks import casalog
    from .mstools import write_history

    _ms = ms( )
    _filler = vlafiller( )
else:
    from taskinit import casac, casalog
    from taskinit import tbtool as table
    from mstools import write_history

    _filler = casac.vlafillertask()
    _ms = casac.ms()

    agentflagger = casac.agentflagger

def importvla(archivefiles,vis,bandname,frequencytol,project,starttime,
              stoptime,applytsys,autocorr,antnamescheme,keepblanks,evlabands):
    i=0
    overwrite=True
    ok = True
    try:
        casalog.origin('importvla')
        if ((type(vis)!=str) | (vis=='') | (os.path.exists(vis))):
            raise ValueError('Need valid visibility file name (bad name or already exists)')
        if (os.path.exists(vis)):
            raise ValueError('Visibility file already exists - remove or rename')
        for archivefile in archivefiles:
            if i>0: overwrite=False
            if ((type(archivefile)==str) & (os.path.exists(archivefile))):
                _filler.fill( msname=vis,inputfile=archivefile,overwrite=overwrite,
                              bandname=bandname,freqtol=frequencytol,project=project,
                              start=starttime,stop=stoptime,applytsys=applytsys,
                              keepautocorr=autocorr,antnamescheme=antnamescheme,
                              keepblanks=keepblanks,evlabands=evlabands )
                i += 1
            else:
                raise ValueError('Archive file not found - please verify the name')
    except Exception as exc:
        msg = "*** Error importing %s to %s: %s" % (archivefiles, vis, exc)
        raise RuntimeError(msg)

    nrows = 0
    try:
        _tb = table()
        ok &=_tb.open(vis)
        nrows =  _tb.nrows()
        _tb.done()
    except Exception:
        msg = "*** Error checking size of visibility file %s: %s" % (vis,instance)
        raise RuntimeError(msg)

    if nrows == 0:
        msg = "*** visibility file is empty: %s" % vis
        raise Exception(msg)

    # Write history
    try:
        param_names = importvla.__code__.co_varnames[:importvla.__code__.co_argcount]
        if is_python3:
            vars = locals( )
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        ok &= write_history(_ms, vis, 'importvla', param_names,param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % instance, 'WARN')

    # write initial flag version
    if ok:
        try:
            _af = agentflagger()
            ok &= _af.open(vis);
            ok &= _af.saveflagversion('Original', comment='Original flags at import into CASA', merge='replace')
            ok &= _af.done();
        except Exception:
            msg = "*** Error writing initial flag version of %s: %s" % (vis, instance)
            raise RuntimeError(msg)

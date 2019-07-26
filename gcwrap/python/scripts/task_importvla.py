from __future__ import absolute_import
from __future__ import print_function
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
    from mstools import write_history

    _filler = casac.vlafillertask()
    _ms = casac.ms()

    agentflagger = casac.agentflagger
    table = tbtool

def importvla(archivefiles,vis,bandname,frequencytol,project,starttime,
              stoptime,applytsys,autocorr,antnamescheme,keepblanks,evlabands):
    i=0
    overwrite=True
    ok = True
    try:
        casalog.origin('importvla')
        if ((type(vis)!=str) | (vis=='') | (os.path.exists(vis))):
            raise Exception('Need valid visibility file name (bad name or already exists)')
        if (os.path.exists(vis)): raise Exception('Visibility file already exists - remove or rename')
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
                raise Exception('Archive file not found - please verify the name')
    except Exception:
        casalog.post("*** Error importing %s to %s" % (archivefiles, vis), 'SEVERE')
        casalog.post("    %s" % instance, 'SEVERE')
	raise

    nrows = 0
    try:
        _tb = table()
        ok &=_tb.open(vis)
        nrows =  _tb.nrows()
	_tb.done()
    except Exception:
        casalog.post("*** Error checking size of visibility file %s: %s" % (vis,instance), 'SEVERE')
        raise

    if nrows == 0:
        msg = "*** visibility file is empty: %s" % vis
        casalog.post(msg, 'SEVERE')
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
    except:
	casalog.post("*** Error writing initial flag version of %s: %s" % (vis, instance), 'SEVERE')
	raise

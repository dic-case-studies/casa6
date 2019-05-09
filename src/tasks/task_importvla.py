from __future__ import absolute_import
from __future__ import print_function
import os

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import vlafiller, ms, agentflagger
    from casatasks import casalog
    from .mstools import write_history

    _ms = ms( )
    _af = agentflagger( )
    _filler = vlafiller( )
else:
    from taskinit import casac, casalog
    from mstools import write_history

    _filler = casac.vlafillertask()
    _af = casac.agentflagger()
    _ms = casa.ms()

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
        print('*** Error importing %s to %s:' % (archivefiles, vis))
        raise

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
    try:
        ok &= _af.open(vis);
        ok &= _af.saveflagversion('Original', comment='Original flags at import into CASA', merge='replace')
        ok &= _af.done();
    except:
        print('*** Error writing initial flag version of %s:' % vis)
        raise

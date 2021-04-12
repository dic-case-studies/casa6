from __future__ import absolute_import

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import miriadfiller
    from casatasks import casalog
    from .mstools import write_history
else:
    from taskinit import *
    from mstools import write_history

    miriadfiller = casac.miriadfiller

def importmiriad (
    mirfile=None,
    vis=None,
    tsys=None,
    spw=None,
    vel=None,
    linecal=None,
    wide=None,
    debug=None,
    ):
    """Convert a Miriad visibility file into a CASA visibility file (MS).
           The conversion of the Miriad visibility format into a measurement set.  This version
           has been tested for both ATNF and CARMA Miriad files.
................          
           Keyword arguments:
        mirfile -- Name of input Miriad visibility file (directory)
               default: none; example: mirfile='mydata.uv'

....   vis      -- Output ms name
               default: mirfile name with suffix replaced by '.ms'

....   tsys   -- Set True to use the system temperature to set the visibility weights
               default: False

....   spw -- specify the spectral windows to use
........ default: all

....   vel -- Velocity system to use: LSRK, LSRD or TOPO
....       default: TOPO for ATCA, LSRK for CARMA

....   linecal -- (CARMA only) apply CARMA linecal on the fly?
....       default: False

....   wide    -- (CARMA only) specify the window averages to use
........ default: all
........ 
....   debug  -- specify level of debug messages (0,1,2,3)
                 default: 0 (=None)

           
        """

    # Python script
    mymf = miriadfiller()
    try:
        try:
            casalog.origin('importmiriad')
            # -----------------------------------------'
            # beginning of importmiriad implementation
            # -----------------------------------------
            mymf.fill(vis,mirfile,tsys,spw,vel,linecal,wide,debug)
        except Exception as e:
            msg = "Failed to import miriad file %s" % mirfile
            raise RuntimeError(msg)

        # Write the args to HISTORY.
        try:
            param_names = importmiriad.__code__.co_varnames[:importmiriad.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            write_history(
                mymf, vis, 'importmiriad', param_names, 
                param_vals, casalog
            )
        except Exception:
            casalog.post("Failed to updated HISTORY", 'WARN')

    finally:
        if (mymf):
            del mymf 
        # -----------------------------------
        # end of importmiriad implementation
        # -----------------------------------

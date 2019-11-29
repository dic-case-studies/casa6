from __future__ import absolute_import
import os
import shutil
import string
import copy
import math

# get is_CASA6, is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from .mstools import write_history
    from casatools import table, ms, mstransformer
    from casatasks import casalog
    from .parallel.parallel_data_helper import ParallelDataHelper

    _tb = table()
else:
    from taskinit import mttool as mstransformer
    from taskinit import mstool as ms
    from taskinit import tbtool, casalog
    from mstools import write_history
    from parallel.parallel_data_helper import ParallelDataHelper

    _tb = tbtool()

def phaseshift(vis=None, 
               outputvis=None,
               keepmms=None,
               field=None,
               spw=None, 
               scan=None, 
               antenna=None, 
               correlation=None,
               timerange=None, 
               intent=None,
               array=None,
               uvrange=None,
               observation=None,
               feed=None,
               datacolumn=None, 
               ):

    """Changes the phase center for either short or large offsets/angles w.r.t. the original

    """

    casalog.origin('phaseshift')
    
    
    # Initiate the helper class    
    pdh = ParallelDataHelper("phaseshift", locals()) 

    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post('%s'%instance,'ERROR')
        return False

    # Input vis is an MMS
    if pdh.isMMSAndNotServer(vis) and keepmms:
        
        if not pdh.validateInputParams():        
            raise Exception('Unable to continue with MMS processing')
                        
        pdh.setupCluster('phaseshift')

        # Execute the jobs
        try:
            pdh.go()
        except Exception as instance:
            casalog.post('%s'%instance,'ERROR')
            return False
                    
        return True


    # Actual task code starts here

    # Create local copies of the MSTransform and ms tools
    mtlocal = mstransformer()
    mslocal = ms()

    try:
                    
        # Gather all the parameters in a dictionary.        
        config = {}
        
        config = pdh.setupParameters(inputms=vis, outputms=outputvis, field=field, 
                    spw=spw, array=array, scan=scan, antenna=antenna, correlation=correlation,
                    uvrange=uvrange,timerange=timerange, intent=intent, observation=observation,
                    feed=feed)
        
        
        # Check if CORRECTED column exists, when requested
        datacolumn = datacolumn.upper()
        if datacolumn == 'CORRECTED':
            _tb.open(vis)
            if 'CORRECTED_DATA' not in _tb.colnames():
                casalog.post('Input CORRECTED_DATA does not exist. Will use DATA','WARN')
                datacolumn = 'DATA'
            _tb.close()
             
        casalog.post('Will use datacolumn = %s'%datacolumn, 'DEBUG')
        config['datacolumn'] = datacolumn
        
        # Call MSTransform framework with phaseshift=True
        config['phaseshift'] = True

        # Configure the tool 
        casalog.post('%s'%config, 'DEBUG1')
        mtlocal.config(config)
        
        # Open the MS, select the data and configure the output
        mtlocal.open()
        
        # Run the tool
        casalog.post('Shift phase center')
        mtlocal.run()        
            
        mtlocal.done()
                    
    except Exception as instance:
        mtlocal.done()
        casalog.post('%s'%instance,'ERROR')
        return False

    # Write history to output MS, not the input ms.
    try:
        param_names = phaseshift.__code__.co_varnames[:phaseshift.__code__.co_argcount]
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        casalog.post('Updating the history in the output', 'DEBUG1')
        write_history(mslocal, outputvis, 'phaseshift', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),'WARN')
        return False

    mslocal = None
    
    return True
 
 

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

def hanningsmooth(vis=None, 
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

    """Hanning smooth frequency channel data to remove Gibbs ringing

    """

    casalog.origin('hanningsmooth')
    
    
    # Initiate the helper class    
    pdh = ParallelDataHelper("hanningsmooth", locals()) 

    # Validate input and output parameters
    pdh.setupIO()

    # Input vis is an MMS
    if pdh.isMMSAndNotServer(vis) and keepmms:
        
        if not pdh.validateInputParams():        
            raise Exception('Unable to continue with MMS processing')
                        
        pdh.setupCluster('hanningsmooth')

        # Execute the jobs
        pdh.go()
        return

    mslocal = ms()

    # Actual task code starts here
    try:    
        mtlocal = mstransformer()

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
        
        # Call MSTransform framework with hanning=True
        config['hanning'] = True

        # Configure the tool 
        casalog.post('%s'%config, 'DEBUG1')
        mtlocal.config(config)
        
        # Open the MS, select the data and configure the output
        mtlocal.open()
        
        # Run the tool
        casalog.post('Apply Hanning smoothing on data')
        mtlocal.run()        
            
    finally:
        mtlocal.done()

    # Write history to output MS, not the input ms.
    try:
        param_names = hanningsmooth.__code__.co_varnames[:hanningsmooth.__code__.co_argcount]
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        casalog.post('Updating the history in the output', 'DEBUG1')
        write_history(mslocal, outputvis, 'hanningsmooth', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),'WARN')

    mslocal = None
 
 

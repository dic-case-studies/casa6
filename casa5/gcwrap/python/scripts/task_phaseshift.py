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
    from casatools import measures as me
    from casatasks import casalog
    from .parallel.parallel_data_helper import ParallelDataHelper
else:
    from mstools import write_history
    from taskinit import tbtool,mstool,mttool,metool
    from taskinit import casalog
    from parallel.parallel_data_helper import ParallelDataHelper

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
               phasecenter=None,
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


    # Create local copies of tools (has to be here, otherwise 
    # ParallelDataHelper has a porblem digest the locals
    if is_CASA6:
        tblocal = table()
        mslocal = ms()
        mtlocal = mstransformer()
        melocal = me()
    else:
        tblocal = tbtool()
        mslocal = mstool()
        mtlocal = mttool()
        melocal = metool()
        
    # Actual task code starts here   
    try:
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
                tblocal.open(vis)
                if 'CORRECTED_DATA' not in tblocal.colnames():
                    casalog.post('Input CORRECTED_DATA does not exist. Will use DATA','WARN')
                    datacolumn = 'DATA'
                tblocal.close()
             
            casalog.post('Will use datacolumn = %s'%datacolumn, 'DEBUG')
            config['datacolumn'] = datacolumn
        
            # Call MSTransform framework with tviphaseshift=True
            config['tviphaseshift'] = True
            tviphaseshift_config = {}
            tviphaseshift_config['phasecenter'] = phasecenter
            config['tviphaseshiftlib'] = dict(tviphaseshift_config)

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
    
        # Update field table
        try:
        
            # Parse phase center string to obtain ra/dec
            dirstr = phasecenter.split(' ')
            try:
                thedir = melocal.direction(dirstr[0], dirstr[1], dirstr[2])
                thenewra_rad = thedir['m0']['value']
                thenewdec_rad = thedir['m1']['value']
            except Exception as instance:
                casalog.post("*** Error \'%s\' when interpreting parameter \'phasecenter\': "
                         % (instance), 'SEVERE')
                return False 
        
            # modify FIELD table 
            tblocal.open(outputvis + '/FIELD', nomodify=False)
            pcol = tblocal.getcol('PHASE_DIR')
            for row in range(0,tblocal.nrows()):
                pcol[0][0][row] = thenewra_rad
                pcol[1][0][row] = thenewdec_rad
            tblocal.putcol('PHASE_DIR', pcol)
            
            
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating FIELD" % (instance),'WARN')
            return False    

    finally:
        if tblocal:
            tblocal.done()
            tblocal = None
        if mslocal:
            mslocal = None
        if mtlocal:
            mtlocal.done()
            mtlocal = None
        if melocal:
            melocal.done()
            melocal = None
    
    return True
 
 

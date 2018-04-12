import os
import shutil
from CASAtools import mstransformer, table
from CASAtools import ms as mstool
from CASAtasks import casalog
from .parallel.parallel_data_helper import ParallelDataHelper
from .mstools import write_history

def cvel2(
    vis,
    outputvis,
    keepmms,
    passall,    # hidden parameter for backwards compatibiliy
    field,
    spw,
    scan,
    antenna,
    correlation,
    timerange,
    intent,
    array,
    uvrange,
    observation,
    feed,
    datacolumn,
    mode,
    nchan,
    start,
    width,
    interpolation,
    phasecenter,
    restfreq,
    outframe,
    veltype,
    hanning,
    ):
    
    """ This task used the MSTransform framework. It needs to use the ParallelDataHelper
        class, implemented in parallel.parallel_data_helper.py. 
    """

    assert outputvis != '', "Must provide output data set name in parameter outputvis."
    assert not os.path.exists(outputvis), "Output MS %s already exists - will not overwrite." % outputvis
    assert not os.path.exists(outputvis+".flagversions"), \
        "The flagversions \"%s.flagversions\" for the output MS already exist. Please delete." % outputvis

    # Initialize the helper class  
    pdh = ParallelDataHelper("cvel2", locals()) 

    casalog.origin('cvel2')
        
    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post('%s'%instance,'ERROR')
        return False

    # Input vis is an MMS
    if pdh.isParallelMS(vis) and keepmms:
        
        status = True   
        
        # Work the heuristics of combinespws=True and the separationaxis of the input             
        retval = pdh.validateInputParams()
        if not retval['status']:
            raise Exception('Unable to continue with MMS processing')
                        
        pdh.setupCluster('cvel2')

        # Execute the jobs
        try:
            status = pdh.go()
        except Exception as instance:
            casalog.post('%s'%instance,'ERROR')
            return status
                           
        return status


    # Create local copy of the MSTransform tool
    mtlocal = mstransformer( )
    tblocal = table( )

    try:
        # Gather all the parameters in a dictionary.        
        config = {}
        
        config = pdh.setupParameters(inputms=vis, outputms=outputvis, field=field, 
                    spw=spw, array=array, scan=scan, antenna=antenna, correlation=correlation,
                    uvrange=uvrange,timerange=timerange, intent=intent, observation=observation,
                    feed=feed)

        config['datacolumn'] = datacolumn
        casalog.post('Will work on datacolumn = %s'%datacolumn.upper())
        
        # In cvel the spws are always combined
        config['combinespws'] = True
        
        # Hanning smoothing
        config['hanning'] = hanning
        
        # Set the regridms parameter in mstransform
        config['regridms'] = True

        if passall == True:
            casalog.post('Parameter passall=True is not supported in cvel2','WARN')
        
        # Reset the defaults depending on the mode
        # Only add non-empty string parameters to config dictionary
        start, width = pdh.defaultRegridParams()
        config['mode'] = mode
        config['nchan'] = nchan
        if start != '':
            config['start'] = start
        if width != '':
            config['width'] = width

        config['interpolation'] = interpolation
        if restfreq != '':
            config['restfreq'] = restfreq
        if outframe != '':
            config['outframe'] = outframe
        if phasecenter != '':
            config['phasecenter'] = phasecenter
        config['veltype'] = veltype
        config['nspw'] = 1
        config['taql'] = "NOT (FLAG_ROW OR ALL(FLAG))"
        
        # Configure the tool and all the parameters        
        casalog.post('%s'%config, 'DEBUG')
        mtlocal.config(config)
        
        # Open the MS, select the data and configure the output
        mtlocal.open()
        
        # Run the tool
        mtlocal.run()        
            
        mtlocal.done()

    except Exception as instance:
        mtlocal.done()
        casalog.post('%s'%instance,'ERROR')
        return False

    # Write history to output MS, not the input ms.
    try:
        mslocal = mstool( )
        vars = locals( )
        param_names = cvel2.__code__.co_varnames[:cvel2.__code__.co_argcount]
        param_vals = [vars[p] for p in param_names]
        write_history(mslocal, outputvis, 'cvel2', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),'WARN')
        return False

    mslocal = None
    
    return True
        

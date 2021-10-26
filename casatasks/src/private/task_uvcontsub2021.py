from __future__ import absolute_import
import os

# get is_CASA6, is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import ms, mstransformer
    from casatasks import casalog
    from .mstools import write_history

else:
    from taskinit import mstool, mttool
    from taskinit import casalog
    from mstools import write_history

    mstransformer = mttool
    ms = mstool


def uvcontsub2021(vis=None, outputvis=None, field=None, spw=None,
                  scan=None, intent=None, array=None, observation=None,
                  datacolumn=None, fitspw=None, fitmethod=None, fitorder=None,
                  writemodel=None):

    """Continuum subtraction in the uv plane"""

    casalog.origin('uvcontsub2021')

    result = {}
    mtlocal = mstransformer()

    try:
        if not outputvis:
            raise ValueError('parameter outputvis: a name is required for the output MS')
        elif os.path.exists(outputvis):
            raise ValueError('Output MS ({}) already exists, refusing to overwrite it'.
                             format(outputvis))

        # Set up the mstransformer config from task parameters
        config = {'uvcontsub': True,
                  'inputms': vis,
                  'outputms': outputvis,
                  'field': field,
                  'spw': spw,
                  'array': array,
                  'scan': scan,
                  'intent': intent,
                  'observation': observation,
                  'datacolumn': datacolumn,
                  'reindex': False
                  }
        # Add uvcontsub-TVI config options
        uvcont_cfg = {'fitspw': fitspw,
                      'fitorder': fitorder,
                      'want_cont': False,
                      'denoising_lib': False,
                      'nthreads': 0,
                      'niter': 0
                      }
        config['uvcontsublib'] = uvcont_cfg

        # Configure the mstransformer tool
        casalog.post('mstransfom config: {}'.format(config), 'DEBUG')
        mtlocal.config(config)

        # Open the MS, select the data and configure the output
        mtlocal.open()

        # Run the tool
        casalog.post('Running continnum subtraction')
        result = mtlocal.run()
        mtlocal.done()

    finally:
        mtlocal.done()

    # Write history to output MS only. Input MS is read only
    try:
        mslocal = ms()
        param_names = uvcontsub2021.__code__.co_varnames[:uvcontsub2021.__code__.co_argcount]
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
            write_history(mslocal, outputvis, 'uvcontsub2021', param_names,
                          param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'{}\' updating HISTORY".format(instance), 'WARN')
    finally:
        mslocal.done()

    return result

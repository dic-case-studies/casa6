import os
from casatools import ms, mstransformer
from casatasks import casalog
from .mstools import write_history


def uvcontsub(vis=None, outputvis=None, field=None, spw=None,
              scan=None, intent=None, array=None, observation=None,
              datacolumn=None, fitspec=None, fitmethod=None, fitorder=None,
              writemodel=None):

    """Continuum subtraction in the uv plane"""

    casalog.origin('uvcontsub')

    casalog.post('The task uvcontsub has been updated. It does not yet support combining '
                 'spws for fitting. If you need this functionality please use '
                 'uvcontsub_old.')

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
        if fitmethod == 'gsl':
            denoising_lib = True
        elif fitmethod == 'casacore':
            denoising_lib = False
        else:
            raise ValueError(f'Unrecognized fit method: {fitmethod}')
        uvcont_cfg = {'fitspec': fitspec,
                      'denoising_lib': denoising_lib,
                      'niter': 1,
                      'fitorder': fitorder,
                      'writemodel': writemodel,
                      'want_cont': False,
                      'nthreads': None
                      }
        config['uvcontsublib'] = uvcont_cfg

        # Configure the mstransformer tool
        casalog.post('mstransfom config: {}'.format(config), 'DEBUG')
        mtlocal.config(config)

        # Open the MS, select the data and configure the output
        mtlocal.open()

        # Run the tool
        casalog.post('Running continnum subtraction')
        mt_result = mtlocal.run()

    finally:
        mtlocal.done()

    # Write history to output MS only. Input MS is read only
    try:
        mslocal = ms()
        param_names = uvcontsub.__code__.co_varnames[:uvcontsub.__code__.co_argcount]
        loc_vars = locals()
        param_vals = [loc_vars[p] for p in param_names]
        write_history(mslocal, outputvis, 'uvcontsub', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'{}\' updating HISTORY".format(instance), 'WARN')
    finally:
        mslocal.done()

    return mt_result['uvcontsub']

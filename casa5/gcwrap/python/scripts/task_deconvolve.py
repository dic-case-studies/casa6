from __future__ import absolute_import
from __future__ import print_function

import os
import shutil
import numpy
import copy
import time

# get is_CASA6 and is_python3, and import other classes
try:
    from casatasks.private.casa_transition import *
except:
    from sys import version_info
    is_python3 = version_info > (3,)
    is_CASA6 = is_python3
if is_CASA6:
    from casatasks import casalog

    from casatasks.private.imagerhelpers.imager_base import PySynthesisImager
    from casatasks.private.imagerhelpers.imager_deconvolver import PyDeconvolver
    from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
else:
    from taskinit import *

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.imager_deconvolver import PyDeconvolver
    from imagerhelpers.input_parameters import ImagerParameters

def deconvolve(
    ####### Data Selection
    #vis,#=''               -> not necessary: need to remove
    #selectdata,            -> not necessary: not used in tclean anywhere?
    #field,#='',               -> testing removal of allselpars
    #spw,#='',                 -> testing removal of allselpars
    #timerange,#='',           -> testing removal of allselpars
    #uvrange,#='',             -> testing removal of allselpars
    #antenna,#='',             -> testing removal of allselpars
    #scan,#='',                -> testing removal of allselpars
    #observation,#='',         -> testing removal of allselpars
    #intent,#='',           -> not necessary: not used in tclean anywhere?
    #datacolumn,#='corrected', -> testing removal of allselpars

    ####### Image definition
    imagename,#='',                                   -> now the first positional argument, indicating the dirty image to start from
    imsize,#=[100,100],                               -> kept: although this can be gleaned from the dirty image, we need to be able to sort out what the user wants in the case this value differs from what is in one or more of the residual images.
    cell,#=['1.0arcsec','1.0arcsec'],                 -> kept: although this can be gleaned from the dirty image, we need to be able to sort out what the user wants in the case this value differs from what is in one or more of the residual images.
    phasecenter,#='J2000 19:59:28.500 +40.44.01.50',  -> kept: although this can be gleaned from the dirty image, we need to be able to sort out what the user wants in the case this value differs from what is in one or more of the residual images.
    stokes,#='I',                                     -> kept: although this can be gleaned from the dirty image, we need to be able to sort out what the user wants in the case this value differs from what is in one or more of the residual images.
    projection,#='SIN',                               -> kept: although this can be gleaned from the dirty image, we need to be able to sort out what the user wants in the case this value differs from what is in one or more of the residual images.
    startmodel,#='',                                  -> kept: allows for a different model name than imagename

    ## Spectral parameters
    #specmode,#='mfs',  -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #reffreq,#='',      -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #nchan,#=1,         -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #start,#='',        -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #width,#='',        -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #outframe,#='LSRK', -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #veltype,#='',      -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #restfreq,#=[''],   -> testing removal of allimpars ? needed in the case that we're working with an image cube ?
    #interpolation,#='',             -> testing removal of allgridpars ? needed in the case that we're working with an image cube ?
    #perchanweightdensity, #='' -> not necessary: part of weightpars which aren't used
    ## 
    ####### Gridding parameters
    #gridder,#='ft',                 -> testing removal of allgridpars ? needed during image to uv gridding ?
    #facets,#=1,                     -> testing removal of allgridpars ? needed during image to uv gridding ?
    #psfphasecenter,#='',                             -> not necessary: used in generating PSF in tclean python code
    #chanchunks,#=1,                 -> testing removal of allgridpars ? needed during image to uv gridding ?

    #wprojplanes,#=1,                -> testing removal of allgridpars ? needed during image to uv gridding ?

    ### PB
    #vptable,#                       -> testing removal of allgridpars ? needed during image to uv gridding ?
    #mosweight, #=True               -> testing removal of allgridpars ? needed during image to uv gridding ?
    #aterm,#=True,                   -> testing removal of allgridpars ? needed during restoration ? 
    #psterm,#=True,                  -> testing removal of allgridpars ? needed during restoration ?
    #wbawp ,#= True,                 -> testing removal of allgridpars ? needed during restoration ?
    #conjbeams ,#= True,             -> testing removal of allgridpars ? needed during restoration ?
    #cfcache ,#= "",                 -> testing removal of allgridpars ? needed during restoration ?
    #usepointing, #=false            -> testing removal of allgridpars
    #computepastep ,#=360.0,         -> testing removal of allgridpars
    #rotatepastep ,#=360.0,          -> testing removal of allgridpars
    #pointingoffsetsigdev ,#=[10.0], -> testing removal of allgridpars

    #pblimit,#=0.01,         -> not necessary: part of allnormpars which aren't used
    #normtype,#='flatnoise', -> not necessary: part of allnormpars which aren't used

    ####### Deconvolution parameters
    deconvolver,#='hogbom',
    scales,#=[],
    nterms,#=1,             -> kept: needed for deconvolver to find multi-term images on disk
    smallscalebias,#=0.0

    ### restoration options
    restoration,
    restoringbeam,#=[],
    # pbcor,

    ##### Outliers
    # outlierfile,#='',     -> not necessary: can manually execute on all outlier images?

    ##### Weighting
    #weighting,#='natural',     -> not necessary: part of weightpars which aren't used
    #robust,#=0.5,              -> not necessary: part of weightpars which aren't used
    #noise,#0.0Jy               -> not necessary: part of weightpars which aren't used
    #npixels,#=0,               -> not necessary: part of weightpars which aren't used
    #uvtaper,#=[],              -> not necessary: part of weightpars which aren't used


    ##### Iteration control
    niter,#=0, 
    gain,#=0.1,
    threshold,#=0.0, 
    nsigma,#=0.0
    #cycleniter,#=0, since we're only doing one minor loop, this can be a mirror of niter
    cyclefactor,#=1.0,
    minpsffraction,#=0.1,
    maxpsffraction,#=0.8,
    interactive,#=False, TODO test with TRUE
    plotReport,#=False
    iterrec,#=False

    ##### (new) Mask parameters
    usemask,#='user',
    mask,#='',
    pbmask,#='',

    ##### automask by multithresh
    sidelobethreshold,#=5.0,
    noisethreshold,#=3.0,
    lownoisethreshold,#=3.0,
    negativethreshold,#=0.0,
    smoothfactor,#=1.0,
    minbeamfrac,#=0.3, 
    cutthreshold,#=0.01,
    growiterations,#=100
    dogrowprune,#=True
    #minpercentchange,#=0.0 -> not necessary: no major loops
    verbose, #=False
    fastnoise): #=False
    """
    Runs the minor cycle only of tclean.
    Most of this code is copied directly from tclean.
    """

    ## Misc

    # restart):#=True,
    #savemodel,#="none",    -> not necessary: should be done in the major cycle

    ####### State parameters
    #parallel):#=False) -> not necessary: no parallelized version of deconvolver-only task

    # used in PySynthesisImager:
    # * allselpars  (initializeImagers) -> not necessary
    # * allgridpars (initializeImagers) -> not necessary
    # * allimpars   (__init__, initializeImagers, estimatememory, hasConverged, runMinorCycleCore)
    # * alldecpars  (initializeDeconvolvers, hasConverged, runMinorCycleCore)
    # * iterpars    (initializeIterationControl)
    #
    # these translate as:
    # * allselpars:  vis, field, spw, scan, timerange, uvrange, antenna, observation, state, datacolumn, savemodel
    # * allgridpars: gridder, aterm, psterm, mterm, wbawp, cfcache, usepointing, dopbcorr, conjbeams, computepastep, rotatepastep,
    #                pointingoffsetsigdev, facets, chanchunks, interpolation, wprojplanes, deconvolver,
    #                vptable, gridfunction, convsupport, truncate, gwidth, jwidth, minweight, clipminmax, imagename
    # * allimpars:   imagename, nchan, imsize, cell, phasecenter, stokes, specmode, start, width, veltype, nterms, restfreq, outframe,
    #                reffreq, sysvel, sysvelframe, projection, restart, startmodel, deconvolver
    # ** ^^notes^^ just need imagename and imsize, I think?
    # * alldecpars:  deconvolver, nterms, scales, smallscalebias, restoringbeam, usemask, mask, pbmask, maskthreshold, maskresolution,
    #                nmask, sidelobethreshold, noisethreshold, lownoisethreshold, negativethreshold,
    #                smoothfactor, minbeamfrac, cutthreshold, growiterations, dogrowprune, minpercentchange, verbose, fastnoise,
    #                interactive, startmodel, nsigma,  imagename
    # * iterpars:    niter, cycleniter, threshold, gain, interactive, cyclefactor, minpsffraction, maxpsffraction, savemodel, nsigma

    # clean input
    inp=locals().copy()
    inp['msname']      = '' # -> no 'vis' parameter for minor cycle only
    inp['cycleniter']  = inp['niter']
    inp['loopgain']    = inp.pop('gain')
    inp['scalebias']   = inp.pop('smallscalebias')

    #####################################################
    #### Construct ImagerParameters and Imager objects
    #####################################################
    
    # make a list of parameters with defaults from tclean
    if is_python3:
        defparm=dict(list(zip(ImagerParameters.__init__.__code__.co_varnames[1:], ImagerParameters.__init__.__defaults__)))
    else:
        defparm=dict(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:], ImagerParameters.__init__.func_defaults))

    ## assign values to the ones passed to deconvolve and if not defined yet in deconvolve...
    ## assign them the default value of the constructor
    bparm={k:  inp[k] if k in inp else defparm[k]  for k in defparm.keys()}

    ## create the parameters list help object
    paramList=ImagerParameters(**bparm)

    ## Setup Imager object
    decon = PyDeconvolver(params=paramList)

    #####################################################
    #### Run the minor cycle
    #####################################################

    iterrecout = False
    isit = 0
    retrec = ''
    try:

        #################################################
        #### Setup
        #################################################

        ## Init minor cycle elements
        # print("initializing deconvolver")
        t0=time.time();
        decon.initializeDeconvolvers()
        ####now is the time to check estimated memory
        decon.estimatememory()
        ## setup iteration controller
        decon.initializeIterationControl()
        if type(iterrec) == dict:
            PyDeconvolver.mergeIterRecords(decon, iterrec)
        t1=time.time();
        casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");

        #################################################
        #### Exec
        #################################################

        if niter > 0:
            ## Set up the internal state of the iterater and automask
            # is this necessary? -> I think so ~bgb200731
            isit = decon.hasConverged()
            decon.updateMask()

            isit = decon.hasConverged() # here in case updateMaskMinor() produces an all-false mask
            if not isit:
                # print("running minor cycle");
                t0=time.time();
                decon.runMinorCycle()
                t1=time.time();
                casalog.post("***Time for minor cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");
                isit = decon.hasConverged() # get the convergence state, to report back to the calling code

            ## Get summary from iterbot
            if type(interactive) != bool:
                retrec=imager.getSummary();

        #################################################
        #### Teardown
        #################################################

        ## Get records from iterbot, to be used in the next call to deconvolve
        iterrecout = decon.getIterRecords()

        # TODO this restoration step is now in at least three different tasks (sdintimaging, tclean, deconvolve). Should it be moved into common code somewhere?
        ## Restore images.
        if restoration==True:  
            t0=time.time();
            decon.restoreImages()
            t1=time.time();
            casalog.post("***Time for restoring images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");

        ##close tools
        decon.deleteTools()

    except Exception as e:
        casalog.post('Exception from deconvolve : ' + str(e), "SEVERE", "deconvolve")
        if decon != None:
            decon.deleteTools() 

        larg = list(e.args)
        larg[0] = 'Exception from deconvolve : ' + str(larg[0])
        e.args = tuple(larg)
        raise

    return iterrecout, isit, retrec
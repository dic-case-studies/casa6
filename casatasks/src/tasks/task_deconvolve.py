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
    print("is_python3: " + str(is_python3) + ", is_CASA6: " + str(is_CASA6))
if is_CASA6:
    from casatasks import casalog

    from casatasks.private.imagerhelpers.imager_base import PySynthesisImager
    from casatasks.private.imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from casatasks.private.imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
else:
    from taskinit import *

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from imagerhelpers.input_parameters import ImagerParameters

def beanvolve(
    ####### Data Selection
    #vis,#=''                                         -> not necessary: need to remove
    #selectdata,                                      -> not necessary: not used in tclean anywhere?
    #field,#='',               -> testing removal of allselpars
    #spw,#='',                 -> testing removal of allselpars
    #timerange,#='',           -> testing removal of allselpars
    #uvrange,#='',             -> testing removal of allselpars
    #antenna,#='',             -> testing removal of allselpars
    #scan,#='',                -> testing removal of allselpars
    #observation,#='',         -> testing removal of allselpars
    #intent,#='',                                     -> not necessary: not used in tclean anywhere?
    #datacolumn,#='corrected', -> testing removal of allselpars

    ####### Image definition
    imagename,#='',                                   -> now the first positional argument, indicating the dirty image to start from
    #imsize,#=[100,100],                              -> not necessary: can be gleaned from the size of the dirty image created in the major cycle from tclean
    #cell,#=['1.0arcsec','1.0arcsec'],                -> not necessary: already set in major cycle of tclean?
    #phasecenter,#='J2000 19:59:28.500 +40.44.01.50', -> not necessary: already set in major cycle of tclean?
    #stokes,#='I',                                    -> not necessary: already set in major cycle of tclean?
    #projection,#='SIN',                              -> not necessary: already set in major cycle of tclean?
    startmodel,#='',                                  -> kept: if the model gets updated then this allows for a choice in where to put the updated model TODO is this true?

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
    #nterms,#=1,        -> testing removal of allimpars
    smallscalebias,#=0.0

    ### restoration options
    restoration,
    restoringbeam,#=[],
    pbcor,

    ##### Outliers
    outlierfile,#='',

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
    fastnoise, #=False

    ## Misc

    restart,#=True,
    savemodel,#="none",

    ####### State parameters
    parallel):#=False)

    # used in PySynthesisImager:
    # * allselpars  (initializeImagers)
    # * allgridpars (initializeImagers)
    # * allimpars   (__init__, initializeImagers, estimatememory, hasConverged, runMinorCycleCore)
    # * alldecpars  (initializeDeconvolvers, hasConverged, runMinorCycleCore)
    # * iterpars    (initializeIterationControl)
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
    inp['msname']      = '' #                   -> no 'vis' parameter for minor cycle only
    inp['cycleniter']  = inp['niter']
    # inp['timestr']   = inp.pop('timerange')   -> testing removal of allselpars
    # inp['uvdist']    = inp.pop('uvrange')     -> testing removal of allselpars
    # inp['obs']       = inp.pop('observation') -> testing removal of allselpars
    # inp['state']     = 'inp.pop('intent')     -> not necessary: not used in tclean anywhere?
    inp['loopgain']    = inp.pop('gain')
    inp['scalebias']   = inp.pop('smallscalebias')

    #####################################################
    #### Construct ImagerParameters and Imager objects
    #####################################################
    
    # make a list of parameters with defaults from tclean
    # TODO this pcube check is now in at least three different tasks (sdintimaging, tclean, deconvolve). Should it be moved into common code somewhere?
    if is_python3:
        defparm=dict(list(zip(ImagerParameters.__init__.__code__.co_varnames[1:], ImagerParameters.__init__.__defaults__)))
    else:
        defparm=dict(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:], ImagerParameters.__init__.func_defaults))

    ## assign values to the ones passed to tclean and if not defined yet in tclean...
    ## assign them the default value of the constructor
    ## stolen directly from tclean
    bparm={k:  inp[k] if k in inp else defparm[k]  for k in defparm.keys()}

    # TODO this pcube check is now in at least two different tasks (tclean, deconvolve). Should it be moved into common code somewhere?
    pcube = False
    if parallel==True and specmode!='mfs':
        pcube=True
        parallel=False

    # TODO this bparm check is now in at least three different tasks (sdintimaging, tclean, deconvolve). Should it be moved into common code somewhere?
    # TODO not necessary if we're removing it along with the rest of the gridding parameters
    ###default mosweight=True is tripping other gridders as they are not
    ###expecting it to be true
    # if(bparm['mosweight']==True and bparm['gridder'].find("mosaic") == -1):
    #     bparm['mosweight']=False

    ## create the parameters list help object
    paramList=ImagerParameters(**bparm)

    ## Setup Imager objects, for different parallelization schemes.
    imager=None
    imagerInst=PySynthesisImager
    if parallel==False and pcube==False:
         imager = PySynthesisImager(params=paramList,activities=PySynthesisImager.minorActivities)
         imagerInst=PySynthesisImager
    elif parallel==True:
         imager = PyParallelContSynthesisImager(params=paramList)
         imagerInst=PyParallelContSynthesisImager
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params=paramList)
         imagerInst=PyParallelCubeSynthesisImager
         # virtualconcat type - changed from virtualmove to virtualcopy 2016-07-20
         #using ia.imageconcat now the name changed to copyvirtual 2019-08-12
         concattype='copyvirtual'
    else:
         print('Invalid parallel combination in doClean.')
         return False

    #####################################################
    #### Run the minor cycle
    #####################################################

    retrec = {}
    try:
        # quick inputs check
        if niter == 0:
            return True

        #################################################
        #### Setup
        #################################################

        ## Init minor cycle elements
        print("initializing imager")
        t0=time.time();
        autoInit=True # TODO testing/discussing initializeAuto
        if autoInit:
            imager.initializeAuto()
        else:
            # imager.initializeImagers()     # -> not needed with activities changes to PySynthesisImager
            # imager.initializeNormalizers() # -> not needed for minor cycle
            # imager.setWeighting()          # -> not needed for minor cycle
            imager.initializeDeconvolvers()
            ####now is the time to check estimated memory
            imager.estimatememory()
            ## setup iteration controller
            imager.initializeIterationControl()
        t1=time.time();
        casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");

        ## Make PSF
        ### Don't need to make the PSF -> The first run
        ### of the major cycle should have done this.
        # if calcpsf==True:
        #     print("building PSF")
        #     imager.makePSF()
        #     if((psfphasecenter != '') and (gridder=='mosaic')):
        #         print("doing with different phasecenter psf")
        #         imager.unlockimages(0)
        #         psfParameters=paramList.getAllPars()
        #         psfParameters['phasecenter']=psfphasecenter
        #         psfParamList=ImagerParameters(**psfParameters)
        #         psfimager=imagerInst(params=psfParamList)
        #         psfimager.initializeImagers()
        #         psfimager.setWeighting()
        #         psfimager.makeImage('psf', psfParameters['imagename']+'.psf')
        #     imager.makePB()

        ## Set up the internal state of the iterater and automask
        # is this necessary? -> I think so ~bgb200731
        isit = imager.hasConverged()
        imager.updateMask()

        #################################################
        #### Exec
        #################################################

        if not imager.hasConverged(): # here in case updateMask() produces an all-false mask
            print("running minor cycle");
            t0=time.time();
            imager.runMinorCycle()
            t1=time.time();
            casalog.post("***Time for minor cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");

        #################################################
        #### Teardown
        #################################################

        ## Get summary from iterbot
        if type(bparm['interactive']) != bool:
            retrec=imager.getSummary();

        # TODO this restoration step is now in at least three different tasks (sdintimaging, tclean, deconvolve). Should it be moved into common code somewhere?
        ## Restore images.
        if restoration==True:  
            t0=time.time();
            imager.restoreImages()
            t1=time.time();
            casalog.post("***Time for restoring images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");
            if pbcor==True:
                t0=time.time();
                imager.pbcorImages()
                t1=time.time();
                casalog.post("***Time for pb-correcting images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");

        ##close tools
        # needs to deletools before concat or lock waits for ever
        imager.deleteTools()
   
        # TODO this concatImages is now in at least three different tasks (sdintimaging, tclean, deconvolve). Should it be moved into common code somewhere?
        if (pcube):
            print("running concatImages ...")
            casalog.post("Running virtualconcat (type=%s) of sub-cubes" % concattype,"INFO2", "task_deconvolve")
            imager.concatImages(type=concattype)

        # TODO this casalog.post is now in at least three different tasks (sdintimaging, tclean, deconvolve). Should it be moved into common code somewhere?        
        # CAS-10721 
        if niter>0 and savemodel != "none":
            casalog.post("Please check the casa log file for a message confirming that the model was saved after the last major cycle. If it doesn't exist, please re-run tclean with niter=0,calcres=False,calcpsf=False in order to trigger a 'predict model' step that obeys the savemodel parameter.","WARN","task_deconvolve")

    except Exception as e:
        casalog.post('Exception from beanvolve : ' + str(e), "SEVERE", "beanvolve")
        if imager != None:
            imager.deleteTools() 

        larg = list(e.args)
        larg[0] = 'Exception from beanvolve : ' + str(larg[0])
        e.args = tuple(larg)
        raise

    return retrec;


if __name__ == "__main__":
    if True:#('beanvolve_args' not in vars()):
        print("importing beanvolve_args")
        exec(open('./beanvolve_args.py').read())

    print("executing beanvolve")
    ret = beanvolve(**(beanvolve_args.copy()))
    globals()['beanvolve_ret'] = ret
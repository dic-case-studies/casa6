from __future__ import absolute_import
from __future__ import print_function

import time
import numpy
import os
import shutil
import re

# get is_CASA6 and is_python3, and import other classes
try:
    from casatasks.private.casa_transition import *
except:
    from sys import version_info
    is_python3 = version_info > (3,)
    is_CASA6 = is_python3
if is_CASA6:
    from casatasks import casalog

    from casatasks.private.imagerhelpers.imager_deconvolver import PyDeconvolver
    from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
else:
    from taskinit import *

    from imagerhelpers.imager_deconvolver import PyDeconvolver
    from imagerhelpers.input_parameters import ImagerParameters
    from imregrid import imregrid

def check_requiredimgs_exist(imagename, deconvolver, nterms):
    # get the list of images to check for
    reqims = []
    if deconvolver == 'mtmfs':
        end = nterms*2-1
        for ttn in range(0, end):
            ttext = ".tt" + str(ttn)
            if ttn != end-1:
                reqims.append(imagename + ".residual" + ttext)
            reqims.append(imagename + ".psf" + ttext)
    else:
        reqims.append(imagename + ".residual")
        reqims.append(imagename + ".psf")

    # find images that exist on disk
    extims = list(filter(lambda im: os.path.exists(im), reqims))

    # verify required images are available
    if len(extims) != len(reqims):
        diffims = list(filter(lambda im: im not in extims, reqims))
        raise RuntimeError("Internal Error: missing one or more of the required images: " + str(diffims))

def check_starmodel_model_collisions(startmodel, imagename, deconvolver):
    # check for startmodel(s)
    startmodels = []
    if type(startmodel) is str:
        if len(startmodel) > 0:
            startmodels = [startmodel]
    else: # if type(startmodel) is list
        startmodels = startmodel

    # verify the existance of startmodel(s), and map to imagename.model
    ttn = 0
    sm_modim_map = []
    for sm in startmodels:
        sm = os.path.normpath(sm)
        smdir = os.path.dirname(sm)
        smbase = os.path.basename(sm)

        # verify startmodel exists
        # Note: this check should be unneccessary (should be done in cpp), but
        # current tests indicate that cpp code does not catch this case.
        if not os.path.exists(sm):
            raise RuntimeError("Internal Error: parameter startmodel set to \"{0}\" but that file does not exist".format(sm))

        # get the path to the destination model
        if ".model" in smbase:
            ext = re.search(r'(\.model.*)', smbase).group(0)
        elif deconvolver == 'mtmfs':
            ext = ".model.tt{0}".format(ttn)
            ttn += 1
        else:
            ext = ".model"
        modim = os.path.join(smdir, imagename+ext)
        sm_modim_map.append([sm, modim])

        # check if both startmodel is set and imagename.model exists
        # Note: this check should be unneccessary (should be done in cpp), but
        # current tests indicate that cpp code does not catch this case.
        if os.path.exists(modim):
            raise RuntimeError("Internal Error: imagename.model already exists! Either parameter startmodel must not be set ('') or imagename.model ({0}) must not exist.".format(modim) +
                               os.linesep+"\tEither unset startmodel or remove {0} to continue".format(modim))

    return sm_modim_map

def deconvolve(
    ####### Data Selection
    imagename,#='',
    startmodel,#='',

    ####### Deconvolution parameters
    deconvolver,#='hogbom',
    scales,#=[],
    nterms,#=1,
    smallscalebias,#=0.0

    ### restoration options
    restoration,#=True,
    restoringbeam,#=[],

    ##### Iteration control
    niter,#=0, 
    gain,#=0.1,
    threshold,#=0.0, 
    nsigma,#=0.0
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
    verbose, #=False
    fastnoise): #=False
    """
    Runs the minor cycle only of tclean.
    Most of this code is copied directly from tclean.
    """

    # discard empty start model strings
    if type(startmodel) is list:
        startmodel = list(filter(lambda v: len(v) > 0, startmodel))

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

    # Make sure that we have all the necessary images and that the startmodel is valid.
    # Note: cpp code should check that .residual and .psf exist, but current tests indicate that it doesn't do that.
    check_requiredimgs_exist(imagename, deconvolver, nterms)
    check_starmodel_model_collisions(startmodel, imagename, deconvolver)

    #####################################################
    #### Run the minor cycle
    #####################################################

    ## Setup Imager object
    decon = PyDeconvolver(params=paramList)

    iterrec = False
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
        t1=time.time();
        casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_deconvolve");

        #################################################
        #### Exec
        #################################################

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
            retrec=decon.getSummary();

        #################################################
        #### Teardown
        #################################################

        ## Get records from iterbot, to be used in the next call to deconvolve
        iterrec = decon.getIterRecords()

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

    return iterrec, isit, retrec
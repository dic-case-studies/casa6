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

    from casatools import image
    from casatasks.private.imagerhelpers.imager_deconvolver import PyDeconvolver
    from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
    # from casatasks import imregrid

    # imregrid = imregrid.imregrid
else:
    from taskinit import *
    image = iatool

    from imagerhelpers.imager_deconvolver import PyDeconvolver
    from imagerhelpers.input_parameters import ImagerParameters
    from imregrid import imregrid

def summaries_diff(sum1, sum2):
    """
    Does the rough comparison betwen two summaries to determine if their
    coordinate systems are probably equal.
    """
    tocheck = ['ndim', 'axisnames', 'axisunits', 'incr', 'refpix', 'refval', 'shape']#, 'tileshape'] TODO anything else?
    diff = []
    for k in tocheck:
        if type(sum1[k]) is numpy.ndarray:
            if not numpy.array_equal(sum1[k], sum2[k]):
                diff.append(k)
        else:
            if sum1[k] != sum2[k]:
                diff.append(k)
    return diff

def differing_coords(*images):
    """
    Compares the coordinate systems between a set of images to make sure they are all the same.
    Each input should be the name/path to a table on the file system.
    Returns the first two images that have differing coordinate systems, if any.
    Returns FALSE if all coordinate systems are equal.
    """

    # can't have differing coordinate systems if the there aren't any images to compare
    if len(images) < 2:
        return []

    # get the coordinate system of the first image to compare to all other images
    _ia = image()
    _ia.open(images[0])
    basesum = _ia.summary(list=False)
    _ia.close()

    # compare to each other image
    # don't actually compare coordinate systmes (cpp method CoordinateSystem::near not exposed, just to a rough check using the summaries)
    for im in images[1:]:
        _ia = image()
        _ia.open(im)
        othersum = _ia.summary(list=False)
        _ia.close()
        sum_diff = summaries_diff(basesum, othersum)
        if len(sum_diff) > 0:
            return [images[0], im, sum_diff]

    # if we've reached this point, then all coordinate systems are equal
    return []

def get_image_info(im):
    _ia = image()
    _ia.open(im)
    csys = _ia.coordsys()
    shape = _ia.shape()
    names = csys.names()
    _ia.close()
    return names, shape, csys.torecord()

def decon_regrid(source, dest, axes, shape, csys):
    _ia = image()
    _ia.open(source)

    source_axes = _ia.coordsys().names()
    if source_axes != axes:
        nl = os.linesep
        iatoolname = "image" if is_CASA6 else "iatool"
        raise RuntimeError("Internal Error: can't regrid image {0} to {1}, axes don't match".format(source, dest) +
                           nl+"\tAxes should be {0}, but instead are {1}".format(axes, source_axes) +
                           nl+"\tMaybe use the task \"imtrans\" or the tool function \"{0}.adddegaxes\" to fix the axes".format(iatoolname))

    _ia.regrid(outfile=dest, shape=shape, csys=csys)
    _ia.close()

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

    #####################################################
    #### Validate the inputs
    #####################################################

    # verify required images are available
    all_inpims = [imagename + v for v in [".residual", ".psf", ".model", ".mask", ".pb"]]
    real_inpims = list(filter(lambda im: os.path.exists(im), all_inpims))
    resim = imagename + ".residual"
    required_images = [resim, imagename+".psf"]
    if (required_images[0] not in real_inpims) or (required_images[1] not in real_inpims):
        raise RuntimeError("Internal Error: missing one of the required images: " + str(required_images))

    # verify the coordinate axes of the residual image
    # all other images will have their coordinate axes confirmed when comparing coordinate systems to the residual
    nl = os.linesep
    resaxes, resshape, rescsys = get_image_info(resim)
    reqaxes = ['Right Ascension', 'Declination', 'Stokes', 'Frequency']
    if (resaxes != reqaxes):
        iatoolname = "image" if is_CASA6 else "iatool"
        raise RuntimeError("Internal Error: the image \"{0}\" has axes {1} when it should have axes {2}".format(resim, imaxes, reqaxes) +
                           nl+"\tMaybe use the task \"imtrans\" or the tool function \"{0}.adddegaxes\" (also available with task \"importfits\") to fix the axes".format(iatoolname))

    # verify images have the correct coordinate systems
    images_diff = differing_coords(*real_inpims)
    if len(images_diff) > 0:
        raise RuntimeError("Internal Error: the images \"{0}\" and \"{1}\" have different coordinate systems in at least the following ways: {2}".format(images_diff[0], images_diff[1], images_diff[2]) +
                           nl+"\tUse the \"imregrid\" task to translate all input images to the same coordinate system.")

    # check for startmodel(s)
    startmodels = []
    if type(startmodel) is str:
        if len(startmodel) > 0:
            startmodels = [startmodel]
    else: # if type(startmodel) is list
        startmodels = startmodel

    # verify the existance of startmodel(s), and map to imagename.model
    ttn = 0
    sms = []
    for sm in startmodels:
        sm = os.path.normpath(sm)
        smdir = os.path.dirname(sm)
        smbase = os.path.basename(sm)

        # verify startmodel exists
        if not os.path.exists(sm):
            raise RuntimeError("Internal Error: parameter startmodel set to \"{0}\" but that file does not exist".format(sm))

        # get the path to the destination model
        if ".model" in smbase:
            ext = re.search(r'(\.model.*)', smbase).group(0)
        elif len(startmodels) > 1:
            ext = ".model.tt{0}".format(ttn)
            ttn += 1
        else:
            ext = ".model"
        modim = os.path.join(smdir, imagename+ext)
        sms.append([sm, modim])

        # check if both startmodel is set and imagename.model exists
        if os.path.exists(modim):
            raise RuntimeError("Internal Error: imagename.model already exists! Either parameter startmodel must not be set ('') or imagename.model ({0}) must not exist.".format(modim) +
                               nl+"\tEither unset startmodel or remove {0} to continue".format(modim))

    # copy, or convert, the startmodel(s) to the mapped imagename.model
    for sm_modim in sms:
        sm = sm_modim[0]    # startmodel .model image path
        modim = sm_modim[1] # imagename.model image path

        # copy startmodel to imagename.model, possibly regridding to the residual image's coordinate system
        startmodel_diff = differing_coords(resim, sm)
        if len(startmodel_diff) > 0:
            casalog.post("Regridding {0} to {1} (issues: {2})".format(sm, modim, startmodel_diff[2]))
            decon_regrid(source=sm, dest=modim, axes=resaxes, shape=resshape, csys=rescsys)
        else:
            casalog.post("Copying {0} to {1}")
            shutil.copytree(sm, modim)

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
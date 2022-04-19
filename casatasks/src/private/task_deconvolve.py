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
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from .cleanhelper import write_tclean_history, get_func_params
    from casatools import synthesisimager
    ia = image( )
else:
    from taskinit import *

    from imagerhelpers.imager_deconvolver import PyDeconvolver
    from imagerhelpers.input_parameters import ImagerParameters
    from imregrid import imregrid
    from parallel.parallel_task_helper import ParallelTaskHelper
    from cleanhelper import write_tclean_history, get_func_params
    synthesisimager=casac.synthesisimager
    ia = iatool( )

try:
    if is_CASA6:
        from casampi.MPIEnvironment import MPIEnvironment
        from casampi import MPIInterface
    else:
        from mpi4casa.MPIEnvironment import MPIEnvironment
        from mpi4casa import MPIInterface
    mpi_available = True
except ImportError:
    mpi_available = False

def check_requiredmask_exists(usemask, mask):
    if usemask != "user":        # don't use the mask parameter
        return
    if type(mask) is type([]):   # mask is an array of values
        return
    if mask == "":               # mask is an empty string -> no file specified
        return
    if "[" in mask:              # mask is a region string
        return
    if not os.path.exists(mask): # mask is a filename string <- only option left
        raise RuntimeError("Internal Error: 'mask' parameter specified as a filename '"+mask+"', but no such file exists")

def check_requiredimgs_exist(imagename, inp):
    # get the list of images to check for
    reqims = []
    if inp['deconvolver'] == 'mtmfs':
        nterms = inp['nterms']
        end = nterms*2-1
        for ttn in range(0, nterms*2-1):
            reqims.append(imagename + ".psf" + ".tt" + str(ttn))
        for ttn in range(0, nterms):
            reqims.append(imagename + ".residual" + ".tt" + str(ttn))
    else:
        reqims.append(imagename + ".residual")
        reqims.append(imagename + ".psf")

    # find images that exist on disk
    allfiles = os.listdir(os.path.dirname(os.path.abspath(imagename)))
    extims = list(filter(lambda im: os.path.exists(im), reqims)) # TODO replace with allfiles

    # verify required images are available
    if len(extims) != len(reqims):
        diffims = list(filter(lambda im: im not in extims, reqims))
        raise RuntimeError("Internal Error: missing one or more of the required images: " + str(diffims))

    # check for .pb image in the case that nsigma > 0
    # see comments on CAS-13144 about casa crashing as to why this check is here
    if (inp['nsigma'] > 0):
        reqpb = ".pb" if (inp['deconvolver'] != 'mtmfs') else ".pb.tt0"
        if (imagename+reqpb not in allfiles):
            raise RuntimeError("The parameter nsigma>0 ("+str(inp['nsigma'])+") requires a "+reqpb+" image to be available.")

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
    interactive,#=False,
    fastnoise,#=True,

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
    verbose): #=False
    """
    Runs the minor cycle only of tclean.
    Most of this code is copied directly from tclean.
    """

    cppparallel=False
    decon=None
    try:

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
        #### Construct ImagerParameters
        #####################################################

        # Make sure that we have all the necessary images and that the startmodel is valid.
        # Note: cpp code should check that .residual and .psf exist, but current tests indicate that it doesn't do that.
        check_requiredmask_exists(usemask, mask)
        check_requiredimgs_exist(imagename, inp)
        check_starmodel_model_collisions(startmodel, imagename, deconvolver)
        
        # make a list of parameters with defaults from tclean
        if is_python3:
            defparm=dict(list(zip(ImagerParameters.__init__.__code__.co_varnames[1:], ImagerParameters.__init__.__defaults__)))
        else:
            defparm=dict(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:], ImagerParameters.__init__.func_defaults))

        ## assign values to the ones passed to deconvolve and if not defined yet in deconvolve...
        ## assign them the default value of the constructor
        bparm={k:  inp[k] if k in inp else defparm[k]  for k in defparm.keys()}

        #=========================================================
        ####set the children to load c++ libraries and applicator
        ### make workers ready for c++ based mpicommands
        if deconvolver != 'mtmfs': # mtmfs isn't currently parallelized
            ia.open(imagename+'.psf')
            isCube=ia.shape()[3] >1 # SynthesisDeconvolver::executeMinorCycle doesn't start mpi without more than one channel
            ia.done()
            if mpi_available and MPIEnvironment.is_mpi_enabled and isCube:
                mint=MPIInterface.MPIInterface()
                cl=mint.getCluster()
                if(is_CASA6):
                    cl._cluster.pgc("from casatools import synthesisimager", False)
                    cl._cluster.pgc("si=synthesisimager()", False)
                else:
                    cl._cluster.pgc("from casac import casac", False)
                    cl._cluster.pgc("si=casac.synthesisimager()", False) 
                cl._cluster.pgc("si.initmpi()", False)
                cppparallel=True
                ###ignore chanchunk
                bparm['chanchunks']=1

        ## create the parameters list help object
        paramList=ImagerParameters(**bparm)

        # Assign cyclethreshold explicitly to threshold
        threshold = threshold if (type(threshold) == str) else (str(threshold*1000)+'mJy')
        paramList.setIterPars({'cyclethreshold': threshold, 'cyclethresholdismutable': False})

        #####################################################
        #### Run the minor cycle
        #####################################################

        iterrec = False
        isit = 0
        retrec = ''

        ## Setup Imager object
        decon = PyDeconvolver(params=paramList)

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
        larg = list(e.args)
        larg[0] = 'Exception from deconvolve : ' + str(larg[0])
        e.args = tuple(larg)
        raise

    finally:
        if decon != None:
            decon.deleteTools()
        if(cppparallel):
            ###release workers back to python mpi control
            si=synthesisimager()
            si.releasempi()

    # Write history at the end, when hopefully all temp files are gone from disk,
    # so they won't be picked up. They need time to disappear on NFS or slow hw.
    # Copied from tclean.
    try:
        params = get_func_params(deconvolve, locals())
        write_tclean_history(imagename, 'deconvolve', params, casalog)
    except Exception as exc:
        casalog.post("Error updating history (logtable): {} ".format(exc),'WARN')

    return { 'iterrec': iterrec, 'isit': isit, 'retrec': retrec }
################################################
# Refactored Clean task
#
# v1.0: 2012.10.05, U.R.V.
#
################################################

from __future__ import absolute_import
from __future__ import print_function

import os
import shutil
import numpy
import copy
import time
# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatasks import casalog

    from casatasks.private.imagerhelpers.imager_base import PySynthesisImager
    from casatasks.private.imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from casatasks.private.imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
    from casatools import table
    from casatools import synthesisimager
else:
    from taskinit import *

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from imagerhelpers.input_parameters import ImagerParameters
    table=casac.table
    synthesisimager=casac.synthesisimager
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

def tclean(
    ####### Data Selection
    vis,#='', 
    selectdata,
    field,#='', 
    spw,#='',
    timerange,#='',
    uvrange,#='',
    antenna,#='',
    scan,#='',
    observation,#='',
    intent,#='',
    datacolumn,#='corrected',

    ####### Image definition
    imagename,#='',
    imsize,#=[100,100],
    cell,#=['1.0arcsec','1.0arcsec'],
    phasecenter,#='J2000 19:59:28.500 +40.44.01.50',
    stokes,#='I',
    projection,#='SIN',
    startmodel,#='',

    ## Spectral parameters
    specmode,#='mfs',
    reffreq,#='',
    nchan,#=1,
    start,#='',
    width,#='',
    outframe,#='LSRK',
    veltype,#='',
    restfreq,#=[''],
#    sysvel,#='',
#    sysvelframe,#='',
    interpolation,#='',
    perchanweightdensity, #=''
    ## 
    ####### Gridding parameters
    gridder,#='ft',
    facets,#=1,
    psfphasecenter,#='',
    chanchunks,#=1,

    wprojplanes,#=1,

    ### PB
    vptable,
    mosweight, #=True
    aterm,#=True,
    psterm,#=True,
    wbawp ,#= True,
    conjbeams ,#= True,
    cfcache ,#= "",
    usepointing, #=false
    computepastep ,#=360.0,
    rotatepastep ,#=360.0,
    pointingoffsetsigdev ,#=[10.0],

    pblimit,#=0.01,
    normtype,#='flatnoise',

    ####### Deconvolution parameters
    deconvolver,#='hogbom',
    scales,#=[],
    nterms,#=1,
    smallscalebias,#=0.0

    ### restoration options
    restoration,
    restoringbeam,#=[],
    pbcor,

    ##### Outliers
    outlierfile,#='',

    ##### Weighting
    weighting,#='natural',
    robust,#=0.5,
    noise,#0.0Jy
    npixels,#=0,
#    uvtaper,#=False,
    uvtaper,#=[],


    ##### Iteration control
    niter,#=0, 
    gain,#=0.1,
    threshold,#=0.0, 
    nsigma,#=0.0
    cycleniter,#=0, 
    cyclefactor,#=1.0,
    minpsffraction,#=0.1,
    maxpsffraction,#=0.8,
    interactive,#=False, 

    ##### (new) Mask parameters
    usemask,#='user',
    mask,#='',
    pbmask,#='',
    # maskthreshold,#='',
    # maskresolution,#='',
    # nmask,#=0,

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
    minpercentchange,#=0.0
    verbose, #=False
    fastnoise, #=False

    ## Misc

    restart,#=True,

    savemodel,#="none",

#    makeimages,#="auto"
    calcres,#=True,
    calcpsf,#=True,

    ####### State parameters
    parallel):#=False):

    #####################################################
    #### Sanity checks and controls
    #####################################################
    
    ### Move these checks elsewhere ? 
    inpparams=locals().copy()
    ###now deal with parameters which are not the same name 
    inpparams['msname']= inpparams.pop('vis')
    inpparams['timestr']= inpparams.pop('timerange')
    inpparams['uvdist']= inpparams.pop('uvrange')
    inpparams['obs']= inpparams.pop('observation')
    inpparams['state']= inpparams.pop('intent')
    inpparams['loopgain']=inpparams.pop('gain')
    inpparams['scalebias']=inpparams.pop('smallscalebias')

    if specmode=='cont':
        specmode='mfs'
        inpparams['specmode']='mfs'
#    if specmode=='mfs' and nterms==1 and deconvolver == "mtmfs":
#        casalog.post( "The MTMFS deconvolution algorithm (deconvolver='mtmfs') needs nterms>1.Please set nterms=2 (or more). ", "WARN", "task_tclean" )
#        return

    if(deconvolver=="mtmfs" and (specmode=='cube' or specmode=='cubedata')):
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') with specmode='cube' is not supported", "WARN", "task_tclean" )
        return

    if(chanchunks!=-1):
        casalog.post( "The parameter chanchunks is no longer used by tclean. It will be removed in CASA 6.2", "WARN", "task_tclean" )

    if((specmode=='cube' or specmode=='cubedata') and parallel==False and mpi_available):
        casalog.post( "Setting parameter parallel=False with specmode='cube' when launching CASA with mpi has no effect.", "WARN", "task_tclean" )
        
    if((specmode=='cube' or specmode=='cubedata') and gridder=='awproject'):
        casalog.post( "The awproject gridder does not currently work with specmode='cube'.", "WARN", "task_tclean" )
        return
      

    #####################################################
    #### Construct ImagerParameters object
    #####################################################

    imager = None
    paramList = None

    # Put all parameters into dictionaries and check them.
    ##make a dictionary of parameters that ImagerParameters take

    if is_python3:
        defparm=dict(list(zip(ImagerParameters.__init__.__code__.co_varnames[1:], ImagerParameters.__init__.__defaults__)))
    else:
        defparm=dict(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:], ImagerParameters.__init__.func_defaults))
        
    ###assign values to the ones passed to tclean and if not defined yet in tclean...
    ###assign them the default value of the constructor
    bparm={k:  inpparams[k] if k in inpparams else defparm[k]  for k in defparm.keys()}

    ###default mosweight=True is tripping other gridders as they are not
    ###expecting it to be true
    if(bparm['mosweight']==True and bparm['gridder'].find("mosaic") == -1):
        bparm['mosweight']=False

    
    # deprecation message
    if usemask=='auto-thresh' or usemask=='auto-thresh2':
        casalog.post(usemask+" is deprecated, will be removed in CASA 5.4.  It is recommended to use auto-multithresh instead", "WARN") 

    #paramList.printParameters()
    
    if len(pointingoffsetsigdev)>0 and pointingoffsetsigdev[0]!=0.0 and usepointing==True and gridder.count('awproj')>1:
        casalog.post("pointingoffsetsigdev will be used for pointing corrections with AWProjection", "WARN") 
#    elif usepointing==True and pointingoffsetsigdev[0] == 0:
#        casalog.post("pointingoffsetsigdev is set to zero which is an unphysical value, will proceed with the native sky pixel resolution instead". "WARN")


    ##pcube may still need to be set to True for some combination of ftmachine etc...
    #=========================================================
    concattype=''
    pcube=False
    if parallel==True and specmode!='mfs':
        pcube=False
        parallel=False
    #=========================================================
    ####set the children to load c++ libraries and applicator
    ### make workers ready for c++ based mpicommands
    cppparallel=False
    if mpi_available and MPIEnvironment.is_mpi_enabled and specmode!='mfs' and not pcube:
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

    # catch non operational case (parallel cube tclean with interative=T)
    if pcube and interactive:
        casalog.post( "Interactive mode is not currently supported with parallel cube CLEANing, please restart by setting interactive=F", "WARN", "task_tclean" )
        return False
    #print('parameters {}'.format(bparm))    
    paramList=ImagerParameters(**bparm)

    ## Setup Imager objects, for different parallelization schemes.
    imagerInst=PySynthesisImager
    if parallel==False and pcube==False:
         imager = PySynthesisImager(params=paramList)
         imagerInst=PySynthesisImager
    elif parallel==True:
         imager = PyParallelContSynthesisImager(params=paramList)
         imagerInst=PySynthesisImager
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params=paramList)
         imagerInst=PyParallelCubeSynthesisImager
         # virtualconcat type - changed from virtualmove to virtualcopy 2016-07-20
         #using ia.imageconcat now the name changed to copyvirtual 2019-08-12
         concattype='copyvirtual'
    else:
         print('Invalid parallel combination in doClean.')
         return False
    
    retrec={}

    try: 
    #if (1):
        #pdb.set_trace()
        ## Init major cycle elements
        t0=time.time();
        imager.initializeImagers()
    
        # Construct the CFCache for AWProject-class of FTMs.  For
        # other choices the following three calls become NoOps.
        # imager.dryGridding();
        # imager.fillCFCache();
        # imager.reloadCFCache();

        imager.initializeNormalizers()
        imager.setWeighting()
        t1=time.time();
        casalog.post("***Time for initializing imager and normalizers: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        ## Init minor cycle elements
        if niter>0 or restoration==True:
            t0=time.time();
            imager.initializeDeconvolvers()
            t1=time.time();
            casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        ####now is the time to check estimated memory
        imager.estimatememory()
            
        if niter>0:
            t0=time.time();
            imager.initializeIterationControl()
            t1=time.time();
            casalog.post("***Time for initializing iteration controller: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
        
        ## Make PSF
        if calcpsf==True:
            t0=time.time();
             
            imager.makePSF()
            if((psfphasecenter != '') and ('mosaic' in gridder)):
                ###for some reason imager keeps the psf open delete it and recreate it afterwards
                imager.deleteTools()
                mytb=table()
                psfname=bparm['imagename']+'.psf.tt0' if(os.path.exists( bparm['imagename']+'.psf.tt0')) else bparm['imagename']+'.psf'
                mytb.open(psfname)
                miscinf=mytb.getkeyword('miscinfo')
                iminf=mytb.getkeyword('imageinfo')
                #print ('miscinfo {} {}'.format(miscinf, iminf))
                mytb.done()
                print("doing with different phasecenter psf")
                imager.unlockimages(0)
                psfParameters=paramList.getAllPars()
                psfParameters['phasecenter']=psfphasecenter
                psfParamList=ImagerParameters(**psfParameters)
                psfimager=imagerInst(params=psfParamList)
                psfimager.initializeImagers()
                psfimager.setWeighting()
                psfimager.makeImage('psf', psfParameters['imagename']+'.psf')
                psfimager.deleteTools()
                mytb.open(psfname, nomodify=False)
                mytb.putkeyword('imageinfo',iminf)
                mytb.putkeyword('miscinfo',miscinf)
                mytb.done()
                imager = PySynthesisImager(params=paramList)
                imager.initializeImagers()
                imager.initializeNormalizers()
                imager.setWeighting()
                ###redo these as we destroyed things for lock issues
                ## Init minor cycle elements
                if niter>0 or restoration==True:
                    imager.initializeDeconvolvers() 
                if niter>0:
                    imager.initializeIterationControl()

            t1=time.time();
            casalog.post("***Time for making PSF: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
        

            imager.makePB()

            t2=time.time();
            casalog.post("***Time for making PB: "+"%.2f"%(t2-t1)+" sec", "INFO3", "task_tclean");

        if niter >=0 : 

            ## Make dirty image
            if calcres==True:
                t0=time.time();
                imager.runMajorCycle()
                t1=time.time();
                casalog.post("***Time for major cycle (calcres=T): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean"); 

            ## In case of no deconvolution iterations....
            if niter==0 and calcres==False:
                if savemodel != "none":
                    imager.predictModel()
            ## Do deconvolution and iterations
            if niter>0 :

                isit = imager.hasConverged()
                t0=time.time();
                imager.updateMask()
                t1=time.time();
                casalog.post("***Time to update mask: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                while ( not isit ):
                    t0=time.time();
                    imager.runMinorCycle()
                    t1=time.time();
                    casalog.post("***Time for minor cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

                    t0=time.time();
                    imager.runMajorCycle()
                    t1=time.time();
                    casalog.post("***Time for major cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                   
                    imager.updateMask()
                    t2=time.time()
                    casalog.post("***Time to update mask: "+"%.2f"%(t2-t1)+" sec", "INFO3", "task_tclean");
                    isit = imager.hasConverged()
                ## Get summary from iterbot
                if type(interactive) != bool:
                    retrec=imager.getSummary();

            ## Restore images.
            if restoration==True:  
                t0=time.time();
                imager.restoreImages()
                t1=time.time();
                casalog.post("***Time for restoring images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                if pbcor==True:
                    t0=time.time();
                    imager.pbcorImages()
                    t1=time.time();
                    casalog.post("***Time for pb-correcting images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                    
        ##close tools
        # needs to deletools before concat or lock waits for ever


        imager.deleteTools()
        if(cppparallel):
            ###release workers back to python mpi control
            si=synthesisimager()
            si.releasempi()
            
        if (pcube):
            print("running concatImages ...")
            casalog.post("Running virtualconcat (type=%s) of sub-cubes" % concattype,"INFO2", "task_tclean")
            imager.concatImages(type=concattype)
        # CAS-10721 
        if niter>0 and savemodel != "none":
            casalog.post("Please check the casa log file for a message confirming that the model was saved after the last major cycle. If it doesn't exist, please re-run tclean with niter=0,calcres=False,calcpsf=False in order to trigger a 'predict model' step that obeys the savemodel parameter.","WARN","task_tclean")


    except Exception as e:
        #print 'Exception : ' + str(e)
        if(cppparallel):
            ###release workers back to python mpi control
            si=synthesisimager()
            si.releasempi()
        casalog.post('Exception from task_tclean : ' + str(e), "SEVERE", "task_tclean")
        if imager != None:
            imager.deleteTools()
            

        larg = list(e.args)
        larg[0] = 'Exception from task_tclean : ' + str(larg[0])
        e.args = tuple(larg)
        raise

    return retrec

##################################################

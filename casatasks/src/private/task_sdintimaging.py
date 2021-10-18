################################################
# single dish + interfermeter join image reconstruction task
#
#
################################################

from __future__ import absolute_import

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
    #from casatasks import imregrid
    from .cleanhelper import write_tclean_history, get_func_params
    from .sdint_helper import *
    from casatools import table
    from casatools import synthesisimager,synthesisutils
else:
    from taskinit import *
    from tasks import *

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from imagerhelpers.input_parameters import ImagerParameters
    from cleanhelper import write_tclean_history, get_func_params
    from sdint_helper import *
    table=casac.table
    synthesisimager=casac.synthesisimager
    synthesisutils=casac.synthesisutils

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

sdintlib = SDINT_helper()
synu = synthesisutils()

# setup functions
def setup_imagerObj(paramList=None):
    """
    setup imaging parameters
    """
    defaultconstructor = False
    if paramList!=None:
        if not isinstance(paramList, ImagerParameters):
            raise RuntimeError("Internal Error: invalid paramList")
    else:
       defaultconstructor = True
       
    if defaultconstructor:
        return PySynthesisImager
    else:
        return PySynthesisImager(params=paramList)


def setup_imager(imagename, specmode,calcres,calpsf,inparams):
    """
     Setup cube imaging for major cycles.
     - Do initialization
     - and run a major cycle
    """
    # create a local copy of input params dict so that it can be 
    # modified
    locparams = copy.deepcopy(inparams)

    # cube imaging setup 
    locparams['imagename']=imagename
    locparams['specmode']='cube'
    locparams['niter']=0
    locparams['deconvolver']='hogbom'

    #casalog.post("local inparams(msname) in setup_imager==",locparams['msname'])
    params = ImagerParameters(**locparams)
    #params = ImagerParameters(msname=self.vis, field=self.field,spw=self.spw,
    #                              imagename=imagename,
    #                              imsize=self.imsize, cell=self.cell, phasecenter=self.phasecenter,
    #                              weighting=self.weighting,
    #                              gridder=self.gridder, pblimit=self.pblimit,wprojplanes=self.wprojplanes,
    #                              specmode='cube',nchan=self.nchan,
    #                              reffreq=self.reffreq, width=self.width, start=self.start, interpolation=self.interpolation,
    #                              interactive=self.interactive,
    #                              deconvolver='hogbom', niter=0,
    #                              wbawp=True)

    ## Major cycle is either PySynthesisImager or PyParallelCubeSynthesisImager
    imagertool = setup_imagerObj(params)

    #self.imagertool = PySynthesisImager(params=params)
    imagertool.initializeImagers()
    imagertool.initializeNormalizers()
    imagertool.setWeighting()
    if 'psfphasecenter' in  locparams:
        psfphasecenter = locparams['psfphasecenter']
    else:
        psfphasecenter = ''

    ## Extra one for psfphasecenter...
    imagerInst=None
    if((psfphasecenter != '') and (gridder=='mosaic')):
        imagerInst = setup_imagerObj()

  
    gridder = locparams['gridder']

    if calpsf == True:
        imagertool.makePSF()
        imagertool.makePB()
        if((psfphasecenter != '') and (gridder=='mosaic')):
            casalog.post("doing with different phasecenter psf", "INFO")
            imagertool.unlockimages(0)
            psfParameters=paramList.getAllPars()
            psfParameters['phasecenter']=psfphasecenter
            psfParamList=ImagerParameters(**psfParameters)
            psfimager=imagerInst(params=psfParamList)
            psfimager.initializeImagers()
            psfimager.setWeighting()
            psfimager.makeImage('psf', psfParameters['imagename']+'.psf')

    # can take out this since niter is fixed to 0
    if locparams['niter'] >=0 :
        ## Make dirty image
        if calcres == True:
            t0=time.time();
            imagertool.runMajorCycle()
            t1=time.time();
            casalog.post("***Time for major cycle (calcres=T): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        ## In case of no deconvolution iterations....
        #if locparams['niter']==0 and calcres==False:
        #    if savemodel != "none":
        #        imagertool.predictModel()

    sdintlib.copy_restoringbeam(fromthis=imagename+'.psf', tothis=imagename+'.residual')
    return imagertool

def setup_deconvolver(imagename,specmode,inparams):
    """
    Cube or MFS minor cycles. 
    """
    #params = ImagerParameters(msname=self.vis, field=self.field,spw=self.spw,
    #                          imagename=imagename,
    #                            imsize=self.imsize, cell=self.cell, phasecenter=self.phasecenter,
    #                              weighting=self.weighting,
    #                              gridder=self.gridder, pblimit=self.pblimit,wprojplanes=self.wprojplanes,
    #                              specmode=self.specmode,nchan=self.nchan,
    #                              reffreq=self.reffreq, width=self.width, start=self.start, interpolation=self.interpolation,
    #                              interactive=self.interactive,
    #                              deconvolver=self.deconvolver, scales=self.scales,nterms=self.nterms,
    #                              niter=self.niter,cycleniter=self.cycleniter, threshold=self.threshold,
    #                              mask=self.mask)
    inparams['imagename']=imagename
    params = ImagerParameters(**inparams)
    deconvolvertool = setup_imagerObj(params)

    ## Why are we initializing these ? 
    deconvolvertool.initializeImagers()
    deconvolvertool.initializeNormalizers()
    deconvolvertool.setWeighting()


        ### These three should be unncessary.  Need a 'makeimage' method for csys generation. 
    deconvolvertool.makePSF() ## Make this to get a coordinate system
    #deconvolvertool.makeImage('psf', imagename+'.psf')
    deconvolvertool.makePB()  ## Make this to turn .weight into .pb maps

        ## Initialize deconvolvers. ( Order is important. This cleans up a leftover tablecache image.... FIX!)
    deconvolvertool.initializeDeconvolvers()
    deconvolvertool.initializeIterationControl() # This needs to be run before runMajorCycle
    deconvolvertool.runMajorCycle() ## Make this to make template residual images.
 
    return deconvolvertool

def setup_sdimaging(template='',output='', inparms=None, sdparms=None):
    """
    Make the SD cube Image and PSF

    Option 1 : Use/Regrid cubes for the observed image and PSF
    Option 2 : Make the SD image and PSF cubes using 'tsdimager's usage of the SD gridder option.

    Currently, only Option 1 is supported. 

    """
    sdintlib = SDINT_helper()
    if 'sdpsf' in sdparms:
        sdpsf = sdparms['sdpsf']
    else:
        raise RuntimeError("Internal Error: missing sdpsf parameter") 

    if 'sdimage' in sdparms:
        sdimage = sdparms['sdimage']
    else:
        raise RuntimeError("Internal Error: missing sdimage parameter") 
    if 'pblimit' in inparms:
        pblimit = inparms['pblimit']

    if sdpsf !="":
        ## check the coordinates of psf with int psf
        sdintlib.checkpsf(sdpsf, template+'.psf') 

    ## Regrid the input SD image and PSF cubes to the target coordinate system. 
    #imregrid(imagename=sdimage, template=template+'.residual',
    #         output=output+'.residual',overwrite=True,axes=[0,1])
    sdintlib.regridimage(imagename=sdimage, template=template+'.residual', outfile=output+'.residual')
    #imregrid(imagename=sdimage, template=template+'.residual',
    #         output=output+'.image',overwrite=True,axes=[0,1])
    sdintlib.regridimage(imagename=sdimage, template=template+'.residual', outfile=output+'.image')

    if sdpsf !="":
        #imregrid(imagename=sdpsf, template=template+'.psf',
        #         output=output+'.psf',overwrite=True,axes=[0,1])
        sdintlib.regridimage(imagename=sdpsf, template=template+'.psf', outfile=output+'.psf')
    else:
        ## Make an internal sdpsf image if the user has not supplied one. 
        casalog.post("Constructing a SD PSF cube by evaluating Gaussians based on the restoring beam information in the regridded SD Image Cube")
        sdintlib.create_sd_psf(sdimage=output+'.residual', sdpsfname=output+'.psf')

    ## Apply the pbmask from the INT image cube, to the SD cubes.
    #TTB: Create *.mask cube  

    sdintlib.addmask(inpimage=output+'.residual', pbimage=template+'.pb', pblimit=pblimit)
    sdintlib.addmask(inpimage=output+'.image', pbimage=template+'.pb', pblimit=pblimit)



def feather_residual(int_cube, sd_cube, joint_cube, applypb, inparm):

    if applypb==True:
        ## Take initial INT_dirty image to flat-sky. 
        sdintlib.modify_with_pb(inpcube=int_cube+'.residual',
                                    pbcube=int_cube+'.pb',
                                    cubewt=int_cube+'.sumwt',
                                chanwt=inparm['chanwt'],
                                    action='div',
                                    pblimit=inparm['pblimit'],
                                    freqdep=True)
            
        ## Feather flat-sky INT dirty image with SD image
    sdintlib.feather_int_sd(sdcube=sd_cube+'.residual',
                                intcube=int_cube+'.residual',
                                jointcube=joint_cube+'.residual', ## output
                                sdgain=inparm['sdgain'],
                                dishdia=inparm['dishdia'],
                                usedata=inparm['usedata'],
                                chanwt = inparm['chanwt'])

    if applypb==True:
        if inparm['specmode'].count('cube')>0:
            ## Multiply the new JOINT dirty image by the frequency-dependent PB. 
            fdep_pb = True
        else:
            ## Multiply new JOINT dirty image by a common PB to get the effect of conjbeams. 
            fdep_pb = False
        sdintlib.modify_with_pb(inpcube=joint_cube+'.residual',
                                pbcube=int_cube+'.pb',
                                cubewt=int_cube+'.sumwt',
                                chanwt=inparm['chanwt'],
                                action='mult',
                                pblimit=inparm['pblimit'],
                                freqdep=fdep_pb)

def deleteTmpFiles():
    if os.path.exists('tmp_intplane'):
       os.system('rm -rf tmp_intplane')
    if os.path.exists('tmp_sdplane'):
       os.system('rm -rf tmp_sdplane')
    if os.path.exists('tmp_jointplane'):
       os.system('rm -rf tmp_jointplane')


def sdintimaging(
    usedata,
    ####### Single dish input data
    sdimage, 
    sdpsf, 
    sdgain, 
    dishdia,
    ####### Interfermeter Data Selection
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
#    chanchunks,#=1,
    perchanweightdensity, #=''
    ## 
    ####### Gridding parameters
    gridder,#='ft',
    facets,#=1,
    psfphasecenter,#='',

    wprojplanes,#=1,

    ### PB
    vptable,
    mosweight, #=True
    aterm,#=True,
    psterm,#=True,
    wbawp ,#= True,
#    conjbeams ,#= True,
    cfcache ,#= "",
    usepointing, #=false
    computepastep ,#=360.0,
    rotatepastep ,#=360.0,
    pointingoffsetsigdev ,#=0.0,

    pblimit,#=0.01,
#    normtype,#='flatnoise',

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
#    outlierfile,#='',    ### RESTRICTION : No support for outlier fields for joint SD-INT imaging. 

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

    #savemodel,#="none",

#    makeimages,#="auto"
    calcres,#=True,
    calcpsf):#=True,

    ####### State parameters
    #parallel):#=False):


    ##################################################
    # copied from SDINT.do_reconstruct 
    #################################################
    int_cube = imagename+'.int.cube'
    sd_cube = imagename+'.sd.cube'
    joint_cube = imagename+'.joint.cube'
    joint_multiterm = imagename+'.joint.multiterm'

    if specmode=='mfs':
        decname = joint_multiterm
    else:
        decname = joint_cube

    #####################################################
    #### Sanity checks and controls
    #####################################################
    
    ### Move these checks elsewhere ? 
    inpparams=locals().copy()
    ###now deal with parameters which are not the same name 
    #casalog.post("current inpparams=",inpparams)
    #casalog.post("inpparams.keys()=",inpparams.keys())
    locvis=inpparams.pop('vis')
    #casalog.post("LOCVIS====",locvis)
    if type(locvis)==list:
        llocvis = [v.lstrip() for v in locvis]
    else:
        llocvis = locvis.lstrip()
    inpparams['msname']=llocvis
    inpparams['timestr']= inpparams.pop('timerange')
    inpparams['uvdist']= inpparams.pop('uvrange')
    inpparams['obs']= inpparams.pop('observation')
    inpparams['state']= inpparams.pop('intent')
    inpparams['loopgain']=inpparams.pop('gain')
    inpparams['scalebias']=inpparams.pop('smallscalebias')

    sdparms={}
    sdparms['sdimage']=inpparams['sdimage']
    sdparms['sdpsf']=inpparams['sdpsf']
    sdparms['sdgain']=inpparams['sdgain']

    if specmode=='cont':
        specmode='mfs'
        inpparams['specmode']='mfs'

    # from sdint
    # automatically decide if pb need to be applied
    if gridder=='mosaic' or gridder=='awproject':
       applypb = True
    else:
       applypb = False
   
    if (deconvolver=="mtmfs") and (specmode!='mfs') and (specmode!='cube' or nterms!=1) and (specmode!='cubedata' or nterms!=1):
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') applies only to specmode='mfs' or specmode='cube' with nterms=1 or specmode='cubedata' with nterms=1.", "WARN", "task_sdintimaging" )
        return
      
    if(deconvolver=="mtmfs" and (specmode=='cube' or specmode=='cubedata') and nterms==1 ):
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') with specmode='cube', nterms=1 is currently not supported. Please use deconvolver='multiscale' instead for cubes.", "WARN", "task_sdintimaging" )
        return

    if(specmode=='mfs' and deconvolver!='mtmfs'):
        casalog.post("Currently, only the multi-term MFS algorithm is supported for specmode=mfs. To make a single plane MFS image (while retaining the frequency dependence for the cube major cycle stage), please pick nterms=1 along with deconvolver=mtmfs. The scales parameter is still usable for multi-scale multi-term deconvolution","WARN","task_sdintimaging")
        return;
        
    if(gridder=='awproject'):
        casalog.post("The awproject gridder is temporarily not supported with cube major cycles. Support will be brought back in a subsequent release.","WARN","task_sdintimaging")
        return;

    if(usedata=='sd'):
        casalog.post("The Single-Dish-Only mode of sdintimaging is better supported via the deconvolve task which supports spectral cube, mfs and multi-term mfs deconvolution in the image domain alone. The deconvolve task is the more appropriate version to use for stand-alone image-domain deconvolution, and will not have the bookkeeping overheads currently present in the sdintimaging task's sd-only mode. Please note that the 'sd' option of the sdintimaging task will be removed in a subsequent release.  Please refer to the task deconvolve documentation for instructions on how to prepare image and psf cubes for the deconvolve task for all these modes.","WARN","task_sdintimaging");


#    if parallel==True:
#        casalog.post("Cube parallelization (all major cycles) is currently not supported via task_sdintimaging. This will be enabled after a cube parallelization rework.")
#        return;

    #####################################################
    #### Construct ImagerParameters object
    #####################################################

    imager = None
    paramList = None
    deconvolvertool = None

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

    ## Two options have been removed from the interface. Hard-code them here.
    bparm['normtype'] = 'flatnoise'  ## Hard-code this since the pbcor steps assume it.
    bparm['conjbeams']=False

    #paramList=ImagerParameters(**bparm)

    #paramList.printParameters()
    
    if len(pointingoffsetsigdev)>0 and pointingoffsetsigdev[0]!=0.0 and usepointing==True and gridder.count('awproj')>1:
        casalog.post("pointingoffsetsigdev will be used for pointing corrections with AWProjection", "WARN") 

    #=========================================================
    ####set the children to load c++ libraries and applicator
    ### make workers ready for c++ based mpicommands
    cppparallel=False
    if mpi_available and MPIEnvironment.is_mpi_enabled:
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

    
    retrec={}

    try: 
        sdintlib = SDINT_helper()
        ## Init major cycle elements
        ### debug (remove it later) 
        casalog.post("INT cube setup ....")
        t0=time.time();
        imager=setup_imager(int_cube, specmode, calcres, calcpsf, bparm) 


        ##imager.initializeImagers()
    
        ##imager.initializeNormalizers()
        ##imager.setWeighting()
        t1=time.time();
        casalog.post("***Time for initializing imager (INT cube) : "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_sdintimaging");

        ## Init minor cycle elements
        if niter>0 or restoration==True:
            ### debug (remove it later) 
            casalog.post("Combined image setup ....")
            t0=time.time();
            deconvolvertool=setup_deconvolver(decname, specmode, bparm )
            #imager.initializeDeconvolvers()
            t1=time.time();
            #casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
            casalog.post("***Time for seting up deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_sdintimaging");

        if usedata!='int':
            ### debug (remove it later) 
            casalog.post("SD cube setup ....")
            setup_sdimaging(template=int_cube, output=sd_cube, inparms=bparm, sdparms=sdparms ) 
            

        ####now is the time to check estimated memory
        # need to move to somewhere below???
        imager.estimatememory()

        ## Do sanity checks on INT and SD cubes
        ### Frequency range of cube, data selection range. mtmfs reffreq.
        ### nchan too small or too large
        ### sumwt : flagged channels in int cubes
        ### sd cube empty channels ? Weight image ? 
        validity, inpparams = sdintlib.check_coords(intres=int_cube+'.residual', intpsf=int_cube+'.psf', 
                                         intwt = int_cube+'.sumwt', 
                                         sdres=sd_cube+'.residual', sdpsf=sd_cube+'.psf',
                                         sdwt = '',
                                         pars=inpparams)

        if validity==False:
            casalog.post('Exiting from the sdintimaging task due to inconsistencies found between the interferometer-only and singledish-only image and psf cubes. Please modify inputs as needed','WARN')
            if imager != None:
                imager.deleteTools()
            if deconvolvertool != None:
                deconvolvertool.deleteTools()
            deleteTmpFiles()
            return

        #### inpparams now has a new parameter "chanwt" with ones and zeros to indicate chans
        ####  that have data from both INT and SD cubes (they are the 'ones'). This is to be used in 
        ####  feathering and in the cube-to-taylor sum and modify_with_pb. 

        #### END MOVE THIS SECTION to setup_imager ....for sdint

        #### SDINT specific feathering....
        ## Feather INT and SD residual images (feather in flat-sky. output has common PB)
        casalog.post("Feathering INT and SD residual images...")
        feather_residual(int_cube, sd_cube, joint_cube, applypb, inpparams)
        sdintlib.feather_int_sd(sdcube=sd_cube+'.psf',
                                intcube=int_cube+'.psf',
                                jointcube=joint_cube+'.psf',
                                sdgain=sdgain,
                                dishdia=dishdia,
                                usedata=usedata,
                                chanwt = inpparams['chanwt'])

        #print("Fitting for cube")
        synu.fitPsfBeam(joint_cube)

        ###############
        ##### Placeholder code for PSF renormalization if needed
        #####  Note : If this is enabled, we'll need to restrict the use of 'faceting' as .sumwt shape changes.
        #sdintlib.calc_renorm(intname=int_cube, jointname=joint_cube)
        #sdintlib.apply_renorm(imname=joint_cube+'.psf', sumwtname=joint_cube+'.sumwt')
        #sdintlib.apply_renorm(imname=joint_cube+'.residual', sumwtname=joint_cube+'.sumwt')
        ###############

        #casalog.post("feather_int_sd DONE")
 
        if specmode=='mfs':
            ## Calculate Spectral PSFs and Taylor Residuals
            casalog.post("Calculate spectral PSFs and Taylor Residuals...")
            sdintlib.cube_to_taylor_sum(cubename=joint_cube+'.psf',
                                        cubewt=int_cube+'.sumwt',
                                        chanwt=inpparams['chanwt'],
                                        mtname=joint_multiterm+'.psf',
                                        nterms=nterms, reffreq=inpparams['reffreq'], dopsf=True)
            sdintlib.cube_to_taylor_sum(cubename=joint_cube+'.residual',
                                        cubewt=int_cube+'.sumwt',
                                        chanwt=inpparams['chanwt'],
                                        mtname=joint_multiterm+'.residual',
                                        nterms=nterms, reffreq=inpparams['reffreq'], dopsf=False)

            #print("Fit for multiterm")
            synu.fitPsfBeam(joint_multiterm,nterms=nterms)

        if niter>0 :
            isit = deconvolvertool.hasConverged()
            deconvolvertool.updateMask()

            while ( not deconvolvertool.hasConverged() ):
 
                t0=time.time();

                #### Print out the peak of the residual image here to check !!! 
#                if specmode=='mfs':
#                    print('Max of joint residual before initminorcycle' + str(imstat(joint_multiterm+'.residual.tt0',verbose=False)['max'][0]))
#                else:
#                    print('Max of joint residual before initminorcycle' + str(imstat(joint_cube+'.residual',verbose=False)['max'][0]))



                deconvolvertool.runMinorCycle()

#                if specmode=='mfs':
#                    print('Max of joint residual after minorcycle' + str(imstat(joint_multiterm+'.residual.tt0',verbose=False)['max'][0]))
#                else:
#                    print('Max of joint residual after minorcycle' + str(imstat(joint_cube+'.residual',verbose=False)['max'][0]))


                t1=time.time();
                casalog.post("***Time for minor cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_sdintimaging");

                ### sdint specific feathering steps HERE
                ## Prepare the joint model cube for INT and SD major cycles
                if specmode=='mfs':
                    ## Convert Taylor model coefficients into a model cube : int_cube.model
                    sdintlib.taylor_model_to_cube(cubename=int_cube, ## output 
                                              mtname=joint_multiterm,  ## input
                                              nterms=nterms, reffreq=inpparams['reffreq'])
                else:
                    ## Copy the joint_model cube to the int_cube.model
                    shutil.rmtree(int_cube+'.model',ignore_errors=True)
                    shutil.copytree(joint_cube+'.model', int_cube+'.model')
                    hasfile=os.path.exists(joint_cube+'.model')
                    #casalog.post("DEBUG: has joint cube .image===",hasfile)

                if applypb==True:
                    ## Take the int_cube.model to flat sky. 
                    if specmode=='cube':
                        ## Divide the model by the frequency-dependent PB to get to flat-sky
                        fdep_pb = True
                    else:
                        ## Divide the model by the common PB to get to get to flat-sky
                        fdep_pb = False
                    sdintlib.modify_with_pb(inpcube=int_cube+'.model',
                                            pbcube=int_cube+'.pb',
                                            cubewt=int_cube+'.sumwt',
                                            chanwt=inpparams['chanwt'],
                                            action='div', pblimit=pblimit,freqdep=fdep_pb)

                if usedata!="int":
                    ## copy the int_cube.model to the sd_cube.model
                    shutil.rmtree(sd_cube+'.model',ignore_errors=True)
                    shutil.copytree(int_cube+'.model', sd_cube+'.model')

                if applypb==True:
                    ## Multiply flat-sky model with freq-dep PB
                    sdintlib.modify_with_pb(inpcube=int_cube+'.model',
                                            pbcube=int_cube+'.pb',
                                            cubewt=int_cube+'.sumwt',
                                            chanwt=inpparams['chanwt'],
                                            action='mult', pblimit=pblimit, freqdep=True)

                ## Major cycle for interferometer data
                t0=time.time();
 #               print('Max of int residual before major cycle' + str(imstat(int_cube+'.residual',verbose=False)['max'][0]))
 #               print('Max of int model before major cycle' + str(imstat(int_cube+'.model',verbose=False)['max'][0]))

                if usedata != "sd":
                    imager.runMajorCycle()
 #               print('Max of int residual after major cycle' + str(imstat(int_cube+'.residual',verbose=False)['max'][0]))
                t1=time.time();
                casalog.post("***Time for major cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

                if usedata!="int":
                    ## Major cycle for Single Dish data (uses the flat sky cube model in sd_cube.model )
                    sdintlib.calc_sd_residual(origcube=sd_cube+'.image',
                                              modelcube=sd_cube+'.model',
                                              residualcube=sd_cube+'.residual',  ## output
                                              psfcube=sd_cube+'.psf')

                ## Feather the residuals
                feather_residual(int_cube, sd_cube, joint_cube, applypb, inpparams )
                ###############
                ##### Placeholder code for PSF renormalization if needed
                #sdintlib.apply_renorm(imname=joint_cube+'.residual', sumwtname=joint_cube+'.sumwt')
                ###############

                if specmode=='mfs':
                    ## Calculate Spectral Taylor Residuals
                    sdintlib.cube_to_taylor_sum(cubename=joint_cube+'.residual',
                                                cubewt=int_cube+'.sumwt',
                                                chanwt=inpparams['chanwt'],
                                                mtname=joint_multiterm+'.residual',
                                                nterms=nterms, reffreq=inpparams['reffreq'], dopsf=False)

#                if specmode=='mfs':
#                    print('Max of residual after feather step ' + str(imstat(joint_multiterm+'.residual.tt0',verbose=False)['max'][0]))
#                else:
#                    print('Max of residual after feather step ' + str(imstat(joint_cube+'.residual',verbose=False)['max'][0]))


                deconvolvertool.updateMask()

                ## Get summary from iterbot
                if type(interactive) != bool:
                    #retrec=imager.getSummary();
                    retrec=deconvolvertool.getSummary();

            ## Restore images.
            if restoration==True:  
                t0=time.time();
                deconvolvertool.restoreImages()
                t1=time.time();
                casalog.post("***Time for restoring images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                if pbcor==True:
                #if applypb==True:
                    t0=time.time();
                    if specmode=='mfs':
                        sdintlib.pbcor(imagename=decname+'.image.tt0' ,  pbimage=decname+'.pb.tt0' , cutoff=pblimit,outfile=decname+'.image.tt0.pbcor')
                    else:
                        sdintlib.pbcor(imagename=joint_cube+'.image' ,  pbimage=int_cube+'.pb' , cutoff=pblimit,outfile=joint_cube+'.image.pbcor')

                    #imager.pbcorImages()
                    t1=time.time();
                    casalog.post("***Time for pb-correcting images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                    
        ##close tools
        # needs to deletools before concat or lock waits for ever
        imager.deleteTools()
        deconvolvertool.deleteTools()
   

    finally:
        if imager != None:
            imager.deleteTools() 
        if(cppparallel):
            ###release workers back to python mpi control
            si=synthesisimager()
            si.releasempi()

        #clean up tmp files
        deleteTmpFiles()

    # Write history at the end, when hopefully all temp files are gone from disk,
    # so they won't be picked up. They need time to disappear on NFS or slow hw.
    # Copied from tclean.
    try:
        params = get_func_params(sdintimaging, locals())
        write_tclean_history(imagename, 'sdintimaging', params, casalog)
    except Exception as exc:
        casalog.post("Error updating history (logtable): {} ".format(exc),'WARN')
       
    return retrec

##################################################

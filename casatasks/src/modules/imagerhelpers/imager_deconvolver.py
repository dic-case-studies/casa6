from __future__ import absolute_import
from __future__ import print_function
import os
import math
import shutil
import string
import time
import re
import copy

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import synthesisdeconvolver, iterbotsink, ctsys, table
    from casatasks import casalog

    ctsys_hostinfo = ctsys.hostinfo
    _tb = table() # TODO is this necessary?
else:
    from taskinit import *

    synthesisdeconvolver = casac.synthesisdeconvolver
    # make it look like the CASA6 version even though it's using the CASA5 named tool not present in CASA6
    iterbotsink = casac.synthesisiterbot

    ctsys_hostinfo = casac.cu.hostinfo

    _tb = tb # TODO is this necessary?
'''
A set of helper functions for deconvolve.

Summary...
    
'''

#############################################
class PyDeconvolver:
    def __init__(self,params):
        ################ Tools
        self.initDefaults()

        # Check all input parameters, after partitioning setup.

        # Imaging/Deconvolution parameters. Same for serial and parallel runs
        self.alldecpars = params.getDecPars()
        self.allimpars = params.getImagePars()
        # Iteration parameters
        self.iterpars = params.getIterPars() ## Or just params.iterpars

        # Not necessary for deconvolver:
        ## self.allselpars = params.getSelPars()
        ## self.allgridpars = params.getGridPars()
        ## self.allnormpars = params.getNormPars()
        ## self.weightpars = params.getWeightPars()

        ## Number of fields ( main + outliers )
        self.NF = len(self.allimpars.keys())
        self.stopMinor = {}  ##[0]*self.NF
        for immod in range(0,self.NF):
            self.stopMinor[str(immod)]=1.0
        ## for debug mode automask incrementation only
        self.ncycle = 0

#############################################
    def initializeDeconvolvers(self):
        for immod in range(0,self.NF):
            self.SDtools.append(synthesisdeconvolver())
            self.SDtools[immod].setupdeconvolution(decpars=self.alldecpars[str(immod)])

#############################################
    def setIterationControl(self, iterbot):
        # not that this might produce strange results if the iteration parameters are different for this iterbot than for this imager_deconvolver
        self.IBtool = iterbot

#############################################
    def initializeIterationControl(self):
        # note that in CASA5 this is casac.synthesisiterbot
        self.IBtool = iterbotsink()
        itbot = self.IBtool.setupiteration(iterpars=self.iterpars)

#############################################
    def estimatememory(self):
        griddermem = 0 # no major cycle memory needed

        deconmem=0
        for immod in range(0,self.NF):
            ims= self.allimpars[str(immod)]['imsize']
            if(type(ims)==int) :
                ims=[ims, ims]
            if(len(ims) ==1):
                ims.append(ims[0])
            #print 'shape', self.allimpars[str(immod)]['imsize'], len(ims) 
            #print "DECON mem usage ", self.SDtools[immod].estimatememory(ims)
            if(len(self.SDtools) > immod):
                if(self.SDtools != None):
                    deconmem+=self.SDtools[immod].estimatememory(ims)

        availmem=ctsys_hostinfo()['memory']['available']
        if((deconmem+griddermem) > 0.8*availmem):
            casalog.post("Memory available "+str(availmem)+" kB is very close to amount of required memory "+str(deconmem+griddermem)+" kB" , "WARN")
        else:
            casalog.post("Memory available "+str(availmem)+" kB and  required memory "+str(deconmem+griddermem)+" kB" , "INFO2")

############################################
    def restoreImages(self):
        for immod in range(0,self.NF):
              self.SDtools[immod].restore()

#############################################
    def pbcorImages(self):
         for immod in range(0,self.NF):
              self.SDtools[immod].pbcor()

#############################################
    def getSummary(self,fignum=1):
        summ = self.IBtool.getiterationsummary()
        # self.plotReport( summ, fignum ) should be called explicitly for deconvolver
        return summ

#############################################
    def deleteDeconvolvers(self):
         for immod in range(0,len(self.SDtools)):
              self.SDtools[immod].done()

    def deleteIterBot(self):
         if self.IBtool != None:
              self.IBtool.done()

    def initDefaults(self):
        # Reset globals/members
         self.NF=1
         self.stopMinor={'0':1.0}  # Flag to call minor cycle for this field or not.
         self.NN=1
         self.SItool=None
         self.SDtools=[]
         self.PStools=[]
         self.IBtool=None
    
#############################################

    def deleteTools(self):
        self.deleteDeconvolvers()
        self.deleteIterBot()
        self.initDefaults()

#############################################

    def hasConverged(self):
        # Merge peak-res info from all fields to decide iteration parameters
         self.IBtool.resetminorcycleinfo() 
         for immod in range(0,self.NF):
              initrec =  self.SDtools[immod].initminorcycle() 
              self.IBtool.mergeinitrecord( initrec );

         # Check with the iteration controller about convergence.
         #print("check convergence")
         stopflag = self.IBtool.cleanComplete()
         #print('Converged : ', stopflag)
         if( stopflag>0 ):
             stopreasons = ['iteration limit', 'threshold', 'force stop','no change in peak residual across two major cycles', 'peak residual increased by more than 3 times from the previous major cycle','peak residual increased by more than 3 times from the minimum reached','zero mask', 'any combination of n-sigma and other valid exit criterion']
             casalog.post("Reached global stopping criterion : " + stopreasons[stopflag-1], "INFO")

             # revert the current automask to the previous one 
             #if self.iterpars['interactive']:
             for immod in range(0,self.NF):
                     if self.alldecpars[str(immod)]['usemask'].count('auto')>0:
                        prevmask = self.allimpars[str(immod)]['imagename']+'.prev.mask'
                        if os.path.isdir(prevmask):
                          # Try to force rmtree even with an error as an nfs mounted disk gives an error 
                          #shutil.rmtree(self.allimpars[str(immod)]['imagename']+'.mask')
                          shutil.rmtree(self.allimpars[str(immod)]['imagename']+'.mask', ignore_errors=True)
                          # For NFS mounted disk it still leave .nfs* file(s) 
                          if os.path.isdir(self.allimpars[str(immod)]['imagename']+'.mask'):
                              import glob
                              if glob.glob(self.allimpars[str(immod)]['imagename']+'.mask/.nfs*'):
                                  for item in os.listdir(prevmask):
                                      src = os.path.join(prevmask,item)
                                      dst = os.path.join(self.allimpars[str(immod)]['imagename']+'.mask',item)
                                      if os.path.isdir(src):
                                          shutil.move(src, dst)
                                      else:
                                          shutil.copy2(src,dst)
                              shutil.rmtree(prevmask)
                          else: 
                              shutil.move(prevmask,self.allimpars[str(immod)]['imagename']+'.mask')
                          casalog.post("[" + str(self.allimpars[str(immod)]['imagename']) + "] : Reverting output mask to one that was last used ", "INFO")

         return (stopflag>0)

#############################################
    def updateMask(self, synthesisImager):
        # Run update mask for this instance
        maskchanged = self.updateMaskMinor()
        
        # Run interactive masking with the major cycle imager
        maskchanged = maskchanged | synthesisImager.runInteractiveGUI2()

        ## Return a flag to say that the mask has changed or not.
        return maskchanged

#############################################
    def updateMaskMinor(self):
        # Setup mask for each field ( input mask, and automask )
        maskchanged = False
        for immod in range(0,self.NF):
            maskchanged = maskchanged | self.SDtools[immod].setupmask()

        ## Return a flag to say that the mask has changed or not.
        return maskchanged

#############################################
    def runMinorCycleCore(self):

        # Set False for release packages. 
        # Only set this to True for testing and debugging automask in parallel mode
        # since in parallel mode, runtime setting of the enviroment variable
        # currently does not work.
        # False = disable always save intermediate images mode
        alwaysSaveIntermediateImages=False

        # Get iteration control parameters
        iterbotrec = self.IBtool.getminorcyclecontrols()
        ##print("Minor Cycle controls : ", iterbotrec)

        self.IBtool.resetminorcycleinfo() 

        #
        # Run minor cycle
        self.ncycle+=1
        for immod in range(0,self.NF):  
            if self.stopMinor[str(immod)]<3 :

                # temporarily disable the check (=> always save the intermediate images
                if alwaysSaveIntermediateImages or ('SAVE_ALL_RESIMS' in os.environ and os.environ['SAVE_ALL_RESIMS']=="true"):
                    resname = self.allimpars[str(immod)]['imagename']+'.residual'
                    tempresname = self.allimpars[str(immod)]['imagename']+'.inputres'+str(self.ncycle)
                    if os.path.isdir(resname):
                        shutil.copytree(resname, tempresname)

                exrec = self.SDtools[immod].executeminorcycle( iterbotrecord = iterbotrec )

                #print('.... iterdone for ', immod, ' : ' , exrec['iterdone'])
                self.IBtool.mergeexecrecord( exrec )
                if alwaysSaveIntermediateImages or ('SAVE_ALL_AUTOMASKS' in os.environ and os.environ['SAVE_ALL_AUTOMASKS']=="true"):
                    maskname = self.allimpars[str(immod)]['imagename']+'.mask'
                    tempmaskname = self.allimpars[str(immod)]['imagename']+'.autothresh'+str(self.ncycle)
                    if os.path.isdir(maskname):
                        shutil.copytree(maskname, tempmaskname)
                
                # Some what duplicated as above but keep a copy of the previous mask
                # for interactive automask to revert to it if the current mask
                # is not used (i.e. reached deconvolution stopping condition).
                #if self.iterpars['interactive'] and self.alldecpars[str(immod)]['usemask']=='auto-thresh':
                if self.alldecpars[str(immod)]['usemask'].count('auto')>0:
                    maskname = self.allimpars[str(immod)]['imagename']+'.mask'
                    prevmaskname=self.allimpars[str(immod)]['imagename']+'.prev.mask'
                    if os.path.isdir(maskname):
                        if os.path.isdir(prevmaskname):
                            shutil.rmtree(prevmaskname)
                        shutil.copytree(maskname, prevmaskname)

#######################################################
#######################################################


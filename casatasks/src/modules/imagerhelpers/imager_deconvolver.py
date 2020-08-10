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
        self.activities=activities

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

        # CFCache params
        self.cfcachepars = params.getCFCachePars()
        ## Number of fields ( main + outliers )
        self.NF = len(self.allimpars.keys())
        self.stopMinor = {}  ##[0]*self.NF
        for immod in range(0,self.NF):
            self.stopMinor[str(immod)]=1.0
        ## Number of nodes. This gets set for parallel runs
        ## It can also be used serially to process the major cycle in pieces.
        self.NN = 1 
        ## for debug mode automask incrementation only
        self.ncycle = 0
#        isvalid = self.checkParameters()
#        if isvalid==False:
#            print('Invalid parameters')

#############################################
    def checkParameters(self):
        return True

#############################################
    def makeCFCache(self,exists):
        pass # convolution functions not needed for deconvolver
        
#############################################
    def initializeImagers(self):
        pass # imager functions not needed for deconvolver

#############################################

    def initializeDeconvolvers(self):
        for immod in range(0,self.NF):
            self.SDtools.append(synthesisdeconvolver())
            self.SDtools[immod].setupdeconvolution(decpars=self.alldecpars[str(immod)])
             

#############################################
    ## Overloaded by ParallelCont
    def initializeNormalizers(self):
        pass # normalizer functions not needed for deconvolver

#############################################

    def initializeIterationControl(self):
        # note that in CASA5 this is casac.synthesisiterbot
        self.IBtool = iterbotsink()
        itbot = self.IBtool.setupiteration(iterpars=self.iterpars)

#############################################
    def estimatememory(self):
        #print "MEMORY usage ", self.SItool.estimatememory(), type(self.SItool.estimatememory())

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
    def deleteImagers(self):
        pass # no SItool for deconvolver

    def deleteDeconvolvers(self):
         for immod in range(0,len(self.SDtools)):
              self.SDtools[immod].done()

    def deleteNormalizers(self):
         pass # no PStools for deconvolver

    def deleteIterBot(self):
         if self.IBtool != None:
              self.IBtool.done()

    def deleteCluster(self):
#         print('no cluster to delete')
        return

    def deleteWorkDir(self):
        # No .workdirectory to delete
        return

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
        self.deleteImagers()
        self.deleteDeconvolvers()
        self.deleteNormalizers()
        self.deleteIterBot()
        self.deleteWorkDir()
        self.initDefaults()
        self.deleteCluster()

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
    def runInteractiveGUI2(self):
        return False # to be run in between major/minor cycle loops

#############################################
    def makePSF(self):
        pass # PSF shoudl have been generated prior to deconvolving

#############################################
    def calcVisAppSens(self):
        pass # no SItool for deconvolver


#############################################
    def runMajorCycle(self):
        raise Exception("Trying to run major cycle in PyDeconvolver. Use PySynthesisImager instead.")

#############################################
    def predictModel(self):
        pass # model should have been created prior to deconvolving

#############################################
    def dryGridding(self):
        pass # no gridding for deconvolver

#############################################
## Overloaded for parallel runs
    def fillCFCache(self):
        pass # convolution functions not needed for deconvolver
                  
#############################################
    def reloadCFCache(self):
        pass # convolution functions not needed for deconvolver

#############################################
    def makePB(self):
        pass # pb should have been generated prior to deconvolving

#############################################
    def makePBCore(self):
        pass # pb core should have been generated prior to deconvolving

#############################################
    def makeSdImage(self):
        pass # Sd image should have been generated prior to deconvolving

#############################################
    def makeSdPSF(self):
        pass # PSF should have been generated prior to deconvolving

#############################################
    def makeSdImageCore(self):
        pass # Sd image should have been generated prior to deconvolving

#############################################
    def makeSdPSFCore(self):
        pass # PSF should have been generated prior to deconvolving

#############################################
    def makeImage(self, imagetype='observed', image='', compleximage='', imagefieldid=0):
        pass # dirty image should have been generated prior to deconvolving

#############################################
## Overloaded for parallel runs
    def setWeighting(self):
        pass # not needed for deconvolver
        
#############################################
## Overloaded for parallel runs
    def makePSFCore(self):
        pass # PSF should have been generated prior to deconvolving

#############################################
## Overloaded for parallel runs
    def runMajorCycleCore(self, lastcycle):
        raise Exception("Trying to run major cycle in PyDeconvolver. Use PySynthesisImager instead.")

#############################################
## Overloaded for parallel runs
    def predictModelCore(self):
        pass # model should have been created prior to deconvolving

#############################################
    def runMinorCycle(self):
        self.runMinorCycleCore()

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

#############################################
    def runMajorMinorLoops(self):
         passraise Exception("Trying to run major cycle loop in PyDeconvolver. Use PySynthesisImager instead.")

#############################################
    def plotReport( self, summ={} ,fignum=1 ):

        if not ( 'summaryminor' in summ and 'summarymajor' in summ and 'threshold' in summ and summ['summaryminor'].shape[0]==6 ):
            print('Cannot make summary plot. Please check contents of the output dictionary from tclean.')
            return summ

        import pylab as pl
        from numpy import max as amax

        # 0 : iteration number (within deconvolver, per cycle)
        # 1 : peak residual
        # 2 : model flux
        # 3 : cyclethreshold
        # 4 : deconvolver id
        # 5 : subimage id (channel, stokes..)

        pl.ioff()

        fig, ax = pl.subplots(nrows=1,ncols=1,num=fignum)
        pl.clf();
        minarr = summ['summaryminor']
        if minarr.size==0:
            casalog.post("Zero iteration: no summary plot is generated.", "WARN")
        else:
            iterlist = minarr[0,:]
            eps=0.0
            peakresstart=[]
            peakresend=[]
            modfluxstart=[]
            modfluxend=[]
            itercountstart=[]
            itercountend=[]
            peakresstart.append( minarr[1,:][0] )
            modfluxstart.append( minarr[2,:][0] )
            itercountstart.append( minarr[0,:][0] + eps )
            peakresend.append( minarr[1,:][0] )
            modfluxend.append( minarr[2,:][0] )
            itercountend.append( minarr[0,:][0] + eps )
            for ii in range(0,len(iterlist)-1):
                if iterlist[ii]==iterlist[ii+1]:
                    peakresend.append( minarr[1,:][ii] )
                    peakresstart.append( minarr[1,:][ii+1] ) 
                    modfluxend.append( minarr[2,:][ii] )
                    modfluxstart.append( minarr[2,:][ii+1] )
                    itercountend.append( iterlist[ii]-eps )
                    itercountstart.append( iterlist[ii]+eps )

            peakresend.append( minarr[1,:][len(iterlist)-1] )
            modfluxend.append( minarr[2,:][len(iterlist)-1] )
            itercountend.append( minarr[0,:][len(iterlist)-1] + eps )

    #        pl.plot( iterlist , minarr[1,:] , 'r.-' , label='peak residual' , linewidth=1.5, markersize=8.0)
    #        pl.plot( iterlist , minarr[2,:] , 'b.-' , label='model flux' )
    #        pl.plot( iterlist , minarr[3,:] , 'g--' , label='cycle threshold' )

            pl.plot( itercountstart , peakresstart , 'r.--' , label='peak residual (start)')
            pl.plot( itercountend , peakresend , 'r.-' , label='peak residual (end)',linewidth=2.5)
            pl.plot( itercountstart , modfluxstart , 'b.--' , label='model flux (start)' )
            pl.plot( itercountend , modfluxend , 'b.-' , label='model flux (end)',linewidth=2.5 )
            pl.plot( iterlist , minarr[3,:] , 'g--' , label='cycle threshold', linewidth=2.5 )
            maxval = amax( minarr[1,:] )
            maxval = max( maxval, amax( minarr[2,:] ) )
            
            bcols = ['b','g','r','y','c']
            minv=1
            niterdone = len(minarr[4,:])
          
            if len(summ['summarymajor'].shape)==1 and summ['summarymajor'].shape[0]>0 :       
                pl.vlines(summ['summarymajor'],0,maxval, label='major cycles', linewidth=2.0)

            pl.hlines( summ['threshold'], 0, summ['iterdone'] , linestyle='dashed' ,label='threshold')
        
            pl.xlabel( 'Iteration Count' )
            pl.ylabel( 'Peak Residual (red), Model Flux (blue)' )

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width, box.height*0.8])

            pl.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),
                      ncol=3, fancybox=True, shadow=True)

            pl.savefig('summaryplot_'+str(fignum)+'.png')
            pl.ion()

        return summ;

    #############################################
    def unlockimages( self, imageid=0 ):
        return False # nothing for deconvolver to unlock

#######################################################
#######################################################


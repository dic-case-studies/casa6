from __future__ import absolute_import
import os
import math
import shutil
import string
import time
import re;
import copy
import pdb
from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import synthesisimager, synthesisnormalizer
    from casatasks import casalog

    from .imager_base import PySynthesisImager
    from .parallel_imager_helper import PyParallelImagerHelper
    synth_imager_name = 'synthesisimager'
    synth_imager_import = 'from casatools import synthesisimager'

else:
    from taskinit import *

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    synthesisimager = casac.synthesisimager
    synthesisnormalizer = casac.synthesisnormalizer
    synth_imager_name = 'casac.synthesisimager'
    synth_imager_import = 'pass'


'''
An implementation of parallel continuum imaging, using synthesisxxxx tools

Datasets are partitioned by row and major cycles are parallelized. 
Gathers and normalization are done before passing the images to a
non-parallel minor cycle. The output model image is them scattered to
all the nodes for the next parallel major cycle.

There are N synthesisimager objects.
There is 1 instance per image field, of the normalizer and deconvolver.
There is 1 iterbot. 
    
'''

#############################################
#############################################
## Parallelize only major cycle.
#############################################
class PyParallelContSynthesisImager(PySynthesisImager):

    def __init__(self,params=None):

         PySynthesisImager.__init__(self,params)

         self.PH = PyParallelImagerHelper()
         self.NN = self.PH.NN
         self.selpars = self.allselpars;
         self.allselpars = self.PH.partitionContDataSelection(self.allselpars)
         # self.allcflist = self.PH.partitionCFCacheList(self.cfcachepars['cflist']);
         # self.allcflist = self.PH.partitionCFCacheList(self.allgridpars['0']);
         self.listOfNodes = self.PH.getNodeList();
         self.coordsyspars = {};
         self.toolsi=None

#############################################
    def resetSaveModelParams(self, params=None):
         mainparams = params.getSelPars()
         for n in self.allselpars: # for all nodes
             for v in self.allselpars[n]: # for all MSes
                 self.allselpars[n][v]['readonly']=mainparams[v]['readonly']
                 self.allselpars[n][v]['usescratch']=mainparams[v]['usescratch']

#############################################
    def initializeImagers(self):
        ### Drygridding, and Coordsys comes from a single imager on MAIN node.
        ### No startmodel confusion. It's created only once and then scattered.
        self.initializeImagers()

        ### Note : Leftover from CAS-9977 
        ### There is a coord system mismatch at scatter/gather, if the MAIN version already
        ###   exists on disk. With startmodel, it's xxx.model.  With aproject, it's xxx.residual.
        ### There is an exception in SIImageStore::openImage to handle this. 
        ### Turn on casalog.filter('DEBUG1') to see the warning message.


#############################################
    def initializeImagersBase(self,thisSelPars,partialSelPars):

        if partialSelPars==False: ## Init only on the zero'th node

            #
            # Use the already-created imager on MAIN node
            #
            ##self.toolsi = synthesisimager()

            #
            # Select data. 
            #
            for mss in sorted( self.selpars.keys() ):
                self.toolsi.selectdata( thisSelPars[mss] )

            # Defineimage. 
            # This makes the global csys. Get csys to distribute to other nodes
            # It also sets 'startmodel' if available (this is later scattered to nodes)
            for fld in range(0,self.NF):
                tmpimpars = copy.deepcopy(self.allimpars[str(fld)])
                #if tmpimpars.has_key('startmodel'):
                #    tmpimpars.pop('startmodel')
                self.toolsi.defineimage( impars=tmpimpars, gridpars = self.allgridpars[str(fld)] )
                fullcoords = self.toolsi.getcsys()
                self.coordsyspars[str(fld)] = fullcoords

            # Modify the coordsys inputs
            for fld in range(0, self.NF):
                self.allimpars[str(fld)]['csys']=self.coordsyspars[str(fld)]['coordsys'].copy()

            # Call the global defineimage again 
            #  (to get around later error of different coordsys latpoles! (CAs-9977)
            #for fld in range(0,self.NF):
            #    self.toolsi.defineimage( impars=self.allimpars[str(fld)], gridpars = self.allgridpars[str(fld)] )
            


        else: ## partialSelPars==True , The actual initialization on all nodes.

            #
            # Start the imagers on all nodes.
            #
            joblist=[]
            for node in self.listOfNodes:
                joblist.append( self.PH.runcmd('{0}; toolsi = {1}()'.format(
                    synth_imager_import, synth_imager_name),
                                               node) );
            self.PH.checkJobs(joblist);

            #
            # Select data.  If partialSelPars is True, use the thisSelPars
            # data structure as a list of partitioned selections.
            #
            joblist=[];
            nodes=self.listOfNodes;#[1];
            for node in nodes:
                for mss in sorted( self.selpars.keys() ):
                    selStr=str(thisSelPars[str(node-1)][mss]);
                    joblist.append( self.PH.runcmd("toolsi.selectdata( "+selStr+")", node) )
            self.PH.checkJobs(joblist);

            #
            # Call defineimage at each node.
            #
            joblist=[];
            for node in nodes:
                ## For each image-field, define imaging parameters
                nimpars = copy.deepcopy(self.allimpars)
                # casalog.post("nimpars = "+str(nimpars))
                ngridpars = copy.deepcopy(self.allgridpars)
                for fld in range(0,self.NF):
                    if self.NN>1:
                        nimpars[str(fld)]['imagename'] = self.PH.getpartimagename( nimpars[str(fld)]['imagename'], node )

                    ## Pop out the startmodel, as it would already have been created on the main node,.
                    tmpimpars = nimpars[str(fld)]
                    if 'startmodel' in tmpimpars:
                        tmpimpars.pop('startmodel')

                    joblist.append( self.PH.runcmd("toolsi.defineimage( impars=" + str( nimpars[str(fld)] ) 
                                                   + ", gridpars=" + str( ngridpars[str(fld)] )   + ")", node ) )
            self.PH.checkJobs(joblist);
        
#############################################

    def initializeImagers(self):

        #---------------------------------------
        #  Check if cfcache exists.
        #
        cfCacheName=''
        if(self.allgridpars['0']['gridder'].startswith('awp')):
            cfCacheName=self.allgridpars['0']['cfcache']
        else:
            self.allgridpars['0']['cfcache']=''
        cfcExists=False
        if(self.allgridpars['0']['gridder'] == 'awproject' or self.allgridpars['0']['gridder'] == 'awprojectft'):
            if (cfCacheName == ''):
                cfCacheName = self.allimpars['0']['imagename'] + '.cf'
                cfCacheName=self.allgridpars['0']['cfcache'] = cfCacheName
 
            cfcExists = (os.path.exists(cfCacheName) and os.path.isdir(cfCacheName));
            if (cfcExists):
                nCFs = len(os.listdir(cfCacheName));
                if (nCFs == 0):
                    casalog.post(cfCacheName + " exists, but is empty.  Attempt is being made to fill it now.","WARN")
                    cfcExists = False;
        # casalog.post("##########################################")
        # casalog.post("CFCACHE = "+cfCacheName,cfcExists)
        # casalog.post("##########################################")

       
        # Start one imager on MAIN node
        self.toolsi = synthesisimager()

        # Init one SI tool ( it records the csys per field in self.coordsyspars )
        self.initializeImagersBase(self.selpars,False);

        # Modify the coordsys inputs
#        for fld in range(0, self.NF):
#            self.allimpars[str(fld)]['csys']=self.coordsyspars[str(fld)]['coordsys'].copy()

        # Dry Gridding on the MAIN node ( i.e. on self.toolsi)
        if (not cfcExists):
            self.dryGridding();

        ##weighting with mosfield=True
        if( (self.weightpars['type']=='briggs')  and (self.weightpars['multifield'])):
            self.toolsi.setweighting(**self.weightpars)
            ###master create the weight density for all fields
            self.toolsi.getweightdensity()
            
        # Clean up the single imager (MAIN node)
        self.toolsi.done()
        self.toolsi = None

        # Do the second round, initializing imagers on ALL nodes
        self.initializeImagersBase(self.allselpars,True);

        # Fill CFCache - it uses all nodes.
        if (not cfcExists):
            self.fillCFCache();
        self.reloadCFCache();

######################################################################################################################################
        #---------------------------------------
        #  4. call setdata() for images on all nodes
        #
        # joblist=[];
        # for node in self.listOfNodes:
        #     ## Send in Selection parameters for all MSs in the list
        #     #### MPIInterface related changes (the -1 in the expression below)
        #     for mss in sorted( (self.allselpars[str(node-1)]).keys() ):
        #         joblist.append( self.PH.runcmd("toolsi.selectdata( "+str(self.allselpars[str(node-1)][mss])+")", node) )
        # self.PH.checkJobs(joblist);

        #---------------------------------------
        #  5. Call defineImage() on all nodes.  This sets up the FTMs.
        #
#         joblist=[];
#         for node in self.listOfNodes:
#             ## For each image-field, define imaging parameters
#             nimpars = copy.deepcopy(self.allimpars)
#             #casalog.post("nimpars = "+str(nimpars))
#             ngridpars = copy.deepcopy(self.allgridpars)
#             for fld in range(0,self.NF):
#                 if self.NN>1:
#                     nimpars[str(fld)]['imagename'] = self.PH.getpath(node) + '/' + nimpars[str(fld)]['imagename']+'.n'+str(node)
# ###                    nimpars[str(fld)]['imagename'] = self.allnormpars[str(fld)]['workdir'] + '/' + nimpars[str(fld)]['imagename']+'.n'+str(node)
# ###                    nimpars[str(fld)]['imagename'] = nimpars[str(fld)]['imagename']+'.n'+str(node)

# #                    ngridpars[str(fld)]['cfcache'] = ngridpars[str(fld)]['cfcache']+'.n'+str(node)
#                     # # Give the same CFCache name to all nodes
#                     ngridpars[str(fld)]['cfcache'] = ngridpars[str(fld)]['cfcache'];

#                 joblist.append( self.PH.runcmd("toolsi.defineimage( impars=" + str( nimpars[str(fld)] ) + ", gridpars=" + str( ngridpars[str(fld)] )   + ")", node ) )
#         self.PH.checkJobs(joblist);

        #---------------------------------------
        #  6. If cfcache does not exist, call fillCFCache()
        #       This will fill the "empty" CFCache in parallel
        #  7. Now call reloadCFCache() on all nodes.
        #     This reloads the latest cfcahce.



        # TRY: Start all over again!
        # self.deleteImagers();

        # joblist=[]

        # for node in self.listOfNodes:
        #     joblist.append( self.PH.runcmd("toolsi = casac.synthesisimager()", node) );
        # self.PH.checkJobs(joblist);

        # joblist=[];
        # nodes=self.listOfNodes;#[1];
        # for node in nodes:
        #     for mss in sorted( (self.allselpars[str(node-1)]).keys() ):
        #         joblist.append( self.PH.runcmd("toolsi.selectdata( "+str(self.allselpars[str(node-1)][mss])+")", node) )
        #             # for mss in sorted( self.selpars.keys() ):
        #             #     joblist.append( self.PH.runcmd("toolsi.selectdata( "+str(self.selpars[mss])+")", node) )
        # self.PH.checkJobs(joblist);

        # joblist=[];
        # for node in self.listOfNodes:
        #     nimpars = copy.deepcopy(self.allimpars)
        #     ngridpars = copy.deepcopy(self.allgridpars)
        #     for fld in range(0,self.NF):
        #         if self.NN>1:
        #             nimpars[str(fld)]['imagename'] = self.PH.getpath(node) + '/' + nimpars[str(fld)]['imagename']+'.n'+str(node)
        #             # # Give the same CFCache name to all nodes
        #             ngridpars[str(fld)]['cfcache'] = ngridpars[str(fld)]['cfcache'];

        #         joblist.append( self.PH.runcmd("toolsi.defineimage( impars=" + str( nimpars[str(fld)] ) + ", gridpars=" + str( ngridpars[str(fld)] )   + ")", node ) )
        # self.PH.checkJobs(joblist);


#############################################


#############################################

    def initializeNormalizers(self):
        for immod in range(0,self.NF):
            self.PStools.append(synthesisnormalizer())
            normpars = copy.deepcopy( self.allnormpars[str(immod)] )
            partnames = []
            if(self.NN>1):
                #### MPIInterface related changes
                #for node in range(0,self.NN):
                for node in self.listOfNodes:
                    partnames.append( self.PH.getpartimagename( self.allimpars[str(immod)]['imagename'], node ) )
                    #onename = self.allimpars[str(immod)]['imagename']+'.n'+str(node)
                    #partnames.append( self.PH.getpath(node) + '/' + onename  )
                    #self.PH.deletepartimages( self.PH.getpath(node), onename ) # To ensure restarts work properly.
                    self.PH.deletepartimages( self.allimpars[str(immod)]['imagename'] ,  node ) # To ensure restarts work properly.
                normpars['partimagenames'] = partnames
            self.PStools[immod].setupnormalizer(normpars=normpars)


#############################################
    def setWeighting(self):

        ## Set weight parameters and accumulate weight density (natural)
        joblist=[];
        if( (self.weightpars['type']=='briggs')  and (self.weightpars['multifield'])):
            ###master created the weight density for all fields
            ##Should have been in  initializeImagersBase_New but it is not being called !
            self.toolsi = synthesisimager()
            for mss in sorted( self.selpars.keys() ):
                self.toolsi.selectdata( self.selpars[mss] )
            for fld in range(0,self.NF):
                self.toolsi.defineimage( impars=self.allimpars[str(fld)], gridpars = self.allgridpars[str(fld)] )
            self.toolsi.setweighting(**self.weightpars)
            ###master create the weight density for all fields
            weightimage=self.toolsi.getweightdensity()
            self.toolsi.done()
            self.toolsi=None
            destWgtim=weightimage+'_moswt'
            shutil.move(weightimage, destWgtim)
            joblist=[];
            for node in self.listOfNodes:
                joblist.append( self.PH.runcmd("toolsi.setweightdensity('"+str(destWgtim)+"')", node ) )
            self.PH.checkJobs( joblist )
            #for node in self.listOfNodes:
            #    ## Set weighting pars
            #   joblist.append( self.PH.runcmd("toolsi.setweighting( **" + str(self.weightpars) + ")", node ) )
            #self.PH.checkJobs( joblist )
            #joblist=[];
            #for node in self.listOfNodes:
            #    joblist.append( self.PH.runcmd("toolsi.getweightdensity()", node ) )
            #self.PH.checkJobs( joblist )
            
            #for immod in range(0,self.NF):
            #    #self.PStools[immod].gatherweightdensity()
             #   self.PStools[immod].scatterweightdensity()
            ## Set weight density for each nodel
            #joblist=[];
            #for node in self.listOfNodes:
             #   joblist.append( self.PH.runcmd("toolsi.setweightdensity()", node ) )
            #self.PH.checkJobs( joblist )
       #### end of multifield or mosweight
        else:
            joblist=[];
            for node in self.listOfNodes:
                ## Set weighting pars
                joblist.append( self.PH.runcmd("toolsi.setweighting( **" + str(self.weightpars) + ")", node ) )
            self.PH.checkJobs( joblist )

            ## If only one field, do the get/gather/set of the weight density.
            if self.NF == 1 and self.allimpars['0']['stokes']=="I":   ## Remove after gridded wts appear for all fields correctly (i.e. new FTM).
                
                if not ( (self.weightpars['type'] ==  'natural') or (self.weightpars['type'] == 'radial'))   :  ## For natural and radial, this array isn't created at all.
                                                                       ## Remove when we switch to new FTM

                    casalog.post("Gathering/Merging/Scattering Weight Density for PSF generation","INFO")

                    joblist=[];
                    for node in self.listOfNodes:
                        joblist.append( self.PH.runcmd("toolsi.getweightdensity()", node ) )
                    self.PH.checkJobs( joblist )

                    ## gather weightdensity and sum and scatter
                    casalog.post("******************************************************")
                    casalog.post(" gather and scatter now ")
                    casalog.post("******************************************************")
                    for immod in range(0,self.NF):
                        self.PStools[immod].gatherweightdensity()
                        self.PStools[immod].scatterweightdensity()

                    ## Set weight density for each nodel
                    joblist=[];
                    for node in self.listOfNodes:
                        joblist.append( self.PH.runcmd("toolsi.setweightdensity()", node ) )
                    self.PH.checkJobs( joblist )



    def deleteImagers(self):
         self.PH.runcmd("toolsi.done()")

    def deleteWorkDir(self):
        ## Delete the contents of the .workdirectory 
        for immod in range(0,self.NF):
            normpars = copy.deepcopy( self.allnormpars[str(immod)] )
            if(self.NN>1):
                for node in self.listOfNodes:
                    self.PH.deletepartimages( self.allimpars[str(immod)]['imagename'],  node ,deldir=True ) 

#        ## Delete the workdirectory
#        casalog.post("Deleting workdirectory : "+self.PH.getworkdir(imagename, node))
#        shutil.rmtree( self.PH.getworkdir(imagename, node) )

    def deleteCluster(self):
        self.PH.takedownCluster()
    
# #############################################
    def dryGridding(self):
        dummy=['']
        self.toolsi.drygridding(dummy)

#    def dryGridding_Old(self):
#        nodes=[1];
#        joblist=[];
#        for node in nodes:
#            dummy=[''];
#            cmd = "toolsi.drygridding("+str(dummy)+")";
#            joblist.append(self.PH.runcmd(cmd,node));
#        self.PH.checkJobs(joblist);

#############################################
    def reloadCFCache(self):
        joblist=[];
        for node in self.listOfNodes:
            cmd = "toolsi.reloadcfcache()";
            casalog.post("reloadCFCache, CMD = {} {}".format(node, cmd))
            joblist.append(self.PH.runcmd(cmd,node));
        self.PH.checkJobs(joblist);
#############################################
    def fillCFCache(self):
        #casalog.post("-----------------------fillCFCache------------------------------------")
        # cflist=[f for f in os.listdir(self.allgridpars['cfcache']) if re.match(r'CFS*', f)];
        # partCFList = 
        if(not str(self.allgridpars['0']['gridder']).startswith("awp")):
            return
         
        allcflist = self.PH.partitionCFCacheList(self.allgridpars['0']);
        cfcPath = "\""+str(self.allgridpars['0']['cfcache'])+"\"";
        ftmname = "\""+str(self.allgridpars['0']['gridder'])+"\"";
        psTermOn = str(self.allgridpars['0']['psterm']);
        aTermOn = str(self.allgridpars['0']['aterm']);
        conjBeams = str(self.allgridpars['0']['conjbeams']);
        #aTermOn = str(True);
        # casalog.post("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        # casalog.post("AllCFList = ",allcflist)
        m = len(allcflist);
        # casalog.post("No. of nodes used: " + m,cfcPath,ftmname)
        # casalog.post("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

        joblist=[];
        for node in self.listOfNodes[:m]:
            # casalog.post("#!$#!%#!$#@$#@$ " + allcflist)
            cmd = "toolsi.fillcfcache("+str(allcflist[node])+","+str(ftmname)+","+str(cfcPath)+","+psTermOn+","+aTermOn+","+conjBeams+")";
            # casalog.post("CMD = " + str(node) +" " + cmd)
            joblist.append(self.PH.runcmd(cmd,node));
        self.PH.checkJobs(joblist);

        # Linear code
        # cfcName = self.allgridpars['0']['cfcache'];
        # cflist=[f for f in os.listdir(cfcName) if re.match(r'CFS*', f)];
        # self.cfcachepars['cflist']=cflist;
        # self.SItool.fillcfcache(**(self.cfcachepars)) ;
#############################################
    def makePSFCore(self):
        ### Make PSFs
        joblist=[]
        #### MPIInterface related changes
        #for node in range(0,self.PH.NN):
        for node in self.listOfNodes:
             joblist.append( self.PH.runcmd("toolsi.makepsf()",node) )
        self.PH.checkJobs( joblist ) # this call blocks until all are done.

#############################################
    def makePBCore(self):
        joblist=[]
        # Only one node needs to make the PB. It reads the freq from the image coordsys
        joblist.append( self.PH.runcmd("toolsi.makepb()",self.listOfNodes[0]) )
        self.PH.checkJobs( joblist )

#############################################

    def runMajorCycleCore(self, lastcycle):
        casalog.post("-----------------------------  Running Parallel Major Cycle ----------------------------","INFO")
        ### Run major cycle
        joblist=[]
        #### MPIInterface related changes
        #for node in range(0,self.PH.NN):
        for node in self.listOfNodes:
             joblist.append( self.PH.runcmd("toolsi.executemajorcycle(controls={'lastcycle':"+str(lastcycle)+"})",node) )
        self.PH.checkJobs( joblist ) # this call blocks until all are done.

#############################################
    def predictModelCore(self):
        joblist=[]
        #### MPIInterface related changes
        #for node in range(0,self.PH.NN):
        for node in self.listOfNodes:
             joblist.append( self.PH.runcmd("toolsi.predictmodel()",node) )
        self.PH.checkJobs( joblist ) # this call blocks until all are done.

    def estimatememory(self):
        joblist=[]
        for node in self.listOfNodes:
            joblist.append( self.PH.runcmd("toolsi.estimatememory()", node) )
        self.PH.checkJobs( joblist )

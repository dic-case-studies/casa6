from __future__ import absolute_import
import os
import math
import shutil
import string
import time
import re
import copy

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import synthesisutils, synthesisimager
    from casatools import image as imageanalysis
    from casatasks import casalog

    from .imager_base import PySynthesisImager
    from .parallel_imager_helper import PyParallelImagerHelper
else:
    from taskinit import *
    from casac import casac

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    synthesisimager = casac.synthesisimager
    synthesisutils = casac.synthesisutils
    imageanalysis = casac.image

'''
An implementation of parallel cube imaging, using synthesisxxxx tools.

Major and minor cycles are parallelized across frequency, by running separate
PySynthesisImagers independently per frequency chunk.
Iteration control is not synchronized,  interactive mask drawing can't be done.
Reference concatenation of all the image products is done at the end.

There are N PySynthesisImager objects, each with their own 
synthesisimager, deconvolvers, normalizers and iterbot. 
   
'''

#############################################
# Parallelize both the major and minor cycle for Cube imaging
# Run a separate instance of PySynthesisImager on each node.
#### ( later, use the live-object interface of ImStore to reference-break the cubes )
#### For nprocesses > nnodes, run the whole 'clean' loop multiple times. 
#############################################
class PyParallelCubeSynthesisImager():

    def __init__(self,params=None):

        self.params=params

        allselpars = params.getSelPars()
        allimagepars = params.getImagePars()
        self.allinimagepars = copy.deepcopy(allimagepars)
        self.allgridpars = params.getGridPars()
        self.allnormpars = params.getNormPars()
        self.weightpars = params.getWeightPars()
        self.decpars = params.getDecPars()
        self.iterpars = params.getIterPars()
        alldataimpars={}
         
        self.PH = PyParallelImagerHelper()
        self.NN = self.PH.NN
        self.NF = len(allimagepars.keys())
        self.listOfNodes = self.PH.getNodeList();
        ## Partition both data and image coords the same way.
        #self.allselpars = self.PH.partitionCubeDataSelection(allselpars)
        #self.allimpars = self.PH.partitionCubeDeconvolution(allimagepars)

        # to define final image coordinates, run selecdata and definemage
        self.SItool = synthesisimager()
        # casalog.post("allselpars="+allselpars)
        origspw={}
        for mss in sorted( allselpars.keys() ): 
#            if(self.allimpars['0']['specmode']=='cubedata'):
#                self.allselpars[mss]['outframe']='Undefined'
            origspw[mss]={'spw':allselpars[mss]['spw']}
            self.SItool.selectdata( allselpars[mss] )
        for fid in sorted( allimagepars.keys() ):
            self.SItool.defineimage( allimagepars[fid], self.allgridpars[fid] )
            # insert coordsys record in imagepars 
            # partionCubeSelection works per field ...
            allimagepars[fid]['csys'] = self.SItool.getcsys()
            if allimagepars[fid]['nchan'] == -1:
                allimagepars[fid]['nchan'] = self.SItool.updatenchan()
            alldataimpars[fid] = self.PH.partitionCubeSelection(allselpars,allimagepars[fid])

        # if there are more nodes than nchan, there would be node(s)
        # that do not get any subcubes and causes an error. So
        # to avoid this, reduces the number of nodes actually uses.
        # Variable nchans among the fields are not supported yet in parallel mode
        fid0nchan = allimagepars['0']['nchan']
        nnodes = len(self.listOfNodes)
        if nnodes > fid0nchan:
           self.modifiedListOfNodes = self.listOfNodes[0:fid0nchan]
           casalog.post("Nchan = "+str(fid0nchan)+", Will use only "+str(len(self.modifiedListOfNodes))+" nodes out of "+str(nnodes), "WARN");
        else:
           self.modifiedListOfNodes = self.listOfNodes[:]

        #casalog.post("********************** " + alldataimpars.keys())
        #for kk in alldataimpars.keys():
        #    casalog.post("KEY : ", kk , " --->", alldataimpars[kk].keys())
            
        # reorganize allselpars and allimpars for partitioned data        
        synu = synthesisutils()
        self.allselpars={}
        self.allimpars={}
        ### casalog.post("self.listOfNodes=",self.listOfNodes)
        # Repack the data/image parameters per node
        #  - internally it stores zero-based node ids
        #  
        #for ipart in self.listOfNodes:
        for ipart in self.modifiedListOfNodes:
            # convert to zero-based indexing for nodes
            nodeidx = str(ipart-1)
            tnode = str(ipart)
            selparsPerNode= {tnode:{}}
            imparsPerNode= {tnode:{}}
            for fid in allimagepars:
                ###restoring original spw selection just to allow weight density to be the same
                ###ultimately should be passed by MPI if done this way
                for mss in origspw.keys():
                    alldataimpars[fid][nodeidx][mss]['spw']=origspw[mss]['spw']
                for ky in alldataimpars[fid][nodeidx]:
###                commenting this as it is resetting the selpars when key is not "msxxx" 
##                    selparsPerNode[tnode]={}
                    if ky.find('ms')==0:
                        # data sel per field
                        selparsPerNode[tnode][ky] = alldataimpars[fid][nodeidx][ky].copy();
                        if alldataimpars[fid][nodeidx][ky]['spw']=='-1':
                            selparsPerNode[tnode][ky]['spw']=''
                        #else:
                        ####using original spw selection for weight calculation
                        #    # remove chan selections (will be adjusted by tuneSelectData)
                        #   newspw=selparsPerNode[tnode][ky]['spw']
                        #  newspwlist = newspw.split(',')
                        #    spwsOnly = ''
                        #    for sp in newspwlist:
                        #        if spwsOnly!='': spwsOnly+=','
                        #        spwsOnly+=sp.split(':')[0]   
                        #       selparsPerNode[tnode][ky]['spw']=spwsOnly

                imparsPerNode[tnode][fid] = allimagepars[fid].copy()
                imparsPerNode[tnode][fid]['csys'] = alldataimpars[fid][nodeidx]['coordsys'].copy()
                imparsPerNode[tnode][fid]['nchan'] = alldataimpars[fid][nodeidx]['nchan']
##                imparsPerNode[tnode][fid]['imagename'] = imparsPerNode[tnode][fid]['imagename'] + '.n'+str(tnode) 
                imparsPerNode[tnode][fid]['imagename'] = self.PH.getpartimagename( imparsPerNode[tnode][fid]['imagename'], ipart )

                # skip this for now (it is not working properly, but should not affect results without this)
                #imparsPerNode[tnode][fid]=synu.updateimpars(imparsPerNode[tnode][fid])
            self.allselpars.update(selparsPerNode)
            self.allimpars.update(imparsPerNode)


            #casalog.post("****** SELPARS in init **********" + self.allselpars)
            #casalog.post("****** SELIMPARS in init **********" + self.allimpars)
        
        joblist=[]
        casa6_import_prefix = ''
        if is_CASA6:
            casa6_import_prefix = 'casatasks.private.'
        cmd_import_pars = ('from {0}imagerhelpers.input_parameters import ImagerParameters'.
                      format(casa6_import_prefix))
        cmd_import_synth = ('from {0}imagerhelpers.imager_base import PySynthesisImager'.
                     format(casa6_import_prefix))
        #### MPIInterface related changes
        #for node in range(0,self.NN):
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd(cmd_import_pars, node) )
            joblist.append( self.PH.runcmd(cmd_import_synth, node) )
        self.PH.checkJobs( joblist )

        self.exitflag={}
        joblist=[]
        #### MPIInterface related changes
        #for node in range(0,self.NN):
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("paramList = ImagerParameters()", node) )
            joblist.append( self.PH.runcmd("paramList.setSelPars("+str(self.allselpars[str(node)])+")", node) )
            joblist.append( self.PH.runcmd("paramList.setImagePars("+str(self.allimpars[str(node)])+")", node) )

            joblist.append( self.PH.runcmd("paramList.setGridPars("+str(self.allgridpars)+")", node) )
            joblist.append( self.PH.runcmd("paramList.setWeightPars("+str(self.weightpars)+")", node) )
            joblist.append( self.PH.runcmd("paramList.setDecPars("+str(self.decpars)+")", node) )
            joblist.append( self.PH.runcmd("paramList.setIterPars("+str(self.iterpars)+")", node) )
            joblist.append( self.PH.runcmd("paramList.setNormPars("+str(self.allnormpars)+")", node) )

            joblist.append( self.PH.runcmd("paramList.checkParameters()", node) )

            joblist.append( self.PH.runcmd("imager = PySynthesisImager(params=paramList)", node) )

            self.exitflag[str(node)] = False

        self.PH.checkJobs( joblist )

    def initializeImagers(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.initializeImagers()", node) )
        self.PH.checkJobs( joblist )

    def initializeDeconvolvers(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.initializeDeconvolvers()", node) )
        self.PH.checkJobs( joblist )

    def initializeNormalizers(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.initializeNormalizers()", node) )
        self.PH.checkJobs( joblist )

    def setWeighting(self):
        ## Set weight parameters and accumulate weight density (natural)
        joblist=[];
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            ## Set weighting pars
            joblist.append( self.PH.runcmd("imager.setWeighting()", node ) )
        self.PH.checkJobs( joblist )


    def initializeIterationControl(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.initializeIterationControl()", node) )
        self.PH.checkJobs( joblist )

    def makePSF(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.makePSF()", node) )
        self.PH.checkJobs( joblist )

    def runMajorMinorLoops(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.runMajorMinorLoops()", node) )
        self.PH.checkJobs( joblist )

    def runMajorCycle(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            if self.exitflag[str(node)]==False:
                joblist.append( self.PH.runcmd("imager.runMajorCycle()", node) )
        self.PH.checkJobs( joblist )

    def runMinorCycle(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            if self.exitflag[str(node)]==False:
                joblist.append( self.PH.runcmd("imager.runMinorCycle()", node) )
        self.PH.checkJobs( joblist )

    ## Merge the results from all pieces. Maintain an 'active' list of nodes...
    def hasConverged(self):

        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            if self.exitflag[str(node)]==False:
                joblist.append( self.PH.runcmd("rest = imager.hasConverged()", node) )
        self.PH.checkJobs( joblist )

#        self.PH.runcmdcheck("rest = imager.hasConverged()")

        retval = True
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            if self.exitflag[str(node)]==False:
                rest = self.PH.pullval("rest", node )
                retval = retval and rest[node]
                self.exitflag[str(node)] = rest[node]
                casalog.post("Node " + str(node) + " converged : " + str(rest[node]) , "INFO")

        return retval

    def updateMask(self):

        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            if self.exitflag[str(node)]==False:
                joblist.append( self.PH.runcmd("maskchanged = imager.updateMask()", node) )
        self.PH.checkJobs( joblist )

#        self.PH.runcmdcheck("maskchanged = imager.updateMask()")

        retval = False
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            if self.exitflag[str(node)]==False:
                rest = self.PH.pullval("maskchanged", node )
                retval = retval or rest[node]
                casalog.post("Node " + str(node) + " maskchanged : " + str(rest[node]) , "INFO")

        return retval

    def predictModel(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.predictmodel()", node) )
        self.PH.checkJobs( joblist )

    def restoreImages(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.restoreImages()", node) )
        self.PH.checkJobs( joblist )

    def pbcorImages(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.pbcorImages()", node) )
        self.PH.checkJobs( joblist )

    def makePB(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.makePB()", node) )
        self.PH.checkJobs( joblist )

    def checkPB(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.checkPB()", node) )
        self.PH.checkJobs( joblist )

    def concatImages(self, type='copyvirtual'):
        import subprocess
        imtypes=['image','psf','model','residual','mask','pb', 'image.pbcor', 'weight', 'sumwt']
        for immod in range(0,self.NF):
            for ext in imtypes:
                subimliststr="'"
                concatimname=self.allinimagepars[str(immod)]['imagename']+'.'+ ext
                distpath = os.getcwd()
                fullconcatimname = distpath+'/'+concatimname
                #for node in self.listOfNodes:
                for node in self.modifiedListOfNodes:
                    #rootimname=self.allinimagepars[str(immod)]['imagename']+'.n'+str(node)
                    #fullimname =  self.PH.getpath(node) + '/' + rootimname 
                    fullimname = self.PH.getpartimagename( self.allinimagepars[str(immod)]['imagename']  , node )
                    if (os.path.exists(fullimname+'.'+ext)):
                        subimliststr+=fullimname+'.'+ext+' '
                subimliststr+="'"
                if subimliststr!="''":
                    # parent images need to be cleaned up for restart=T
                    if self.allinimagepars[str(immod)]['restart'] and os.path.exists(fullconcatimname):
                        try:
                            casalog.post("Cleaning up the existing "+fullconcatimname,"DEBUG")
                            shutil.rmtree(fullconcatimname)
                        except:
                            casalog.post("Cleaning up the existing file named "+fullconcatimname,"DEBUG")
                            os.remove(fullconcatimname)
                    # set tempclose = false to avoid a long accessing issue
                    #cmd = 'imageconcat inimages='+subimliststr+' outimage='+"'"+fullconcatimname+"'"+' type='+type+' tempclose=false'      
                    #ret=os.system(cmd)
                    #if ret!=0:
                    #    casalog.post("concatenation of "+concatimname+" failed","WARN")
                    iatool=imageanalysis()
                    concattool = iatool.imageconcat(outfile=fullconcatimname, mode=type, infiles=subimliststr.strip("'"), axis=-1, tempclose=False, overwrite=True)
                    if(len(concattool.shape())==0):
                        casalog.post("concatenation of "+concatimname+" failed","WARN")
                    concattool.done()
                    

    def getSummary(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("summ = imager.getSummary("+str(node)+")", node) )
        self.PH.checkJobs( joblist )

        fullsumm={}
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
             summ = self.PH.pullval("summ", node )
             fullsumm["node"+str(node)] = summ

        return fullsumm

    def deleteTools(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.deleteTools()", node) )
        self.PH.checkJobs( joblist )

    def estimatememory(self):
        joblist=[]
        #for node in self.listOfNodes:
        for node in self.modifiedListOfNodes:
            joblist.append( self.PH.runcmd("imager.estimatememory()", node) )
        self.PH.checkJobs( joblist )
#############################################

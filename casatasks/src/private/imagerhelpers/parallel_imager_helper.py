from __future__ import absolute_import
import os
import math
import shutil
import string
import time
import re;
import copy

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import synthesisutils
    from casatasks import casalog
else:
    from taskinit import *

    synthesisutils = casac.synthesisutils

'''
A set of helper functions for the tasks  tclean

Summary...
    
'''

#############################################
###  Parallel Imager Helper.
#############################################
#casalog.post('Using clustermanager from MPIInterface', 'WARN')
try:
    if is_CASA6:
        from casampi.MPIInterface import MPIInterface as mpi_clustermanager
        mpi_available = True
    else:
        from mpi4casa.MPIInterface import MPIInterface as mpi_clustermanager
        mpi_available = True
except ImportError:
    mpi_available = False
    
class PyParallelImagerHelper():

    def __init__(self):

        ############### Cluster Info
         self.CL=None
         self.sc=None
         self.nodeList=None;
         # Initialize cluster, and partitioning.
        ############### Number of nodes to parallelize on

         # self.nodeList gets filled by setupCluster()
         self.NN = self.setupCluster()

    def getNodeList(self):
        return self.nodeList;

#############################################
    def chunkify(self,lst,n):
        return [ lst[i::n] for i in range(n) ]

    def partitionCFCacheList(self,gridPars):

        cfcName = gridPars['cfcache'];
        cflist=[];
        if (not (cfcName == '')):
            cflist=[f for f in os.listdir(cfcName) if re.match(r'CFS*', f)];
        nCF = len(cflist);
        nProcs=len(self.nodeList);
        
        if (nProcs > nCF):
            n=nCF;
        else:
            n=nProcs;
        if (nCF > 0):
            casalog.post("########################################################");
            casalog.post("nCF = " + str(nCF) + " nProcs = " + str(n) + " NodeList=" + str(self.nodeList));
            casalog.post("########################################################");
        xx=self.chunkify(cflist,n);
        allcfs={};
        for i in range(n):
            allcfs[i+1]=xx[i];

        return allcfs;
#############################################
# The above version works better (better balanced chunking).
# Keeping the code below in the file sometime, just in case...(SB).
    # def partitionCFCacheList(self,gridPars):

    #     cfcName = gridPars['cfcache'];
    #     cflist=[];
    #     if (not (cfcName == '')):
    #         cflist=[f for f in os.listdir(cfcName) if re.match(r'CFS*', f)];

    #     nCF = len(cflist);
    #     nProcs=len(self.nodeList);
        
    #     casalog.post("########################################################")
    #     casalog.post("nCF = " + nCF + " nProcs = " + nProcs + " NodeList=" + self.nodeList)
    #     casalog.post("########################################################")

    #     #n0=int(nCF/self.NN);
    #     n0=int(float(nCF)/nProcs+0.5);
    #     if (nProcs >= nCF):
    #         n0 = 1;
    #     allcfs = {};
    #     nUsed=0; i=1;
    #     while (nUsed < nCF):
    #         m = nUsed+n0;
    #         if (m > nCF): 
    #     	m=nCF;
    #         allcfs[i]=cflist[nUsed:m];
    #         nUsed = m;
    #         if (i >= nProcs):
    #             break;
    #         i=i+1;
    #     if (nUsed < nCF):
    #         allcfs[nProcs].append(cflist[i]);
    #     return allcfs;
            
#############################################
## Very rudimentary partitioning - only for tests. The actual code needs to go here.
    def partitionContDataSelection(self,oneselpars={}):

        synu = synthesisutils()
        allselpars =  synu.contdatapartition( oneselpars , self.NN )
        synu.done()

        casalog.post('Partitioned Selection : {}'.format(allselpars))
        return allselpars

#############################################
## Very rudimentary partitioning - only for tests. The actual code needs to go here.
    def partitionCubeDataSelection(self,oneselpars={}):

        synu = synthesisutils()
        allselpars =  synu.cubedatapartition( oneselpars , self.NN )
        synu.done()

        casalog.post('Partitioned Selection : {}'.format(allselpars))
        return allselpars

#############################################
    def partitionCubeDeconvolution(self,impars={}):

        synu = synthesisutils()
        allimpars =  synu.cubeimagepartition( impars , self.NN )
        synu.done()

        casalog.post('ImSplit : {}'.format(allimpars))
        return allimpars

#############################################
    def partitionCubeSelection(self, oneselpars={}, oneimpars={}):
        incsys = oneimpars['csys']
        nchan = oneimpars['nchan']
        synu = synthesisutils()
        allpars = synu.cubedataimagepartition(oneselpars, incsys, self.NN, nchan)
        synu.done()

        # casalog.post("Cube Data/Im partitioned selection : {}".format(allpars))
        return allpars

#############################################
    def setupCluster(self):
        # Initialize cluster

        # * Terminal: Client logs + Server logs
        # * casapy-<timestamp>.log: Client logs
        # * casapy-<timestamp>.log-server-<rank>-host-<hostname>-pid-<pid>: Server logs 
        mpi_clustermanager.set_log_mode('redirect');

        self.sc=mpi_clustermanager.getCluster()
        self.sc.set_log_level('DEBUG')

        self.CL=self.sc._cluster
        self.nodeList = self.CL.get_engines();
        numproc=len(self.CL.get_engines())
        numprocperhost=len(self.nodeList)/len(self.nodeList) if (len(self.nodeList) >0 ) else 1

        owd=os.getcwd()
        self.CL.pgc('import os')
        self.CL.pgc('from numpy import array,int32')
        self.CL.pgc('os.chdir("'+owd+'")')
        os.chdir(owd)
        casalog.post("setupCluster, Setting up {} engines.".format(numproc))
        return numproc

#############################################
    def takedownCluster(self):
        # Check that all nodes have returned, before stopping the cluster
         self.checkJobs()
         casalog.post('Ending use of cluster, but not closing it. Call clustermanager.stop_cluster() to close it if needed.')
#         self.sc.stop_cluster()
         self.CL=None
         self.sc=None

#############################################
    # This is a blocking call that will wait until jobs are done.
    def checkJobs(self,joblist=[]):
        #### MPIInterface related changes
        numcpu = len(self.nodeList)
        
        if len(joblist)==0:
             joblist = list(range(numcpu))
             #for k in range(numcpu):
             for k in self.nodeList:
                 joblist[k-1] = self.CL.odo('casalog.post("node '+str(k)+' has completed its job")', k)

        casalog.post('checkJobs. Blocking for nodes to finish')
        over=False
        while(not over):
            overone=True
            time.sleep(1)
            for k in range(len(joblist)):
                try:
                    overone =  self.CL.check_job(joblist[k],False) and overone
                except Exception:
                     raise
            over = overone
        casalog.post('...done')

#############################################
    def runcmd(self, cmdstr="", node=-1):
         if node >= 0:
              return self.CL.odo( cmdstr , node)
         else:
              self.CL.pgc( cmdstr )

#############################################
    def runcmdcheck(self, cmdstr):
         joblist=[]
         #### MPIInterface related changes
         #for node in range(0,self.NN):
         for node in self.nodeList:
              joblist.append( self.CL.odo( cmdstr, node ) )
         self.checkJobs( joblist )

#############################################
    def pullval(self, varname="", node=0):
         return self.CL.pull( varname , node )

#############################################
    def pushval(self, varname="", val=None, node=0):
         return self.CL.push( varname , val, node )

#############################################
    def getpath(self, node=0):
        enginepath = self.sc.get_engine_store(node)
        if enginepath==None:
            return ""
        else:
            return enginepath
#############################################
#    def deletepartimages(self, dirname, imname):
#        namelist = shutil.fnmatch.filter( os.listdir(dirname), imname+".*" )
#        #casalog.post("Deleting : " +  namelist + ' from ' + dirname +  ' starting with ' + imname)
#        for aname in namelist:
#            shutil.rmtree( dirname + "/" + aname )
#############################################
    def deletepartimages(self, imagename, node, deldir=False):
        namelist = shutil.fnmatch.filter( os.listdir(self.getworkdir(imagename, node)), "*" )
        #casalog.post("Deleting : " + namelist + ' from ' + self.getworkdir(imagename, node)  + ' starting with ' + imagename)
        for aname in namelist:
              shutil.rmtree( os.path.join(self.getworkdir(imagename, node), aname) )
        if deldir==True:
            #casalog.post("Deleting workdirectory : "+self.getworkdir(imagename, node))
            shutil.rmtree( self.getworkdir(imagename, node) )

#############################################
    def getworkdir(self, imagename, nodeid):
        workdir = ''
        workdir = os.path.join(self.getpath(nodeid), imagename + '.workdirectory')

        if( not os.path.exists(workdir) ):
            os.mkdir( workdir )

        return workdir
                                    
#############################################
    def getpartimagename(self, imagename, nodeid):
        """
        For imagename = 'imaging_subdir/foo_img', it produces something like:
        'imaging_subdir/foo_img.workdirectory/foo_img.n5.gridwt' (where n5 is the node idx)

        :param imagename: imagename as passed to the tclean task
        :param nodeid: id of MPI node

        :returns: (full path) name of a part/sub-image for nodeid, produced by concatenating
        the working directory, the image basename and the node id as a string.
        """
        # don't include subdirs again here - the workdir is already inside the subdir(s)
        image_basename = os.path.basename(imagename)
        return os.path.join(self.getworkdir(imagename,nodeid), image_basename + '.n' +
                            str(nodeid))


#############################################



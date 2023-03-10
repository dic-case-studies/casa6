#!/usr/bin/env python
from __future__ import absolute_import
import os
import sys
import copy
import shutil

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from .. import partitionhelper as ph
    from casatools import table as tbtool
    from casatools import ms as mstool
    from casatasks import casalog
    from casatasks.private.parallel.rflag_post_proc import combine_rflag_subreport, is_rflag_report
    from casatasks.private.parallel.rflag_post_proc import finalize_agg_rflag_thresholds
else:
    from parallel.rflag_post_proc import combine_rflag_subreport, is_rflag_report
    from parallel.rflag_post_proc import finalize_agg_rflag_thresholds
    import partitionhelper as ph
    from taskinit import *

# string.find (python 2) vs str_instance.find
if is_python3:
    def strfind(str_instance, a):
        return str_instance.find(a)
else:
    def strfind(str_instance, a):
        return string.find(str_instance,a)

# common function to use to get a dictionary values iterator
if is_python3:
    def locitervalues(adict):
        return adict.values()
else:
    def locitervalues(adict):
        return adict.itervalues()

# To handle thread-based Tier-2 parallelization
import threading
if not is_python3:
    import thread

# jagonzal (CAS-4106): Properly report all the exceptions and errors in the cluster framework
import traceback

# jagonzal (Migration to MPI)
try:
    if is_CASA6:
        from casampi.MPIEnvironment import MPIEnvironment
        from casampi.MPICommandClient import MPICommandClient
        mpi_available = True
    else:
        from mpi4casa.MPIEnvironment import MPIEnvironment
        from mpi4casa.MPICommandClient import MPICommandClient
        mpi_available = True
except ImportError:
    mpi_available = False

class JobData:
    """
    This class incapsulates a single job.  The commandName is the name
    of the task to be executed.  The jobInfo is a dictionary of all
    parameters that need to be handled.
    """
    class CommandInfo:

        def __init__(self, commandName, commandInfo, returnVariable):
            self.commandName = commandName
            self.commandInfo = commandInfo
            self.returnVariable = returnVariable

        def getReturnVariable(self):
            return self.returnVariable
        
        def getCommandLine(self):
            firstArgument = True
            output = "%s = %s(" % (self.returnVariable, self.commandName)
            for (arg,value) in self.commandInfo.items():
                if firstArgument:
                    firstArgument = False
                else:
                    output += ', '
                if isinstance(value, str):
                    output += ("%s = '%s'" % (arg, value))
                else:
                    output += ("%s = " % arg) + str(value)
            output += ')'
            return output
    
    
    def __init__(self, commandName, commandInfo = {}):
        self._commandList = []
        self.status  = 'new'
        self.addCommand(commandName, commandInfo)
        self._returnValues = None
            

    def addCommand(self, commandName, commandInfo):
        """
        Add an additional command to this Job to be exectued after
        previous Jobs.
        """
        rtnVar = "returnVar%d" % len(self._commandList)
        self._commandList.append(JobData.CommandInfo(commandName,
                                                     commandInfo,
                                                     rtnVar))
    def getCommandLine(self):
        """
        This method will return the command line(s) to be executed on the
        remote engine.  It is usually only needed for debugging or for
        the JobQueueManager.
        """
        output = ''
        for idx in range(len(self._commandList)):
            if idx > 0:
                output += '; '
            output += self._commandList[idx].getCommandLine()
        return output

    def getCommandNames(self):
        """
        This method will return a list of command names that are associated
        with this job.
        """
        return [command.commandName for command in self._commandList]
    

    def getCommandArguments(self, commandName = None):
        """
        This method will return the command arguments associated with a
        particular job.
           * If commandName is not none the arguments for the command with
             that name are returned.
           * Otherwise a dictionary (with keys being the commandName and
             the value being the dictionary of arguments) is returned.
           * If there is only a single command the arguments for that
             command are returned as a dictionary.
        """
        returnValue = {}
        for command in self._commandList:
            if commandName is None or commandName == command.commandName:
                returnValue[command.commandName] = command.commandInfo
                                                   
        if len(returnValue) == 1:
            return list(returnValue.values())[0]
        return returnValue
    
    def getReturnVariableList(self):
        return [ci.returnVariable for ci in self._commandList]

    def setReturnValues(self, valueList):
        self._returnValues = valueList

    def getReturnValues(self):
        if self._returnValues is not None:
            if len(self._returnValues) == 1:
                return self._returnValues[0]
        return self._returnValues

class ParallelTaskHelper:
    """
    This is the extension of the TaskHelper to allow for parallel
    operation.  For simple tasks all that should be required to make
    a task parallel is to use this rather than the TaskHelper method
    above
    """

    __bypass_parallel_processing = 0
    __async_mode = False
    __multithreading = False    
    
    def __init__(self, task_name, args = {}):
        self._arg = dict(args)
        self._arguser = {}
        self._taskName = task_name
        self._executionList = []
        self._jobQueue = None
        # Cache the initial inputs
        self.__originalParams = args
        # jagonzal: Add reference to cluster object
        self._cluster = None
        self._mpi_cluster = False
        self._command_request_id_list = None
        if not mpi_available or not MPIEnvironment.is_mpi_enabled:
            self.__bypass_parallel_processing = 1
        if (self.__bypass_parallel_processing == 0):
            self._mpi_cluster = True
            self._command_request_id_list = []
            self._cluster = MPICommandClient()
        # jagonzal: To inhibit return values consolidation
        self._consolidateOutput = True
        # jagonzal (CAS-4287): Add a cluster-less mode to by-pass parallel processing for MMSs as requested
        # This is actually a dict, with key=vis and value= the 'success' field of the cmd.
        # (exception: for tasks with parameter outputvis (like partition), key=outputvis)
        self._sequential_return_list = {}
        
    def override_arg(self,arg,value):
        self._arguser[arg] = value

    def initialize(self):
        """
        This is the setup portion.
        Currently it:
           * Finds the full path for the input vis.
           * Initialize the MPICommandClient
        """
        self._arg['vis'] = os.path.abspath(self._arg['vis'])
        
        # jagonzal (Migration to MPI)
        if self._mpi_cluster:
            self._cluster.start_services()
            
    def getNumberOfServers(self):
        """
        Return the number of engines (iPython cluster) or the number of servers (MPI cluster)
        """
        if (mpi_available and self.__bypass_parallel_processing == 0):
            return len(MPIEnvironment.mpi_server_rank_list())
        else:
            return None

    def generateJobs(self):
        """
        This is the method which generates all of the actual jobs to be
        done.  The default is to assume the input vis is a reference ms and
        build one job for each referenced ms.
        """
        
        casalog.origin("ParallelTaskHelper")
        
        try:
            msTool = mstool()
            if not msTool.open(self._arg['vis']):
                raise ValueError("Unable to open MS %s," % self._arg['vis'])
            if not msTool.ismultims():
                raise ValueError("MS is not a MultiMS, simple parallelization failed")

            subMs_idx = 0
            for subMS in msTool.getreferencedtables():
                localArgs = copy.deepcopy(self._arg)
                localArgs['vis'] = subMS
                
                for key in self._arguser:
                    localArgs[key] = self._arguser[key][subMs_idx]
                subMs_idx += 1
                
                if self._mpi_cluster:
                    self._executionList.append([self._taskName + '()',localArgs])
                else:
                    self._executionList.append(JobData(self._taskName,localArgs))
                
            msTool.close()
            return True
        except Exception as instance:
            casalog.post("Error handling MMS %s: %s" % (self._arg['vis'],instance),"WARN","generateJobs")
            msTool.close()
            return False


    def executeJobs(self):
        
        casalog.origin("ParallelTaskHelper")
        
        # jagonzal (CAS-4287): Add a cluster-less mode to by-pass parallel processing for MMSs as requested 
        if (self.__bypass_parallel_processing == 1):
            for job in self._executionList:
                parameters = job.getCommandArguments()
                try:
                    if is_CASA6:
                        gvars = globals( )
                        try:
                            exec("from casatasks import *; " + job.getCommandLine(),gvars)
                        except Exception as exc:
                            casalog.post("exec in parallel_task_helper.executeJobs failed: {}'".format(exc))
                            raise

                        # jagonzal: Special case for partition
                        # The 'True' values emulate the command_response['successful'] that
                        # we'd get in parallel runs from other MPI processes.
                        if 'outputvis' in parameters:
                            self._sequential_return_list[parameters['outputvis']] = True
                        else:
                            self._sequential_return_list[parameters['vis']] = gvars['returnVar0'] or True
                    else:
                        exec("from taskinit import *; from tasks import *; " + job.getCommandLine())

                        # jagonzal: Special case for partition
                        if ('outputvis' in parameters):
                            self._sequential_return_list[parameters['outputvis']] = True
                        else:
                            self._sequential_return_list[parameters['vis']] = returnVar0 or True

                except Exception as instance:
                    str_instance = str(instance)
                    if (strfind(str_instance, "NullSelection") == 0):
                        casalog.post("Error running task sequentially %s: %s" % (job.getCommandLine(),str_instance),"WARN","executeJobs")
                        traceback.print_tb(sys.exc_info()[2])
                    else:
                        casalog.post("Ignoring NullSelection error from %s" % (parameters['vis']),"INFO","executeJobs")
            self._executionList = []
        else:
            for job in self._executionList:
                command_request_id = self._cluster.push_command_request(job[0],False,None,job[1])
                self._command_request_id_list.append(command_request_id[0])


    def postExecution(self):   

        casalog.origin("ParallelTaskHelper")

        ret_list = {}
        if (self.__bypass_parallel_processing==1):
            ret_list = self._sequential_return_list
            self._sequential_return_list = {}        
        elif (self._cluster != None):
            # jagonzal (CAS-7631): Support for thread-based Tier-2 parallelization
            if ParallelTaskHelper.getMultithreadingMode():
                event = self._cluster.get_command_response_event(self._command_request_id_list)
                ParallelTaskWorker.releaseTaskLock()
                event.wait()
                ParallelTaskWorker.acquireTaskLock()
            # Get command response
            command_response_list =  self._cluster.get_command_response(self._command_request_id_list,True,True)
            # Format list in the form of vis dict
            ret_list = {}
            for command_response in command_response_list:
                vis = command_response['parameters']['vis']
                if 'uvcontsub' in command_response['command']:
                    # One more particular case, similar as in 'executeJob' for partition.
                    # The design of these lists and how they are used in different ways in
                    # tasks uvcontsub, setjy, flagdata, etc. is evil
                    # uvcontsub expects a 'success' True/False value for every subMS rather
                    # than the return value of the subMS uvcontsub.
                    ret_list[vis] = command_response['successful']
                else:
                    ret_list[vis] = command_response['ret']
        else:
            return None

        ret = ret_list
        if self._consolidateOutput:
            ret = ParallelTaskHelper.consolidateResults(ret_list,self._taskName)

        return ret


    @staticmethod
    def consolidateResults(ret_list,taskname):
        if isinstance(list(ret_list.values())[0],bool):
            retval = True
            for subMs in ret_list:
                if not ret_list[subMs]:
                    casalog.post("%s failed for sub-MS %s" % (taskname,subMs),"WARN","consolidateResults")
                    retval = False
            return retval
        elif any(isinstance(v,dict) for v in locitervalues(ret_list)):
            ret_dict = {}
            for _key, subMS_dict in ret_list.items():
                casalog.post(" ***** consolidateResults, subMS: {0}".format(subMS_dict),
                             "WARN", "consolidateResults")
                if isinstance(subMS_dict, dict):
                    try:
                        ret_dict = ParallelTaskHelper.combine_dictionaries(subMS_dict, ret_dict)
                    except Exception as instance:
                        casalog.post("Error post processing MMS results {0}: {1}".format(
                            subMS_dict, instance), 'WARN', 'consolidateResults')
                        raise
            return ParallelTaskHelper.finalize_consolidate_results(ret_dict)


    @staticmethod
    def combine_dictionaries(dict_list,ret_dict):
        """
        Combines a flagging (sub-)report dictionary dict_list (from a subMS) into an overall
        report dictionary (ret_dict).
        """
        for key, item in dict_list.items():
            if isinstance(item, dict):
                if key in ret_dict:
                    if is_rflag_report(item):
                        ret_dict[key] = combine_rflag_subreport(item, ret_dict[key])
                    else:
                        ret_dict[key] = ParallelTaskHelper.combine_dictionaries(item,ret_dict[key])
                else:
                    ret_dict[key] = ParallelTaskHelper.combine_dictionaries(item,{})
            else:
                if key in ret_dict:
                    # the 'nreport' field should not be summed - it's an index
                    if not isinstance(ret_dict[key],str) and 'nreport' != key:
                        # This is a good default for all reports that have flag counters
                        ret_dict[key] += item
                else:
                    ret_dict[key] = item

        return ret_dict


    @staticmethod
    def finalize_consolidate_results(ret):
        """ Applies final step to the items of the report dictionary.
        For now only needs specific processing to finalize the aggregation of the RFlag
        thresholds (freqdev/timedev) vectors. """

        for key, item in ret.items():
            if isinstance(item, dict) and is_rflag_report(item):
                ret[key] = finalize_agg_rflag_thresholds(item)

        return ret


    @staticmethod
    def getResult(command_request_id_list,taskname):
        
        # Access MPICommandClietn singleton instance
        client = MPICommandClient()
        
        # Get response list
        command_response_list =  client.get_command_response(command_request_id_list,True,True)
                
        # Format list in the form of vis dict
        ret_list = {}
        for command_response in command_response_list:
            vis = command_response['parameters']['vis']
            ret_list[vis] = command_response['ret']
            
        # Consolidate results and return
        ret = ParallelTaskHelper.consolidateResults(ret_list,taskname)
        
        return ret                    


    def go(self):

        casalog.origin("ParallelTaskHelper")

        self.initialize()
        if (self.generateJobs()):
            self.executeJobs()

            if ParallelTaskHelper.__async_mode:
                res_list = [] if self._command_request_id_list is None else list(self._command_request_id_list)
                return res_list
            else:
                try:
                    retVar = self.postExecution()
                except Exception as instance:
                    casalog.post("Error post processing MMS results %s: %s" % (self._arg['vis'],instance),"WARN","go")
                    traceback.print_tb(sys.exc_info()[2])
                    return False
        else:
            retVar = False

        # Restore casalog origin
        casalog.origin(self._taskName)

        return retVar

    @staticmethod
    def getReferencedMSs(vis):
        
        msTool = mstool()
        if not msTool.open(vis):
            raise ValueError("Unable to open MS %s." % vis)

        if not msTool.ismultims():
            raise ValueError("MS %s is not a reference MS." % vis)

        rtnValue = msTool.getreferencedtables()
        if not isinstance(rtnValue, list):
            rtnValue = [rtnValue]
      
        msTool.close()
        return rtnValue


    @staticmethod
    def restoreSubtableAgreement(vis, mastersubms='', subtables=[]):
        """
        Tidy up the MMS vis by replacing the subtables of all SubMSs
        by the subtables from the SubMS given by "mastersubms".
        If specified, only the subtables in the list "subtables"
        are replaced, otherwise all.
        If "mastersubms" is not given, the first SubMS of the MMS
        will be used as master.
        """

        msTool = mstool();
        msTool.open(vis)
        theSubMSs = msTool.getreferencedtables()
        msTool.close()

        tbTool = tbtool( );
        
        if mastersubms=='':
            tbTool.open(vis)
            myKeyw = tbTool.getkeywords()
            tbTool.close()
            mastersubms=os.path.dirname(myKeyw['ANTENNA'].split(' ')[1]) #assume ANTENNA is present

        mastersubms = os.path.abspath(mastersubms)
            
        theSubTables = ph.getSubtables(mastersubms)

        if subtables==[]:
            subtables=theSubTables
        else:
            for s in subtables:
                if not (s in theSubTables):
                    raise ValueError( s+' is not a subtable of '+ mastersubms )

        origpath = os.getcwd()      
        masterbase = os.path.basename(mastersubms)
        
        for r in theSubMSs:
            rbase = os.path.basename(r)
            if not rbase==masterbase:
                for s in subtables:
                    theSubTab = r+'/'+s
                    if os.path.islink(theSubTab): # don't copy over links
                        if(os.path.basename(os.path.dirname(os.path.realpath(theSubTab)))!=masterbase):
                            # the mastersubms has changed: make new link
                            os.chdir(r)
                            shutil.rmtree(s, ignore_errors=True)
                            os.symlink('../'+masterbase+'/'+s, s)
                            os.chdir(origpath)
                    else:    
                        shutil.rmtree(theSubTab, ignore_errors=True)
                        shutil.copytree(mastersubms+'/'+s, theSubTab)

        return True

    @staticmethod
    def bypassParallelProcessing(switch=1):
        """
        # jagonzal (CAS-4287): Add a cluster-less mode to by-pass parallel processing for MMSs as requested 
        switch=1 => Process each sub-Ms sequentially
        switch=2 => Process the MMS as a normal MS
        """        
        ParallelTaskHelper.__bypass_parallel_processing = switch
        
    @staticmethod
    def getBypassParallelProcessing():
        """
        # jagonzal (CAS-4287): Add a cluster-less mode to by-pass parallel processing for MMSs as requested 
        switch=1 => Process each sub-Ms sequentially
        switch=2 => Process the MMS as a normal MS
        """        
        return ParallelTaskHelper.__bypass_parallel_processing        
    
    @staticmethod
    def setAsyncMode(async_mode=False):     
        ParallelTaskHelper.__async_mode = async_mode
        
    @staticmethod
    def getAsyncMode():
        return ParallelTaskHelper.__async_mode    
    
    @staticmethod
    def setMultithreadingMode(multithreading=False):     
        ParallelTaskHelper.__multithreading = multithreading
        
    @staticmethod
    def getMultithreadingMode():
        return ParallelTaskHelper.__multithreading
    
    @staticmethod
    def isParallelMS(vis):
        """
        This method will let us know if we can do the simple form
        of parallelization by invoking on many referenced mss.
        """
        
        # jagonzal (CAS-4287): Add a cluster-less mode to by-pass parallel processing for MMSs as requested 
        if (ParallelTaskHelper.__bypass_parallel_processing == 2):
            return False
        
        msTool = mstool()
        if not msTool.open(vis):
            raise ValueError( "Unable to open MS %s," % vis)
        rtnVal = msTool.ismultims() and \
                 isinstance(msTool.getreferencedtables(), list)

        msTool.close()
        return rtnVal
    
    @staticmethod
    def findAbsPath(input):
        if isinstance(input,str):
            return os.path.abspath(input)

        if isinstance(input, list):
            rtnValue = []
            for file_i in input:
                rtnValue.append(os.path.abspath(file_i))
            return rtnValue

        # Your on your own, don't know what to do
        return input

    @staticmethod
    def isMPIEnabled():
        return MPIEnvironment.is_mpi_enabled if mpi_available else False

    @staticmethod
    def isMPIClient():
        return MPIEnvironment.is_mpi_client if mpi_available else False

    @staticmethod
    def listToCasaString(inputList):
        """
        This Method will take a list of integers and try to express them as a 
        compact set using the CASA notation.
        """
        if inputList is None or len(inputList) == 0:
            return ''
        
        def selectionString(rangeStart, rangeEnd):
            if rangeStart == rangeEnd:
                return str(rangeStart)
            return "%d~%d" % (rangeStart, rangeEnd)
    
        inputList.sort()
        compactStrings = []
        rangeStart = inputList[0]
        lastValue = inputList[0]
        for val in inputList[1:]:
            if val > lastValue + 1:
                compactStrings.append(selectionString(rangeStart,lastValue))
                rangeStart = val
            lastValue = val
        compactStrings.append(selectionString(rangeStart,lastValue))

        return ','.join([a for a in compactStrings])
    

class ParallelTaskWorker:
    
    # Initialize task lock
    __task_lock = threading.Lock()
    
    def __init__(self, cmd):
        
        self.__cmd = compile(cmd,"ParallelTaskWorker", "eval")
        self.__state = "initialized"
        self.__res = None        
        self.__thread = None
        self.__environment = self.getEnvironment()
        self.__formatted_traceback = None        
        self.__completion_event = threading.Event()  

    def getEnvironment(self):
        try:
            # casampi should not depend on globals (casashell). And CASA6/casashell doesn't
            # anyway have init_tasks:update_params. Keep going w/o globals
            import casampi
            return {}
        except ImportError:
            stack=inspect.stack()
            for stack_level in range(len(stack)):
                frame_globals=sys._getframe(stack_level).f_globals
                if 'update_params' in frame_globals:
                    return dict(frame_globals)

            raise Exception("CASA top level environment not found")
        
    def start(self):
        
        # Initialize completion event
        self.__completion_event.clear()        
               
        # Spawn thread
        if is_python3:
            self.__thread = threading.Thread(target=self.runCmd, args=(), kwargs=())
            self.__thread.setDaemon(True)
            self.__thread.start()
        else:
            self.__thread = thread.start_new_thread(self.runCmd, ())

        # Mark state as running
        self.__state = "running"        

    def runCmd(self):
        
        # Acquire lock
        ParallelTaskWorker.acquireTaskLock()
        
        # Update environment with globals from calling context
        globals().update(self.__environment)
        
        # Run compiled command
        try:
            self.__res = eval(self.__cmd)
            # Mark state as successful
            self.__state = "successful"
            # Release task lock
            ParallelTaskWorker.releaseTaskLock()            
        except Exception as instance:
            # Mark state as failed
            self.__state = "failed"
            # Release task lock if necessary
            if ParallelTaskWorker.checkTaskLock():ParallelTaskWorker.releaseTaskLock()
            # Post error message
            self.__formatted_traceback = traceback.format_exc()
            casalog.post("Exception executing command '%s': %s" 
                         % (self.__cmd,self.__formatted_traceback),
                         "SEVERE","ParallelTaskWorker::runCmd")
        
        # Send completion event signal
        self.__completion_event.set()
        
    def getResult(self):
        
        if self.__state == "running":
            # Wait until completion event signal is received
            self.__completion_event.wait()
            
            
        if self.__state == "initialized":
            casalog.post("Worker not started",
                         "WARN","ParallelTaskWorker::getResult")
        elif self.__state == "successful":
            return self.__res            
        elif self.__state == "failed":
            casalog.post("Exception executing command '%s': %s" 
                         % (self.__cmd,self.__formatted_traceback),
                         "SEVERE","ParallelTaskWorker::runCmd")                        
    
    @staticmethod
    def acquireTaskLock():
        
        ParallelTaskWorker.__task_lock.acquire()
        
    @staticmethod
    def releaseTaskLock():
        
        ParallelTaskWorker.__task_lock.release()
        
    @staticmethod
    def checkTaskLock():
        
        return ParallelTaskWorker.__task_lock.locked()        
          
          

#from .compare import *
#from .testhelper import *
#from .extractcasascript import main
#from .testhelpers import TestHelpers
#import imagehelpers.imagetesthelpers
#from imagehelpers import imagetesthelpers

import os
import sys
import time
from functools import wraps
import fnmatch
import logging
import filecmp
import unittest
import pickle
import numbers
import operator
import subprocess
import numpy


_casa5 = False
_casa6 = False
_importmpi = False
__bypass_parallel_processing = 0

# https://stackoverflow.com/questions/52580105/exception-similar-to-modulenotfounderror-in-python-2-7
try:
    ModuleNotFoundError
except NameError:
    ModuleNotFoundError = ImportError

def import_casamods():
    try:
        # CASA 6
        logging.debug("Importing CASAtools")
        import casatools
        logging.debug("Importing CASAtasks")
        #try:
        #    import casatasks
        #    from casatasks import casalog
        #except (ImportError, ModuleNotFoundError):
        #    pass

        try:
            from casampi.MPIEnvironment import MPIEnvironment
            _importmpi = True
            if not MPIEnvironment.is_mpi_enabled:
                __bypass_parallel_processing = 1
        except ImportError:
            print("MPIEnvironment not Enabled")

        _casa6 = True

    except (ImportError, ModuleNotFoundError):
        # CASA 5
        logging.debug("Import casa6 errors. Trying casa5...")
        from __main__ import default
        from taskinit import tbtool, mstool, iatool
        from casa_stack_manip import stack_find, find_casa

        try:
            from mpi4casa.MPIEnvironment import MPIEnvironment
            _importmpi = True
            if not MPIEnvironment.is_mpi_enabled:
                __bypass_parallel_processing = 1
        except ImportError:
            print("MPIEnvironment not Enabled")

        casa = find_casa()
        if casa.has_key('state') and casa['state'].has_key('init_version') and casa['state']['init_version'] > 0:
            casaglobals=True
            casac = stack_find("casac")
            #casalog = stack_find("casalog")
        _casa5 = True

_casa6tools = set([
    "agentflagger", "atcafiller", "atmosphere", "calanalysis", "calibrater", "coercetype", "componentlist", "config", "constants", "coordsys", "ctuser", "functional", "image",
    "imagemetadata", "imagepol", "imager", "iterbotsink", "logsink", "measures", "miriadfiller", "ms", "msmetadata", "mstransformer", "platform", "quanta", "regionmanager", "sakura",
    "sdm", "simulator", "singledishms", "spectralline", "synthesisdeconvolver", "synthesisimager", "synthesisimstore", "synthesisnormalizer", "synthesisutils", "table", "typecheck", "utils",
    "vlafiller", "vpmanager"])


_casa6tasks = set([
    "accor", "accum", "applycal", "asdmsummary", "bandpass", "blcal", "calstat", "clearcal", "clearstat", "concat", "conjugatevis", "cvel", "cvel2",
    "delmod" ,"exportasdm", "exportfits", "exportuvfits", "feather", "fixplanets", "fixvis", "flagcmd", "flagdata", "flagmanager", "fluxscale", "ft", "gaincal",
    "gencal", "hanningsmooth", "imcollapse", "imcontsub", "imdev", "imfit", "imhead", "imhistory", "immath", "immoments", "impbcor", "importasap", "importasdm",
    "importatca", "importfits", "importfitsidi", "importgmrt", "importmiriad", "importnro", "importuvfits", "importvla", "impv", "imrebin", "imreframe",
    "imregrid", "imsmooth", "imstat", "imsubimage", "imtrans", "imval", "initweights", "listcal", "listfits", "listhistory", "listobs", "listpartition",
    "listsdm", "listvis", "makemask", "mstransform", "partition", "polcal", 'polfromgain', "predictcomp", "rerefant", "rmfit", "rmtables", "sdbaseline", "sdcal",
    "sdfit", "sdfixscan", "sdgaincal", "sdimaging", "sdsmooth", "setjy", "simalma", "simanalyze", "simobserve", "slsearch", "smoothcal", "specfit",
    "specflux", "specsmooth", "splattotable", "split", "spxfit", "statwt", "tclean", "uvcontsub", "uvmodelfit", "uvsub", "virtualconcat", "vishead", "visstat", "widebandpbcor","deconvolve"])

_miscellaneous_tasks = set(['wvrgcal','plotms'])

 
############################################################################################
##################################       General Functions       ###########################
############################################################################################


def getNumberOfServers( __bypass_parallel_processing ):
    """
    Return the number of engines (iPython cluster) or the number of servers (MPI cluster)
    """
    import_casamods()
    if (__bypass_parallel_processing == 0) and (_importmpi):
        return len(MPIEnvironment.mpi_server_rank_list()) 
    else:
        return None

def add_to_dict(self, output=None, dataset="TestData", status=False, **kwargs):
    '''
        This function adds key value pairs to a provided dictionary. Any additional keys and values can be added as keyword arguments to this function
        @param output: This is the dictionary that the key-value pairs will be appended to
        @param filename: This is the name of the test script file
        @param dataset: This is the name of the dataset used when executing this test case
        @return: Nothing is returned, the output dict is modified by this function
    '''
    import inspect
    import multiprocessing
    frame = inspect.stack()[1]
    #print(frame)
    module = inspect.getmodule(frame[0])
    filename = module.__file__
    #print(filename)
    testcase = unittest.TestCase.id(self)
    #print(testcase)
    test_split = testcase.split('.')
    test_case = test_split[-1]
    #taskname = test_split[1].split('_')[0]
    #print(taskname)
    if (sys.version_info > (3, 3)):
        try:
            casapath = os.environ.get('CASAPATH').split()[0] + '/bin/'
        except:
            casapath = ''
        rerun = "{}python3 {} {}.{}".format(casapath,filename, test_split[1], test_split[2])
       
    else:
        filename = "{}.py".format(filename.split('.')[0])
        casapath = os.environ.get('CASAPATH').split()[0]
        rerun = "{}/bin/casa -c {}/lib/python2.7/runUnitTest.py {}".format(casapath,casapath, filename.split('.')[0])
    current_case = None
    func_calls = []
    values = {key:kwargs[key] for key in kwargs}
    with open(filename, 'r') as file:
        for line in file:
            #print(line)
            line = line.strip()
            if line.startswith('def test_'):
                if line.split()[1][:-7].endswith(test_case):
                    current_case = test_case
                else:
                    current_case = None
                    #_casa6tasks, _miscellaneous_tasks
            for task in _casa6tasks.union(_miscellaneous_tasks):
                if current_case == test_case:
                    if "{}(".format(task) in line:
                        #print(line)
                        taskname = line.split("(")[0]
                        ## Optional: Can print first index of casa function call but does not print the string name if it's an assigned object
                        ## Attempt to Get Dataset from casa task call
                        if dataset== "TestData":
                            import re
                            dataset = re.search('(?<=\().+?(?=\,)',line).group()
                            if len(dataset) == 0 or dataset is None:
                                dataset = "TestData"
                        params = line.split(',')[1::]
                        #print(params)
                        while ')' not in list(line):
                            line = next(file)
                            new_line = [x.strip() for x in line.split(',')]
                            for i in new_line:
                                params.append(i)
                            params = list(filter(lambda a: a != '', params))
                        call = "{}({},{}".format(taskname, dataset, ','.join(params))
                        #print(call)
                        func_calls.append(call)
                        #print(func_calls)
    values['runtime'] = -1.0
    #This is a temp error value
    values['status'] = status
    if test_case not in output.keys():
        output[test_case]= {}
    for key in values.keys():
        if test_case in output.keys():
            #print("output[test_case].keys(): {}".format(output[test_case].keys()))
            if key in output[test_case].keys():
                values[key] = output[test_case][key].append(values[key])
            else:
                output[test_case][key] = [values[key]]
    #output[test_case] = values
    output[test_case]['taskcall'] = func_calls
    output[test_case]['rerun'] = rerun
    output[test_case]['description'] = unittest.TestCase.shortDescription(self)
    output[test_case]['images'] = [ ]
    if getNumberOfServers(__bypass_parallel_processing) == None:
        output[test_case]['Serial Mode'] = "MPI Environment Not Enabled"
    else: 
        output[test_case]['MPI Mode'] = "{}".format(str(str(getNumberOfServers(__bypass_parallel_processing)) + " MPI Servers + 1 Client"))

    output[test_case]['Number of processors available in the system'] = "{}".format(str(multiprocessing.cpu_count()))


    #print("Test Case: {}".format(test_case))
    #print("{} : {}".format(test_case,output[test_case]))

def to_pickle(input_dict, picklefile):
    '''
        Add a new dictionary into the existing pickle file
        @param input_dict: The dictionary object to add to the pickle file
        @param picklefile: The picklefile containing a dictionary to be appended to
        @return: Nothing is returned by this function
    '''
    pickle_read = open(picklefile, 'rb')
    pickle_dict = pickle.load(pickle_read)
    # Make sure that the pickle file contains a dictionary
    if type(pickle_dict) != type({}):
        logging.warning('The pickle file is not a dictionary')
    # Add to the dictionary in the pickle file
    for item in list(input_dict.keys()):
        pickle_dict[item] = input_dict[item]
    # Re-write the pickle file with the new dictionary
    with open(picklefile, 'wb') as fout:
        pickle.dump(pickle_dict, fout)

def generate_weblog(task,dictionary,show_passed = True):
    """Generate Test Summary Weblog
    Example:
        generate_weblog("taskname", dictionary, show_passed)
    """
    import_casamods()
    from .weblog import Weblog
    Weblog(task, dictionary).generate_weblog(show_passed = show_passed)

############################################################################################
##################################       Decorators       ##################################
############################################################################################

#import casatestutils
#@casatestutils.skipIfMissingModule
def skipIfMissingModule(required_module, strict=False):
    '''
    Decorator: skip test if specified module is not avaliable
    Example:
        @casatestutils.skipIfMissingModule('astropy')
        def test_test(self):
    '''
    try:
        __import__(required_module)
        flag = True
    except ImportError:
        flag = False
    def deco(function):
        if not _casa6:
            return deco
        def wrapper(self, *args, **kwargs):
            if not flag:
                # If there is a strict flag run the tests as normal
                print(sys.argv)
                if strict:
                    function(self)
                else:
                    # Module ImportError and no strict flag
                    self.skipTest("ModuleNotFoundError: No module named '{}'".format(required_module))
            else:
                function(self)
        return wrapper
    return deco

#import casatestutils
#@casatestutils.time_execution
def time_execution(out_dict):
    import_casamods()
    def time_decorator(function):
        '''
        Decorator: time execution of test
        Example:
            @casatestutils.time_execution
            def test_test(self):
        '''
        @wraps(function)
        def function_timer(*args, **kwargs):
            failed = False
            result = None
            out_dict[function.__name__] = {}
            out_dict[function.__name__]['description'] = "An unexpected exception was raised"
            out_dict[function.__name__]['taskcall'] = ["An unexpected exception was raised by the test case"]
            out_dict[function.__name__]['rerun'] = "An unexpected exception was raised by the test case"
            t0 = time.time()
            try:
                result = function(*args, **kwargs)
            except Exception as e:
                failed = True
                t1 = time.time()
                out_dict[function.__name__]['runtime'] = t1-t0
                #casalog.post("Total time running {}: {} seconds".format(function.__name__, str(t1-t0)))
                out_dict[function.__name__]['status'] = False
                out_dict[function.__name__]['Failure Message'] = e
                raise
            t1 = time.time()
            #print ("Total time running %s: %s seconds" % (function.__name__, str(t1-t0)))
            #casalog.post("Total time running {}: {} seconds".format(function.__name__, str(t1-t0)))
            #print('======================================================')
            #print(function.__name__)
            out_dict[function.__name__]['runtime'] = t1-t0
            out_dict[function.__name__]['status'] = True
            return result
        return function_timer
    return time_decorator

def cpu_usage(out_dict):
    def cpu_decorator(function):
        @wraps(function)
        def function_usage(*args, **kwargs):
            #Temp Fix : CASA 5 Doesnt Have psutil by default
            try:
                import psutil
                use_psutil = True
            except ImportError:
                use_psutil = False
            if use_psutil:
                process = psutil.Process(os.getpid())
                snapshot1 = process.memory_info()
                open_files1 = process.open_files()
                num_file_descriptors1 = process.num_fds()
                #print ("Function: {}, {} MBs".format(function.__name__, megs1))
                #print ("Function: {}, Open Files: {}".format(function.__name__, open_files1))
                #print ("Function: {}, num_file_descriptors: {}".format(function.__name__, num_file_descriptors1))
                result = function(*args, **kwargs)
                process = psutil.Process(os.getpid())
                snapshot2 = process.memory_info()
                open_files2 = process.open_files()
                num_file_descriptors2 = process.num_fds()
                #print ("Function: {}, {} MBs".format(function.__name__, megs2))
                #print ("Function: {}, Open Files: {}".format(function.__name__, open_files2))
                #print ("Function: {}, num_file_descriptors: {}".format(function.__name__, num_file_descriptors2))
                #print('{:.2f} MB\n'.format(process.memory_info().rss / 1024 / 1024))
                #print ("Total Mem Info { }: {:.2f} MB".format(function.__name__,(process.memory_info().rss) / 1024 / 1024 ))
                out_dict[function.__name__]['cpu_usage'] = {"number of file descriptors opened" :  num_file_descriptors2 - num_file_descriptors1,
                                                            "Open files" : open_files2,
                                                            "Pre Memory Snapshot (bytes)" : snapshot1,
                                                            "Post Memory Snapshot (bytes)" : snapshot2
                                                           }
            else:
                #TODO: Add methods to get mem snapshots when psutils is not available
                result = function(*args, **kwargs)
                out_dict[function.__name__]['cpu_usage'] = {"number of file descriptors opened" : "Unknown",
                                                            "Open files" : "Unknown",
                                                            "Pre Memory Snapshot (bytes)" : "Unknown",
                                                            "Post Memory Snapshot (bytes)" : "Unknown"
                                                           }
            return result
        return function_usage
    return cpu_decorator

def peak_mem(out_dict):
    #TODO: https://pytracemalloc.readthedocs.io/examples.html
    ### NOTE: Only for python3.4+
    def mem_decorator(function):
        @wraps(function)
        def function_mem(*args, **kwargs):
            if sys.version_info > (3, 3):
                import tracemalloc
                tracemalloc.clear_traces()
                tracemalloc.start()
                snapshot1 = tracemalloc.take_snapshot() # Snapshot of traces of memory blocks allocated by Python.
                result = function(*args, **kwargs)
                snapshot2 = tracemalloc.take_snapshot()
                peak_traced_memory = ("{} MiB".format(tracemalloc.get_traced_memory()[1] / 1024 /1024)) #Get the current size and peak size of memory blocks traced by the tracemalloc module as a tuple: (current: int, peak: int)
                tracemalloc.stop()
                top_stats = snapshot2.compare_to(snapshot1, 'lineno') # Compute the differences with an old snapshot.
                out_dict[function.__name__]['peakmem'] = peak_traced_memory
                out_dict[function.__name__]['memleaks'] = top_stats[:10] #
            else:
                result = function(*args, **kwargs)
                out_dict[function.__name__]['peakmem'] = "Unknown"
                out_dict[function.__name__]['memleaks'] = "Unknown" #
            return result
        return function_mem
    return mem_decorator

def mem_use_deco(out_dict):
    def mem_decorator(function):
        @wraps(function)
        def function_mem(*args, **kwargs):
            out = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())], stdout=subprocess.PIPE).communicate()[0].split(b'\n')
            vsz_index = out[0].split().index(b'RSS')
            out_start = float(out[1].split()[vsz_index]) / 1024
            result = function(*args, **kwargs)
            out = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())], stdout=subprocess.PIPE).communicate()[0].split(b'\n')
            vsz_index = out[0].split().index(b'RSS')
            out_end = float(out[1].split()[vsz_index]) / 1024
            out_dict[function.__name__]['Mem Use'] = "{} MiB".format(out_end-out_start)
            return result
        return function_mem
    return mem_decorator

def stats_dict(out_dict):
    def stats_decorator(function):
        @time_execution(out_dict)
        #@cpu_usage(out_dict)
        #@peakmem(out_dict)
        @mem_use_deco(out_dict)
        @wraps(function)
        def all_wrapped(*args, **kwargs):
            return function(*args, **kwargs)
        return all_wrapped
    return stats_decorator

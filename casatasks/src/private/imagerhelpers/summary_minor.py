import copy
import numpy as np
import os

from casatasks.private.casa_transition import is_CASA6
from casampi.MPIEnvironment import MPIEnvironment
if is_CASA6:
    from casatasks import casalog
else:
    from taskinit import *

class SummaryMinor(dict):
    """A return object from tclean's minor cycles. Gathers the information together
    in a way that makes it easier to query for exactly what you want.

    The structure for this dictionary is:
        {
            channel id: {
                polarity id: {
                    summary key: {
                        cycle: value
                    }
                }
            }
        }

    Examples:
    
        1. To get the number of iterations done on the channel 5, first polarity, during the second minor cycle:
            chan5iters = summ_min[5][summ_min.pol0]['iterDone'][1]
    
        2. To get the number of available channels, and the ids of those channels:
            nchans = len(summ_min)
            avail_chans = summ_min.keys()

        3. To get the available minor cycle summary statistics:
            summaryKeys = summ_min.getRowDescriptions()

    There are also a few additional attributes/methods:
        chan0, pol0: the ids of the first channel/polarity
        nCycles, nChan, nPol, nFields: number of minor cycles/channels/polarities/outlier fields logged
        fieldIds: ids of the available outlier fields
        getMatrix(fieldId): returns the original numpy.ndarray matrix
        getAll(summary_key): returns the list of all values for the given summary key
        getDict(calc_iterdone_deltas, keep_startvals): to get the iterDone stat for iterations across all channels
        getOutlierField(fieldid): to get a new SummaryMinor instance the contents of an outlier field other than field 0

    Extends the python dictionary interface (try the command "help(dict)" for more information on builtin python dicts).
    
    Note that this class will be saved as a plain python dict when saved with methods like pickle.dump() or numpy.save().
    This is done to prevent issues when loading later, when this class might not be available to the python interpretter."""
    #                           0           1          2            3              4           5       6      7                  8                9               10                11 "No Mask"      12           13         14          15         16         17       18
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "mapperId", "chan", "pol", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "mpiServer", "peakMem", "runtime", "immod", "stopCode"]
    rowDescriptions13683    = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "immod",    "chan"]
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "mapperId", "cycleStartIters", "masksum", "mpiServer", "peakMem", "runtime", "immod", "stopCode", "chan", "pol"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def __init__(self, summaryminor_matrix, summaryminor_dict = None):
        self.summaryminor_matrix = summaryminor_matrix
        self.singleFieldMatrix = summaryminor_matrix
        self.fieldIds = SummaryMinor._getFieldIds(summaryminor_matrix)
        if len(self.fieldIds) > 1:
            self.singleFieldMatrix = SummaryMinor._getSingleFieldMatrix(summaryminor_matrix, self.fieldIds[0])
        if (len(self.fieldIds) > 1):
            for fieldId in self.fieldIds:
                testMatrix = SummaryMinor._getSingleFieldMatrix(summaryminor_matrix, fieldId)
                sm = SummaryMinor(testMatrix)

        if (summaryminor_dict == None):
            summaryminor_dict = SummaryMinor.indexMinorCycleSummaryBySubimage(self.singleFieldMatrix)
        self.summaryminor_dict = summaryminor_dict
        percycleiters_dict = SummaryMinor._getPerCycleDict(copy.deepcopy(summaryminor_dict))
        self.update(percycleiters_dict)

        self.chan0  = None
        self.pol0   = None
        self.nCycles = 0
        self.nChan   = 0
        self.nPol    = 0
        self.nFields = len( self.fieldIds )
        if len(self.keys()) > 0:
            self.chan0 = list( self.keys() )[0]
            self.nChan = len(  self.keys() )
            if len(self[self.chan0].keys()) > 0:
                self.pol0    = list( self[self.chan0].keys() )[0]
                self.nCycles = len(  self[self.chan0][self.pol0]["iterDone"])
                self.nPol    = len(  self[self.chan0].keys() )

    def _getFieldIds(matrix):
        """Get a sorted list of available outlier field ids in the given matrix"""
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "immod" in availRows:
            return [0]
        immodIdx = availRows.index("immod")
        nrows = matrix.shape[0]
        ncols = matrix.shape[1]
        fieldIds = sorted(np.unique(matrix[immodIdx,:]).tolist())
        return fieldIds

    def _getSingleFieldMatrix(matrixIn, fieldId):
        """Create a new matrix to hold all the values of the given matrix, but only for the given outlier field id"""
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "immod" in availRows:
            return matrixIn
        immodIdx = availRows.index("immod")
        nrowsIn = matrixIn.shape[0]
        ncolsIn = matrixIn.shape[1]
        nrowsOut = nrowsIn
        ncolsOut = matrixIn[immodIdx,:].tolist().count(fieldId)

        matrixOut = np.zeros((nrowsOut, ncolsOut))
        colOut = 0
        maxColOut = 0
        maxRowOut = 0
        for colIn in range(ncolsIn):
            if matrixIn[immodIdx,colIn] != fieldId:
                continue
            for rowIn in range(nrowsIn):
                rowOut = rowIn
                matrixOut[rowOut,colOut] = matrixIn[rowIn,colIn]
                maxRowOut = max(rowOut, maxRowOut)
            maxColOut = colOut
            colOut += 1

        return matrixOut

    def getMatrix(self, fieldId=-1):
        """Returns the original numpy.ndarray matrix.
        fieldId: None for the original matrix, -1 for the matrix used to populate this dictionary, or choose another outlier fields
        Index 0: row (see this.rowDescriptionsOldOrder/rowDescriptions13683)
        Index 1: values for all the minor cycles and outlier fields"""
        if (fieldId == None):
            return self.summaryminor_matrix
        if (fieldId == -1):
            return self.singleFieldMatrix
        return SummaryMinor._getSingleFieldMatrix(self.summaryminor_matrix, fieldId)

    def getAll(self, summary_key):
        """Return the numpy matrix of all values for the given key (eg "chan")"""
        idx = self.getRowDescriptionsOldOrder().index(summary_key)
        return self.getMatrix()[idx,:]

    def useSmallSummaryminor(ignored_parameter=None):
        """Temporary CAS-13683 workaround"""
        if ('USE_SMALL_SUMMARYMINOR' in os.environ):
            uss = os.environ['USE_SMALL_SUMMARYMINOR'].lower()
            if (uss == "true"):
                return True
        return False

    def _getRowDescriptionsOldOrder(useSmallSummaryminor):
        """Temporary CAS-13683 workaround"""
        if (useSmallSummaryminor):
            return SummaryMinor.rowDescriptions13683
        return SummaryMinor.rowDescriptionsOldOrder

    def getRowDescriptionsOldOrder():
        return SummaryMinor._getRowDescriptionsOldOrder(SummaryMinor.useSmallSummaryminor())

    def _getRowDescriptions(useSmallSummaryminor):
        """Temporary CAS-13683 workaround"""
        ret = SummaryMinor.rowDescriptions
        availRows = SummaryMinor._getRowDescriptionsOldOrder(useSmallSummaryminor)
        ret = list(filter(lambda x: x in availRows, ret))
        return ret

    def getRowDescriptions():
        return SummaryMinor._getRowDescriptions(SummaryMinor.useSmallSummaryminor())

    def _getRowStartDescs(useSmallSummaryminor):
        """Temporary CAS-13683 workaround"""
        ret = SummaryMinor.rowStartDescs
        availRows = SummaryMinor._getRowDescriptionsOldOrder(useSmallSummaryminor)
        ret = list(filter(lambda x: x in availRows, ret))
        return ret

    def getRowStartDescs():
        return SummaryMinor._getRowStartDescs(SummaryMinor.useSmallSummaryminor())

    def indexMinorCycleSummaryBySubimage(matrix):
        """Re-indexes matrix from [row,column] to [channel,polarity,row,cycle]."""
        # get some properties of the summary_minor matrix
        nrows = matrix.shape[0]
        ncols = matrix.shape[1]
        uss = SummaryMinor.useSmallSummaryminor() # Temporary CAS-13683 workaround
        import sys
        oldChanIdx = SummaryMinor._getRowDescriptionsOldOrder(uss).index("chan")
        newChanIdx = SummaryMinor._getRowDescriptions(uss).index("chan")
        if not uss:
            oldPolIdx  = SummaryMinor._getRowDescriptionsOldOrder(uss).index("pol")
            newPolIdx  = SummaryMinor._getRowDescriptions(uss).index("pol")
        chans = list(np.sort(np.unique(matrix[oldChanIdx])))
        chans = [int(x) for x in chans]
        if uss:
            pols = [0]
        else:
            pols = list(np.sort(np.unique(matrix[oldPolIdx])))
            pols = [int(x) for x in pols]
        ncycles = 0
        if len(chans) > 0 and len(pols) > 0:
            ncycles = int( ncols / (len(chans)*len(pols)) )
            if uss and MPIEnvironment.is_mpi_enabled:
                # This is necessary because we may have an odd number of "channels" due to each process getting only a subchunk.
                # Example:
                #     Process 1 gets polarities 0-1, process 2 gets polarity 2
                #     Each of them assigns channel id = chan + pol * nsubpols
                #     Process 1 assigns channel ids [0,2], Process 2 assigns channel id 0.
                # This hack is not needed when not using a small summary minor because we have the extra knowledge of the polarities, instead of mapping polarities + channels onto chunks.
                chanslist = matrix[oldChanIdx].tolist()
                for chan in chans:
                    singlechan_occurances = list(filter(lambda x: x == chan, chanslist))
                    ncycles = max(ncycles, len(singlechan_occurances))

        # ret is the return dictionary[chans][pols][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/polarity/row
        ret = {desc:[0]*ncycles for desc in SummaryMinor._getRowDescriptions(uss)}
        # channel and polarity information is in the indexing, don't need to add it to the return dict
        for desc in ["chan", "pol"]:
            if desc in ret:
                del ret[desc]
        ret = {pol:copy.deepcopy(ret) for pol in pols}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure

        # reindex based on subimage index (aka chan/pol index)
        for desc in SummaryMinor._getRowDescriptions(uss):
            if desc in ["chan", "pol"]:
                # channel and polarity information is in the indexing, don't need to add it to the return dict
                continue
            oldRowIdx = SummaryMinor._getRowDescriptionsOldOrder(uss).index(desc)
            for colIdx in range(ncols):
                chan = int(matrix[oldChanIdx][colIdx])
                if (uss):
                    pol = 0
                else:
                    pol = int(matrix[oldPolIdx][colIdx])
                val = matrix[oldRowIdx][colIdx]
                cummulativeCol = int(cummulativeCnt[chan][pol][desc][0]) # const 0: cummulativeCnt doesn't make use of 'cycle' index from copied ret structure
                ret[chan][pol][desc][cummulativeCol] = val
                cummulativeCnt[chan][pol][desc][0] += 1

        return ret

    def _getPerCycleDict(summaryminor_dict, calc_iterdone_deltas=None, keep_startvals=None):
        calc_iterdone_deltas = True if (calc_iterdone_deltas == None) else calc_iterdone_deltas
        keep_startvals       = True if (keep_startvals == None)       else keep_startvals
        ret = summaryminor_dict
        uss = SummaryMinor.useSmallSummaryminor() # Temporary CAS-13683 workaround
        availRows = SummaryMinor._getRowDescriptionsOldOrder(uss)

        if (calc_iterdone_deltas) and ("startIterDone" in availRows):
            for chan in ret:
                for pol in ret[chan]:
                    for cyc in range(len(ret[chan][pol]["startIterDone"])):
                        ret[chan][pol]["iterDone"][cyc] -= ret[chan][pol]["startIterDone"][cyc]
        if not keep_startvals:
            for chan in ret:
                for pol in ret[chan]:
                    for desc in SummaryMinor._getRowStartDescs(uss):
                        del ret[chan][pol][desc]

        return ret

    def getOutlierField(self, fieldId):
        """Get a new SummaryMinor instance for the given fieldId.

        Example:
        nchan_tot = 0
        for field_id in summ_min.fieldIds:
            tmp_summ_min = summ_min.getOutlierField(field_id)
            nchan_tot += tmp_summ_min.nChan"""
        singleFieldMatrix = SummaryMinor._getSingleFieldMatrix(self.summaryminor_matrix, fieldId)
        return SummaryMinor(singleFieldMatrix)

    def getDict(self, calc_iterdone_deltas=None, keep_startvals=None):
        """Computes the per-minor-cycle values for iterDone.

        calc_iterdone_deltas: replaces the original "iterDone" value with the iterations done per cycle
        keep_startvals: don't toss out the start* statistics
        """
        ret = SummaryMinor(self.singleFieldMatrix, self.summaryminor_dict)
        ret = SummaryMinor._getPerCycleDict(ret, calc_iterdone_deltas, keep_startvals)
        return ret

    # pickle as a standard dictionary instead of this class, in case we unpickle and this class isn't available (aka not on the importable PYTHONPATH path)
    def __reduce__(self):
        #       class  init args  state  list vals  dict vals
        return (dict,  tuple(),   None,  None,      self.items().__iter__())

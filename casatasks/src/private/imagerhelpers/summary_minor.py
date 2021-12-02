import copy
import numpy as np
import os

from casatasks.private.casa_transition import is_CASA6
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
    
        1. To get the number of iterations done on the channel 5, polarity 0 during the second minor cycle:
            chan5iters = summMin[5][0]['iterDone'][1]
    
        2. To get the number of available channels, and the ids of those channels:
            nchans = len(summMin)
            avail_chans = summMin.keys()

        3. To get the available minor cycle summary statistics:
            summaryKeys = summMin.rowDescriptions

    There are also a few additional methods:
        getMatrix(): returns the original numpy.ndarray matrix
        getDict(calc_iterdone_deltas, keep_startvals): to get the iterDone stat for iterations across all channels

    Extends the python dictionary interface (try the command "help(dict)" for more information on builtin python dicts).
    
    Note that this class will be saved as a plain python dict when saved with methods like pickle.dump() or numpy.save().
    This is done to prevent issues when loading later, when this class might not be available to the python interpretter."""
    #                           0           1          2            3              4           5       6      7                  8                9               10                11 "No Mask"      12           13         14          15         16         17
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "mapperId", "chan", "pol", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "mpiServer", "peakMem", "runtime", "stopCode"]
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "mapperId", "cycleStartIters", "masksum", "mpiServer", "peakMem", "runtime", "stopCode", "chan", "pol"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def __init__(self, summaryminor_matrix, summaryminor_dict = None):
        self.summaryminor_matrix = summaryminor_matrix
        if (summaryminor_dict == None):
            summaryminor_dict = SummaryMinor.indexMinorCycleSummaryBySubimage(self.summaryminor_matrix)
        self.summaryminor_dict = summaryminor_dict
        percycleiters_dict = SummaryMinor._getPerCycleDict(copy.deepcopy(summaryminor_dict), nrows=summaryminor_matrix.shape[0])
        self.update(percycleiters_dict)

    def getMatrix(self):
        """Returns the original numpy.ndarray matrix.
        Index 0: row (see this.rowDescriptionsOldOrder)
        Index 1: values for all the minor cycles"""
        return self.summaryminor_matrix

    def useSmallSummaryminor():
        """Temporary CAS-13683 workaround"""
        if ('USE_SMALL_SUMMARYMINOR' in os.environ):
            uss = os.environ['USE_SMALL_SUMMARYMINOR'].lower()
            if (uss == "true"):
                return True
        return False

    def getRowDescriptionsOldOrder(nrows):
        """Temporary CAS-13683 workaround"""
        return SummaryMinor.rowDescriptionsOldOrder[0:nrows]

    def getRowDescriptions(nrows):
        """Temporary CAS-13683 workaround"""
        ret = SummaryMinor.rowDescriptions
        rowDescriptionsOldOrder = SummaryMinor.getRowDescriptionsOldOrder(nrows)
        ret = list(filter(lambda x: x in rowDescriptionsOldOrder, ret))
        return ret

    def getRowStartDescs(nrows):
        """Temporary CAS-13683 workaround"""
        ret = SummaryMinor.rowStartDescs
        rowDescriptionsOldOrder = SummaryMinor.getRowDescriptionsOldOrder(nrows)
        ret = list(filter(lambda x: x in rowDescriptionsOldOrder, ret))
        return ret

    def indexMinorCycleSummaryBySubimage(summaryminor):
        """Re-indexes summaryminor from [row,column] to [channel,polarity,row,cycle]."""
        # get some properties of the summaryminor matrix
        nrows = summaryminor.shape[0]
        ncols = summaryminor.shape[1]
        uss = SummaryMinor.useSmallSummaryminor() # Temporary CAS-13683 workaround
        import sys
        oldChanIdx = SummaryMinor.getRowDescriptionsOldOrder(nrows).index("chan")
        newChanIdx = SummaryMinor.getRowDescriptions(nrows).index("chan")
        if not uss:
            oldPolIdx  = SummaryMinor.getRowDescriptionsOldOrder(nrows).index("pol")
            newPolIdx  = SummaryMinor.getRowDescriptions(nrows).index("pol")
        chans = list(np.sort(np.unique(summaryminor[oldChanIdx])))
        chans = [int(x) for x in chans]
        if uss:
            pols = [0]
        else:
            pols = list(np.sort(np.unique(summaryminor[oldPolIdx])))
            pols = [int(x) for x in pols]
        ncycles = int( ncols / (len(chans)*len(pols)) )

        # ret is the return dictionary[chans][pols][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/polarity/row
        ret = {desc:[0]*ncycles for desc in SummaryMinor.getRowDescriptions(nrows)}
        # channel and polarity information is in the indexing, don't need to add it to the return dict
        for desc in ["chan", "pol"]:
            if desc in ret:
                del ret[desc]
        ret = {pol:copy.deepcopy(ret) for pol in pols}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure

        # reindex based on subimage index (aka chan/pol index)
        for desc in SummaryMinor.getRowDescriptions(nrows):
            if desc in ["chan", "pol"]:
                # channel and polarity information is in the indexing, don't need to add it to the return dict
                continue
            oldRowIdx = SummaryMinor.getRowDescriptionsOldOrder(nrows).index(desc)
            for colIdx in range(ncols):
                chan = int(summaryminor[oldChanIdx][colIdx])
                if (uss):
                    pol = 0
                else:
                    pol = int(summaryminor[oldPolIdx][colIdx])
                val = summaryminor[oldRowIdx][colIdx]
                cummulativeCol = int(cummulativeCnt[chan][pol][desc][0]) # cummulativeCnt doesn't make use of last index
                ret[chan][pol][desc][cummulativeCol] = val
                cummulativeCnt[chan][pol][desc][0] += 1

        return ret

    def _getPerCycleDict(summaryminor_dict, calc_iterdone_deltas=None, keep_startvals=None, nrows=0):
        calc_iterdone_deltas = True if (calc_iterdone_deltas == None) else calc_iterdone_deltas
        keep_startvals       = True if (keep_startvals == None)       else keep_startvals
        ret = summaryminor_dict

        # Temporary CAS-13683 workaround
        uss = SummaryMinor.useSmallSummaryminor()
        rowDescriptionsOldOrder = SummaryMinor.getRowDescriptionsOldOrder(nrows)

        if (calc_iterdone_deltas) and ("startIterDone" in rowDescriptionsOldOrder):
            for chan in ret:
                for pol in ret[chan]:
                    for cyc in range(len(ret[chan][pol]["startIterDone"])):
                        ret[chan][pol]["iterDone"][cyc] -= ret[chan][pol]["startIterDone"][cyc]
        if not keep_startvals:
            for chan in ret:
                for pol in ret[chan]:
                    for desc in SummaryMinor.getRowStartDescs(nrows):
                        del ret[chan][pol][desc]

        return ret

    def getMatrix(self):
        """Returns the numpy.ndarray representation of the minor cycle summary."""
        return self.summaryminor_matrix

    def getDict(self, calc_iterdone_deltas=None, keep_startvals=None):
        """Computes the per-minor-cycle values for iterDone.

        calc_iterdone_deltas: replaces the original "iterDone" value with the iterations done per cycle
        keep_startvals: don't toss out the start* statistics
        """
        ret = SummaryMinor(self.summaryminor_matrix, self.summaryminor_dict)
        ret = SummaryMinor._getPerCycleDict(ret, calc_iterdone_deltas, keep_startvals)
        return ret

    # pickle as a standard dictionary instead of this class, in case we unpickle and this class isn't available (aka not on the importable PYTHONPATH path)
    def __reduce__(self):
        #       class  init args  state  list vals  dict vals
        return (dict,  tuple(),   None,  None,      self.items().__iter__())

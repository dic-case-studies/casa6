import copy
import numpy as np

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
    
        1. To get the number of iterations done on the channel 5 during the first minor cycle:
            chan5iters = ret[5][0][0]['iterDone']
    
        2. To get the number of available channels, and the ids of those channels:
            nchans = len(ret)
            avail_chans = ret.keys()

        3. To get the available minor cycle summary statistics:
            statkeys = ret.rowDescriptions

    There are also a few additional methods:
        getMatrix(): returns the original numpy.ndarray matrix
        getDict(calc_iterdone_deltas, keep_startvals): to get the iterDone stat for iterations across all channels

    Extends the python dictionary interface (try the command "help(dict)" for more information on builtin python dicts)."""
    #                           0           1          2            3              4           5       6      7                  8                9               10                11 "No Mask"      12           13         14          15         16         17
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "mapperId", "chan", "pol", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "mpiServer", "peakMem", "runtime", "stopCode"]
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "mapperId", "cycleStartIters", "masksum", "mpiServer", "peakMem", "runtime", "stopCode", "chan", "pol"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def __init__(self, summaryminor_matrix, summaryminor_dict = None):
        self.summaryminor_matrix = summaryminor_matrix
        if (summaryminor_dict == None):
            summaryminor_dict = SummaryMinor.indexMinorCycleSummaryBySubimage(self.summaryminor_matrix)
        self.summaryminor_dict = summaryminor_dict
        percycleiters_dict = SummaryMinor._getPerCycleDict(copy.deepcopy(summaryminor_dict))
        self.update(percycleiters_dict)

    def getMatrix(self):
        """Returns the original numpy.ndarray matrix.
        Index 0: row (see this.rowDescriptionsOldOrder)
        Index 1: values for all the minor cycles"""
        return self.summaryminor_matrix

    def indexMinorCycleSummaryBySubimage(summaryminor):
        """Re-indexes summaryminor from [row,column] to [channel,polarity,row,cycle]."""
        # get some properties of the summaryminor matrix
        oldChanIdx = SummaryMinor.rowDescriptionsOldOrder.index("chan")
        oldPolIdx  = SummaryMinor.rowDescriptionsOldOrder.index("pol")
        newChanIdx = SummaryMinor.rowDescriptionsOldOrder.index("chan")
        newPolIdx  = SummaryMinor.rowDescriptionsOldOrder.index("pol")
        nrows = summaryminor.shape[0]
        ncols = summaryminor.shape[1]
        chans = list(np.sort(np.unique(summaryminor[oldChanIdx])))
        chans = [int(x) for x in chans]
        pols = list(np.sort(np.unique(summaryminor[oldPolIdx])))
        pols = [int(x) for x in pols]
        ncycles = int( ncols / (len(chans)*len(pols)) )

        # ret is the return dictionary[chans][pols][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/polarity/row
        ret = {desc:[0]*ncycles for desc in SummaryMinor.rowDescriptions}
        # channel and polarity information is in the indexing, don't need to add it to the return dict
        for desc in ["chan", "pol"]:
            del ret[desc]
        ret = {pol:copy.deepcopy(ret) for pol in pols}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure

        # reindex based on subimage index (aka chan/pol index)
        for desc in SummaryMinor.rowDescriptions:
            if desc in ["chan", "pol"]:
                # channel and polarity information is in the indexing, don't need to add it to the return dict
                continue
            oldRowIdx = SummaryMinor.rowDescriptionsOldOrder.index(desc)
            for colIdx in range(ncols):
                chan = int(summaryminor[oldChanIdx][colIdx])
                pol = int(summaryminor[oldPolIdx][colIdx])
                val = summaryminor[oldRowIdx][colIdx]
                cummulativeCol = int(cummulativeCnt[chan][pol][desc][0]) # cummulativeCnt doesn't make use of last index
                ret[chan][pol][desc][cummulativeCol] = val
                cummulativeCnt[chan][pol][desc][0] += 1

        return ret

    def _getPerCycleDict(summaryminor_dict, calc_iterdone_deltas=None, keep_startvals=None):
        calc_iterdone_deltas = True if (calc_iterdone_deltas == None) else calc_iterdone_deltas
        keep_startvals       = True if (keep_startvals == None)       else keep_startvals
        ret = summaryminor_dict

        if calc_iterdone_deltas:
            for chan in ret:
                for pol in ret[chan]:
                    for cyc in range(len(ret[chan][pol]["startIterDone"])):
                        ret[chan][pol]["iterDone"][cyc] -= ret[chan][pol]["startIterDone"][cyc]
        if not keep_startvals:
            for chan in ret:
                for pol in ret[chan]:
                    for desc in SummaryMinor.rowStartDescs:
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

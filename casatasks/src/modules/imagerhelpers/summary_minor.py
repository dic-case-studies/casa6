import copy
import numpy as np

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
else:
    from taskinit import *

'''
A return object from tclean's minor cycles. Gathers the information together
in a way that makes it easier to query for exactly what you want.
'''

class SummaryMinor:
    #                           0           1          2            3              4           5       6      7                  8                9               10                11 "No Mask"      12           13         14
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "mapperId", "chan", "pol", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "stopCode"]
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "mapperId", "cycleStartIters", "masksum", "stopCode", "chan", "pol"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def __init__(self, summaryminor_matrix):
        self.summaryminor_matrix = summaryminor_matrix
        self.summaryminor_reorg = None

    def getMatrix(self):
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

    def getDict(self, calc_iterdone_deltas=True, keep_startvals=False):
        if self.summaryminor_reorg == None:
            self.summaryminor_reorg = SummaryMinor.indexMinorCycleSummaryBySubimage(self.summaryminor_matrix)
        ret = self.summaryminor_reorg

        if (calc_iterdone_deltas) or (not keep_startvals):
            ret = copy.deepcopy(ret)
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

    def __str__(self):
        diffsDict = self.getDict()
        ret = "SummaryMinor obj {\n"
        for chan in diffsDict:
            ret += "    "+str(chan)+": "+str(diffsDict[chan])+",\n"
        ret += "}"
        return ret
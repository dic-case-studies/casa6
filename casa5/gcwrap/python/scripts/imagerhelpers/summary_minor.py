import copy
import numpy as np
from collections import OrderedDict

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
    #                           0           1          2            3              4           5       6      7                  8                9               10                11                    12               13
    rowDescriptionsPreOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "mapperId", "chan", "pol", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNoMask", "peakResNoMask", "stopCode"]
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNoMask", "peakResNoMask", "cycleThresh", "mapperId", "cycleStartIters", "stopCode", "chan", "pol"]
    diffRows                = [("startIterDone", "iterDone"), ("startPeakRes", "peakRes"), ("startModelFlux", "modelFlux"), ("startPeakResNoMask", "peakResNoMask")]

    def __init__(self, summaryminor_matrix):
        self.summaryminor_matrix = summaryminor_matrix
        self.summaryminor_reorg = None
        self.summaryminor_reorg_justdiffs = None

    def getMatrix(self):
        return self.summaryminor_matrix

    def indexMinorCycleSummaryBySubimage(summaryminor):
        """Re-indexes summaryminor from [row,column] to [channel,polarity,row,cycle]."""
        # get some properties of the summaryminor matrix
        nrows = summaryminor.shape[0]
        ncols = summaryminor.shape[1]
        chans = list(np.sort(np.unique(summaryminor[5])))
        chans = [int(x) for x in chans]
        pols = list(np.sort(np.unique(summaryminor[6])))
        pols = [int(x) for x in pols]
        ncycles = int( ncols / (len(chans)*len(pols)) )

        # get the row reordering indices
        oldRowIdxs = []
        for desc in SummaryMinor.rowDescriptions:
            oldRowIdxs.append(SummaryMinor.rowDescriptionsPreOrder.index(desc))

        # reindex based on subimage index (aka chan/pol index)
        # ret is the return dictionary[chans][pols][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/polarity/row
        ret = OrderedDict()
        for desc in SummaryMinor.rowDescriptions:
            if desc in ["chan", "pol"]:
                # channel and polarity information is in the indexing, don't need to add it to the return dict
                continue
            ret[desc] = [0]*ncycles
        ret = {pol:copy.deepcopy(ret) for pol in pols}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure
        for rowIdx in range(nrows):
            if rowIdx in [5, 6]:
                # channel and polarity information is in the indexing, don't need to add it to the return dict
                continue
            oldRowIdx = oldRowIdxs[rowIdx]
            desc = SummaryMinor.rowDescriptions[rowIdx]
            for colIdx in range(ncols):
                chan = int(summaryminor[5][colIdx])
                pol = int(summaryminor[6][colIdx])
                val = summaryminor[oldRowIdx][colIdx]
                cummulativeCol = int(cummulativeCnt[chan][pol][desc][0]) # cummulativeCol doesn't make use of last index
                ret[chan][pol][desc][cummulativeCol] = val
                cummulativeCnt[chan][pol][desc][0] += 1

        return ret

    def getDict(self, just_the_diffs=True):
        if self.summaryminor_reorg == None:
            self.summaryminor_reorg = SummaryMinor.indexMinorCycleSummaryBySubimage(self.summaryminor_matrix)
        ret = self.summaryminor_reorg

        if calc_diffs:
            if self.summaryminor_reorg_justdiffs == None:
                ret = copy.deepcopy(ret)
                for chan in ret:
                    for pol in ret[chan]:
                        for (startDesc, desc) in diffRows:
                            for cyc in ret[chan][pol][startDesc]:
                                ret[chan][pol][desc][cyc] -= ret[chan][pol][startDesc][cyc]
                            del ret[chan][pol][startDesc]
                self.summaryminor_reorg_justdiffs = ret
            ret = self.summaryminor_reorg_justdiffs

        return ret

    def __str__(self):
        return "SummaryMinor obj "+str(self.getDict())
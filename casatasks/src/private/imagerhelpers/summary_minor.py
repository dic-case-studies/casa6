import copy
import numpy as np
import os

from casatasks import casalog

class SummaryMinor:
    """Gathers the information together from the tclean return value in a way that makes it easier to query for exactly what you want.

    The structure for this nested dictionary is:
        {
            multi-field id: {
                "chans": [chan0, chan1, ...],
                "pols":  [pol0, pol1, ...],
                "ncycs": number of major/minor cycles
                channel id: {
                    polarity id: {
                        summary key: {
                            cycle: value
                        }
                    }
                }
            },
            "matrix": original numpy ndarray matrix returned from the cpp code
        }

    Examples:
    
        1. To get the number of iterations done on the first field, fifth channel, first polarity, during the middle minor cycle:
            field0 = summ_min[0]
            chan   = field0['chans'][5] # channel index doesn't necessarily start at 0
            pol    = field0['pols'][0]  # polarity index doesn't necessarily start at 0
            cycle  = field0['ncycs']/2
            chan5iters = field0[field][chan][pol]['iterDone'][cycle]
    
        2. To get the number of available channels, and the ids of those channels:
            nchans = len(summ_min[0]["chans"])
            avail_chans = summ_min[0]["chans"]

        3. To get the available minor cycle summary statistics:
            field0 = summ_min[0]
            summaryKeys = field0[field0['chans'][0]][field0['pols'][0]].keys()
    """
    #                           0           1          2            3              4          5       6      7                  8                9               10                11 "No Mask"      12           13         14           15         16         17              18
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "deconId", "chan", "pol", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "mpiServer", "peakMem", "runtime", "multifieldId", "stopCode"]
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "deconId", "cycleStartIters", "masksum", "mpiServer", "peakMem", "runtime", "multifieldId", "stopCode", "chan", "pol"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def convertMatrix(summaryminor_matrix, calc_iterdone_deltas=None, keep_startvals=None):
        ret = {}

        # edge case: no iterations
        if summaryminor_matrix.shape[1] == 0:
            return { 0: {
                'chans': [],
                'pols': [],
                'ncycs': 0
            }}

        # get individual dictionaries for each field id
        field_ids = SummaryMinor._getFieldIds(summaryminor_matrix)
        if len(field_ids) > 1:
            for fieldId in field_ids:
                singleFieldMatrix = SummaryMinor._getSingleFieldMatrix(summaryminor_matrix, field_ids[0])
                ret[fieldId] = SummaryMinor._convertSingleFieldMatrix(singleFieldMatrix, calc_iterdone_deltas, keep_startvals)
        elif len(field_ids) == 1:
            ret[field_ids[0]] = SummaryMinor._convertSingleFieldMatrix(summaryminor_matrix, calc_iterdone_deltas, keep_startvals)
        else:
            casalog.post("Failed to parse summary minor matrix from tclean. No multifield ids were found.", "SEVERE")
            return summaryminor_matrix

        # insert convenience information
        for field_id in field_ids:
            chan_ids = list( ret[field_id].keys() )
            pol_ids  = []
            nCycles = 0

            if len(chan_ids) > 0:
                chan0  = chan_ids[0]
                pol_ids = list( ret[field_id][chan0].keys() )
                if len(pol_ids) > 0:
                    pol0    = pol_ids[0]
                    nCycles = len( ret[field_id][chan0][pol0]['iterDone'] )

            ret[field_id]['chans'] = chan_ids
            ret[field_id]['pols']  = pol_ids
            ret[field_id]['ncycs'] = nCycles

        return ret

    def _convertSingleFieldMatrix(single_field_matrix, calc_iterdone_deltas=None, keep_startvals=None):
        summaryminor_dict = SummaryMinor.indexMinorCycleSummaryBySubimage(single_field_matrix)
        percycleiters_dict = SummaryMinor._getPerCycleDict(copy.deepcopy(summaryminor_dict), calc_iterdone_deltas, keep_startvals)
        return percycleiters_dict

    def _getFieldIds(matrix):
        """ Get a sorted list of available outlier field ids in the given matrix """
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "multifieldId" in availRows:
            return [0]
        multifieldIdx = availRows.index("multifieldId")
        nrows = matrix.shape[0]
        ncols = matrix.shape[1]
        fieldIds = sorted(np.unique(matrix[multifieldIdx,:]).tolist())
        fieldIds = list(map(lambda x: int(x), fieldIds))
        return fieldIds

    def _getSingleFieldMatrix(matrixIn, fieldId):
        """ Create a new matrix to hold all the values of the given matrix, but only for the given outlier field id """
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "multifieldId" in availRows:
            return matrixIn
        multifieldIdx = availRows.index("multifieldId")
        nrowsIn = matrixIn.shape[0]
        ncolsIn = matrixIn.shape[1]
        nrowsOut = nrowsIn
        ncolsOut = matrixIn[multifieldIdx,:].tolist().count(fieldId)

        matrixOut = np.zeros((nrowsOut, ncolsOut))
        colOut = 0
        maxColOut = 0
        maxRowOut = 0
        for colIn in range(ncolsIn):
            if matrixIn[multifieldIdx,colIn] != fieldId:
                continue
            for rowIn in range(nrowsIn):
                rowOut = rowIn
                matrixOut[rowOut,colOut] = matrixIn[rowIn,colIn]
                maxRowOut = max(rowOut, maxRowOut)
            maxColOut = colOut
            colOut += 1

        return matrixOut

    def getRowDescriptionsOldOrder():
        return SummaryMinor.rowDescriptionsOldOrder

    def getRowDescriptions():
        return SummaryMinor.rowDescriptions

    def getRowStartDescs():
        return SummaryMinor.rowStartDescs

    def indexMinorCycleSummaryBySubimage(matrix):
        """ Re-indexes matrix from [row,column] to [channel,polarity,row,cycle]. 

        Param matrix: the original matrix to convert.
        """
        # get some properties of the summary_minor matrix
        nrows = matrix.shape[0]
        ncols = matrix.shape[1]
        import sys
        oldChanIdx = SummaryMinor.getRowDescriptionsOldOrder().index("chan")
        newChanIdx = SummaryMinor.getRowDescriptions().index("chan")
        oldPolIdx  = SummaryMinor.getRowDescriptionsOldOrder().index("pol")
        newPolIdx  = SummaryMinor.getRowDescriptions().index("pol")
        chans = list(np.sort(np.unique(matrix[oldChanIdx])))
        chans = [int(x) for x in chans]
        pols = list(np.sort(np.unique(matrix[oldPolIdx])))
        pols = [int(x) for x in pols]
        ncycles = 0
        if len(chans) > 0 and len(pols) > 0:
            ncycles = int( ncols / (len(chans)*len(pols)) )

        # ret is the return dictionary[chans][pols][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/polarity/row
        ret = {desc:[0]*ncycles for desc in SummaryMinor.getRowDescriptions()}
        # channel and polarity information is in the indexing, don't need to add it to the return dict
        for desc in ["chan", "pol"]:
            if desc in ret:
                del ret[desc]
        ret = {pol:copy.deepcopy(ret) for pol in pols}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure

        # reindex based on subimage index (aka chan/pol index)
        for desc in SummaryMinor.getRowDescriptions():
            if desc in ["chan", "pol"]:
                # channel and polarity information is in the indexing, don't need to add it to the return dict
                continue
            oldRowIdx = SummaryMinor.getRowDescriptionsOldOrder().index(desc)
            for colIdx in range(ncols):
                chan = int(matrix[oldChanIdx][colIdx])
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
        availRows = SummaryMinor.getRowDescriptionsOldOrder()

        if (calc_iterdone_deltas) and ("startIterDone" in availRows):
            for chan in ret:
                for pol in ret[chan]:
                    for cyc in range(len(ret[chan][pol]["startIterDone"])):
                        ret[chan][pol]["iterDone"][cyc] -= ret[chan][pol]["startIterDone"][cyc]
        if not keep_startvals:
            for chan in ret:
                for pol in ret[chan]:
                    for desc in SummaryMinor.getRowStartDescs():
                        del ret[chan][pol][desc]

        return ret